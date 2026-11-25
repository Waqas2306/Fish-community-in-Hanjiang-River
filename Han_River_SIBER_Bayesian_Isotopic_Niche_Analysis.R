# RUN IN RSTUDIO: open this R file and click Source.
# Input file(s): Han_River_SIBER_Bayesian_Isotopic_Niche_R_Dataset(1).xlsx
# Worksheet tabs read automatically: SIBER_Input; Group_Lookup; Pairwise_Comparisons; Centroid_Movement_Pairs
# Use the corrected files beside this script. See WHICH_EXCEL_TO_USE.md.

# Updated for the supplied 21 November 2026 master and corrected data package.
.hr_paths <- vapply(sys.frames(), function(f) {
  p <- f$ofile
  if (is.character(p) && length(p) == 1L && nzchar(p)) p else ""
}, character(1))
.hr_paths <- .hr_paths[nzchar(.hr_paths)]
.hr_args <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.hr_here <- if (length(.hr_paths)) dirname(normalizePath(tail(.hr_paths, 1L), mustWork = TRUE)) else if (length(.hr_args)) dirname(normalizePath(.hr_args[1L], mustWork = TRUE)) else getwd()
if (!file.exists(file.path(.hr_here, "Han_River_Runtime.R")) && requireNamespace("rstudioapi", quietly = TRUE)) {
  .hr_editor <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
  if (nzchar(.hr_editor)) .hr_here <- dirname(.hr_editor)
}
if (!file.exists(file.path(.hr_here, "Han_River_Runtime.R"))) stop("Extract the whole package, then open START_HERE.R and click Source.")
options(hanriver.project = normalizePath(.hr_here, winslash = "/", mustWork = TRUE))
source(file.path(.hr_here, "Han_River_Runtime.R"), local = TRUE, encoding = "UTF-8")

# =============================================================================
# LOWER HAN RIVER FISH STABLE ISOTOPES
# BAYESIAN ISOTOPIC NICHE ANALYSIS USING SIBER
# FINAL READY-TO-USE R SCRIPT
#
# PURPOSE
#   Compares the six reach–campaign isotope groups using:
#     1. Bayesian Standard Ellipse Area (SEA.B);
#     2. pairwise 95% Bayesian ellipse overlap;
#     3. posterior centroid location;
#     4. pairwise centroid distance; and
#     5. wet-to-dry centroid movement within each reach.
#
# INPUT WORKBOOK
#   Han_River_SIBER_Bayesian_Isotopic_Niche_R_Dataset(1).xlsx
#
# REQUIRED SYSTEM SOFTWARE
#   JAGS 4.x must be installed separately. The R package rjags is only the
#   interface to JAGS.
#
# IMPORTANT DESIGN
#   SIBER community = river reach:
#       1 = S1; 2 = S2; 3 = S3
#
#   SIBER group = campaign nested within reach:
#       1 = Wet_2023; 2 = Dry_2023
#
#   The resulting SIBER codes are:
#       1.1 HC wet; 1.2 HC dry
#       2.1 QJ wet; 2.2 QJ dry
#       3.1 SY wet; 3.2 SY dry
#
# HOW TO RUN
#   1. Install standalone JAGS before running this script.
#   2. Install the required R packages.
#   3. Place the Excel workbook in the same folder as this script, Downloads,
#      Desktop, or Documents. Otherwise select it when prompted.
#   4. Click Source once.
#   5. Do not click Source again while JAGS or overlap calculations are running.
# =============================================================================


# =============================================================================
# 1. USER SETTINGS
# =============================================================================

FORCE_REFIT <- isTRUE(getOption("hanriver.force_refit", FALSE))

SEED <- 20260730L

# Publication site labels. HC, QJ, and SY remain internal data codes.
REACH_LABELS <- c(
  HC = "S1",
  QJ = "S2",
  SY = "S3"
)

REACH_LABEL_LEVELS <- unname(
  REACH_LABELS[c("HC", "QJ", "SY")]
)

# Publication typography in points.
AXIS_TITLE_SIZE <- 16
AXIS_TEXT_SIZE <- 14
FACET_TEXT_SIZE <- 14
LEGEND_TEXT_SIZE <- 13
LEGEND_TITLE_SIZE <- 13

# Official SIBER-style MCMC settings
N_ITER <- 20000L
N_BURNIN <- 1000L
N_THIN <- 10L
N_CHAINS <- 2L

# Bayesian ellipse-overlap calculations are computationally intensive.
OVERLAP_DRAWS <- 500L
OVERLAP_INTERVAL <- 0.95
ELLIPSE_POINTS <- 180L

# Equal-tailed posterior interval
LOWER_PROB <- 0.025
UPPER_PROB <- 0.975

set.seed(SEED)
graphics.off()


# =============================================================================
# 2. CHECK AND LOAD PACKAGES
# =============================================================================

required_packages <- c(
  "SIBER",
  "rjags",
  "coda",
  "readxl",
  "dplyr",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    paste0(
      "Missing R package(s): ",
      paste(missing_packages, collapse = ", "),
      "\nInstall the missing package(s), restart RStudio, and run again."
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(SIBER)
  library(rjags)
  library(coda)
  library(readxl)
  library(dplyr)
  library(ggplot2)
})


# =============================================================================
# 3. VERIFY THE STANDALONE JAGS INSTALLATION
# =============================================================================

jags_version <- tryCatch(
  {
    as.character(
      rjags::jags.version()
    )
  },
  error = function(e) {
    NA_character_
  }
)

if (
  length(jags_version) == 0L ||
  is.na(jags_version) ||
  !nzchar(jags_version)
) {
  stop(
    paste0(
      "JAGS could not be detected.\n",
      "Install the standalone JAGS 4.x software, restart RStudio, ",
      "and run this script again."
    ),
    call. = FALSE
  )
}

cat(
  "\nJAGS detected: version ",
  jags_version,
  "\nSIBER package version: ",
  as.character(
    utils::packageVersion("SIBER")
  ),
  "\n",
  sep = ""
)


# =============================================================================
# 4. LOCATE THE INPUT WORKBOOK
# =============================================================================

INPUT_FILE <- hr_input("siber")
hr_validate_workbook(INPUT_FILE, "siber")


cat(
  "\nInput workbook:\n",
  INPUT_FILE,
  "\n",
  sep = ""
)


# =============================================================================
# 5. CREATE OUTPUT FOLDERS AND STATUS FILES
# =============================================================================

RESULTS_DIR <- hr_results("siber")

MODEL_DIR <- file.path(
  RESULTS_DIR,
  "Model"
)

JAGS_DIR <- file.path(
  MODEL_DIR,
  "Raw_JAGS_output"
)

TABLE_DIR <- file.path(
  RESULTS_DIR,
  "Tables"
)

FIGURE_DIR <- file.path(
  RESULTS_DIR,
  "Figures"
)

LOG_DIR <- file.path(
  RESULTS_DIR,
  "Logs"
)

for (directory_path in c(
  RESULTS_DIR,
  MODEL_DIR,
  JAGS_DIR,
  TABLE_DIR,
  FIGURE_DIR,
  LOG_DIR
)) {
  dir.create(
    directory_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

RUNNING_FILE <- file.path(
  LOG_DIR,
  "ANALYSIS_RUNNING.txt"
)

COMPLETED_FILE <- file.path(
  LOG_DIR,
  "ANALYSIS_COMPLETED.txt"
)

ERROR_FILE <- file.path(
  LOG_DIR,
  "ANALYSIS_ERROR.txt"
)

unlink(
  c(
    COMPLETED_FILE,
    ERROR_FILE
  )
)

writeLines(
  c(
    "SIBER Bayesian isotopic niche analysis is running.",
    paste("Started:", Sys.time()),
    paste("Input workbook:", INPUT_FILE),
    paste("SIBER version:", utils::packageVersion("SIBER")),
    paste("JAGS version:", jags_version),
    paste("Overlap posterior draws:", OVERLAP_DRAWS)
  ),
  RUNNING_FILE
)


# =============================================================================
# 6. READ AND VALIDATE WORKBOOK SHEETS
# =============================================================================

available_sheets <- readxl::excel_sheets(
  INPUT_FILE
)

required_sheets <- c(
  "SIBER_Input",
  "Group_Lookup",
  "Pairwise_Comparisons",
  "Centroid_Movement_Pairs"
)

missing_sheets <- setdiff(
  required_sheets,
  available_sheets
)

if (length(missing_sheets) > 0L) {
  stop(
    paste0(
      "The workbook is missing worksheet(s): ",
      paste(missing_sheets, collapse = ", ")
    ),
    call. = FALSE
  )
}

siber_data <- hr_read_excel(
  INPUT_FILE,
  sheet = "SIBER_Input"
)

siber_data <- as.data.frame(
  siber_data,
  stringsAsFactors = FALSE
)

if (!identical(
  names(siber_data),
  c(
    "iso1",
    "iso2",
    "group",
    "community"
  )
)) {
  stop(
    paste0(
      "SIBER_Input must contain exactly these four columns in this order:\n",
      "iso1, iso2, group, community"
    ),
    call. = FALSE
  )
}

for (column_name in c(
  "iso1",
  "iso2",
  "group",
  "community"
)) {

  siber_data[[column_name]] <- suppressWarnings(
    as.numeric(
      siber_data[[column_name]]
    )
  )
}

if (anyNA(siber_data)) {
  stop(
    "Missing or non-numeric values were detected in SIBER_Input.",
    call. = FALSE
  )
}

if (
  any(
    siber_data$group != round(
      siber_data$group
    )
  ) ||
  any(
    siber_data$community != round(
      siber_data$community
    )
  )
) {
  stop(
    "The group and community variables must contain integer codes.",
    call. = FALSE
  )
}

siber_data$group <- as.integer(
  siber_data$group
)

siber_data$community <- as.integer(
  siber_data$community
)

group_lookup <- hr_read_excel(
  INPUT_FILE,
  sheet = "Group_Lookup"
)

group_lookup <- as.data.frame(
  group_lookup,
  stringsAsFactors = FALSE
)

pairwise_comparisons <- hr_read_excel(
  INPUT_FILE,
  sheet = "Pairwise_Comparisons"
)

pairwise_comparisons <- as.data.frame(
  pairwise_comparisons,
  stringsAsFactors = FALSE
)

centroid_movement_pairs <- hr_read_excel(
  INPUT_FILE,
  sheet = "Centroid_Movement_Pairs"
)

centroid_movement_pairs <- as.data.frame(
  centroid_movement_pairs,
  stringsAsFactors = FALSE
)

required_lookup_columns <- c(
  "SIBER_code",
  "community",
  "group",
  "Reach_code",
  "Reach_name",
  "Campaign",
  "Campaign_label",
  "Group_label",
  "n_fish",
  "species_richness"
)

missing_lookup_columns <- setdiff(
  required_lookup_columns,
  names(group_lookup)
)

if (length(missing_lookup_columns) > 0L) {
  stop(
    paste0(
      "Group_Lookup is missing column(s): ",
      paste(missing_lookup_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

required_pair_columns <- c(
  "Pair_ID",
  "Ellipse_1_code",
  "Ellipse_2_code",
  "Comparison_label",
  "Comparison_type",
  "Manuscript_priority"
)

missing_pair_columns <- setdiff(
  required_pair_columns,
  names(pairwise_comparisons)
)

if (length(missing_pair_columns) > 0L) {
  stop(
    paste0(
      "Pairwise_Comparisons is missing column(s): ",
      paste(missing_pair_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

required_movement_columns <- c(
  "Reach_code",
  "Reach_name",
  "Wet_code",
  "Dry_code"
)

missing_movement_columns <- setdiff(
  required_movement_columns,
  names(centroid_movement_pairs)
)

if (length(missing_movement_columns) > 0L) {
  stop(
    paste0(
      "Centroid_Movement_Pairs is missing column(s): ",
      paste(missing_movement_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

# Force all SIBER identifiers to character form.
group_lookup$SIBER_code <- as.character(
  group_lookup$SIBER_code
)

pairwise_comparisons$Ellipse_1_code <- as.character(
  pairwise_comparisons$Ellipse_1_code
)

pairwise_comparisons$Ellipse_2_code <- as.character(
  pairwise_comparisons$Ellipse_2_code
)

centroid_movement_pairs$Wet_code <- as.character(
  centroid_movement_pairs$Wet_code
)

centroid_movement_pairs$Dry_code <- as.character(
  centroid_movement_pairs$Dry_code
)

# Replace workbook display names with the publication labels everywhere while
# preserving the original reach codes used by SIBER and the statistical logic.
group_lookup$Reach_name <- unname(
  REACH_LABELS[as.character(group_lookup$Reach_code)]
)

centroid_movement_pairs$Reach_name <- unname(
  REACH_LABELS[as.character(centroid_movement_pairs$Reach_code)]
)

if (anyNA(group_lookup$Reach_name) || anyNA(centroid_movement_pairs$Reach_name)) {
  stop(
    "A reach code could not be mapped to S1, S2, or S3.",
    call. = FALSE
  )
}

group_lookup$Group_label <- paste(
  group_lookup$Reach_name,
  group_lookup$Campaign_label,
  sep = " – "
)

expected_codes <- c(
  "1.1",
  "1.2",
  "2.1",
  "2.2",
  "3.1",
  "3.2"
)

observed_codes <- sort(
  unique(
    paste(
      siber_data$community,
      siber_data$group,
      sep = "."
    )
  )
)

if (!identical(
  observed_codes,
  expected_codes
)) {
  stop(
    paste0(
      "The six expected SIBER groups were not found.\n",
      "Observed codes: ",
      paste(observed_codes, collapse = ", ")
    ),
    call. = FALSE
  )
}

group_sample_sizes <- as.data.frame(
  table(
    siber_data$community,
    siber_data$group
  ),
  stringsAsFactors = FALSE
)

names(group_sample_sizes) <- c(
  "community",
  "group",
  "n_fish"
)

group_sample_sizes <- group_sample_sizes[
  group_sample_sizes$n_fish > 0,
  ,
  drop = FALSE
]

if (
  nrow(group_sample_sizes) != 6L ||
  min(group_sample_sizes$n_fish) < 5L
) {
  stop(
    paste0(
      "Each of the six SIBER groups must contain at least five observations."
    ),
    call. = FALSE
  )
}

group_labels <- stats::setNames(
  as.character(
    group_lookup$Group_label
  ),
  group_lookup$SIBER_code
)

short_group_labels <- stats::setNames(
  paste0(
    group_lookup$Reach_name,
    " – ",
    ifelse(
      group_lookup$Campaign == "Wet_2023",
      "Wet",
      "Dry"
    )
  ),
  group_lookup$SIBER_code
)

pairwise_comparisons$Comparison_label <- paste(
  unname(short_group_labels[pairwise_comparisons$Ellipse_1_code]),
  "vs",
  unname(short_group_labels[pairwise_comparisons$Ellipse_2_code])
)

cat(
  "\nData validation completed.",
  "\nObservations: ",
  nrow(siber_data),
  "\nGroups: ",
  length(observed_codes),
  "\nMinimum group sample size: ",
  min(group_sample_sizes$n_fish),
  "\n",
  sep = ""
)

write.csv(
  group_sample_sizes,
  file.path(
    TABLE_DIR,
    "Table_S1_group_sample_sizes.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# =============================================================================
# 7. CREATE THE SIBER OBJECT AND MAXIMUM-LIKELIHOOD CHECKS
# =============================================================================

siber_object <- SIBER::createSiberObject(
  siber_data
)

group_metrics_ml <- SIBER::groupMetricsML(
  siber_object
)

ml_metric_table <- data.frame(
  SIBER_code = colnames(
    group_metrics_ml
  ),
  TA_ML = as.numeric(
    group_metrics_ml[
      "TA",
      ,
      drop = TRUE
    ]
  ),
  SEA_ML = as.numeric(
    group_metrics_ml[
      "SEA",
      ,
      drop = TRUE
    ]
  ),
  SEAc_ML = as.numeric(
    group_metrics_ml[
      "SEAc",
      ,
      drop = TRUE
    ]
  ),
  stringsAsFactors = FALSE
)

ml_metric_table <- merge(
  group_lookup[
    ,
    required_lookup_columns,
    drop = FALSE
  ],
  ml_metric_table,
  by = "SIBER_code",
  all.x = TRUE,
  sort = FALSE
)

write.csv(
  ml_metric_table,
  file.path(
    TABLE_DIR,
    "Table_1_maximum_likelihood_niche_metrics.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# =============================================================================
# 8. FIT OR LOAD THE BAYESIAN SIBER MODEL
# =============================================================================

MODEL_FILE <- file.path(
  MODEL_DIR,
  "SIBER_ellipse_posteriors.rds"
)

parms <- list()
parms$n.iter <- N_ITER
parms$n.burnin <- N_BURNIN
parms$n.thin <- N_THIN
parms$n.chains <- N_CHAINS
parms$save.output <- TRUE
parms$save.dir <- JAGS_DIR

priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

is_valid_siber_posterior <- function(
  object
) {

  if (
    !is.list(object) ||
    length(object) != 6L
  ) {
    return(
      FALSE
    )
  }

  if (!all(
    expected_codes %in% names(object)
  )) {
    return(
      FALSE
    )
  }

  matrix_checks <- vapply(
    expected_codes,
    function(code) {

      current_matrix <- object[[code]]

      is.matrix(current_matrix) &&
      ncol(current_matrix) >= 6L &&
      nrow(current_matrix) > 100L &&
      all(
        is.finite(
          current_matrix
        )
      )
    },
    FUN.VALUE = logical(1)
  )

  all(
    matrix_checks
  )
}

valid_saved_model <- FALSE

if (
  file.exists(MODEL_FILE) &&
  !isTRUE(FORCE_REFIT)
) {

  loaded_model <- try(
    readRDS(MODEL_FILE),
    silent = TRUE
  )

  valid_saved_model <- (
    !inherits(loaded_model, "try-error") &&
    is_valid_siber_posterior(
      loaded_model
    )
  )
}

if (valid_saved_model) {

  ellipses_posterior <- loaded_model

  cat(
    "\nLoaded the existing SIBER posterior object:\n",
    MODEL_FILE,
    "\n",
    sep = ""
  )

} else {

  if (file.exists(MODEL_FILE)) {
    unlink(
      MODEL_FILE
    )
  }

  if (isTRUE(FORCE_REFIT)) {

    old_jags_files <- list.files(
      JAGS_DIR,
      pattern = "^jags_output_.*[.]RData$",
      full.names = TRUE
    )

    if (length(old_jags_files) > 0L) {
      unlink(
        old_jags_files
      )
    }
  }

  cat(
    "\nBeginning SIBER Bayesian ellipse fitting.",
    "\nSix group-specific bivariate models will be fitted using JAGS.",
    "\nDo not run the script again while fitting is in progress.\n",
    sep = ""
  )

  set.seed(
    SEED
  )

  ellipses_posterior <- SIBER::siberMVN(
    siber_object,
    parms,
    priors
  )

  if (!is_valid_siber_posterior(
    ellipses_posterior
  )) {
    stop(
      "The object returned by siberMVN() did not pass validation.",
      call. = FALSE
    )
  }

  saveRDS(
    ellipses_posterior,
    MODEL_FILE
  )
}

posterior_draw_count <- min(
  vapply(
    ellipses_posterior,
    nrow,
    FUN.VALUE = integer(1)
  )
)

cat(
  "\nBayesian ellipse posteriors available.",
  "\nPosterior draws per group: ",
  posterior_draw_count,
  "\n",
  sep = ""
)


# =============================================================================
# 9. MCMC CONVERGENCE AND EFFECTIVE SAMPLE-SIZE DIAGNOSTICS
#
# SIBER can save each group's original mcmc.list object. These files preserve
# the separate chains required for Gelman-Rubin and effective-size diagnostics.
# =============================================================================

diagnostic_parameter_rows <- list()
diagnostic_group_rows <- list()

diagnostic_parameter_index <- 1L
diagnostic_group_index <- 1L

for (lookup_index in seq_len(
  nrow(group_lookup)
)) {

  community_index <- as.integer(
    group_lookup$community[
      lookup_index
    ]
  )

  group_index <- as.integer(
    group_lookup$group[
      lookup_index
    ]
  )

  current_code <- as.character(
    group_lookup$SIBER_code[
      lookup_index
    ]
  )

  current_file <- file.path(
    JAGS_DIR,
    paste0(
      "jags_output_community",
      community_index,
      "_group",
      group_index,
      ".RData"
    )
  )

  if (!file.exists(
    current_file
  )) {

    diagnostic_group_rows[[diagnostic_group_index]] <- data.frame(
      SIBER_code = current_code,
      Group_label = group_labels[[current_code]],
      Maximum_Rhat = NA_real_,
      Maximum_Rhat_upper_CI = NA_real_,
      Minimum_effective_sample_size = NA_real_,
      Diagnostic_status = "Raw chain file not available",
      stringsAsFactors = FALSE
    )

    diagnostic_group_index <- diagnostic_group_index + 1L
    next
  }

  diagnostic_environment <- new.env(
    parent = emptyenv()
  )

  loaded_names <- load(
    current_file,
    envir = diagnostic_environment
  )

  if (!("output" %in% loaded_names)) {

    diagnostic_group_rows[[diagnostic_group_index]] <- data.frame(
      SIBER_code = current_code,
      Group_label = group_labels[[current_code]],
      Maximum_Rhat = NA_real_,
      Maximum_Rhat_upper_CI = NA_real_,
      Minimum_effective_sample_size = NA_real_,
      Diagnostic_status = "No object named output in chain file",
      stringsAsFactors = FALSE
    )

    diagnostic_group_index <- diagnostic_group_index + 1L
    next
  }

  current_output <- diagnostic_environment$output

  gelman_result <- try(
    coda::gelman.diag(
      current_output,
      autoburnin = FALSE,
      multivariate = FALSE
    ),
    silent = TRUE
  )

  effective_sizes <- try(
    coda::effectiveSize(
      current_output
    ),
    silent = TRUE
  )

  if (
    inherits(gelman_result, "try-error") ||
    inherits(effective_sizes, "try-error")
  ) {

    diagnostic_group_rows[[diagnostic_group_index]] <- data.frame(
      SIBER_code = current_code,
      Group_label = group_labels[[current_code]],
      Maximum_Rhat = NA_real_,
      Maximum_Rhat_upper_CI = NA_real_,
      Minimum_effective_sample_size = NA_real_,
      Diagnostic_status = "Diagnostic calculation failed",
      stringsAsFactors = FALSE
    )

    diagnostic_group_index <- diagnostic_group_index + 1L
    next
  }

  psrf_matrix <- gelman_result$psrf

  common_parameters <- intersect(
    rownames(
      psrf_matrix
    ),
    names(
      effective_sizes
    )
  )

  if (length(common_parameters) == 0L) {

    diagnostic_group_rows[[diagnostic_group_index]] <- data.frame(
      SIBER_code = current_code,
      Group_label = group_labels[[current_code]],
      Maximum_Rhat = NA_real_,
      Maximum_Rhat_upper_CI = NA_real_,
      Minimum_effective_sample_size = NA_real_,
      Diagnostic_status = "No common diagnostic parameter names",
      stringsAsFactors = FALSE
    )

    diagnostic_group_index <- diagnostic_group_index + 1L
    next
  }

  current_parameter_table <- data.frame(
    SIBER_code = current_code,
    Group_label = group_labels[[current_code]],
    Parameter = common_parameters,
    Rhat = as.numeric(
      psrf_matrix[
        common_parameters,
        "Point est."
      ]
    ),
    Rhat_upper_CI = as.numeric(
      psrf_matrix[
        common_parameters,
        "Upper C.I."
      ]
    ),
    Effective_sample_size = as.numeric(
      effective_sizes[
        common_parameters
      ]
    ),
    stringsAsFactors = FALSE
  )

  diagnostic_parameter_rows[[diagnostic_parameter_index]] <- current_parameter_table
  diagnostic_parameter_index <- diagnostic_parameter_index + 1L

  diagnostic_group_rows[[diagnostic_group_index]] <- data.frame(
    SIBER_code = current_code,
    Group_label = group_labels[[current_code]],
    Maximum_Rhat = max(
      current_parameter_table$Rhat,
      na.rm = TRUE
    ),
    Maximum_Rhat_upper_CI = max(
      current_parameter_table$Rhat_upper_CI,
      na.rm = TRUE
    ),
    Minimum_effective_sample_size = min(
      current_parameter_table$Effective_sample_size,
      na.rm = TRUE
    ),
    Diagnostic_status = "Calculated",
    stringsAsFactors = FALSE
  )

  diagnostic_group_index <- diagnostic_group_index + 1L
}

diagnostic_parameter_table <- if (
  length(diagnostic_parameter_rows) > 0L
) {
  dplyr::bind_rows(
    diagnostic_parameter_rows
  )
} else {
  data.frame()
}

diagnostic_group_table <- dplyr::bind_rows(
  diagnostic_group_rows
)

write.csv(
  diagnostic_group_table,
  file.path(
    TABLE_DIR,
    "Table_2_group_MCMC_diagnostics.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

if (nrow(
  diagnostic_parameter_table
) > 0L) {

  write.csv(
    diagnostic_parameter_table,
    file.path(
      TABLE_DIR,
      "Table_S2_parameter_MCMC_diagnostics.csv"
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}


# =============================================================================
# 10. BAYESIAN STANDARD ELLIPSE AREA
# =============================================================================

SEA_B <- SIBER::siberEllipses(
  ellipses_posterior
)

colnames(
  SEA_B
) <- names(
  ellipses_posterior
)

summarise_numeric <- function(
  values
) {

  values <- as.numeric(
    values
  )

  values <- values[
    is.finite(
      values
    )
  ]

  if (length(values) == 0L) {
    return(
      c(
        Mean = NA_real_,
        Median = NA_real_,
        Lower_95 = NA_real_,
        Upper_95 = NA_real_
      )
    )
  }

  c(
    Mean = mean(
      values
    ),
    Median = stats::median(
      values
    ),
    Lower_95 = unname(
      stats::quantile(
        values,
        probs = LOWER_PROB
      )
    ),
    Upper_95 = unname(
      stats::quantile(
        values,
        probs = UPPER_PROB
      )
    )
  )
}

sea_summary_rows <- vector(
  mode = "list",
  length = length(
    expected_codes
  )
)

for (code_index in seq_along(
  expected_codes
)) {

  current_code <- expected_codes[
    code_index
  ]

  current_summary <- summarise_numeric(
    SEA_B[
      ,
      current_code
    ]
  )

  sea_summary_rows[[code_index]] <- data.frame(
    SIBER_code = current_code,
    Group_label = group_labels[[current_code]],
    Posterior_mean_SEA = current_summary[["Mean"]],
    Posterior_median_SEA = current_summary[["Median"]],
    Lower_95 = current_summary[["Lower_95"]],
    Upper_95 = current_summary[["Upper_95"]],
    stringsAsFactors = FALSE
  )
}

sea_summary <- dplyr::bind_rows(
  sea_summary_rows
)

sea_summary <- merge(
  sea_summary,
  ml_metric_table[
    ,
    c(
      "SIBER_code",
      "Reach_code",
      "Reach_name",
      "Campaign",
      "Campaign_label",
      "n_fish",
      "species_richness",
      "SEA_ML",
      "SEAc_ML"
    ),
    drop = FALSE
  ],
  by = "SIBER_code",
  all.x = TRUE,
  sort = FALSE
)

sea_summary <- sea_summary[
  match(
    expected_codes,
    sea_summary$SIBER_code
  ),
  ,
  drop = FALSE
]

write.csv(
  sea_summary,
  file.path(
    TABLE_DIR,
    "Table_3_Bayesian_SEA_summary.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# Pairwise posterior comparisons of niche width
sea_pairwise_rows <- vector(
  mode = "list",
  length = nrow(
    pairwise_comparisons
  )
)

for (pair_index in seq_len(
  nrow(
    pairwise_comparisons
  )
)) {

  code_1 <- pairwise_comparisons$Ellipse_1_code[
    pair_index
  ]

  code_2 <- pairwise_comparisons$Ellipse_2_code[
    pair_index
  ]

  n_draws <- min(
    nrow(
      SEA_B
    ),
    length(
      SEA_B[
        ,
        code_1
      ]
    ),
    length(
      SEA_B[
        ,
        code_2
      ]
    )
  )

  difference_draws <- (
    SEA_B[
      seq_len(
        n_draws
      ),
      code_2
    ] -
    SEA_B[
      seq_len(
        n_draws
      ),
      code_1
    ]
  )

  ratio_draws <- (
    SEA_B[
      seq_len(
        n_draws
      ),
      code_2
    ] /
    SEA_B[
      seq_len(
        n_draws
      ),
      code_1
    ]
  )

  difference_summary <- summarise_numeric(
    difference_draws
  )

  ratio_summary <- summarise_numeric(
    ratio_draws
  )

  sea_pairwise_rows[[pair_index]] <- data.frame(
    Pair_ID = pairwise_comparisons$Pair_ID[
      pair_index
    ],
    Ellipse_1_code = code_1,
    Ellipse_2_code = code_2,
    Comparison_label = pairwise_comparisons$Comparison_label[
      pair_index
    ],
    Comparison_type = pairwise_comparisons$Comparison_type[
      pair_index
    ],
    Manuscript_priority = pairwise_comparisons$Manuscript_priority[
      pair_index
    ],
    Median_SEA_difference_2_minus_1 = difference_summary[["Median"]],
    Lower_95_difference = difference_summary[["Lower_95"]],
    Upper_95_difference = difference_summary[["Upper_95"]],
    Probability_SEA_2_greater_1 = mean(
      difference_draws > 0
    ),
    Median_SEA_ratio_2_over_1 = ratio_summary[["Median"]],
    Lower_95_ratio = ratio_summary[["Lower_95"]],
    Upper_95_ratio = ratio_summary[["Upper_95"]],
    stringsAsFactors = FALSE
  )
}

sea_pairwise_summary <- dplyr::bind_rows(
  sea_pairwise_rows
)

write.csv(
  sea_pairwise_summary,
  file.path(
    TABLE_DIR,
    "Table_4_pairwise_SEA_comparisons.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# =============================================================================
# 11. POSTERIOR CENTROID LOCATIONS
#
# In the object returned by siberMVN(), columns 1:4 contain the covariance
# matrix and columns 5:6 contain the back-transformed isotope means.
# =============================================================================

centroid_summary_rows <- vector(
  mode = "list",
  length = length(
    expected_codes
  )
)

for (code_index in seq_along(
  expected_codes
)) {

  current_code <- expected_codes[
    code_index
  ]

  current_posterior <- ellipses_posterior[[
    current_code
  ]]

  d13C_summary <- summarise_numeric(
    current_posterior[
      ,
      5
    ]
  )

  d15N_summary <- summarise_numeric(
    current_posterior[
      ,
      6
    ]
  )

  centroid_summary_rows[[code_index]] <- data.frame(
    SIBER_code = current_code,
    Group_label = group_labels[[current_code]],
    Posterior_mean_d13C = d13C_summary[["Mean"]],
    Posterior_median_d13C = d13C_summary[["Median"]],
    Lower_95_d13C = d13C_summary[["Lower_95"]],
    Upper_95_d13C = d13C_summary[["Upper_95"]],
    Posterior_mean_d15N = d15N_summary[["Mean"]],
    Posterior_median_d15N = d15N_summary[["Median"]],
    Lower_95_d15N = d15N_summary[["Lower_95"]],
    Upper_95_d15N = d15N_summary[["Upper_95"]],
    stringsAsFactors = FALSE
  )
}

centroid_summary <- dplyr::bind_rows(
  centroid_summary_rows
)

centroid_summary <- merge(
  centroid_summary,
  group_lookup[
    ,
    c(
      "SIBER_code",
      "Reach_code",
      "Reach_name",
      "Campaign",
      "Campaign_label",
      "n_fish",
      "mean_d13C",
      "mean_d15N"
    ),
    drop = FALSE
  ],
  by = "SIBER_code",
  all.x = TRUE,
  sort = FALSE
)

centroid_summary <- centroid_summary[
  match(
    expected_codes,
    centroid_summary$SIBER_code
  ),
  ,
  drop = FALSE
]

write.csv(
  centroid_summary,
  file.path(
    TABLE_DIR,
    "Table_5_posterior_centroid_locations.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# =============================================================================
# 12. PAIRWISE CENTROID DIFFERENCES AND DISTANCES
# =============================================================================

wrap_angle_radians <- function(
  angle
) {
  atan2(
    sin(
      angle
    ),
    cos(
      angle
    )
  )
}

summarise_angle_degrees <- function(
  angle_radians
) {

  circular_mean <- atan2(
    mean(
      sin(
        angle_radians
      )
    ),
    mean(
      cos(
        angle_radians
      )
    )
  )

  relative_angles <- wrap_angle_radians(
    angle_radians -
    circular_mean
  )

  lower_relative <- unname(
    stats::quantile(
      relative_angles,
      probs = LOWER_PROB
    )
  )

  upper_relative <- unname(
    stats::quantile(
      relative_angles,
      probs = UPPER_PROB
    )
  )

  resultant_length <- sqrt(
    mean(
      cos(
        angle_radians
      )
    ) ^ 2 +
    mean(
      sin(
        angle_radians
      )
    ) ^ 2
  )

  c(
    Circular_mean_degrees = circular_mean * 180 / pi,
    Lower_95_degrees = wrap_angle_radians(
      circular_mean +
      lower_relative
    ) * 180 / pi,
    Upper_95_degrees = wrap_angle_radians(
      circular_mean +
      upper_relative
    ) * 180 / pi,
    Resultant_length = resultant_length
  )
}

centroid_pairwise_rows <- vector(
  mode = "list",
  length = nrow(
    pairwise_comparisons
  )
)

for (pair_index in seq_len(
  nrow(
    pairwise_comparisons
  )
)) {

  code_1 <- pairwise_comparisons$Ellipse_1_code[
    pair_index
  ]

  code_2 <- pairwise_comparisons$Ellipse_2_code[
    pair_index
  ]

  posterior_1 <- ellipses_posterior[[
    code_1
  ]]

  posterior_2 <- ellipses_posterior[[
    code_2
  ]]

  n_draws <- min(
    nrow(
      posterior_1
    ),
    nrow(
      posterior_2
    )
  )

  delta_d13C <- (
    posterior_2[
      seq_len(
        n_draws
      ),
      5
    ] -
    posterior_1[
      seq_len(
        n_draws
      ),
      5
    ]
  )

  delta_d15N <- (
    posterior_2[
      seq_len(
        n_draws
      ),
      6
    ] -
    posterior_1[
      seq_len(
        n_draws
      ),
      6
    ]
  )

  distance_draws <- sqrt(
    delta_d13C ^ 2 +
    delta_d15N ^ 2
  )

  angle_draws <- atan2(
    delta_d15N,
    delta_d13C
  )

  delta_d13C_summary <- summarise_numeric(
    delta_d13C
  )

  delta_d15N_summary <- summarise_numeric(
    delta_d15N
  )

  distance_summary <- summarise_numeric(
    distance_draws
  )

  angle_summary <- summarise_angle_degrees(
    angle_draws
  )

  centroid_pairwise_rows[[pair_index]] <- data.frame(
    Pair_ID = pairwise_comparisons$Pair_ID[
      pair_index
    ],
    Ellipse_1_code = code_1,
    Ellipse_2_code = code_2,
    Comparison_label = pairwise_comparisons$Comparison_label[
      pair_index
    ],
    Comparison_type = pairwise_comparisons$Comparison_type[
      pair_index
    ],
    Manuscript_priority = pairwise_comparisons$Manuscript_priority[
      pair_index
    ],
    Median_delta_d13C_2_minus_1 = delta_d13C_summary[["Median"]],
    Lower_95_delta_d13C = delta_d13C_summary[["Lower_95"]],
    Upper_95_delta_d13C = delta_d13C_summary[["Upper_95"]],
    Probability_delta_d13C_positive = mean(
      delta_d13C > 0
    ),
    Median_delta_d15N_2_minus_1 = delta_d15N_summary[["Median"]],
    Lower_95_delta_d15N = delta_d15N_summary[["Lower_95"]],
    Upper_95_delta_d15N = delta_d15N_summary[["Upper_95"]],
    Probability_delta_d15N_positive = mean(
      delta_d15N > 0
    ),
    Median_centroid_distance = distance_summary[["Median"]],
    Lower_95_distance = distance_summary[["Lower_95"]],
    Upper_95_distance = distance_summary[["Upper_95"]],
    Mean_direction_degrees = angle_summary[["Circular_mean_degrees"]],
    Lower_95_direction_degrees = angle_summary[["Lower_95_degrees"]],
    Upper_95_direction_degrees = angle_summary[["Upper_95_degrees"]],
    Direction_resultant_length = angle_summary[["Resultant_length"]],
    stringsAsFactors = FALSE
  )
}

centroid_pairwise_summary <- dplyr::bind_rows(
  centroid_pairwise_rows
)

write.csv(
  centroid_pairwise_summary,
  file.path(
    TABLE_DIR,
    "Table_6_pairwise_centroid_differences.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# Wet-to-dry movement within each reach
centroid_movement_summary <- merge(
  centroid_movement_pairs[
    ,
    c(
      "Reach_code",
      "Reach_name",
      "Wet_code",
      "Dry_code"
    ),
    drop = FALSE
  ],
  centroid_pairwise_summary,
  by.x = c(
    "Wet_code",
    "Dry_code"
  ),
  by.y = c(
    "Ellipse_1_code",
    "Ellipse_2_code"
  ),
  all.x = TRUE,
  sort = FALSE
)

write.csv(
  centroid_movement_summary,
  file.path(
    TABLE_DIR,
    "Table_7_wet_to_dry_centroid_movement.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# =============================================================================
# 13. BAYESIAN AND MAXIMUM-LIKELIHOOD ELLIPSE OVERLAP
#
# Overlap is calculated for 95% prediction ellipses. This is different from
# SEA.B, which describes the approximately 40% standard ellipse.
# =============================================================================

overlap_draws_list <- vector(
  mode = "list",
  length = nrow(
    pairwise_comparisons
  )
)

overlap_summary_rows <- vector(
  mode = "list",
  length = nrow(
    pairwise_comparisons
  )
)

overlap_ml_rows <- vector(
  mode = "list",
  length = nrow(
    pairwise_comparisons
  )
)

cat(
  "\nBeginning Bayesian ellipse-overlap calculations.",
  "\nComparisons: ",
  nrow(
    pairwise_comparisons
  ),
  "\nRequested posterior draws per comparison: ",
  OVERLAP_DRAWS,
  "\n",
  sep = ""
)

for (pair_index in seq_len(
  nrow(
    pairwise_comparisons
  )
)) {

  code_1 <- pairwise_comparisons$Ellipse_1_code[
    pair_index
  ]

  code_2 <- pairwise_comparisons$Ellipse_2_code[
    pair_index
  ]

  draws_to_use <- min(
    OVERLAP_DRAWS,
    nrow(
      ellipses_posterior[[
        code_1
      ]]
    ),
    nrow(
      ellipses_posterior[[
        code_2
      ]]
    )
  )

  cat(
    "Overlap ",
    pair_index,
    " / ",
    nrow(
      pairwise_comparisons
    ),
    ": ",
    code_1,
    " vs ",
    code_2,
    "\n",
    sep = ""
  )

  current_overlap <- SIBER::bayesianOverlap(
    ellipse1 = code_1,
    ellipse2 = code_2,
    ellipses.posterior = ellipses_posterior,
    draws = draws_to_use,
    p.interval = OVERLAP_INTERVAL,
    n = ELLIPSE_POINTS,
    do.plot = FALSE
  )

  current_overlap <- as.data.frame(
    current_overlap,
    stringsAsFactors = FALSE
  )

  required_overlap_columns <- c(
    "area1",
    "area2",
    "overlap"
  )

  if (!all(
    required_overlap_columns %in% names(
      current_overlap
    )
  )) {
    stop(
      paste0(
        "Unexpected bayesianOverlap() output for ",
        code_1,
        " vs ",
        code_2,
        "."
      ),
      call. = FALSE
    )
  }

  current_overlap$Pair_ID <- pairwise_comparisons$Pair_ID[
    pair_index
  ]

  current_overlap$Ellipse_1_code <- code_1
  current_overlap$Ellipse_2_code <- code_2

  current_overlap$Overlap_percent_ellipse_1 <- (
    current_overlap$overlap /
    current_overlap$area1 *
    100
  )

  current_overlap$Overlap_percent_ellipse_2 <- (
    current_overlap$overlap /
    current_overlap$area2 *
    100
  )

  current_overlap$Overlap_percent_smaller_ellipse <- (
    current_overlap$overlap /
    pmin(
      current_overlap$area1,
      current_overlap$area2
    ) *
    100
  )

  current_overlap$Jaccard_overlap_percent <- (
    current_overlap$overlap /
    (
      current_overlap$area1 +
      current_overlap$area2 -
      current_overlap$overlap
    ) *
    100
  )

  overlap_draws_list[[pair_index]] <- current_overlap

  overlap_area_summary <- summarise_numeric(
    current_overlap$overlap
  )

  pct_1_summary <- summarise_numeric(
    current_overlap$Overlap_percent_ellipse_1
  )

  pct_2_summary <- summarise_numeric(
    current_overlap$Overlap_percent_ellipse_2
  )

  pct_smaller_summary <- summarise_numeric(
    current_overlap$Overlap_percent_smaller_ellipse
  )

  jaccard_summary <- summarise_numeric(
    current_overlap$Jaccard_overlap_percent
  )

  overlap_summary_rows[[pair_index]] <- data.frame(
    Pair_ID = pairwise_comparisons$Pair_ID[
      pair_index
    ],
    Ellipse_1_code = code_1,
    Ellipse_2_code = code_2,
    Comparison_label = pairwise_comparisons$Comparison_label[
      pair_index
    ],
    Comparison_type = pairwise_comparisons$Comparison_type[
      pair_index
    ],
    Manuscript_priority = pairwise_comparisons$Manuscript_priority[
      pair_index
    ],
    Posterior_draws_used = draws_to_use,
    Median_overlap_area = overlap_area_summary[["Median"]],
    Lower_95_overlap_area = overlap_area_summary[["Lower_95"]],
    Upper_95_overlap_area = overlap_area_summary[["Upper_95"]],
    Median_percent_ellipse_1 = pct_1_summary[["Median"]],
    Lower_95_percent_ellipse_1 = pct_1_summary[["Lower_95"]],
    Upper_95_percent_ellipse_1 = pct_1_summary[["Upper_95"]],
    Median_percent_ellipse_2 = pct_2_summary[["Median"]],
    Lower_95_percent_ellipse_2 = pct_2_summary[["Lower_95"]],
    Upper_95_percent_ellipse_2 = pct_2_summary[["Upper_95"]],
    Median_percent_smaller_ellipse = pct_smaller_summary[["Median"]],
    Lower_95_percent_smaller_ellipse = pct_smaller_summary[["Lower_95"]],
    Upper_95_percent_smaller_ellipse = pct_smaller_summary[["Upper_95"]],
    Median_Jaccard_percent = jaccard_summary[["Median"]],
    Lower_95_Jaccard_percent = jaccard_summary[["Lower_95"]],
    Upper_95_Jaccard_percent = jaccard_summary[["Upper_95"]],
    stringsAsFactors = FALSE
  )

  current_ml_overlap <- SIBER::maxLikOverlap(
    ellipse1 = code_1,
    ellipse2 = code_2,
    siber.object = siber_object,
    p.interval = OVERLAP_INTERVAL,
    n = ELLIPSE_POINTS,
    do.plot = FALSE
  )

  current_ml_overlap <- as.numeric(
    current_ml_overlap
  )

  if (length(
    current_ml_overlap
  ) != 3L) {
    stop(
      "Unexpected maxLikOverlap() output.",
      call. = FALSE
    )
  }

  ml_area_1 <- current_ml_overlap[1]
  ml_area_2 <- current_ml_overlap[2]
  ml_overlap_area <- current_ml_overlap[3]

  overlap_ml_rows[[pair_index]] <- data.frame(
    Pair_ID = pairwise_comparisons$Pair_ID[
      pair_index
    ],
    Ellipse_1_code = code_1,
    Ellipse_2_code = code_2,
    ML_area_ellipse_1 = ml_area_1,
    ML_area_ellipse_2 = ml_area_2,
    ML_overlap_area = ml_overlap_area,
    ML_percent_ellipse_1 = (
      ml_overlap_area /
      ml_area_1 *
      100
    ),
    ML_percent_ellipse_2 = (
      ml_overlap_area /
      ml_area_2 *
      100
    ),
    ML_percent_smaller_ellipse = (
      ml_overlap_area /
      min(
        ml_area_1,
        ml_area_2
      ) *
      100
    ),
    ML_Jaccard_percent = (
      ml_overlap_area /
      (
        ml_area_1 +
        ml_area_2 -
        ml_overlap_area
      ) *
      100
    ),
    stringsAsFactors = FALSE
  )
}

overlap_draws <- dplyr::bind_rows(
  overlap_draws_list
)

overlap_summary <- dplyr::bind_rows(
  overlap_summary_rows
)

overlap_ml_summary <- dplyr::bind_rows(
  overlap_ml_rows
)

saveRDS(
  overlap_draws,
  file.path(
    TABLE_DIR,
    "Bayesian_overlap_draws.rds"
  )
)

write.csv(
  overlap_summary,
  file.path(
    TABLE_DIR,
    "Table_8_Bayesian_95pct_ellipse_overlap.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  overlap_ml_summary,
  file.path(
    TABLE_DIR,
    "Table_S3_ML_95pct_ellipse_overlap.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# =============================================================================
# 14. CREATE A FULL PAIRWISE OVERLAP MATRIX
# =============================================================================

heatmap_group_order <- expected_codes

overlap_matrix <- matrix(
  NA_real_,
  nrow = length(
    heatmap_group_order
  ),
  ncol = length(
    heatmap_group_order
  ),
  dimnames = list(
    heatmap_group_order,
    heatmap_group_order
  )
)

diag(
  overlap_matrix
) <- 100

for (row_index in seq_len(
  nrow(
    overlap_summary
  )
)) {

  code_1 <- overlap_summary$Ellipse_1_code[
    row_index
  ]

  code_2 <- overlap_summary$Ellipse_2_code[
    row_index
  ]

  overlap_value <- overlap_summary$Median_Jaccard_percent[
    row_index
  ]

  overlap_matrix[
    code_1,
    code_2
  ] <- overlap_value

  overlap_matrix[
    code_2,
    code_1
  ] <- overlap_value
}

overlap_matrix_long <- expand.grid(
  Row_code = heatmap_group_order,
  Column_code = heatmap_group_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

overlap_matrix_long$Median_Jaccard_percent <- mapply(
  function(
    row_code,
    column_code
  ) {

    overlap_matrix[
      row_code,
      column_code
    ]
  },
  overlap_matrix_long$Row_code,
  overlap_matrix_long$Column_code
)

overlap_matrix_long$Row_label <- short_group_labels[
  overlap_matrix_long$Row_code
]

overlap_matrix_long$Column_label <- short_group_labels[
  overlap_matrix_long$Column_code
]

write.csv(
  overlap_matrix_long,
  file.path(
    TABLE_DIR,
    "Table_9_pairwise_overlap_matrix.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# =============================================================================
# 15. OPTIONAL COMBINED EXCEL RESULTS WORKBOOK
# =============================================================================

if (
  requireNamespace(
    "writexl",
    quietly = TRUE
  )
) {

  writexl::write_xlsx(
    list(
      ML_Niche_Metrics = ml_metric_table,
      MCMC_Diagnostics = diagnostic_group_table,
      Bayesian_SEA = sea_summary,
      Pairwise_SEA = sea_pairwise_summary,
      Centroid_Locations = centroid_summary,
      Pairwise_Centroids = centroid_pairwise_summary,
      Seasonal_Movement = centroid_movement_summary,
      Bayesian_Overlap = overlap_summary,
      ML_Overlap = overlap_ml_summary,
      Overlap_Matrix = overlap_matrix_long
    ),
    path = file.path(
      RESULTS_DIR,
      "SIBER_Bayesian_isotopic_niche_results.xlsx"
    )
  )
}


# =============================================================================
# 16. PUBLICATION FIGURES
# =============================================================================

COL_WET <- "#0072B2"
COL_DRY <- "#D55E00"
COL_SEAC <- "#222222"

theme_publication <- function(
  base_size = 14
) {

  ggplot2::theme_classic(
    base_size = base_size,
    base_family = "sans"
  ) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(
        size = AXIS_TEXT_SIZE,
        colour = "black"
      ),
      axis.title = ggplot2::element_text(
        size = AXIS_TITLE_SIZE
      ),
      strip.background = ggplot2::element_rect(
        fill = "grey94",
        colour = "black",
        linewidth = 0.4
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = FACET_TEXT_SIZE
      ),
      legend.position = "top",
      legend.title = ggplot2::element_text(
        face = "bold",
        size = LEGEND_TITLE_SIZE
      ),
      legend.text = ggplot2::element_text(
        size = LEGEND_TEXT_SIZE
      ),
      plot.margin = ggplot2::margin(
        8,
        10,
        8,
        8
      )
    ) + hr_typography()
}

save_publication_figure <- function(
  plot_object,
  filename,
  width,
  height
) {

  hr_ggsave(
    filename = file.path(
      FIGURE_DIR,
      paste0(
        filename,
        ".png"
      )
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  hr_ggsave(
    filename = file.path(
      FIGURE_DIR,
      paste0(
        filename,
        ".tiff"
      )
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )

  pdf_saved <- tryCatch(
    {
      hr_ggsave(
        filename = file.path(
          FIGURE_DIR,
          paste0(
            filename,
            ".pdf"
          )
        ),
        plot = plot_object,
        device = grDevices::cairo_pdf,
        width = width,
        height = height,
        units = "in",
        bg = "white"
      )
      TRUE
    },
    error = function(e) {
      FALSE
    }
  )

  if (!pdf_saved) {

    hr_ggsave(
      filename = file.path(
        FIGURE_DIR,
        paste0(
          filename,
          ".pdf"
        )
      ),
      plot = plot_object,
      device = "pdf",
      width = width,
      height = height,
      units = "in",
      bg = "white"
    )
  }
}


# -----------------------------------------------------------------------------
# Figure 1. Raw isotope data with ML-corrected standard ellipses
# -----------------------------------------------------------------------------

plot_data <- data.frame(
  d13C = siber_data$iso1,
  d15N = siber_data$iso2,
  group = siber_data$group,
  community = siber_data$community,
  stringsAsFactors = FALSE
)

plot_data$SIBER_code <- paste(
  plot_data$community,
  plot_data$group,
  sep = "."
)

plot_data <- merge(
  plot_data,
  group_lookup[
    ,
    c(
      "SIBER_code",
      "Reach_name",
      "Campaign_label"
    ),
    drop = FALSE
  ],
  by = "SIBER_code",
  all.x = TRUE,
  sort = FALSE
)

create_ml_ellipse_coordinates <- function(
  code
) {

  code_parts <- as.integer(
    strsplit(
      code,
      ".",
      fixed = TRUE
    )[[1]]
  )

  community_index <- code_parts[1]
  group_index <- code_parts[2]

  current_mu <- as.numeric(
    siber_object$ML.mu[[
      community_index
    ]][
      ,
      ,
      group_index
    ]
  )

  current_covariance <- siber_object$ML.cov[[
    community_index
  ]][
    ,
    ,
    group_index
  ]

  current_sample_size <- as.numeric(
    siber_object$sample.sizes[
      community_index,
      as.character(
        group_index
      )
    ]
  )

  current_coordinates <- SIBER::addEllipse(
    mu = current_mu,
    sigma = current_covariance,
    m = current_sample_size,
    n = ELLIPSE_POINTS,
    p.interval = NULL,
    ci.mean = FALSE,
    small.sample = TRUE,
    do.plot = FALSE
  )

  current_lookup <- group_lookup[
    group_lookup$SIBER_code == code,
    ,
    drop = FALSE
  ]

  data.frame(
    d13C = current_coordinates[
      ,
      1
    ],
    d15N = current_coordinates[
      ,
      2
    ],
    SIBER_code = code,
    Reach_name = current_lookup$Reach_name[1],
    Campaign_label = current_lookup$Campaign_label[1],
    stringsAsFactors = FALSE
  )
}

ellipse_coordinates <- dplyr::bind_rows(
  lapply(
    expected_codes,
    create_ml_ellipse_coordinates
  )
)

figure_1 <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = d13C,
    y = d15N,
    colour = Campaign_label,
    shape = Campaign_label
  )
) +
  ggplot2::geom_point(
    size = 1.8,
    alpha = 0.58
  ) +
  ggplot2::geom_path(
    data = ellipse_coordinates,
    mapping = ggplot2::aes(
      x = d13C,
      y = d15N,
      colour = Campaign_label,
      group = SIBER_code
    ),
    inherit.aes = FALSE,
    linewidth = 1.0
  ) +
  ggplot2::geom_point(
    data = centroid_summary,
    mapping = ggplot2::aes(
      x = Posterior_median_d13C,
      y = Posterior_median_d15N,
      fill = Campaign_label
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 3.2,
    stroke = 0.65,
    colour = "black"
  ) +
  ggplot2::facet_wrap(
    ~ Reach_name,
    nrow = 1
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "Wet campaign" = COL_WET,
      "Dry campaign" = COL_DRY
    )
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Wet campaign" = COL_WET,
      "Dry campaign" = COL_DRY
    ),
    guide = "none"
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "Wet campaign" = 16,
      "Dry campaign" = 17
    )
  ) +
  ggplot2::labs(
    x = expression(
      delta^13*C~("\u2030")
    ),
    y = expression(
      delta^15*N~("\u2030")
    ),
    colour = "Campaign",
    shape = "Campaign",
    caption = paste(
      "Lines show maximum-likelihood small-sample-corrected standard ellipses;",
      "large points show posterior median centroids."
    )
  ) +
  theme_publication(
    base_size = 14
  )

save_publication_figure(
  figure_1,
  "Figure_1_isotope_space_standard_ellipses",
  width = 7.6,
  height = 4.6
)


# -----------------------------------------------------------------------------
# Figure 2. Bayesian standard ellipse area
# -----------------------------------------------------------------------------

sea_plot_data <- sea_summary |>
  dplyr::mutate(
    Group_label_short = factor(
      short_group_labels[
        SIBER_code
      ],
      levels = short_group_labels[
        expected_codes
      ]
    )
  )

figure_2 <- ggplot2::ggplot(
  sea_plot_data,
  ggplot2::aes(
    x = Group_label_short,
    y = Posterior_median_SEA,
    colour = Campaign_label,
    shape = Campaign_label
  )
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = Lower_95,
      ymax = Upper_95
    ),
    width = 0,
    linewidth = 0.65
  ) +
  ggplot2::geom_point(
    size = 3
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      y = SEAc_ML
    ),
    shape = 4,
    size = 2.7,
    stroke = 0.9,
    colour = COL_SEAC
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "Wet campaign" = COL_WET,
      "Dry campaign" = COL_DRY
    )
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "Wet campaign" = 16,
      "Dry campaign" = 17
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = expression(
      "Standard ellipse area ("
      * "\u2030"^2 *
      ")"
    ),
    colour = "Campaign",
    shape = "Campaign",
    caption = paste(
      "Points and intervals are posterior medians and 95% credible intervals;",
      "black crosses show maximum-likelihood SEAc."
    )
  ) +
  theme_publication(
    base_size = 14
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1
    )
  )

save_publication_figure(
  figure_2,
  "Figure_2_Bayesian_standard_ellipse_area",
  width = 7.1,
  height = 4.8
)


# -----------------------------------------------------------------------------
# Figure 3. Median Bayesian 95% ellipse-overlap heatmap
# -----------------------------------------------------------------------------

overlap_heatmap_data <- overlap_matrix_long |>
  dplyr::mutate(
    Row_label = factor(
      Row_label,
      levels = rev(
        short_group_labels[
          expected_codes
        ]
      )
    ),
    Column_label = factor(
      Column_label,
      levels = short_group_labels[
        expected_codes
      ]
    )
  )

figure_3 <- ggplot2::ggplot(
  overlap_heatmap_data,
  ggplot2::aes(
    x = Column_label,
    y = Row_label,
    fill = Median_Jaccard_percent
  )
) +
  ggplot2::geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf(
        "%.1f",
        Median_Jaccard_percent
      )
    ),
    size = 4.25
  ) +
  ggplot2::scale_fill_gradient(
    low = "white",
    high = "#2166AC",
    limits = c(
      0,
      100
    ),
    name = "Jaccard\noverlap (%)"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    caption = "Cells show posterior median Jaccard overlap of 95% prediction ellipses."
  ) +
  ggplot2::coord_fixed() +
  theme_publication(
    base_size = 14
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 40,
      hjust = 1
    ),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank()
  )

save_publication_figure(
  figure_3,
  "Figure_3_pairwise_ellipse_overlap_heatmap",
  width = 6.8,
  height = 5.9
)


# -----------------------------------------------------------------------------
# Figure 4. Posterior median wet-to-dry centroid movement
# -----------------------------------------------------------------------------

movement_plot_data <- centroid_movement_summary |>
  dplyr::mutate(
    Wet_d13C = centroid_summary$Posterior_median_d13C[
      match(
        Wet_code,
        centroid_summary$SIBER_code
      )
    ],
    Wet_d15N = centroid_summary$Posterior_median_d15N[
      match(
        Wet_code,
        centroid_summary$SIBER_code
      )
    ],
    Dry_d13C = centroid_summary$Posterior_median_d13C[
      match(
        Dry_code,
        centroid_summary$SIBER_code
      )
    ],
    Dry_d15N = centroid_summary$Posterior_median_d15N[
      match(
        Dry_code,
        centroid_summary$SIBER_code
      )
    ]
  )

figure_4 <- ggplot2::ggplot(
  movement_plot_data
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = Wet_d13C,
      y = Wet_d15N,
      xend = Dry_d13C,
      yend = Dry_d15N,
      linetype = Reach_name
    ),
    linewidth = 0.9,
    arrow = grid::arrow(
      length = grid::unit(
        0.18,
        "cm"
      ),
      type = "closed"
    )
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = Wet_d13C,
      y = Wet_d15N
    ),
    shape = 21,
    size = 3.2,
    stroke = 0.6,
    fill = COL_WET,
    colour = "black"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = Dry_d13C,
      y = Dry_d15N
    ),
    shape = 21,
    size = 3.2,
    stroke = 0.6,
    fill = COL_DRY,
    colour = "black"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      x = Dry_d13C,
      y = Dry_d15N,
      label = Reach_name
    ),
    nudge_x = 0.12,
    nudge_y = 0.10,
    size = 4.25,
    hjust = 0
  ) +
  ggplot2::labs(
    x = expression(
      delta^13*C~("\u2030")
    ),
    y = expression(
      delta^15*N~("\u2030")
    ),
    linetype = "Reach",
    caption = paste(
      "Arrows connect posterior median wet centroids (blue)",
      "to dry centroids (orange)."
    )
  ) +
  theme_publication(
    base_size = 14
  )

save_publication_figure(
  figure_4,
  "Figure_4_wet_to_dry_centroid_movement",
  width = 6.5,
  height = 5.2
)


# -----------------------------------------------------------------------------
# Figure 5. Wet-to-dry centroid movement distances
# -----------------------------------------------------------------------------

movement_distance_data <- centroid_movement_summary |>
  dplyr::mutate(
    Reach_name = factor(
      Reach_name,
      levels = REACH_LABEL_LEVELS
    )
  )

figure_5 <- ggplot2::ggplot(
  movement_distance_data,
  ggplot2::aes(
    x = Reach_name,
    y = Median_centroid_distance
  )
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = Lower_95_distance,
      ymax = Upper_95_distance
    ),
    width = 0,
    linewidth = 0.65
  ) +
  ggplot2::geom_point(
    size = 3,
    shape = 21,
    fill = "white",
    colour = "black"
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Wet-to-dry centroid distance (‰)",
    caption = "Points and intervals show posterior medians and 95% credible intervals."
  ) +
  theme_publication(
    base_size = 14
  )

save_publication_figure(
  figure_5,
  "Figure_5_centroid_movement_distance",
  width = 6.2,
  height = 4.3
)


# -----------------------------------------------------------------------------
# Figure S1. Pairwise probability that ellipse 2 has a larger SEA than ellipse 1
# -----------------------------------------------------------------------------

sea_probability_plot_data <- sea_pairwise_summary |>
  dplyr::filter(
    Manuscript_priority == "Primary"
  ) |>
  dplyr::mutate(
    Comparison_label = factor(
      Comparison_label,
      levels = rev(
        unique(
          Comparison_label
        )
      )
    )
  )

figure_s1 <- ggplot2::ggplot(
  sea_probability_plot_data,
  ggplot2::aes(
    x = Probability_SEA_2_greater_1,
    y = Comparison_label
  )
) +
  ggplot2::geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey40"
  ) +
  ggplot2::geom_point(
    size = 2.8
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 0.25
    )
  ) +
  ggplot2::labs(
    x = "Posterior probability that ellipse 2 has larger SEA",
    y = NULL
  ) +
  theme_publication(
    base_size = 14
  )

save_publication_figure(
  figure_s1,
  "Figure_S1_pairwise_SEA_probabilities",
  width = 7.4,
  height = 5.4
)


# =============================================================================
# 17. SAVE SESSION INFORMATION AND COMPLETION LOG
# =============================================================================

model_settings <- data.frame(
  Setting = c(
    "SIBER version",
    "JAGS version",
    "n.iter",
    "n.burnin setting",
    "n.thin",
    "n.chains",
    "Posterior draws per ellipse",
    "Overlap posterior draws",
    "Overlap prediction interval",
    "Ellipse boundary points"
  ),
  Value = c(
    as.character(
      utils::packageVersion("SIBER")
    ),
    jags_version,
    N_ITER,
    N_BURNIN,
    N_THIN,
    N_CHAINS,
    posterior_draw_count,
    OVERLAP_DRAWS,
    OVERLAP_INTERVAL,
    ELLIPSE_POINTS
  ),
  stringsAsFactors = FALSE
)

write.csv(
  model_settings,
  file.path(
    TABLE_DIR,
    "Table_S4_model_settings.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    LOG_DIR,
    "R_session_information.txt"
  )
)

unlink(
  RUNNING_FILE
)

writeLines(
  c(
    "SIBER Bayesian isotopic niche analysis completed successfully.",
    paste("Completed:", Sys.time()),
    paste("Input workbook:", INPUT_FILE),
    paste("Model file:", MODEL_FILE),
    paste("Results folder:", RESULTS_DIR),
    paste("Posterior draws per group:", posterior_draw_count),
    paste("Overlap draws per comparison:", min(OVERLAP_DRAWS, posterior_draw_count))
  ),
  COMPLETED_FILE
)

cat(
  "\n============================================================\n",
  "SIBER BAYESIAN ISOTOPIC NICHE ANALYSIS COMPLETED\n",
  "============================================================\n",
  "\nResults folder:\n",
  RESULTS_DIR,
  "\n\nMain tables:\n",
  "Table_1_maximum_likelihood_niche_metrics.csv\n",
  "Table_2_group_MCMC_diagnostics.csv\n",
  "Table_3_Bayesian_SEA_summary.csv\n",
  "Table_4_pairwise_SEA_comparisons.csv\n",
  "Table_5_posterior_centroid_locations.csv\n",
  "Table_6_pairwise_centroid_differences.csv\n",
  "Table_7_wet_to_dry_centroid_movement.csv\n",
  "Table_8_Bayesian_95pct_ellipse_overlap.csv\n",
  "Table_9_pairwise_overlap_matrix.csv\n",
  "\nMain figures:\n",
  "Figure_1_isotope_space_standard_ellipses\n",
  "Figure_2_Bayesian_standard_ellipse_area\n",
  "Figure_3_pairwise_ellipse_overlap_heatmap\n",
  "Figure_4_wet_to_dry_centroid_movement\n",
  "Figure_5_centroid_movement_distance\n",
  "Figure_S1_pairwise_SEA_probabilities\n",
  "============================================================\n",
  sep = ""
)

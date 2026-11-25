# RUN IN RSTUDIO: open this R file and click Source.
# Input file(s): Han_River_Within_Species_Turnover_Decomposition_R_Dataset(1).xlsx
# Worksheet tabs read automatically: Individual_Data
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
# WITHIN-SPECIES VERSUS SPECIES-TURNOVER DECOMPOSITION
# FINAL READY-TO-USE R SCRIPT
#
# PURPOSE
#   Determines whether wet-to-dry changes in community-weighted mean d13C and
#   d15N were caused by:
#     1. isotope shifts within species;
#     2. changes in the relative representation of species shared by both
#        campaigns; and
#     3. gains or losses of species between campaigns.
#
# METHOD
#   Exact midpoint (Kitagawa/Oaxaca-type) decomposition:
#
#     Community mean: C_t = sum_i p_it * z_it
#
#     Dry - wet change:
#       Delta C = sum_i pbar_i * Delta z_i
#               + sum_i Delta p_i * zbar_i
#
#   The first term is the within-species effect. The composition term is split
#   into shared-species reweighting and wet-only/dry-only turnover.
#
# INPUT WORKBOOK
#   Han_River_Within_Species_Turnover_Decomposition_R_Dataset(1).xlsx
#   Worksheet used by R: Individual_Data
#
# OUTPUT
#   Tables, bootstrap uncertainty summaries, publication figures, an optional
#   Excel results workbook, and validation logs.
#
# IMPORTANT INTERPRETATION
#   "sample_count" weights describe the composition of the analysed fish
#   sample. They represent true abundance only when sampling effort and
#   catchability are sufficiently comparable. "equal_species" is included as
#   a sensitivity analysis.
#
# HOW TO RUN
#   1. Save this script.
#   2. Put the input workbook in the same folder, Downloads, Desktop, or
#      Documents. Otherwise select it in the file-selection window.
#   3. Click Source once.
# =============================================================================


# =============================================================================
# 1. USER SETTINGS
# =============================================================================

N_BOOT <- 2000L
# Use 2000 for the final analysis.
# For a quick test only, temporarily set N_BOOT <- 100L.

SEED <- 20260730L

# Publication site labels. HC, QJ, and SY remain the internal data codes;
# every figure displays only the corresponding numbered site.
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

WEIGHTING_MODES <- c(
  "sample_count",
  "equal_species"
)

set.seed(SEED)
graphics.off()


# =============================================================================
# 2. CHECK AND LOAD PACKAGES
# =============================================================================

required_packages <- c(
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
  library(readxl)
  library(dplyr)
  library(ggplot2)
})


# =============================================================================
# 3. LOCATE THE INPUT WORKBOOK
# =============================================================================

INPUT_FILE <- hr_input("decomposition")
hr_validate_workbook(INPUT_FILE, "decomposition")


cat(
  "\nInput workbook:\n",
  INPUT_FILE,
  "\n",
  sep = ""
)


# =============================================================================
# 4. CREATE OUTPUT FOLDERS AND STATUS FILES
# =============================================================================

RESULTS_DIR <- hr_results("decomposition")

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
    "Within-species versus turnover decomposition is running.",
    paste("Started:", Sys.time()),
    paste("Input workbook:", INPUT_FILE),
    paste("Bootstrap replicates:", N_BOOT)
  ),
  RUNNING_FILE
)


# =============================================================================
# 5. READ AND VALIDATE INDIVIDUAL-LEVEL DATA
# =============================================================================

available_sheets <- readxl::excel_sheets(
  INPUT_FILE
)

if (!("Individual_Data" %in% available_sheets)) {
  stop(
    paste0(
      'Worksheet "Individual_Data" was not found.\n',
      "Available worksheets: ",
      paste(available_sheets, collapse = ", ")
    ),
    call. = FALSE
  )
}

fish_data <- hr_read_excel(
  INPUT_FILE,
  sheet = "Individual_Data"
)

fish_data <- as.data.frame(
  fish_data,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "Fish_id",
  "Reach_code",
  "Campaign",
  "Species_code",
  "Feeding_guild",
  "d13C_permil",
  "d15N_permil"
)

missing_columns <- setdiff(
  required_columns,
  names(fish_data)
)

if (length(missing_columns) > 0L) {
  stop(
    paste0(
      "The Individual_Data worksheet is missing column(s): ",
      paste(missing_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

fish_data$d13C_permil <- suppressWarnings(
  as.numeric(
    fish_data$d13C_permil
  )
)

fish_data$d15N_permil <- suppressWarnings(
  as.numeric(
    fish_data$d15N_permil
  )
)

if (anyNA(fish_data[required_columns])) {
  stop(
    "Missing or non-numeric values were found in required variables.",
    call. = FALSE
  )
}

expected_reaches <- c(
  "HC",
  "QJ",
  "SY"
)

expected_campaigns <- c(
  "Wet_2023",
  "Dry_2023"
)

unexpected_reaches <- setdiff(
  unique(
    fish_data$Reach_code
  ),
  expected_reaches
)

unexpected_campaigns <- setdiff(
  unique(
    fish_data$Campaign
  ),
  expected_campaigns
)

if (length(unexpected_reaches) > 0L) {
  stop(
    paste0(
      "Unexpected Reach_code value(s): ",
      paste(unexpected_reaches, collapse = ", ")
    ),
    call. = FALSE
  )
}

if (length(unexpected_campaigns) > 0L) {
  stop(
    paste0(
      "Unexpected Campaign value(s): ",
      paste(unexpected_campaigns, collapse = ", ")
    ),
    call. = FALSE
  )
}

fish_data <- fish_data |>
  dplyr::mutate(
    Reach_code = factor(
      as.character(Reach_code),
      levels = expected_reaches
    ),

    Campaign = factor(
      as.character(Campaign),
      levels = expected_campaigns
    ),

    Species_code = factor(
      as.character(Species_code)
    ),

    Feeding_guild = factor(
      as.character(Feeding_guild)
    )
  )

group_counts <- fish_data |>
  dplyr::count(
    Reach_code,
    Campaign,
    name = "n"
  )

if (nrow(group_counts) != 6L) {
  stop(
    "All three reaches must be represented in both campaigns.",
    call. = FALSE
  )
}

data_summary <- data.frame(
  statistic = c(
    "Individual observations",
    "Unique species",
    "Reach-campaign groups",
    "Wet observations",
    "Dry observations"
  ),

  value = c(
    nrow(fish_data),
    nlevels(fish_data$Species_code),
    nrow(group_counts),
    sum(fish_data$Campaign == "Wet_2023"),
    sum(fish_data$Campaign == "Dry_2023")
  )
)

write.csv(
  data_summary,
  file.path(
    TABLE_DIR,
    "Table_S1_data_summary.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat(
  "\nData validation completed.",
  "\nObservations: ", nrow(fish_data),
  "\nUnique species: ", nlevels(fish_data$Species_code),
  "\n",
  sep = ""
)


# =============================================================================
# 6. ANALYSIS FUNCTIONS
# =============================================================================

make_species_summary <- function(
  individual_data
) {

  species_summary <- individual_data |>
    dplyr::group_by(
      Reach_code,
      Campaign,
      Species_code,
      Feeding_guild
    ) |>
    dplyr::summarise(
      n_individuals = dplyr::n(),

      mean_d13C = mean(
        d13C_permil,
        na.rm = TRUE
      ),

      mean_d15N = mean(
        d15N_permil,
        na.rm = TRUE
      ),

      .groups = "drop"
    ) |>
    dplyr::group_by(
      Reach_code,
      Campaign
    ) |>
    dplyr::mutate(
      total_individuals = sum(
        n_individuals
      ),

      species_richness = dplyr::n(),

      p_sample_count = (
        n_individuals /
        total_individuals
      ),

      p_equal_species = (
        1 /
        species_richness
      )
    ) |>
    dplyr::ungroup()

  species_summary
}


decompose_one_reach <- function(
  species_summary,
  reach_code,
  response_column,
  response_label,
  weighting
) {

  weight_column <- if (
    weighting == "sample_count"
  ) {
    "p_sample_count"
  } else {
    "p_equal_species"
  }

  wet_data <- species_summary[
    species_summary$Reach_code == reach_code &
    species_summary$Campaign == "Wet_2023",
    c(
      "Species_code",
      "Feeding_guild",
      "n_individuals",
      weight_column,
      response_column
    ),
    drop = FALSE
  ]

  dry_data <- species_summary[
    species_summary$Reach_code == reach_code &
    species_summary$Campaign == "Dry_2023",
    c(
      "Species_code",
      "Feeding_guild",
      "n_individuals",
      weight_column,
      response_column
    ),
    drop = FALSE
  ]

  names(wet_data) <- c(
    "Species_code",
    "Feeding_guild_wet",
    "n_wet",
    "p_wet",
    "z_wet"
  )

  names(dry_data) <- c(
    "Species_code",
    "Feeding_guild_dry",
    "n_dry",
    "p_dry",
    "z_dry"
  )

  merged_data <- merge(
    wet_data,
    dry_data,
    by = "Species_code",
    all = TRUE,
    sort = TRUE
  )

  merged_data$Feeding_guild <- ifelse(
    !is.na(
      merged_data$Feeding_guild_wet
    ),
    as.character(
      merged_data$Feeding_guild_wet
    ),
    as.character(
      merged_data$Feeding_guild_dry
    )
  )

  merged_data$n_wet[
    is.na(
      merged_data$n_wet
    )
  ] <- 0

  merged_data$n_dry[
    is.na(
      merged_data$n_dry
    )
  ] <- 0

  merged_data$p_wet[
    is.na(
      merged_data$p_wet
    )
  ] <- 0

  merged_data$p_dry[
    is.na(
      merged_data$p_dry
    )
  ] <- 0

  merged_data$Status <- ifelse(
    merged_data$n_wet > 0 &
    merged_data$n_dry > 0,
    "Shared",
    ifelse(
      merged_data$n_wet > 0,
      "Wet_only",
      "Dry_only"
    )
  )

  wet_mean_for_sum <- merged_data$z_wet
  wet_mean_for_sum[
    is.na(
      wet_mean_for_sum
    )
  ] <- 0

  dry_mean_for_sum <- merged_data$z_dry
  dry_mean_for_sum[
    is.na(
      dry_mean_for_sum
    )
  ] <- 0

  wet_community_mean <- sum(
    merged_data$p_wet *
    wet_mean_for_sum
  )

  dry_community_mean <- sum(
    merged_data$p_dry *
    dry_mean_for_sum
  )

  reference_mean <- mean(
    c(
      wet_community_mean,
      dry_community_mean
    )
  )

  merged_data$p_bar <- (
    merged_data$p_wet +
    merged_data$p_dry
  ) /
    2

  merged_data$delta_p <- (
    merged_data$p_dry -
    merged_data$p_wet
  )

  merged_data$z_mid <- ifelse(
    merged_data$Status == "Shared",
    (
      merged_data$z_wet +
      merged_data$z_dry
    ) /
      2,
    ifelse(
      merged_data$Status == "Wet_only",
      merged_data$z_wet,
      merged_data$z_dry
    )
  )

  merged_data$delta_z <- ifelse(
    merged_data$Status == "Shared",
    (
      merged_data$z_dry -
      merged_data$z_wet
    ),
    0
  )

  merged_data$Within_species <- (
    merged_data$p_bar *
    merged_data$delta_z
  )

  merged_data$Shared_species_reweighting <- ifelse(
    merged_data$Status == "Shared",
    (
      merged_data$delta_p *
      (
        merged_data$z_mid -
        reference_mean
      )
    ),
    0
  )

  merged_data$Gains_losses_turnover <- ifelse(
    merged_data$Status != "Shared",
    (
      merged_data$delta_p *
      (
        merged_data$z_mid -
        reference_mean
      )
    ),
    0
  )

  merged_data$Composition_total <- (
    merged_data$Shared_species_reweighting +
    merged_data$Gains_losses_turnover
  )

  merged_data$Total_contribution <- (
    merged_data$Within_species +
    merged_data$Composition_total
  )

  total_change <- (
    dry_community_mean -
    wet_community_mean
  )

  within_total <- sum(
    merged_data$Within_species
  )

  shared_reweighting_total <- sum(
    merged_data$Shared_species_reweighting
  )

  turnover_total <- sum(
    merged_data$Gains_losses_turnover
  )

  composition_total <- (
    shared_reweighting_total +
    turnover_total
  )

  reconstructed_change <- (
    within_total +
    composition_total
  )

  closure_error <- (
    total_change -
    reconstructed_change
  )

  absolute_denominator <- sum(
    abs(
      c(
        within_total,
        shared_reweighting_total,
        turnover_total
      )
    )
  )

  if (
    !is.finite(
      absolute_denominator
    ) ||
    absolute_denominator == 0
  ) {

    within_share <- 0
    reweighting_share <- 0
    turnover_share <- 0

  } else {

    within_share <- (
      abs(
        within_total
      ) /
      absolute_denominator
    )

    reweighting_share <- (
      abs(
        shared_reweighting_total
      ) /
      absolute_denominator
    )

    turnover_share <- (
      abs(
        turnover_total
      ) /
      absolute_denominator
    )
  }

  dominant_values <- abs(
    c(
      "Within-species" = within_total,
      "Shared-species reweighting" = shared_reweighting_total,
      "Gains/losses turnover" = turnover_total
    )
  )

  dominant_process <- names(
    dominant_values
  )[
    which.max(
      dominant_values
    )
  ]

  summary_row <- data.frame(
    Weighting = weighting,
    Reach_code = as.character(
      reach_code
    ),
    Response = response_label,
    Wet_community_mean = wet_community_mean,
    Dry_community_mean = dry_community_mean,
    Total_change = total_change,
    Within_species = within_total,
    Shared_species_reweighting = shared_reweighting_total,
    Gains_losses_turnover = turnover_total,
    Composition_total = composition_total,
    Reconstructed_change = reconstructed_change,
    Closure_error = closure_error,
    Dominant_process = dominant_process,
    Absolute_within_share = within_share,
    Absolute_reweighting_share = reweighting_share,
    Absolute_turnover_share = turnover_share,
    stringsAsFactors = FALSE
  )

  species_rows <- data.frame(
    Weighting = weighting,
    Reach_code = as.character(
      reach_code
    ),
    Response = response_label,
    Species_code = as.character(
      merged_data$Species_code
    ),
    Feeding_guild = merged_data$Feeding_guild,
    Status = merged_data$Status,
    n_wet = merged_data$n_wet,
    n_dry = merged_data$n_dry,
    p_wet = merged_data$p_wet,
    p_dry = merged_data$p_dry,
    z_wet = merged_data$z_wet,
    z_dry = merged_data$z_dry,
    p_bar = merged_data$p_bar,
    delta_p = merged_data$delta_p,
    z_mid = merged_data$z_mid,
    delta_z = merged_data$delta_z,
    Reference_mean = reference_mean,
    Within_species = merged_data$Within_species,
    Shared_species_reweighting = merged_data$Shared_species_reweighting,
    Gains_losses_turnover = merged_data$Gains_losses_turnover,
    Composition_total = merged_data$Composition_total,
    Total_contribution = merged_data$Total_contribution,
    stringsAsFactors = FALSE
  )

  list(
    summary = summary_row,
    species = species_rows
  )
}


run_all_decompositions <- function(
  species_summary,
  include_species = TRUE
) {

  summary_results <- list()
  species_results <- list()

  summary_index <- 1L
  species_index <- 1L

  for (weighting in WEIGHTING_MODES) {

    for (reach_code in expected_reaches) {

      for (response_label in c(
        "d13C",
        "d15N"
      )) {

        response_column <- if (
          response_label == "d13C"
        ) {
          "mean_d13C"
        } else {
          "mean_d15N"
        }

        current_result <- decompose_one_reach(
          species_summary = species_summary,
          reach_code = reach_code,
          response_column = response_column,
          response_label = response_label,
          weighting = weighting
        )

        summary_results[[summary_index]] <- current_result$summary
        summary_index <- summary_index + 1L

        if (isTRUE(include_species)) {

          species_results[[species_index]] <- current_result$species
          species_index <- species_index + 1L
        }
      }
    }
  }

  combined_summary <- dplyr::bind_rows(
    summary_results
  )

  if (isTRUE(include_species)) {

    combined_species <- dplyr::bind_rows(
      species_results
    )

  } else {

    combined_species <- NULL
  }

  list(
    summary = combined_summary,
    species = combined_species
  )
}


# =============================================================================
# 7. OBSERVED DECOMPOSITION
# =============================================================================

observed_species_summary <- make_species_summary(
  fish_data
)

observed_results <- run_all_decompositions(
  species_summary = observed_species_summary,
  include_species = TRUE
)

observed_decomposition <- observed_results$summary
species_contributions <- observed_results$species

maximum_observed_closure_error <- max(
  abs(
    observed_decomposition$Closure_error
  )
)

if (
  !is.finite(
    maximum_observed_closure_error
  ) ||
  maximum_observed_closure_error > 1e-8
) {
  stop(
    paste0(
      "Observed decomposition did not close exactly. ",
      "Maximum absolute error = ",
      maximum_observed_closure_error
    ),
    call. = FALSE
  )
}

write.csv(
  observed_decomposition,
  file.path(
    TABLE_DIR,
    "Table_1_observed_reach_decomposition.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  species_contributions,
  file.path(
    TABLE_DIR,
    "Table_4_species_level_contributions.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  observed_species_summary,
  file.path(
    TABLE_DIR,
    "Table_S2_species_campaign_summary.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat(
  "\nObserved decomposition completed.",
  "\nMaximum closure error: ",
  format(
    maximum_observed_closure_error,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)


# =============================================================================
# 8. NONPARAMETRIC BOOTSTRAP
#
# Individuals are resampled with replacement within each reach x campaign
# group, preserving the observed sample size of every group.
# =============================================================================

group_identifier <- interaction(
  fish_data$Reach_code,
  fish_data$Campaign,
  drop = TRUE,
  lex.order = TRUE
)

group_indices <- split(
  seq_len(
    nrow(fish_data)
  ),
  group_identifier
)

bootstrap_results <- vector(
  mode = "list",
  length = N_BOOT
)

cat(
  "\nBeginning ",
  N_BOOT,
  " bootstrap replicates.\n",
  sep = ""
)

for (bootstrap_index in seq_len(N_BOOT)) {

  sampled_indices <- unlist(
    lapply(
      group_indices,
      function(index_vector) {

        sample(
          index_vector,
          size = length(
            index_vector
          ),
          replace = TRUE
        )
      }
    ),
    use.names = FALSE
  )

  bootstrap_data <- fish_data[
    sampled_indices,
    ,
    drop = FALSE
  ]

  bootstrap_species_summary <- make_species_summary(
    bootstrap_data
  )

  bootstrap_decomposition <- run_all_decompositions(
    species_summary = bootstrap_species_summary,
    include_species = FALSE
  )$summary

  bootstrap_decomposition$Bootstrap_replicate <- bootstrap_index

  bootstrap_results[[bootstrap_index]] <- bootstrap_decomposition

  if (
    bootstrap_index == 1L ||
    bootstrap_index %% 100L == 0L ||
    bootstrap_index == N_BOOT
  ) {
    cat(
      "Bootstrap replicate ",
      bootstrap_index,
      " / ",
      N_BOOT,
      " completed.\n",
      sep = ""
    )
  }
}

bootstrap_draws <- dplyr::bind_rows(
  bootstrap_results
)

maximum_bootstrap_closure_error <- max(
  abs(
    bootstrap_draws$Closure_error
  )
)

if (
  !is.finite(
    maximum_bootstrap_closure_error
  ) ||
  maximum_bootstrap_closure_error > 1e-8
) {
  stop(
    paste0(
      "At least one bootstrap decomposition failed the closure check. ",
      "Maximum absolute error = ",
      maximum_bootstrap_closure_error
    ),
    call. = FALSE
  )
}

saveRDS(
  bootstrap_draws,
  file.path(
    TABLE_DIR,
    "bootstrap_decomposition_draws.rds"
  )
)


# =============================================================================
# 9. BOOTSTRAP UNCERTAINTY SUMMARIES
# =============================================================================

component_columns <- c(
  "Total_change",
  "Within_species",
  "Shared_species_reweighting",
  "Gains_losses_turnover",
  "Composition_total"
)

bootstrap_long_list <- vector(
  mode = "list",
  length = length(
    component_columns
  )
)

for (
  component_index in seq_along(
    component_columns
  )
) {

  component_name <- component_columns[
    component_index
  ]

  bootstrap_long_list[[component_index]] <- data.frame(
    Weighting = bootstrap_draws$Weighting,
    Reach_code = bootstrap_draws$Reach_code,
    Response = bootstrap_draws$Response,
    Bootstrap_replicate = bootstrap_draws$Bootstrap_replicate,
    Component = component_name,
    Value = bootstrap_draws[[component_name]],
    stringsAsFactors = FALSE
  )
}

bootstrap_long <- dplyr::bind_rows(
  bootstrap_long_list
)

bootstrap_component_intervals <- bootstrap_long |>
  dplyr::group_by(
    Weighting,
    Reach_code,
    Response,
    Component
  ) |>
  dplyr::summarise(
    Bootstrap_mean = mean(
      Value,
      na.rm = TRUE
    ),

    Bootstrap_median = stats::median(
      Value,
      na.rm = TRUE
    ),

    Lower_95 = unname(
      stats::quantile(
        Value,
        probs = 0.025,
        na.rm = TRUE
      )
    ),

    Upper_95 = unname(
      stats::quantile(
        Value,
        probs = 0.975,
        na.rm = TRUE
      )
    ),

    Probability_positive = mean(
      Value > 0,
      na.rm = TRUE
    ),

    Probability_negative = mean(
      Value < 0,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

observed_long_list <- vector(
  mode = "list",
  length = length(
    component_columns
  )
)

for (
  component_index in seq_along(
    component_columns
  )
) {

  component_name <- component_columns[
    component_index
  ]

  observed_long_list[[component_index]] <- data.frame(
    Weighting = observed_decomposition$Weighting,
    Reach_code = observed_decomposition$Reach_code,
    Response = observed_decomposition$Response,
    Component = component_name,
    Observed_value = observed_decomposition[[component_name]],
    stringsAsFactors = FALSE
  )
}

observed_long <- dplyr::bind_rows(
  observed_long_list
)

bootstrap_component_intervals <- merge(
  bootstrap_component_intervals,
  observed_long,
  by = c(
    "Weighting",
    "Reach_code",
    "Response",
    "Component"
  ),
  all.x = TRUE,
  sort = FALSE
)

bootstrap_component_intervals$Credible_direction_95 <- (
  bootstrap_component_intervals$Lower_95 > 0 |
  bootstrap_component_intervals$Upper_95 < 0
)

write.csv(
  bootstrap_component_intervals,
  file.path(
    TABLE_DIR,
    "Table_2_bootstrap_component_intervals.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# Dominant process probability
dominant_matrix <- cbind(
  abs(
    bootstrap_draws$Within_species
  ),
  abs(
    bootstrap_draws$Shared_species_reweighting
  ),
  abs(
    bootstrap_draws$Gains_losses_turnover
  )
)

dominant_names <- c(
  "Within-species",
  "Shared-species reweighting",
  "Gains/losses turnover"
)

bootstrap_draws$Dominant_process <- dominant_names[
  max.col(
    dominant_matrix,
    ties.method = "first"
  )
]

dominant_process_probabilities <- bootstrap_draws |>
  dplyr::count(
    Weighting,
    Reach_code,
    Response,
    Dominant_process,
    name = "Bootstrap_count"
  ) |>
  dplyr::group_by(
    Weighting,
    Reach_code,
    Response
  ) |>
  dplyr::mutate(
    Probability_dominant = (
      Bootstrap_count /
      sum(
        Bootstrap_count
      )
    )
  ) |>
  dplyr::ungroup()

write.csv(
  dominant_process_probabilities,
  file.path(
    TABLE_DIR,
    "Table_3_dominant_process_probabilities.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# Validation table
validation_summary <- data.frame(
  Check = c(
    "Observed maximum absolute closure error",
    "Bootstrap maximum absolute closure error",
    "Bootstrap replicates",
    "Weighting modes",
    "Reach-response combinations per weighting"
  ),

  Value = c(
    maximum_observed_closure_error,
    maximum_bootstrap_closure_error,
    N_BOOT,
    length(
      WEIGHTING_MODES
    ),
    length(
      expected_reaches
    ) * 2
  ),

  Criterion = c(
    "< 1e-8",
    "< 1e-8",
    "User setting",
    "2",
    "6"
  ),

  stringsAsFactors = FALSE
)

write.csv(
  validation_summary,
  file.path(
    TABLE_DIR,
    "Table_5_decomposition_validation.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# =============================================================================
# 10. OPTIONAL EXCEL RESULTS WORKBOOK
# =============================================================================

if (
  requireNamespace(
    "writexl",
    quietly = TRUE
  )
) {

  writexl::write_xlsx(
    list(
      Observed_Decomposition = observed_decomposition,
      Bootstrap_Intervals = bootstrap_component_intervals,
      Dominant_Probabilities = dominant_process_probabilities,
      Species_Contributions = species_contributions,
      Validation = validation_summary
    ),
    path = file.path(
      RESULTS_DIR,
      "Within_species_turnover_decomposition_results.xlsx"
    )
  )
}


# =============================================================================
# 11. PUBLICATION FIGURES
# =============================================================================

reach_labels <- REACH_LABELS

component_labels <- c(
  "Within_species" = "Within-species",
  "Shared_species_reweighting" = "Shared-species reweighting",
  "Gains_losses_turnover" = "Species gains/losses",
  "Composition_total" = "Composition total",
  "Total_change" = "Total change"
)

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

  pdf_success <- tryCatch(
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

  if (!pdf_success) {

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


# Figure 1: detailed decomposition, sample-count weighting
figure_1_data <- bootstrap_component_intervals |>
  dplyr::filter(
    Weighting == "sample_count",
    Component %in% c(
      "Within_species",
      "Shared_species_reweighting",
      "Gains_losses_turnover"
    )
  ) |>
  dplyr::mutate(
    Reach = factor(
      reach_labels[
        Reach_code
      ],
      levels = REACH_LABEL_LEVELS
    ),

    Response_label = dplyr::recode(
      Response,
      d13C = "d13C",
      d15N = "d15N"
    ),

    Component_label = factor(
      component_labels[
        Component
      ],
      levels = c(
        "Within-species",
        "Shared-species reweighting",
        "Species gains/losses"
      )
    )
  )

figure_1 <- ggplot2::ggplot(
  figure_1_data,
  ggplot2::aes(
    x = Component_label,
    y = Observed_value
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey40"
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
    ggplot2::aes(
      fill = Component_label
    ),
    shape = 21,
    size = 3,
    stroke = 0.55,
    colour = "black"
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(
      Response_label
    ),
    cols = ggplot2::vars(
      Reach
    ),
    scales = "free_y"
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Within-species" = "#0072B2",
      "Shared-species reweighting" = "#E69F00",
      "Species gains/losses" = "#009E73"
    ),
    guide = "none"
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Contribution to dry - wet change (per mil)",
    caption = paste(
      "Points show the observed exact decomposition;",
      "intervals show 95% individual-level bootstrap intervals."
    )
  ) +
  theme_publication(
    base_size = 14
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 30,
      hjust = 1
    )
  )

save_publication_figure(
  figure_1,
  "Figure_1_count_weighted_decomposition_components",
  width = 8.2,
  height = 5.7
)


# Figure 2: within-species versus total composition
figure_2_data <- bootstrap_component_intervals |>
  dplyr::filter(
    Weighting == "sample_count",
    Component %in% c(
      "Within_species",
      "Composition_total",
      "Total_change"
    )
  ) |>
  dplyr::mutate(
    Reach = factor(
      reach_labels[
        Reach_code
      ],
      levels = REACH_LABEL_LEVELS
    ),

    Component_label = factor(
      component_labels[
        Component
      ],
      levels = c(
        "Within-species",
        "Composition total",
        "Total change"
      )
    )
  )

figure_2 <- ggplot2::ggplot(
  figure_2_data,
  ggplot2::aes(
    x = Reach,
    y = Observed_value,
    colour = Component_label,
    shape = Component_label
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = Lower_95,
      ymax = Upper_95
    ),
    width = 0,
    linewidth = 0.6,
    position = ggplot2::position_dodge(
      width = 0.38
    )
  ) +
  ggplot2::geom_point(
    size = 2.8,
    position = ggplot2::position_dodge(
      width = 0.38
    )
  ) +
  ggplot2::facet_wrap(
    ~ Response,
    ncol = 1,
    scales = "free_y"
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "Within-species" = "#0072B2",
      "Composition total" = "#D55E00",
      "Total change" = "#222222"
    )
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "Within-species" = 16,
      "Composition total" = 17,
      "Total change" = 18
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Dry - wet change (per mil)",
    colour = "Component",
    shape = "Component",
    caption = "Intervals are 95% individual-level bootstrap intervals."
  ) +
  theme_publication(
    base_size = 14
  )

save_publication_figure(
  figure_2,
  "Figure_2_within_species_vs_composition",
  width = 7.0,
  height = 5.8
)


# Figure 3: sensitivity to the weighting scheme
figure_3_data <- observed_long |>
  dplyr::filter(
    Component %in% c(
      "Within_species",
      "Composition_total",
      "Total_change"
    )
  ) |>
  dplyr::mutate(
    Reach = factor(
      reach_labels[
        Reach_code
      ],
      levels = REACH_LABEL_LEVELS
    ),

    Weighting_label = dplyr::recode(
      Weighting,
      sample_count = "Sample-count weighting",
      equal_species = "Equal-species weighting"
    ),

    Component_label = factor(
      component_labels[
        Component
      ],
      levels = c(
        "Within-species",
        "Composition total",
        "Total change"
      )
    )
  )

figure_3 <- ggplot2::ggplot(
  figure_3_data,
  ggplot2::aes(
    x = Reach,
    y = Observed_value,
    colour = Weighting_label,
    shape = Weighting_label
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  ggplot2::geom_point(
    size = 2.6,
    position = ggplot2::position_dodge(
      width = 0.35
    )
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(
      Response
    ),
    cols = ggplot2::vars(
      Component_label
    ),
    scales = "free_y"
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "Sample-count weighting" = "#0072B2",
      "Equal-species weighting" = "#D55E00"
    )
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "Sample-count weighting" = 16,
      "Equal-species weighting" = 17
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Dry - wet change (per mil)",
    colour = "Weighting",
    shape = "Weighting"
  ) +
  theme_publication(
    base_size = 14
  )

save_publication_figure(
  figure_3,
  "Figure_3_weighting_sensitivity",
  width = 8.2,
  height = 5.6
)


# Supplementary species-level contribution figure
species_d13C <- species_contributions[
  species_contributions$Weighting == "sample_count" &
  species_contributions$Response == "d13C",
  c(
    "Reach_code",
    "Species_code",
    "Feeding_guild",
    "Status",
    "Total_contribution"
  )
]

names(species_d13C)[
  names(species_d13C) == "Total_contribution"
] <- "d13C_total_contribution"

species_d15N <- species_contributions[
  species_contributions$Weighting == "sample_count" &
  species_contributions$Response == "d15N",
  c(
    "Reach_code",
    "Species_code",
    "Total_contribution"
  )
]

names(species_d15N)[
  names(species_d15N) == "Total_contribution"
] <- "d15N_total_contribution"

species_plot_data <- merge(
  species_d13C,
  species_d15N,
  by = c(
    "Reach_code",
    "Species_code"
  ),
  all = TRUE,
  sort = FALSE
)

species_plot_data$Reach <- factor(
  reach_labels[
    species_plot_data$Reach_code
  ],
  levels = REACH_LABEL_LEVELS
)

figure_s1 <- ggplot2::ggplot(
  species_plot_data,
  ggplot2::aes(
    x = d13C_total_contribution,
    y = d15N_total_contribution,
    shape = Status
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey45"
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey45"
  ) +
  ggplot2::geom_point(
    size = 2.2,
    alpha = 0.75
  ) +
  ggplot2::facet_wrap(
    ~ Reach,
    nrow = 1
  ) +
  ggplot2::labs(
    x = "Species contribution to d13C change (per mil)",
    y = "Species contribution to d15N change (per mil)",
    shape = "Species status"
  ) +
  theme_publication(
    base_size = 14
  )

save_publication_figure(
  figure_s1,
  "Figure_S1_species_level_total_contributions",
  width = 7.5,
  height = 4.5
)


# =============================================================================
# 12. SAVE SESSION INFORMATION AND COMPLETION LOG
# =============================================================================

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
    "Within-species versus turnover decomposition completed successfully.",
    paste("Completed:", Sys.time()),
    paste("Input workbook:", INPUT_FILE),
    paste("Bootstrap replicates:", N_BOOT),
    paste("Observed maximum closure error:", maximum_observed_closure_error),
    paste("Bootstrap maximum closure error:", maximum_bootstrap_closure_error),
    paste("Results folder:", RESULTS_DIR)
  ),
  COMPLETED_FILE
)

cat(
  "\n============================================================\n",
  "WITHIN-SPECIES VERSUS TURNOVER DECOMPOSITION COMPLETED\n",
  "============================================================\n",
  "Observed and bootstrap decompositions closed exactly.\n",
  "\nResults folder:\n",
  RESULTS_DIR,
  "\n\nMain tables:\n",
  "Table_1_observed_reach_decomposition.csv\n",
  "Table_2_bootstrap_component_intervals.csv\n",
  "Table_3_dominant_process_probabilities.csv\n",
  "Table_4_species_level_contributions.csv\n",
  "Table_5_decomposition_validation.csv\n",
  "\nMain figures:\n",
  "Figure_1_count_weighted_decomposition_components\n",
  "Figure_2_within_species_vs_composition\n",
  "Figure_3_weighting_sensitivity\n",
  "Figure_S1_species_level_total_contributions\n",
  "============================================================\n",
  sep = ""
)

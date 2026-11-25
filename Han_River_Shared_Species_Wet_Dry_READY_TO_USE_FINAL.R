# RUN IN RSTUDIO: open this R file and click Source.
# Input file(s): Han_River_Shared_Species_Wet_Dry_R_Dataset(1).xlsx
# Worksheet tabs read automatically: Shared_Species_Long
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
# SHARED-SPECIES WET-DRY BAYESIAN MULTIVARIATE COMPARISON
# FINAL READY-TO-USE WINDOWS SCRIPT
#
# PURPOSE
#   Tests wet-dry changes in d13C and d15N for species occurring in both
#   campaigns within the same reach.
#
# MODEL
#   d13C_z ~ Reach_code * Campaign + Feeding_guild
#            + (1 + Campaign | p | Pair_id)
#   d15N_z ~ Reach_code * Campaign + Feeding_guild
#            + (1 + Campaign | p | Pair_id)
#
#   Pair_id is the reach-specific species identifier.
#   Residual correlation between d13C and d15N is estimated.
#
# IMPORTANT WINDOWS STABILITY SETTING
#   Four chains are run sequentially with cores = 1. This avoids the parallel
#   RStan worker stall previously observed on this computer.
#
# INPUT WORKBOOK
#   Han_River_Shared_Species_Wet_Dry_R_Dataset(1).xlsx
#   Worksheet: Shared_Species_Long
#
# HOW TO RUN
#   1. Stop all earlier Background Jobs.
#   2. Select Session -> Restart R.
#   3. Open this script and click Source once.
#   4. Select the Excel workbook when requested.
#   5. Do not click Source again while sampling is running.
# =============================================================================


# =============================================================================
# 1. USER SETTINGS
# =============================================================================

ANALYSIS_SET <- "all"
# Allowed values:
#   "all"  = every reach-species pair present in both campaigns
#   "min2" = at least 2 fish in each campaign
#   "min3" = at least 3 fish in each campaign

FORCE_REFIT <- isTRUE(getOption("hanriver.force_refit", FALSE))

SEED <- 20260730L
CHAINS <- 4L
ITER <- 4000L
WARMUP <- 2000L
CORES <- 1L
ADAPT_DELTA <- 0.995
MAX_TREEDEPTH <- 15L

set.seed(SEED)
options(mc.cores = CORES)
graphics.off()


# =============================================================================
# 2. PACKAGE CHECK
# =============================================================================

required_packages <- c(
  "brms",
  "posterior",
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
  library(brms)
  library(posterior)
  library(readxl)
  library(dplyr)
  library(ggplot2)
})


# =============================================================================
# 3. MAIN ANALYSIS WRAPPER
# =============================================================================

run_shared_species_analysis <- function() {

  # ---------------------------------------------------------------------------
  # 3.1 Validate user settings
  # ---------------------------------------------------------------------------

  if (!(ANALYSIS_SET %in% c("all", "min2", "min3"))) {
    stop(
      'ANALYSIS_SET must be "all", "min2", or "min3".',
      call. = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # 3.2 Locate the workbook
  # ---------------------------------------------------------------------------

  INPUT_FILE <- hr_input("shared")
  hr_validate_workbook(INPUT_FILE, "shared")


  cat(
    "\nInput workbook:\n",
    INPUT_FILE,
    "\n",
    sep = ""
  )


  # ---------------------------------------------------------------------------
  # 3.3 Create output folders
  # ---------------------------------------------------------------------------

  RESULTS_DIR <- hr_results("shared")

  MODEL_DIR <- file.path(
    RESULTS_DIR,
    "Model"
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

  running_file <- file.path(
    LOG_DIR,
    "ANALYSIS_RUNNING.txt"
  )

  completed_file <- file.path(
    LOG_DIR,
    "ANALYSIS_COMPLETED.txt"
  )

  error_file <- file.path(
    LOG_DIR,
    "ANALYSIS_ERROR.txt"
  )

  unlink(
    c(
      completed_file,
      error_file
    )
  )

  writeLines(
    c(
      "Shared-species wet-dry Bayesian analysis is running.",
      paste("Started:", Sys.time()),
      paste("Analysis set:", ANALYSIS_SET),
      paste("Input file:", INPUT_FILE)
    ),
    running_file
  )


  # ---------------------------------------------------------------------------
  # 3.4 Read and validate the data
  # ---------------------------------------------------------------------------

  available_sheets <- readxl::excel_sheets(
    INPUT_FILE
  )

  if (!("Shared_Species_Long" %in% available_sheets)) {
    stop(
      paste0(
        'Worksheet "Shared_Species_Long" was not found.\n',
        "Available worksheets: ",
        paste(available_sheets, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  shared_data <- hr_read_excel(
    INPUT_FILE,
    sheet = "Shared_Species_Long"
  )

  shared_data <- as.data.frame(
    shared_data,
    stringsAsFactors = FALSE
  )

  required_columns <- c(
    "Fish_id",
    "Pair_id",
    "Reach_code",
    "Campaign",
    "Species_code",
    "Feeding_guild",
    "d13C_permil",
    "d15N_permil",
    "Eligible_min2",
    "Eligible_min3"
  )

  missing_columns <- setdiff(
    required_columns,
    names(shared_data)
  )

  if (length(missing_columns) > 0L) {
    stop(
      paste0(
        "The input worksheet is missing column(s): ",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (ANALYSIS_SET == "min2") {
    shared_data <- shared_data[
      as.logical(shared_data$Eligible_min2),
      ,
      drop = FALSE
    ]
  }

  if (ANALYSIS_SET == "min3") {
    shared_data <- shared_data[
      as.logical(shared_data$Eligible_min3),
      ,
      drop = FALSE
    ]
  }

  if (nrow(shared_data) == 0L) {
    stop(
      "No observations remained after applying the analysis-set filter.",
      call. = FALSE
    )
  }

  numeric_columns <- c(
    "d13C_permil",
    "d15N_permil"
  )

  for (column_name in numeric_columns) {
    shared_data[[column_name]] <- suppressWarnings(
      as.numeric(
        shared_data[[column_name]]
      )
    )
  }

  if (anyNA(shared_data[required_columns])) {
    stop(
      "Missing values were found in required analysis variables.",
      call. = FALSE
    )
  }

  shared_data <- shared_data |>
    dplyr::mutate(
      Reach_code = factor(
        as.character(Reach_code),
        levels = c("HC", "QJ", "SY")
      ),

      Campaign = factor(
        as.character(Campaign),
        levels = c("Wet_2023", "Dry_2023")
      ),

      Feeding_guild = factor(
        as.character(Feeding_guild),
        levels = c(
          "Omnivorous",
          "Carnivorous",
          "Herbivorous",
          "Planktivorous"
        )
      ),

      Pair_id = factor(
        as.character(Pair_id)
      ),

      Species_code = factor(
        as.character(Species_code)
      )
    )

  if (
    anyNA(shared_data$Reach_code) ||
    anyNA(shared_data$Campaign) ||
    anyNA(shared_data$Feeding_guild)
  ) {
    stop(
      "Unexpected reach, campaign, or feeding-guild values were detected.",
      call. = FALSE
    )
  }

  campaign_coverage <- shared_data |>
    dplyr::distinct(
      Pair_id,
      Reach_code,
      Campaign
    ) |>
    dplyr::count(
      Pair_id,
      Reach_code,
      name = "n_campaigns"
    )

  invalid_pairs <- campaign_coverage |>
    dplyr::filter(
      n_campaigns != 2L
    )

  if (nrow(invalid_pairs) > 0L) {
    stop(
      paste0(
        nrow(invalid_pairs),
        " Pair_id value(s) were not represented in both campaigns."
      ),
      call. = FALSE
    )
  }

  d13C_center <- mean(
    shared_data$d13C_permil
  )

  d13C_scale <- stats::sd(
    shared_data$d13C_permil
  )

  d15N_center <- mean(
    shared_data$d15N_permil
  )

  d15N_scale <- stats::sd(
    shared_data$d15N_permil
  )

  if (
    !is.finite(d13C_scale) ||
    !is.finite(d15N_scale) ||
    d13C_scale <= 0 ||
    d15N_scale <= 0
  ) {
    stop(
      "The isotope standard deviations must be finite and greater than zero.",
      call. = FALSE
    )
  }

  shared_data$d13C_z <- (
    shared_data$d13C_permil -
    d13C_center
  ) /
    d13C_scale

  shared_data$d15N_z <- (
    shared_data$d15N_permil -
    d15N_center
  ) /
    d15N_scale

  scaling_constants <- data.frame(
    response = c("d13C", "d15N"),
    center = c(d13C_center, d15N_center),
    scale = c(d13C_scale, d15N_scale)
  )

  write.csv(
    scaling_constants,
    file.path(
      TABLE_DIR,
      "Table_S1_scaling_constants.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    shared_data,
    file.path(
      TABLE_DIR,
      "analysis_data_used.csv"
    ),
    row.names = FALSE
  )

  data_summary <- data.frame(
    statistic = c(
      "Observations",
      "Reach-species pairs",
      "Unique species",
      "Reaches",
      "Campaigns"
    ),

    value = c(
      nrow(shared_data),
      nlevels(shared_data$Pair_id),
      nlevels(shared_data$Species_code),
      nlevels(droplevels(shared_data$Reach_code)),
      nlevels(droplevels(shared_data$Campaign))
    )
  )

  write.csv(
    data_summary,
    file.path(
      TABLE_DIR,
      "Table_0_data_summary.csv"
    ),
    row.names = FALSE
  )

  cat(
    "\nData validation completed.",
    "\nObservations: ", nrow(shared_data),
    "\nReach-species pairs: ", nlevels(shared_data$Pair_id),
    "\nUnique species: ", nlevels(shared_data$Species_code),
    "\n",
    sep = ""
  )


  # ---------------------------------------------------------------------------
  # 3.5 Define the Bayesian multivariate model
  # ---------------------------------------------------------------------------

  formula_d13C <- brms::bf(
    d13C_z ~
      Reach_code * Campaign +
      Feeding_guild +
      (1 + Campaign | p | Pair_id),
    family = gaussian()
  )

  formula_d15N <- brms::bf(
    d15N_z ~
      Reach_code * Campaign +
      Feeding_guild +
      (1 + Campaign | p | Pair_id),
    family = gaussian()
  )

  model_formula <- formula_d13C +
    formula_d15N +
    brms::set_rescor(TRUE)


  # ---------------------------------------------------------------------------
  # 3.6 Obtain the valid parameter structure and define response-specific priors
  #
  # IMPORTANT FIX:
  #   get_prior() is used here. It is the documented formula-based interface.
  #   Calling default_prior(formula = ...) caused the earlier
  #   "argument object is missing" error in the installed brms version.
  # ---------------------------------------------------------------------------

  available_prior_structure <- brms::get_prior(
    formula = model_formula,
    data = shared_data
  )

  write.csv(
    as.data.frame(available_prior_structure),
    file.path(
      TABLE_DIR,
      "Table_S2_available_prior_structure.csv"
    ),
    row.names = FALSE
  )

  response_ids <- unique(
    as.character(
      available_prior_structure$resp[
        !is.na(available_prior_structure$resp) &
        nzchar(
          as.character(
            available_prior_structure$resp
          )
        )
      ]
    )
  )

  response_d13C <- response_ids[
    grepl(
      "13",
      response_ids,
      ignore.case = TRUE
    )
  ][1]

  response_d15N <- response_ids[
    grepl(
      "15",
      response_ids,
      ignore.case = TRUE
    )
  ][1]

  if (
    length(response_d13C) == 0L ||
    length(response_d15N) == 0L ||
    is.na(response_d13C) ||
    is.na(response_d15N)
  ) {
    stop(
      paste0(
        "The internal response names could not be identified.\n",
        "Detected response identifiers: ",
        paste(response_ids, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  model_priors <- c(
    brms::set_prior(
      "normal(0, 0.75)",
      class = "b",
      resp = response_d13C
    ),

    brms::set_prior(
      "normal(0, 1)",
      class = "Intercept",
      resp = response_d13C
    ),

    brms::set_prior(
      "normal(0, 0.75)",
      class = "b",
      resp = response_d15N
    ),

    brms::set_prior(
      "normal(0, 1)",
      class = "Intercept",
      resp = response_d15N
    )
  )

  validated_priors <- brms::validate_prior(
    prior = model_priors,
    formula = model_formula,
    data = shared_data
  )

  write.csv(
    as.data.frame(validated_priors),
    file.path(
      TABLE_DIR,
      "Table_S3_validated_priors.csv"
    ),
    row.names = FALSE
  )

  cat(
    "\nPrior validation completed successfully.",
    "\nInternal d13C response: ", response_d13C,
    "\nInternal d15N response: ", response_d15N,
    "\n",
    sep = ""
  )


  # ---------------------------------------------------------------------------
  # 3.7 Fit or load the model
  # ---------------------------------------------------------------------------

  MODEL_FILE <- file.path(
    MODEL_DIR,
    paste0(
      "Shared_species_Bayesian_MV_",
      ANALYSIS_SET,
      ".rds"
    )
  )

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
      inherits(loaded_model, "brmsfit")
    )
  }

  if (valid_saved_model) {

    fit_shared <- loaded_model

    cat(
      "\nLoaded the existing fitted model:\n",
      MODEL_FILE,
      "\n",
      sep = ""
    )

  } else {

    if (file.exists(MODEL_FILE)) {
      unlink(MODEL_FILE)
    }

    cat(
      "\nBeginning Bayesian sampling.",
      "\nFour chains will run sequentially with cores = 1.",
      "\nDo not click Source again while the model is running.\n",
      sep = ""
    )

    fit_shared <- brms::brm(
      formula = model_formula,
      data = shared_data,
      prior = model_priors,
      backend = "rstan",
      chains = CHAINS,
      cores = CORES,
      iter = ITER,
      warmup = WARMUP,
      seed = SEED,
      control = list(
        adapt_delta = ADAPT_DELTA,
        max_treedepth = MAX_TREEDEPTH
      ),
      refresh = 100,
      save_pars = brms::save_pars(
        all = TRUE
      )
    )

    saveRDS(
      fit_shared,
      MODEL_FILE
    )
  }


  # ---------------------------------------------------------------------------
  # 3.8 Complete model summary
  # ---------------------------------------------------------------------------

  model_summary <- summary(
    fit_shared
  )

  writeLines(
    capture.output(
      model_summary
    ),
    file.path(
      TABLE_DIR,
      "Table_S4_complete_model_summary.txt"
    )
  )


  # ---------------------------------------------------------------------------
  # 3.9 Posterior draws and fixed effects
  # ---------------------------------------------------------------------------

  posterior_draws <- posterior::as_draws_df(
    fit_shared
  )

  fixed_matrix <- brms::fixef(
    fit_shared,
    summary = TRUE,
    probs = c(0.025, 0.975)
  )

  fixed_effects <- data.frame(
    parameter = rownames(fixed_matrix),
    estimate = fixed_matrix[, "Estimate"],
    posterior_sd = fixed_matrix[, "Est.Error"],
    lower_95 = fixed_matrix[, "Q2.5"],
    upper_95 = fixed_matrix[, "Q97.5"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  fixed_effects$probability_positive <- vapply(
    fixed_effects$parameter,
    function(parameter_name) {

      draw_name <- paste0(
        "b_",
        parameter_name
      )

      if (!(draw_name %in% names(posterior_draws))) {
        return(NA_real_)
      }

      values <- as.numeric(
        posterior_draws[[draw_name]]
      )

      mean(
        values > 0,
        na.rm = TRUE
      )
    },
    FUN.VALUE = numeric(1)
  )

  fixed_effects$probability_negative <- vapply(
    fixed_effects$parameter,
    function(parameter_name) {

      draw_name <- paste0(
        "b_",
        parameter_name
      )

      if (!(draw_name %in% names(posterior_draws))) {
        return(NA_real_)
      }

      values <- as.numeric(
        posterior_draws[[draw_name]]
      )

      mean(
        values < 0,
        na.rm = TRUE
      )
    },
    FUN.VALUE = numeric(1)
  )

  fixed_effects$probability_direction <- pmax(
    fixed_effects$probability_positive,
    fixed_effects$probability_negative,
    na.rm = TRUE
  )

  fixed_effects$credible_95 <- (
    fixed_effects$lower_95 > 0 |
    fixed_effects$upper_95 < 0
  )

  write.csv(
    fixed_effects,
    file.path(
      TABLE_DIR,
      "Table_1_fixed_effects_and_posterior_probabilities.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.10 Variance components
  # ---------------------------------------------------------------------------

  variance_components <- brms::VarCorr(
    fit_shared,
    summary = TRUE,
    probs = c(0.025, 0.975)
  )

  writeLines(
    capture.output(
      print(
        variance_components
      )
    ),
    file.path(
      TABLE_DIR,
      "Table_2_pair_level_variance_components.txt"
    )
  )


  # ---------------------------------------------------------------------------
  # 3.11 Bayesian R-squared
  # ---------------------------------------------------------------------------

  r2_d13C <- as.data.frame(
    brms::bayes_R2(
      fit_shared,
      resp = response_d13C,
      summary = TRUE,
      probs = c(0.025, 0.975)
    )
  )

  r2_d13C$response <- "d13C"

  r2_d15N <- as.data.frame(
    brms::bayes_R2(
      fit_shared,
      resp = response_d15N,
      summary = TRUE,
      probs = c(0.025, 0.975)
    )
  )

  r2_d15N$response <- "d15N"

  bayesian_r2 <- dplyr::bind_rows(
    r2_d13C,
    r2_d15N
  )

  write.csv(
    bayesian_r2,
    file.path(
      TABLE_DIR,
      "Table_3_Bayesian_R2.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.12 Residual correlation
  # ---------------------------------------------------------------------------

  residual_parameter_names <- grep(
    "rescor",
    names(posterior_draws),
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(residual_parameter_names) > 0L) {

    residual_parameter <- residual_parameter_names[1]

    residual_values <- as.numeric(
      posterior_draws[[residual_parameter]]
    )

    residual_summary <- data.frame(
      parameter = residual_parameter,

      posterior_mean = mean(
        residual_values,
        na.rm = TRUE
      ),

      posterior_median = stats::median(
        residual_values,
        na.rm = TRUE
      ),

      posterior_sd = stats::sd(
        residual_values,
        na.rm = TRUE
      ),

      lower_95 = unname(
        stats::quantile(
          residual_values,
          probs = 0.025,
          na.rm = TRUE
        )
      ),

      upper_95 = unname(
        stats::quantile(
          residual_values,
          probs = 0.975,
          na.rm = TRUE
        )
      ),

      probability_positive = mean(
        residual_values > 0,
        na.rm = TRUE
      ),

      probability_negative = mean(
        residual_values < 0,
        na.rm = TRUE
      )
    )

    write.csv(
      residual_summary,
      file.path(
        TABLE_DIR,
        "Table_4_residual_correlation.csv"
      ),
      row.names = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # 3.13 Convergence and sampler diagnostics
  # ---------------------------------------------------------------------------

  parameter_diagnostics <- posterior::summarise_draws(
    posterior_draws
  )

  finite_rhat <- parameter_diagnostics$rhat[
    is.finite(
      parameter_diagnostics$rhat
    )
  ]

  finite_bulk_ess <- parameter_diagnostics$ess_bulk[
    is.finite(
      parameter_diagnostics$ess_bulk
    )
  ]

  finite_tail_ess <- parameter_diagnostics$ess_tail[
    is.finite(
      parameter_diagnostics$ess_tail
    )
  ]

  nuts_diagnostics <- brms::nuts_params(
    fit_shared
  )

  divergences <- sum(
    nuts_diagnostics$Parameter == "divergent__" &
    nuts_diagnostics$Value == 1
  )

  treedepth_hits <- sum(
    nuts_diagnostics$Parameter == "treedepth__" &
    nuts_diagnostics$Value >= MAX_TREEDEPTH
  )

  diagnostic_summary <- data.frame(
    diagnostic = c(
      "Maximum Rhat",
      "Minimum bulk ESS",
      "Minimum tail ESS",
      "Divergent transitions",
      "Maximum-treedepth hits"
    ),

    value = c(
      max(
        finite_rhat
      ),

      min(
        finite_bulk_ess
      ),

      min(
        finite_tail_ess
      ),

      divergences,

      treedepth_hits
    )
  )

  write.csv(
    diagnostic_summary,
    file.path(
      TABLE_DIR,
      "Table_5_convergence_diagnostics.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.14 Helper functions for predictions and posterior contrasts
  # ---------------------------------------------------------------------------

  as_prediction_matrix <- function(
    prediction_object
  ) {

    object_dimensions <- dim(
      prediction_object
    )

    if (length(object_dimensions) == 2L) {
      return(
        prediction_object
      )
    }

    if (
      length(object_dimensions) == 3L &&
      object_dimensions[3] == 1L
    ) {
      return(
        prediction_object[, , 1]
      )
    }

    stop(
      paste0(
        "Unexpected posterior prediction dimensions: ",
        paste(object_dimensions, collapse = " x ")
      ),
      call. = FALSE
    )
  }

  summarise_change_draws <- function(
    change_draws,
    response,
    comparison,
    reach_code = NA_character_,
    pair_id = NA_character_,
    species_code = NA_character_,
    feeding_guild = NA_character_
  ) {

    data.frame(
      response = response,
      comparison = comparison,
      Reach_code = reach_code,
      Pair_id = pair_id,
      Species_code = species_code,
      Feeding_guild = feeding_guild,

      posterior_mean = mean(
        change_draws,
        na.rm = TRUE
      ),

      posterior_median = stats::median(
        change_draws,
        na.rm = TRUE
      ),

      lower_95 = unname(
        stats::quantile(
          change_draws,
          probs = 0.025,
          na.rm = TRUE
        )
      ),

      upper_95 = unname(
        stats::quantile(
          change_draws,
          probs = 0.975,
          na.rm = TRUE
        )
      ),

      probability_positive = mean(
        change_draws > 0,
        na.rm = TRUE
      ),

      probability_negative = mean(
        change_draws < 0,
        na.rm = TRUE
      ),

      stringsAsFactors = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # 3.15 Reach-specific adjusted dry-minus-wet contrasts
  # ---------------------------------------------------------------------------

  guild_weights <- prop.table(
    table(
      shared_data$Feeding_guild
    )
  )

  reach_prediction_grid <- expand.grid(
    Reach_code = levels(
      shared_data$Reach_code
    ),

    Campaign = levels(
      shared_data$Campaign
    ),

    Feeding_guild = levels(
      shared_data$Feeding_guild
    ),

    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  reach_prediction_grid <- reach_prediction_grid |>
    dplyr::mutate(
      Reach_code = factor(
        Reach_code,
        levels = levels(
          shared_data$Reach_code
        )
      ),

      Campaign = factor(
        Campaign,
        levels = levels(
          shared_data$Campaign
        )
      ),

      Feeding_guild = factor(
        Feeding_guild,
        levels = levels(
          shared_data$Feeding_guild
        )
      ),

      Pair_id = factor(
        levels(
          shared_data$Pair_id
        )[1],
        levels = levels(
          shared_data$Pair_id
        )
      )
    )

  calculate_reach_changes <- function(
    response_id,
    raw_center,
    raw_scale,
    response_label
  ) {

    expected_draws <- brms::posterior_epred(
      fit_shared,
      newdata = reach_prediction_grid,
      resp = response_id,
      re_formula = NA,
      ndraws = 2000
    )

    expected_draws <- as_prediction_matrix(
      expected_draws
    )

    reach_results <- list()

    for (reach_code in levels(shared_data$Reach_code)) {

      campaign_draws <- list()

      for (campaign_name in levels(shared_data$Campaign)) {

        selected_columns <- which(
          reach_prediction_grid$Reach_code == reach_code &
          reach_prediction_grid$Campaign == campaign_name
        )

        selected_guilds <- as.character(
          reach_prediction_grid$Feeding_guild[
            selected_columns
          ]
        )

        selected_weights <- as.numeric(
          guild_weights[
            selected_guilds
          ]
        )

        selected_weights <- selected_weights /
          sum(
            selected_weights
          )

        standardized_average <- as.numeric(
          expected_draws[
            ,
            selected_columns,
            drop = FALSE
          ] %*%
          selected_weights
        )

        campaign_draws[[campaign_name]] <- (
          standardized_average *
          raw_scale
        ) +
          raw_center
      }

      dry_minus_wet <- (
        campaign_draws[["Dry_2023"]] -
        campaign_draws[["Wet_2023"]]
      )

      reach_results[[reach_code]] <- summarise_change_draws(
        change_draws = dry_minus_wet,
        response = response_label,
        comparison = "Dry - wet",
        reach_code = reach_code
      )
    }

    dplyr::bind_rows(
      reach_results
    )
  }

  reach_changes <- dplyr::bind_rows(
    calculate_reach_changes(
      response_id = response_d13C,
      raw_center = d13C_center,
      raw_scale = d13C_scale,
      response_label = "d13C"
    ),

    calculate_reach_changes(
      response_id = response_d15N,
      raw_center = d15N_center,
      raw_scale = d15N_scale,
      response_label = "d15N"
    )
  )

  reach_changes$probability_direction <- pmax(
    reach_changes$probability_positive,
    reach_changes$probability_negative
  )

  reach_changes$credible_95 <- (
    reach_changes$lower_95 > 0 |
    reach_changes$upper_95 < 0
  )

  write.csv(
    reach_changes,
    file.path(
      TABLE_DIR,
      "Table_6_reach_specific_dry_minus_wet_contrasts_permil.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.16 Reach-species-pair-specific dry-minus-wet changes
  # ---------------------------------------------------------------------------

  pair_information <- shared_data |>
    dplyr::distinct(
      Pair_id,
      Reach_code,
      Species_code,
      Feeding_guild
    ) |>
    dplyr::arrange(
      Pair_id
    )

  pair_prediction_grid <- dplyr::bind_rows(
    pair_information |>
      dplyr::mutate(
        Campaign = "Wet_2023"
      ),

    pair_information |>
      dplyr::mutate(
        Campaign = "Dry_2023"
      )
  ) |>
    dplyr::mutate(
      Reach_code = factor(
        as.character(Reach_code),
        levels = levels(
          shared_data$Reach_code
        )
      ),

      Campaign = factor(
        Campaign,
        levels = levels(
          shared_data$Campaign
        )
      ),

      Feeding_guild = factor(
        as.character(Feeding_guild),
        levels = levels(
          shared_data$Feeding_guild
        )
      ),

      Pair_id = factor(
        as.character(Pair_id),
        levels = levels(
          shared_data$Pair_id
        )
      ),

      Species_code = factor(
        as.character(Species_code),
        levels = levels(
          shared_data$Species_code
        )
      )
    ) |>
    dplyr::arrange(
      Pair_id,
      Campaign
    )

  calculate_pair_changes <- function(
    response_id,
    raw_scale,
    response_label
  ) {

    pair_expected_draws <- brms::posterior_epred(
      fit_shared,
      newdata = pair_prediction_grid,
      resp = response_id,
      re_formula = NULL,
      allow_new_levels = FALSE,
      ndraws = 2000
    )

    pair_expected_draws <- as_prediction_matrix(
      pair_expected_draws
    )

    pair_results <- vector(
      mode = "list",
      length = nrow(pair_information)
    )

    pair_names <- as.character(
      pair_information$Pair_id
    )

    for (pair_index in seq_along(pair_names)) {

      current_pair <- pair_names[pair_index]

      wet_column <- which(
        pair_prediction_grid$Pair_id == current_pair &
        pair_prediction_grid$Campaign == "Wet_2023"
      )

      dry_column <- which(
        pair_prediction_grid$Pair_id == current_pair &
        pair_prediction_grid$Campaign == "Dry_2023"
      )

      if (
        length(wet_column) != 1L ||
        length(dry_column) != 1L
      ) {
        stop(
          paste0(
            "Prediction rows could not be identified for Pair_id: ",
            current_pair
          ),
          call. = FALSE
        )
      }

      standardized_change <- (
        pair_expected_draws[, dry_column] -
        pair_expected_draws[, wet_column]
      )

      raw_change <- standardized_change *
        raw_scale

      pair_results[[pair_index]] <- summarise_change_draws(
        change_draws = raw_change,
        response = response_label,
        comparison = "Dry - wet",
        reach_code = as.character(
          pair_information$Reach_code[pair_index]
        ),
        pair_id = current_pair,
        species_code = as.character(
          pair_information$Species_code[pair_index]
        ),
        feeding_guild = as.character(
          pair_information$Feeding_guild[pair_index]
        )
      )
    }

    dplyr::bind_rows(
      pair_results
    )
  }

  pair_changes <- dplyr::bind_rows(
    calculate_pair_changes(
      response_id = response_d13C,
      raw_scale = d13C_scale,
      response_label = "d13C"
    ),

    calculate_pair_changes(
      response_id = response_d15N,
      raw_scale = d15N_scale,
      response_label = "d15N"
    )
  )

  pair_changes$probability_direction <- pmax(
    pair_changes$probability_positive,
    pair_changes$probability_negative
  )

  pair_changes$credible_95 <- (
    pair_changes$lower_95 > 0 |
    pair_changes$upper_95 < 0
  )

  write.csv(
    pair_changes,
    file.path(
      TABLE_DIR,
      "Table_7_pair_specific_dry_minus_wet_changes_permil.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.17 Publication figures
  # ---------------------------------------------------------------------------

  colour_d13C <- "#0072B2"
  colour_d15N <- "#D55E00"

  reach_labels <- c(
    HC = "S1",
    QJ = "S2",
    SY = "S3"
  )

  theme_publication <- function(base_size = 14) {

    ggplot2::theme_classic(
      base_size = base_size,
      base_family = "sans"
    ) +
      ggplot2::theme(
        axis.text = ggplot2::element_text(
          colour = "black"
        ),
        strip.background = ggplot2::element_rect(
          fill = "grey94",
          colour = "black",
          linewidth = 0.4
        ),
        strip.text = ggplot2::element_text(
          face = "bold"
        ),
        legend.position = "top",
        plot.margin = ggplot2::margin(
          8,
          10,
          8,
          8
        )
      ) + hr_typography()
}

  save_figure <- function(
    plot_object,
    filename,
    width,
    height
  ) {

    hr_ggsave(
      filename = file.path(
        FIGURE_DIR,
        paste0(filename, ".png")
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
        paste0(filename, ".tiff")
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
            paste0(filename, ".pdf")
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
          paste0(filename, ".pdf")
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

  reach_plot_data <- reach_changes |>
    dplyr::mutate(
      Reach = dplyr::recode(
        Reach_code,
        HC = "S1",
        QJ = "S2",
        SY = "S3"
      ),

      Reach = factor(
        Reach,
        levels = c(
          "S1",
          "S2",
          "S3"
        )
      ),

      response_label = dplyr::recode(
        response,
        d13C = "d13C",
        d15N = "d15N"
      )
    )

  figure_1 <- ggplot2::ggplot(
    reach_plot_data,
    ggplot2::aes(
      x = Reach,
      y = posterior_median
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.45,
      colour = "grey35"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = lower_95,
        ymax = upper_95
      ),
      width = 0,
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = response_label
      ),
      shape = 21,
      size = 3.1,
      stroke = 0.6,
      colour = "black"
    ) +
    ggplot2::facet_wrap(
      ~ response_label,
      ncol = 1,
      scales = "free_y"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        d13C = colour_d13C,
        d15N = colour_d15N
      ),
      guide = "none"
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Posterior dry - wet change (per mil)",
      caption = "Points are posterior medians; intervals are 95% credible intervals."
    ) +
    theme_publication(
      base_size = 14
    )

  save_figure(
    figure_1,
    "Figure_1_reach_specific_dry_minus_wet_changes",
    width = 6.8,
    height = 6.0
  )

  pair_wide <- merge(
    pair_changes[
      pair_changes$response == "d13C",
      c(
        "Pair_id",
        "Reach_code",
        "Species_code",
        "Feeding_guild",
        "posterior_median"
      )
    ],

    pair_changes[
      pair_changes$response == "d15N",
      c(
        "Pair_id",
        "posterior_median"
      )
    ],

    by = "Pair_id",
    suffixes = c(
      "_d13C",
      "_d15N"
    )
  )

  pair_wide$Reach <- factor(
    reach_labels[
      pair_wide$Reach_code
    ],
    levels = c(
      "S1",
      "S2",
      "S3"
    )
  )

  figure_2 <- ggplot2::ggplot(
    pair_wide,
    ggplot2::aes(
      x = posterior_median_d13C,
      y = posterior_median_d15N,
      shape = Feeding_guild
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
      size = 2.4,
      alpha = 0.8
    ) +
    ggplot2::facet_wrap(
      ~ Reach,
      nrow = 1
    ) +
    ggplot2::labs(
      x = "Pair-specific dry - wet d13C change (per mil)",
      y = "Pair-specific dry - wet d15N change (per mil)",
      shape = "Feeding guild"
    ) +
    theme_publication(
      base_size = 14
    )

  save_figure(
    figure_2,
    "Figure_2_pair_specific_bivariate_changes",
    width = 7.4,
    height = 4.4
  )


  # ---------------------------------------------------------------------------
  # 3.18 Posterior predictive density checks
  # ---------------------------------------------------------------------------

  create_ppc_density <- function(
    response_id,
    observed_values,
    response_label,
    filename
  ) {

    set.seed(
      SEED
    )

    replicated_values <- brms::posterior_predict(
      fit_shared,
      resp = response_id,
      ndraws = 30
    )

    replicated_values <- as_prediction_matrix(
      replicated_values
    )

    value_range <- range(
      c(
        observed_values,
        as.numeric(replicated_values)
      ),
      finite = TRUE
    )

    observed_density <- stats::density(
      observed_values,
      from = value_range[1],
      to = value_range[2],
      n = 512
    )

    observed_density_data <- data.frame(
      x = observed_density$x,
      density = observed_density$y
    )

    replicated_density_data <- dplyr::bind_rows(
      lapply(
        seq_len(
          nrow(replicated_values)
        ),
        function(draw_index) {

          current_density <- stats::density(
            replicated_values[draw_index, ],
            from = value_range[1],
            to = value_range[2],
            n = 512
          )

          data.frame(
            x = current_density$x,
            density = current_density$y,
            draw = draw_index
          )
        }
      )
    )

    ppc_plot <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = replicated_density_data,
        mapping = ggplot2::aes(
          x = x,
          y = density,
          group = draw
        ),
        linewidth = 0.35,
        alpha = 0.35,
        colour = "#9ECAE1"
      ) +
      ggplot2::geom_line(
        data = observed_density_data,
        mapping = ggplot2::aes(
          x = x,
          y = density
        ),
        linewidth = 1.0,
        colour = "#162A43"
      ) +
      ggplot2::labs(
        x = paste(
          "Standardized",
          response_label
        ),
        y = "Density",
        caption = "Dark line: observed distribution; pale lines: posterior predictive distributions."
      ) +
      theme_publication(
        base_size = 14
      )

    save_figure(
      ppc_plot,
      filename,
      width = 6.6,
      height = 4.6
    )
  }

  create_ppc_density(
    response_id = response_d13C,
    observed_values = shared_data$d13C_z,
    response_label = "d13C",
    filename = "Figure_S1_PPC_density_d13C"
  )

  create_ppc_density(
    response_id = response_d15N,
    observed_values = shared_data$d15N_z,
    response_label = "d15N",
    filename = "Figure_S2_PPC_density_d15N"
  )


  # ---------------------------------------------------------------------------
  # 3.19 Save session information and completion status
  # ---------------------------------------------------------------------------

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
    running_file
  )

  writeLines(
    c(
      "Shared-species wet-dry Bayesian analysis completed successfully.",
      paste("Completed:", Sys.time()),
      paste("Analysis set:", ANALYSIS_SET),
      paste("Input file:", INPUT_FILE),
      paste("Model file:", MODEL_FILE),
      paste("Results folder:", RESULTS_DIR)
    ),
    completed_file
  )

  cat(
    "\n============================================================\n",
    "SHARED-SPECIES WET-DRY ANALYSIS COMPLETED\n",
    "============================================================\n",
    "Four chains were run sequentially with cores = 1.\n",
    "The priors were validated using brms::get_prior().\n",
    "\nResults folder:\n",
    RESULTS_DIR,
    "\n\nMain outputs:\n",
    "Table_1_fixed_effects_and_posterior_probabilities.csv\n",
    "Table_2_pair_level_variance_components.txt\n",
    "Table_3_Bayesian_R2.csv\n",
    "Table_4_residual_correlation.csv\n",
    "Table_5_convergence_diagnostics.csv\n",
    "Table_6_reach_specific_dry_minus_wet_contrasts_permil.csv\n",
    "Table_7_pair_specific_dry_minus_wet_changes_permil.csv\n",
    "Figure_1_reach_specific_dry_minus_wet_changes\n",
    "Figure_2_pair_specific_bivariate_changes\n",
    "============================================================\n",
    sep = ""
  )

  invisible(
    list(
      fit = fit_shared,
      reach_changes = reach_changes,
      pair_changes = pair_changes,
      results_directory = RESULTS_DIR
    )
  )
}


# =============================================================================
# 4. RUN THE ANALYSIS AND SAVE A USEFUL ERROR LOG IF IT FAILS
# =============================================================================

analysis_result <- tryCatch(
  run_shared_species_analysis(),

  error = function(e) {

    error_message <- paste0(
      "Shared-species analysis stopped with an error.\n",
      "Time: ",
      Sys.time(),
      "\n\n",
      conditionMessage(e),
      "\n"
    )

    message(
      "\n",
      error_message
    )

    # Save an error file in the current working directory as a fallback.
    fallback_error_file <- file.path(
      getwd(),
      "SHARED_SPECIES_ANALYSIS_ERROR.txt"
    )

    writeLines(
      error_message,
      fallback_error_file
    )

    stop(
      conditionMessage(e),
      call. = FALSE
    )
  }
)

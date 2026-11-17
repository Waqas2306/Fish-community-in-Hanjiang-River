# =============================================================================
# LOWER HAN RIVER FISH STABLE ISOTOPES
# BAYESIAN TROPHIC-POSITION ANALYSIS
# FINAL READY-TO-USE WINDOWS SCRIPT
#
# PURPOSE
#   Calculates and compares fish trophic position among reaches and campaigns,
#   while propagating uncertainty in:
#     1. fish d15N estimated by a hierarchical Bayesian model;
#     2. baseline d15N; and
#     3. the nitrogen trophic discrimination factor (Delta15N).
#
# TROPHIC-POSITION EQUATION
#   TP = lambda + (d15N_fish - d15N_baseline) / Delta15N
#
# IMPORTANT PUBLICATION NOTE
#   The supplied workbook contains provisional herbivorous-fish baseline
#   proxies because no measured primary-consumer baseline was supplied.
#   The workflow will run, but absolute trophic-position values should only be
#   reported as final ecological estimates after the baseline is independently
#   justified or replaced with measured primary-consumer data.
#
# INPUT WORKBOOK
#   Han_River_Trophic_Position_R_Dataset.xlsx
#
# PRIMARY MODEL
#   d15N_z ~ Reach_code * Campaign + Feeding_guild
#            + (1 + Campaign | p | Species_code)
#
# WINDOWS STABILITY
#   Four chains run sequentially with cores = 1.
#
# HOW TO RUN
#   1. Stop earlier Background Jobs.
#   2. Select Session -> Restart R.
#   3. Open this script and click Source once.
#   4. Select the workbook when requested.
#   5. Do not click Source again while sampling is running.
# =============================================================================


# =============================================================================
# 1. USER SETTINGS
# =============================================================================

FORCE_REFIT <- FALSE

SEED <- 20260731L
CHAINS <- 4L
CORES <- 1L
ITER <- 4000L
WARMUP <- 2000L
ADAPT_DELTA <- 0.995
MAX_TREEDEPTH <- 15L

MONTE_CARLO_DRAWS <- 4000L

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
      "\nInstall them, restart RStudio, and run the script again."
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

run_trophic_position_analysis <- function() {

  # ---------------------------------------------------------------------------
  # 3.1 Locate the workbook
  # ---------------------------------------------------------------------------

  workbook_pattern <- "^Han_River_Trophic_Position_R_Dataset.*[.]xlsx$"

  user_home <- path.expand("~")
  windows_user <- Sys.getenv("USERPROFILE")

  search_folders <- unique(
    c(
      getwd(),
      file.path(windows_user, "Downloads"),
      file.path(user_home, "Downloads"),
      file.path(windows_user, "Desktop"),
      file.path(user_home, "Desktop"),
      file.path(windows_user, "Documents"),
      file.path(user_home, "Documents")
    )
  )

  search_folders <- search_folders[
    nzchar(search_folders) &
    dir.exists(search_folders)
  ]

  candidate_files <- character(0)

  for (folder_path in search_folders) {

    detected_files <- list.files(
      path = folder_path,
      pattern = workbook_pattern,
      full.names = TRUE,
      recursive = FALSE,
      ignore.case = TRUE
    )

    candidate_files <- c(
      candidate_files,
      detected_files
    )
  }

  candidate_files <- unique(
    candidate_files[
      file.exists(candidate_files)
    ]
  )

  if (length(candidate_files) > 1L) {

    candidate_information <- file.info(
      candidate_files
    )

    candidate_files <- candidate_files[
      order(
        candidate_information$mtime,
        decreasing = TRUE,
        na.last = TRUE
      )
    ]
  }

  if (length(candidate_files) == 0L) {

    message(
      "\nSelect the trophic-position Excel workbook."
    )

    selected_file <- tryCatch(
      file.choose(),
      error = function(e) {
        ""
      }
    )

    if (
      !nzchar(selected_file) ||
      !file.exists(selected_file)
    ) {
      stop(
        "No valid Excel workbook was selected.",
        call. = FALSE
      )
    }

    candidate_files <- selected_file
  }

  INPUT_FILE <- normalizePath(
    candidate_files[1],
    winslash = "/",
    mustWork = TRUE
  )

  cat(
    "\nInput workbook:\n",
    INPUT_FILE,
    "\n",
    sep = ""
  )


  # ---------------------------------------------------------------------------
  # 3.2 Create output folders
  # ---------------------------------------------------------------------------

  RESULTS_DIR <- file.path(
    dirname(INPUT_FILE),
    "Trophic_position_results"
  )

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

  warning_file <- file.path(
    LOG_DIR,
    "PUBLICATION_BASELINE_WARNING.txt"
  )

  unlink(
    c(
      completed_file,
      error_file,
      warning_file
    )
  )

  writeLines(
    c(
      "Bayesian trophic-position analysis is running.",
      paste("Started:", Sys.time()),
      paste("Input workbook:", INPUT_FILE)
    ),
    running_file
  )


  # ---------------------------------------------------------------------------
  # 3.3 Read and validate workbook sheets
  # ---------------------------------------------------------------------------

  available_sheets <- readxl::excel_sheets(
    INPUT_FILE
  )

  required_sheets <- c(
    "Fish_Data",
    "Baseline_Inputs",
    "TDF_Inputs"
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

  fish_data <- readxl::read_excel(
    INPUT_FILE,
    sheet = "Fish_Data"
  )

  baseline_inputs <- readxl::read_excel(
    INPUT_FILE,
    sheet = "Baseline_Inputs"
  )

  tdf_inputs <- readxl::read_excel(
    INPUT_FILE,
    sheet = "TDF_Inputs"
  )

  fish_data <- as.data.frame(
    fish_data,
    stringsAsFactors = FALSE
  )

  baseline_inputs <- as.data.frame(
    baseline_inputs,
    stringsAsFactors = FALSE
  )

  tdf_inputs <- as.data.frame(
    tdf_inputs,
    stringsAsFactors = FALSE
  )

  fish_required <- c(
    "Fish_id",
    "Reach_code",
    "Campaign",
    "Species_code",
    "Feeding_guild",
    "d15N_permil"
  )

  baseline_required <- c(
    "Scenario_ID",
    "Scenario_label",
    "Baseline_type",
    "Reach_code",
    "Campaign",
    "Baseline_mean_d15N",
    "Baseline_SD_d15N",
    "Baseline_TP_lambda",
    "n_baseline_species",
    "Active",
    "Primary_scenario",
    "Publication_status"
  )

  tdf_required <- c(
    "TDF_ID",
    "TDF_label",
    "Delta15N_mean",
    "Delta15N_SD",
    "Active",
    "Primary_TDF"
  )

  missing_fish_columns <- setdiff(
    fish_required,
    names(fish_data)
  )

  missing_baseline_columns <- setdiff(
    baseline_required,
    names(baseline_inputs)
  )

  missing_tdf_columns <- setdiff(
    tdf_required,
    names(tdf_inputs)
  )

  if (length(missing_fish_columns) > 0L) {
    stop(
      paste0(
        "Fish_Data is missing: ",
        paste(missing_fish_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (length(missing_baseline_columns) > 0L) {
    stop(
      paste0(
        "Baseline_Inputs is missing: ",
        paste(missing_baseline_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (length(missing_tdf_columns) > 0L) {
    stop(
      paste0(
        "TDF_Inputs is missing: ",
        paste(missing_tdf_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  fish_data$d15N_permil <- suppressWarnings(
    as.numeric(
      fish_data$d15N_permil
    )
  )

  if (anyNA(
    fish_data[
      ,
      fish_required
    ]
  )) {
    stop(
      "Missing or invalid values were found in Fish_Data.",
      call. = FALSE
    )
  }

  logical_from_excel <- function(
    values
  ) {

    if (is.logical(
      values
    )) {
      return(
        values
      )
    }

    values_character <- toupper(
      trimws(
        as.character(
          values
        )
      )
    )

    values_character %in% c(
      "TRUE",
      "T",
      "1",
      "YES",
      "Y"
    )
  }

  baseline_inputs$Active <- logical_from_excel(
    baseline_inputs$Active
  )

  baseline_inputs$Primary_scenario <- logical_from_excel(
    baseline_inputs$Primary_scenario
  )

  tdf_inputs$Active <- logical_from_excel(
    tdf_inputs$Active
  )

  tdf_inputs$Primary_TDF <- logical_from_excel(
    tdf_inputs$Primary_TDF
  )

  baseline_inputs$Baseline_mean_d15N <- suppressWarnings(
    as.numeric(
      baseline_inputs$Baseline_mean_d15N
    )
  )

  baseline_inputs$Baseline_SD_d15N <- suppressWarnings(
    as.numeric(
      baseline_inputs$Baseline_SD_d15N
    )
  )

  baseline_inputs$Baseline_TP_lambda <- suppressWarnings(
    as.numeric(
      baseline_inputs$Baseline_TP_lambda
    )
  )

  baseline_inputs$n_baseline_species <- suppressWarnings(
    as.numeric(
      baseline_inputs$n_baseline_species
    )
  )

  tdf_inputs$Delta15N_mean <- suppressWarnings(
    as.numeric(
      tdf_inputs$Delta15N_mean
    )
  )

  tdf_inputs$Delta15N_SD <- suppressWarnings(
    as.numeric(
      tdf_inputs$Delta15N_SD
    )
  )

  active_baselines <- baseline_inputs[
    baseline_inputs$Active,
    ,
    drop = FALSE
  ]

  active_tdfs <- tdf_inputs[
    tdf_inputs$Active,
    ,
    drop = FALSE
  ]

  if (nrow(
    active_baselines
  ) == 0L) {
    stop(
      "No active baseline scenarios were found.",
      call. = FALSE
    )
  }

  if (nrow(
    active_tdfs
  ) == 0L) {
    stop(
      "No active TDF scenarios were found.",
      call. = FALSE
    )
  }

  if (anyNA(
    active_baselines[
      ,
      c(
        "Baseline_mean_d15N",
        "Baseline_SD_d15N",
        "Baseline_TP_lambda",
        "n_baseline_species"
      )
    ]
  )) {
    stop(
      paste0(
        "At least one active baseline row has missing numerical values. ",
        "Complete or deactivate that row."
      ),
      call. = FALSE
    )
  }

  if (any(
    active_baselines$n_baseline_species <= 0
  )) {
    stop(
      "n_baseline_species must be greater than zero for active baselines.",
      call. = FALSE
    )
  }

  if (anyNA(
    active_tdfs[
      ,
      c(
        "Delta15N_mean",
        "Delta15N_SD"
      )
    ]
  )) {
    stop(
      "At least one active TDF row has missing numerical values.",
      call. = FALSE
    )
  }

  if (any(
    active_tdfs$Delta15N_mean <= 0
  )) {
    stop(
      "Active Delta15N means must be greater than zero.",
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

  expected_guilds <- c(
    "Omnivorous",
    "Carnivorous",
    "Herbivorous",
    "Planktivorous"
  )

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
        as.character(Feeding_guild),
        levels = expected_guilds
      )
    )

  if (
    anyNA(
      fish_data$Reach_code
    ) ||
    anyNA(
      fish_data$Campaign
    ) ||
    anyNA(
      fish_data$Feeding_guild
    )
  ) {
    stop(
      "Unexpected reach, campaign, or feeding-guild values were found.",
      call. = FALSE
    )
  }

  baseline_coverage <- active_baselines |>
    dplyr::distinct(
      Scenario_ID,
      Reach_code,
      Campaign
    ) |>
    dplyr::count(
      Scenario_ID,
      name = "n_groups"
    )

  incomplete_scenarios <- baseline_coverage |>
    dplyr::filter(
      n_groups != 6L
    )

  if (nrow(
    incomplete_scenarios
  ) > 0L) {
    stop(
      paste0(
        "Every active baseline scenario must contain all six reach-campaign rows. ",
        "Incomplete scenario(s): ",
        paste(
          incomplete_scenarios$Scenario_ID,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }

  primary_baseline_ids <- unique(
    active_baselines$Scenario_ID[
      active_baselines$Primary_scenario
    ]
  )

  primary_tdf_ids <- unique(
    active_tdfs$TDF_ID[
      active_tdfs$Primary_TDF
    ]
  )

  if (length(
    primary_baseline_ids
  ) != 1L) {
    stop(
      "Exactly one active baseline scenario must be flagged Primary_scenario=TRUE.",
      call. = FALSE
    )
  }

  if (length(
    primary_tdf_ids
  ) != 1L) {
    stop(
      "Exactly one active TDF must be flagged Primary_TDF=TRUE.",
      call. = FALSE
    )
  }

  primary_baseline_id <- primary_baseline_ids[1]
  primary_tdf_id <- primary_tdf_ids[1]

  primary_baseline_rows <- active_baselines[
    active_baselines$Scenario_ID == primary_baseline_id,
    ,
    drop = FALSE
  ]

  primary_baseline_type <- unique(
    primary_baseline_rows$Baseline_type
  )

  primary_is_provisional <- any(
    grepl(
      "Provisional",
      primary_baseline_type,
      ignore.case = TRUE
    )
  )

  if (primary_is_provisional) {

    writeLines(
      c(
        "PUBLICATION BASELINE WARNING",
        "",
        paste("Primary baseline scenario:", primary_baseline_id),
        paste("Baseline type:", paste(primary_baseline_type, collapse = "; ")),
        "",
        paste(
          "Absolute trophic-position values are conditional on a provisional",
          "herbivorous-fish proxy. Replace or independently justify the baseline",
          "before treating the estimates as publication-ready absolute trophic positions."
        )
      ),
      warning_file
    )
  }

  cat(
    "\nData validation completed.",
    "\nFish observations: ",
    nrow(
      fish_data
    ),
    "\nActive baseline scenarios: ",
    length(
      unique(
        active_baselines$Scenario_ID
      )
    ),
    "\nActive TDF scenarios: ",
    nrow(
      active_tdfs
    ),
    "\nPrimary baseline: ",
    primary_baseline_id,
    "\nPrimary TDF: ",
    primary_tdf_id,
    "\n",
    sep = ""
  )


  # ---------------------------------------------------------------------------
  # 3.4 Standardize d15N
  # ---------------------------------------------------------------------------

  d15N_center <- mean(
    fish_data$d15N_permil
  )

  d15N_scale <- stats::sd(
    fish_data$d15N_permil
  )

  if (
    !is.finite(
      d15N_scale
    ) ||
    d15N_scale <= 0
  ) {
    stop(
      "The d15N standard deviation must be greater than zero.",
      call. = FALSE
    )
  }

  fish_data$d15N_z <- (
    fish_data$d15N_permil -
    d15N_center
  ) /
    d15N_scale

  write.csv(
    data.frame(
      response = "d15N",
      center = d15N_center,
      scale = d15N_scale
    ),
    file.path(
      TABLE_DIR,
      "Table_S1_scaling_constant.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.5 Define and validate the hierarchical Bayesian model
  # ---------------------------------------------------------------------------

  model_formula <- brms::bf(
    d15N_z ~
      Reach_code * Campaign +
      Feeding_guild +
      (1 + Campaign | p | Species_code),
    family = gaussian()
  )

  model_priors <- c(
    brms::set_prior(
      "normal(0, 0.75)",
      class = "b"
    ),
    brms::set_prior(
      "normal(0, 1)",
      class = "Intercept"
    ),
    brms::set_prior(
      "exponential(1)",
      class = "sd"
    ),
    brms::set_prior(
      "exponential(1)",
      class = "sigma"
    )
  )

  validated_priors <- brms::validate_prior(
    prior = model_priors,
    formula = model_formula,
    data = fish_data
  )

  write.csv(
    as.data.frame(
      validated_priors
    ),
    file.path(
      TABLE_DIR,
      "Table_S2_validated_priors.csv"
    ),
    row.names = FALSE
  )

  cat(
    "\nPrior validation completed successfully.\n"
  )


  # ---------------------------------------------------------------------------
  # 3.6 Fit or load the Bayesian model
  # ---------------------------------------------------------------------------

  model_file <- file.path(
    MODEL_DIR,
    "Fish_d15N_hierarchical_model.rds"
  )

  valid_saved_model <- FALSE

  if (
    file.exists(
      model_file
    ) &&
    !isTRUE(
      FORCE_REFIT
    )
  ) {

    loaded_model <- try(
      readRDS(
        model_file
      ),
      silent = TRUE
    )

    valid_saved_model <- (
      !inherits(
        loaded_model,
        "try-error"
      ) &&
      inherits(
        loaded_model,
        "brmsfit"
      )
    )
  }

  if (valid_saved_model) {

    fit_d15N <- loaded_model

    cat(
      "\nLoaded existing model:\n",
      model_file,
      "\n",
      sep = ""
    )

  } else {

    if (file.exists(
      model_file
    )) {
      unlink(
        model_file
      )
    }

    cat(
      "\nBeginning Bayesian sampling.",
      "\nFour chains will run sequentially with cores = 1.",
      "\nDo not run the script again while sampling is active.\n",
      sep = ""
    )

    fit_d15N <- brms::brm(
      formula = model_formula,
      data = fish_data,
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
      fit_d15N,
      model_file
    )
  }


  # ---------------------------------------------------------------------------
  # 3.7 Model diagnostics
  # ---------------------------------------------------------------------------

  posterior_draws <- posterior::as_draws_df(
    fit_d15N
  )

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
    fit_d15N
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
    Diagnostic = c(
      "Maximum Rhat",
      "Minimum bulk ESS",
      "Minimum tail ESS",
      "Divergent transitions",
      "Maximum-treedepth hits"
    ),

    Value = c(
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
      "Table_1_model_diagnostics.csv"
    ),
    row.names = FALSE
  )

  writeLines(
    capture.output(
      summary(
        fit_d15N
      )
    ),
    file.path(
      TABLE_DIR,
      "Table_S3_complete_model_summary.txt"
    )
  )


  # ---------------------------------------------------------------------------
  # 3.8 Population-level adjusted d15N posterior draws
  # ---------------------------------------------------------------------------

  guild_weights <- prop.table(
    table(
      fish_data$Feeding_guild
    )
  )

  prediction_grid <- expand.grid(
    Reach_code = levels(
      fish_data$Reach_code
    ),
    Campaign = levels(
      fish_data$Campaign
    ),
    Feeding_guild = levels(
      fish_data$Feeding_guild
    ),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  prediction_grid <- prediction_grid |>
    dplyr::mutate(
      Reach_code = factor(
        Reach_code,
        levels = levels(
          fish_data$Reach_code
        )
      ),

      Campaign = factor(
        Campaign,
        levels = levels(
          fish_data$Campaign
        )
      ),

      Feeding_guild = factor(
        Feeding_guild,
        levels = levels(
          fish_data$Feeding_guild
        )
      ),

      Species_code = factor(
        levels(
          fish_data$Species_code
        )[1],
        levels = levels(
          fish_data$Species_code
        )
      )
    )

  epred_standardized <- brms::posterior_epred(
    fit_d15N,
    newdata = prediction_grid,
    re_formula = NA,
    ndraws = MONTE_CARLO_DRAWS
  )

  if (
    length(
      dim(
        epred_standardized
      )
    ) == 3L
  ) {
    epred_standardized <- epred_standardized[
      ,
      ,
      1
    ]
  }

  if (
    !is.matrix(
      epred_standardized
    )
  ) {
    epred_standardized <- as.matrix(
      epred_standardized
    )
  }

  group_keys <- expand.grid(
    Reach_code = levels(
      fish_data$Reach_code
    ),
    Campaign = levels(
      fish_data$Campaign
    ),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  adjusted_d15N_draws <- matrix(
    NA_real_,
    nrow = nrow(
      epred_standardized
    ),
    ncol = nrow(
      group_keys
    )
  )

  colnames(
    adjusted_d15N_draws
  ) <- paste(
    group_keys$Reach_code,
    group_keys$Campaign,
    sep = "__"
  )

  for (group_index in seq_len(
    nrow(
      group_keys
    )
  )) {

    reach_code <- group_keys$Reach_code[
      group_index
    ]

    campaign_name <- group_keys$Campaign[
      group_index
    ]

    selected_columns <- which(
      prediction_grid$Reach_code == reach_code &
      prediction_grid$Campaign == campaign_name
    )

    selected_guilds <- as.character(
      prediction_grid$Feeding_guild[
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

    adjusted_standardized <- as.numeric(
      epred_standardized[
        ,
        selected_columns,
        drop = FALSE
      ] %*%
      selected_weights
    )

    adjusted_d15N_draws[
      ,
      group_index
    ] <- (
      adjusted_standardized *
      d15N_scale
    ) +
      d15N_center
  }


  # ---------------------------------------------------------------------------
  # 3.9 Helper functions for sensitivity propagation
  # ---------------------------------------------------------------------------

  summarise_draws <- function(
    values
  ) {

    values <- as.numeric(
      values
    )

    data.frame(
      Posterior_mean = mean(
        values,
        na.rm = TRUE
      ),

      Posterior_median = stats::median(
        values,
        na.rm = TRUE
      ),

      Lower_95 = unname(
        stats::quantile(
          values,
          probs = 0.025,
          na.rm = TRUE
        )
      ),

      Upper_95 = unname(
        stats::quantile(
          values,
          probs = 0.975,
          na.rm = TRUE
        )
      ),

      Probability_positive = mean(
        values > 0,
        na.rm = TRUE
      ),

      Probability_negative = mean(
        values < 0,
        na.rm = TRUE
      )
    )
  }


  baseline_draw_key <- function(
    scenario_id,
    reach_code,
    campaign_name
  ) {

    if (scenario_id == "Reach_Pooled_Herbivore_Proxy") {
      return(
        paste(
          scenario_id,
          reach_code,
          sep = "__"
        )
      )
    }

    if (scenario_id == "Global_Herbivore_Proxy") {
      return(
        scenario_id
      )
    }

    paste(
      scenario_id,
      reach_code,
      campaign_name,
      sep = "__"
    )
  }


  draw_positive_normal <- function(
    n,
    mean_value,
    sd_value,
    minimum_value = 0.50
  ) {

    draws <- stats::rnorm(
      n = n,
      mean = mean_value,
      sd = sd_value
    )

    invalid <- (
      !is.finite(
        draws
      ) |
      draws <= minimum_value
    )

    while (any(
      invalid
    )) {

      draws[
        invalid
      ] <- stats::rnorm(
        n = sum(
          invalid
        ),
        mean = mean_value,
        sd = sd_value
      )

      invalid <- (
        !is.finite(
          draws
        ) |
        draws <= minimum_value
      )
    }

    draws
  }


  # ---------------------------------------------------------------------------
  # 3.10 Calculate trophic-position draws for all active scenarios
  # ---------------------------------------------------------------------------

  scenario_ids <- unique(
    active_baselines$Scenario_ID
  )

  tdf_ids <- unique(
    active_tdfs$TDF_ID
  )

  tp_draw_storage <- list()
  group_summary_rows <- list()
  group_summary_index <- 1L

  set.seed(
    SEED + 10L
  )

  for (scenario_id in scenario_ids) {

    scenario_rows <- active_baselines[
      active_baselines$Scenario_ID == scenario_id,
      ,
      drop = FALSE
    ]

    baseline_draw_cache <- list()

    for (group_index in seq_len(
      nrow(
        group_keys
      )
    )) {

      reach_code <- group_keys$Reach_code[
        group_index
      ]

      campaign_name <- group_keys$Campaign[
        group_index
      ]

      baseline_row <- scenario_rows[
        scenario_rows$Reach_code == reach_code &
        scenario_rows$Campaign == campaign_name,
        ,
        drop = FALSE
      ]

      if (nrow(
        baseline_row
      ) != 1L) {
        stop(
          paste0(
            "Expected one baseline row for ",
            scenario_id,
            ", ",
            reach_code,
            ", ",
            campaign_name,
            "."
          ),
          call. = FALSE
        )
      }

      cache_key <- baseline_draw_key(
        scenario_id,
        reach_code,
        campaign_name
      )

      if (!(cache_key %in% names(
        baseline_draw_cache
      ))) {

        baseline_se <- (
          baseline_row$Baseline_SD_d15N /
          sqrt(
            baseline_row$n_baseline_species
          )
        )

        baseline_draw_cache[[
          cache_key
        ]] <- stats::rnorm(
          n = nrow(
            adjusted_d15N_draws
          ),
          mean = baseline_row$Baseline_mean_d15N,
          sd = baseline_se
        )
      }
    }

    for (tdf_id in tdf_ids) {

      tdf_row <- active_tdfs[
        active_tdfs$TDF_ID == tdf_id,
        ,
        drop = FALSE
      ]

      if (nrow(
        tdf_row
      ) != 1L) {
        stop(
          paste0(
            "Expected exactly one active row for TDF_ID ",
            tdf_id,
            "."
          ),
          call. = FALSE
        )
      }

      tdf_draws <- draw_positive_normal(
        n = nrow(
          adjusted_d15N_draws
        ),
        mean_value = tdf_row$Delta15N_mean,
        sd_value = tdf_row$Delta15N_SD,
        minimum_value = 0.50
      )

      current_tp_matrix <- matrix(
        NA_real_,
        nrow = nrow(
          adjusted_d15N_draws
        ),
        ncol = ncol(
          adjusted_d15N_draws
        )
      )

      colnames(
        current_tp_matrix
      ) <- colnames(
        adjusted_d15N_draws
      )

      for (group_index in seq_len(
        nrow(
          group_keys
        )
      )) {

        reach_code <- group_keys$Reach_code[
          group_index
        ]

        campaign_name <- group_keys$Campaign[
          group_index
        ]

        baseline_row <- scenario_rows[
          scenario_rows$Reach_code == reach_code &
          scenario_rows$Campaign == campaign_name,
          ,
          drop = FALSE
        ]

        cache_key <- baseline_draw_key(
          scenario_id,
          reach_code,
          campaign_name
        )

        baseline_draws <- baseline_draw_cache[[
          cache_key
        ]]

        tp_draws <- (
          baseline_row$Baseline_TP_lambda +
          (
            adjusted_d15N_draws[
              ,
              group_index
            ] -
            baseline_draws
          ) /
          tdf_draws
        )

        current_tp_matrix[
          ,
          group_index
        ] <- tp_draws

        current_summary <- summarise_draws(
          tp_draws
        )

        group_summary_rows[[
          group_summary_index
        ]] <- data.frame(
          Scenario_ID = scenario_id,
          Scenario_label = baseline_row$Scenario_label,
          Baseline_type = baseline_row$Baseline_type,
          TDF_ID = tdf_id,
          TDF_label = tdf_row$TDF_label,
          Reach_code = reach_code,
          Campaign = campaign_name,
          Baseline_mean_d15N = baseline_row$Baseline_mean_d15N,
          Baseline_SE_d15N = (
            baseline_row$Baseline_SD_d15N /
            sqrt(
              baseline_row$n_baseline_species
            )
          ),
          Baseline_TP_lambda = baseline_row$Baseline_TP_lambda,
          Delta15N_mean = tdf_row$Delta15N_mean,
          Delta15N_SD = tdf_row$Delta15N_SD,
          current_summary,
          stringsAsFactors = FALSE
        )

        group_summary_index <- group_summary_index + 1L
      }

      storage_key <- paste(
        scenario_id,
        tdf_id,
        sep = "__"
      )

      tp_draw_storage[[
        storage_key
      ]] <- current_tp_matrix
    }
  }

  all_group_tp_summary <- dplyr::bind_rows(
    group_summary_rows
  )

  write.csv(
    all_group_tp_summary,
    file.path(
      TABLE_DIR,
      "Table_2_all_sensitivity_group_TP.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.11 Primary-scenario group estimates and contrasts
  # ---------------------------------------------------------------------------

  primary_key <- paste(
    primary_baseline_id,
    primary_tdf_id,
    sep = "__"
  )

  primary_tp_matrix <- tp_draw_storage[[
    primary_key
  ]]

  primary_group_summary <- all_group_tp_summary |>
    dplyr::filter(
      Scenario_ID == primary_baseline_id,
      TDF_ID == primary_tdf_id
    ) |>
    dplyr::mutate(
      Reach_name = dplyr::recode(
        Reach_code,
        HC = "Hanchuan",
        QJ = "Qianjiang",
        SY = "Shayang"
      ),

      Campaign_label = dplyr::recode(
        Campaign,
        Wet_2023 = "Wet campaign",
        Dry_2023 = "Dry campaign"
      )
    )

  write.csv(
    primary_group_summary,
    file.path(
      TABLE_DIR,
      "Table_3_primary_group_trophic_positions.csv"
    ),
    row.names = FALSE
  )


  contrast_rows <- list()
  contrast_index <- 1L

  # Within-reach dry minus wet
  for (reach_code in expected_reaches) {

    wet_key <- paste(
      reach_code,
      "Wet_2023",
      sep = "__"
    )

    dry_key <- paste(
      reach_code,
      "Dry_2023",
      sep = "__"
    )

    contrast_draws <- (
      primary_tp_matrix[
        ,
        dry_key
      ] -
      primary_tp_matrix[
        ,
        wet_key
      ]
    )

    current_summary <- summarise_draws(
      contrast_draws
    )

    contrast_rows[[
      contrast_index
    ]] <- data.frame(
      Contrast_type = "Within-reach seasonal",
      Reach_or_campaign = reach_code,
      Group_1 = wet_key,
      Group_2 = dry_key,
      Contrast = "Dry - wet",
      current_summary,
      stringsAsFactors = FALSE
    )

    contrast_index <- contrast_index + 1L
  }

  seasonal_contrasts <- dplyr::bind_rows(
    contrast_rows
  )

  write.csv(
    seasonal_contrasts,
    file.path(
      TABLE_DIR,
      "Table_4_primary_wet_to_dry_TP_contrasts.csv"
    ),
    row.names = FALSE
  )


  # Among-reach contrasts within each campaign
  reach_pairs <- combn(
    expected_reaches,
    2,
    simplify = FALSE
  )

  reach_contrast_rows <- list()
  reach_contrast_index <- 1L

  for (campaign_name in expected_campaigns) {

    for (reach_pair in reach_pairs) {

      group_1 <- paste(
        reach_pair[1],
        campaign_name,
        sep = "__"
      )

      group_2 <- paste(
        reach_pair[2],
        campaign_name,
        sep = "__"
      )

      contrast_draws <- (
        primary_tp_matrix[
          ,
          group_2
        ] -
        primary_tp_matrix[
          ,
          group_1
        ]
      )

      current_summary <- summarise_draws(
        contrast_draws
      )

      reach_contrast_rows[[
        reach_contrast_index
      ]] <- data.frame(
        Contrast_type = "Among-reach same campaign",
        Reach_or_campaign = campaign_name,
        Group_1 = group_1,
        Group_2 = group_2,
        Contrast = paste0(
          reach_pair[2],
          " - ",
          reach_pair[1]
        ),
        current_summary,
        stringsAsFactors = FALSE
      )

      reach_contrast_index <- reach_contrast_index + 1L
    }
  }

  among_reach_contrasts <- dplyr::bind_rows(
    reach_contrast_rows
  )

  write.csv(
    among_reach_contrasts,
    file.path(
      TABLE_DIR,
      "Table_5_primary_among_reach_TP_contrasts.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.12 Sensitivity of seasonal change
  # ---------------------------------------------------------------------------

  sensitivity_contrast_rows <- list()
  sensitivity_index <- 1L

  for (storage_key in names(
    tp_draw_storage
  )) {

    key_parts <- strsplit(
      storage_key,
      "__",
      fixed = TRUE
    )[[1]]

    scenario_id <- key_parts[1]
    tdf_id <- paste(
      key_parts[
        -1
      ],
      collapse = "__"
    )

    current_matrix <- tp_draw_storage[[
      storage_key
    ]]

    for (reach_code in expected_reaches) {

      wet_key <- paste(
        reach_code,
        "Wet_2023",
        sep = "__"
      )

      dry_key <- paste(
        reach_code,
        "Dry_2023",
        sep = "__"
      )

      contrast_draws <- (
        current_matrix[
          ,
          dry_key
        ] -
        current_matrix[
          ,
          wet_key
        ]
      )

      current_summary <- summarise_draws(
        contrast_draws
      )

      sensitivity_contrast_rows[[
        sensitivity_index
      ]] <- data.frame(
        Scenario_ID = scenario_id,
        TDF_ID = tdf_id,
        Reach_code = reach_code,
        Contrast = "Dry - wet",
        current_summary,
        stringsAsFactors = FALSE
      )

      sensitivity_index <- sensitivity_index + 1L
    }
  }

  sensitivity_seasonal_contrasts <- dplyr::bind_rows(
    sensitivity_contrast_rows
  )

  write.csv(
    sensitivity_seasonal_contrasts,
    file.path(
      TABLE_DIR,
      "Table_6_sensitivity_wet_to_dry_TP_contrasts.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.13 Descriptive individual TP under the primary assumptions
  # ---------------------------------------------------------------------------

  primary_baseline_lookup <- primary_baseline_rows[
    ,
    c(
      "Reach_code",
      "Campaign",
      "Baseline_mean_d15N",
      "Baseline_TP_lambda"
    ),
    drop = FALSE
  ]

  primary_tdf_row <- active_tdfs[
    active_tdfs$TDF_ID == primary_tdf_id,
    ,
    drop = FALSE
  ]

  individual_tp <- merge(
    fish_data,
    primary_baseline_lookup,
    by = c(
      "Reach_code",
      "Campaign"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  individual_tp$TP_point_estimate <- (
    individual_tp$Baseline_TP_lambda +
    (
      individual_tp$d15N_permil -
      individual_tp$Baseline_mean_d15N
    ) /
    primary_tdf_row$Delta15N_mean
  )

  individual_tp$Primary_baseline_scenario <- primary_baseline_id
  individual_tp$Primary_TDF <- primary_tdf_id

  write.csv(
    individual_tp,
    file.path(
      TABLE_DIR,
      "Table_7_individual_TP_point_estimates.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.14 Bayesian R-squared and variance components
  # ---------------------------------------------------------------------------

  bayesian_r2 <- as.data.frame(
    brms::bayes_R2(
      fit_d15N,
      summary = TRUE,
      probs = c(
        0.025,
        0.975
      )
    )
  )

  write.csv(
    bayesian_r2,
    file.path(
      TABLE_DIR,
      "Table_S4_Bayesian_R2.csv"
    ),
    row.names = FALSE
  )

  writeLines(
    capture.output(
      brms::VarCorr(
        fit_d15N,
        summary = TRUE,
        probs = c(
          0.025,
          0.975
        )
      )
    ),
    file.path(
      TABLE_DIR,
      "Table_S5_species_variance_components.txt"
    )
  )


  # ---------------------------------------------------------------------------
  # 3.15 Publication figures
  # ---------------------------------------------------------------------------

  reach_labels <- c(
    HC = "Hanchuan",
    QJ = "Qianjiang",
    SY = "Shayang"
  )

  campaign_colours <- c(
    "Wet campaign" = "#0072B2",
    "Dry campaign" = "#D55E00"
  )

  theme_publication <- function(
    base_size = 10
  ) {

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
      )
  }


  save_figure <- function(
    plot_object,
    filename,
    width,
    height
  ) {

    ggplot2::ggsave(
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

    ggplot2::ggsave(
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
        ggplot2::ggsave(
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

      ggplot2::ggsave(
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


  # Figure 1: observed individual TP distributions
  individual_plot_data <- individual_tp |>
    dplyr::mutate(
      Reach_name = factor(
        reach_labels[
          as.character(
            Reach_code
          )
        ],
        levels = c(
          "Hanchuan",
          "Qianjiang",
          "Shayang"
        )
      ),

      Campaign_label = factor(
        dplyr::recode(
          as.character(
            Campaign
          ),
          Wet_2023 = "Wet campaign",
          Dry_2023 = "Dry campaign"
        ),
        levels = c(
          "Wet campaign",
          "Dry campaign"
        )
      )
    )

  figure_1 <- ggplot2::ggplot(
    individual_plot_data,
    ggplot2::aes(
      x = Campaign_label,
      y = TP_point_estimate,
      fill = Campaign_label
    )
  ) +
    ggplot2::geom_violin(
      trim = FALSE,
      alpha = 0.45,
      colour = "grey30",
      linewidth = 0.45
    ) +
    ggplot2::geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      alpha = 0.75,
      linewidth = 0.45
    ) +
    ggplot2::facet_wrap(
      ~ Reach_name,
      nrow = 1
    ) +
    ggplot2::scale_fill_manual(
      values = campaign_colours,
      guide = "none"
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Individual trophic-position estimate",
      caption = paste(
        "Point estimates use the primary baseline and TDF assumptions;",
        "absolute values remain conditional on baseline validity."
      )
    ) +
    theme_publication(
      base_size = 9.5
    )

  save_figure(
    figure_1,
    "Figure_1_observed_individual_trophic_positions",
    width = 7.5,
    height = 4.8
  )


  # Figure 2: adjusted group TP estimates
  group_plot_data <- primary_group_summary |>
    dplyr::mutate(
      Reach_name = factor(
        Reach_name,
        levels = c(
          "Hanchuan",
          "Qianjiang",
          "Shayang"
        )
      ),

      Campaign_label = factor(
        Campaign_label,
        levels = c(
          "Wet campaign",
          "Dry campaign"
        )
      )
    )

  figure_2 <- ggplot2::ggplot(
    group_plot_data,
    ggplot2::aes(
      x = Reach_name,
      y = Posterior_median,
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
      linewidth = 0.65,
      position = ggplot2::position_dodge(
        width = 0.38
      )
    ) +
    ggplot2::geom_point(
      size = 3,
      position = ggplot2::position_dodge(
        width = 0.38
      )
    ) +
    ggplot2::scale_colour_manual(
      values = campaign_colours
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "Wet campaign" = 16,
        "Dry campaign" = 17
      )
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Adjusted trophic position",
      colour = "Campaign",
      shape = "Campaign",
      caption = "Points and intervals are posterior medians and 95% credible intervals."
    ) +
    theme_publication(
      base_size = 10
    )

  save_figure(
    figure_2,
    "Figure_2_adjusted_group_trophic_positions",
    width = 6.8,
    height = 4.7
  )


  # Figure 3: seasonal contrasts
  seasonal_plot_data <- seasonal_contrasts |>
    dplyr::mutate(
      Reach_name = factor(
        reach_labels[
          Reach_or_campaign
        ],
        levels = c(
          "Hanchuan",
          "Qianjiang",
          "Shayang"
        )
      )
    )

  figure_3 <- ggplot2::ggplot(
    seasonal_plot_data,
    ggplot2::aes(
      x = Reach_name,
      y = Posterior_median
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.45,
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
      size = 3,
      shape = 21,
      fill = "white",
      colour = "black"
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Dry - wet trophic-position change",
      caption = "Points and intervals are posterior medians and 95% credible intervals."
    ) +
    theme_publication(
      base_size = 10
    )

  save_figure(
    figure_3,
    "Figure_3_wet_to_dry_trophic_position_changes",
    width = 6.3,
    height = 4.4
  )


  # Figure 4: baseline sensitivity under primary TDF
  baseline_sensitivity_data <- all_group_tp_summary |>
    dplyr::filter(
      TDF_ID == primary_tdf_id
    ) |>
    dplyr::mutate(
      Group_label = paste0(
        Reach_code,
        "-",
        ifelse(
          Campaign == "Wet_2023",
          "Wet",
          "Dry"
        )
      ),

      Group_label = factor(
        Group_label,
        levels = c(
          "HC-Wet",
          "HC-Dry",
          "QJ-Wet",
          "QJ-Dry",
          "SY-Wet",
          "SY-Dry"
        )
      )
    )

  figure_4 <- ggplot2::ggplot(
    baseline_sensitivity_data,
    ggplot2::aes(
      x = Group_label,
      y = Posterior_median,
      group = Scenario_ID,
      colour = Scenario_ID
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      size = 2.5
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Median trophic position",
      colour = "Baseline scenario",
      caption = paste(
        "Sensitivity at the primary TDF;",
        "differences among lines quantify baseline dependence."
      )
    ) +
    theme_publication(
      base_size = 9.5
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 35,
        hjust = 1
      )
    )

  save_figure(
    figure_4,
    "Figure_4_baseline_sensitivity",
    width = 7.5,
    height = 4.9
  )


  # Figure 5: TDF sensitivity under primary baseline
  tdf_sensitivity_data <- all_group_tp_summary |>
    dplyr::filter(
      Scenario_ID == primary_baseline_id
    ) |>
    dplyr::mutate(
      Group_label = paste0(
        Reach_code,
        "-",
        ifelse(
          Campaign == "Wet_2023",
          "Wet",
          "Dry"
        )
      ),

      Group_label = factor(
        Group_label,
        levels = c(
          "HC-Wet",
          "HC-Dry",
          "QJ-Wet",
          "QJ-Dry",
          "SY-Wet",
          "SY-Dry"
        )
      )
    )

  figure_5 <- ggplot2::ggplot(
    tdf_sensitivity_data,
    ggplot2::aes(
      x = Group_label,
      y = Posterior_median,
      group = TDF_ID,
      colour = TDF_ID
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      size = 2.5
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Median trophic position",
      colour = "TDF scenario",
      caption = paste(
        "Sensitivity under the primary baseline;",
        "larger Delta15N values compress trophic-position differences."
      )
    ) +
    theme_publication(
      base_size = 9.5
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 35,
        hjust = 1
      )
    )

  save_figure(
    figure_5,
    "Figure_5_discrimination_factor_sensitivity",
    width = 7.5,
    height = 4.9
  )


  # Posterior predictive check for d15N model
  set.seed(
    SEED + 20L
  )

  pp_check_plot <- brms::pp_check(
    fit_d15N,
    type = "dens_overlay",
    ndraws = 30
  ) +
    ggplot2::labs(
      x = "Standardized d15N",
      y = "Density"
    ) +
    theme_publication(
      base_size = 10
    )

  save_figure(
    pp_check_plot,
    "Figure_S1_posterior_predictive_check_d15N",
    width = 6.5,
    height = 4.6
  )


  # ---------------------------------------------------------------------------
  # 3.16 Optional combined Excel output
  # ---------------------------------------------------------------------------

  if (
    requireNamespace(
      "writexl",
      quietly = TRUE
    )
  ) {

    writexl::write_xlsx(
      list(
        Model_Diagnostics = diagnostic_summary,
        All_Sensitivity_Group_TP = all_group_tp_summary,
        Primary_Group_TP = primary_group_summary,
        Seasonal_Contrasts = seasonal_contrasts,
        Among_Reach_Contrasts = among_reach_contrasts,
        Sensitivity_Seasonal = sensitivity_seasonal_contrasts
      ),
      path = file.path(
        RESULTS_DIR,
        "Trophic_position_results.xlsx"
      )
    )
  }


  # ---------------------------------------------------------------------------
  # 3.17 Save session information and completion status
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
      "Bayesian trophic-position analysis completed successfully.",
      paste("Completed:", Sys.time()),
      paste("Input workbook:", INPUT_FILE),
      paste("Model file:", model_file),
      paste("Primary baseline:", primary_baseline_id),
      paste("Primary TDF:", primary_tdf_id),
      paste("Results folder:", RESULTS_DIR),
      paste("Primary baseline provisional:", primary_is_provisional)
    ),
    completed_file
  )

  cat(
    "\n============================================================\n",
    "BAYESIAN TROPHIC-POSITION ANALYSIS COMPLETED\n",
    "============================================================\n",
    "\nResults folder:\n",
    RESULTS_DIR,
    "\n\nMain tables:\n",
    "Table_1_model_diagnostics.csv\n",
    "Table_2_all_sensitivity_group_TP.csv\n",
    "Table_3_primary_group_trophic_positions.csv\n",
    "Table_4_primary_wet_to_dry_TP_contrasts.csv\n",
    "Table_5_primary_among_reach_TP_contrasts.csv\n",
    "Table_6_sensitivity_wet_to_dry_TP_contrasts.csv\n",
    "Table_7_individual_TP_point_estimates.csv\n",
    "\nMain figures:\n",
    "Figure_1_observed_individual_trophic_positions\n",
    "Figure_2_adjusted_group_trophic_positions\n",
    "Figure_3_wet_to_dry_trophic_position_changes\n",
    "Figure_4_baseline_sensitivity\n",
    "Figure_5_discrimination_factor_sensitivity\n",
    "Figure_S1_posterior_predictive_check_d15N\n",
    "\nPublication baseline warning: ",
    ifelse(
      primary_is_provisional,
      "YES — read Logs/PUBLICATION_BASELINE_WARNING.txt",
      "No"
    ),
    "\n============================================================\n",
    sep = ""
  )

  invisible(
    list(
      fit = fit_d15N,
      primary_group_summary = primary_group_summary,
      seasonal_contrasts = seasonal_contrasts,
      sensitivity_summary = all_group_tp_summary,
      results_directory = RESULTS_DIR
    )
  )
}


# =============================================================================
# 4. RUN ANALYSIS AND SAVE A CLEAR ERROR LOG
# =============================================================================

analysis_result <- tryCatch(
  run_trophic_position_analysis(),

  error = function(e) {

    error_message <- paste0(
      "Bayesian trophic-position analysis stopped with an error.\n",
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

    fallback_error_file <- file.path(
      getwd(),
      "TROPHIC_POSITION_ANALYSIS_ERROR.txt"
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

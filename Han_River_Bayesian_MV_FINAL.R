# RUN IN RSTUDIO: open this R file and click Source.
# Input file(s): Bayesian Multivariate analysis.xlsx
# Worksheet tabs read automatically: Model_Data
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
# HAN RIVER FISH STABLE ISOTOPES
# FINAL BAYESIAN MULTIVARIATE HIERARCHICAL ANALYSIS
#
# Create a NEW R script, paste this entire code, save it, and click Source.
# The package launcher runs this module in an isolated environment.
# =============================================================================

graphics.off()
set.seed(20260729)

# ---- File locations and model settings ---------------------------------------

DATA_FILE <- hr_input("bayesian")
DATA_SHEET <- "Model_Data"
OUTPUT_DIR <- hr_results("main_model")

SEED <- 20260729L
CHAINS <- 4L
CORES <- 1L
ITER <- 5000L
WARMUP <- 2500L
ADAPT_DELTA <- 0.995
MAX_TREEDEPTH <- 15L

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Packages ---------------------------------------------------------------

required_packages <- c(
  "readxl", "dplyr", "ggplot2", "brms",
  "rstan", "posterior", "pkgbuild"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE,
          FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0L) {
  stop(
    paste0(
      "Missing package(s): ", paste(missing_packages, collapse = ", "),
      "\nInstall them, restart RStudio, and run this script again."
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(brms)
  library(rstan)
  library(posterior)
})

options(mc.cores = CORES, rstan.auto_write = TRUE)

if (!isTRUE(pkgbuild::has_build_tools(debug = FALSE))) {
  stop(
    "Rtools is not available to R. Run pkgbuild::has_build_tools(debug = TRUE).",
    call. = FALSE
  )
}

cat("\nRTOOLS CHECK: PASSED\n")

# ---- Import and validate the workbook ---------------------------------------

if (!file.exists(DATA_FILE)) {
  stop(paste0("Excel file not found:\n", DATA_FILE), call. = FALSE)
}

available_sheets <- readxl::excel_sheets(DATA_FILE)

if (!(DATA_SHEET %in% available_sheets)) {
  stop(
    paste0(
      "Worksheet '", DATA_SHEET, "' not found.\nAvailable sheets: ",
      paste(available_sheets, collapse = ", ")
    ),
    call. = FALSE
  )
}

exact_data <- hr_read_excel(DATA_FILE, sheet = DATA_SHEET)

required_columns <- c(
  "Fish_id", "Reach_code", "Campaign", "Species_code",
  "Feeding_guild", "d13C_permil", "d15N_permil"
)

missing_columns <- setdiff(required_columns, names(exact_data))

if (length(missing_columns) > 0L) {
  stop(
    paste0("Missing column(s): ", paste(missing_columns, collapse = ", ")),
    call. = FALSE
  )
}

fish_model <- exact_data |>
  dplyr::transmute(
    Fish_id = trimws(as.character(Fish_id)),
    Reach_code = trimws(as.character(Reach_code)),
    Campaign = trimws(as.character(Campaign)),
    Species_code = trimws(as.character(Species_code)),
    Feeding_guild = trimws(as.character(Feeding_guild)),
    d13C_permil = suppressWarnings(as.numeric(d13C_permil)),
    d15N_permil = suppressWarnings(as.numeric(d15N_permil))
  )

missing_summary <- data.frame(
  variable = names(fish_model),
  missing_values = vapply(fish_model, function(x) sum(is.na(x)), numeric(1))
)

write.csv(
  missing_summary,
  file.path(OUTPUT_DIR, "missing_value_summary.csv"),
  row.names = FALSE
)

if (any(missing_summary$missing_values > 0L)) {
  stop(
    "Missing or non-numeric values detected. See missing_value_summary.csv.",
    call. = FALSE
  )
}

duplicate_ids <- fish_model |>
  dplyr::count(Fish_id, name = "number_of_records") |>
  dplyr::filter(number_of_records > 1L)

if (nrow(duplicate_ids) > 0L) {
  write.csv(
    duplicate_ids,
    file.path(OUTPUT_DIR, "duplicate_fish_ids.csv"),
    row.names = FALSE
  )
  stop("Duplicate Fish_id values detected.", call. = FALSE)
}

expected_reaches <- c("HC", "QJ", "SY")
expected_campaigns <- c("Wet_2023", "Dry_2023")
expected_guilds <- c(
  "Omnivorous", "Carnivorous", "Herbivorous", "Planktivorous"
)

unexpected_reaches <- setdiff(unique(fish_model$Reach_code), expected_reaches)
unexpected_campaigns <- setdiff(unique(fish_model$Campaign), expected_campaigns)
unexpected_guilds <- setdiff(unique(fish_model$Feeding_guild), expected_guilds)

if (length(unexpected_reaches) > 0L) {
  stop(
    paste0("Unexpected Reach_code: ", paste(unexpected_reaches, collapse = ", ")),
    call. = FALSE
  )
}

if (length(unexpected_campaigns) > 0L) {
  stop(
    paste0("Unexpected Campaign: ", paste(unexpected_campaigns, collapse = ", ")),
    call. = FALSE
  )
}

if (length(unexpected_guilds) > 0L) {
  stop(
    paste0("Unexpected Feeding_guild: ", paste(unexpected_guilds, collapse = ", ")),
    call. = FALSE
  )
}

guild_consistency <- fish_model |>
  dplyr::distinct(Species_code, Feeding_guild) |>
  dplyr::count(Species_code, name = "number_of_guilds") |>
  dplyr::filter(number_of_guilds > 1L)

if (nrow(guild_consistency) > 0L) {
  write.csv(
    guild_consistency,
    file.path(OUTPUT_DIR, "species_with_multiple_feeding_guilds.csv"),
    row.names = FALSE
  )
  stop("At least one species has multiple feeding-guild classifications.",
       call. = FALSE)
}

if (nrow(fish_model) != 609L) {
  stop(
    paste0("Expected 609 observations, but found ", nrow(fish_model), "."),
    call. = FALSE
  )
}

# ---- Factors and standardized responses -------------------------------------

fish_model <- fish_model |>
  dplyr::mutate(
    Reach_code = factor(Reach_code, levels = expected_reaches),
    Campaign = factor(Campaign, levels = expected_campaigns),
    Species_code = factor(Species_code),
    Feeding_guild = factor(Feeding_guild, levels = expected_guilds)
  )

factor_columns <- c(
  "Reach_code", "Campaign", "Species_code", "Feeding_guild"
)

factor_na_check <- colSums(is.na(fish_model[factor_columns]))

if (any(factor_na_check > 0L)) {
  stop(
    paste0(
      "Factor conversion created missing values:\n",
      paste(names(factor_na_check), factor_na_check,
            sep = " = ", collapse = "\n")
    ),
    call. = FALSE
  )
}

d13C_mean <- mean(fish_model$d13C_permil)
d13C_sd <- stats::sd(fish_model$d13C_permil)
d15N_mean <- mean(fish_model$d15N_permil)
d15N_sd <- stats::sd(fish_model$d15N_permil)

if (!is.finite(d13C_sd) || d13C_sd <= 0 ||
    !is.finite(d15N_sd) || d15N_sd <= 0) {
  stop("The isotope responses cannot be standardized.", call. = FALSE)
}

fish_model <- fish_model |>
  dplyr::mutate(
    d13C_z = (d13C_permil - d13C_mean) / d13C_sd,
    d15N_z = (d15N_permil - d15N_mean) / d15N_sd
  )

model_variables <- c(
  "d13C_z", "d15N_z", "Reach_code", "Campaign",
  "Species_code", "Feeding_guild"
)

if (!all(stats::complete.cases(fish_model[model_variables]))) {
  stop("Some model observations are incomplete.", call. = FALSE)
}

write.csv(
  fish_model,
  file.path(OUTPUT_DIR, "model_data_used.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    response = c("d13C_permil", "d15N_permil"),
    mean = c(d13C_mean, d15N_mean),
    standard_deviation = c(d13C_sd, d15N_sd)
  ),
  file.path(OUTPUT_DIR, "response_scaling_parameters.csv"),
  row.names = FALSE
)

cat(
  "\nDATA CHECK: PASSED",
  "\nFish observations: ", nrow(fish_model),
  "\nSpecies: ", nlevels(fish_model$Species_code),
  "\n"
)

# ---- Observed-data figure ----------------------------------------------------

fish_plot <- fish_model |>
  dplyr::mutate(
    Reach_label = dplyr::recode(
      as.character(Reach_code),
      HC = "S1", QJ = "S2", SY = "S3"
    ),
    Campaign_label = dplyr::recode(
      as.character(Campaign),
      Wet_2023 = "Wet campaign", Dry_2023 = "Dry campaign"
    )
  )

observed_plot <- ggplot2::ggplot(
  fish_plot,
  ggplot2::aes(
    x = d13C_permil, y = d15N_permil,
    shape = Campaign_label, colour = Campaign_label
  )
) +
  ggplot2::geom_point(alpha = 0.70, size = 2) +
  ggplot2::facet_wrap(~ Reach_label, nrow = 1) +
  ggplot2::labs(
    title = "Fish stable-isotope distributions",
    subtitle = "Wet- and dry-campaign observations in the lower Han River",
    x = expression(delta^13*C~("\u2030")),
    y = expression(delta^15*N~("\u2030")),
    shape = "Campaign",
    colour = "Campaign"
  ) +
  ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank()
  )

print(observed_plot)

hr_ggsave(companions = TRUE, 
  file.path(OUTPUT_DIR, "Figure_1_observed_isotopes.png"),
  observed_plot,
  width = 10, height = 5, units = "in", dpi = 400
)

# ---- Bayesian multivariate hierarchical model -------------------------------

formula_d13C <- brms::bf(
  d13C_z ~ Reach_code * Campaign + Feeding_guild +
    (1 + Campaign | p | Species_code),
  family = stats::gaussian()
)

formula_d15N <- brms::bf(
  d15N_z ~ Reach_code * Campaign + Feeding_guild +
    (1 + Campaign | p | Species_code),
  family = stats::gaussian()
)

model_formula <- formula_d13C + formula_d15N + brms::set_rescor(TRUE)

# ---- Response names and proper priors ----------------------------------------

prior_structure <- brms::get_prior(
  formula = model_formula,
  data = fish_model
)

write.csv(
  prior_structure,
  file.path(OUTPUT_DIR, "accepted_prior_structure.csv"),
  row.names = FALSE
)

response_ids <- unique(
  prior_structure$resp[
    !is.na(prior_structure$resp) & nzchar(prior_structure$resp)
  ]
)

resp_d13C <- response_ids[
  grepl("13", response_ids, ignore.case = TRUE)
][1]

resp_d15N <- response_ids[
  grepl("15", response_ids, ignore.case = TRUE)
][1]

if (length(resp_d13C) == 0L || length(resp_d15N) == 0L ||
    is.na(resp_d13C) || is.na(resp_d15N)) {
  stop(
    paste0(
      "Internal response names not identified. Detected: ",
      paste(response_ids, collapse = ", ")
    ),
    call. = FALSE
  )
}

model_priors <- c(
  brms::set_prior(
    "normal(0, 0.75)", class = "b", resp = resp_d13C
  ),
  brms::set_prior(
    "normal(0, 1)", class = "Intercept", resp = resp_d13C
  ),
  brms::set_prior(
    "normal(0, 0.75)", class = "b", resp = resp_d15N
  ),
  brms::set_prior(
    "normal(0, 1)", class = "Intercept", resp = resp_d15N
  )
)

validated_priors <- brms::validate_prior(
  prior = model_priors,
  formula = model_formula,
  data = fish_model
)

write.csv(
  as.data.frame(validated_priors),
  file.path(OUTPUT_DIR, "validated_priors.csv"),
  row.names = FALSE
)

cat(
  "\nPRIOR VALIDATION: PASSED",
  "\nd13C response: ", resp_d13C,
  "\nd15N response: ", resp_d15N,
  "\n"
)

# ---- Fit and cache the model -------------------------------------------------

model_cache_base <- file.path(
  OUTPUT_DIR,
  "Han_River_Bayesian_MV_cache"
)

model_file <- file.path(
  OUTPUT_DIR,
  "Han_River_Bayesian_MV_final.rds"
)

cat(
  "\n============================================================\n",
  "MODEL COMPILATION AND SAMPLING ARE STARTING\n",
  "Do not close RStudio while the four chains are running.\n",
  "============================================================\n",
  sep = ""
)

fit_mv <- brms::brm(
  formula = model_formula,
  data = fish_model,
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
  save_pars = brms::save_pars(all = TRUE),
  file = model_cache_base,
  file_refit = if (isTRUE(getOption("hanriver.force_refit", FALSE))) "always" else "on_change",
  refresh = 100
)

if (!inherits(fit_mv, "brmsfit")) {
  stop("The fitted object is not a valid brmsfit model.", call. = FALSE)
}

saveRDS(fit_mv, model_file)

if (!file.exists(model_file)) {
  stop("The fitted model could not be saved.", call. = FALSE)
}

model_summary <- summary(fit_mv)
print(model_summary)

writeLines(
  capture.output(model_summary),
  file.path(OUTPUT_DIR, "model_summary.txt")
)

# ---- Fixed effects -----------------------------------------------------------

fixed_effect_matrix <- brms::fixef(
  fit_mv,
  summary = TRUE,
  probs = c(0.025, 0.975)
)

fixed_effects <- data.frame(
  parameter = rownames(fixed_effect_matrix),
  estimate = fixed_effect_matrix[, "Estimate"],
  posterior_sd = fixed_effect_matrix[, "Est.Error"],
  lower_95 = fixed_effect_matrix[, "Q2.5"],
  upper_95 = fixed_effect_matrix[, "Q97.5"],
  row.names = NULL,
  stringsAsFactors = FALSE
)

write.csv(
  fixed_effects,
  file.path(OUTPUT_DIR, "fixed_effects_standardized.csv"),
  row.names = FALSE
)

# ---- Convergence and sampler diagnostics ------------------------------------

posterior_draws <- posterior::as_draws_df(fit_mv)
parameter_diagnostics <- posterior::summarise_draws(posterior_draws)

write.csv(
  parameter_diagnostics,
  file.path(OUTPUT_DIR, "parameter_diagnostics.csv"),
  row.names = FALSE
)

finite_rhat <- parameter_diagnostics$rhat[
  is.finite(parameter_diagnostics$rhat)
]
finite_bulk_ess <- parameter_diagnostics$ess_bulk[
  is.finite(parameter_diagnostics$ess_bulk)
]
finite_tail_ess <- parameter_diagnostics$ess_tail[
  is.finite(parameter_diagnostics$ess_tail)
]

nuts_diagnostics <- brms::nuts_params(fit_mv)

divergent_transitions <- sum(
  nuts_diagnostics$Parameter == "divergent__" &
    nuts_diagnostics$Value == 1
)

treedepth_hits <- sum(
  nuts_diagnostics$Parameter == "treedepth__" &
    nuts_diagnostics$Value >= MAX_TREEDEPTH
)

diagnostic_overview <- data.frame(
  diagnostic = c(
    "Number of observations",
    "Number of species",
    "Maximum Rhat",
    "Minimum bulk ESS",
    "Minimum tail ESS",
    "Divergent transitions",
    "Maximum-treedepth hits"
  ),
  value = c(
    nrow(fish_model),
    nlevels(fish_model$Species_code),
    if (length(finite_rhat) > 0L) max(finite_rhat) else NA_real_,
    if (length(finite_bulk_ess) > 0L) min(finite_bulk_ess) else NA_real_,
    if (length(finite_tail_ess) > 0L) min(finite_tail_ess) else NA_real_,
    divergent_transitions,
    treedepth_hits
  )
)

print(diagnostic_overview)

write.csv(
  diagnostic_overview,
  file.path(OUTPUT_DIR, "diagnostic_overview.csv"),
  row.names = FALSE
)

# ---- Residual correlation and Bayesian R-squared ----------------------------

all_parameter_summary <- as.data.frame(
  brms::posterior_summary(
    fit_mv,
    probs = c(0.025, 0.975)
  )
)

all_parameter_summary$parameter <- rownames(all_parameter_summary)
rownames(all_parameter_summary) <- NULL

residual_correlation <- all_parameter_summary[
  grepl("rescor", all_parameter_summary$parameter, ignore.case = TRUE),
  ,
  drop = FALSE
]

write.csv(
  residual_correlation,
  file.path(OUTPUT_DIR, "residual_correlation.csv"),
  row.names = FALSE
)

r2_d13C <- as.data.frame(
  brms::bayes_R2(fit_mv, resp = resp_d13C, summary = TRUE)
)
r2_d13C$response <- "d13C"

r2_d15N <- as.data.frame(
  brms::bayes_R2(fit_mv, resp = resp_d15N, summary = TRUE)
)
r2_d15N$response <- "d15N"

bayesian_r2 <- dplyr::bind_rows(r2_d13C, r2_d15N)

write.csv(
  bayesian_r2,
  file.path(OUTPUT_DIR, "Bayesian_R2.csv"),
  row.names = FALSE
)

# ---- Fixed-effect forest plot ------------------------------------------------

fixed_effects$response <- ifelse(
  grepl(paste0("^", resp_d13C, "_"), fixed_effects$parameter),
  "δ13C",
  ifelse(
    grepl(paste0("^", resp_d15N, "_"), fixed_effects$parameter),
    "δ15N",
    NA_character_
  )
)

fixed_effects$term <- fixed_effects$parameter
fixed_effects$term <- sub(
  paste0("^", resp_d13C, "_"), "", fixed_effects$term
)
fixed_effects$term <- sub(
  paste0("^", resp_d15N, "_"), "", fixed_effects$term
)

fixed_effects$term <- gsub(
  "Reach_codeQJ", "S2 vs S1",
  fixed_effects$term, fixed = TRUE
)
fixed_effects$term <- gsub(
  "Reach_codeSY", "S3 vs S1",
  fixed_effects$term, fixed = TRUE
)
fixed_effects$term <- gsub(
  "CampaignDry_2023", "Dry vs wet campaign",
  fixed_effects$term, fixed = TRUE
)
fixed_effects$term <- gsub(
  "Feeding_guildCarnivorous", "Carnivorous vs omnivorous",
  fixed_effects$term, fixed = TRUE
)
fixed_effects$term <- gsub(
  "Feeding_guildHerbivorous", "Herbivorous vs omnivorous",
  fixed_effects$term, fixed = TRUE
)
fixed_effects$term <- gsub(
  "Feeding_guildPlanktivorous", "Planktivorous vs omnivorous",
  fixed_effects$term, fixed = TRUE
)
fixed_effects$term <- gsub(
  ":", " × ", fixed_effects$term, fixed = TRUE
)

fixed_effect_plot_data <- fixed_effects[
  !grepl("Intercept", fixed_effects$parameter, ignore.case = TRUE) &
    !is.na(fixed_effects$response),
  ,
  drop = FALSE
]

if (nrow(fixed_effect_plot_data) > 0L) {
  fixed_effect_plot_data$term <- factor(
    fixed_effect_plot_data$term,
    levels = rev(unique(fixed_effect_plot_data$term))
  )

  forest_plot <- ggplot2::ggplot(
    fixed_effect_plot_data,
    ggplot2::aes(x = estimate, y = term)
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = lower_95, xend = upper_95,
        y = term, yend = term
      ),
      linewidth = 0.8
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::facet_wrap(
      ~ response,
      ncol = 1,
      scales = "free_y"
    ) +
    ggplot2::labs(
      title = "Posterior population-level effects",
      subtitle = "Points are posterior means; lines are 95% credible intervals",
      x = "Standardized posterior effect",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    )

  print(forest_plot)

  hr_ggsave(companions = TRUE, 
    file.path(OUTPUT_DIR, "Figure_2_fixed_effect_forest.png"),
    forest_plot,
    width = 10, height = 8, units = "in", dpi = 400
  )
}

# ---- Posterior predictive checks ---------------------------------------------

save_plot <- function(plot_object, filename, width = 8, height = 5) {
  print(plot_object)
  hr_ggsave(companions = TRUE, 
    file.path(OUTPUT_DIR, filename),
    plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 400
  )
}

pp_density_d13C <- brms::pp_check(
  fit_mv, resp = resp_d13C,
  type = "dens_overlay", ndraws = 50
) +
  ggplot2::labs(
    title = "Posterior predictive check for δ13C",
    x = "Standardized δ13C", y = "Density"
  ) +
  ggplot2::theme_bw(base_size = 14)

pp_density_d15N <- brms::pp_check(
  fit_mv, resp = resp_d15N,
  type = "dens_overlay", ndraws = 50
) +
  ggplot2::labs(
    title = "Posterior predictive check for δ15N",
    x = "Standardized δ15N", y = "Density"
  ) +
  ggplot2::theme_bw(base_size = 14)

pp_mean_d13C <- brms::pp_check(
  fit_mv, resp = resp_d13C,
  type = "stat", stat = "mean", ndraws = 500
) +
  ggplot2::labs(title = "Posterior predictive check: mean δ13C") +
  ggplot2::theme_bw(base_size = 14)

pp_mean_d15N <- brms::pp_check(
  fit_mv, resp = resp_d15N,
  type = "stat", stat = "mean", ndraws = 500
) +
  ggplot2::labs(title = "Posterior predictive check: mean δ15N") +
  ggplot2::theme_bw(base_size = 14)

pp_sd_d13C <- brms::pp_check(
  fit_mv, resp = resp_d13C,
  type = "stat", stat = "sd", ndraws = 500
) +
  ggplot2::labs(title = "Posterior predictive check: SD of δ13C") +
  ggplot2::theme_bw(base_size = 14)

pp_sd_d15N <- brms::pp_check(
  fit_mv, resp = resp_d15N,
  type = "stat", stat = "sd", ndraws = 500
) +
  ggplot2::labs(title = "Posterior predictive check: SD of δ15N") +
  ggplot2::theme_bw(base_size = 14)

save_plot(pp_density_d13C, "Figure_3_PPC_density_d13C.png")
save_plot(pp_density_d15N, "Figure_4_PPC_density_d15N.png")
save_plot(pp_mean_d13C, "Figure_5_PPC_mean_d13C.png")
save_plot(pp_mean_d15N, "Figure_6_PPC_mean_d15N.png")
save_plot(pp_sd_d13C, "Figure_7_PPC_SD_d13C.png")
save_plot(pp_sd_d15N, "Figure_8_PPC_SD_d15N.png")

# ---- Save session information and completion note ----------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "session_information.txt")
)

completion_note <- c(
  "HAN RIVER BAYESIAN MULTIVARIATE MODEL COMPLETED",
  "",
  paste0("Data file: ", DATA_FILE),
  paste0("Model file: ", model_file),
  paste0("Output folder: ", OUTPUT_DIR),
  "",
  "CHECK BEFORE INTERPRETATION:",
  "Maximum Rhat should be below 1.01.",
  "Divergent transitions should equal zero.",
  "Maximum-treedepth hits should equal zero.",
  "Inspect all posterior predictive checks.",
  "",
  "MAIN OUTPUTS:",
  "model_summary.txt",
  "fixed_effects_standardized.csv",
  "diagnostic_overview.csv",
  "residual_correlation.csv",
  "Bayesian_R2.csv"
)

writeLines(
  completion_note,
  file.path(OUTPUT_DIR, "READ_ME_FIRST.txt")
)

cat(
  "\n============================================================\n",
  "HAN RIVER BAYESIAN MULTIVARIATE MODEL COMPLETED\n",
  "============================================================\n",
  "fit_mv exists: ", exists("fit_mv"),
  "\nValid brmsfit object: ", inherits(fit_mv, "brmsfit"),
  "\nSaved model exists: ", file.exists(model_file),
  "\nResults folder:\n", OUTPUT_DIR,
  "\n\nOpen READ_ME_FIRST.txt first.\n",
  "============================================================\n",
  sep = ""
)

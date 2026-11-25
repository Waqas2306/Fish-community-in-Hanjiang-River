# RUN IN RSTUDIO: open this R file and click Source.
# Input file(s): No direct Excel input
# Worksheet tabs read automatically: not applicable
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
# HAN RIVER BAYESIAN MULTIVARIATE MODEL
# PUBLICATION-QUALITY VISUALIZATION SCRIPT
#
# PURPOSE
#   Creates professional figures from the already fitted Bayesian multivariate
#   hierarchical model. This script DOES NOT refit the model.
#
# REQUIRED FILES
#   Results/<current-input-signature>/main_model/
#       Han_River_Bayesian_MV_final.rds
#       model_data_used.csv
#
# OUTPUT
#   Results/<current-input-signature>/main_model/
#       Figures_publication/
#
# HOW TO RUN
#   1. Create a new blank R script.
#   2. Paste this entire script.
#   3. Save it as Han_River_publication_figures.R.
#   4. Click Source once.
# =============================================================================


# =============================================================================
# 1. SETTINGS
# =============================================================================

graphics.off()

RESULTS_DIR <- hr_results("main_model")

FIGURE_DIR <- file.path(
  RESULTS_DIR,
  "Figures_publication"
)

MODEL_FILE <- file.path(
  RESULTS_DIR,
  "Han_River_Bayesian_MV_final.rds"
)

DATA_FILE <- file.path(
  RESULTS_DIR,
  "model_data_used.csv"
)

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

SEED <- 20260729L

# Journal figures normally place titles in the manuscript caption.
# Change to TRUE when you want titles shown inside the figures.
SHOW_TITLES <- FALSE

# "sans" is portable and avoids missing-font errors.
FONT_FAMILY <- "sans"

set.seed(SEED)


# =============================================================================
# 2. CHECK AND LOAD PACKAGES
# =============================================================================

required_packages <- c(
  "brms",
  "posterior",
  "ggplot2",
  "dplyr",
  "tidyr"
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
      "\nInstall them, restart RStudio, and run this script again."
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(brms)
  library(posterior)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})


# =============================================================================
# 3. LOAD THE SAVED MODEL AND MODEL DATA
# =============================================================================

hr_ensure_main_model()

if (!file.exists(MODEL_FILE)) {
  stop(
    paste0(
      "The fitted model was not found:\n",
      MODEL_FILE
    ),
    call. = FALSE
  )
}

if (!file.exists(DATA_FILE)) {
  stop(
    paste0(
      "The model-data file was not found:\n",
      DATA_FILE
    ),
    call. = FALSE
  )
}

fit_mv <- readRDS(
  MODEL_FILE
)

if (!inherits(fit_mv, "brmsfit")) {
  stop(
    "The saved object is not a valid brmsfit model.",
    call. = FALSE
  )
}

fish_model <- read.csv(
  DATA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# =============================================================================
# 4. VALIDATE AND RESTORE THE MODEL DATA
# =============================================================================

required_columns <- c(
  "Fish_id",
  "Reach_code",
  "Campaign",
  "Species_code",
  "Feeding_guild",
  "d13C_permil",
  "d15N_permil",
  "d13C_z",
  "d15N_z"
)

missing_columns <- setdiff(
  required_columns,
  names(fish_model)
)

if (length(missing_columns) > 0L) {
  stop(
    paste0(
      "The model data are missing these column(s): ",
      paste(missing_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

fish_model <- fish_model |>
  dplyr::mutate(
    Reach_code = factor(
      as.character(Reach_code),
      levels = c("HC", "QJ", "SY")
    ),

    Campaign = factor(
      as.character(Campaign),
      levels = c("Wet_2023", "Dry_2023")
    ),

    Species_code = factor(
      as.character(Species_code)
    ),

    Feeding_guild = factor(
      as.character(Feeding_guild),
      levels = c(
        "Omnivorous",
        "Carnivorous",
        "Herbivorous",
        "Planktivorous"
      )
    )
  )

if (anyNA(fish_model[required_columns])) {
  stop(
    "Missing values were detected in required plotting variables.",
    call. = FALSE
  )
}

if (nrow(fish_model) != 609L) {
  warning(
    paste0(
      "The plotting data contain ",
      nrow(fish_model),
      " rows rather than the expected 609 rows."
    ),
    call. = FALSE
  )
}


# =============================================================================
# 5. EXTRACT POSTERIOR DRAWS AND IDENTIFY RESPONSE NAMES
# =============================================================================

posterior_draws <- posterior::as_draws_df(
  fit_mv
)

fixed_parameter_names_all <- grep(
  "^b_",
  names(posterior_draws),
  value = TRUE
)

response_candidates <- unique(
  sub(
    "^b_([^_]+)_.*$",
    "\\1",
    fixed_parameter_names_all
  )
)

resp_d13C <- response_candidates[
  grepl(
    "13",
    response_candidates,
    ignore.case = TRUE
  )
][1]

resp_d15N <- response_candidates[
  grepl(
    "15",
    response_candidates,
    ignore.case = TRUE
  )
][1]

if (
  length(resp_d13C) == 0L ||
  length(resp_d15N) == 0L ||
  is.na(resp_d13C) ||
  is.na(resp_d15N)
) {
  stop(
    paste0(
      "The internal response names could not be identified.\n",
      "Detected response candidates: ",
      paste(response_candidates, collapse = ", ")
    ),
    call. = FALSE
  )
}

cat(
  "\nModel and data loaded successfully.",
  "\nInternal response names:",
  "\n  δ13C = ", resp_d13C,
  "\n  δ15N = ", resp_d15N,
  "\n",
  sep = ""
)


# =============================================================================
# 6. COLOUR-BLIND-SAFE PALETTE
# =============================================================================

COL_DRY <- "#D55E00"
COL_WET <- "#0072B2"

COL_OBSERVED <- "#162A43"
COL_PREDICTED <- "#9ECAE1"

COL_INTERVAL_95 <- "#9A9A9A"
COL_INTERVAL_50 <- "#222222"
COL_REFERENCE <- "#666666"

COL_CARBON <- "#0072B2"
COL_NITROGEN <- "#D55E00"


# =============================================================================
# 7. PUBLICATION THEME
# =============================================================================

theme_publication <- function(base_size = 14) {

  ggplot2::theme_classic(
    base_size = base_size,
    base_family = FONT_FAMILY
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 2,
        hjust = 0
      ),

      plot.subtitle = ggplot2::element_text(
        size = base_size,
        hjust = 0,
        margin = ggplot2::margin(b = 8)
      ),

      plot.caption = ggplot2::element_text(
        size = base_size - 1,
        colour = "grey30",
        hjust = 0,
        margin = ggplot2::margin(t = 7)
      ),

      axis.title = ggplot2::element_text(
        size = base_size + 1
      ),

      axis.text = ggplot2::element_text(
        size = base_size,
        colour = "black"
      ),

      axis.line = ggplot2::element_line(
        linewidth = 0.45,
        colour = "black"
      ),

      axis.ticks = ggplot2::element_line(
        linewidth = 0.4,
        colour = "black"
      ),

      strip.background = ggplot2::element_rect(
        fill = "grey94",
        colour = "black",
        linewidth = 0.45
      ),

      strip.text = ggplot2::element_text(
        face = "bold",
        size = base_size,
        colour = "black",
        margin = ggplot2::margin(
          t = 5,
          r = 5,
          b = 5,
          l = 5
        )
      ),

      legend.position = "top",
      legend.justification = "left",

      legend.title = ggplot2::element_text(
        face = "bold"
      ),

      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),

      panel.spacing = grid::unit(
        7,
        "pt"
      ),

      plot.margin = ggplot2::margin(
        t = 8,
        r = 10,
        b = 8,
        l = 8
      )
    ) + hr_typography()
}


# =============================================================================
# 8. OPTIONAL TITLES
# =============================================================================

add_optional_title <- function(
  plot_object,
  title_text,
  subtitle_text = NULL
) {

  if (isTRUE(SHOW_TITLES)) {
    plot_object <- plot_object +
      ggplot2::labs(
        title = title_text,
        subtitle = subtitle_text
      )
  }

  plot_object
}


# =============================================================================
# 9. SAVE EACH FIGURE AS PNG, TIFF AND PDF
# =============================================================================

save_publication_figure <- function(
  plot_object,
  filename,
  width,
  height
) {

  png_file <- file.path(
    FIGURE_DIR,
    paste0(filename, ".png")
  )

  tiff_file <- file.path(
    FIGURE_DIR,
    paste0(filename, ".tiff")
  )

  pdf_file <- file.path(
    FIGURE_DIR,
    paste0(filename, ".pdf")
  )

  hr_ggsave(
    filename = png_file,
    plot = plot_object,
    device = "png",
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  hr_ggsave(
    filename = tiff_file,
    plot = plot_object,
    device = "tiff",
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
        filename = pdf_file,
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
      filename = pdf_file,
      plot = plot_object,
      device = "pdf",
      width = width,
      height = height,
      units = "in",
      bg = "white"
    )
  }

  message(
    "Saved: ",
    filename
  )
}


# =============================================================================
# 10. HELPER: CONVERT BRMS PREDICTIONS TO A MATRIX
# =============================================================================

as_prediction_matrix <- function(prediction_object) {

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
      "Unexpected prediction dimensions: ",
      paste(object_dimensions, collapse = " × ")
    ),
    call. = FALSE
  )
}


# =============================================================================
# FIGURE 1. OBSERVED δ13C–δ15N DISTRIBUTIONS
# =============================================================================

observed_plot_data <- fish_model |>
  dplyr::mutate(
    Reach = dplyr::recode(
      as.character(Reach_code),
      HC = "S1",
      QJ = "S2",
      SY = "S3"
    ),

    Campaign_label = dplyr::recode(
      as.character(Campaign),
      Wet_2023 = "Wet campaign",
      Dry_2023 = "Dry campaign"
    )
  )

figure_1 <- ggplot2::ggplot(
  observed_plot_data,
  ggplot2::aes(
    x = d13C_permil,
    y = d15N_permil,
    fill = Campaign_label,
    shape = Campaign_label
  )
) +
  ggplot2::geom_point(
    size = 2.05,
    stroke = 0.35,
    colour = "black",
    alpha = 0.72
  ) +
  ggplot2::facet_wrap(
    ~ Reach,
    nrow = 1
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Dry campaign" = COL_DRY,
      "Wet campaign" = COL_WET
    ),
    breaks = c(
      "Dry campaign",
      "Wet campaign"
    )
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "Dry campaign" = 24,
      "Wet campaign" = 21
    ),
    breaks = c(
      "Dry campaign",
      "Wet campaign"
    )
  ) +
  ggplot2::labs(
    x = expression(delta^13*C~("\u2030")),
    y = expression(delta^15*N~("\u2030")),
    fill = "Campaign",
    shape = "Campaign"
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      override.aes = list(
        size = 3,
        alpha = 1
      )
    ),
    shape = ggplot2::guide_legend(
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  ) +
  theme_publication(
    base_size = 14
  )

figure_1 <- add_optional_title(
  figure_1,
  title_text = "Fish stable-isotope distributions",
  subtitle_text = paste(
    "Wet- and dry-campaign observations",
    "in the lower Han River"
  )
)

print(
  figure_1
)

save_publication_figure(
  plot_object = figure_1,
  filename = "Figure_1_observed_isotope_distributions",
  width = 7.2,
  height = 4.4
)


# =============================================================================
# FIGURE 2. POSTERIOR FIXED-EFFECT FOREST PLOT
# =============================================================================

fixed_parameter_names <- fixed_parameter_names_all[
  !grepl(
    "_Intercept$",
    fixed_parameter_names_all
  )
]

clean_effect_label <- function(x) {

  x <- gsub(
    "Reach_codeQJ",
    "S2 vs S1",
    x,
    fixed = TRUE
  )

  x <- gsub(
    "Reach_codeSY",
    "S3 vs S1",
    x,
    fixed = TRUE
  )

  x <- gsub(
    "CampaignDry_2023",
    "Dry vs wet campaign",
    x,
    fixed = TRUE
  )

  x <- gsub(
    "Feeding_guildCarnivorous",
    "Carnivorous vs omnivorous",
    x,
    fixed = TRUE
  )

  x <- gsub(
    "Feeding_guildHerbivorous",
    "Herbivorous vs omnivorous",
    x,
    fixed = TRUE
  )

  x <- gsub(
    "Feeding_guildPlanktivorous",
    "Planktivorous vs omnivorous",
    x,
    fixed = TRUE
  )

  x <- gsub(
    ":",
    " × ",
    x,
    fixed = TRUE
  )

  x
}

summarise_fixed_parameter <- function(parameter_name) {

  values <- as.numeric(
    posterior_draws[[parameter_name]]
  )

  if (
    grepl(
      paste0("^b_", resp_d13C, "_"),
      parameter_name
    )
  ) {
    response_name <- "δ13C"
  } else if (
    grepl(
      paste0("^b_", resp_d15N, "_"),
      parameter_name
    )
  ) {
    response_name <- "δ15N"
  } else {
    response_name <- NA_character_
  }

  raw_term <- parameter_name

  raw_term <- sub(
    paste0("^b_", resp_d13C, "_"),
    "",
    raw_term
  )

  raw_term <- sub(
    paste0("^b_", resp_d15N, "_"),
    "",
    raw_term
  )

  data.frame(
    parameter = parameter_name,

    response = response_name,

    term = clean_effect_label(
      raw_term
    ),

    estimate = stats::median(
      values,
      na.rm = TRUE
    ),

    lower_95 = unname(
      stats::quantile(
        values,
        probs = 0.025,
        na.rm = TRUE
      )
    ),

    upper_95 = unname(
      stats::quantile(
        values,
        probs = 0.975,
        na.rm = TRUE
      )
    ),

    lower_50 = unname(
      stats::quantile(
        values,
        probs = 0.25,
        na.rm = TRUE
      )
    ),

    upper_50 = unname(
      stats::quantile(
        values,
        probs = 0.75,
        na.rm = TRUE
      )
    ),

    probability_direction = max(
      mean(
        values > 0,
        na.rm = TRUE
      ),
      mean(
        values < 0,
        na.rm = TRUE
      )
    ),

    stringsAsFactors = FALSE
  )
}

fixed_plot_data <- dplyr::bind_rows(
  lapply(
    fixed_parameter_names,
    summarise_fixed_parameter
  )
)

fixed_plot_data <- fixed_plot_data |>
  dplyr::filter(
    !is.na(response)
  ) |>
  dplyr::mutate(
    credible_95 = (
      lower_95 > 0 |
      upper_95 < 0
    )
  )

term_levels <- rev(
  unique(
    fixed_plot_data$term
  )
)

fixed_plot_data$term <- factor(
  fixed_plot_data$term,
  levels = term_levels
)

figure_2 <- ggplot2::ggplot(
  fixed_plot_data,
  ggplot2::aes(
    x = estimate,
    y = term
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = COL_REFERENCE
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = lower_95,
      xend = upper_95,
      y = term,
      yend = term
    ),
    linewidth = 0.65,
    colour = COL_INTERVAL_95
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = lower_50,
      xend = upper_50,
      y = term,
      yend = term
    ),
    linewidth = 2.1,
    colour = COL_INTERVAL_50,
    lineend = "round"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      fill = credible_95
    ),
    shape = 21,
    size = 2.8,
    stroke = 0.55,
    colour = "black"
  ) +
  ggplot2::facet_wrap(
    ~ response,
    ncol = 1,
    scales = "free_y"
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "TRUE" = "#222222",
      "FALSE" = "white"
    ),
    labels = c(
      "FALSE" = "95% CrI includes zero",
      "TRUE" = "95% CrI excludes zero"
    ),
    name = NULL
  ) +
  ggplot2::labs(
    x = "Standardized posterior effect",
    y = NULL,
    caption = paste(
      "Thin lines show 95% credible intervals;",
      "thick lines show 50% credible intervals."
    )
  ) +
  theme_publication(
    base_size = 14
  ) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_line(
      colour = "grey90",
      linewidth = 0.3
    ),
    axis.line.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )

figure_2 <- add_optional_title(
  figure_2,
  title_text = "Posterior population-level effects",
  subtitle_text = paste(
    "Posterior medians with 50% and",
    "95% credible intervals"
  )
)

print(
  figure_2
)

save_publication_figure(
  plot_object = figure_2,
  filename = "Figure_2_posterior_fixed_effects",
  width = 7.2,
  height = 6.5
)

write.csv(
  fixed_plot_data,
  file.path(
    FIGURE_DIR,
    "Figure_2_fixed_effect_plot_data.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# FIGURE 3. ADJUSTED REACH × CAMPAIGN ESTIMATES
# =============================================================================

guild_weights <- prop.table(
  table(
    fish_model$Feeding_guild
  )
)

prediction_grid <- expand.grid(
  Reach_code = levels(
    fish_model$Reach_code
  ),

  Campaign = levels(
    fish_model$Campaign
  ),

  Feeding_guild = levels(
    fish_model$Feeding_guild
  ),

  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

prediction_grid <- prediction_grid |>
  dplyr::mutate(
    Reach_code = factor(
      Reach_code,
      levels = levels(
        fish_model$Reach_code
      )
    ),

    Campaign = factor(
      Campaign,
      levels = levels(
        fish_model$Campaign
      )
    ),

    Feeding_guild = factor(
      Feeding_guild,
      levels = levels(
        fish_model$Feeding_guild
      )
    ),

    Species_code = factor(
      levels(
        fish_model$Species_code
      )[1],
      levels = levels(
        fish_model$Species_code
      )
    )
  )

summarise_adjusted_predictions <- function(
  response_id,
  raw_mean,
  raw_sd,
  response_label
) {

  predictions <- brms::posterior_epred(
    fit_mv,
    newdata = prediction_grid,
    resp = response_id,
    re_formula = NA,
    ndraws = 2000
  )

  predictions <- as_prediction_matrix(
    predictions
  )

  combinations <- prediction_grid |>
    dplyr::distinct(
      Reach_code,
      Campaign
    )

  summaries <- vector(
    mode = "list",
    length = nrow(combinations)
  )

  for (
    combination_index in seq_len(
      nrow(combinations)
    )
  ) {

    current_reach <- combinations$Reach_code[
      combination_index
    ]

    current_campaign <- combinations$Campaign[
      combination_index
    ]

    column_indices <- which(
      prediction_grid$Reach_code == current_reach &
      prediction_grid$Campaign == current_campaign
    )

    current_guilds <- as.character(
      prediction_grid$Feeding_guild[
        column_indices
      ]
    )

    current_weights <- as.numeric(
      guild_weights[
        current_guilds
      ]
    )

    if (
      anyNA(current_weights) ||
      sum(current_weights) <= 0
    ) {
      stop(
        "Feeding-guild weights could not be calculated.",
        call. = FALSE
      )
    }

    current_weights <- current_weights /
      sum(current_weights)

    averaged_draws <- as.numeric(
      predictions[
        ,
        column_indices,
        drop = FALSE
      ] %*%
      current_weights
    )

    raw_draws <- (
      averaged_draws *
      raw_sd
    ) +
      raw_mean

    summaries[[combination_index]] <- data.frame(
      Reach_code = as.character(
        current_reach
      ),

      Campaign = as.character(
        current_campaign
      ),

      response = response_label,

      estimate = stats::median(
        raw_draws,
        na.rm = TRUE
      ),

      lower_95 = unname(
        stats::quantile(
          raw_draws,
          probs = 0.025,
          na.rm = TRUE
        )
      ),

      upper_95 = unname(
        stats::quantile(
          raw_draws,
          probs = 0.975,
          na.rm = TRUE
        )
      ),

      lower_50 = unname(
        stats::quantile(
          raw_draws,
          probs = 0.25,
          na.rm = TRUE
        )
      ),

      upper_50 = unname(
        stats::quantile(
          raw_draws,
          probs = 0.75,
          na.rm = TRUE
        )
      ),

      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(
    summaries
  )
}

d13C_mean <- mean(
  fish_model$d13C_permil
)

d13C_sd <- stats::sd(
  fish_model$d13C_permil
)

d15N_mean <- mean(
  fish_model$d15N_permil
)

d15N_sd <- stats::sd(
  fish_model$d15N_permil
)

adjusted_d13C <- summarise_adjusted_predictions(
  response_id = resp_d13C,
  raw_mean = d13C_mean,
  raw_sd = d13C_sd,
  response_label = "δ13C"
)

adjusted_d15N <- summarise_adjusted_predictions(
  response_id = resp_d15N,
  raw_mean = d15N_mean,
  raw_sd = d15N_sd,
  response_label = "δ15N"
)

adjusted_plot_data <- dplyr::bind_rows(
  adjusted_d13C,
  adjusted_d15N
) |>
  dplyr::mutate(
    Reach = dplyr::recode(
      Reach_code,
      HC = "S1",
      QJ = "S2",
      SY = "S3"
    ),

    Campaign_label = dplyr::recode(
      Campaign,
      Wet_2023 = "Wet campaign",
      Dry_2023 = "Dry campaign"
    ),

    Reach = factor(
      Reach,
      levels = c(
        "S1",
        "S2",
        "S3"
      )
    )
  )

figure_3 <- ggplot2::ggplot(
  adjusted_plot_data,
  ggplot2::aes(
    x = Reach,
    y = estimate,
    colour = Campaign_label,
    shape = Campaign_label
  )
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = lower_95,
      ymax = upper_95
    ),
    width = 0,
    linewidth = 0.6,
    position = ggplot2::position_dodge(
      width = 0.34
    )
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = lower_50,
      ymax = upper_50
    ),
    width = 0,
    linewidth = 2.0,
    position = ggplot2::position_dodge(
      width = 0.34
    )
  ) +
  ggplot2::geom_point(
    size = 3,
    stroke = 0.7,
    position = ggplot2::position_dodge(
      width = 0.34
    )
  ) +
  ggplot2::facet_wrap(
    ~ response,
    ncol = 1,
    scales = "free_y"
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "Dry campaign" = COL_DRY,
      "Wet campaign" = COL_WET
    ),
    breaks = c(
      "Dry campaign",
      "Wet campaign"
    )
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "Dry campaign" = 17,
      "Wet campaign" = 16
    ),
    breaks = c(
      "Dry campaign",
      "Wet campaign"
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Adjusted isotope value (‰)",
    colour = "Campaign",
    shape = "Campaign",
    caption = paste(
      "Population-level posterior estimates averaged over the observed",
      "feeding-guild distribution. Thin and thick intervals represent",
      "95% and 50% credible intervals, respectively."
    )
  ) +
  theme_publication(
    base_size = 14
  ) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(
      colour = "grey90",
      linewidth = 0.3
    )
  )

figure_3 <- add_optional_title(
  figure_3,
  title_text = "Adjusted reach and campaign estimates",
  subtitle_text = paste(
    "Population-level posterior estimates",
    "in original isotope units"
  )
)

print(
  figure_3
)

save_publication_figure(
  plot_object = figure_3,
  filename = "Figure_3_adjusted_reach_campaign_estimates",
  width = 7.0,
  height = 6.2
)

write.csv(
  adjusted_plot_data,
  file.path(
    FIGURE_DIR,
    "Figure_3_adjusted_estimates.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 11. HELPER: CREATE A MANUAL DENSITY PPC
# =============================================================================

create_density_ppc <- function(
  response_id,
  observed_values,
  x_axis_label,
  title_text,
  subtitle_text,
  filename
) {

  yrep <- brms::posterior_predict(
    fit_mv,
    resp = response_id,
    ndraws = 35
  )

  yrep <- as_prediction_matrix(
    yrep
  )

  combined_range <- range(
    c(
      observed_values,
      as.numeric(yrep)
    ),
    finite = TRUE
  )

  range_padding <- diff(
    combined_range
  ) * 0.04

  if (
    !is.finite(range_padding) ||
    range_padding <= 0
  ) {
    range_padding <- 0.1
  }

  density_from <- combined_range[1] -
    range_padding

  density_to <- combined_range[2] +
    range_padding

  observed_density <- stats::density(
    observed_values,
    from = density_from,
    to = density_to,
    n = 512
  )

  observed_density_data <- data.frame(
    x = observed_density$x,
    density = observed_density$y,
    series = "Observed"
  )

  replicated_density_list <- lapply(
    seq_len(
      nrow(yrep)
    ),
    function(draw_number) {

      current_density <- stats::density(
        yrep[draw_number, ],
        from = density_from,
        to = density_to,
        n = 512
      )

      data.frame(
        x = current_density$x,
        density = current_density$y,
        draw = draw_number,
        series = "Posterior predictive"
      )
    }
  )

  replicated_density_data <- dplyr::bind_rows(
    replicated_density_list
  )

  plot_object <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = replicated_density_data,
      mapping = ggplot2::aes(
        x = x,
        y = density,
        group = draw,
        colour = series
      ),
      linewidth = 0.35,
      alpha = 0.55
    ) +
    ggplot2::geom_line(
      data = observed_density_data,
      mapping = ggplot2::aes(
        x = x,
        y = density,
        colour = series
      ),
      linewidth = 1.05
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "Observed" = COL_OBSERVED,
        "Posterior predictive" = COL_PREDICTED
      ),
      breaks = c(
        "Observed",
        "Posterior predictive"
      ),
      name = NULL
    ) +
    ggplot2::labs(
      x = x_axis_label,
      y = "Density"
    ) +
    theme_publication(
      base_size = 14
    )

  plot_object <- add_optional_title(
    plot_object,
    title_text = title_text,
    subtitle_text = subtitle_text
  )

  print(
    plot_object
  )

  save_publication_figure(
    plot_object = plot_object,
    filename = filename,
    width = 6.8,
    height = 4.8
  )

  invisible(
    yrep
  )
}


# =============================================================================
# FIGURE 4A. POSTERIOR PREDICTIVE DENSITY FOR δ13C
# =============================================================================

set.seed(
  SEED
)

density_yrep_d13C <- create_density_ppc(
  response_id = resp_d13C,
  observed_values = fish_model$d13C_z,
  x_axis_label = "Standardized δ13C",
  title_text = "Posterior predictive check for δ13C",
  subtitle_text = paste(
    "Observed distribution compared with",
    "35 posterior predictive datasets"
  ),
  filename = "Figure_4a_PPC_density_d13C"
)


# =============================================================================
# FIGURE 4B. POSTERIOR PREDICTIVE DENSITY FOR δ15N
# =============================================================================

set.seed(
  SEED + 1L
)

density_yrep_d15N <- create_density_ppc(
  response_id = resp_d15N,
  observed_values = fish_model$d15N_z,
  x_axis_label = "Standardized δ15N",
  title_text = "Posterior predictive check for δ15N",
  subtitle_text = paste(
    "Observed distribution compared with",
    "35 posterior predictive datasets"
  ),
  filename = "Figure_4b_PPC_density_d15N"
)


# =============================================================================
# FIGURE 5. COMBINED PPCs FOR MEAN AND STANDARD DEVIATION
# =============================================================================

set.seed(
  SEED + 2L
)

moment_yrep_d13C <- brms::posterior_predict(
  fit_mv,
  resp = resp_d13C,
  ndraws = 1000
)

moment_yrep_d13C <- as_prediction_matrix(
  moment_yrep_d13C
)

moment_yrep_d15N <- brms::posterior_predict(
  fit_mv,
  resp = resp_d15N,
  ndraws = 1000
)

moment_yrep_d15N <- as_prediction_matrix(
  moment_yrep_d15N
)

moment_plot_data <- dplyr::bind_rows(
  data.frame(
    response = "δ13C",
    statistic = "Mean",
    value = rowMeans(
      moment_yrep_d13C
    ),
    observed = mean(
      fish_model$d13C_z
    )
  ),

  data.frame(
    response = "δ13C",
    statistic = "Standard deviation",
    value = apply(
      moment_yrep_d13C,
      1,
      stats::sd
    ),
    observed = stats::sd(
      fish_model$d13C_z
    )
  ),

  data.frame(
    response = "δ15N",
    statistic = "Mean",
    value = rowMeans(
      moment_yrep_d15N
    ),
    observed = mean(
      fish_model$d15N_z
    )
  ),

  data.frame(
    response = "δ15N",
    statistic = "Standard deviation",
    value = apply(
      moment_yrep_d15N,
      1,
      stats::sd
    ),
    observed = stats::sd(
      fish_model$d15N_z
    )
  )
)

figure_5 <- ggplot2::ggplot(
  moment_plot_data,
  ggplot2::aes(
    x = value
  )
) +
  ggplot2::geom_histogram(
    ggplot2::aes(
      y = ggplot2::after_stat(
        density
      )
    ),
    bins = 30,
    fill = COL_PREDICTED,
    colour = "white",
    linewidth = 0.25
  ) +
  ggplot2::geom_density(
    colour = "#4A90C2",
    linewidth = 0.65
  ) +
  ggplot2::geom_vline(
    ggplot2::aes(
      xintercept = observed
    ),
    colour = COL_OBSERVED,
    linewidth = 0.95
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(
      response
    ),
    cols = ggplot2::vars(
      statistic
    ),
    scales = "free_x"
  ) +
  ggplot2::labs(
    x = "Posterior predictive statistic",
    y = "Density",
    caption = paste(
      "Vertical lines represent the observed statistic;",
      "distributions represent 1,000 posterior predictive datasets."
    )
  ) +
  theme_publication(
    base_size = 14
  ) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(
      colour = "grey92",
      linewidth = 0.3
    )
  )

figure_5 <- add_optional_title(
  figure_5,
  title_text = "Posterior predictive checks for location and spread",
  subtitle_text = paste(
    "Observed means and standard deviations compared",
    "with posterior predictive datasets"
  )
)

print(
  figure_5
)

save_publication_figure(
  plot_object = figure_5,
  filename = "Figure_5_PPC_mean_and_standard_deviation",
  width = 7.2,
  height = 5.8
)

write.csv(
  moment_plot_data,
  file.path(
    FIGURE_DIR,
    "Figure_5_PPC_moment_data.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# FIGURE 6. POSTERIOR RESIDUAL CORRELATION
# =============================================================================

residual_correlation_parameters <- grep(
  "rescor",
  names(posterior_draws),
  value = TRUE,
  ignore.case = TRUE
)

if (length(residual_correlation_parameters) > 0L) {

  residual_parameter <- residual_correlation_parameters[1]

  residual_draws <- as.numeric(
    posterior_draws[[residual_parameter]]
  )

  residual_density <- stats::density(
    residual_draws,
    from = -1,
    to = 1,
    n = 1024
  )

  residual_density_data <- data.frame(
    correlation = residual_density$x,
    density = residual_density$y
  )

  residual_median <- stats::median(
    residual_draws,
    na.rm = TRUE
  )

  residual_lower_95 <- unname(
    stats::quantile(
      residual_draws,
      probs = 0.025,
      na.rm = TRUE
    )
  )

  residual_upper_95 <- unname(
    stats::quantile(
      residual_draws,
      probs = 0.975,
      na.rm = TRUE
    )
  )

  residual_lower_50 <- unname(
    stats::quantile(
      residual_draws,
      probs = 0.25,
      na.rm = TRUE
    )
  )

  residual_upper_50 <- unname(
    stats::quantile(
      residual_draws,
      probs = 0.75,
      na.rm = TRUE
    )
  )

  interval_y <- max(
    residual_density_data$density
  ) * 0.06

  annotation_y <- max(
    residual_density_data$density
  ) * 0.92

  figure_6 <- ggplot2::ggplot(
    residual_density_data,
    ggplot2::aes(
      x = correlation,
      y = density
    )
  ) +
    ggplot2::geom_area(
      fill = COL_CARBON,
      alpha = 0.22
    ) +
    ggplot2::geom_line(
      colour = COL_CARBON,
      linewidth = 0.9
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.45,
      colour = COL_REFERENCE
    ) +
    ggplot2::annotate(
      geom = "segment",
      x = residual_lower_95,
      xend = residual_upper_95,
      y = interval_y,
      yend = interval_y,
      linewidth = 1.0,
      colour = COL_INTERVAL_95
    ) +
    ggplot2::annotate(
      geom = "segment",
      x = residual_lower_50,
      xend = residual_upper_50,
      y = interval_y,
      yend = interval_y,
      linewidth = 3.2,
      lineend = "round",
      colour = COL_INTERVAL_50
    ) +
    ggplot2::annotate(
      geom = "point",
      x = residual_median,
      y = interval_y,
      shape = 21,
      size = 2.8,
      stroke = 0.55,
      fill = "white",
      colour = "black"
    ) +
    ggplot2::annotate(
      geom = "text",
      x = -0.96,
      y = annotation_y,
      hjust = 0,
      vjust = 1,
      family = FONT_FAMILY,
      size = 4.25,
      label = sprintf(
        paste0(
          "Median = %.2f\n",
          "95%% CrI [%.2f, %.2f]"
        ),
        residual_median,
        residual_lower_95,
        residual_upper_95
      )
    ) +
    ggplot2::scale_x_continuous(
      limits = c(-1, 1),
      breaks = seq(
        -1,
        1,
        by = 0.25
      )
    ) +
    ggplot2::labs(
      x = expression(
        rho[residual]
      ),
      y = "Posterior density",
      caption = paste(
        "Residual correlation after accounting for reach, campaign,",
        "feeding guild and species-level effects."
      )
    ) +
    theme_publication(
      base_size = 14
    )

  figure_6 <- add_optional_title(
    figure_6,
    title_text = "Residual correlation between δ13C and δ15N",
    subtitle_text = paste(
      "Posterior distribution of the",
      "multivariate residual correlation"
    )
  )

  print(
    figure_6
  )

  save_publication_figure(
    plot_object = figure_6,
    filename = "Figure_6_posterior_residual_correlation",
    width = 6.5,
    height = 4.3
  )

  write.csv(
    data.frame(
      parameter = residual_parameter,
      median = residual_median,
      lower_95 = residual_lower_95,
      upper_95 = residual_upper_95,
      lower_50 = residual_lower_50,
      upper_50 = residual_upper_50
    ),
    file.path(
      FIGURE_DIR,
      "Figure_6_residual_correlation_summary.csv"
    ),
    row.names = FALSE
  )

} else {

  warning(
    "No residual-correlation parameter was found; Figure 6 was skipped.",
    call. = FALSE
  )
}


# =============================================================================
# FIGURE 7. OBSERVED VERSUS FITTED VALUES
# =============================================================================

create_observed_fitted_data <- function(
  response_id,
  observed_raw,
  raw_mean,
  raw_sd,
  response_label
) {

  fitted_draws <- brms::posterior_epred(
    fit_mv,
    resp = response_id,
    ndraws = 1000
  )

  fitted_draws <- as_prediction_matrix(
    fitted_draws
  )

  fitted_standardized <- apply(
    fitted_draws,
    2,
    stats::median,
    na.rm = TRUE
  )

  fitted_raw <- (
    fitted_standardized *
    raw_sd
  ) +
    raw_mean

  data.frame(
    observed = observed_raw,
    fitted = fitted_raw,
    response = response_label,
    Campaign = as.character(
      fish_model$Campaign
    ),
    stringsAsFactors = FALSE
  )
}

observed_fitted_d13C <- create_observed_fitted_data(
  response_id = resp_d13C,
  observed_raw = fish_model$d13C_permil,
  raw_mean = d13C_mean,
  raw_sd = d13C_sd,
  response_label = "δ13C"
)

observed_fitted_d15N <- create_observed_fitted_data(
  response_id = resp_d15N,
  observed_raw = fish_model$d15N_permil,
  raw_mean = d15N_mean,
  raw_sd = d15N_sd,
  response_label = "δ15N"
)

observed_fitted_data <- dplyr::bind_rows(
  observed_fitted_d13C,
  observed_fitted_d15N
) |>
  dplyr::mutate(
    Campaign_label = dplyr::recode(
      Campaign,
      Wet_2023 = "Wet campaign",
      Dry_2023 = "Dry campaign"
    )
  )

figure_7 <- ggplot2::ggplot(
  observed_fitted_data,
  ggplot2::aes(
    x = observed,
    y = fitted,
    colour = Campaign_label,
    shape = Campaign_label
  )
) +
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.55,
    colour = COL_REFERENCE
  ) +
  ggplot2::geom_point(
    size = 1.8,
    alpha = 0.68
  ) +
  ggplot2::facet_wrap(
    ~ response,
    scales = "free"
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "Dry campaign" = COL_DRY,
      "Wet campaign" = COL_WET
    ),
    breaks = c(
      "Dry campaign",
      "Wet campaign"
    )
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "Dry campaign" = 17,
      "Wet campaign" = 16
    ),
    breaks = c(
      "Dry campaign",
      "Wet campaign"
    )
  ) +
  ggplot2::labs(
    x = "Observed isotope value (‰)",
    y = "Posterior median fitted value (‰)",
    colour = "Campaign",
    shape = "Campaign",
    caption = "The dashed line represents perfect agreement between observed and fitted values."
  ) +
  theme_publication(
    base_size = 14
  ) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(
      colour = "grey92",
      linewidth = 0.3
    )
  )

figure_7 <- add_optional_title(
  figure_7,
  title_text = "Observed and fitted isotope values",
  subtitle_text = "Posterior median fitted values including species-level effects"
)

print(
  figure_7
)

save_publication_figure(
  plot_object = figure_7,
  filename = "Figure_7_observed_vs_fitted",
  width = 7.2,
  height = 4.2
)

write.csv(
  observed_fitted_data,
  file.path(
    FIGURE_DIR,
    "Figure_7_observed_fitted_data.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 12. SAVE SESSION INFORMATION
# =============================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    FIGURE_DIR,
    "visualization_session_information.txt"
  )
)


# =============================================================================
# 13. COMPLETION MESSAGE
# =============================================================================

cat(
  "\n============================================================\n",
  "PUBLICATION-QUALITY VISUALIZATIONS COMPLETED\n",
  "============================================================\n",
  "The Bayesian model was not refitted.\n",
  "\nFigures were saved as:\n",
  "  - 600-dpi PNG\n",
  "  - 600-dpi LZW-compressed TIFF\n",
  "  - vector PDF\n",
  "\nOutput folder:\n",
  FIGURE_DIR,
  "\n\nCreated figures:\n",
  "  Figure 1: Observed isotope distributions\n",
  "  Figure 2: Posterior fixed effects\n",
  "  Figure 3: Adjusted reach × campaign estimates\n",
  "  Figure 4a: δ13C density posterior predictive check\n",
  "  Figure 4b: δ15N density posterior predictive check\n",
  "  Figure 5: Mean and SD posterior predictive checks\n",
  "  Figure 6: Posterior residual correlation\n",
  "  Figure 7: Observed versus fitted values\n",
  "============================================================\n",
  sep = ""
)

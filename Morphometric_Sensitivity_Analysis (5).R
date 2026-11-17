# Morphometric sensitivity analysis for lower Han River fish isotopes
#
# Purpose
# -------
# Test whether wet-to-dry changes in mean fish size at the reach-species level
# explain apparent within-species changes in mean d13C or d15N.
#
# IMPORTANT INTERPRETATION
# Individual isotope and morphometric records cannot be linked. Therefore,
# this is a reach-species-pair sensitivity analysis, not an individual-level
# covariate adjustment.
#
# Before running:
# 1. Keep Morphometric_Sensitivity_Analysis_Data.xlsx anywhere on your computer.
#    If it is not in RStudio's working directory, the script will ask you to
#    select it using a file-selection window.
# 2. Install missing packages if needed, for example:
#    install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "broom",
#                       "car", "lmtest", "sandwich", "writexl"))
# 3. Source the complete script. Results are saved beside the selected workbook;
#    in Windows, the completed output folder opens automatically.


# 1. Setup -----------------------------------------------------------------

required_packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "broom",
  "car", "lmtest", "sandwich", "writexl"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install these packages before running the analysis: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# broom, car, lmtest, sandwich, and writexl are called explicitly with the
# package::function form below. Keeping them off the search path prevents
# function-name conflicts such as car::recode masking dplyr::recode.

options(contrasts = c("contr.treatment", "contr.poly"))

# The script first looks for the workbook in the current working directory.
# If it is not there, an interactive file-selection window opens so that the
# workbook can be selected from Downloads, Documents, Desktop, or elsewhere.
default_input_file <- "Morphometric_Sensitivity_Analysis_Data.xlsx"

if (file.exists(default_input_file)) {
  input_file <- default_input_file
} else if (interactive()) {
  message(
    "The workbook was not found in the current working directory. ",
    "Please select Morphometric_Sensitivity_Analysis_Data.xlsx."
  )
  input_file <- file.choose()
} else {
  stop(
    "Input workbook not found. Place ", default_input_file,
    " in the working directory or set input_file to its full path."
  )
}

if (!file.exists(input_file)) {
  stop("Selected input file does not exist: ", input_file)
}

input_file <- normalizePath(input_file, mustWork = TRUE)

# Save every output beside the selected Excel workbook. On Windows, the output
# folder opens automatically when the analysis finishes.
output_dir <- file.path(dirname(input_file), "morphometric_sensitivity_results")
figure_dir <- file.path(output_dir, "publication_figures")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)


# 2. Helper functions -------------------------------------------------------

first_nonmissing <- function(x) {
  keep <- !is.na(x) & trimws(as.character(x)) != ""
  if (!any(keep)) return(NA_character_)
  as.character(x[which(keep)[1]])
}

mode_nonmissing <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & trimws(x) != ""]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

sample_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x)
}

correlation_row <- function(data, x, y, analysis_set, response_label) {
  dat <- data %>%
    filter(is.finite(.data[[x]]), is.finite(.data[[y]]))

  test <- cor.test(dat[[x]], dat[[y]], method = "pearson")

  tibble(
    analysis_set = analysis_set,
    response = response_label,
    predictor = "Dry-minus-wet change in mean total length (mm)",
    n_pairs = nrow(dat),
    pearson_r = unname(test$estimate),
    conf_low = unname(test$conf.int[1]),
    conf_high = unname(test$conf.int[2]),
    p_value = test$p.value
  )
}

prepare_model_data <- function(data) {
  out <- data %>%
    mutate(
      Reach_code = factor(Reach_code),
      Feeding_guild = factor(Feeding_guild),
      Delta_total_length_10mm = Delta_total_length_mm / 10,
      Delta_total_length_z = as.numeric(scale(Delta_total_length_mm))
    )

  if ("HC" %in% levels(out$Reach_code)) {
    out$Reach_code <- relevel(out$Reach_code, ref = "HC")
  }
  if ("Omnivorous" %in% levels(out$Feeding_guild)) {
    out$Feeding_guild <- relevel(out$Feeding_guild, ref = "Omnivorous")
  }
  out
}

tidy_coeftest <- function(test_matrix, model_name, variance_type) {
  mat <- unclass(test_matrix)
  out <- tibble(
    term = rownames(mat),
    estimate = mat[, 1],
    std_error = mat[, 2],
    statistic = mat[, 3],
    p_value = mat[, 4]
  )
  out %>%
    mutate(
      conf_low = estimate - qnorm(0.975) * std_error,
      conf_high = estimate + qnorm(0.975) * std_error,
      model = model_name,
      variance_type = variance_type,
      .before = 1
    )
}

vif_table <- function(model, model_name) {
  result <- car::vif(model)

  if (is.matrix(result)) {
    tibble(
      model = model_name,
      term = rownames(result),
      GVIF = result[, "GVIF"],
      df = result[, "Df"],
      GVIF_adjusted = result[, "GVIF^(1/(2*Df))"]
    )
  } else {
    tibble(
      model = model_name,
      term = names(result),
      GVIF = as.numeric(result),
      df = 1,
      GVIF_adjusted = sqrt(as.numeric(result))
    )
  }
}

diagnostic_row <- function(fit, model_name) {
  # Use `fit` for the lm object because the output table also has a column
  # named `model`. This avoids tibble's sequential evaluation masking the
  # fitted model object after the first column is created.
  residual_values <- residuals(fit)
  cook_values <- cooks.distance(fit)
  shapiro_p <- if (length(residual_values) >= 3 && length(residual_values) <= 5000) {
    shapiro.test(residual_values)$p.value
  } else {
    NA_real_
  }

  glance_row <- broom::glance(fit)

  tibble(
    model = model_name,
    n_pairs = nobs(fit),
    r_squared = glance_row$r.squared,
    adjusted_r_squared = glance_row$adj.r.squared,
    residual_sd = sigma(fit),
    AIC = AIC(fit),
    breusch_pagan_p = lmtest::bptest(fit)$p.value,
    shapiro_wilk_p = shapiro_p,
    max_cooks_distance = max(cook_values, na.rm = TRUE),
    cooks_cutoff_4_over_n = 4 / nobs(fit),
    influential_n = sum(cook_values > 4 / nobs(fit), na.rm = TRUE)
  )
}

format_p_value <- function(p_value) {
  if (!is.finite(p_value)) return("p = NA")
  if (p_value < 0.001) return("p < 0.001")
  paste0("p = ", formatC(p_value, format = "f", digits = 3))
}

publication_theme <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.35),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      axis.title = element_text(color = "black", size = base_size + 0.5),
      axis.text = element_text(color = "black", size = base_size),
      strip.background = element_blank(),
      strip.text = element_text(color = "black", face = "bold", size = base_size + 0.5),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 0.5),
      legend.key = element_blank(),
      panel.spacing = grid::unit(8, "pt"),
      plot.margin = margin(6, 8, 6, 6)
    )
}

save_publication_figure <- function(plot, file_stem, width_mm, height_mm) {
  # Scatterplots with points, regression bands, lines, and labels are treated
  # as combination artwork. TIFF and PNG are exported at 600 dpi, exceeding
  # Elsevier's general 500-dpi requirement for combination artwork. A vector
  # PDF is also supplied for lossless scaling and editing.
  tiff_file <- file.path(figure_dir, paste0(file_stem, ".tif"))
  png_file <- file.path(figure_dir, paste0(file_stem, ".png"))
  pdf_file <- file.path(figure_dir, paste0(file_stem, ".pdf"))

  ggsave(
    filename = tiff_file,
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600,
    device = "tiff",
    compression = "lzw",
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    filename = png_file,
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600,
    device = "png",
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    filename = pdf_file,
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = grDevices::pdf,
    useDingbats = FALSE,
    bg = "white",
    limitsize = FALSE
  )

  invisible(c(TIFF = tiff_file, PNG = png_file, PDF = pdf_file))
}

prepare_publication_plot_data <- function(data) {
  complete_data <- data %>%
    filter(
      Morphometric_complete,
      is.finite(Delta_total_length_mm),
      is.finite(Delta_d13C_permil),
      is.finite(Delta_d15N_permil)
    )

  bind_rows(
    complete_data %>% mutate(Analysis_subset = "All complete pairs"),
    complete_data %>%
      filter(Count_matched_both_campaigns) %>%
      mutate(Analysis_subset = "Count-matched pairs")
  ) %>%
    pivot_longer(
      cols = c(Delta_d13C_permil, Delta_d15N_permil),
      names_to = "Isotope",
      values_to = "Delta_isotope_permil"
    ) %>%
    mutate(
      Analysis_subset = factor(
        Analysis_subset,
        levels = c("All complete pairs", "Count-matched pairs")
      ),
      Isotope = dplyr::recode(
        Isotope,
        Delta_d13C_permil = "delta^{13}*C",
        Delta_d15N_permil = "delta^{15}*N"
      ),
      Isotope = factor(Isotope, levels = c("delta^{13}*C", "delta^{15}*N"))
    )
}

publication_panel_statistics <- function(plot_data) {
  stats <- plot_data %>%
    group_by(Isotope, Analysis_subset) %>%
    summarise(
      n_pairs = n(),
      pearson_r = cor(
        Delta_total_length_mm,
        Delta_isotope_permil,
        method = "pearson"
      ),
      p_value = cor.test(
        Delta_total_length_mm,
        Delta_isotope_permil,
        method = "pearson"
      )$p.value,
      .groups = "drop"
    ) %>%
    arrange(Isotope, Analysis_subset) %>%
    mutate(
      panel_label = LETTERS[seq_len(n())],
      statistics_label = paste0(
        "r = ", sprintf("%.2f", pearson_r),
        "; ", vapply(p_value, format_p_value, character(1)),
        "; n = ", n_pairs
      )
    )

  stats
}

make_publication_plot <- function(plot_data, panel_stats, y_axis_label) {
  ggplot(
    plot_data,
    aes(
      x = Delta_total_length_mm,
      y = Delta_isotope_permil,
      color = Reach_code,
      shape = Feeding_guild
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey75") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey75") +
    geom_smooth(
      data = plot_data,
      mapping = aes(
        x = Delta_total_length_mm,
        y = Delta_isotope_permil,
        group = 1
      ),
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      color = "black",
      fill = "grey82",
      alpha = 0.55,
      linewidth = 0.65,
      inherit.aes = FALSE
    ) +
    geom_point(size = 2.0, stroke = 0.55, alpha = 0.88) +
    geom_text(
      data = panel_stats,
      mapping = aes(x = -Inf, y = Inf, label = panel_label),
      inherit.aes = FALSE,
      hjust = -0.15,
      vjust = 1.25,
      family = "sans",
      fontface = "bold",
      size = 3.4
    ) +
    geom_text(
      data = panel_stats,
      mapping = aes(x = Inf, y = Inf, label = statistics_label),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.25,
      family = "sans",
      size = 2.8
    ) +
    scale_color_manual(
      values = c(HC = "#0072B2", QJ = "#D55E00", SY = "#009E73"),
      drop = FALSE
    ) +
    labs(
      x = "Change in mean total length (dry - wet, mm)",
      y = y_axis_label,
      color = "Reach",
      shape = "Feeding guild"
    ) +
    guides(
      color = guide_legend(order = 1, nrow = 1),
      shape = guide_legend(order = 2, nrow = 1)
    ) +
    coord_cartesian(clip = "off") +
    publication_theme(base_size = 9)
}


# 3. Read and clean the prepared workbook ----------------------------------

fish_isotopes <- read_excel(input_file, sheet = "Fish_Isotopes") %>%
  filter(
    Reach_code %in% c("HC", "QJ", "SY"),
    Campaign %in% c("Wet_2023", "Dry_2023"),
    is.finite(d13C_permil),
    is.finite(d15N_permil)
  )

# Do not remove rows flagged only because the original trophic-level column
# was below 1. Those flags concern the old trophic-position calculation, not
# the measured d13C or d15N values used here.

fish_morphometrics <- read_excel(input_file, sheet = "Fish_Morphometrics") %>%
  filter(
    Principal_isotope_analysis == "Yes",
    Reach_code %in% c("HC", "QJ", "SY"),
    Campaign %in% c("Wet_2023", "Dry_2023")
  )


# 4. Build reach-species-campaign summaries --------------------------------

isotope_summary <- fish_isotopes %>%
  group_by(Reach_code, Reach, Scientific_name, Campaign) %>%
  summarise(
    Common_name = first_nonmissing(Common_name_clean),
    Feeding_guild = mode_nonmissing(Feeding_guild),
    Isotope_n = n(),
    Mean_d13C_permil = mean(d13C_permil, na.rm = TRUE),
    Mean_d15N_permil = mean(d15N_permil, na.rm = TRUE),
    Isotope_QC_flag_n = sum(QC_status != "OK", na.rm = TRUE),
    .groups = "drop"
  )

morphometric_summary <- fish_morphometrics %>%
  group_by(Reach_code, Reach, Scientific_name, Campaign) %>%
  summarise(
    Morphometric_n = n(),
    Mean_total_length_mm = mean(Total_length_mm, na.rm = TRUE),
    SD_total_length_mm = sample_sd(Total_length_mm),
    Mean_body_length_mm = mean(Body_length_mm, na.rm = TRUE),
    Mean_body_weight_g = mean(Body_weight_g, na.rm = TRUE),
    SD_body_weight_g = sample_sd(Body_weight_g),
    .groups = "drop"
  ) %>%
  mutate(
    across(
      c(
        Mean_total_length_mm, SD_total_length_mm, Mean_body_length_mm,
        Mean_body_weight_g, SD_body_weight_g
      ),
      ~ ifelse(is.nan(.x), NA_real_, .x)
    )
  )

shared_pair_keys <- isotope_summary %>%
  distinct(Reach_code, Reach, Scientific_name, Campaign) %>%
  count(Reach_code, Reach, Scientific_name, name = "Campaign_n") %>%
  filter(Campaign_n == 2) %>%
  select(-Campaign_n)

campaign_summary <- isotope_summary %>%
  semi_join(
    shared_pair_keys,
    by = c("Reach_code", "Reach", "Scientific_name")
  ) %>%
  left_join(
    morphometric_summary,
    by = c("Reach_code", "Reach", "Scientific_name", "Campaign")
  ) %>%
  mutate(Pair_ID = paste(Reach_code, Scientific_name, sep = "|"), .before = 1) %>%
  arrange(Reach_code, Scientific_name, Campaign)


# 5. Create one row per shared reach-species pair ---------------------------

pair_metadata <- campaign_summary %>%
  group_by(Pair_ID, Reach_code, Reach, Scientific_name) %>%
  summarise(
    Common_name = first_nonmissing(Common_name),
    Feeding_guild = mode_nonmissing(Feeding_guild),
    .groups = "drop"
  )

isotope_wide <- campaign_summary %>%
  select(
    Pair_ID, Reach_code, Reach, Scientific_name, Campaign,
    Isotope_n, Mean_d13C_permil, Mean_d15N_permil, Isotope_QC_flag_n
  ) %>%
  pivot_wider(
    names_from = Campaign,
    values_from = c(
      Isotope_n, Mean_d13C_permil, Mean_d15N_permil, Isotope_QC_flag_n
    ),
    names_glue = "{.value}_{Campaign}"
  )

morphometric_wide <- campaign_summary %>%
  select(
    Pair_ID, Reach_code, Reach, Scientific_name, Campaign,
    Morphometric_n, Mean_total_length_mm, Mean_body_length_mm,
    Mean_body_weight_g
  ) %>%
  pivot_wider(
    names_from = Campaign,
    values_from = c(
      Morphometric_n, Mean_total_length_mm, Mean_body_length_mm,
      Mean_body_weight_g
    ),
    names_glue = "{.value}_{Campaign}"
  )

morphometric_pairs <- pair_metadata %>%
  left_join(
    isotope_wide,
    by = c("Pair_ID", "Reach_code", "Reach", "Scientific_name")
  ) %>%
  left_join(
    morphometric_wide,
    by = c("Pair_ID", "Reach_code", "Reach", "Scientific_name")
  ) %>%
  mutate(
    Delta_d13C_permil = Mean_d13C_permil_Dry_2023 - Mean_d13C_permil_Wet_2023,
    Delta_d15N_permil = Mean_d15N_permil_Dry_2023 - Mean_d15N_permil_Wet_2023,
    Delta_total_length_mm = Mean_total_length_mm_Dry_2023 - Mean_total_length_mm_Wet_2023,
    Delta_body_length_mm = Mean_body_length_mm_Dry_2023 - Mean_body_length_mm_Wet_2023,
    Delta_body_weight_g = Mean_body_weight_g_Dry_2023 - Mean_body_weight_g_Wet_2023,
    Morphometric_complete = if_all(
      c(
        Mean_total_length_mm_Wet_2023,
        Mean_total_length_mm_Dry_2023,
        Mean_body_weight_g_Wet_2023,
        Mean_body_weight_g_Dry_2023
      ),
      ~ !is.na(.x)
    ),
    Count_matched_both_campaigns =
      Morphometric_complete &
      Isotope_n_Wet_2023 == Morphometric_n_Wet_2023 &
      Isotope_n_Dry_2023 == Morphometric_n_Dry_2023,
    Analysis_set = case_when(
      Count_matched_both_campaigns ~ "Count-matched",
      Morphometric_complete ~ "Morphometric-complete",
      TRUE ~ "Morphometric-incomplete"
    )
  ) %>%
  arrange(Reach_code, Scientific_name)


# 6. Verify the intended analysis population -------------------------------

verification <- tibble(
  check = c(
    "Shared reach-species pairs",
    "Pairs with length and weight in both campaigns",
    "Pairs with isotope and morphometric counts matched in both campaigns"
  ),
  observed = c(
    nrow(morphometric_pairs),
    sum(morphometric_pairs$Morphometric_complete),
    sum(morphometric_pairs$Count_matched_both_campaigns)
  ),
  expected = c(85, 81, 51)
) %>%
  mutate(status = if_else(observed == expected, "PASS", "CHECK"))

print(verification)

if (!all(verification$status == "PASS")) {
  stop(
    "The reconstructed analysis population does not match 85/81/51. ",
    "Check whether the input workbook or filtering rules changed."
  )
}


# 7. Descriptive correlations ----------------------------------------------

complete_pairs <- morphometric_pairs %>% filter(Morphometric_complete)
count_matched_pairs <- morphometric_pairs %>% filter(Count_matched_both_campaigns)

correlation_results <- bind_rows(
  correlation_row(
    complete_pairs,
    "Delta_d13C_permil",
    "Delta_total_length_mm",
    "All morphometric-complete pairs",
    "Dry-minus-wet change in mean d13C (permil)"
  ),
  correlation_row(
    count_matched_pairs,
    "Delta_d13C_permil",
    "Delta_total_length_mm",
    "Count-matched subset",
    "Dry-minus-wet change in mean d13C (permil)"
  ),
  correlation_row(
    complete_pairs,
    "Delta_d15N_permil",
    "Delta_total_length_mm",
    "All morphometric-complete pairs",
    "Dry-minus-wet change in mean d15N (permil)"
  ),
  correlation_row(
    count_matched_pairs,
    "Delta_d15N_permil",
    "Delta_total_length_mm",
    "Count-matched subset",
    "Dry-minus-wet change in mean d15N (permil)"
  )
)

print(correlation_results)


# 8. Pair-level linear models ----------------------------------------------

# The primary coefficient is scaled per 10-mm dry-minus-wet length change.
# Reach and feeding guild are categorical adjustment variables.

model_data_all <- prepare_model_data(complete_pairs)
model_data_matched <- prepare_model_data(count_matched_pairs)

model_specs <- list(
  list(
    name = "d13C_all_complete",
    response = "Delta_d13C_permil",
    data = model_data_all,
    fit = lm(
      Delta_d13C_permil ~ Delta_total_length_10mm + Reach_code + Feeding_guild,
      data = model_data_all
    )
  ),
  list(
    name = "d15N_all_complete",
    response = "Delta_d15N_permil",
    data = model_data_all,
    fit = lm(
      Delta_d15N_permil ~ Delta_total_length_10mm + Reach_code + Feeding_guild,
      data = model_data_all
    )
  ),
  list(
    name = "d13C_count_matched",
    response = "Delta_d13C_permil",
    data = model_data_matched,
    fit = lm(
      Delta_d13C_permil ~ Delta_total_length_10mm + Reach_code + Feeding_guild,
      data = model_data_matched
    )
  ),
  list(
    name = "d15N_count_matched",
    response = "Delta_d15N_permil",
    data = model_data_matched,
    fit = lm(
      Delta_d15N_permil ~ Delta_total_length_10mm + Reach_code + Feeding_guild,
      data = model_data_matched
    )
  )
)

model_coefficients_ols <- bind_rows(lapply(model_specs, function(spec) {
  broom::tidy(spec$fit, conf.int = TRUE) %>%
    mutate(
      model = spec$name,
      variance_type = "Ordinary least squares",
      .before = 1
    )
}))

model_coefficients_hc3 <- bind_rows(lapply(model_specs, function(spec) {
  robust_test <- lmtest::coeftest(
    spec$fit,
    vcov. = sandwich::vcovHC(spec$fit, type = "HC3")
  )
  tidy_coeftest(robust_test, spec$name, "HC3 heteroskedasticity-consistent")
}))

# Species may appear in more than one reach. Cluster-robust standard errors
# provide a sensitivity analysis for that repeated-species structure.
model_coefficients_species_cluster <- bind_rows(lapply(model_specs, function(spec) {
  clustered_test <- lmtest::coeftest(
    spec$fit,
    vcov. = sandwich::vcovCL(
      spec$fit,
      cluster = spec$data$Scientific_name,
      type = "HC1"
    )
  )
  tidy_coeftest(clustered_test, spec$name, "Clustered by scientific name")
}))

model_diagnostics <- bind_rows(lapply(model_specs, function(spec) {
  diagnostic_row(spec$fit, spec$name)
}))

model_vif <- bind_rows(lapply(model_specs, function(spec) {
  vif_table(spec$fit, spec$name)
}))


# 9. Influence sensitivity --------------------------------------------------

# Refit after excluding observations with Cook's distance > 4/n. This is a
# sensitivity check only; it should not silently replace the primary models.

influence_sensitivity <- bind_rows(lapply(model_specs, function(spec) {
  cook <- cooks.distance(spec$fit)
  cutoff <- 4 / nobs(spec$fit)
  keep <- cook <= cutoff

  refit <- lm(formula(spec$fit), data = spec$data[keep, , drop = FALSE])

  broom::tidy(refit, conf.int = TRUE) %>%
    filter(term == "Delta_total_length_10mm") %>%
    mutate(
      model = spec$name,
      original_n = nobs(spec$fit),
      retained_n = nobs(refit),
      excluded_n = original_n - retained_n,
      cooks_cutoff = cutoff,
      .before = 1
    )
}))


# 10. Figures ---------------------------------------------------------------

publication_plot_data <- prepare_publication_plot_data(morphometric_pairs)
overview_statistics <- publication_panel_statistics(publication_plot_data)

# Main four-panel overview: full Elsevier double-column width (190 mm).
# Rows show isotopes; columns show all complete and count-matched analyses.
plot_overview <- make_publication_plot(
  publication_plot_data,
  overview_statistics,
  "Change in mean isotope value (dry - wet, permil)"
) +
  facet_grid(
    rows = vars(Isotope),
    cols = vars(Analysis_subset),
    scales = "free_y",
    labeller = labeller(Isotope = label_parsed)
  )

save_publication_figure(
  plot = plot_overview,
  file_stem = "FIG1_Morphometric_Sensitivity_Overview",
  width_mm = 190,
  height_mm = 145
)

# Two-panel d13C figure for supplementary use or separate submission.
plot_d13C_data <- publication_plot_data %>%
  filter(Isotope == "delta^{13}*C") %>%
  droplevels()
plot_d13C_statistics <- publication_panel_statistics(plot_d13C_data)

plot_d13C <- make_publication_plot(
  plot_d13C_data,
  plot_d13C_statistics,
  expression(Delta*delta^{13}*C~"(dry - wet, \u2030)")
) +
  facet_wrap(vars(Analysis_subset), nrow = 1)

save_publication_figure(
  plot = plot_d13C,
  file_stem = "FIGS1_Morphometric_Sensitivity_d13C",
  width_mm = 190,
  height_mm = 85
)

# Two-panel d15N figure for supplementary use or separate submission.
plot_d15N_data <- publication_plot_data %>%
  filter(Isotope == "delta^{15}*N") %>%
  droplevels()
plot_d15N_statistics <- publication_panel_statistics(plot_d15N_data)

plot_d15N <- make_publication_plot(
  plot_d15N_data,
  plot_d15N_statistics,
  expression(Delta*delta^{15}*N~"(dry - wet, \u2030)")
) +
  facet_wrap(vars(Analysis_subset), nrow = 1)

save_publication_figure(
  plot = plot_d15N,
  file_stem = "FIGS2_Morphometric_Sensitivity_d15N",
  width_mm = 190,
  height_mm = 85
)

figure_captions <- c(
  paste0(
    "Figure 1. Morphometric sensitivity of seasonal stable-isotope shifts. ",
    "Panels A and B show dry-minus-wet changes in mean d13C, and panels C and D ",
    "show changes in mean d15N. The left column includes all 81 reach-species pairs ",
    "with morphometric observations in both campaigns; the right column includes the ",
    "51 pairs with matching isotope and morphometric counts in both campaigns. Solid ",
    "lines are ordinary least-squares fits and grey bands are 95% confidence intervals. ",
    "Colors identify reaches and symbols identify feeding guilds. Each point represents ",
    "a reach-species pair, not an individual fish."
  ),
  paste0(
    "Figure S1. Relationship between dry-minus-wet changes in mean total length and ",
    "mean d13C for all morphometric-complete reach-species pairs (A) and the count-matched ",
    "subset (B). Solid lines are ordinary least-squares fits; grey bands are 95% confidence ",
    "intervals. Each point represents a reach-species pair."
  ),
  paste0(
    "Figure S2. Relationship between dry-minus-wet changes in mean total length and ",
    "mean d15N for all morphometric-complete reach-species pairs (A) and the count-matched ",
    "subset (B). Solid lines are ordinary least-squares fits; grey bands are 95% confidence ",
    "intervals. Each point represents a reach-species pair."
  )
)

writeLines(
  figure_captions,
  con = file.path(figure_dir, "Figure_captions.txt")
)

writeLines(
  c(
    "Publication figure specifications",
    "- Main overview size: 190 x 145 mm (Elsevier double-column width)",
    "- Supplementary figure size: 190 x 85 mm",
    "- TIFF: RGB, LZW-compressed, 600 dpi",
    "- PNG: RGB, 600 dpi",
    "- PDF: vector, fonts retained as standard sans serif",
    "- Palette: color-vision-accessible; shapes also distinguish feeding guilds",
    "- Figure titles/captions are supplied separately, not embedded in the artwork"
  ),
  con = file.path(figure_dir, "Artwork_specifications.txt")
)


# 11. Export results --------------------------------------------------------

write.csv(
  campaign_summary,
  file.path(output_dir, "Campaign_summary_170_rows.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  morphometric_pairs,
  file.path(output_dir, "Morphometric_pairs_85_rows.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  correlation_results,
  file.path(output_dir, "Correlation_results.csv"),
  row.names = FALSE
)

write.csv(
  overview_statistics,
  file.path(output_dir, "Figure_panel_statistics.csv"),
  row.names = FALSE
)

write.csv(
  model_coefficients_hc3,
  file.path(output_dir, "Model_coefficients_HC3.csv"),
  row.names = FALSE
)

writexl::write_xlsx(
  list(
    Verification = verification,
    Correlations = correlation_results,
    Figure_statistics = overview_statistics,
    Model_OLS = model_coefficients_ols,
    Model_HC3 = model_coefficients_hc3,
    Model_species_cluster = model_coefficients_species_cluster,
    Model_diagnostics = model_diagnostics,
    Model_VIF = model_vif,
    Influence_sensitivity = influence_sensitivity,
    Morphometric_pairs = morphometric_pairs,
    Campaign_summary = campaign_summary
  ),
  path = file.path(output_dir, "Morphometric_Sensitivity_Results.xlsx")
)


# 12. Create a concise manuscript-results draft -----------------------------

r_c13_all <- correlation_results %>%
  filter(
    analysis_set == "All morphometric-complete pairs",
    grepl("d13C", response, fixed = TRUE)
  )
r_c13_matched <- correlation_results %>%
  filter(
    analysis_set == "Count-matched subset",
    grepl("d13C", response, fixed = TRUE)
  )
r_n15_all <- correlation_results %>%
  filter(
    analysis_set == "All morphometric-complete pairs",
    grepl("d15N", response, fixed = TRUE)
  )
r_n15_matched <- correlation_results %>%
  filter(
    analysis_set == "Count-matched subset",
    grepl("d15N", response, fixed = TRUE)
  )

manuscript_text <- c(
  paste0(
    "Morphometric sensitivity analysis was conducted at the reach-species-pair level because ",
    "individual isotope and morphometric records could not be linked. Among 85 shared ",
    "reach-species pairs, 81 had length and weight observations in both campaigns and 51 had ",
    "matching isotope and morphometric counts in both campaigns."
  ),
  paste0(
    "Dry-minus-wet change in mean total length was weakly correlated with change in mean d13C ",
    "(r = ", sprintf("%.2f", r_c13_all$pearson_r), ", n = ", r_c13_all$n_pairs,
    ") and remained negligible in the count-matched subset (r = ",
    sprintf("%.2f", r_c13_matched$pearson_r), ", n = ", r_c13_matched$n_pairs, ")."
  ),
  paste0(
    "The corresponding correlation was more positive for d15N (r = ",
    sprintf("%.2f", r_n15_all$pearson_r), ", n = ", r_n15_all$n_pairs,
    ") but weakened in the count-matched subset (r = ",
    sprintf("%.2f", r_n15_matched$pearson_r), ", n = ", r_n15_matched$n_pairs, ")."
  ),
  paste0(
    "These results indicate that size differences are unlikely to explain the d13C shift, ",
    "whereas d15N change may be partly associated with ontogenetic differences. The analysis ",
    "does not constitute an individual-level size adjustment."
  )
)

writeLines(
  manuscript_text,
  con = file.path(output_dir, "Manuscript_results_draft.txt")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "R_session_info.txt")
)

cat("\nAnalysis complete. Outputs were written to:\n")
cat(normalizePath(output_dir), "\n\n")
print(correlation_results)

if (interactive()) {
  print(plot_overview)
}

# Open the completed output folder automatically in Windows/RStudio.
if (interactive() && identical(.Platform$OS.type, "windows")) {
  try(shell.exec(normalizePath(output_dir)), silent = TRUE)
}

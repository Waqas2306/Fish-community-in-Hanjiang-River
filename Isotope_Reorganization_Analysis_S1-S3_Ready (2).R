#!/usr/bin/env Rscript

# Reviewer-focused stable-isotope analysis for RStudio
# ====================================================
#
# This script reads the supplied Han River Excel workbook and produces:
#   1. Basal-resource seasonal-shift tables.
#   2. Baseline-relative fish d13C difference-in-differences.
#   3. Baseline-corrected trophic position with TDF sensitivity.
#   4. Species-level size-structure sensitivity.
#   5. Publication-ready PNG, PDF, and 600-dpi TIFF figures.
#   6. A results workbook and CSV tables.
#
# In RStudio, place this script and the input workbook in the same folder,
# edit INPUT_FILE below if necessary, and click Source.

options(stringsAsFactors = FALSE, scipen = 999)

# ------------------------------- User settings -------------------------------
# Resolve paths relative to this script, not RStudio's current working folder.
# This lets you extract the package anywhere and simply click Source.
get_script_directory <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
  }

  for (frame_number in rev(seq_len(sys.nframe()))) {
    source_file <- tryCatch(sys.frame(frame_number)$ofile, error = function(e) NULL)
    if (!is.null(source_file) && length(source_file) == 1 && nzchar(source_file)) {
      return(dirname(normalizePath(source_file, mustWork = FALSE)))
    }
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_file <- tryCatch(
      rstudioapi::getActiveDocumentContext()$path,
      error = function(e) ""
    )
    if (nzchar(active_file)) {
      return(dirname(normalizePath(active_file, mustWork = FALSE)))
    }
  }

  normalizePath(getwd(), mustWork = FALSE)
}

SCRIPT_DIRECTORY <- get_script_directory()
WORKBOOK_NAME <- "eb571c79-7e89-41f9-9915-0f393ddaedad.xlsx"
INPUT_FILE <- file.path(SCRIPT_DIRECTORY, WORKBOOK_NAME)
OUTPUT_DIR <- file.path(SCRIPT_DIRECTORY, "isotope_analysis_R_results")
BOOTSTRAP_DRAWS <- 5000
RANDOM_SEED <- 20260816
INSTALL_MISSING_PACKAGES <- TRUE

# Command-line use is also supported:
# Rscript run_isotope_reorganization_analysis.R input.xlsx output_folder 5000
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[1])) INPUT_FILE <- args[1]
if (length(args) >= 2 && nzchar(args[2])) OUTPUT_DIR <- args[2]
if (length(args) >= 3 && nzchar(args[3])) BOOTSTRAP_DRAWS <- as.integer(args[3])

required_packages <- c(
  "readxl", "dplyr", "tidyr", "purrr", "tibble", "readr",
  "ggplot2", "patchwork", "sandwich", "openxlsx"
)

missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  if (!INSTALL_MISSING_PACKAGES) {
    stop("Install these packages before running: ", paste(missing_packages, collapse = ", "))
  }
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(sandwich)
  library(openxlsx)
})

if (!file.exists(INPUT_FILE)) {
  workbook_candidates <- list.files(
    SCRIPT_DIRECTORY,
    pattern = "\\.xlsx$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  workbook_candidates <- workbook_candidates[
    !startsWith(basename(workbook_candidates), "~$")
  ]

  if (length(workbook_candidates) == 1) {
    INPUT_FILE <- workbook_candidates[1]
    message("Using the Excel workbook found beside the script: ", basename(INPUT_FILE))
  } else if (interactive()) {
    message(
      "The workbook was not found beside the script. ",
      "Please select Morphometric_Sensitivity_Analysis_Data(1).xlsx or the supplied analysis workbook."
    )
    selected_file <- tryCatch(file.choose(), error = function(e) "")
    if (nzchar(selected_file) && file.exists(selected_file)) INPUT_FILE <- selected_file
  }
}

if (!file.exists(INPUT_FILE)) {
  stop(
    "Input workbook not found. Extract the complete package so the R script and ",
    "Excel workbook are in the same folder, or run the script again and select ",
    "the workbook when prompted. Expected location: ",
    normalizePath(INPUT_FILE, mustWork = FALSE)
  )
}
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(RANDOM_SEED)

CAMPAIGNS <- c("Wet_2023", "Dry_2023")

# The source workbook uses HC, QJ, and SY. Convert them immediately after
# import so all result tables, models, and figures use the reporting codes
# S1, S2, and S3.
SOURCE_REACH_TO_SITE <- c(HC = "S1", QJ = "S2", SY = "S3")
REACH_ORDER <- c("S1", "S2", "S3")
REACH_NAMES <- c(S1 = "Hanchuan", S2 = "Qianjiang", S3 = "Shayang")
REACH_LABELS <- c(
  S1 = "S1 (Hanchuan)",
  S2 = "S2 (Qianjiang)",
  S3 = "S3 (Shayang)"
)
PATHWAY_ORDER <- c("Phytoplankton", "POM", "Gastropod")
SEASON_COLORS <- c(Wet_2023 = "#0072B2", Dry_2023 = "#D55E00")
REACH_COLORS <- c(S1 = "#0072B2", S2 = "#E69F00", S3 = "#009E73")
PATHWAY_COLORS <- c(
  "Raw fish" = "#333333", Phytoplankton = "#009E73",
  POM = "#CC79A7", Gastropod = "#E69F00"
)

# ------------------------------- Functions -----------------------------------
clean_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.finite(x)]
}

standardize_site_code <- function(x) {
  x <- as.character(x)
  converted <- unname(SOURCE_REACH_TO_SITE[x])
  converted[is.na(converted)] <- x[is.na(converted)]
  converted
}

validate_site_mapping <- function(data, source_name) {
  expected_reach <- unname(REACH_NAMES[as.character(data$Reach_code)])
  mismatch <- is.na(expected_reach) |
    is.na(data$Reach) |
    trimws(as.character(data$Reach)) != expected_reach

  if (any(mismatch)) {
    bad_pairs <- unique(
      paste(data$Reach_code[mismatch], data$Reach[mismatch], sep = " = ")
    )
    stop(
      "Unexpected site-code mapping in ", source_name, ": ",
      paste(bad_pairs, collapse = ", "),
      ". Expected S1 = Hanchuan, S2 = Qianjiang, and S3 = Shayang."
    )
  }

  data
}

bootstrap_means <- function(x, B = BOOTSTRAP_DRAWS) {
  x <- clean_numeric(x)
  if (length(x) == 0) return(rep(NA_real_, B))
  replicate(B, mean(sample(x, size = length(x), replace = TRUE)))
}

bootstrap_ci <- function(draws, level = 0.95) {
  alpha <- (1 - level) / 2
  as.numeric(quantile(draws, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE))
}

mean_summary <- function(x) {
  x <- clean_numeric(x)
  if (length(x) == 0) {
    return(tibble(n = 0L, mean = NA_real_, sd = NA_real_, se = NA_real_,
                  ci_low = NA_real_, ci_high = NA_real_))
  }
  boot <- bootstrap_means(x)
  ci <- bootstrap_ci(boot)
  tibble(
    n = length(x),
    mean = mean(x),
    sd = ifelse(length(x) > 1, sd(x), NA_real_),
    se = ifelse(length(x) > 1, sd(x) / sqrt(length(x)), NA_real_),
    ci_low = ci[1],
    ci_high = ci[2]
  )
}

select_pathway <- function(data, pathway) {
  if (pathway == "Phytoplankton") {
    return(filter(data, Resource_class == "Phytoplankton"))
  }
  if (pathway == "POM") {
    return(filter(data, Resource_class == "POM/detritus"))
  }
  if (pathway == "Gastropod") {
    return(filter(data, Resource_type_original == "Rusty-patina Ring-ribbed Snail"))
  }
  stop("Unknown pathway: ", pathway)
}

save_submission_plot <- function(plot_object, stem, width, height) {
  ggsave(file.path(OUTPUT_DIR, paste0(stem, ".png")), plot_object,
         width = width, height = height, units = "in", dpi = 600, bg = "white")
  ggsave(file.path(OUTPUT_DIR, paste0(stem, ".pdf")), plot_object,
         width = width, height = height, units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(OUTPUT_DIR, paste0(stem, ".tiff")), plot_object,
         width = width, height = height, units = "in", dpi = 600,
         device = "tiff", compression = "lzw", bg = "white")
}

theme_submission <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.35),
      panel.border = element_blank(),
      axis.line = element_line(color = "#555555", linewidth = 0.4),
      axis.title = element_text(face = "bold", size = base_size + 2),
      axis.text = element_text(size = base_size),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      strip.background = element_rect(fill = "#F3F5F6", color = NA),
      strip.text = element_text(face = "bold", size = base_size + 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.5)
    )
}

# ------------------------------- Read data ------------------------------------
fish <- read_excel(INPUT_FILE, sheet = "Fish_Isotopes") %>%
  mutate(
    Reach_code = standardize_site_code(Reach_code),
    d13C_permil = as.numeric(d13C_permil),
    d15N_permil = as.numeric(d15N_permil)
  ) %>%
  filter(
    Reach_code %in% REACH_ORDER,
    Campaign %in% CAMPAIGNS,
    is.finite(d13C_permil),
    is.finite(d15N_permil)
  ) %>%
  validate_site_mapping("Fish_Isotopes")

morph <- read_excel(INPUT_FILE, sheet = "Fish_Morphometrics") %>%
  mutate(Reach_code = standardize_site_code(Reach_code)) %>%
  filter(
    Reach_code %in% REACH_ORDER,
    Campaign %in% CAMPAIGNS,
    QC_status == "OK"
  ) %>%
  validate_site_mapping("Fish_Morphometrics")

resources_raw <- read_excel(INPUT_FILE, sheet = "Resource_Isotopes") %>%
  mutate(
    Reach_code = standardize_site_code(Reach_code),
    Campaign_as_supplied = Campaign
  ) %>%
  filter(Reach_code %in% REACH_ORDER) %>%
  validate_site_mapping("Resource_Isotopes")

# Source code 2306 is Wet_2023 and 2311 is Dry_2023. In the supplied file,
# three Hanchuan gastropod rows with code 2306 are labelled Wet_2024,
# Wet_2025, and Wet_2026. Standardize them for analysis without changing the
# source workbook, and retain the recoding in an audit table.
resources <- resources_raw %>%
  mutate(
    Original_time_code = as.numeric(Original_time_code),
    Campaign = case_when(
      Original_time_code == 2306 ~ "Wet_2023",
      Original_time_code == 2311 ~ "Dry_2023",
      TRUE ~ Campaign_as_supplied
    )
  )

campaign_audit <- resources %>%
  filter(Campaign != Campaign_as_supplied) %>%
  transmute(
    Resource_record_id, Reach_code, Reach, Original_time_code,
    Campaign_as_supplied, Campaign_analysis = Campaign,
    Resource_type_original, d15N_permil, d13C_permil,
    Audit_action = "Recoded from source code 2306 to Wet_2023"
  )

resources <- resources %>%
  filter(Campaign %in% CAMPAIGNS) %>%
  mutate(
    d13C_permil = as.numeric(d13C_permil),
    d15N_permil = as.numeric(d15N_permil)
  )

if (nrow(fish) != 609) {
  warning("Expected 609 fish isotope rows; found ", nrow(fish), ". Check the input workbook.")
}

# ------------------------- Fish and resource summaries ------------------------
fish_summary <- map_dfr(REACH_ORDER, function(reach) {
  map_dfr(CAMPAIGNS, function(campaign) {
    group <- filter(fish, Reach_code == reach, Campaign == campaign)
    bind_rows(
      mean_summary(group$d13C_permil) %>% mutate(Isotope = "d13C"),
      mean_summary(group$d15N_permil) %>% mutate(Isotope = "d15N")
    ) %>%
      mutate(Reach_code = reach, Reach = REACH_LABELS[[reach]], Campaign = campaign, .before = 1)
  })
})

resource_group_summary <- map_dfr(REACH_ORDER, function(reach) {
  map_dfr(PATHWAY_ORDER, function(pathway) {
    selected <- resources %>% filter(Reach_code == reach) %>% select_pathway(pathway)
    map_dfr(CAMPAIGNS, function(campaign) {
      group <- filter(selected, Campaign == campaign)
      bind_rows(
        mean_summary(group$d13C_permil) %>% mutate(Isotope = "d13C"),
        mean_summary(group$d15N_permil) %>% mutate(Isotope = "d15N")
      ) %>%
        mutate(
          Reach_code = reach, Reach = REACH_LABELS[[reach]],
          Pathway = pathway, Campaign = campaign, .before = 1
        )
    })
  })
})

resource_seasonal_shifts <- map_dfr(REACH_ORDER, function(reach) {
  map_dfr(PATHWAY_ORDER, function(pathway) {
    selected <- resources %>% filter(Reach_code == reach) %>% select_pathway(pathway)
    imap_dfr(c(d13C = "d13C_permil", d15N = "d15N_permil"), function(variable, isotope_label) {
      wet <- selected %>% filter(Campaign == "Wet_2023") %>% pull(all_of(variable)) %>% clean_numeric()
      dry <- selected %>% filter(Campaign == "Dry_2023") %>% pull(all_of(variable)) %>% clean_numeric()
      if (length(wet) > 0 && length(dry) > 0) {
        delta_draws <- bootstrap_means(dry) - bootstrap_means(wet)
        ci <- bootstrap_ci(delta_draws)
        tibble(
          Reach_code = reach, Reach = REACH_LABELS[[reach]], Pathway = pathway,
          Isotope = isotope_label,
          Wet_n = length(wet), Dry_n = length(dry),
          Wet_mean = mean(wet), Dry_mean = mean(dry),
          Delta_dry_minus_wet = mean(dry) - mean(wet),
          CI_low = ci[1], CI_high = ci[2],
          Status = "Observed in both campaigns"
        )
      } else {
        tibble(
          Reach_code = reach, Reach = REACH_LABELS[[reach]], Pathway = pathway,
          Isotope = isotope_label,
          Wet_n = length(wet), Dry_n = length(dry),
          Wet_mean = ifelse(length(wet) > 0, mean(wet), NA_real_),
          Dry_mean = ifelse(length(dry) > 0, mean(dry), NA_real_),
          Delta_dry_minus_wet = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
          Status = "Not estimable: pathway missing from one campaign"
        )
      }
    })
  })
})

# ----------------------- Baseline-relative fish d13C --------------------------
baseline_sensitivity_d13C <- map_dfr(REACH_ORDER, function(reach) {
  fish_wet <- fish %>% filter(Reach_code == reach, Campaign == "Wet_2023") %>% pull(d13C_permil)
  fish_dry <- fish %>% filter(Reach_code == reach, Campaign == "Dry_2023") %>% pull(d13C_permil)
  fish_wet_draws <- bootstrap_means(fish_wet)
  fish_dry_draws <- bootstrap_means(fish_dry)
  fish_delta_draws <- fish_dry_draws - fish_wet_draws
  fish_delta <- mean(fish_dry) - mean(fish_wet)
  raw_ci <- bootstrap_ci(fish_delta_draws)

  raw_row <- tibble(
    Reach_code = reach, Reach = REACH_LABELS[[reach]], Baseline_pathway = "Raw fish",
    Fish_wet_n = length(fish_wet), Fish_dry_n = length(fish_dry),
    Fish_wet_mean = mean(fish_wet), Fish_dry_mean = mean(fish_dry), Fish_delta = fish_delta,
    Baseline_wet_n = NA_integer_, Baseline_dry_n = NA_integer_,
    Baseline_wet_mean = NA_real_, Baseline_dry_mean = NA_real_, Baseline_delta = NA_real_,
    Relative_delta = fish_delta, CI_low = raw_ci[1], CI_high = raw_ci[2],
    Direction_vs_raw = "Reference", Status = "Unadjusted fish isotope-space change"
  )

  adjusted_rows <- map_dfr(PATHWAY_ORDER, function(pathway) {
    selected <- resources %>% filter(Reach_code == reach) %>% select_pathway(pathway)
    base_wet <- selected %>% filter(Campaign == "Wet_2023") %>% pull(d13C_permil) %>% clean_numeric()
    base_dry <- selected %>% filter(Campaign == "Dry_2023") %>% pull(d13C_permil) %>% clean_numeric()

    if (length(base_wet) > 0 && length(base_dry) > 0) {
      base_delta_draws <- bootstrap_means(base_dry) - bootstrap_means(base_wet)
      relative_draws <- fish_delta_draws - base_delta_draws
      relative_delta <- fish_delta - (mean(base_dry) - mean(base_wet))
      ci <- bootstrap_ci(relative_draws)
      direction <- ifelse(sign(relative_delta) == sign(fish_delta), "Same direction", "Direction reversed")
      tibble(
        Reach_code = reach, Reach = REACH_LABELS[[reach]], Baseline_pathway = pathway,
        Fish_wet_n = length(fish_wet), Fish_dry_n = length(fish_dry),
        Fish_wet_mean = mean(fish_wet), Fish_dry_mean = mean(fish_dry), Fish_delta = fish_delta,
        Baseline_wet_n = length(base_wet), Baseline_dry_n = length(base_dry),
        Baseline_wet_mean = mean(base_wet), Baseline_dry_mean = mean(base_dry),
        Baseline_delta = mean(base_dry) - mean(base_wet),
        Relative_delta = relative_delta, CI_low = ci[1], CI_high = ci[2],
        Direction_vs_raw = direction, Status = "Observed baseline in both campaigns"
      )
    } else {
      missing_note <- ifelse(length(base_wet) == 0, "Not estimable: wet baseline absent",
                             "Not estimable: dry baseline absent")
      tibble(
        Reach_code = reach, Reach = REACH_LABELS[[reach]], Baseline_pathway = pathway,
        Fish_wet_n = length(fish_wet), Fish_dry_n = length(fish_dry),
        Fish_wet_mean = mean(fish_wet), Fish_dry_mean = mean(fish_dry), Fish_delta = fish_delta,
        Baseline_wet_n = length(base_wet), Baseline_dry_n = length(base_dry),
        Baseline_wet_mean = ifelse(length(base_wet) > 0, mean(base_wet), NA_real_),
        Baseline_dry_mean = ifelse(length(base_dry) > 0, mean(base_dry), NA_real_),
        Baseline_delta = NA_real_, Relative_delta = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
        Direction_vs_raw = "Not estimable", Status = missing_note
      )
    }
  })
  bind_rows(raw_row, adjusted_rows)
})

# -------------------------- Trophic-position analysis -------------------------
target_gastropod_d15N <- function(reach, campaign) {
  resources %>%
    filter(
      Reach_code == reach, Campaign == campaign,
      Resource_type_original == "Rusty-patina Ring-ribbed Snail"
    ) %>%
    pull(d15N_permil) %>%
    clean_numeric()
}

TDF_VALUES <- c(2.9, 3.4, 3.9)
tp_store <- list()
tp_group_rows <- list()

for (reach in REACH_ORDER) {
  for (campaign in CAMPAIGNS) {
    fish_values <- fish %>% filter(Reach_code == reach, Campaign == campaign) %>% pull(d15N_permil)
    baseline_campaign <- campaign
    baseline_values <- target_gastropod_d15N(reach, baseline_campaign)
    if (length(baseline_values) == 0) {
      stop("Target gastropod baseline absent for ", reach, " ", campaign)
    }
    baseline_use <- "Observed"
    fish_draws <- bootstrap_means(fish_values)
    baseline_draws <- bootstrap_means(baseline_values)
    tp_draws <- 2 + (fish_draws - baseline_draws) / 3.4
    ci <- bootstrap_ci(tp_draws)

    tp_group_rows[[paste(reach, campaign, sep = "|")]] <- tibble(
      Reach_code = reach, Reach = REACH_LABELS[[reach]], Campaign = campaign,
      Fish_n = length(fish_values), Fish_mean_d15N = mean(fish_values),
      Baseline_campaign = baseline_campaign, Baseline_n = length(baseline_values),
      Baseline_mean_d15N = mean(baseline_values), Baseline_use = baseline_use,
      TP_2.9 = 2 + (mean(fish_values) - mean(baseline_values)) / 2.9,
      TP_3.4 = 2 + (mean(fish_values) - mean(baseline_values)) / 3.4,
      TP_3.9 = 2 + (mean(fish_values) - mean(baseline_values)) / 3.9,
      TP_3.4_CI_low = ci[1], TP_3.4_CI_high = ci[2]
    )
    tp_store[[paste(reach, campaign, sep = "|")]] <- list(
      fish = fish_values, baseline = baseline_values,
      fish_draws = fish_draws, baseline_draws = baseline_draws
    )
  }
}
tp_groups <- bind_rows(tp_group_rows)

tp_changes <- map_dfr(REACH_ORDER, function(reach) {
  wet <- tp_store[[paste(reach, "Wet_2023", sep = "|")]]
  dry <- tp_store[[paste(reach, "Dry_2023", sep = "|")]]
  fish_delta_draws <- dry$fish_draws - wet$fish_draws
  fish_delta <- mean(dry$fish) - mean(wet$fish)
  baseline_delta_draws <- dry$baseline_draws - wet$baseline_draws
  baseline_delta <- mean(dry$baseline) - mean(wet$baseline)
  status <- "Primary: observed site- and campaign-specific baselines"
  corrected_draws <- fish_delta_draws - baseline_delta_draws
  tp_change_draws <- corrected_draws / 3.4
  ci <- bootstrap_ci(tp_change_draws)
  corrected_delta <- fish_delta - baseline_delta
  tibble(
    Reach_code = reach, Reach = REACH_LABELS[[reach]],
    Fish_delta_d15N = fish_delta, Baseline_delta_d15N = baseline_delta,
    Corrected_delta_d15N = corrected_delta,
    Delta_TP_2.9 = corrected_delta / 2.9,
    Delta_TP_3.4 = corrected_delta / 3.4,
    Delta_TP_3.9 = corrected_delta / 3.9,
    Delta_TP_3.4_CI_low = ci[1], Delta_TP_3.4_CI_high = ci[2],
    Status = status
  )
})

# Individual gastropod-referenced coordinates. These are baseline-referenced
# isotope-space values and are not a carbon-source mixing model.
gastropod_baselines <- resources %>%
  filter(Resource_type_original == "Rusty-patina Ring-ribbed Snail") %>%
  group_by(Reach_code, Campaign) %>%
  summarise(
    Gastropod_n = n(),
    Gastropod_mean_d13C = mean(d13C_permil),
    Gastropod_mean_d15N = mean(d15N_permil),
    .groups = "drop"
  )

baseline_normalized_fish <- fish %>%
  select(
    Fish_isotope_record_id, Reach_code, Reach, Campaign, Scientific_name,
    Feeding_guild, d13C_permil, d15N_permil
  ) %>%
  left_join(gastropod_baselines, by = c("Reach_code", "Campaign")) %>%
  mutate(
    d13C_relative_to_gastropod = d13C_permil - Gastropod_mean_d13C,
    TP_TDF_3.4 = 2 + (d15N_permil - Gastropod_mean_d15N) / 3.4
  )

if (any(!is.finite(baseline_normalized_fish$d13C_relative_to_gastropod)) ||
    any(!is.finite(baseline_normalized_fish$TP_TDF_3.4))) {
  stop("At least one fish group lacks an observed target-gastropod baseline.")
}

# ---------------------------- Size sensitivity --------------------------------
fish_by_species <- fish %>%
  group_by(Reach_code, Campaign, Scientific_name) %>%
  summarise(
    Isotope_n = n(), Mean_d13C = mean(d13C_permil), Mean_d15N = mean(d15N_permil),
    Feeding_guild = first(na.omit(Feeding_guild)), .groups = "drop"
  )

morph_by_species <- morph %>%
  filter(Principal_isotope_analysis == "Yes") %>%
  group_by(Reach_code, Campaign, Scientific_name) %>%
  summarise(
    Morphometric_n = n(),
    Mean_total_length_mm = mean(Total_length_mm, na.rm = TRUE),
    Mean_body_length_mm = mean(Body_length_mm, na.rm = TRUE),
    Mean_body_weight_g = mean(Body_weight_g, na.rm = TRUE),
    Mean_log_weight = mean(log(Body_weight_g), na.rm = TRUE),
    .groups = "drop"
  )

combined_species <- inner_join(
  fish_by_species, morph_by_species,
  by = c("Reach_code", "Campaign", "Scientific_name")
)

wet_species <- combined_species %>%
  filter(Campaign == "Wet_2023") %>%
  select(-Campaign) %>%
  rename_with(~ paste0(.x, "_Wet"), -c(Reach_code, Scientific_name))

dry_species <- combined_species %>%
  filter(Campaign == "Dry_2023") %>%
  select(-Campaign) %>%
  rename_with(~ paste0(.x, "_Dry"), -c(Reach_code, Scientific_name))

size_pairs <- inner_join(wet_species, dry_species, by = c("Reach_code", "Scientific_name")) %>%
  mutate(
    Counts_match_both = Isotope_n_Wet == Morphometric_n_Wet &
      Isotope_n_Dry == Morphometric_n_Dry,
    Minimum_n = pmin(Isotope_n_Wet, Isotope_n_Dry, Morphometric_n_Wet, Morphometric_n_Dry),
    Reach = REACH_LABELS[Reach_code],
    Feeding_guild = coalesce(Feeding_guild_Wet, Feeding_guild_Dry),
    Delta_d13C = Mean_d13C_Dry - Mean_d13C_Wet,
    Delta_d15N = Mean_d15N_Dry - Mean_d15N_Wet,
    Delta_mean_total_length_mm = Mean_total_length_mm_Dry - Mean_total_length_mm_Wet,
    Delta_log_mean_total_length = log(Mean_total_length_mm_Dry) - log(Mean_total_length_mm_Wet),
    Delta_mean_log_weight = Mean_log_weight_Dry - Mean_log_weight_Wet
  ) %>%
  filter(Counts_match_both, Minimum_n >= 3) %>%
  select(
    Reach_code, Reach, Scientific_name, Feeding_guild,
    Isotope_n_Wet, Isotope_n_Dry, Morphometric_n_Wet, Morphometric_n_Dry, Minimum_n,
    Mean_d13C_Wet, Mean_d13C_Dry, Delta_d13C,
    Mean_d15N_Wet, Mean_d15N_Dry, Delta_d15N,
    Mean_total_length_mm_Wet, Mean_total_length_mm_Dry, Delta_mean_total_length_mm,
    Delta_log_mean_total_length, Mean_body_weight_g_Wet, Mean_body_weight_g_Dry,
    Delta_mean_log_weight
  ) %>%
  arrange(Reach_code, Scientific_name)

correlation_table <- crossing(
  Size_predictor = c("Delta_log_mean_total_length", "Delta_mean_log_weight"),
  Isotope_change = c("Delta_d13C", "Delta_d15N")
) %>%
  pmap_dfr(function(Size_predictor, Isotope_change) {
    x <- size_pairs[[Size_predictor]]
    y <- size_pairs[[Isotope_change]]
    keep <- is.finite(x) & is.finite(y)
    pearson <- cor.test(x[keep], y[keep], method = "pearson")
    spearman <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
    tibble(
      Size_predictor = Size_predictor, Isotope_change = Isotope_change,
      n = sum(keep), Pearson_r = unname(pearson$estimate), Pearson_p = pearson$p.value,
      Spearman_rho = unname(spearman$estimate), Spearman_p = spearman$p.value
    )
  })

fit_cluster_model <- function(data, outcome, predictor) {
  model_data <- data %>%
    select(all_of(c(outcome, predictor, "Reach_code", "Feeding_guild", "Scientific_name"))) %>%
    drop_na()
  formula <- reformulate(c(predictor, "Reach_code", "Feeding_guild"), response = outcome)
  fit <- lm(formula, data = model_data)
  vcov_cluster <- sandwich::vcovCL(
    fit, cluster = model_data$Scientific_name, type = "HC1", cadjust = TRUE
  )
  estimate <- coef(fit)
  se <- sqrt(diag(vcov_cluster))
  clusters <- n_distinct(model_data$Scientific_name)
  df_cluster <- clusters - 1
  t_value <- estimate / se
  p_value <- 2 * pt(abs(t_value), df = df_cluster, lower.tail = FALSE)
  critical <- qt(0.975, df = df_cluster)
  tibble(
    Outcome = outcome, Size_predictor = predictor, Term = names(estimate),
    Estimate = unname(estimate), Cluster_robust_SE = unname(se),
    t = unname(t_value), df = df_cluster, p = unname(p_value),
    CI_low = unname(estimate - critical * se), CI_high = unname(estimate + critical * se),
    n_pairs = nrow(model_data), n_species_clusters = clusters,
    Model = paste0(outcome, " ~ ", predictor,
                   " + Reach_code + Feeding_guild; SE clustered by species")
  )
}

size_models <- crossing(
  Size_predictor = c("Delta_log_mean_total_length", "Delta_mean_log_weight"),
  Outcome = c("Delta_d13C", "Delta_d15N")
) %>%
  pmap_dfr(function(Size_predictor, Outcome) {
    fit_cluster_model(size_pairs, Outcome, Size_predictor)
  })

# ----------------------------- Save result tables -----------------------------
result_tables <- list(
  fish_group_summary = fish_summary,
  resource_group_summary = resource_group_summary,
  resource_seasonal_shifts = resource_seasonal_shifts,
  baseline_sensitivity_d13C = baseline_sensitivity_d13C,
  trophic_position_groups = tp_groups,
  trophic_position_changes = tp_changes,
  baseline_normalized_fish = baseline_normalized_fish,
  resource_campaign_audit = campaign_audit,
  size_strict_pairs = size_pairs,
  size_correlations = correlation_table,
  size_cluster_robust_models = size_models
)

iwalk(result_tables, function(data, name) {
  write_excel_csv(data, file.path(OUTPUT_DIR, paste0(name, ".csv")))
})

# R-generated results workbook
wb <- createWorkbook(creator = "Reproducible R analysis")
header_style <- createStyle(
  fgFill = "#1F6F78", fontColour = "#FFFFFF", textDecoration = "bold",
  halign = "center", valign = "center", wrapText = TRUE
)
title_style <- createStyle(
  fgFill = "#173F4F", fontColour = "#FFFFFF", textDecoration = "bold",
  fontSize = 14, valign = "center"
)

addWorksheet(wb, "README", gridLines = FALSE)
writeData(wb, "README", "RStudio isotope-space reorganization analysis", startRow = 1, startCol = 1)
addStyle(wb, "README", title_style, rows = 1, cols = 1:4, gridExpand = TRUE)
readme_rows <- data.frame(
  Item = c("Input file", "Fish isotope records", "2023 resource records", "Strict size pairs",
           "Bootstrap draws", "Random seed", "Primary terminology", "Trophic inference",
           "Baseline taxon", "Campaign recoding"),
  Value = c(
    basename(INPUT_FILE), nrow(fish), nrow(resources), nrow(size_pairs),
    BOOTSTRAP_DRAWS, RANDOM_SEED,
    "Use isotope-space reorganization for d13C displacement and ellipse geometry",
    "Restrict to baseline-corrected d15N trophic position",
    "Rusty-patina ring-ribbed snail (Sinotaia quadrata; synonym Bellamya aeruginosa)",
    paste(nrow(campaign_audit), "Hanchuan gastropod rows recoded to Wet_2023 from source code 2306")
  )
)
writeDataTable(wb, "README", readme_rows, startRow = 3, tableStyle = "TableStyleMedium2")
setColWidths(wb, "README", cols = 1:2, widths = c(28, 90))

sheet_names <- c(
  "Fish Summary", "Resource Summary", "Resource Shifts", "Baseline Sensitivity",
  "TP Groups", "TP Changes", "Normalized Fish", "Campaign Audit",
  "Size Pairs", "Size Correlations", "Size Models"
)
for (i in seq_along(result_tables)) {
  sheet <- sheet_names[i]
  data <- result_tables[[i]]
  addWorksheet(wb, sheet, gridLines = FALSE)
  writeDataTable(wb, sheet, data, startRow = 1, tableStyle = "TableStyleMedium2")
  addStyle(wb, sheet, header_style, rows = 1, cols = seq_len(ncol(data)), gridExpand = TRUE)
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(wb, sheet, cols = seq_len(ncol(data)), widths = "auto")
}
saveWorkbook(
  wb, file.path(OUTPUT_DIR, "Isotope_Reorganization_Analysis_from_R.xlsx"),
  overwrite = TRUE
)

# ------------------------------ Publication figures ---------------------------
available_combos <- resource_group_summary %>%
  filter(n > 0) %>%
  distinct(Reach_code, Pathway, Campaign) %>%
  count(Reach_code, Pathway, name = "campaign_n") %>%
  filter(campaign_n == 2) %>%
  arrange(factor(Reach_code, REACH_ORDER), factor(Pathway, PATHWAY_ORDER))

combo_levels <- paste(available_combos$Reach_code, available_combos$Pathway, sep = "\n")
figure1_data <- resource_group_summary %>%
  semi_join(available_combos, by = c("Reach_code", "Pathway")) %>%
  mutate(
    Campaign = factor(Campaign, levels = CAMPAIGNS, labels = c("Wet", "Dry")),
    Combo = factor(paste(Reach_code, Pathway, sep = "\n"), levels = combo_levels),
    Isotope = factor(Isotope, levels = c("d13C", "d15N"),
                     labels = c("Basal-resource δ13C (‰)", "Basal-resource δ15N (‰)"))
  )

pd <- position_dodge(width = 0.36)
figure1 <- ggplot(figure1_data, aes(Combo, mean, color = Campaign)) +
  geom_line(aes(group = interaction(Reach_code, Pathway)), color = "#A0A0A0",
            linewidth = 0.45, position = pd) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.10,
                linewidth = 0.45, position = pd) +
  geom_point(size = 2.4, position = pd) +
  facet_wrap(~ Isotope, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c(Wet = SEASON_COLORS[["Wet_2023"]],
                                Dry = SEASON_COLORS[["Dry_2023"]])) +
  labs(
    title = "Seasonal movement of recurrent basal-resource pathways",
    x = NULL, y = NULL
  ) +
  theme_submission(11) +
  theme(axis.text.x = element_text(size = 10.5), legend.position = "top")
save_submission_plot(figure1, "Figure_1_Basal_resource_shifts", 7.2, 7.0)

figure2_data <- baseline_sensitivity_d13C %>%
  filter(is.finite(Relative_delta)) %>%
  mutate(
    Reach = factor(Reach, levels = REACH_LABELS[REACH_ORDER]),
    Baseline_pathway = factor(Baseline_pathway,
                              levels = rev(c("Raw fish", "Phytoplankton", "POM", "Gastropod")))
  )

figure2 <- ggplot(figure2_data, aes(Relative_delta, Baseline_pathway, color = Baseline_pathway)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, color = "#555555") +
  geom_segment(aes(x = CI_low, xend = CI_high, yend = Baseline_pathway), linewidth = 0.55) +
  geom_point(size = 2.6) +
  facet_wrap(~ Reach, nrow = 1) +
  scale_color_manual(values = PATHWAY_COLORS, guide = "none") +
  labs(
    title = "Carbon-shift inference depends on the basal pathway",
    x = "Dry–wet fish δ13C change relative to baseline (‰)", y = "Carbon reference"
  ) +
  theme_submission(11)
save_submission_plot(figure2, "Figure_2_Baseline_relative_d13C", 7.4, 3.9)

figure3_data <- tp_groups %>%
  mutate(
    Reach = factor(Reach, levels = REACH_LABELS[REACH_ORDER]),
    Campaign = factor(Campaign, levels = CAMPAIGNS, labels = c("Wet", "Dry"))
  )
pd_tp <- position_dodge(width = 0.35)
figure3 <- ggplot(figure3_data, aes(Reach, TP_3.4, color = Campaign)) +
  geom_errorbar(aes(ymin = TP_3.4_CI_low, ymax = TP_3.4_CI_high),
                width = 0.10, linewidth = 0.55, position = pd_tp) +
  geom_point(size = 2.8, position = pd_tp) +
  scale_color_manual(values = c(Wet = SEASON_COLORS[["Wet_2023"]],
                                Dry = SEASON_COLORS[["Dry_2023"]])) +
  labs(
    title = "Baseline correction reveals reach-specific trophic-position changes",
    x = NULL, y = "Baseline-corrected trophic position"
  ) +
  theme_submission(11) +
  theme(legend.position = "bottom")
save_submission_plot(figure3, "Figure_3_Trophic_position", 6.2, 4.4)

make_size_panel <- function(outcome, y_label, panel_label) {
  rho_test <- suppressWarnings(cor.test(
    size_pairs$Delta_log_mean_total_length, size_pairs[[outcome]],
    method = "spearman", exact = FALSE
  ))
  ggplot(size_pairs, aes(Delta_log_mean_total_length, .data[[outcome]], color = Reach_code)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "#777777") +
    geom_smooth(method = "lm", se = TRUE, color = "#333333", fill = "#BDBDBD",
                alpha = 0.25, linewidth = 0.6) +
    geom_point(size = 2.2, alpha = 0.88) +
    scale_color_manual(values = REACH_COLORS, labels = REACH_LABELS) +
    annotate(
      "text", x = -Inf, y = Inf,
      label = sprintf("Spearman ρ = %.2f\np = %.3f", unname(rho_test$estimate), rho_test$p.value),
      hjust = -0.05, vjust = 1.15, size = 3.6
    ) +
    labs(
      title = panel_label,
      x = "Dry–wet change in log mean total length", y = y_label
    ) +
    theme_submission(11)
}

size_panel_1 <- make_size_panel("Delta_d13C", "Dry–wet mean δ13C (‰)", "A")
size_panel_2 <- make_size_panel("Delta_d15N", "Dry–wet mean δ15N (‰)", "B") +
  theme(legend.position = "none")
figure4 <- (size_panel_1 + size_panel_2) +
  plot_annotation(title = "Species-level isotope shifts show limited size dependence") &
  theme(plot.title = element_text(face = "bold", size = 13))
save_submission_plot(figure4, "Figure_4_Size_sensitivity", 7.2, 3.7)

figure5_data <- baseline_normalized_fish %>%
  mutate(
    Reach = factor(Reach, levels = REACH_LABELS[REACH_ORDER]),
    Campaign = factor(Campaign, levels = CAMPAIGNS, labels = c("Wet", "Dry"))
  )

figure5 <- ggplot(
  figure5_data,
  aes(d13C_relative_to_gastropod, TP_TDF_3.4, color = Campaign)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "#777777") +
  geom_point(size = 0.8, alpha = 0.22) +
  stat_ellipse(type = "norm", level = 0.95, linewidth = 0.7) +
  stat_summary(fun = mean, geom = "point", size = 2.8, shape = 21, fill = "white") +
  facet_wrap(~ Reach, nrow = 1) +
  scale_color_manual(values = c(Wet = SEASON_COLORS[["Wet_2023"]],
                                Dry = SEASON_COLORS[["Dry_2023"]])) +
  labs(
    title = "Gastropod-referenced fish isotope space",
    x = "Fish δ13C relative to target gastropod (‰)",
    y = "Gastropod-baseline trophic position (TDF = 3.4‰)"
  ) +
  theme_submission(11) +
  theme(legend.position = "top")
save_submission_plot(figure5, "Figure_5_Gastropod_referenced_isotope_space", 8.2, 3.8)

# ------------------------------ Captions and log ------------------------------
captions <- c(
  "Figure 1. Seasonal movement of recurrent basal-resource pathways. Points show mean stable-isotope values and error bars show nonparametric 95% bootstrap confidence intervals. Lines connect wet- and dry-season observations within the same reach and basal pathway. Only pathways observed in both 2023 campaigns are displayed.",
  "",
  "Figure 2. Sensitivity of community fish carbon-isotope change to baseline selection. Points show dry-minus-wet change in mean fish δ13C before adjustment or after subtracting the seasonal change in phytoplankton, POM, or target-gastropod δ13C. Error bars are 95% bootstrap confidence intervals that propagate uncertainty in fish and basal-resource means.",
  "",
  "Figure 3. Baseline-corrected trophic position by reach and campaign. Trophic position was calculated as TP = 2 + (δ15Nfish − δ15Nbaseline)/3.4 using observed site- and campaign-specific target-gastropod baselines. Error bars are 95% bootstrap confidence intervals.",
  "",
  "Figure 4. Species-level sensitivity of isotope shifts to changes in size structure. Each point represents a reach–species pair with matching isotope and morphometric sample counts and at least three observations in both campaigns. Lines show unadjusted ordinary least-squares relationships with 95% confidence bands. Formal models adjust for reach and feeding guild and use species-clustered standard errors.",
  "",
  "Figure 5. Gastropod-referenced fish isotope space. Individual fish are shown after subtracting the site- and campaign-specific mean target-gastropod δ13C and converting δ15N to trophic position with a 3.4‰ trophic discrimination factor. Large points are group centroids and outlines are normal-theory 95% data ellipses. This is a baseline-referenced isotope-space display, not a carbon-source mixing model."
)
writeLines(captions, file.path(OUTPUT_DIR, "Figure_captions.txt"), useBytes = TRUE)

run_summary <- c(
  paste("Input file:", normalizePath(INPUT_FILE)),
  paste("Fish isotope records:", nrow(fish)),
  paste("Matched 2023 resource records:", nrow(resources)),
  paste("Campaign labels standardized:", nrow(campaign_audit)),
  paste("Strict size pairs:", nrow(size_pairs)),
  paste("Bootstrap draws:", BOOTSTRAP_DRAWS),
  paste("Random seed:", RANDOM_SEED),
  "Site mapping: S1 = Hanchuan; S2 = Qianjiang; S3 = Shayang",
  "Figure axis titles: 13 pt bold; axis tick labels: 11 pt",
  "",
  "Primary interpretation:",
  "Use isotope-space reorganization for d13C displacement and ellipse geometry.",
  "Restrict trophic inference to baseline-corrected d15N trophic position.",
  "All three reaches use observed wet and dry target-gastropod baselines."
)
writeLines(run_summary, file.path(OUTPUT_DIR, "Run_Summary.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(OUTPUT_DIR, "R_sessionInfo.txt"))

message("Analysis complete. Results written to: ", normalizePath(OUTPUT_DIR))

# FIGURE 8 REVIEWER REVISION 2026-11-01 | R1-R7 / A1-A4
# =============================================================================
# LOWER HAN RIVER FISH STUDY
# ENVIRONMENTAL INTEGRATION, RESOURCE BASELINES, AND REVISED TROPHIC POSITION
# MASTER-DATA-ALIGNED, DO-INCLUSIVE FINAL R SCRIPT — ALL 609 FISH
#
# This script uses the authoritative matched campaign windows in the master
# dataset:
#   20–27 June 2023     = Wet_2023
#   14–23 November 2023 = Dry_2023
#
# PRIMARY MODULES
#   1. Environmental variable transformation and correlation screening
#   2. Transect-level environmental PCA and standardized heat map
#   3. Environmental PERMANOVA and multivariate dispersion tests
#   4. Univariate Reach × Campaign environmental models
#   5. Resource-isotope summaries and isotope-space visualization
#   6. Revised community-mean trophic position using measured primary-consumer
#      baselines, observed-mean points, fixed TDF and bootstrap uncertainty
#   7. Exact matrix tests linking environmental dissimilarity to:
#        - Bray-Curtis fish-community dissimilarity
#        - Jaccard fish-community dissimilarity
#        - isotope-centroid distance
#        - standard-ellipse-area difference
#        - proportional ellipse-area difference
#        - isotope-niche non-overlap
#        - revised trophic-position difference
#   8. Procrustes concordance among environmental, community, and isotope
#      configurations
#
# ANALYTICAL STATUS
#   - Environmental PCA/PERMANOVA: non-Bayesian.
#   - Environmental–biological Mantel tests: exact permutation analyses.
#   - Revised trophic position: nonparametric bootstrap at fixed TDF, not Bayesian.
#
# IMPORTANT INFERENCE LIMITS
#   - Environmental analyses use 18 transect observations.
#   - Biological correspondence analyses use six aggregated reach–campaign
#     groups, so exact matrix tests are exploratory despite complete
#     exploratory effect estimates and two restricted-null sensitivities.
#     Three paired reaches do not establish indicator complementarity.
#   - Fish-level environmental regressions are intentionally not performed.
#   - All six reach–campaign trophic baselines use observed primary-consumer
#     measurements from the same reach and campaign. HC wet has three
#     measured wet-campaign S. quadrata; cross-season proxies are prohibited.
#   - TP interpretation is conditional on pathway representation and uncalibrated
#     tissue integration times. Bulk-zooplankton comparisons are exploratory.
#   - June 2024 resource-isotope records are retained in the input workbook for
#     provenance but excluded from matched 2023 resource summaries and figures.
#
# INPUT WORKBOOK
#   Han_River_Figure8_Analysis_Ready.xlsx
#
# HOW TO RUN
#   1. Put this R script and the Excel workbook in the same folder.
#   2. Restart RStudio.
#   3. Open this script and click Source once.
#   4. Select the workbook if R cannot locate it automatically.
# =============================================================================


# =============================================================================
# 1. USER SETTINGS
# =============================================================================

SEED <- 20260731L

# Publication labels. HC, QJ, and SY remain internal analysis codes.
REACH_LABELS <- c(
  HC = "S1",
  QJ = "S2",
  SY = "S3"
)

REACH_LABEL_LEVELS <- unname(
  REACH_LABELS[c("HC", "QJ", "SY")]
)

# Publication typography in points.
AXIS_TITLE_SIZE <- 14
AXIS_TEXT_SIZE <- 12
FACET_TEXT_SIZE <- 12
LEGEND_TEXT_SIZE <- 11
LEGEND_TITLE_SIZE <- 11

N_PERMUTATIONS_RANDOM <- 9999L
N_TP_DRAWS <- 5000L

TDF_MEAN <- 3.4  # Fixed in the primary analysis.
FIXED_TDF_VALUES <- c(2.9, 3.4, 3.9)
TROPHIC_LEVEL_BASELINE <- 2

CORRELATION_CUTOFF <- 0.90

# DO is scientifically central to this analysis. If correlation screening
# identifies a redundant pair (for example, DO versus water temperature), DO
# is retained and the later-listed correlated variable is removed.
CORRELATION_PRIORITY_VARIABLES <- c(
  "DO_mg_L"
)

MATCHED_RESOURCE_CAMPAIGNS <- c(
  "Wet_2023",
  "Dry_2023"
)

PRIMARY_ENVIRONMENTAL_VARIABLES <- c(
  "Depth_m",
  "Transparency_m",
  "Flow_velocity_m_s",
  "Water_temperature_C",
  "DO_mg_L",
  "Conductivity_uS_cm",
  "pH",
  "Water_TN_mg_L",
  "Water_TP_mg_L",
  "NH4_N_mg_L",
  "NO3_N_mg_L",
  "Phosphate_P_mg_L",
  "Chlorophyll_a_ug_L",
  "TSS_mg_L",
  "COD_mg_L"
)

LOG1P_VARIABLES <- c(
  "Flow_velocity_m_s",
  "Water_TN_mg_L",
  "Water_TP_mg_L",
  "NH4_N_mg_L",
  "NO3_N_mg_L",
  "Phosphate_P_mg_L",
  "Chlorophyll_a_ug_L",
  "TSS_mg_L",
  "COD_mg_L"
)

ENVIRONMENTAL_DOMAINS <- list(
  Hydrology_Habitat = c(
    "Depth_m",
    "Transparency_m",
    "Flow_velocity_m_s"
  ),

  Thermal_Oxygen = c(
    "Water_temperature_C",
    "DO_mg_L",
    "Conductivity_uS_cm",
    "pH"
  ),

  Nutrient_Productivity = c(
    "Water_TN_mg_L",
    "Water_TP_mg_L",
    "NH4_N_mg_L",
    "NO3_N_mg_L",
    "Phosphate_P_mg_L",
    "Chlorophyll_a_ug_L",
    "TSS_mg_L",
    "COD_mg_L"
  )
)

SELECTED_ENVIRONMENTAL_PLOT_VARIABLES <- c(
  "Water_temperature_C",
  "DO_mg_L",
  "Flow_velocity_m_s",
  "Water_TN_mg_L",
  "Water_TP_mg_L",
  "Chlorophyll_a_ug_L"
)

set.seed(SEED)
graphics.off()


# =============================================================================
# 2. PACKAGE CHECK
# =============================================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggrepel",
  "vegan",
  "permute",
  "pheatmap",
  "car",
  "writexl"
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
      "\nInstall them with:\ninstall.packages(c(",
      paste(
        paste0('"', missing_packages, '"'),
        collapse = ", "
      ),
      "))\nThen restart RStudio and run the script again."
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(vegan)
  library(permute)
  library(pheatmap)
  library(car)
  library(writexl)
})


# =============================================================================
# 3. GENERAL HELPERS
# =============================================================================

safe_numeric <- function(
  x,
  variable_name
) {

  output <- suppressWarnings(
    as.numeric(
      x
    )
  )

  if (anyNA(
    output
  )) {
    stop(
      paste0(
        "Variable ",
        variable_name,
        " contains missing or non-numeric values."
      ),
      call. = FALSE
    )
  }

  output
}


standard_error <- function(
  x
) {

  x <- x[
    is.finite(
      x
    )
  ]

  if (length(
    x
  ) < 2L) {
    return(
      NA_real_
    )
  }

  stats::sd(
    x
  ) /
    sqrt(
      length(
        x
      )
    )
}


quantile_025 <- function(
  x
) {

  as.numeric(
    stats::quantile(
      x,
      probs = 0.025,
      na.rm = TRUE,
      names = FALSE
    )
  )
}


quantile_975 <- function(
  x
) {

  as.numeric(
    stats::quantile(
      x,
      probs = 0.975,
      na.rm = TRUE,
      names = FALSE
    )
  )
}


screen_correlated_variables <- function(
  data,
  candidate_variables,
  cutoff = 0.90,
  priority_variables = character(0)
) {

  candidate_variables <- candidate_variables[
    candidate_variables %in% names(
      data
    )
  ]

  priority_variables <- priority_variables[
    priority_variables %in% candidate_variables
  ]

  candidate_variables <- unique(
    c(
      priority_variables,
      candidate_variables
    )
  )

  zero_variance <- candidate_variables[
    vapply(
      candidate_variables,
      function(variable_name) {

        current_values <- data[[
          variable_name
        ]]

        stats::sd(
          current_values,
          na.rm = TRUE
        ) == 0
      },
      FUN.VALUE = logical(1)
    )
  ]

  remaining_variables <- setdiff(
    candidate_variables,
    zero_variance
  )

  retained_variables <- character(0)

  screening_rows <- list()
  screening_index <- 1L

  for (variable_name in remaining_variables) {

    if (length(
      retained_variables
    ) == 0L) {

      retained_variables <- c(
        retained_variables,
        variable_name
      )

      screening_rows[[
        screening_index
      ]] <- data.frame(
        Variable = variable_name,
        Status = "Retained",
        Most_correlated_retained_variable = NA_character_,
        Absolute_correlation = NA_real_,
        Reason = "First retained variable",
        stringsAsFactors = FALSE
      )

      screening_index <- screening_index + 1L
      next
    }

    current_correlations <- vapply(
      retained_variables,
      function(retained_name) {

        abs(
          stats::cor(
            data[[
              variable_name
            ]],
            data[[
              retained_name
            ]],
            use = "pairwise.complete.obs"
          )
        )
      },
      FUN.VALUE = numeric(1)
    )

    strongest_index <- which.max(
      current_correlations
    )

    strongest_correlation <- current_correlations[
      strongest_index
    ]

    strongest_variable <- retained_variables[
      strongest_index
    ]

    if (
      is.finite(
        strongest_correlation
      ) &&
      strongest_correlation > cutoff
    ) {

      screening_rows[[
        screening_index
      ]] <- data.frame(
        Variable = variable_name,
        Status = "Removed",
        Most_correlated_retained_variable = strongest_variable,
        Absolute_correlation = strongest_correlation,
        Reason = paste0(
          "|r| > ",
          cutoff
        ),
        stringsAsFactors = FALSE
      )

    } else {

      retained_variables <- c(
        retained_variables,
        variable_name
      )

      screening_rows[[
        screening_index
      ]] <- data.frame(
        Variable = variable_name,
        Status = "Retained",
        Most_correlated_retained_variable = strongest_variable,
        Absolute_correlation = strongest_correlation,
        Reason = paste0(
          "|r| <= ",
          cutoff
        ),
        stringsAsFactors = FALSE
      )
    }

    screening_index <- screening_index + 1L
  }

  if (length(
    zero_variance
  ) > 0L) {

    zero_rows <- lapply(
      zero_variance,
      function(variable_name) {

        data.frame(
          Variable = variable_name,
          Status = "Removed",
          Most_correlated_retained_variable = NA_character_,
          Absolute_correlation = NA_real_,
          Reason = "Zero variance",
          stringsAsFactors = FALSE
        )
      }
    )

    screening_rows <- c(
      screening_rows,
      zero_rows
    )
  }

  list(
    retained = retained_variables,
    screening = dplyr::bind_rows(
      screening_rows
    )
  )
}


all_permutations <- function(
  values
) {

  if (length(
    values
  ) == 1L) {

    return(
      matrix(
        values,
        nrow = 1L
      )
    )
  }

  permutation_blocks <- lapply(
    seq_along(
      values
    ),
    function(index) {

      remaining_permutations <- all_permutations(
        values[
          -index
        ]
      )

      cbind(
        values[
          index
        ],
        remaining_permutations
      )
    }
  )

  do.call(
    rbind,
    permutation_blocks
  )
}


build_symmetric_matrix <- function(
  pairwise_data,
  value_column,
  sample_order
) {

  output_matrix <- matrix(
    0,
    nrow = length(
      sample_order
    ),
    ncol = length(
      sample_order
    ),
    dimnames = list(
      sample_order,
      sample_order
    )
  )

  for (row_index in seq_len(
    nrow(
      pairwise_data
    )
  )) {

    sample_1 <- as.character(
      pairwise_data$Sample_1[
        row_index
      ]
    )

    sample_2 <- as.character(
      pairwise_data$Sample_2[
        row_index
      ]
    )

    if (
      !(sample_1 %in% sample_order) ||
      !(sample_2 %in% sample_order)
    ) {

      stop(
        paste0(
          "Unknown biological sample in pairwise table: ",
          sample_1,
          " or ",
          sample_2
        ),
        call. = FALSE
      )
    }

    current_value <- as.numeric(
      pairwise_data[[
        value_column
      ]][
        row_index
      ]
    )

    output_matrix[
      sample_1,
      sample_2
    ] <- current_value

    output_matrix[
      sample_2,
      sample_1
    ] <- current_value
  }

  output_matrix
}


# Reviewer 3 matrix calculations use r3_matrix_tests below.

save_ggplot <- function(
  plot_object,
  output_folder,
  filename,
  width,
  height
) {

  ggplot2::ggsave(
    filename = file.path(
      output_folder,
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
      output_folder,
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
          output_folder,
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
        output_folder,
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


theme_publication <- function(
  base_size = 10
) {

  ggplot2::theme_classic(
    base_size = base_size,
    base_family = "sans"
  ) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),

      panel.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),

      panel.border = ggplot2::element_rect(
        fill = NA,
        colour = "#262626",
        linewidth = 0.55
      ),

      panel.grid.major.y = ggplot2::element_line(
        colour = "#E6E6E6",
        linewidth = 0.30
      ),

      panel.grid.major.x = ggplot2::element_blank(),

      panel.grid.minor = ggplot2::element_blank(),

      axis.line = ggplot2::element_blank(),

      axis.ticks = ggplot2::element_line(
        colour = "#262626",
        linewidth = 0.45
      ),

      axis.ticks.length = grid::unit(
        1.6,
        "mm"
      ),

      axis.title = ggplot2::element_text(
        face = "bold",
        colour = "#1A1A1A",
        size = AXIS_TITLE_SIZE
      ),

      axis.text = ggplot2::element_text(
        colour = "#1A1A1A",
        size = AXIS_TEXT_SIZE
      ),

      strip.background = ggplot2::element_rect(
        fill = "#F2F2F2",
        colour = "#262626",
        linewidth = 0.50
      ),

      strip.text = ggplot2::element_text(
        face = "bold",
        colour = "#1A1A1A",
        size = FACET_TEXT_SIZE
      ),

      legend.position = "top",

      legend.direction = "horizontal",

      legend.box = "vertical",

      legend.justification = "left",

      legend.title = ggplot2::element_text(
        face = "bold",
        size = LEGEND_TITLE_SIZE
      ),

      legend.text = ggplot2::element_text(
        size = LEGEND_TEXT_SIZE
      ),

      legend.key = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),

      legend.spacing.x = grid::unit(
        2.2,
        "mm"
      ),

      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0,
        colour = "#1A1A1A",
        size = base_size + 1
      ),

      plot.subtitle = ggplot2::element_text(
        hjust = 0,
        colour = "#4D4D4D",
        size = base_size
      ),

      plot.caption = ggplot2::element_text(
        hjust = 0,
        colour = "#4D4D4D",
        size = base_size - 1
      ),

      plot.margin = ggplot2::margin(
        7,
        9,
        7,
        7
      )
    )
}


# Figure 8 primary TP and reference-sensitivity module.
# Primary point estimates use observed means and fixed TDF; bootstrap medians
# are retained only as diagnostics. No historical TP/QC field filters fish.
analyse_figure8_tp <- function(fish, resources, baselines,
                               draws = N_TP_DRAWS, seed = SEED + 100L) {
  keys <- c("HC|Wet_2023", "HC|Dry_2023", "QJ|Wet_2023", "QJ|Dry_2023",
            "SY|Wet_2023", "SY|Dry_2023")
  sample_ids <- setNames(gsub("|", "__", keys, fixed = TRUE), keys)
  target <- "Rusty-patina Ring-ribbed Snail"
  if (length(draws) != 1L || !is.finite(draws) || draws < 1000L)
    stop("N_TP_DRAWS must be a single integer of at least 1000.")
  if (!is.finite(TDF_MEAN) || TDF_MEAN <= 0)
    stop("The primary fixed discrimination factor must be positive.")
  if (nrow(fish) != 609L || anyDuplicated(fish$Fish_isotope_record_id) ||
      any(!is.finite(fish$d15N_permil)))
    stop("Figure 8 requires all 609 unique fish with finite measured delta15N.")
  if (anyDuplicated(resources$Resource_record_id))
    stop("Resource_record_id must be unique; baseline matching is otherwise ambiguous.")
  fish$Fish_group_key <- paste(fish$Reach_code, fish$Campaign, sep = "|")
  if (!setequal(unique(fish$Fish_group_key), keys) ||
      nrow(baselines) != 6L || anyDuplicated(baselines$Fish_group_key) ||
      !setequal(baselines$Fish_group_key, keys))
    stop("Fish and TP_Baselines must contain exactly the six expected groups.")

  fish_draws <- matrix(NA_real_, draws, 6L, dimnames = list(NULL, keys))
  baseline_draws <- fish_draws
  group_rows <- baseline_rows <- list()
  set.seed(seed)
  for (i in seq_along(keys)) {
    key <- keys[i]
    parts <- strsplit(key, "|", fixed = TRUE)[[1L]]
    b <- baselines[match(key, baselines$Fish_group_key), , drop = FALSE]
    ids <- trimws(strsplit(as.character(b$Resource_record_ids), ";", fixed = TRUE)[[1L]])
    if (any(!nzchar(ids)) || anyDuplicated(ids) ||
        length(ids) != as.integer(b$Baseline_n) ||
        anyNA(match(ids, resources$Resource_record_id)))
      stop(paste("Missing, duplicated or miscounted baseline IDs for", key))
    r <- resources[match(ids, resources$Resource_record_id), , drop = FALSE]
    if (any(r$Reach_code != parts[1L]) || any(r$Campaign != parts[2L]) ||
        any(r$Resource_type_original != target) || any(r$Resource_class != "Gastropod") ||
        any(!is.finite(r$d15N_permil)))
      stop(paste("Baseline must be the target gastropod from the same reach and campaign:", key))
    if (as.character(b$Fish_site_code) != parts[1L] ||
        as.character(b$Fish_season) != parts[2L] ||
        as.character(b$Source_site_code) != parts[1L] ||
        as.character(b$Source_season) != parts[2L] ||
        toupper(as.character(b$Baseline_use_type)) != "OBSERVED" ||
        as.character(b$Site_match) != "PASS" ||
        as.character(b$Season_rule_check) != "PASS")
      stop(paste("Baseline metadata must confirm an observed same-site same-campaign reference:", key))
    observed_b <- mean(r$d15N_permil)
    if (!is.finite(as.numeric(b$Mean_d15N_permil)) ||
        abs(as.numeric(b$Mean_d15N_permil) - observed_b) > 1e-8 ||
        abs(as.numeric(b$SD_d15N_permil) - stats::sd(r$d15N_permil)) > 1e-8)
      stop(paste("Baseline summaries disagree with their underlying measurements:", key))
    f <- fish$d15N_permil[fish$Fish_group_key == key]
    if (!length(f)) stop(paste("No measured fish for", key))
    # One reference mean per group per replicate; never independent baselines
    # assigned to individual fish in a community-mean uncertainty calculation.
    fish_draws[, i] <- replicate(draws, mean(sample(f, length(f), replace = TRUE)))
    baseline_draws[, i] <- replicate(draws,
      mean(sample(r$d15N_permil, nrow(r), replace = TRUE)))
    group_rows[[i]] <- data.frame(
      Fish_group_key = key, Sample_ID = unname(sample_ids[key]),
      Group_label = paste(unname(REACH_LABELS[parts[1L]]),
                          if (parts[2L] == "Wet_2023") "Wet" else "Dry", sep = " - "),
      Reach_code = parts[1L], Campaign = parts[2L],
      Fish_n = length(f), Baseline_n = nrow(r), Baseline_use_type = "OBSERVED",
      Fish_mean_d15N = mean(f), Baseline_mean_d15N = observed_b,
      Fixed_TDF = TDF_MEAN,
      Point_TP = TROPHIC_LEVEL_BASELINE + (mean(f) - observed_b) / TDF_MEAN,
      stringsAsFactors = FALSE)
    baseline_rows[[i]] <- data.frame(
      Fish_group_key = key, Baseline_use_type = "OBSERVED",
      Resource_record_id = r$Resource_record_id,
      Source_reach = r$Reach_code, Source_campaign = r$Campaign,
      Resource_type = r$Resource_type_original,
      Baseline_d15N_permil = r$d15N_permil, Site_match = "PASS",
      Campaign_match = "PASS", Taxon_match = "PASS", stringsAsFactors = FALSE)
  }
  groups <- do.call(rbind, group_rows)
  rownames(groups) <- NULL
  tp_draws <- TROPHIC_LEVEL_BASELINE + (fish_draws - baseline_draws) / TDF_MEAN
  groups$Lower_95_TP <- apply(tp_draws, 2L, quantile_025)
  groups$Upper_95_TP <- apply(tp_draws, 2L, quantile_975)
  groups$Bootstrap_median_TP <- apply(tp_draws, 2L, stats::median)
  groups$Bootstrap_SD_TP <- apply(tp_draws, 2L, stats::sd)

  contrast_rows <- lapply(c("HC", "QJ", "SY"), function(reach) {
    w <- paste0(reach, "|Wet_2023"); d <- paste0(reach, "|Dry_2023")
    wi <- match(w, keys); di <- match(d, keys)
    delta <- tp_draws[, di] - tp_draws[, wi]
    data.frame(Reach_code = reach, Wet_group = w, Dry_group = d,
      Fish_delta_d15N = groups$Fish_mean_d15N[di] - groups$Fish_mean_d15N[wi],
      Baseline_delta_d15N = groups$Baseline_mean_d15N[di] - groups$Baseline_mean_d15N[wi],
      Point_Dry_minus_Wet_TP = groups$Point_TP[di] - groups$Point_TP[wi],
      Lower_95 = quantile_025(delta), Upper_95 = quantile_975(delta),
      Bootstrap_median_Dry_minus_Wet_TP = stats::median(delta),
      Bootstrap_fraction_Dry_greater_Wet = mean(delta > 0),
      stringsAsFactors = FALSE)
  })
  contrasts <- do.call(rbind, contrast_rows)
  pairs <- utils::combn(seq_along(keys), 2L)
  pairwise <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(j) {
    i1 <- pairs[1L, j]; i2 <- pairs[2L, j]
    delta <- tp_draws[, i2] - tp_draws[, i1]
    point <- groups$Point_TP[i2] - groups$Point_TP[i1]
    data.frame(Sample_1 = unname(sample_ids[keys[i1]]),
      Sample_2 = unname(sample_ids[keys[i2]]),
      Point_TP_difference_2_minus_1 = point,
      Absolute_point_TP_difference = abs(point),
      Lower_95_TP_difference = quantile_025(delta),
      Upper_95_TP_difference = quantile_975(delta),
      Bootstrap_median_TP_difference = stats::median(delta), stringsAsFactors = FALSE)
  }))

  tdf_groups <- do.call(rbind, lapply(FIXED_TDF_VALUES, function(tdf) {
    values <- TROPHIC_LEVEL_BASELINE + (fish_draws - baseline_draws) / tdf
    data.frame(Fish_group_key = keys, Sample_ID = unname(sample_ids[keys]),
      Fixed_TDF = tdf,
      Point_TP = TROPHIC_LEVEL_BASELINE +
        (groups$Fish_mean_d15N - groups$Baseline_mean_d15N) / tdf,
      Lower_95 = apply(values, 2L, quantile_025),
      Upper_95 = apply(values, 2L, quantile_975), row.names = NULL)
  }))
  tdf_contrasts <- do.call(rbind, lapply(FIXED_TDF_VALUES, function(tdf) {
    data.frame(Reach_code = contrasts$Reach_code, Fixed_TDF = tdf,
      Point_Dry_minus_Wet_TP = contrasts$Point_Dry_minus_Wet_TP * TDF_MEAN / tdf,
      Lower_95 = contrasts$Lower_95 * TDF_MEAN / tdf,
      Upper_95 = contrasts$Upper_95 * TDF_MEAN / tdf)
  }))

  eligible_resources <- resources[resources$Campaign %in% MATCHED_RESOURCE_CAMPAIGNS &
    resources$Resource_class %in% c("Gastropod", "Zooplankton", "Bivalve"), , drop = FALSE]
  inv_groups <- unique(eligible_resources[c("Reach_code", "Resource_type_original", "Resource_class")])
  inventory <- do.call(rbind, lapply(seq_len(nrow(inv_groups)), function(i) {
    rr <- eligible_resources[eligible_resources$Reach_code == inv_groups$Reach_code[i] &
      eligible_resources$Resource_type_original == inv_groups$Resource_type_original[i], , drop = FALSE]
    wn <- sum(rr$Campaign == "Wet_2023"); dn <- sum(rr$Campaign == "Dry_2023")
    status <- if (!wn || !dn) "Not paired across campaigns" else if (
      inv_groups$Resource_type_original[i] == target) "Primary gastropod reference" else
      "Exploratory only: taxonomic composition and feeding status unverified"
    data.frame(Reach_code = inv_groups$Reach_code[i],
      Resource_type = inv_groups$Resource_type_original[i],
      Resource_class = inv_groups$Resource_class[i], Wet_n = wn, Dry_n = dn,
      Assessment = status, stringsAsFactors = FALSE)
  }))

  # Reference sensitivity is separate from fixed-factor sensitivity.
  # Bulk zooplankton is NOT silently promoted to a verified TP=2 baseline.
  sensitivity_rows <- list(); si <- 1L
  set.seed(seed + 1000L)
  for (reach in c("HC", "QJ", "SY")) {
    wi <- match(paste0(reach, "|Wet_2023"), keys)
    di <- match(paste0(reach, "|Dry_2023"), keys)
    fish_delta <- groups$Fish_mean_d15N[di] - groups$Fish_mean_d15N[wi]
    for (reference in c(target, "Zooplankton")) {
      rw <- resources$d15N_permil[resources$Reach_code == reach &
        resources$Campaign == "Wet_2023" & resources$Resource_type_original == reference]
      rd <- resources$d15N_permil[resources$Reach_code == reach &
        resources$Campaign == "Dry_2023" & resources$Resource_type_original == reference]
      if (!length(rw) || !length(rd)) next
      if (reference == target) {
        delta_draws <- tp_draws[, di] - tp_draws[, wi]
      } else {
        rw_draws <- replicate(draws, mean(sample(rw, length(rw), replace = TRUE)))
        rd_draws <- replicate(draws, mean(sample(rd, length(rd), replace = TRUE)))
        delta_draws <- ((fish_draws[, di] - fish_draws[, wi]) - (rd_draws - rw_draws)) / TDF_MEAN
      }
      baseline_delta <- mean(rd) - mean(rw)
      sensitivity_rows[[si]] <- data.frame(Reach_code = reach,
        Reference = if (reference == target) "S. quadrata" else "Bulk zooplankton",
        Wet_n = length(rw), Dry_n = length(rd),
        Wet_reference_mean_d15N = mean(rw), Dry_reference_mean_d15N = mean(rd),
        Fish_delta_d15N = fish_delta, Reference_delta_d15N = baseline_delta,
        Corrected_delta_d15N = fish_delta - baseline_delta,
        Scaled_Dry_minus_Wet = (fish_delta - baseline_delta) / TDF_MEAN,
        Lower_95 = quantile_025(delta_draws), Upper_95 = quantile_975(delta_draws),
        Interpretation = if (reference == target)
          "Gastropod-referenced TP; conditional on pathway representation and tissue integration" else
          "Exploratory scaled isotope contrast; bulk zooplankton taxonomy and feeding status unverified",
        stringsAsFactors = FALSE)
      si <- si + 1L
    }
  }
  ref_sensitivity <- do.call(rbind, sensitivity_rows)
  s8d_I <- data.frame(Reach_code = contrasts$Reach_code,
    Fish_delta_d15N = contrasts$Fish_delta_d15N,
    Baseline_delta_d15N = contrasts$Baseline_delta_d15N,
    Corrected_delta_d15N = contrasts$Fish_delta_d15N - contrasts$Baseline_delta_d15N,
    Delta_TP_2.9 = contrasts$Point_Dry_minus_Wet_TP * TDF_MEAN / 2.9,
    Delta_TP_3.4 = contrasts$Point_Dry_minus_Wet_TP,
    Delta_TP_3.9 = contrasts$Point_Dry_minus_Wet_TP * TDF_MEAN / 3.9,
    Lower_95_TP_3.4 = contrasts$Lower_95, Upper_95_TP_3.4 = contrasts$Upper_95)
  s8d_II <- data.frame(Reach_code = groups$Reach_code, Campaign = groups$Campaign,
    Fish_n = groups$Fish_n, Fish_mean_d15N = groups$Fish_mean_d15N,
    Baseline_campaign = groups$Campaign, Baseline_n = groups$Baseline_n,
    Baseline_mean_d15N = groups$Baseline_mean_d15N,
    TP_2.9 = TROPHIC_LEVEL_BASELINE + (groups$Fish_mean_d15N - groups$Baseline_mean_d15N) / 2.9,
    TP_3.4 = groups$Point_TP,
    TP_3.9 = TROPHIC_LEVEL_BASELINE + (groups$Fish_mean_d15N - groups$Baseline_mean_d15N) / 3.9,
    Lower_95_TP_3.4 = groups$Lower_95_TP, Upper_95_TP_3.4 = groups$Upper_95_TP)
  list(groups = groups, contrasts = contrasts, pairwise = pairwise,
    baseline_raw = do.call(rbind, baseline_rows), inventory = inventory,
    tdf_groups = tdf_groups, tdf_contrasts = tdf_contrasts,
    reference_sensitivity = ref_sensitivity, s8d_I = s8d_I, s8d_II = s8d_II,
    fish_draws = fish_draws, baseline_draws = baseline_draws, tp_draws = tp_draws)
}

save_combined_figure8 <- function(panel_a, panel_b, folder) {
  draw_panels <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(
      2L, 1L, heights = grid::unit(c(1, 1.65), "null"))))
    print(panel_a + ggplot2::labs(tag = "(a)"), newpage = FALSE,
      vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 1L))
    print(panel_b + ggplot2::labs(tag = "(b)"), newpage = FALSE,
      vp = grid::viewport(layout.pos.row = 2L, layout.pos.col = 1L))
    grid::popViewport()
  }
  save_device <- function(kind) {
    filename <- file.path(folder, paste0("Figure_8_gastropod_TP_environmental_concordance.", kind))
    if (kind == "png") {
      grDevices::png(filename, width = 9, height = 12, units = "in", res = 600, bg = "white")
    } else if (capabilities("cairo")) {
      grDevices::cairo_pdf(filename, width = 9, height = 12, bg = "white")
    } else {
      grDevices::pdf(filename, width = 9, height = 12, useDingbats = FALSE)
    }
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_panels()
  }
  save_device("png")
  save_device("pdf")
  invisible(NULL)
}


# =============================================================================
# 4. MAIN ANALYSIS
# =============================================================================
# =============================================================================
# 4. MAIN ANALYSIS
# =============================================================================

# =============================================================================
# REVIEWER 3: EXPLORATORY COMPARISON AND RESTRICTED NULL SENSITIVITIES
# The main scientific result is descriptive. Neither null is guaranteed by an
# observational river survey. Exact enumeration is conditional on the stated
# exchangeability assumption and does not increase biological replication.
# =============================================================================
R3_DISPLAY_SCHEME <- "Whole_reach_pairs"  # Reference output, not a selected significant null.
R3_SHOW_P_IN_FIGURES <- FALSE

r3_all_perms <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], r3_all_perms(x[-i]))
  }))
}

r3_make_designs <- function(group_order) {
  parts <- strsplit(as.character(group_order), "__", fixed = TRUE)
  if (length(parts) != 6L || any(lengths(parts) != 2L) || anyDuplicated(group_order))
    stop("Reviewer-3 design requires six unique Reach__Campaign identifiers.")
  reach <- vapply(parts, `[`, character(1), 1L)
  campaign <- vapply(parts, `[`, character(1), 2L)
  reaches <- unique(reach)
  campaigns <- unique(campaign)
  if (length(reaches) != 3L || length(campaigns) != 2L ||
      any(table(reach, campaign) != 1L))
    stop("Reviewer-3 design requires exactly three reaches, each with two campaigns.")
  rp <- r3_all_perms(seq_along(reaches))
  whole <- t(vapply(seq_len(nrow(rp)), function(k) {
    source_reach <- reaches[rp[k, match(reach, reaches)]]
    match(paste(source_reach, campaign, sep = "__"), group_order)
  }, integer(length(group_order))))
  flips <- as.matrix(expand.grid(rep(list(0:1), length(reaches))))
  within <- t(vapply(seq_len(nrow(flips)), function(k) {
    p <- seq_along(group_order)
    for (i in seq_along(reaches)) {
      ix <- which(reach == reaches[i])
      if (flips[k, i] == 1L) p[ix] <- rev(ix)
    }
    p
  }, integer(length(group_order))))
  designs <- list(Whole_reach_pairs = whole, Within_reach_swaps = within)
  for (nm in names(designs)) {
    p <- designs[[nm]]
    stopifnot(!anyNA(p), ncol(p) == 6L,
              nrow(unique(p)) == nrow(p),
              any(apply(p, 1L, function(z) all(z == seq_along(group_order)))))
    if (nm == "Whole_reach_pairs") {
      stopifnot(nrow(p) == 6L)
      for (k in seq_len(nrow(p))) {
        stopifnot(all(campaign[p[k, ]] == campaign))
        for (r in reaches) stopifnot(length(unique(reach[p[k, reach == r]])) == 1L)
      }
    } else {
      stopifnot(nrow(p) == 8L)
      for (k in seq_len(nrow(p))) stopifnot(all(reach[p[k, ]] == reach))
    }
  }
  designs
}

r3_nonidentity <- function(p) {
  p[!apply(p, 1L, function(z) all(z == seq_len(ncol(p)))), , drop = FALSE]
}

r3_null_hypotheses <- function() {
  data.frame(
    Permutation_scheme = c("Whole_reach_pairs", "Within_reach_swaps"),
    Null_hypothesis = c(
      "Correspondence of complete wet-dry reach trajectories between matrices is exchangeable, retaining campaign order.",
      "Within each reach, correspondence of the two campaign observations between matrices is exchangeable, retaining reach identity."),
    Required_assumption = c(
      "The three whole reach trajectories are exchangeable under the null; spatial structure is not guaranteed to satisfy this.",
      "Campaign observations are exchangeable within reach under the null; unmodelled shared seasonality can violate this."),
    Upper_tail_alternative = "Observed matrix correlation is unusually large and positive under the stated null.",
    Arrangements_including_observed = c(6L, 8L),
    Nonidentity_permutations = c(5L, 7L),
    Minimum_attainable_P = c(1/6, 1/8),
    Interpretation = "Conditional exploratory sensitivity; not a test of complementarity or independence.",
    stringsAsFactors = FALSE
  )
}

r3_permutation_table <- function(designs) {
  do.call(rbind, lapply(names(designs), function(nm) {
    p <- designs[[nm]]
    result <- data.frame(Permutation_scheme = nm, Arrangement = seq_len(nrow(p)),
      Observed_order = apply(p, 1L, function(z) all(z == seq_len(ncol(p)))))
    p <- as.data.frame(p); names(p) <- paste0("Index_", seq_len(ncol(p)))
    cbind(result, p)
  }))
}

r3_mantel_core <- function(a, b, p, method = "spearman") {
  stopifnot(identical(dim(a), dim(b)), nrow(a) == 6L, ncol(p) == 6L,
            all(is.finite(a)), all(is.finite(b)),
            max(abs(a-t(a))) < 1e-10, max(abs(b-t(b))) < 1e-10)
  if (!is.null(rownames(a)) && !is.null(rownames(b)) &&
      !identical(rownames(a), rownames(b))) stop("Matrix group order mismatch.")
  ix <- lower.tri(a)
  observed <- as.numeric(stats::cor(a[ix], b[ix], method = method))
  permuted <- vapply(seq_len(nrow(p)), function(i) {
    x <- a[p[i, ], p[i, ], drop = FALSE]
    as.numeric(stats::cor(x[ix], b[ix], method = method))
  }, numeric(1))
  if (!is.finite(observed) || any(!is.finite(permuted)))
    stop("Undefined matrix correlation; check for constant or missing distances.")
  # p includes the observed identity exactly once. Do not add a second +1.
  list(statistic = observed,
       upper_P = mean(permuted >= observed - 1e-12),
       two_P = mean(abs(permuted) >= abs(observed) - 1e-12))
}

r3_matrix_tests <- function(left, right, methods, designs,
                           labels = c("Community_metric", "Trophic_reorganization_metric")) {
  rows <- list(); k <- 1L
  for (scheme in names(designs)) for (method in methods) {
    for (ln in names(left)) for (rn in names(right)) {
      p <- designs[[scheme]]
      z <- r3_mantel_core(left[[ln]], right[[rn]], p, method)
      d <- data.frame(Left_metric = ln, Right_metric = rn,
        Correlation_method = method, Matrix_correlation = z$statistic,
        Permutation_scheme = scheme, Exact_one_sided_P = z$upper_P,
        Exact_two_sided_P = z$two_P, Permutations = nrow(p)-1L,
        Arrangements_including_observed = nrow(p), Minimum_attainable_P = 1/nrow(p),
        Direction = if (z$statistic > 0) "Positive" else if (z$statistic < 0) "Negative" else "Zero",
        Inference_status = "Conditional exploratory sensitivity; no complementarity test",
        stringsAsFactors = FALSE)
      names(d)[1:2] <- labels
      rows[[k]] <- d; k <- k+1L
    }
  }
  out <- do.call(rbind, rows)
  out$FDR_adjusted_one_sided_P <- NA_real_
  out$FDR_adjusted_two_sided_P <- NA_real_
  for (i in split(seq_len(nrow(out)), interaction(out$Permutation_scheme,
                                                out$Correlation_method, drop = TRUE))) {
    out$FDR_adjusted_one_sided_P[i] <- stats::p.adjust(out$Exact_one_sided_P[i], "BH")
    out$FDR_adjusted_two_sided_P[i] <- stats::p.adjust(out$Exact_two_sided_P[i], "BH")
  }
  rownames(out) <- NULL
  out
}

r3_seasonal_rankings <- function(seasonal_data) {
  fields <- intersect(c("Bray_Curtis_turnover", "Jaccard_turnover",
    "Median_centroid_distance", "Absolute_SEA_change", "Absolute_log_SEA_ratio",
    "Niche_nonoverlap_proportion"), names(seasonal_data))
  do.call(rbind, lapply(fields, function(v) {
    has_interval <- v == "Median_centroid_distance" &&
      all(c("Lower_95_centroid_distance", "Upper_95_centroid_distance") %in% names(seasonal_data))
    data.frame(Metric = v, Reach_code = seasonal_data$Reach_code,
      Point_estimate = seasonal_data[[v]],
      Point_rank_1_largest = rank(-seasonal_data[[v]], ties.method = "min"),
      Lower_95_if_supplied = if (has_interval) seasonal_data$Lower_95_centroid_distance else NA_real_,
      Upper_95_if_supplied = if (has_interval) seasonal_data$Upper_95_centroid_distance else NA_real_,
      Rank_inference = "Descriptive point ranks only; no rank-difference test or posterior rank probability.",
      stringsAsFactors = FALSE)
  }))
}

r3_write_design_files <- function(folder, designs) {
  write.csv(r3_null_hypotheses(), file.path(folder, "Reviewer3_null_hypotheses.csv"), row.names = FALSE)
  write.csv(r3_permutation_table(designs), file.path(folder, "Reviewer3_permutation_indices.csv"), row.names = FALSE)
  writeLines(c(
    "Exploratory comparison of three reaches sampled in two campaigns.",
    "Main plots report observed effect estimates, not evidence of indicator complementarity.",
    "Whole-reach-pair sensitivity: 6 arrangements, including observed; minimum exact P = 1/6.",
    "Within-reach campaign-swap sensitivity: 8 arrangements, including observed; minimum exact P = 1/8.",
    "These are distinct conditional nulls. Do not select a scheme or sidedness by its P value.",
    "Upper-tail alternative: unusually large positive association; absolute-correlation two-sided sensitivity is also reported.",
    "BH adjustment is separate within each metric family, correlation method and permutation scheme.",
    "Neither non-significance nor different point ranks establishes complementarity or independence.",
    "Full posterior draws and independent reach/campaign replication are unavailable for rank validation."),
    file.path(folder, "Reviewer3_interpretation.txt"))
}


run_environmental_integration <- function() {

  # ---------------------------------------------------------------------------
  # 4.1 Locate input workbook
  # ---------------------------------------------------------------------------

  workbook_pattern <- "^Han_River_Figure8_Analysis_Ready.*[.]xlsx$"

  # RStudio Source and Rscript both locate the paired workbook beside the script.
  script_directory <- getwd()
  source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(source_file) && nzchar(source_file)) {
    script_directory <- dirname(normalizePath(source_file, mustWork = FALSE))
  } else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(active_path)) script_directory <- dirname(active_path)
  } else {
    file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_argument)) script_directory <- dirname(sub("^--file=", "", file_argument[1L]))
  }

  user_home <- path.expand(
    "~"
  )

  windows_user <- Sys.getenv(
    "USERPROFILE"
  )

  search_folders <- unique(
    c(
      script_directory,
      getwd(),
      file.path(
        windows_user,
        "Downloads"
      ),
      file.path(
        user_home,
        "Downloads"
      ),
      file.path(
        windows_user,
        "Desktop"
      ),
      file.path(
        user_home,
        "Desktop"
      ),
      file.path(
        windows_user,
        "Documents"
      ),
      file.path(
        user_home,
        "Documents"
      )
    )
  )

  search_folders <- search_folders[
    nzchar(
      search_folders
    ) &
      dir.exists(
        search_folders
      )
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
      file.exists(
        candidate_files
      )
    ]
  )

  if (length(
    candidate_files
  ) > 1L) {

    file_information <- file.info(
      candidate_files
    )

    candidate_files <- candidate_files[
      order(
        file_information$mtime,
        decreasing = TRUE,
        na.last = TRUE
      )
    ]
  }

  if (length(
    candidate_files
  ) == 0L) {

    message(
      "\nSelect the environmental-integration Excel workbook."
    )

    selected_file <- tryCatch(
      file.choose(),

      error = function(e) {

        ""
      }
    )

    if (
      !nzchar(
        selected_file
      ) ||
        !file.exists(
          selected_file
        )
    ) {

      stop(
        "No valid Excel workbook was selected.",
        call. = FALSE
      )
    }

    candidate_files <- selected_file
  }

  # Prefer the revised workbook if earlier versions are in the same folder.
  preferred_r3 <- grepl("_Reviewer3[.]xlsx$", basename(candidate_files), ignore.case = TRUE)
  candidate_files <- c(candidate_files[preferred_r3], candidate_files[!preferred_r3])

  INPUT_FILE <- normalizePath(
    candidate_files[
      1
    ],
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
  # 4.2 Output directories and status files
  # ---------------------------------------------------------------------------

  RESULTS_DIR <- file.path(
    dirname(
      INPUT_FILE
    ),
    "Han_River_Figure8_results_Reviewer3_ALL_609"
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
      "Environmental integration and revised trophic-position analysis is running.",
      paste(
        "Started:",
        Sys.time()
      ),
      paste(
        "Input workbook:",
        INPUT_FILE
      )
    ),
    RUNNING_FILE
  )


  # ---------------------------------------------------------------------------
  # 4.3 Read workbook sheets
  # ---------------------------------------------------------------------------

  required_sheets <- c(
    "Environmental_Transects",
    "Environmental_Group_Means",
    "Biological_Group_Metadata",
    "Biological_Pairwise",
    "Resource_Isotopes",
    "TP_Baselines",
    "Fish_Isotopes"
  )

  available_sheets <- readxl::excel_sheets(
    INPUT_FILE
  )

  missing_sheets <- setdiff(
    required_sheets,
    available_sheets
  )

  if (length(
    missing_sheets
  ) > 0L) {

    stop(
      paste0(
        "Workbook is missing worksheet(s): ",
        paste(
          missing_sheets,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }

  environmental_data <- readxl::read_excel(
    INPUT_FILE,
    sheet = "Environmental_Transects"
  ) |>
    as.data.frame(
      stringsAsFactors = FALSE
    )

  environmental_group_means <- readxl::read_excel(
    INPUT_FILE,
    sheet = "Environmental_Group_Means"
  ) |>
    as.data.frame(
      stringsAsFactors = FALSE
    )

  biological_groups <- readxl::read_excel(
    INPUT_FILE,
    sheet = "Biological_Group_Metadata"
  ) |>
    as.data.frame(
      stringsAsFactors = FALSE
    )

  biological_pairwise <- readxl::read_excel(
    INPUT_FILE,
    sheet = "Biological_Pairwise"
  ) |>
    as.data.frame(
      stringsAsFactors = FALSE
    )

  resource_isotopes <- readxl::read_excel(
    INPUT_FILE,
    sheet = "Resource_Isotopes"
  ) |>
    as.data.frame(
      stringsAsFactors = FALSE
    )

  tp_baselines <- readxl::read_excel(
    INPUT_FILE,
    sheet = "TP_Baselines"
  ) |>
    as.data.frame(
      stringsAsFactors = FALSE
    )

  fish_isotopes <- readxl::read_excel(
    INPUT_FILE,
    sheet = "Fish_Isotopes"
  ) |>
    as.data.frame(
      stringsAsFactors = FALSE
    )


  # ---------------------------------------------------------------------------
  # 4.4 Input validation
  # ---------------------------------------------------------------------------

  required_environmental_columns <- c(
    "Environmental_sample_ID",
    "Biological_group_ID",
    "Campaign",
    "Reach_code",
    "Reach_name",
    "Transect_position",
    PRIMARY_ENVIRONMENTAL_VARIABLES
  )

  required_pairwise_columns <- c(
    "Pair_ID",
    "Sample_1",
    "Sample_2",
    "Bray_Curtis",
    "Jaccard",
    "Median_centroid_distance",
    "Absolute_SEA_difference",
    "Absolute_log_SEA_ratio",
    "Niche_nonoverlap_proportion"
  )

  required_group_columns <- c(
    "Sample_ID",
    "Reach_code",
    "Campaign",
    "Posterior_median_d13C",
    "Posterior_median_d15N"
  )

  required_resource_columns <- c(
    "Resource_record_id",
    "Reach_code",
    "Reach",
    "Campaign",
    "Trophic_role",
    "Resource_class",
    "Baseline_candidate",
    "d15N_permil",
    "d13C_permil"
  )

  required_baseline_columns <- c(
    "Fish_group_key",
    "Baseline_use_type",
    "Baseline_n",
    "Resource_record_ids"
  )

  required_fish_columns <- c(
    "Fish_isotope_record_id",
    "Reach_code",
    "Reach",
    "Campaign",
    "Feeding_guild",
    "Species_code",
    "d15N_permil",
    "d13C_permil",
    "TP_primary_include"
  )

  validation_list <- list(
    Environmental_Transects = setdiff(
      required_environmental_columns,
      names(
        environmental_data
      )
    ),

    Biological_Pairwise = setdiff(
      required_pairwise_columns,
      names(
        biological_pairwise
      )
    ),

    Biological_Group_Metadata = setdiff(
      required_group_columns,
      names(
        biological_groups
      )
    ),

    Resource_Isotopes = setdiff(
      required_resource_columns,
      names(
        resource_isotopes
      )
    ),

    TP_Baselines = setdiff(
      required_baseline_columns,
      names(
        tp_baselines
      )
    ),

    Fish_Isotopes = setdiff(
      required_fish_columns,
      names(
        fish_isotopes
      )
    )
  )

  validation_failures <- names(
    validation_list
  )[
    vapply(
      validation_list,
      length,
      FUN.VALUE = integer(1)
    ) > 0L
  ]

  if (length(
    validation_failures
  ) > 0L) {

    failure_messages <- vapply(
      validation_failures,
      function(sheet_name) {

        paste0(
          sheet_name,
          ": ",
          paste(
            validation_list[[
              sheet_name
            ]],
            collapse = ", "
          )
        )
      },
      FUN.VALUE = character(1)
    )

    stop(
      paste0(
        "Missing required columns:\n",
        paste(
          failure_messages,
          collapse = "\n"
        )
      ),
      call. = FALSE
    )
  }

  # Standardize all user-facing reach fields immediately after validation.
  # Internal HC/QJ/SY codes are retained for joins, models, and matrix order.
  environmental_data$Reach_name <- unname(
    REACH_LABELS[as.character(environmental_data$Reach_code)]
  )

  environmental_group_means$Reach_name <- unname(
    REACH_LABELS[as.character(environmental_group_means$Reach_code)]
  )

  biological_groups$Reach_name <- unname(
    REACH_LABELS[as.character(biological_groups$Reach_code)]
  )

  resource_isotopes$Reach <- unname(
    REACH_LABELS[as.character(resource_isotopes$Reach_code)]
  )

  fish_isotopes$Reach <- unname(
    REACH_LABELS[as.character(fish_isotopes$Reach_code)]
  )

  if (
    anyNA(environmental_data$Reach_name) ||
      anyNA(environmental_group_means$Reach_name) ||
      anyNA(biological_groups$Reach_name) ||
      anyNA(resource_isotopes$Reach) ||
      anyNA(fish_isotopes$Reach)
  ) {
    stop(
      "A reach code could not be mapped to S1, S2, or S3.",
      call. = FALSE
    )
  }

  if (nrow(
    environmental_data
  ) != 18L) {

    stop(
      paste0(
        "Environmental_Transects must contain 18 rows; observed ",
        nrow(
          environmental_data
        ),
        "."
      ),
      call. = FALSE
    )
  }

  for (variable_name in PRIMARY_ENVIRONMENTAL_VARIABLES) {

    environmental_data[[
      variable_name
    ]] <- safe_numeric(
      environmental_data[[
        variable_name
      ]],
      variable_name
    )

    environmental_group_means[[
      variable_name
    ]] <- safe_numeric(
      environmental_group_means[[
        variable_name
      ]],
      variable_name
    )
  }

  if (
    any(!is.finite(environmental_data$DO_mg_L)) ||
      any(environmental_data$DO_mg_L < 0)
  ) {

    stop(
      "DO_mg_L must contain 18 finite, non-negative master-dataset values.",
      call. = FALSE
    )
  }

  dissolved_oxygen_audit <- data.frame(
    Check = c(
      "Complete DO observations",
      "Minimum DO (mg/L)",
      "Maximum DO (mg/L)",
      "Multivariate screening priority"
    ),
    Observed = c(
      sum(is.finite(environmental_data$DO_mg_L)),
      min(environmental_data$DO_mg_L),
      max(environmental_data$DO_mg_L),
      "Enabled"
    ),
    Expected = c(
      18,
      "Master dataset",
      "Master dataset",
      "DO_mg_L retained"
    ),
    Status = "PASS",
    stringsAsFactors = FALSE
  )

  numeric_pairwise_columns <- c(
    "Bray_Curtis",
    "Jaccard",
    "Median_centroid_distance",
    "Absolute_SEA_difference",
    "Absolute_log_SEA_ratio",
    "Niche_nonoverlap_proportion"
  )

  for (variable_name in numeric_pairwise_columns) {

    biological_pairwise[[
      variable_name
    ]] <- safe_numeric(
      biological_pairwise[[
        variable_name
      ]],
      variable_name
    )
  }

  biological_groups$Posterior_median_d13C <- safe_numeric(
    biological_groups$Posterior_median_d13C,
    "Posterior_median_d13C"
  )

  biological_groups$Posterior_median_d15N <- safe_numeric(
    biological_groups$Posterior_median_d15N,
    "Posterior_median_d15N"
  )

  resource_isotopes$d15N_permil <- safe_numeric(
    resource_isotopes$d15N_permil,
    "resource d15N_permil"
  )

  resource_isotopes$d13C_permil <- safe_numeric(
    resource_isotopes$d13C_permil,
    "resource d13C_permil"
  )

  fish_isotopes$d15N_permil <- safe_numeric(
    fish_isotopes$d15N_permil,
    "fish d15N_permil"
  )

  fish_isotopes$d13C_permil <- safe_numeric(
    fish_isotopes$d13C_permil,
    "fish d13C_permil"
  )

  # -------------------------------------------------------------------------
  # Corrected primary trophic-position inclusion rule
  # -------------------------------------------------------------------------
  # All 609 fish are retained. Historical trophic-level values and original
  # QC flags are provenance fields only and must never determine inclusion.
  # No observation is excluded because a calculated TP is below 1 or because
  # the value appears biologically unexpected.

  fish_isotopes$TP_primary_include <- as.logical(
    fish_isotopes$TP_primary_include
  )

  if (nrow(
    fish_isotopes
  ) != 609L) {

    stop(
      paste0(
        "The corrected primary TP analysis requires all 609 fish; observed ",
        nrow(
          fish_isotopes
        ),
        "."
      ),
      call. = FALSE
    )
  }

  if (anyDuplicated(
    fish_isotopes$Fish_isotope_record_id
  )) {

    stop(
      "Fish_isotope_record_id contains duplicates.",
      call. = FALSE
    )
  }

  if (
    anyNA(
      fish_isotopes$TP_primary_include
    ) ||
      !all(
        fish_isotopes$TP_primary_include
      )
  ) {

    stop(
      paste0(
        "Every fish must have TP_primary_include = TRUE in the corrected ",
        "primary analysis. Outcome-dependent TP exclusions are prohibited."
      ),
      call. = FALSE
    )
  }

  group_order <- c(
    "HC__Wet_2023",
    "HC__Dry_2023",
    "QJ__Wet_2023",
    "QJ__Dry_2023",
    "SY__Wet_2023",
    "SY__Dry_2023"
  )

  if (
    !setequal(
      environmental_group_means$Biological_group_ID,
      group_order
    )
  ) {

    stop(
      "Environmental_Group_Means does not contain the six expected group identifiers.",
      call. = FALSE
    )
  }

  if (
    !setequal(
      biological_groups$Sample_ID,
      group_order
    )
  ) {

    stop(
      "Biological_Group_Metadata does not contain the six expected group identifiers.",
      call. = FALSE
    )
  }

  if (nrow(
    biological_pairwise
  ) != 15L) {

    stop(
      "Biological_Pairwise must contain all 15 pairwise comparisons.",
      call. = FALSE
    )
  }

  # Every reach–campaign baseline is an observed baseline. In particular,
  # S1 Wet_2023 follows the same observed-baseline workflow as the
  # other five groups.

  tp_baselines$Baseline_use_type <- toupper(
    trimws(
      as.character(
        tp_baselines$Baseline_use_type
      )
    )
  )

  expected_baseline_groups <- c(
    "HC|Wet_2023",
    "HC|Dry_2023",
    "QJ|Wet_2023",
    "QJ|Dry_2023",
    "SY|Wet_2023",
    "SY|Dry_2023"
  )

  if (
    nrow(
      tp_baselines
    ) != length(
      expected_baseline_groups
    ) ||
      !setequal(
        tp_baselines$Fish_group_key,
        expected_baseline_groups
      )
  ) {

    stop(
      "TP_Baselines must contain exactly the six expected reach–campaign groups.",
      call. = FALSE
    )
  }

  if (
    anyNA(
      tp_baselines$Baseline_use_type
    ) ||
      !all(
        tp_baselines$Baseline_use_type == "OBSERVED"
      )
  ) {

    stop(
      paste0(
        "All six TP_Baselines rows, including HC|Wet_2023, must have ",
        "Baseline_use_type = OBSERVED."
      ),
      call. = FALSE
    )
  }

  environmental_data <- environmental_data |>
    dplyr::mutate(
      Reach_code = factor(
        Reach_code,
        levels = c(
          "HC",
          "QJ",
          "SY"
        )
      ),

      Campaign = factor(
        Campaign,
        levels = c(
          "Wet_2023",
          "Dry_2023"
        )
      ),

      Transect_position = factor(
        Transect_position,
        levels = c(
          "Left",
          "Mid",
          "Right"
        )
      )
    ) |>
    dplyr::arrange(
      Reach_code,
      Campaign,
      Transect_position
    )

  environmental_group_means <- environmental_group_means |>
    dplyr::mutate(
      Biological_group_ID = factor(
        Biological_group_ID,
        levels = group_order
      )
    ) |>
    dplyr::arrange(
      Biological_group_ID
    )

  biological_groups <- biological_groups |>
    dplyr::mutate(
      Sample_ID = factor(
        Sample_ID,
        levels = group_order
      )
    ) |>
    dplyr::arrange(
      Sample_ID
    )


  # ---------------------------------------------------------------------------
  # 4.5 Transform environmental variables
  # ---------------------------------------------------------------------------

  transformed_environmental_data <- environmental_data

  transformed_group_means <- environmental_group_means

  transformation_summary <- data.frame(
    Variable = PRIMARY_ENVIRONMENTAL_VARIABLES,
    Transformation = ifelse(
      PRIMARY_ENVIRONMENTAL_VARIABLES %in% LOG1P_VARIABLES,
      "log1p",
      "None"
    ),
    stringsAsFactors = FALSE
  )

  for (variable_name in LOG1P_VARIABLES) {

    if (any(
      transformed_environmental_data[[
        variable_name
      ]] < 0
    )) {

      stop(
        paste0(
          "Cannot apply log1p to negative values in ",
          variable_name,
          "."
        ),
        call. = FALSE
      )
    }

    transformed_environmental_data[[
      variable_name
    ]] <- log1p(
      transformed_environmental_data[[
        variable_name
      ]]
    )

    transformed_group_means[[
      variable_name
    ]] <- log1p(
      transformed_group_means[[
        variable_name
      ]]
    )
  }

  variable_screening <- screen_correlated_variables(
    transformed_environmental_data,
    PRIMARY_ENVIRONMENTAL_VARIABLES,
    cutoff = CORRELATION_CUTOFF,
    priority_variables = CORRELATION_PRIORITY_VARIABLES
  )

  retained_environmental_variables <- variable_screening$retained

  if (!all(
    CORRELATION_PRIORITY_VARIABLES %in%
      retained_environmental_variables
  )) {

    stop(
      paste0(
        "Priority environmental variable(s) were not retained: ",
        paste(
          setdiff(
            CORRELATION_PRIORITY_VARIABLES,
            retained_environmental_variables
          ),
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }

  if (length(
    retained_environmental_variables
  ) < 3L) {

    stop(
      "Fewer than three environmental variables remained after screening.",
      call. = FALSE
    )
  }

  variable_screening_table <- dplyr::left_join(
    transformation_summary,
    variable_screening$screening,
    by = "Variable"
  )

  write.csv(
    variable_screening_table,
    file.path(
      TABLE_DIR,
      "Table_1_environmental_variable_screening.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    dissolved_oxygen_audit,
    file.path(
      TABLE_DIR,
      "Table_1A_dissolved_oxygen_audit.csv"
    ),
    row.names = FALSE
  )

  cat(
    "\nEnvironmental variable screening completed.",
    "\nCandidate variables: ",
    length(
      PRIMARY_ENVIRONMENTAL_VARIABLES
    ),
    "\nRetained variables: ",
    length(
      retained_environmental_variables
    ),
    "\nRetained: ",
    paste(
      retained_environmental_variables,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )


  # ---------------------------------------------------------------------------
  # 4.6 Environmental PCA
  # ---------------------------------------------------------------------------

  environmental_numeric_matrix <- as.matrix(
    transformed_environmental_data[
      ,
      retained_environmental_variables,
      drop = FALSE
    ]
  )

  environmental_scaled_matrix <- scale(
    environmental_numeric_matrix,
    center = TRUE,
    scale = TRUE
  )

  environmental_pca <- stats::prcomp(
    environmental_scaled_matrix,
    center = FALSE,
    scale. = FALSE
  )

  pca_variance <- (
    environmental_pca$sdev^2
  ) /
    sum(
      environmental_pca$sdev^2
    )

  pca_variance_table <- data.frame(
    Principal_component = paste0(
      "PC",
      seq_along(
        pca_variance
      )
    ),
    Eigenvalue = environmental_pca$sdev^2,
    Variance_explained = pca_variance,
    Cumulative_variance = cumsum(
      pca_variance
    ),
    stringsAsFactors = FALSE
  )

  pca_site_scores <- data.frame(
    Environmental_sample_ID = transformed_environmental_data$Environmental_sample_ID,
    Biological_group_ID = transformed_environmental_data$Biological_group_ID,
    Reach_code = transformed_environmental_data$Reach_code,
    Reach_name = transformed_environmental_data$Reach_name,
    Campaign = transformed_environmental_data$Campaign,
    Campaign_label = transformed_environmental_data$Campaign_label,
    Transect_position = transformed_environmental_data$Transect_position,
    environmental_pca$x,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  pca_loadings <- data.frame(
    Variable = rownames(
      environmental_pca$rotation
    ),
    environmental_pca$rotation,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  write.csv(
    pca_variance_table,
    file.path(
      TABLE_DIR,
      "Table_2_environmental_PCA_variance.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    pca_site_scores,
    file.path(
      TABLE_DIR,
      "Table_S1_environmental_PCA_site_scores.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    pca_loadings,
    file.path(
      TABLE_DIR,
      "Table_S2_environmental_PCA_loadings.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 4.7 Environmental PCA figure
  # ---------------------------------------------------------------------------

  # Okabe-Ito-inspired, colourblind-friendly colours.
  reach_colours <- c(
    "HC" = "#0072B2",
    "QJ" = "#D55E00",
    "SY" = "#009E73"
  )

  reach_shapes <- c(
    "HC" = 21,
    "QJ" = 22,
    "SY" = 24
  )

  campaign_shapes <- c(
    "Wet_2023" = 21,
    "Dry_2023" = 24
  )

  environmental_variable_labels <- c(
    Depth_m = "Depth",
    Transparency_m = "Transparency",
    Flow_velocity_m_s = "Flow",
    Water_temperature_C = "Temperature",
    DO_mg_L = "DO",
    Conductivity_uS_cm = "Conductivity",
    pH = "pH",
    Water_TN_mg_L = "TN",
    Water_TP_mg_L = "TP",
    NH4_N_mg_L = expression(NH[4]^"+" * "-N"),
    NO3_N_mg_L = expression(NO[3]^"-" * "-N"),
    Phosphate_P_mg_L = expression(PO[4]^"3-" * "-P"),
    Chlorophyll_a_ug_L = "Chl-a",
    TSS_mg_L = "TSS",
    COD_mg_L = "COD"
  )

  environmental_variable_text <- c(
    Depth_m = "Depth",
    Transparency_m = "Transparency",
    Flow_velocity_m_s = "Flow",
    Water_temperature_C = "Temperature",
    DO_mg_L = "DO",
    Conductivity_uS_cm = "Conductivity",
    pH = "pH",
    Water_TN_mg_L = "TN",
    Water_TP_mg_L = "TP",
    NH4_N_mg_L = "NH4-N",
    NO3_N_mg_L = "NO3-N",
    Phosphate_P_mg_L = "PO4-P",
    Chlorophyll_a_ug_L = "Chl-a",
    TSS_mg_L = "TSS",
    COD_mg_L = "COD"
  )

  pca_plot_data <- pca_site_scores |>
    dplyr::mutate(
      Reach_code = factor(
        Reach_code,
        levels = c(
          "HC",
          "QJ",
          "SY"
        )
      ),

      Campaign = factor(
        Campaign,
        levels = c(
          "Wet_2023",
          "Dry_2023"
        )
      )
    )

  pca_group_centroids <- pca_plot_data |>
    dplyr::group_by(
      Reach_code,
      Campaign
    ) |>
    dplyr::summarise(
      Centroid_PC1 = mean(
        PC1
      ),
      Centroid_PC2 = mean(
        PC2
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Group_label = paste0(
        unname(
          REACH_LABELS[as.character(Reach_code)]
        ),
        "–",
        ifelse(
          Campaign == "Wet_2023",
          "Wet",
          "Dry"
        )
      )
    )

  pca_segment_data <- pca_plot_data |>
    dplyr::left_join(
      pca_group_centroids,
      by = c(
        "Reach_code",
        "Campaign"
      )
    )

  loading_scale <- 0.75 *
    min(
      diff(
        range(
          pca_plot_data$PC1
        )
      ) /
        max(
          abs(
            pca_loadings$PC1
          )
        ),

      diff(
        range(
          pca_plot_data$PC2
        )
      ) /
        max(
          abs(
            pca_loadings$PC2
          )
        )
    )

  loading_plot_data <- pca_loadings |>
    dplyr::mutate(
      Contribution_PC12 = PC1^2 + PC2^2
    ) |>
    dplyr::slice_max(
      order_by = Contribution_PC12,
      n = min(
        8L,
        nrow(
          pca_loadings
        )
      ),
      with_ties = FALSE
    ) |>
    dplyr::mutate(
      PC1_plot = PC1 * loading_scale,
      PC2_plot = PC2 * loading_scale,
      Variable_label = dplyr::coalesce(
        unname(
          environmental_variable_text[
            Variable
          ]
        ),
        Variable
      )
    )

  figure_1 <- ggplot2::ggplot(
    pca_plot_data,
    ggplot2::aes(
      x = PC1,
      y = PC2
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "#BDBDBD",
      linewidth = 0.35
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      colour = "#BDBDBD",
      linewidth = 0.35
    ) +
    ggplot2::geom_segment(
      data = pca_segment_data,
      ggplot2::aes(
        xend = Centroid_PC1,
        yend = Centroid_PC2,
        colour = Reach_code
      ),
      linewidth = 0.45,
      alpha = 0.38,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = Reach_code,
        shape = Campaign
      ),
      size = 2.7,
      colour = "black",
      stroke = 0.55,
      alpha = 0.58
    ) +
    ggplot2::geom_point(
      data = pca_group_centroids,
      ggplot2::aes(
        x = Centroid_PC1,
        y = Centroid_PC2,
        fill = Reach_code,
        shape = Campaign
      ),
      size = 4.4,
      colour = "black",
      stroke = 0.85
    ) +
    ggrepel::geom_text_repel(
      data = pca_group_centroids,
      ggplot2::aes(
        x = Centroid_PC1,
        y = Centroid_PC2,
        label = Group_label,
        colour = Reach_code
      ),
      size = 3.0,
      fontface = "bold",
      show.legend = FALSE,
      min.segment.length = 0,
      max.overlaps = Inf
    ) +
    ggplot2::geom_segment(
      data = loading_plot_data,
      ggplot2::aes(
        x = 0,
        y = 0,
        xend = PC1_plot,
        yend = PC2_plot
      ),
      inherit.aes = FALSE,
      arrow = ggplot2::arrow(
        length = grid::unit(
          0.15,
          "cm"
        ),
        type = "closed"
      ),
      colour = "#4D4D4D",
      linewidth = 0.55
    ) +
    ggrepel::geom_text_repel(
      data = loading_plot_data,
      ggplot2::aes(
        x = PC1_plot,
        y = PC2_plot,
        label = Variable_label
      ),
      inherit.aes = FALSE,
      colour = "#333333",
      size = 2.8,
      fontface = "bold",
      box.padding = 0.25,
      point.padding = 0.15,
      max.overlaps = Inf
    ) +
    ggplot2::scale_fill_manual(
      values = reach_colours,
      labels = REACH_LABEL_LEVELS,
      name = "Reach"
    ) +
    ggplot2::scale_colour_manual(
      values = reach_colours,
      guide = "none"
    ) +
    ggplot2::scale_shape_manual(
      values = campaign_shapes,
      labels = c(
        "Wet campaign",
        "Dry campaign"
      ),
      name = "Campaign"
    ) +
    ggplot2::labs(
      x = sprintf(
        "PC1 (%.1f%%)",
        100 *
          pca_variance[
            1
          ]
      ),
      y = sprintf(
        "PC2 (%.1f%%)",
        100 *
          pca_variance[
            2
          ]
      )
    ) +
    ggplot2::coord_equal() +
    theme_publication(
      base_size = 9.6
    ) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          shape = 21,
          size = 3.4
        )
      ),
      shape = ggplot2::guide_legend(
        order = 2,
        override.aes = list(
          fill = "white",
          size = 3.4
        )
      )
    )

  save_ggplot(
    figure_1,
    FIGURE_DIR,
    "Figure_1_environmental_PCA",
    width = 7.5,
    height = 6.2
  )


  # ---------------------------------------------------------------------------
  # 4.8 Standardized group-level environmental heat map
  # ---------------------------------------------------------------------------

  transformed_group_means <- transformed_group_means |>
    dplyr::mutate(
      Biological_group_ID = as.character(
        Biological_group_ID
      )
    ) |>
    dplyr::arrange(
      match(
        Biological_group_ID,
        group_order
      )
    )

  group_environment_matrix <- as.matrix(
    transformed_group_means[
      ,
      retained_environmental_variables,
      drop = FALSE
    ]
  )

  rownames(
    group_environment_matrix
  ) <- transformed_group_means$Biological_group_ID

  group_environment_scaled <- scale(
    group_environment_matrix,
    center = TRUE,
    scale = TRUE
  )

  ordered_heatmap_variables <- unique(
    c(
      intersect(
        ENVIRONMENTAL_DOMAINS$Hydrology_Habitat,
        retained_environmental_variables
      ),

      intersect(
        ENVIRONMENTAL_DOMAINS$Thermal_Oxygen,
        retained_environmental_variables
      ),

      intersect(
        ENVIRONMENTAL_DOMAINS$Nutrient_Productivity,
        retained_environmental_variables
      ),

      setdiff(
        retained_environmental_variables,
        unlist(
          ENVIRONMENTAL_DOMAINS,
          use.names = FALSE
        )
      )
    )
  )

  group_environment_scaled <- group_environment_scaled[
    ,
    ordered_heatmap_variables,
    drop = FALSE
  ]

  heatmap_row_labels <- c(
    "HC__Wet_2023" = "S1 — Wet",
    "HC__Dry_2023" = "S1 — Dry",
    "QJ__Wet_2023" = "S2 — Wet",
    "QJ__Dry_2023" = "S2 — Dry",
    "SY__Wet_2023" = "S3 — Wet",
    "SY__Dry_2023" = "S3 — Dry"
  )

  heatmap_column_labels <- dplyr::coalesce(
    unname(
      environmental_variable_text[
        ordered_heatmap_variables
      ]
    ),
    ordered_heatmap_variables
  )

  heatmap_annotation <- data.frame(
    Reach = unname(
      REACH_LABELS[as.character(transformed_group_means$Reach_code)]
    ),
    Campaign = ifelse(
      transformed_group_means$Campaign == "Wet_2023",
      "Wet",
      "Dry"
    ),
    row.names = transformed_group_means$Biological_group_ID,
    stringsAsFactors = FALSE
  )

  heatmap_annotation_colours <- list(
    Reach = c(
      "S1" = "#0072B2",
      "S2" = "#D55E00",
      "S3" = "#009E73"
    ),

    Campaign = c(
      "Wet" = "#56B4E9",
      "Dry" = "#CC79A7"
    )
  )

  heatmap_colours <- grDevices::colorRampPalette(
    c(
      "#313695",
      "#74ADD1",
      "#F7F7F7",
      "#F46D43",
      "#A50026"
    )
  )(
    121
  )

  heatmap_limit <- max(
    2.5,
    ceiling(
      10 *
        max(
          abs(
            group_environment_scaled
          ),
          na.rm = TRUE
        )
    ) /
      10
  )

  heatmap_breaks <- seq(
    -heatmap_limit,
    heatmap_limit,
    length.out = 122
  )

  hydrology_count <- sum(
    ordered_heatmap_variables %in%
      ENVIRONMENTAL_DOMAINS$Hydrology_Habitat
  )

  thermal_count <- sum(
    ordered_heatmap_variables %in%
      ENVIRONMENTAL_DOMAINS$Thermal_Oxygen
  )

  heatmap_gaps <- c(
    hydrology_count,
    hydrology_count + thermal_count
  )

  heatmap_gaps <- heatmap_gaps[
    heatmap_gaps > 0 &
      heatmap_gaps < length(
        ordered_heatmap_variables
      )
  ]

  save_environmental_heatmap <- function(
    filename,
    device_type
  ) {

    if (device_type == "png") {

      grDevices::png(
        filename,
        width = 5600,
        height = 3500,
        res = 600,
        bg = "white"
      )
    }

    if (device_type == "tiff") {

      grDevices::tiff(
        filename,
        width = 5600,
        height = 3500,
        res = 600,
        compression = "lzw",
        bg = "white"
      )
    }

    if (device_type == "pdf") {

      grDevices::pdf(
        filename,
        width = 9.3,
        height = 5.8,
        useDingbats = FALSE,
        bg = "white"
      )
    }

    pheatmap::pheatmap(
      mat = group_environment_scaled,
      color = heatmap_colours,
      breaks = heatmap_breaks,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      gaps_col = heatmap_gaps,
      labels_row = unname(
        heatmap_row_labels[
          rownames(
            group_environment_scaled
          )
        ]
      ),
      labels_col = heatmap_column_labels,
      annotation_row = heatmap_annotation,
      annotation_colors = heatmap_annotation_colours,
      annotation_names_row = TRUE,
      border_color = "white",
      legend_breaks = c(
        -2,
        -1,
        0,
        1,
        2
      ),
      legend_labels = c(
        "−2",
        "−1",
        "0",
        "1",
        "2"
      ),
      fontsize = LEGEND_TEXT_SIZE,
      fontsize_row = AXIS_TEXT_SIZE,
      fontsize_col = AXIS_TEXT_SIZE,
      angle_col = 45,
      treeheight_row = 0,
      treeheight_col = 0,
      main = ""
    )

    grDevices::dev.off()
  }

  save_environmental_heatmap(
    file.path(
      FIGURE_DIR,
      "Figure_2_standardized_environmental_heatmap.png"
    ),
    "png"
  )

  save_environmental_heatmap(
    file.path(
      FIGURE_DIR,
      "Figure_2_standardized_environmental_heatmap.tiff"
    ),
    "tiff"
  )

  save_environmental_heatmap(
    file.path(
      FIGURE_DIR,
      "Figure_2_standardized_environmental_heatmap.pdf"
    ),
    "pdf"
  )


  # ---------------------------------------------------------------------------
  # 4.9 Environmental PERMANOVA
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # 4.9 Environmental PERMANOVA
  # ---------------------------------------------------------------------------

  environmental_distance <- stats::dist(
    environmental_scaled_matrix,
    method = "euclidean"
  )

  environmental_metadata <- transformed_environmental_data |>
    dplyr::select(
      Reach_code,
      Campaign,
      Transect_position
    ) |>
    as.data.frame(
      stringsAsFactors = FALSE
    )

  set.seed(
    SEED
  )

  overall_permanova <- vegan::adonis2(
    environmental_distance ~
      Reach_code *
      Campaign +
      Transect_position,
    data = environmental_metadata,
    permutations = N_PERMUTATIONS_RANDOM,
    by = "margin"
  )

  overall_permanova_table <- data.frame(
    Term = rownames(
      overall_permanova
    ),
    overall_permanova,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  campaign_blocks <- interaction(
    environmental_metadata$Reach_code,
    environmental_metadata$Transect_position,
    drop = TRUE
  )

  reach_blocks <- interaction(
    environmental_metadata$Campaign,
    environmental_metadata$Transect_position,
    drop = TRUE
  )

  transect_blocks <- interaction(
    environmental_metadata$Reach_code,
    environmental_metadata$Campaign,
    drop = TRUE
  )

  set.seed(
    SEED + 1L
  )

  campaign_blocked <- vegan::adonis2(
    environmental_distance ~ Campaign,
    data = environmental_metadata,
    permutations = N_PERMUTATIONS_RANDOM,
    strata = campaign_blocks
  )

  set.seed(
    SEED + 2L
  )

  reach_blocked <- vegan::adonis2(
    environmental_distance ~ Reach_code,
    data = environmental_metadata,
    permutations = N_PERMUTATIONS_RANDOM,
    strata = reach_blocks
  )

  set.seed(
    SEED + 3L
  )

  transect_blocked <- vegan::adonis2(
    environmental_distance ~ Transect_position,
    data = environmental_metadata,
    permutations = N_PERMUTATIONS_RANDOM,
    strata = transect_blocks
  )

  blocked_permanova_table <- dplyr::bind_rows(
    data.frame(
      Test = "Campaign within Reach × Transect",
      Term = rownames(
        campaign_blocked
      ),
      campaign_blocked,
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),

    data.frame(
      Test = "Reach within Campaign × Transect",
      Term = rownames(
        reach_blocked
      ),
      reach_blocked,
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),

    data.frame(
      Test = "Transect within Reach × Campaign",
      Term = rownames(
        transect_blocked
      ),
      transect_blocked,
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  )

  write.csv(
    overall_permanova_table,
    file.path(
      TABLE_DIR,
      "Table_3_environmental_PERMANOVA.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    blocked_permanova_table,
    file.path(
      TABLE_DIR,
      "Table_S3_blocked_environmental_PERMANOVA.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 4.10 Environmental multivariate dispersion
  # ---------------------------------------------------------------------------

  dispersion_groupings <- list(
    Reach = environmental_metadata$Reach_code,
    Campaign = environmental_metadata$Campaign,
    Reach_Campaign = interaction(
      environmental_metadata$Reach_code,
      environmental_metadata$Campaign,
      drop = TRUE
    )
  )

  dispersion_rows <- list()
  dispersion_index <- 1L

  for (grouping_name in names(
    dispersion_groupings
  )) {

    current_fit <- vegan::betadisper(
      environmental_distance,
      group = dispersion_groupings[[
        grouping_name
      ]],
      type = "median"
    )

    set.seed(
      SEED + 20L + dispersion_index
    )

    current_test <- vegan::permutest(
      current_fit,
      permutations = N_PERMUTATIONS_RANDOM
    )

    current_anova <- current_test$tab

    dispersion_rows[[
      dispersion_index
    ]] <- data.frame(
      Grouping = grouping_name,
      F_value = as.numeric(
        current_anova[
          1,
          "F"
        ]
      ),
      Permutation_P = as.numeric(
        current_anova[
          1,
          "Pr(>F)"
        ]
      ),
      Group_count = length(
        unique(
          dispersion_groupings[[
            grouping_name
          ]]
        )
      ),
      stringsAsFactors = FALSE
    )

    dispersion_index <- dispersion_index + 1L
  }

  environmental_dispersion_table <- dplyr::bind_rows(
    dispersion_rows
  )

  write.csv(
    environmental_dispersion_table,
    file.path(
      TABLE_DIR,
      "Table_4_environmental_PERMDISP.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 4.11 Univariate environmental models
  # ---------------------------------------------------------------------------

  old_contrasts <- options(
    "contrasts"
  )

  options(
    contrasts = c(
      "contr.sum",
      "contr.poly"
    )
  )

  on.exit(
    options(
      old_contrasts
    ),
    add = TRUE
  )

  univariate_rows <- list()
  univariate_index <- 1L

  for (variable_name in PRIMARY_ENVIRONMENTAL_VARIABLES) {

    current_data <- transformed_environmental_data |>
      dplyr::select(
        Reach_code,
        Campaign,
        Transect_position,
        dplyr::all_of(
          variable_name
        )
      )

    model_formula <- stats::as.formula(
      paste0(
        variable_name,
        " ~ Reach_code * Campaign + Transect_position"
      )
    )

    current_model <- stats::lm(
      model_formula,
      data = current_data
    )

    current_anova <- car::Anova(
      current_model,
      type = 3
    )

    current_anova_data <- data.frame(
      Term = rownames(
        current_anova
      ),
      current_anova,
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    current_anova_data <- current_anova_data |>
      dplyr::filter(
        Term %in% c(
          "Reach_code",
          "Campaign",
          "Transect_position",
          "Reach_code:Campaign"
        )
      )

    for (row_index in seq_len(
      nrow(
        current_anova_data
      )
    )) {

      univariate_rows[[
        univariate_index
      ]] <- data.frame(
        Variable = variable_name,
        Transformation = ifelse(
          variable_name %in% LOG1P_VARIABLES,
          "log1p",
          "None"
        ),
        Term = current_anova_data$Term[
          row_index
        ],
        Df = current_anova_data$Df[
          row_index
        ],
        F_value = current_anova_data$`F value`[
          row_index
        ],
        P_value = current_anova_data$`Pr(>F)`[
          row_index
        ],
        stringsAsFactors = FALSE
      )

      univariate_index <- univariate_index + 1L
    }
  }

  univariate_environmental_tests <- dplyr::bind_rows(
    univariate_rows
  ) |>
    dplyr::group_by(
      Term
    ) |>
    dplyr::mutate(
      FDR_adjusted_P = stats::p.adjust(
        P_value,
        method = "BH"
      )
    ) |>
    dplyr::ungroup()

  write.csv(
    univariate_environmental_tests,
    file.path(
      TABLE_DIR,
      "Table_5_environmental_univariate_tests.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 4.12 Environmental group summaries
  # ---------------------------------------------------------------------------

  environmental_group_summary <- environmental_data |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(
        PRIMARY_ENVIRONMENTAL_VARIABLES
      ),
      names_to = "Variable",
      values_to = "Value"
    ) |>
    dplyr::group_by(
      Reach_code,
      Reach_name,
      Campaign,
      Campaign_label,
      Variable
    ) |>
    dplyr::summarise(
      n = dplyr::n(),
      Mean = mean(
        Value
      ),
      SD = stats::sd(
        Value
      ),
      SE = standard_error(
        Value
      ),
      .groups = "drop"
    )

  write.csv(
    environmental_group_summary,
    file.path(
      TABLE_DIR,
      "Table_6_environmental_group_summary.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 4.13 Selected environmental variables figure
  # ---------------------------------------------------------------------------

  # Plain Unicode labels are used here. This prevents ggplot2 from
  # trying to parse ordinary text such as "Water temperature (°C)" as R code.
  selected_environmental_labels <- c(
    Water_temperature_C = "Water temperature (°C)",
    DO_mg_L = "Dissolved oxygen (mg L⁻¹)",
    Flow_velocity_m_s = "Flow velocity (m s⁻¹)",
    Water_TN_mg_L = "Total nitrogen (mg L⁻¹)",
    Water_TP_mg_L = "Total phosphorus (mg L⁻¹)",
    Chlorophyll_a_ug_L = "Chlorophyll-a (µg L⁻¹)"
  )

  selected_plot_summary <- environmental_group_summary |>
    dplyr::filter(
      Variable %in% SELECTED_ENVIRONMENTAL_PLOT_VARIABLES
    ) |>
    dplyr::mutate(
      Campaign = factor(
        Campaign,
        levels = c(
          "Wet_2023",
          "Dry_2023"
        ),
        labels = c(
          "Wet",
          "Dry"
        )
      ),

      Reach_code = factor(
        Reach_code,
        levels = c(
          "HC",
          "QJ",
          "SY"
        )
      ),

      Variable_label = factor(
        Variable,
        levels = SELECTED_ENVIRONMENTAL_PLOT_VARIABLES,
        labels = unname(
          selected_environmental_labels[
            SELECTED_ENVIRONMENTAL_PLOT_VARIABLES
          ]
        )
      )
    )

  selected_plot_raw <- environmental_data |>
    dplyr::select(
      Reach_code,
      Campaign,
      dplyr::all_of(
        SELECTED_ENVIRONMENTAL_PLOT_VARIABLES
      )
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(
        SELECTED_ENVIRONMENTAL_PLOT_VARIABLES
      ),
      names_to = "Variable",
      values_to = "Value"
    ) |>
    dplyr::mutate(
      Campaign = factor(
        Campaign,
        levels = c(
          "Wet_2023",
          "Dry_2023"
        ),
        labels = c(
          "Wet",
          "Dry"
        )
      ),

      Reach_code = factor(
        Reach_code,
        levels = c(
          "HC",
          "QJ",
          "SY"
        )
      ),

      Variable_label = factor(
        Variable,
        levels = SELECTED_ENVIRONMENTAL_PLOT_VARIABLES,
        labels = unname(
          selected_environmental_labels[
            SELECTED_ENVIRONMENTAL_PLOT_VARIABLES
          ]
        )
      )
    )

  environmental_dodge <- ggplot2::position_dodge(
    width = 0.58
  )

  environmental_jitter_dodge <- ggplot2::position_jitterdodge(
    jitter.width = 0.055,
    jitter.height = 0,
    dodge.width = 0.58,
    seed = SEED
  )

  figure_3 <- ggplot2::ggplot(
    selected_plot_summary,
    ggplot2::aes(
      x = Campaign,
      y = Mean,
      fill = Reach_code,
      shape = Reach_code,
      group = Reach_code
    )
  ) +
    ggplot2::geom_point(
      data = selected_plot_raw,
      ggplot2::aes(
        y = Value
      ),
      position = environmental_jitter_dodge,
      size = 1.8,
      colour = "black",
      stroke = 0.35,
      alpha = 0.28,
      inherit.aes = TRUE
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = Mean - SE,
        ymax = Mean + SE,
        colour = Reach_code
      ),
      position = environmental_dodge,
      width = 0.10,
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      position = environmental_dodge,
      size = 3.6,
      colour = "black",
      stroke = 0.70
    ) +
    ggplot2::facet_wrap(
      ~ Variable_label,
      scales = "free_y",
      ncol = 3,
      labeller = ggplot2::label_wrap_gen(
        width = 32
      )
    ) +
    ggplot2::scale_colour_manual(
      values = reach_colours,
      guide = "none"
    ) +
    ggplot2::scale_fill_manual(
      values = reach_colours,
      labels = REACH_LABEL_LEVELS,
      name = "Reach"
    ) +
    ggplot2::scale_shape_manual(
      values = reach_shapes,
      labels = REACH_LABEL_LEVELS,
      name = "Reach"
    ) +
    ggplot2::scale_x_discrete(
      expand = ggplot2::expansion(
        add = 0.45
      )
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0.06,
          0.12
        )
      )
    ) +
    ggplot2::labs(
      x = "Campaign",
      y = "Observed values and group mean ± SE"
    ) +
    theme_publication(
      base_size = 9.2
    ) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "top"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          size = 3.5
        )
      ),
      shape = "none"
    )

  save_ggplot(
    figure_3,
    FIGURE_DIR,
    "Figure_3_selected_environmental_variables",
    width = 8.2,
    height = 6.2
  )


  # ---------------------------------------------------------------------------
  # 4.14 Resource-isotope summaries
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # 4.14 Resource-isotope summaries
  # ---------------------------------------------------------------------------

  resource_campaign_audit <- resource_isotopes |>
    dplyr::count(
      Campaign,
      name = "Resource_records"
    ) |>
    dplyr::mutate(
      Included_in_matched_2023_analysis = Campaign %in%
        MATCHED_RESOURCE_CAMPAIGNS,
      Reason = dplyr::if_else(
        Included_in_matched_2023_analysis,
        "Matched to the 2023 fish campaigns",
        "Retained for provenance but excluded from matched 2023 summaries"
      )
    )

  write.csv(
    resource_campaign_audit,
    file.path(
      TABLE_DIR,
      "Table_7A_resource_campaign_audit.csv"
    ),
    row.names = FALSE
  )

  matched_resource_isotopes <- resource_isotopes |>
    dplyr::filter(
      Campaign %in% MATCHED_RESOURCE_CAMPAIGNS
    )

  if (nrow(matched_resource_isotopes) == 0L) {

    stop(
      "No resource-isotope records matched the 2023 campaigns.",
      call. = FALSE
    )
  }

  resource_isotope_summary <- matched_resource_isotopes |>
    dplyr::group_by(
      Reach_code,
      Reach,
      Campaign,
      Trophic_role,
      Resource_class
    ) |>
    dplyr::summarise(
      n = dplyr::n(),
      Mean_d13C = mean(
        d13C_permil
      ),
      SD_d13C = stats::sd(
        d13C_permil
      ),
      SE_d13C = standard_error(
        d13C_permil
      ),
      Mean_d15N = mean(
        d15N_permil
      ),
      SD_d15N = stats::sd(
        d15N_permil
      ),
      SE_d15N = standard_error(
        d15N_permil
      ),
      .groups = "drop"
    )

  write.csv(
    resource_isotope_summary,
    file.path(
      TABLE_DIR,
      "Table_7_resource_isotope_summary.csv"
    ),
    row.names = FALSE
  )

  trophic_role_colours <- c(
    "Producer" = "#009E73",
    "Basal material" = "#0072B2",
    "Primary consumer" = "#E69F00",
    "Consumer/prey" = "#CC79A7",
    "Other" = "#666666"
  )

  resource_class_shapes <- c(
    "Phytoplankton" = 16,
    "Plant" = 15,
    "Periphyton" = 17,
    "POM/detritus" = 18,
    "Gastropod" = 1,
    "Bivalve" = 0,
    "Zooplankton" = 2,
    "Invertebrate prey" = 5,
    "Other" = 4
  )

  resource_panel_levels <- c(
    "S1 — Wet",
    "S1 — Dry",
    "S2 — Wet",
    "S2 — Dry",
    "S3 — Wet",
    "S3 — Dry"
  )

  resource_plot_data <- matched_resource_isotopes |>
    dplyr::mutate(
      Trophic_role_plot = dplyr::if_else(
        is.na(
          Trophic_role
        ) |
          !Trophic_role %in% names(
            trophic_role_colours
          ),
        "Other",
        Trophic_role
      ),

      Resource_class_plot = dplyr::if_else(
        is.na(
          Resource_class
        ) |
          !Resource_class %in% names(
            resource_class_shapes
          ),
        "Other",
        Resource_class
      ),

      Panel = factor(
        paste(
          Reach,
          ifelse(
            Campaign == "Wet_2023",
            "Wet",
            "Dry"
          ),
          sep = " — "
        ),
        levels = resource_panel_levels
      )
    )

  resource_centroid_plot_data <- resource_isotope_summary |>
    dplyr::mutate(
      Trophic_role_plot = dplyr::if_else(
        is.na(
          Trophic_role
        ) |
          !Trophic_role %in% names(
            trophic_role_colours
          ),
        "Other",
        Trophic_role
      ),

      Resource_class_plot = dplyr::if_else(
        is.na(
          Resource_class
        ) |
          !Resource_class %in% names(
            resource_class_shapes
          ),
        "Other",
        Resource_class
      ),

      Panel = factor(
        paste(
          Reach,
          ifelse(
            Campaign == "Wet_2023",
            "Wet",
            "Dry"
          ),
          sep = " — "
        ),
        levels = resource_panel_levels
      )
    )

  fish_centroid_plot_data <- biological_groups |>
    dplyr::mutate(
      Sample_ID = as.character(
        Sample_ID
      ),

      Panel = factor(
        paste(
          Reach_name,
          ifelse(
            Campaign == "Wet_2023",
            "Wet",
            "Dry"
          ),
          sep = " — "
        ),
        levels = resource_panel_levels
      )
    )

  figure_4 <- ggplot2::ggplot(
    resource_plot_data,
    ggplot2::aes(
      x = d13C_permil,
      y = d15N_permil
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        colour = Trophic_role_plot,
        shape = Resource_class_plot
      ),
      size = 1.8,
      alpha = 0.35,
      stroke = 0.45
    ) +
    ggplot2::geom_segment(
      data = resource_centroid_plot_data,
      ggplot2::aes(
        x = Mean_d13C - SE_d13C,
        xend = Mean_d13C + SE_d13C,
        y = Mean_d15N,
        yend = Mean_d15N,
        colour = Trophic_role_plot
      ),
      inherit.aes = FALSE,
      linewidth = 0.65,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = resource_centroid_plot_data,
      ggplot2::aes(
        x = Mean_d13C,
        xend = Mean_d13C,
        y = Mean_d15N - SE_d15N,
        yend = Mean_d15N + SE_d15N,
        colour = Trophic_role_plot
      ),
      inherit.aes = FALSE,
      linewidth = 0.65,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = resource_centroid_plot_data,
      ggplot2::aes(
        x = Mean_d13C,
        y = Mean_d15N,
        colour = Trophic_role_plot,
        shape = Resource_class_plot
      ),
      inherit.aes = FALSE,
      size = 3.1,
      stroke = 0.95
    ) +
    ggplot2::geom_point(
      data = fish_centroid_plot_data,
      ggplot2::aes(
        x = Posterior_median_d13C,
        y = Posterior_median_d15N
      ),
      inherit.aes = FALSE,
      shape = 8,
      size = 4.0,
      stroke = 1.0,
      colour = "black"
    ) +
    ggplot2::facet_wrap(
      ~ Panel,
      ncol = 3,
      drop = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = trophic_role_colours,
      name = "Trophic role"
    ) +
    ggplot2::scale_shape_manual(
      values = resource_class_shapes,
      name = "Resource class"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0.05,
          0.08
        )
      )
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0.06,
          0.10
        )
      )
    ) +
    ggplot2::labs(
      x = expression(
        delta^13 *
          C ~ ("\u2030")
      ),
      y = expression(
        delta^15 *
          N ~ ("\u2030")
      )
    ) +
    theme_publication(
      base_size = 8.9
    ) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      legend.position = "top"
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        order = 1,
        nrow = 1,
        override.aes = list(
          shape = 16,
          size = 3,
          alpha = 1
        )
      ),
      shape = ggplot2::guide_legend(
        order = 2,
        nrow = 2,
        override.aes = list(
          colour = "#333333",
          size = 2.7,
          alpha = 1
        )
      )
    )

  save_ggplot(
    figure_4,
    FIGURE_DIR,
    "Figure_4_resource_isotope_space",
    width = 10.2,
    height = 7.1
  )


  # ---------------------------------------------------------------------------
  # 4.15 Gastropod-referenced TP and separate reference-sensitivity analysis
  # ---------------------------------------------------------------------------
  fish_tp_data <- fish_isotopes
  tp <- analyse_figure8_tp(fish_tp_data, resource_isotopes, tp_baselines)
  revised_tp_group_summary <- tp$groups
  revised_tp_seasonal_contrasts <- tp$contrasts
  revised_tp_pairwise <- tp$pairwise
  revised_tp_pairwise$Pair_ID <- vapply(seq_len(nrow(revised_tp_pairwise)), function(i) {
    x <- revised_tp_pairwise[i, ]
    hit <- biological_pairwise[(biological_pairwise$Sample_1 == x$Sample_1 &
      biological_pairwise$Sample_2 == x$Sample_2) |
      (biological_pairwise$Sample_1 == x$Sample_2 & biological_pairwise$Sample_2 == x$Sample_1), ]
    if (nrow(hit) != 1L) stop("TP pair must map to exactly one biological pair.")
    as.character(hit$Pair_ID)
  }, character(1))
  revised_tp_pairwise <- revised_tp_pairwise[order(revised_tp_pairwise$Pair_ID), ]
  tp_group_sample_sizes <- data.frame(Reach_code = tp$groups$Reach_code,
    Campaign = tp$groups$Campaign, Fish_n = tp$groups$Fish_n)
  tp_inclusion_audit <- data.frame(
    Check = c("Fish available", "Fish retained in primary TP analysis",
      "Fish excluded from primary TP analysis", "Outcome-dependent TP filters applied"),
    Value = c(nrow(fish_isotopes), nrow(fish_tp_data), 0L, 0L),
    Primary_analysis_rule = c("Expected 609", "Expected 609", "Expected 0", "Prohibited"))
  tp_exports <- list(
    Table_S4a_TP_inclusion_audit = tp_inclusion_audit,
    Table_S4b_TP_group_sample_sizes_all_609 = tp_group_sample_sizes,
    Table_S4_primary_consumer_baseline_values = tp$baseline_raw,
    Table_8_revised_trophic_position_group_summary = revised_tp_group_summary,
    Table_9_revised_trophic_position_seasonal_contrasts = revised_tp_seasonal_contrasts,
    Table_S5_revised_trophic_position_pairwise = revised_tp_pairwise,
    Table_S6_revised_TP_TDF_sensitivity = tp$tdf_groups,
    Table_S6b_TP_TDF_contrast_sensitivity = tp$tdf_contrasts,
    Table_S8d_I_gastropod_TP_changes = tp$s8d_I,
    Table_S8d_II_gastropod_TP_groups = tp$s8d_II,
    Table_S8g_reference_sensitivity = tp$reference_sensitivity,
    Table_S8g_consumer_availability = tp$inventory)
  for (output_name in names(tp_exports)) {
    write.csv(tp_exports[[output_name]], file.path(TABLE_DIR, paste0(output_name, ".csv")),
      row.names = FALSE, fileEncoding = "UTF-8")
  }

  # ---------------------------------------------------------------------------
  # 4.16 Revised trophic-position figure
  # ---------------------------------------------------------------------------

  tp_plot_data <- revised_tp_group_summary |>
    dplyr::mutate(
      Campaign = factor(
        Campaign,
        levels = c(
          "Wet_2023",
          "Dry_2023"
        ),
        labels = c(
          "Wet",
          "Dry"
        )
      ),

      Reach_code = factor(
        Reach_code,
        levels = c(
          "HC",
          "QJ",
          "SY"
        )
      )
    )

  figure_5 <- ggplot2::ggplot(
    tp_plot_data,
    ggplot2::aes(
      x = Campaign,
      y = Point_TP,
      group = Reach_code,
      colour = Reach_code
    )
  ) +
    ggplot2::geom_hline(
      yintercept = TROPHIC_LEVEL_BASELINE,
      linetype = "dashed",
      colour = "#777777",
      linewidth = 0.45
    ) +
    ggplot2::geom_line(
      linewidth = 0.85,
      alpha = 0.90
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = Lower_95_TP,
        ymax = Upper_95_TP
      ),
      width = 0.08,
      linewidth = 0.65
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = Reach_code
      ),
      shape = 21,
      size = 4.1,
      colour = "black",
      stroke = 0.75
    ) +
    ggplot2::scale_colour_manual(
      values = reach_colours,
      labels = REACH_LABEL_LEVELS,
      name = "Reach"
    ) +
    ggplot2::scale_fill_manual(
      values = reach_colours,
      guide = "none"
    ) +
    ggplot2::scale_x_discrete(
      expand = ggplot2::expansion(
        add = 0.35
      )
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0.06,
          0.14
        )
      )
    ) +
    ggplot2::labs(
      x = "Campaign",
      y = "Gastropod-referenced fish TP"
    ) +
    theme_publication(
      base_size = 10
    ) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          linewidth = 1.0
        )
      )
    )

  save_ggplot(
    figure_5,
    FIGURE_DIR,
    "Figure_5_revised_trophic_position",
    width = 6.8,
    height = 5.2
  )


  # ---------------------------------------------------------------------------
  # 4.17 Exact environmental–biological matrix tests
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # 4.17 Exact environmental–biological matrix tests
  # ---------------------------------------------------------------------------

  # Reviewer 3: no unrestricted group-label permutations are used here.
  r3_designs <- r3_make_designs(group_order)
  exact_permutation_matrix <- r3_nonidentity(r3_designs[[R3_DISPLAY_SCHEME]])
  r3_write_design_files(TABLE_DIR, r3_designs)

  # Overall environmental distance.
  overall_group_environment_matrix <- group_environment_scaled

  environmental_distance_matrices <- list(
    Overall_Environment = as.matrix(
      stats::dist(
        overall_group_environment_matrix,
        method = "euclidean"
      )
    )
  )

  domain_screening_rows <- list()
  domain_screening_index <- 1L

  for (domain_name in names(
    ENVIRONMENTAL_DOMAINS
  )) {

    domain_candidates <- ENVIRONMENTAL_DOMAINS[[
      domain_name
    ]]

    domain_screening <- screen_correlated_variables(
      transformed_environmental_data,
      domain_candidates,
      cutoff = CORRELATION_CUTOFF
    )

    domain_variables <- domain_screening$retained

    if (length(
      domain_variables
    ) < 1L) {

      next
    }

    domain_group_matrix <- as.matrix(
      transformed_group_means[
        ,
        domain_variables,
        drop = FALSE
      ]
    )

    rownames(
      domain_group_matrix
    ) <- as.character(
      transformed_group_means$Biological_group_ID
    )

    domain_group_scaled <- scale(
      domain_group_matrix,
      center = TRUE,
      scale = TRUE
    )

    environmental_distance_matrices[[
      domain_name
    ]] <- as.matrix(
      stats::dist(
        domain_group_scaled,
        method = "euclidean"
      )
    )

    current_screening <- domain_screening$screening |>
      dplyr::mutate(
        Domain = domain_name
      )

    domain_screening_rows[[
      domain_screening_index
    ]] <- current_screening

    domain_screening_index <- domain_screening_index + 1L
  }

  domain_screening_table <- dplyr::bind_rows(
    domain_screening_rows
  )

  write.csv(
    domain_screening_table,
    file.path(
      TABLE_DIR,
      "Table_S7_environmental_domain_screening.csv"
    ),
    row.names = FALSE
  )

  biological_distance_matrices <- list(
    Bray_Curtis = build_symmetric_matrix(
      biological_pairwise,
      "Bray_Curtis",
      group_order
    ),

    Jaccard = build_symmetric_matrix(
      biological_pairwise,
      "Jaccard",
      group_order
    ),

    Centroid_Distance = build_symmetric_matrix(
      biological_pairwise,
      "Median_centroid_distance",
      group_order
    ),

    Absolute_SEA_Difference = build_symmetric_matrix(
      biological_pairwise,
      "Absolute_SEA_difference",
      group_order
    ),

    Absolute_Log_SEA_Ratio = build_symmetric_matrix(
      biological_pairwise,
      "Absolute_log_SEA_ratio",
      group_order
    ),

    Niche_Nonoverlap = build_symmetric_matrix(
      biological_pairwise,
      "Niche_nonoverlap_proportion",
      group_order
    ),

    Revised_TP_Difference = build_symmetric_matrix(
      revised_tp_pairwise,
      "Absolute_point_TP_difference",
      group_order
    )
  )

  paired_environmental_biological_results <- r3_matrix_tests(
    environmental_distance_matrices, biological_distance_matrices, "spearman", r3_designs,
    labels = c("Environmental_metric", "Biological_metric"))
  names(paired_environmental_biological_results)[
    names(paired_environmental_biological_results) == "Matrix_correlation"] <- "Mantel_r"
  exact_environmental_biological_results <- subset(paired_environmental_biological_results,
    Permutation_scheme == R3_DISPLAY_SCHEME)
  write.csv(paired_environmental_biological_results, file.path(TABLE_DIR,
    "Reviewer3_all_paired_environmental_Mantel.csv"), row.names = FALSE)

  overall_environmental_linkage <- exact_environmental_biological_results |>
    dplyr::filter(
      Environmental_metric == "Overall_Environment"
    ) |>
    dplyr::mutate(
      FDR_within_overall_tests = stats::p.adjust(
        Exact_one_sided_P,
        method = "BH"
      )
    )

  domain_environmental_linkage <- exact_environmental_biological_results |>
    dplyr::filter(
      Environmental_metric != "Overall_Environment"
    )

  write.csv(
    overall_environmental_linkage,
    file.path(
      TABLE_DIR,
      "Table_10_exact_environmental_biological_Mantel_results.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    domain_environmental_linkage,
    file.path(
      TABLE_DIR,
      "Table_S8_domain_specific_Mantel_results.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 4.18 Exact environmental–biological Mantel heat-map figure
  # ---------------------------------------------------------------------------

  biological_metric_labels <- c(
    Bray_Curtis = "Bray–Curtis community",
    Jaccard = "Jaccard community",
    Centroid_Distance = "Isotope-centroid distance",
    Absolute_SEA_Difference = "|SEA difference|",
    Absolute_Log_SEA_Ratio = "|log SEA ratio|",
    Niche_Nonoverlap = "Niche non-overlap",
    Revised_TP_Difference = "Gastropod-referenced TP difference"
  )

  environmental_metric_labels <- c(
    Overall_Environment = "Overall\nenvironment",
    Hydrology_Habitat = "Hydrology–\nhabitat",
    Thermal_Oxygen = "Thermal–\noxygen",
    Nutrient_Productivity = "Nutrient–\nproductivity"
  )

  mantel_plot_data <- exact_environmental_biological_results |>
    dplyr::mutate(
      Environmental_label = factor(
        environmental_metric_labels[
          Environmental_metric
        ],
        levels = unname(
          environmental_metric_labels
        )
      ),

      Biological_label = factor(
        biological_metric_labels[
          Biological_metric
        ],
        levels = rev(
          unname(
            biological_metric_labels
          )
        )
      ),

      Annotation = if (R3_SHOW_P_IN_FIGURES) sprintf(
        "\u03c1 = %.2f\nP(pair) = %.3f", Mantel_r, Exact_one_sided_P
      ) else sprintf("\u03c1 = %.2f", Mantel_r),

      Text_colour = ifelse(
        abs(
          Mantel_r
        ) >= 0.58,
        "white",
        "black"
      )
    )

  figure_6 <- ggplot2::ggplot(
    mantel_plot_data,
    ggplot2::aes(
      x = Environmental_label,
      y = Biological_label,
      fill = Mantel_r
    )
  ) +
    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.85
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = Annotation,
        colour = Text_colour
      ),
      size = 3.0,
      fontface = "bold",
      lineheight = 0.95,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0,
      limits = c(
        -1,
        1
      ),
      breaks = c(
        -1,
        -0.5,
        0,
        0.5,
        1
      ),
      name = "Spearman\nMantel \u03c1"
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL
    ) +
    ggplot2::coord_fixed(
      ratio = 0.82
    ) +
    theme_publication(
      base_size = 9.2
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        face = "bold",
        lineheight = 0.95
      ),
      axis.text.y = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "top"
    )

  save_ggplot(
    figure_6,
    FIGURE_DIR,
    "Figure_6_environmental_biological_Mantel_heatmap",
    width = 8.2,
    height = 6.6
  )

  save_combined_figure8(figure_5, figure_6, FIGURE_DIR)



  # ---------------------------------------------------------------------------
  # 4.19 Procrustes concordance
  # ---------------------------------------------------------------------------

  group_environment_pca <- stats::prcomp(
    group_environment_scaled,
    center = FALSE,
    scale. = FALSE
  )

  environmental_group_scores_2d <- group_environment_pca$x[
    ,
    1:2,
    drop = FALSE
  ]

  rownames(
    environmental_group_scores_2d
  ) <- group_order

  bray_pcoa <- stats::cmdscale(
    stats::as.dist(
      biological_distance_matrices$Bray_Curtis
    ),
    k = 2,
    eig = TRUE,
    add = TRUE
  )

  jaccard_pcoa <- stats::cmdscale(
    stats::as.dist(
      biological_distance_matrices$Jaccard
    ),
    k = 2,
    eig = TRUE,
    add = TRUE
  )

  isotope_centroid_coordinates <- biological_groups |>
    dplyr::mutate(
      Sample_ID = as.character(
        Sample_ID
      )
    ) |>
    dplyr::arrange(
      match(
        Sample_ID,
        group_order
      )
    ) |>
    dplyr::select(
      Posterior_median_d13C,
      Posterior_median_d15N
    ) |>
    as.matrix()

  rownames(
    isotope_centroid_coordinates
  ) <- group_order

  procrustes_targets <- list(
    Bray_Curtis_community = bray_pcoa$points,
    Jaccard_community = jaccard_pcoa$points,
    Isotope_centroids = isotope_centroid_coordinates
  )

  procrustes_summary_rows <- list()
  procrustes_plot_rows <- list()

  procrustes_index <- 1L

  for (comparison_name in names(
    procrustes_targets
  )) {

    current_target <- procrustes_targets[[
      comparison_name
    ]]

    current_target <- current_target[
      group_order,
      ,
      drop = FALSE
    ]

    current_procrustes <- vegan::procrustes(
      X = environmental_group_scores_2d,
      Y = current_target,
      symmetric = TRUE
    )

    set.seed(
      SEED + 300L + procrustes_index
    )

    current_protest <- vegan::protest(
      X = environmental_group_scores_2d,
      Y = current_target,
      permutations = exact_permutation_matrix,
      symmetric = TRUE
    )

    procrustes_summary_rows[[
      procrustes_index
    ]] <- data.frame(
      Comparison = comparison_name,
      Procrustes_correlation = as.numeric(
        current_protest$t0
      ),
      Exact_P = as.numeric(
        current_protest$signif
      ),
      Permutations = nrow(
        exact_permutation_matrix
      ),
      stringsAsFactors = FALSE
    )

    procrustes_plot_rows[[
      procrustes_index
    ]] <- data.frame(
      Sample_ID = group_order,
      Environmental_axis_1 = current_procrustes$X[
        ,
        1
      ],
      Environmental_axis_2 = current_procrustes$X[
        ,
        2
      ],
      Biological_axis_1 = current_procrustes$Yrot[
        ,
        1
      ],
      Biological_axis_2 = current_procrustes$Yrot[
        ,
        2
      ],
      Comparison = comparison_name,
      Reach_code = sub(
        "__.*$",
        "",
        group_order
      ),
      Campaign = ifelse(
        grepl(
          "Wet",
          group_order
        ),
        "Wet_2023",
        "Dry_2023"
      ),
      stringsAsFactors = FALSE
    )

    procrustes_index <- procrustes_index + 1L
  }

  procrustes_summary <- dplyr::bind_rows(procrustes_summary_rows)
  procrustes_summary$Permutation_scheme <- R3_DISPLAY_SCHEME
  procrustes_summary$Arrangements_including_observed <- nrow(r3_designs[[R3_DISPLAY_SCHEME]])
  procrustes_summary$Inference_status <- "Conditional exploratory sensitivity"
  procrustes_secondary <- lapply(names(procrustes_targets), function(nm) {
    target <- procrustes_targets[[nm]][group_order, , drop = FALSE]
    fit <- vegan::protest(environmental_group_scores_2d, target,
      permutations = r3_nonidentity(r3_designs[["Within_reach_swaps"]]), symmetric = TRUE)
    data.frame(Comparison = nm, Procrustes_correlation = as.numeric(fit$t0),
      Exact_P = as.numeric(fit$signif), Permutations = 7L,
      Permutation_scheme = "Within_reach_swaps", Arrangements_including_observed = 8L,
      Inference_status = "Conditional exploratory sensitivity", stringsAsFactors = FALSE)
  })
  procrustes_paired_results <- dplyr::bind_rows(procrustes_summary,
    dplyr::bind_rows(procrustes_secondary)) |>
    dplyr::group_by(Permutation_scheme) |>
    dplyr::mutate(FDR_adjusted_P = stats::p.adjust(Exact_P, "BH")) |>
    dplyr::ungroup()
  write.csv(procrustes_paired_results, file.path(TABLE_DIR,
    "Reviewer3_paired_Procrustes_sensitivities.csv"), row.names = FALSE)

  procrustes_plot_data <- dplyr::bind_rows(
    procrustes_plot_rows
  )

  write.csv(
    procrustes_summary,
    file.path(
      TABLE_DIR,
      "Table_11_Procrustes_concordance.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    procrustes_plot_data,
    file.path(
      TABLE_DIR,
      "Table_S9_Procrustes_coordinates.csv"
    ),
    row.names = FALSE
  )

  procrustes_facet_labels <- procrustes_summary |>
    dplyr::mutate(
      Comparison_title = dplyr::recode(
        Comparison,
        Bray_Curtis_community = "Environment vs Bray–Curtis community",
        Jaccard_community = "Environment vs Jaccard community",
        Isotope_centroids = "Environment vs isotope centroids"
      ),

      Facet_label = if (R3_SHOW_P_IN_FIGURES) sprintf(
        "%s\nr = %.2f, P(pair) = %.3f", Comparison_title, Procrustes_correlation, Exact_P
      ) else sprintf("%s\nr = %.2f", Comparison_title, Procrustes_correlation)
    ) |>
    dplyr::select(
      Comparison,
      Facet_label
    )

  procrustes_plot_data <- procrustes_plot_data |>
    dplyr::left_join(
      procrustes_facet_labels,
      by = "Comparison"
    ) |>
    dplyr::mutate(
      Reach_code = factor(
        Reach_code,
        levels = c(
          "HC",
          "QJ",
          "SY"
        )
      ),

      Campaign = factor(
        Campaign,
        levels = c(
          "Wet_2023",
          "Dry_2023"
        )
      ),

      Sample_label = dplyr::recode(
        Sample_ID,
        HC__Wet_2023 = "S1–Wet",
        HC__Dry_2023 = "S1–Dry",
        QJ__Wet_2023 = "S2–Wet",
        QJ__Dry_2023 = "S2–Dry",
        SY__Wet_2023 = "S3–Wet",
        SY__Dry_2023 = "S3–Dry"
      )
    )

  figure_7 <- ggplot2::ggplot(
    procrustes_plot_data
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "#D9D9D9",
      linewidth = 0.35
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      colour = "#D9D9D9",
      linewidth = 0.35
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = Environmental_axis_1,
        y = Environmental_axis_2,
        xend = Biological_axis_1,
        yend = Biological_axis_2,
        colour = Reach_code,
        linetype = Campaign
      ),
      arrow = ggplot2::arrow(
        length = grid::unit(
          0.14,
          "cm"
        ),
        type = "closed"
      ),
      linewidth = 0.75,
      alpha = 0.82
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = Environmental_axis_1,
        y = Environmental_axis_2,
        fill = Reach_code,
        shape = Campaign
      ),
      size = 3.7,
      colour = "black",
      stroke = 0.70
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = Biological_axis_1,
        y = Biological_axis_2,
        colour = Reach_code
      ),
      shape = 4,
      size = 3.4,
      stroke = 1.05
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(
        x = Biological_axis_1,
        y = Biological_axis_2,
        label = Sample_label,
        colour = Reach_code
      ),
      size = 2.7,
      fontface = "bold",
      min.segment.length = 0,
      show.legend = FALSE,
      max.overlaps = Inf
    ) +
    ggplot2::facet_wrap(
      ~ Facet_label,
      ncol = 3
    ) +
    ggplot2::scale_colour_manual(
      values = reach_colours,
      labels = REACH_LABEL_LEVELS,
      name = "Reach"
    ) +
    ggplot2::scale_fill_manual(
      values = reach_colours,
      guide = "none"
    ) +
    ggplot2::scale_shape_manual(
      values = campaign_shapes,
      labels = c(
        "Wet campaign",
        "Dry campaign"
      ),
      name = "Campaign"
    ) +
    ggplot2::scale_linetype_manual(
      values = c(
        "Wet_2023" = "solid",
        "Dry_2023" = "dashed"
      ),
      labels = c(
        "Wet campaign",
        "Dry campaign"
      ),
      name = "Campaign"
    ) +
    ggplot2::labs(
      x = "Procrustes axis 1",
      y = "Procrustes axis 2"
    ) +
    ggplot2::coord_equal() +
    theme_publication(
      base_size = 8.7
    ) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      legend.position = "top",
      strip.text = ggplot2::element_text(
        size = FACET_TEXT_SIZE,
        lineheight = 0.95
      )
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          linetype = 0,
          shape = 16,
          size = 3
        )
      ),
      shape = ggplot2::guide_legend(
        order = 2,
        override.aes = list(
          fill = "white",
          size = 3.2
        )
      ),
      linetype = "none"
    )

  save_ggplot(
    figure_7,
    FIGURE_DIR,
    "Figure_7_Procrustes_environment_biology",
    width = 11.2,
    height = 5.2
  )


  # ---------------------------------------------------------------------------
  # Manuscript-ready figure captions
  # ---------------------------------------------------------------------------

  figure_captions <- data.frame(
    Figure = c(
      "Figure 1",
      "Figure 2",
      "Figure 3",
      "Figure 4",
      "Figure 5",
      "Figure 6",
      "Figure 7"
    ),

    Filename = c(
      "Figure_1_environmental_PCA",
      "Figure_2_standardized_environmental_heatmap",
      "Figure_3_selected_environmental_variables",
      "Figure_4_resource_isotope_space",
      "Figure_5_revised_trophic_position",
      "Figure_6_environmental_biological_Mantel_heatmap",
      "Figure_7_Procrustes_environment_biology"
    ),

    Caption = c(
      paste(
        "Principal component analysis of matched environmental observations",
        "from three transect positions in each reach and campaign.",
        "Small symbols represent individual transects, large symbols represent",
        "reach–campaign centroids, and arrows show the strongest environmental loadings."
      ),

      paste(
        "Standardized environmental conditions for the six reach–campaign groups.",
        "Variables are ordered by hydrology–habitat, thermal–oxygen, and",
        "nutrient–productivity domains. Values are column-wise z scores."
      ),

      paste(
        "Spatial and seasonal variation in selected environmental indicators.",
        "Faint symbols represent individual transects and large symbols show",
        "reach–campaign means ± SE."
      ),

      paste(
        "Stable-isotope distributions of basal resources and primary consumers",
        "by reach and campaign. Small symbols are individual observations,",
        "large symbols and crosshairs are resource-class means ± SE, and black",
        "stars indicate posterior median fish-community isotope centroids."
      ),

      paste(
        "Gastropod-referenced community-mean fish trophic position using observed",
        "same-reach and same-campaign S. quadrata baselines and fixed TDF 3.4 per mil.",
        "Points are calculated from observed means; bars are 95% percentile bootstrap intervals.",
        "Estimates are conditional on baseline pathway representation and tissue integration."
      ),

      paste(
        "Exact Spearman Mantel associations between environmental dissimilarity",
        "and fish-community or isotope-based response distances.",
        "Cells report observed correlations from six summaries in three paired reaches.",
        "Two conditional restricted-null sensitivities are reported separately;",
        "non-significance does not demonstrate complementarity or independence."
      ),

      paste(
        "Procrustes concordance between environmental and biological configurations.",
        "Filled symbols denote environmental configurations, crosses denote rotated",
        "biological configurations, and arrows represent residual mismatch.",
        "Facet headings report descriptive Procrustes correlations; conditional paired-null P values are tabulated separately."
      )
    ),

    stringsAsFactors = FALSE
  )

  write.csv(
    figure_captions,
    file.path(
      TABLE_DIR,
      "Table_S11_manuscript_figure_captions.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 4.20 Matrix validation
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # 4.20 Matrix validation
  # ---------------------------------------------------------------------------

  matrix_validation_rows <- list()
  matrix_validation_index <- 1L

  all_matrix_lists <- c(
    environmental_distance_matrices,
    biological_distance_matrices
  )

  for (matrix_name in names(
    all_matrix_lists
  )) {

    current_matrix <- all_matrix_lists[[
      matrix_name
    ]]

    matrix_validation_rows[[
      matrix_validation_index
    ]] <- data.frame(
      Matrix = matrix_name,
      Rows = nrow(
        current_matrix
      ),
      Columns = ncol(
        current_matrix
      ),
      Maximum_symmetry_error = max(
        abs(
          current_matrix -
            t(
              current_matrix
            )
        )
      ),
      Maximum_absolute_diagonal = max(
        abs(
          diag(
            current_matrix
          )
        )
      ),
      Minimum_value = min(
        current_matrix
      ),
      Maximum_value = max(
        current_matrix
      ),
      stringsAsFactors = FALSE
    )

    matrix_validation_index <- matrix_validation_index + 1L
  }

  matrix_validation_table <- dplyr::bind_rows(
    matrix_validation_rows
  )

  if (
    any(
      matrix_validation_table$Maximum_symmetry_error > 1e-10
    ) ||
      any(
        matrix_validation_table$Maximum_absolute_diagonal > 1e-10
      )
  ) {

    stop(
      "At least one environmental or biological distance matrix failed validation.",
      call. = FALSE
    )
  }

  write.csv(
    matrix_validation_table,
    file.path(
      TABLE_DIR,
      "Table_S10_matrix_validation.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 4.21 Design and inference limitations
  # ---------------------------------------------------------------------------

  design_limitations <- data.frame(
    Item = c(
      "Environmental survey mapping",
      "Environmental replication",
      "Biological group count",
      "Environmental PCA",
      "Environmental PERMANOVA",
      "Matrix tests",
      "Matrix-test interpretation",
      "Revised trophic-position method",
      "Primary TP sample size",
      "TP exclusion rule",
      "Historical trophic-level and QC fields",
      "Baseline treatment",
      "Fish-level environmental regression",
      "Resource isotope analysis"
    ),

    Value = c(
      "20–27 June 2023 = Wet_2023; 14–23 November 2023 = Dry_2023",
      "Three transects per Reach × Campaign group; 18 observations",
      "Six aggregated Reach × Campaign groups",
      "Transect-level, transformed, standardized, and correlation-screened with DO priority-retained",
      "Reach × Campaign plus transect position; supplemented by blocked tests",
      "Observed Spearman correlations; whole reach pairs (6 arrangements) and within-reach swaps (8 arrangements)",
      "Exploratory comparison of three paired reaches; no demonstrated complementarity or independence",
      "Observed-mean gastropod-referenced TP; fixed TDF 3.4; 5000 within-group bootstrap replicates",
      "All 609 fish retained",
      "No exclusion based on calculated TP, including TP < 1",
      "Retained only for provenance and never used for primary-analysis inclusion",
      "Six observed same-reach same-campaign S. quadrata baselines; cross-season proxies prohibited",
      "Not performed to avoid pseudoreplication",
      "Matched 2023 campaigns only; June 2024 records retained for provenance but excluded"
    ),

    stringsAsFactors = FALSE
  )

  write.csv(
    design_limitations,
    file.path(
      TABLE_DIR,
      "Table_12_design_and_inference_limitations.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  design_limitations <- rbind(design_limitations, data.frame(
    Item = c("Primary fixed discrimination factor", "Reference-choice sensitivity",
      "Tissue integration", "TP matrix distance"),
    Value = c("3.4 per mil; fixed-factor checks at 2.9, 3.4 and 3.9 per mil",
      "Paired bulk zooplankton at QJ/SY is exploratory: taxonomy and feeding status unverified",
      "Fish muscle and gastropod foot muscle incorporation windows were not calibrated",
      "Absolute difference of observed-mean point estimates, consistent with Figure 8a"),
    stringsAsFactors = FALSE))
  write.csv(design_limitations,
    file.path(TABLE_DIR, "Table_12_design_and_inference_limitations.csv"), row.names = FALSE)
  writeLines(c(
    "Figure 8. Gastropod-referenced trophic position and environmental-biota concordance.",
    "(a) Points use observed fish and S. quadrata means at fixed TDF 3.4 per mil; bars are 95% percentile bootstrap intervals from 5000 within-group resamples.",
    "(b) Observed Spearman matrix correlations among six reach-campaign summaries. The TP distance uses the same points as panel (a). Conditional whole-reach-pair and within-reach-swap sensitivities are reported separately; the main panel does not establish indicator complementarity.",
    "TP estimates are conditional on the gastropod pathway reference and uncalibrated tissue integration; bulk-zooplankton diagnostics are reported separately in Table S8g."),
    file.path(TABLE_DIR, "Figure_8_caption.txt"))

  # 4.22 Combined Excel results workbook
  # ---------------------------------------------------------------------------

  writexl::write_xlsx(
    list(
      Environmental_Variable_Screening = variable_screening_table,
      Environmental_PCA_Variance = pca_variance_table,
      Environmental_PERMANOVA = overall_permanova_table,
      Environmental_Blocked_PERMANOVA = blocked_permanova_table,
      Environmental_PERMDISP = environmental_dispersion_table,
      Environmental_Univariate_Tests = univariate_environmental_tests,
      Environmental_Group_Summary = environmental_group_summary,
      Dissolved_Oxygen_Audit = dissolved_oxygen_audit,
      Resource_Campaign_Audit = resource_campaign_audit,
      Resource_Isotope_Summary = resource_isotope_summary,
      Revised_TP_Group_Summary = revised_tp_group_summary,
      Revised_TP_Seasonal_Contrasts = revised_tp_seasonal_contrasts,
      Revised_TP_Pairwise = revised_tp_pairwise,
      Table_S8d_I = tp$s8d_I,
      Table_S8d_II = tp$s8d_II,
      Table_S8g_Reference_Sensitivity = tp$reference_sensitivity,
      Consumer_Reference_Availability = tp$inventory,
      Baseline_Raw_Validation = tp$baseline_raw,
      Fixed_TDF_Group_Sensitivity = tp$tdf_groups,
      Fixed_TDF_Contrast_Sensitivity = tp$tdf_contrasts,
      TP_Inclusion_Audit = tp_inclusion_audit,
      TP_Group_Sample_Sizes = tp_group_sample_sizes,
      Environmental_Biological_Mantel = overall_environmental_linkage,
      Domain_Mantel = domain_environmental_linkage,
      Paired_Mantel_All = paired_environmental_biological_results,
      Permutation_Design = r3_permutation_table(r3_designs),
      Null_Hypotheses = r3_null_hypotheses(),
      Paired_Procrustes = procrustes_paired_results,
      Procrustes_Concordance = procrustes_summary,
      Matrix_Validation = matrix_validation_table,
      Design_Limitations = design_limitations
    ),
    path = file.path(
      RESULTS_DIR,
      "Han_River_Figure8_results_Reviewer3_ALL_609.xlsx"
    )
  )


  # ---------------------------------------------------------------------------
  # 4.23 Session information and completion log
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
    RUNNING_FILE
  )

  writeLines(
    c(
      "Environmental integration and revised trophic-position analysis completed successfully.",
      paste(
        "Completed:",
        Sys.time()
      ),
      paste(
        "Input workbook:",
        INPUT_FILE
      ),
      paste(
        "Environmental transects:",
        nrow(
          environmental_data
        )
      ),
      paste(
        "Retained environmental variables:",
        length(
          retained_environmental_variables
        )
      ),
      paste(
        "Biological groups:",
        length(
          group_order
        )
      ),
      paste(
        "Exact matrix permutations:",
        nrow(
          exact_permutation_matrix
        )
      ),
      paste(
        "Fish retained in primary TP analysis:",
        nrow(
          fish_tp_data
        )
      ),
      paste(
        "Fish excluded from primary TP analysis:",
        nrow(
          fish_isotopes
        ) - nrow(
          fish_tp_data
        )
      ),
      paste(
        "TP fixed-factor bootstrap replicates:",
        N_TP_DRAWS
      ),
      paste(
        "Results folder:",
        RESULTS_DIR
      )
    ),
    COMPLETED_FILE
  )

  cat(
    "\n============================================================\n",
    "ENVIRONMENTAL INTEGRATION AND REVISED TP ANALYSIS COMPLETED — ECOLOGICAL INDICATORS FIGURES\n",
    "============================================================\n",
    "\nResults folder:\n",
    RESULTS_DIR,
    "\n\nMain tables:\n",
    "Table_1_environmental_variable_screening.csv\n",
    "Table_2_environmental_PCA_variance.csv\n",
    "Table_3_environmental_PERMANOVA.csv\n",
    "Table_4_environmental_PERMDISP.csv\n",
    "Table_5_environmental_univariate_tests.csv\n",
    "Table_6_environmental_group_summary.csv\n",
    "Table_7_resource_isotope_summary.csv\n",
    "Table_8_revised_trophic_position_group_summary.csv\n",
    "Table_9_revised_trophic_position_seasonal_contrasts.csv\n",
    "Table_10_exact_environmental_biological_Mantel_results.csv\n",
    "Table_11_Procrustes_concordance.csv\n",
    "Table_12_design_and_inference_limitations.csv\n",
    "\nMain figures:\n",
    "Figure_1_environmental_PCA\n",
    "Figure_2_standardized_environmental_heatmap\n",
    "Figure_3_selected_environmental_variables\n",
    "Figure_4_resource_isotope_space\n",
    "Figure_5_revised_trophic_position\n",
    "Figure_6_environmental_biological_Mantel_heatmap\n",
    "Figure_7_Procrustes_environment_biology\n",
    "\nImportant:\n",
    "- Environmental–biological matrix tests are exact but exploratory.\n",
    "- Revised TP uses all 609 fish with no outcome-dependent exclusions.\n",
    "- Historical TP and original QC fields are ignored for inclusion.\n",
    "- Gastropod-referenced TP uses fixed-factor nonparametric bootstrap intervals.\n",
    "- All six target-gastropod baselines are observed and campaign-matched.\n",
    "- Bulk zooplankton is an exploratory reference, not a validated TP baseline.\n",
    "- Figure_8_gastropod_TP_environmental_concordance combines panels (a) and (b).\n",
    "============================================================\n",
    sep = ""
  )

  invisible(
    list(
      environmental_pca = environmental_pca,
      environmental_permanova = overall_permanova_table,
      revised_tp_group_summary = revised_tp_group_summary,
      environmental_biological_mantel = overall_environmental_linkage,
      procrustes_summary = procrustes_summary,
      results_directory = RESULTS_DIR
    )
  )
}


# =============================================================================
# 5. RUN ANALYSIS WITH CLEAR ERROR LOGGING
# =============================================================================

analysis_result <- tryCatch(
  run_environmental_integration(),

  error = function(e) {

    error_message <- paste0(
      "Environmental integration and revised trophic-position analysis stopped with an error.\n",
      "Time: ",
      Sys.time(),
      "\n\n",
      conditionMessage(
        e
      ),
      "\n"
    )

    message(
      "\n",
      error_message
    )

    fallback_error_file <- file.path(
      getwd(),
      "ENVIRONMENTAL_INTEGRATION_REVISED_TP_ERROR.txt"
    )

    writeLines(
      error_message,
      fallback_error_file
    )

    stop(
      conditionMessage(
        e
      ),
      call. = FALSE
    )
  }
)

# =============================================================================
# LOWER HAN RIVER FISH STABLE ISOTOPES
# LINKING FISH-COMMUNITY TURNOVER TO TROPHIC REORGANIZATION
# FINAL READY-TO-USE R SCRIPT — IMPROVED PUBLICATION VISUALS
#
# THIS IS NOT A BAYESIAN MODEL.
#
# The analysis links community dissimilarity from the fish-community analysis
# to posterior-median isotope-niche metrics from the preceding SIBER analysis.
#
# PRIMARY QUESTION
#   Does greater fish-community turnover correspond to:
#     1. greater isotope-centroid separation;
#     2. larger differences in standard ellipse area;
#     3. larger proportional differences in ellipse area; or
#     4. lower isotope-ellipse overlap?
#
# PRIMARY INFERENCE
#   Exploratory Spearman matrix correlations compare six group summaries.
#   Whole wet-dry reach pairs (6 arrangements) and within-reach campaign swaps
#   (8 arrangements) are separate conditional null sensitivities. The main
#   figures show effect estimates. Non-significance does not demonstrate
#   complementarity, independence, or equivalence.
#
# IMPORTANT LIMITATIONS
#   - Only six aggregated reach-campaign groups are available.
#   - The 15 pairwise observations are not independent; ordinary regression
#     P-values are therefore not used.
#   - The wet-dry comparison among reaches has n = 3 and is descriptive only.
#   - Isotope outcomes are posterior medians imported from SIBER; their full
#     posterior uncertainty is not propagated in this linking analysis.
#
# INPUT WORKBOOK
#   Han_River_Community_Turnover_Trophic_Reorganization_R_Dataset.xlsx
#
# HOW TO RUN
#   1. Restart RStudio.
#   2. Open this script and click Source once.
#   3. Select the Excel workbook if requested.
# =============================================================================


# =============================================================================
# 1. USER SETTINGS
# =============================================================================

PRIMARY_METHOD <- "spearman"
SENSITIVITY_METHOD <- "pearson"
SEED <- 20260731L

set.seed(SEED)
graphics.off()


# =============================================================================
# 2. PACKAGE CHECK
# =============================================================================

required_packages <- c(
  "vegan",
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
  library(vegan)
  library(readxl)
  library(dplyr)
  library(ggplot2)
})


# =============================================================================
# 3. MAIN ANALYSIS
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


run_linkage_analysis <- function() {

  # ---------------------------------------------------------------------------
  # 3.1 Locate the workbook
  # ---------------------------------------------------------------------------

  workbook_pattern <- paste0(
    "^Han_River_Community_Turnover_Trophic_",
    "Reorganization_R_Dataset.*[.]xlsx$"
  )

  user_home <- path.expand("~")
  windows_user <- Sys.getenv("USERPROFILE")

  search_folders <- unique(
    c(
      getwd(),
      { src <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
        if (!is.null(src) && nzchar(src)) dirname(normalizePath(src, mustWork = FALSE)) else
        if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
          p <- rstudioapi::getActiveDocumentContext()$path
          if (nzchar(p)) dirname(p) else getwd()
        } else {
          a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
          if (length(a)) dirname(sub("^--file=", "", a[1L])) else getwd()
        }
      },
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

  if (length(candidate_files) == 0L) {

    message(
      "\nSelect the community-turnover linkage workbook."
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

  # Prefer the revised workbook if earlier versions are in the same folder.
  preferred_r3 <- grepl("_Reviewer3[.]xlsx$", basename(candidate_files), ignore.case = TRUE)
  candidate_files <- c(candidate_files[preferred_r3], candidate_files[!preferred_r3])

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
  # 3.2 Output folders and status files
  # ---------------------------------------------------------------------------

  RESULTS_DIR <- file.path(
    dirname(INPUT_FILE),
    "Community_turnover_results_Reviewer3"
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
      "Community-turnover linkage analysis is running.",
      paste("Started:", Sys.time()),
      paste("Input workbook:", INPUT_FILE)
    ),
    RUNNING_FILE
  )


  # ---------------------------------------------------------------------------
  # 3.3 Read workbook
  # ---------------------------------------------------------------------------

  available_sheets <- readxl::excel_sheets(
    INPUT_FILE
  )

  required_sheets <- c(
    "Group_Metadata",
    "Pairwise_Linkage",
    "Seasonal_Linkage",
    "Expected_Mantel_Audit"
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

  group_metadata <- readxl::read_excel(
    INPUT_FILE,
    sheet = "Group_Metadata"
  )

  pairwise_data <- readxl::read_excel(
    INPUT_FILE,
    sheet = "Pairwise_Linkage"
  )

  seasonal_data <- readxl::read_excel(
    INPUT_FILE,
    sheet = "Seasonal_Linkage"
  )

  expected_audit <- readxl::read_excel(
    INPUT_FILE,
    sheet = "Expected_Mantel_Audit"
  )

  group_metadata <- as.data.frame(
    group_metadata,
    stringsAsFactors = FALSE
  )

  pairwise_data <- as.data.frame(
    pairwise_data,
    stringsAsFactors = FALSE
  )

  seasonal_data <- as.data.frame(
    seasonal_data,
    stringsAsFactors = FALSE
  )

  expected_audit <- as.data.frame(
    expected_audit,
    stringsAsFactors = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.4 Validate input columns
  # ---------------------------------------------------------------------------

  required_pairwise_columns <- c(
    "Pair_ID",
    "Sample_1",
    "Sample_2",
    "Comparison_type",
    "Manuscript_priority",
    "Bray_Curtis",
    "Jaccard",
    "Median_centroid_distance",
    "Absolute_delta_d13C",
    "Absolute_delta_d15N",
    "Absolute_SEA_difference",
    "Absolute_log_SEA_ratio",
    "Niche_nonoverlap_proportion"
  )

  required_seasonal_columns <- c(
    "Reach_code",
    "Reach_name",
    "Bray_Curtis_turnover",
    "Jaccard_turnover",
    "Median_centroid_distance",
    "Absolute_delta_d13C",
    "Absolute_delta_d15N",
    "Absolute_SEA_change",
    "Absolute_log_SEA_ratio",
    "Niche_nonoverlap_proportion"
  )

  required_group_columns <- c(
    "Sample_ID",
    "SIBER_code",
    "Reach_code",
    "Campaign",
    "Posterior_median_d13C",
    "Posterior_median_d15N",
    "Posterior_median_SEA"
  )

  missing_pairwise_columns <- setdiff(
    required_pairwise_columns,
    names(pairwise_data)
  )

  missing_seasonal_columns <- setdiff(
    required_seasonal_columns,
    names(seasonal_data)
  )

  missing_group_columns <- setdiff(
    required_group_columns,
    names(group_metadata)
  )

  if (length(missing_pairwise_columns) > 0L) {
    stop(
      paste0(
        "Pairwise_Linkage is missing column(s): ",
        paste(missing_pairwise_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (length(missing_seasonal_columns) > 0L) {
    stop(
      paste0(
        "Seasonal_Linkage is missing column(s): ",
        paste(missing_seasonal_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (length(missing_group_columns) > 0L) {
    stop(
      paste0(
        "Group_Metadata is missing column(s): ",
        paste(missing_group_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  numeric_pairwise_columns <- c(
    "Bray_Curtis",
    "Jaccard",
    "Median_centroid_distance",
    "Absolute_delta_d13C",
    "Absolute_delta_d15N",
    "Absolute_SEA_difference",
    "Absolute_log_SEA_ratio",
    "Niche_nonoverlap_proportion"
  )

  numeric_seasonal_columns <- c(
    "Bray_Curtis_turnover",
    "Jaccard_turnover",
    "Median_centroid_distance",
    "Absolute_delta_d13C",
    "Absolute_delta_d15N",
    "Absolute_SEA_change",
    "Absolute_log_SEA_ratio",
    "Niche_nonoverlap_proportion"
  )

  for (column_name in numeric_pairwise_columns) {

    pairwise_data[[column_name]] <- suppressWarnings(
      as.numeric(
        pairwise_data[[column_name]]
      )
    )
  }

  for (column_name in numeric_seasonal_columns) {

    seasonal_data[[column_name]] <- suppressWarnings(
      as.numeric(
        seasonal_data[[column_name]]
      )
    )
  }

  if (anyNA(
    pairwise_data[
      ,
      c(
        "Sample_1",
        "Sample_2",
        numeric_pairwise_columns
      )
    ]
  )) {
    stop(
      "Pairwise_Linkage contains missing or invalid values.",
      call. = FALSE
    )
  }

  if (anyNA(
    seasonal_data[
      ,
      c(
        "Reach_code",
        numeric_seasonal_columns
      )
    ]
  )) {
    stop(
      "Seasonal_Linkage contains missing or invalid values.",
      call. = FALSE
    )
  }

  group_order <- as.character(
    group_metadata$Sample_ID
  )

  if (
    length(
      group_order
    ) != 6L ||
    anyDuplicated(
      group_order
    )
  ) {
    stop(
      "Group_Metadata must contain six unique Sample_ID values.",
      call. = FALSE
    )
  }

  expected_pair_count <- choose(
    length(
      group_order
    ),
    2
  )

  if (
    nrow(
      pairwise_data
    ) != expected_pair_count ||
    anyDuplicated(
      pairwise_data$Pair_ID
    )
  ) {
    stop(
      "Pairwise_Linkage must contain all 15 unique pairwise comparisons.",
      call. = FALSE
    )
  }

  if (nrow(
    seasonal_data
  ) != 3L) {
    stop(
      "Seasonal_Linkage must contain exactly three reach-level wet-dry comparisons.",
      call. = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # 3.5 Build matrices
  # ---------------------------------------------------------------------------

  build_symmetric_matrix <- function(
    data,
    value_column,
    sample_order
  ) {

    n_groups <- length(
      sample_order
    )

    output_matrix <- matrix(
      0,
      nrow = n_groups,
      ncol = n_groups,
      dimnames = list(
        sample_order,
        sample_order
      )
    )

    for (row_index in seq_len(
      nrow(
        data
      )
    )) {

      sample_1 <- as.character(
        data$Sample_1[
          row_index
        ]
      )

      sample_2 <- as.character(
        data$Sample_2[
          row_index
        ]
      )

      if (
        !(sample_1 %in% sample_order) ||
        !(sample_2 %in% sample_order)
      ) {
        stop(
          paste0(
            "Unknown sample identifier in pairwise data: ",
            sample_1,
            " or ",
            sample_2
          ),
          call. = FALSE
        )
      }

      current_value <- as.numeric(
        data[[
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


  matrix_list <- list(
    Bray_Curtis = build_symmetric_matrix(
      pairwise_data,
      "Bray_Curtis",
      group_order
    ),

    Jaccard = build_symmetric_matrix(
      pairwise_data,
      "Jaccard",
      group_order
    ),

    Centroid_Distance = build_symmetric_matrix(
      pairwise_data,
      "Median_centroid_distance",
      group_order
    ),

    Absolute_SEA_Difference = build_symmetric_matrix(
      pairwise_data,
      "Absolute_SEA_difference",
      group_order
    ),

    Absolute_Log_SEA_Ratio = build_symmetric_matrix(
      pairwise_data,
      "Absolute_log_SEA_ratio",
      group_order
    ),

    Niche_Nonoverlap = build_symmetric_matrix(
      pairwise_data,
      "Niche_nonoverlap_proportion",
      group_order
    )
  )


  matrix_validation_rows <- lapply(
    names(
      matrix_list
    ),
    function(matrix_name) {

      current_matrix <- matrix_list[[
        matrix_name
      ]]

      data.frame(
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
    }
  )

  matrix_validation <- dplyr::bind_rows(
    matrix_validation_rows
  )

  if (
    any(
      matrix_validation$Maximum_symmetry_error > 1e-12
    ) ||
    any(
      matrix_validation$Maximum_absolute_diagonal > 1e-12
    )
  ) {
    stop(
      "At least one linkage matrix failed symmetry or diagonal validation.",
      call. = FALSE
    )
  }

  write.csv(
    matrix_validation,
    file.path(
      TABLE_DIR,
      "Table_S1_matrix_validation.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.6 Reviewer-3 paired null sensitivities
  # ---------------------------------------------------------------------------
  n_groups <- length(group_order)
  r3_designs <- r3_make_designs(group_order)
  exact_permutation_matrix <- r3_nonidentity(r3_designs[[R3_DISPLAY_SCHEME]])
  r3_write_design_files(TABLE_DIR, r3_designs)
  community_metric_names <- c("Bray_Curtis", "Jaccard")
  outcome_metric_names <- c("Centroid_Distance", "Absolute_SEA_Difference",
                            "Absolute_Log_SEA_Ratio", "Niche_Nonoverlap")
  paired_mantel_results <- r3_matrix_tests(
    matrix_list[community_metric_names], matrix_list[outcome_metric_names],
    c(PRIMARY_METHOD, SENSITIVITY_METHOD), r3_designs)
  primary_mantel_results <- subset(paired_mantel_results,
    Permutation_scheme == R3_DISPLAY_SCHEME & Correlation_method == PRIMARY_METHOD)
  sensitivity_mantel_results <- subset(paired_mantel_results,
    Permutation_scheme == R3_DISPLAY_SCHEME & Correlation_method == SENSITIVITY_METHOD)
  write.csv(primary_mantel_results, file.path(TABLE_DIR,
    "Table_1_exact_Spearman_Mantel_results.csv"), row.names = FALSE)
  write.csv(sensitivity_mantel_results, file.path(TABLE_DIR,
    "Table_S2_exact_Pearson_Mantel_sensitivity.csv"), row.names = FALSE)
  write.csv(paired_mantel_results, file.path(TABLE_DIR,
    "Reviewer3_all_paired_Mantel_results.csv"), row.names = FALSE)
  point_rankings <- r3_seasonal_rankings(seasonal_data)
  write.csv(point_rankings, file.path(TABLE_DIR,
    "Reviewer3_seasonal_point_rankings.csv"), row.names = FALSE)

  # The earlier reviewer-2 catch sensitivity must use the same restricted nulls.
  catch_paired_results <- data.frame()
  if ("Catch_Pairwise_R3" %in% available_sheets) {
    catch_pairs <- as.data.frame(readxl::read_excel(INPUT_FILE, sheet = "Catch_Pairwise_R3"))
    if (nrow(catch_pairs) != 15L || anyDuplicated(catch_pairs$Pair_ID))
      stop("Catch_Pairwise_R3 requires 15 unique pairwise records.")
    catch_matrices <- setNames(lapply(community_metric_names, function(v) {
      build_symmetric_matrix(catch_pairs, v, group_order)
    }), community_metric_names)
    catch_paired_results <- r3_matrix_tests(catch_matrices,
      matrix_list[outcome_metric_names], PRIMARY_METHOD, r3_designs)
    write.csv(catch_paired_results, file.path(TABLE_DIR,
      "Reviewer3_Table_S10d_catch_Mantel.csv"), row.names = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 3.8 Verify exact results against the workbook audit
  # ---------------------------------------------------------------------------

  expected_audit$Spearman_Mantel_rho <- as.numeric(
    expected_audit$Spearman_Mantel_rho
  )

  expected_audit$Exact_one_sided_P <- as.numeric(
    expected_audit$Exact_one_sided_P
  )

  expected_audit$Exact_two_sided_P <- as.numeric(
    expected_audit$Exact_two_sided_P
  )

  audit_comparison <- merge(
    primary_mantel_results,
    expected_audit[
      ,
      c(
        "Community_metric",
        "Trophic_reorganization_metric",
        "Spearman_Mantel_rho",
        "Exact_one_sided_P",
        "Exact_two_sided_P"
      ),
      drop = FALSE
    ],
    by = c(
      "Community_metric",
      "Trophic_reorganization_metric"
    ),
    suffixes = c(
      "_R",
      "_Workbook"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  audit_comparison$Absolute_rho_error <- abs(
    audit_comparison$Matrix_correlation -
    audit_comparison$Spearman_Mantel_rho
  )

  audit_comparison$Absolute_one_sided_P_error <- abs(
    audit_comparison$Exact_one_sided_P_R -
    audit_comparison$Exact_one_sided_P_Workbook
  )

  audit_comparison$Absolute_two_sided_P_error <- abs(
    audit_comparison$Exact_two_sided_P_R -
    audit_comparison$Exact_two_sided_P_Workbook
  )

  maximum_audit_error <- max(
    c(
      audit_comparison$Absolute_rho_error,
      audit_comparison$Absolute_one_sided_P_error,
      audit_comparison$Absolute_two_sided_P_error
    ),
    na.rm = TRUE
  )

  if (
    !is.finite(
      maximum_audit_error
    ) ||
    maximum_audit_error > 1e-10
  ) {
    stop(
      paste0(
        "Mantel-result audit failed. Maximum absolute error = ",
        maximum_audit_error
      ),
      call. = FALSE
    )
  }

  write.csv(
    audit_comparison,
    file.path(
      TABLE_DIR,
      "Table_S3_exact_result_audit.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.9 Pairwise linked metrics
  # ---------------------------------------------------------------------------

  write.csv(
    pairwise_data,
    file.path(
      TABLE_DIR,
      "Table_2_linked_pairwise_metrics.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.10 Descriptive seasonal associations among three reaches
  # ---------------------------------------------------------------------------

  seasonal_predictors <- c(
    "Bray_Curtis_turnover",
    "Jaccard_turnover"
  )

  seasonal_outcomes <- c(
    "Median_centroid_distance",
    "Absolute_delta_d13C",
    "Absolute_delta_d15N",
    "Absolute_SEA_change",
    "Absolute_log_SEA_ratio",
    "Niche_nonoverlap_proportion"
  )

  seasonal_association_rows <- list()
  seasonal_association_index <- 1L

  for (predictor_name in seasonal_predictors) {

    for (outcome_name in seasonal_outcomes) {

      current_rho <- stats::cor(
        seasonal_data[[
          predictor_name
        ]],
        seasonal_data[[
          outcome_name
        ]],
        method = "spearman",
        use = "complete.obs"
      )

      current_pearson <- stats::cor(
        seasonal_data[[
          predictor_name
        ]],
        seasonal_data[[
          outcome_name
        ]],
        method = "pearson",
        use = "complete.obs"
      )

      seasonal_association_rows[[
        seasonal_association_index
      ]] <- data.frame(
        Community_metric = predictor_name,
        Trophic_reorganization_metric = outcome_name,
        Spearman_rho = current_rho,
        Pearson_r = current_pearson,
        Number_of_reaches = nrow(
          seasonal_data
        ),
        Inference_status = "Descriptive only; no inferential P value",
        stringsAsFactors = FALSE
      )

      seasonal_association_index <- seasonal_association_index + 1L
    }
  }

  seasonal_associations <- dplyr::bind_rows(
    seasonal_association_rows
  )

  write.csv(
    seasonal_data,
    file.path(
      TABLE_DIR,
      "Table_3_within_reach_seasonal_linkage.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    seasonal_associations,
    file.path(
      TABLE_DIR,
      "Table_4_descriptive_seasonal_associations.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.11 Design and inference limitations
  # ---------------------------------------------------------------------------

  design_limitations <- data.frame(
    Item = c(
      "Analysis type",
      "Reach-campaign groups",
      "Unique pairwise distances",
      "Exact permutations",
      "Primary test",
      "Pairwise regression inference",
      "Seasonal reach-level inference",
      "SIBER posterior uncertainty",
      "Recommended interpretation"
    ),

    Value = c(
      "Non-Bayesian matrix-association analysis",
      "6",
      "15",
      "Whole reach pairs: 6 including observed; within-reach swaps: 8 including observed",
      "Exact Spearman Mantel test",
      "Not used because pairwise distances are non-independent",
      "Descriptive only because n = 3",
      "Not fully propagated; posterior medians are used",
      "Exploratory effect estimates; no demonstrated complementarity or independence"
    ),

    stringsAsFactors = FALSE
  )

  write.csv(
    design_limitations,
    file.path(
      TABLE_DIR,
      "Table_5_design_and_inference_limitations.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.12 Publication figures — improved colourblind-friendly design
  # ---------------------------------------------------------------------------

  # Metric colours
  COL_BRAY <- "#1F78B4"
  COL_JACCARD <- "#E66101"

  # Comparison-type colours
  comparison_colours <- c(
    "Within-reach wet–dry" = "#0072B2",
    "Among-reach, same campaign" = "#E69F00",
    "Cross-reach, cross-campaign" = "#CC79A7"
  )

  comparison_shapes <- c(
    "Within-reach wet–dry" = 21,
    "Among-reach, same campaign" = 22,
    "Cross-reach, cross-campaign" = 24
  )

  # Reach colours
  reach_colours <- c(
    "Hanchuan" = "#0072B2",
    "Qianjiang" = "#D55E00",
    "Shayang" = "#009E73"
  )

  reach_shapes <- c(
    "Hanchuan" = 21,
    "Qianjiang" = 22,
    "Shayang" = 24
  )

  # Diverging heat-map colours
  COL_HEAT_LOW <- "#2166AC"
  COL_HEAT_MID <- "#F7F7F7"
  COL_HEAT_HIGH <- "#B2182B"

  comparison_label_function <- function(
    comparison_type
  ) {

    dplyr::recode(
      comparison_type,
      "Within-reach wet-dry" = "Within-reach wet–dry",
      "Among-reach same campaign" = "Among-reach, same campaign",
      "Cross-reach cross-campaign" = "Cross-reach, cross-campaign"
    )
  }


  theme_publication <- function(
    base_size = 10
  ) {

    ggplot2::theme_bw(
      base_size = base_size,
      base_family = "sans"
    ) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(
          colour = "#E5E5E5",
          linewidth = 0.35
        ),

        panel.grid.minor = ggplot2::element_blank(),

        panel.border = ggplot2::element_rect(
          colour = "black",
          fill = NA,
          linewidth = 0.6
        ),

        strip.background = ggplot2::element_rect(
          fill = "#F0F0F0",
          colour = "black",
          linewidth = 0.5
        ),

        strip.text = ggplot2::element_text(
          face = "bold",
          colour = "#222222"
        ),

        axis.title = ggplot2::element_text(
          face = "bold",
          colour = "#222222"
        ),

        axis.text = ggplot2::element_text(
          colour = "#222222"
        ),

        legend.position = "top",

        legend.title = ggplot2::element_text(
          face = "bold"
        ),

        legend.key = ggplot2::element_rect(
          fill = "white",
          colour = NA
        ),

        plot.title = ggplot2::element_text(
          face = "bold",
          colour = "#222222"
        ),

        plot.caption = ggplot2::element_text(
          colour = "#444444",
          size = base_size - 1
        ),

        plot.margin = ggplot2::margin(
          9,
          12,
          9,
          9
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


  # ---------------------------------------------------------------------------
  # Figure 1. Community turnover versus isotope-centroid separation
  # ---------------------------------------------------------------------------

  centroid_scatter_data <- dplyr::bind_rows(
    data.frame(
      Pair_ID = pairwise_data$Pair_ID,
      Community_metric = "Bray–Curtis",
      Community_dissimilarity = pairwise_data$Bray_Curtis,
      Outcome = pairwise_data$Median_centroid_distance,
      Comparison_type = pairwise_data$Comparison_type,
      Manuscript_priority = pairwise_data$Manuscript_priority,
      stringsAsFactors = FALSE
    ),

    data.frame(
      Pair_ID = pairwise_data$Pair_ID,
      Community_metric = "Jaccard",
      Community_dissimilarity = pairwise_data$Jaccard,
      Outcome = pairwise_data$Median_centroid_distance,
      Comparison_type = pairwise_data$Comparison_type,
      Manuscript_priority = pairwise_data$Manuscript_priority,
      stringsAsFactors = FALSE
    )
  ) |>
    dplyr::mutate(
      Community_metric = factor(
        Community_metric,
        levels = c(
          "Bray–Curtis",
          "Jaccard"
        )
      ),

      Comparison_label = factor(
        comparison_label_function(
          Comparison_type
        ),
        levels = names(
          comparison_colours
        )
      )
    )

  figure_1 <- ggplot2::ggplot(
    centroid_scatter_data,
    ggplot2::aes(
      x = Community_dissimilarity,
      y = Outcome
    )
  ) +
    ggplot2::geom_smooth(
      ggplot2::aes(
        colour = Community_metric,
        group = Community_metric
      ),
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.9,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = Comparison_label,
        shape = Comparison_label
      ),
      size = 3.3,
      stroke = 0.65,
      colour = "black",
      alpha = 0.95
    ) +
    ggplot2::facet_wrap(
      ~ Community_metric,
      scales = "free_x"
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "Bray–Curtis" = COL_BRAY,
        "Jaccard" = COL_JACCARD
      )
    ) +
    ggplot2::scale_fill_manual(
      values = comparison_colours,
      name = "Comparison"
    ) +
    ggplot2::scale_shape_manual(
      values = comparison_shapes,
      name = "Comparison"
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
          0.05,
          0.10
        )
      )
    ) +
    ggplot2::labs(
      x = "Fish-community dissimilarity",
      y = "Posterior median isotope-centroid distance (‰)",
      caption = paste(
        "Lines are descriptive visual summaries;",
        "restricted permutations are conditional exploratory sensitivities."
      )
    ) +
    theme_publication(
      base_size = 10
    )

  save_figure(
    figure_1,
    "Figure_1_community_turnover_vs_centroid_distance",
    width = 7.8,
    height = 5.0
  )


  # ---------------------------------------------------------------------------
  # Figure 2. Community turnover versus isotope-niche non-overlap
  # ---------------------------------------------------------------------------

  nonoverlap_scatter_data <- dplyr::bind_rows(
    data.frame(
      Pair_ID = pairwise_data$Pair_ID,
      Community_metric = "Bray–Curtis",
      Community_dissimilarity = pairwise_data$Bray_Curtis,
      Outcome = pairwise_data$Niche_nonoverlap_proportion,
      Comparison_type = pairwise_data$Comparison_type,
      stringsAsFactors = FALSE
    ),

    data.frame(
      Pair_ID = pairwise_data$Pair_ID,
      Community_metric = "Jaccard",
      Community_dissimilarity = pairwise_data$Jaccard,
      Outcome = pairwise_data$Niche_nonoverlap_proportion,
      Comparison_type = pairwise_data$Comparison_type,
      stringsAsFactors = FALSE
    )
  ) |>
    dplyr::mutate(
      Community_metric = factor(
        Community_metric,
        levels = c(
          "Bray–Curtis",
          "Jaccard"
        )
      ),

      Comparison_label = factor(
        comparison_label_function(
          Comparison_type
        ),
        levels = names(
          comparison_colours
        )
      )
    )

  figure_2 <- ggplot2::ggplot(
    nonoverlap_scatter_data,
    ggplot2::aes(
      x = Community_dissimilarity,
      y = Outcome
    )
  ) +
    ggplot2::geom_smooth(
      ggplot2::aes(
        colour = Community_metric,
        group = Community_metric
      ),
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.9,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = Comparison_label,
        shape = Comparison_label
      ),
      size = 3.3,
      stroke = 0.65,
      colour = "black",
      alpha = 0.95
    ) +
    ggplot2::facet_wrap(
      ~ Community_metric,
      scales = "free_x"
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "Bray–Curtis" = COL_BRAY,
        "Jaccard" = COL_JACCARD
      )
    ) +
    ggplot2::scale_fill_manual(
      values = comparison_colours,
      name = "Comparison"
    ) +
    ggplot2::scale_shape_manual(
      values = comparison_shapes,
      name = "Comparison"
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
      limits = c(
        0,
        1
      ),
      breaks = seq(
        0,
        1,
        by = 0.2
      ),
      expand = ggplot2::expansion(
        mult = c(
          0.02,
          0.04
        )
      )
    ) +
    ggplot2::labs(
      x = "Fish-community dissimilarity",
      y = "Isotope-niche non-overlap proportion",
      caption = paste(
        "Niche non-overlap = 1 − posterior median Jaccard overlap",
        "of the 95% prediction ellipses."
      )
    ) +
    theme_publication(
      base_size = 10
    )

  save_figure(
    figure_2,
    "Figure_2_community_turnover_vs_niche_nonoverlap",
    width = 7.8,
    height = 5.0
  )


  # ---------------------------------------------------------------------------
  # Figure 3. Exact Mantel association heat map
  # ---------------------------------------------------------------------------

  outcome_labels <- c(
    Centroid_Distance = "Centroid distance",
    Absolute_SEA_Difference = "|SEA difference|",
    Absolute_Log_SEA_Ratio = "|log SEA ratio|",
    Niche_Nonoverlap = "Niche non-overlap"
  )

  community_labels <- c(
    Bray_Curtis = "Bray–Curtis",
    Jaccard = "Jaccard"
  )

  mantel_plot_data <- primary_mantel_results |>
    dplyr::mutate(
      Community_label = factor(
        community_labels[
          Community_metric
        ],
        levels = c(
          "Bray–Curtis",
          "Jaccard"
        )
      ),

      Outcome_label = factor(
        outcome_labels[
          Trophic_reorganization_metric
        ],
        levels = rev(
          c(
            "Centroid distance",
            "|SEA difference|",
            "|log SEA ratio|",
            "Niche non-overlap"
          )
        )
      ),

      Annotation = if (R3_SHOW_P_IN_FIGURES) sprintf(
        "ρ = %.2f\nP(pair) = %.3f", Matrix_correlation, Exact_one_sided_P
      ) else sprintf("ρ = %.2f", Matrix_correlation),

      Text_colour = ifelse(
        abs(
          Matrix_correlation
        ) >= 0.60,
        "white",
        "black"
      )
    )

  figure_3 <- ggplot2::ggplot(
    mantel_plot_data,
    ggplot2::aes(
      x = Community_label,
      y = Outcome_label,
      fill = Matrix_correlation
    )
  ) +
    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.8
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = Annotation,
        colour = Text_colour
      ),
      size = 3.4,
      fontface = "bold",
      lineheight = 0.95,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_gradient2(
      low = COL_HEAT_LOW,
      mid = COL_HEAT_MID,
      high = COL_HEAT_HIGH,
      midpoint = 0,
      limits = c(
        -1,
        1
      ),
      name = "Spearman\nMantel ρ"
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      caption = paste(
        "Observed correlations from six summaries in three paired reaches.",
        "Conditional restricted-null P values are reported in the sensitivity tables."
      )
    ) +
    ggplot2::coord_fixed() +
    theme_publication(
      base_size = 10
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank()
    )

  save_figure(
    figure_3,
    "Figure_3_exact_Mantel_association_heatmap",
    width = 6.9,
    height = 5.2
  )


  # ---------------------------------------------------------------------------
  # Figure 4. Seasonal community turnover versus centroid movement
  # ---------------------------------------------------------------------------

  seasonal_centroid_data <- dplyr::bind_rows(
    data.frame(
      Reach_name = seasonal_data$Reach_name,
      Community_metric = "Bray–Curtis",
      Community_turnover = seasonal_data$Bray_Curtis_turnover,
      Outcome = seasonal_data$Median_centroid_distance,
      stringsAsFactors = FALSE
    ),

    data.frame(
      Reach_name = seasonal_data$Reach_name,
      Community_metric = "Jaccard",
      Community_turnover = seasonal_data$Jaccard_turnover,
      Outcome = seasonal_data$Median_centroid_distance,
      stringsAsFactors = FALSE
    )
  ) |>
    dplyr::mutate(
      Community_metric = factor(
        Community_metric,
        levels = c(
          "Bray–Curtis",
          "Jaccard"
        )
      ),

      Reach_name = factor(
        Reach_name,
        levels = names(
          reach_colours
        )
      )
    )

  figure_4 <- ggplot2::ggplot(
    seasonal_centroid_data,
    ggplot2::aes(
      x = Community_turnover,
      y = Outcome
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = Reach_name,
        shape = Reach_name
      ),
      size = 4.0,
      stroke = 0.75,
      colour = "black"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = Reach_name,
        colour = Reach_name
      ),
      nudge_y = 0.045,
      size = 3.2,
      fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(
      ~ Community_metric,
      scales = "free_x"
    ) +
    ggplot2::scale_fill_manual(
      values = reach_colours,
      name = "Reach"
    ) +
    ggplot2::scale_colour_manual(
      values = reach_colours,
      guide = "none"
    ) +
    ggplot2::scale_shape_manual(
      values = reach_shapes,
      name = "Reach"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0.10,
          0.10
        )
      )
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0.08,
          0.18
        )
      )
    ) +
    ggplot2::labs(
      x = "Wet–dry fish-community dissimilarity",
      y = "Wet–dry isotope-centroid movement (‰)",
      caption = "Descriptive comparison only because three reaches were available."
    ) +
    theme_publication(
      base_size = 10
    )

  save_figure(
    figure_4,
    "Figure_4_seasonal_turnover_vs_centroid_movement",
    width = 7.6,
    height = 4.9
  )


  # ---------------------------------------------------------------------------
  # Figure 5. Seasonal community turnover versus niche-area change
  # ---------------------------------------------------------------------------

  seasonal_sea_data <- dplyr::bind_rows(
    data.frame(
      Reach_name = seasonal_data$Reach_name,
      Community_metric = "Bray–Curtis",
      Community_turnover = seasonal_data$Bray_Curtis_turnover,
      Outcome = seasonal_data$Absolute_SEA_change,
      stringsAsFactors = FALSE
    ),

    data.frame(
      Reach_name = seasonal_data$Reach_name,
      Community_metric = "Jaccard",
      Community_turnover = seasonal_data$Jaccard_turnover,
      Outcome = seasonal_data$Absolute_SEA_change,
      stringsAsFactors = FALSE
    )
  ) |>
    dplyr::mutate(
      Community_metric = factor(
        Community_metric,
        levels = c(
          "Bray–Curtis",
          "Jaccard"
        )
      ),

      Reach_name = factor(
        Reach_name,
        levels = names(
          reach_colours
        )
      )
    )

  figure_5 <- ggplot2::ggplot(
    seasonal_sea_data,
    ggplot2::aes(
      x = Community_turnover,
      y = Outcome
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = Reach_name,
        shape = Reach_name
      ),
      size = 4.0,
      stroke = 0.75,
      colour = "black"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = Reach_name,
        colour = Reach_name
      ),
      nudge_y = 0.45,
      size = 3.2,
      fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(
      ~ Community_metric,
      scales = "free_x"
    ) +
    ggplot2::scale_fill_manual(
      values = reach_colours,
      name = "Reach"
    ) +
    ggplot2::scale_colour_manual(
      values = reach_colours,
      guide = "none"
    ) +
    ggplot2::scale_shape_manual(
      values = reach_shapes,
      name = "Reach"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0.10,
          0.10
        )
      )
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0.06,
          0.18
        )
      )
    ) +
    ggplot2::labs(
      x = "Wet–dry fish-community dissimilarity",
      y = expression(
        "Absolute wet–dry SEA change ("
        * "\u2030"^2 *
        ")"
      ),
      caption = "Descriptive comparison only because three reaches were available."
    ) +
    theme_publication(
      base_size = 10
    )

  save_figure(
    figure_5,
    "Figure_5_seasonal_turnover_vs_niche_area_change",
    width = 7.6,
    height = 4.9
  )


  # ---------------------------------------------------------------------------
  # Figure S1. Community turnover versus proportional SEA difference
  # ---------------------------------------------------------------------------

  log_ratio_scatter_data <- dplyr::bind_rows(
    data.frame(
      Pair_ID = pairwise_data$Pair_ID,
      Community_metric = "Bray–Curtis",
      Community_dissimilarity = pairwise_data$Bray_Curtis,
      Outcome = pairwise_data$Absolute_log_SEA_ratio,
      Comparison_type = pairwise_data$Comparison_type,
      stringsAsFactors = FALSE
    ),

    data.frame(
      Pair_ID = pairwise_data$Pair_ID,
      Community_metric = "Jaccard",
      Community_dissimilarity = pairwise_data$Jaccard,
      Outcome = pairwise_data$Absolute_log_SEA_ratio,
      Comparison_type = pairwise_data$Comparison_type,
      stringsAsFactors = FALSE
    )
  ) |>
    dplyr::mutate(
      Community_metric = factor(
        Community_metric,
        levels = c(
          "Bray–Curtis",
          "Jaccard"
        )
      ),

      Comparison_label = factor(
        comparison_label_function(
          Comparison_type
        ),
        levels = names(
          comparison_colours
        )
      )
    )

  figure_s1 <- ggplot2::ggplot(
    log_ratio_scatter_data,
    ggplot2::aes(
      x = Community_dissimilarity,
      y = Outcome
    )
  ) +
    ggplot2::geom_smooth(
      ggplot2::aes(
        colour = Community_metric,
        group = Community_metric
      ),
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.9,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = Comparison_label,
        shape = Comparison_label
      ),
      size = 3.3,
      stroke = 0.65,
      colour = "black",
      alpha = 0.95
    ) +
    ggplot2::facet_wrap(
      ~ Community_metric,
      scales = "free_x"
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "Bray–Curtis" = COL_BRAY,
        "Jaccard" = COL_JACCARD
      )
    ) +
    ggplot2::scale_fill_manual(
      values = comparison_colours,
      name = "Comparison"
    ) +
    ggplot2::scale_shape_manual(
      values = comparison_shapes,
      name = "Comparison"
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
          0.05,
          0.10
        )
      )
    ) +
    ggplot2::labs(
      x = "Fish-community dissimilarity",
      y = "|log posterior median SEA ratio|",
      caption = paste(
        "Lines are descriptive visual summaries;",
        "exact matrix tests are reported in Table 1."
      )
    ) +
    theme_publication(
      base_size = 10
    )

  save_figure(
    figure_s1,
    "Figure_S1_community_turnover_vs_log_SEA_ratio",
    width = 7.8,
    height = 5.0
  )


  # ---------------------------------------------------------------------------
  # 3.13 Optional combined Excel result workbook
  # ---------------------------------------------------------------------------

  if (
    requireNamespace(
      "writexl",
      quietly = TRUE
    )
  ) {

    writexl::write_xlsx(
      list(
        Exact_Spearman_Mantel = primary_mantel_results,
        Exact_Pearson_Sensitivity = sensitivity_mantel_results,
        Paired_Mantel_All = paired_mantel_results,
        Catch_Mantel_R3 = catch_paired_results,
        Permutation_Design = r3_permutation_table(r3_designs),
        Null_Hypotheses = r3_null_hypotheses(),
        Seasonal_Point_Rankings = point_rankings,
        Linked_Pairwise_Metrics = pairwise_data,
        Seasonal_Linkage = seasonal_data,
        Seasonal_Associations = seasonal_associations,
        Matrix_Validation = matrix_validation,
        Exact_Result_Audit = audit_comparison,
        Design_Limitations = design_limitations
      ),
      path = file.path(
        RESULTS_DIR,
        "Community_turnover_trophic_reorganization_results.xlsx"
      )
    )
  }


  # ---------------------------------------------------------------------------
  # 3.14 Save session and completion log
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
      "Community-turnover linkage analysis completed successfully.",
      paste("Completed:", Sys.time()),
      paste("Input workbook:", INPUT_FILE),
      paste("Groups:", n_groups),
      paste("Pairwise distances:", nrow(pairwise_data)),
      paste("Exact non-identity permutations:", nrow(exact_permutation_matrix)),
      paste("Maximum audit error:", maximum_audit_error),
      paste("Results folder:", RESULTS_DIR)
    ),
    COMPLETED_FILE
  )

  cat(
    "\n============================================================\n",
    "COMMUNITY TURNOVER–TROPHIC REORGANIZATION ANALYSIS COMPLETED — IMPROVED VISUALS\n",
    "============================================================\n",
    "This was an exploratory matrix comparison with conditional paired-null sensitivities.\n",
    "\nResults folder:\n",
    RESULTS_DIR,
    "\n\nMain tables:\n",
    "Table_1_exact_Spearman_Mantel_results.csv\n",
    "Table_2_linked_pairwise_metrics.csv\n",
    "Table_3_within_reach_seasonal_linkage.csv\n",
    "Table_4_descriptive_seasonal_associations.csv\n",
    "Table_5_design_and_inference_limitations.csv\n",
    "\nMain figures:\n",
    "Figure_1_community_turnover_vs_centroid_distance\n",
    "Figure_2_community_turnover_vs_niche_nonoverlap\n",
    "Figure_3_exact_Mantel_association_heatmap\n",
    "Figure_4_seasonal_turnover_vs_centroid_movement\n",
    "Figure_5_seasonal_turnover_vs_niche_area_change\n",
    "Figure_S1_community_turnover_vs_log_SEA_ratio\n",
    "\nImportant: seasonal reach-level associations are descriptive only (n = 3).\n",
    "============================================================\n",
    sep = ""
  )

  invisible(
    list(
      exact_spearman_mantel = primary_mantel_results,
      exact_pearson_sensitivity = sensitivity_mantel_results,
      seasonal_associations = seasonal_associations,
      results_directory = RESULTS_DIR
    )
  )
}


# =============================================================================
# 4. RUN AND SAVE A CLEAR ERROR LOG
# =============================================================================

analysis_result <- tryCatch(
  run_linkage_analysis(),

  error = function(e) {

    error_message <- paste0(
      "Community-turnover linkage analysis stopped with an error.\n",
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
      "COMMUNITY_TURNOVER_TROPHIC_LINKAGE_ERROR.txt"
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

# RUN IN RSTUDIO: open this R file and click Source.
# Input file(s): Reviewer4_Community_Counts (1).csv; Reviewer4_CPUE_Event_Records (1).csv
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

# ============================================================================
# Reviewer comment 4 - complete RStudio analysis for Tables S5b, S5c and S9
# ============================================================================
# No Python and no additional R packages are required.
# R >= 3.6 is required; a current R 4.x installation is recommended.
#
# Put these THREE files in the same folder:
#   Reviewer4_Analysis_RStudio.R
#   Reviewer4_Community_Counts (1).csv
#   Reviewer4_CPUE_Event_Records (1).csv
# Open this script in RStudio and click Source to run the whole script.
# Alternatively: source("Reviewer4_Analysis_RStudio.R", encoding = "UTF-8")
#
# Scope:
#   - S3e: exports the header correction map; editing headings needs no analysis.
#   - S5b: within-reach wet-to-dry Bray-Curtis and Jaccard dissimilarities.
#   - S5c: separate ONE-FACTOR community PERMANOVA models, not an interaction.
#   - S9a: recorded event groups, catch, inferred net-hours and pooled CPUE CIs.
#   - S9b: wet-to-dry percentage changes with conditional bootstrap CIs.
#   - Ready-to-copy manuscript wording, audit files, settings and sessionInfo.
# This is not the entire manuscript's brms/SIBER/environmental workflow.
#
# IMPORTANT DATA SCOPE
# The CPUE input contains reconstructed catch-record groups, not an independently
# verified register of physical deployments. Do not call their counts verified
# deployment counts without checking the original field log and zero catches.
# The bootstrap is conditional on these groups and their reported effort.
# R's random-number generator differs from NumPy's; use the intervals generated
# by THIS run throughout the tables and manuscript. Do not mix R/Python CIs.
#
# Method references:
# https://vegandevs.github.io/vegan/reference/adonis.html
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/Random.html
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/sample.html
# ============================================================================

# ---------------------------- EDITABLE SETTINGS -----------------------------
INPUT_DIR <- hr_root()
OUTPUT_DIR <- hr_results("reviewer4")
COMMUNITY_FILENAME <- "Reviewer4_Community_Counts (1).csv"
EVENT_FILENAME <- "Reviewer4_CPUE_Event_Records (1).csv"
BOOTSTRAP_REPLICATES <- 5000L
BOOTSTRAP_SEED <- 20230826L
CONFIDENCE_LEVEL <- 0.95
CHECK_WITH_VEGAN_IF_INSTALLED <- TRUE
SAVE_ANALYSIS_OBJECTS <- TRUE
RUN_CONDITIONAL_CPUE <- TRUE

# These are the isotope/environmental windows stated in the current manuscript.
# They are used for AUDITING source dates, not for filtering out catch records.
STATED_CAMPAIGN_WINDOWS <- data.frame(
  Campaign = c("Wet_2023", "Dry_2023"),
  Start_date = as.Date(c("2023-06-20", "2023-11-14")),
  End_date = as.Date(c("2023-06-27", "2023-11-23")),
  stringsAsFactors = FALSE
)

# ------------------------------ MAIN FUNCTION -------------------------------
run_reviewer4_analysis <- function(
    input_dir = "", output_dir = "",
    community_filename = "Reviewer4_Community_Counts (1).csv",
    event_filename = "Reviewer4_CPUE_Event_Records (1).csv",
    n_boot = 5000L, seed = 20230826L, ci_level = 0.95,
    check_vegan = TRUE, save_objects = TRUE,
    run_conditional_cpue = TRUE,
    stated_windows = STATED_CAMPAIGN_WINDOWS) {

  if (getRversion() < "3.6.0") stop("Please use R 3.6 or newer.")
  if (length(n_boot) != 1L || !is.finite(n_boot) || n_boot < 1000 || n_boot != as.integer(n_boot)) {
    stop("n_boot must be an integer of at least 1000; the study setting is 5000.")
  }
  if (length(seed) != 1L || !is.finite(seed) || seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("seed must be one nonnegative integer within R's integer range.")
  }
  if (length(ci_level) != 1L || !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be between zero and one.")
  }
  n_boot <- as.integer(n_boot)
  seed <- as.integer(seed)
  old_options <- options(digits = 17)
  on.exit(options(old_options), add = TRUE)

  # Local utilities ----------------------------------------------------------
  read_input <- function(path) {
    x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                         fileEncoding = "UTF-8-BOM", na.strings = c("", "NA"))
    names(x) <- sub("^\ufeff", "", names(x))
    hr_normalize_site_names(x)
  }
  write_output <- function(x, filename) {
    utils::write.csv(x, file.path(output_dir, filename), row.names = FALSE,
                     na = "NA", fileEncoding = "UTF-8")
  }
  required_columns <- function(x, columns, label) {
    missing <- setdiff(columns, names(x))
    if (length(missing)) stop(label, " is missing: ", paste(missing, collapse = ", "))
  }
  numeric_column <- function(x, label) {
    result <- suppressWarnings(as.numeric(x))
    if (anyNA(result) || any(!is.finite(result))) stop(label, " has missing or nonnumeric values.")
    result
  }
  format_p <- function(x) ifelse(x < 0.001, "< 0.001", paste0("= ", sprintf("%.3f", x)))
  format_number <- function(x, digits = 3L) {
    ifelse(is.na(x), "Not estimable", formatC(x, digits = digits, format = "f"))
  }
  interval_text <- function(low, high, digits = 3L) {
    ifelse(is.na(low) | is.na(high), "Not estimable",
           paste0("[", format_number(low, digits), ", ", format_number(high, digits), "]"))
  }
  site_name <- function(code) {
    map <- c(HC = "S1", QJ = "S2", SY = "S3")
    result <- unname(map[code])
    result[is.na(result)] <- code[is.na(result)]
    result
  }

  # Locate both inputs without changing the user's working directory ----------
  files_present <- function(directory) {
    length(directory) == 1L && !is.na(directory) && nzchar(directory) &&
      all(file.exists(file.path(directory, c(community_filename, event_filename))))
  }
  if (nzchar(input_dir)) {
    if (!files_present(input_dir)) stop("Both input CSV files must be present in INPUT_DIR: ", input_dir)
  } else {
    source_paths <- vapply(sys.frames(), function(frame) {
      value <- frame$ofile
      if (is.character(value) && length(value) == 1L) value else ""
    }, character(1L))
    source_paths <- source_paths[nzchar(source_paths)]
    command_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    command_file <- sub("^--file=", "", command_file)
    candidates <- unique(c(dirname(rev(source_paths)), dirname(command_file), getwd()))
    if (requireNamespace("rstudioapi", quietly = TRUE)) {
      editor_path <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
      if (nzchar(editor_path)) candidates <- unique(c(candidates, dirname(editor_path)))
    }
    found <- candidates[vapply(candidates, files_present, logical(1L))]
    if (!length(found) && interactive()) {
      message("Select Reviewer4_Community_Counts (1).csv; both inputs must be in its folder.")
      chosen <- file.choose()
      if (basename(chosen) == community_filename && files_present(dirname(chosen))) found <- dirname(chosen)
    }
    if (!length(found)) stop("Cannot locate both CSV files. Put them beside the script or set INPUT_DIR.")
    input_dir <- found[1L]
  }
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  if (!nzchar(output_dir)) output_dir <- file.path(input_dir, "Reviewer4_R_Results")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Could not create the output folder.")
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  community_path <- file.path(input_dir, community_filename)
  event_path <- file.path(input_dir, event_filename)
  message("Inputs: ", input_dir)
  message("Outputs: ", output_dir)

  # 1. Validate and prepare the six community samples -------------------------
  community <- read_input(community_path)
  required_columns(community, "Sample_ID", "Community input")
  if (nrow(community) != 6L || anyNA(community$Sample_ID) || anyDuplicated(community$Sample_ID)) {
    stop("This design requires six unique groups: three reaches, each wet and dry.")
  }
  split_ids <- strsplit(as.character(community$Sample_ID), "__", fixed = TRUE)
  if (any(lengths(split_ids) != 2L)) stop("Sample_ID must look like HC__Wet_2023.")
  metadata <- data.frame(
    Sample_ID = as.character(community$Sample_ID),
    Reach = vapply(split_ids, `[`, character(1L), 1L),
    Campaign = vapply(split_ids, `[`, character(1L), 2L),
    stringsAsFactors = FALSE
  )
  reaches <- sort(unique(metadata$Reach))
  campaigns <- c("Wet_2023", "Dry_2023")
  if (length(reaches) != 3L || !setequal(metadata$Campaign, campaigns)) {
    stop("Expected three reaches and Wet_2023/Dry_2023 campaigns.")
  }
  for (reach in reaches) {
    if (!identical(sort(metadata$Campaign[metadata$Reach == reach]), sort(campaigns))) {
      stop("Each reach must have exactly one wet and one dry community sample.")
    }
  }
  species_columns <- setdiff(names(community), "Sample_ID")
  if (!length(species_columns)) stop("No species-count columns were found.")
  counts <- as.matrix(as.data.frame(lapply(community[species_columns], function(x) {
    numeric_column(x, "Species count")
  }), check.names = FALSE))
  storage.mode(counts) <- "double"
  if (any(counts < 0) || any(abs(counts - round(counts)) > 1e-8) || any(rowSums(counts) <= 0)) {
    stop("Species counts must be nonnegative integers and each sample must be nonempty.")
  }
  rownames(counts) <- metadata$Sample_ID
  colnames(counts) <- species_columns
  metadata$Fish_record_n <- rowSums(counts)
  metadata$Species_richness <- rowSums(counts > 0)
  relative_counts <- counts / rowSums(counts)
  occurrence <- counts > 0
  n <- nrow(counts)
  D_bray <- D_jaccard <- matrix(0, n, n, dimnames = list(metadata$Sample_ID, metadata$Sample_ID))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      D_bray[i, j] <- sum(abs(relative_counts[i, ] - relative_counts[j, ])) / 2
      D_jaccard[i, j] <- 1 - sum(occurrence[i, ] & occurrence[j, ]) /
        sum(occurrence[i, ] | occurrence[j, ])
    }
  }
  if (max(abs(D_bray - t(D_bray))) > 1e-12 || max(abs(D_jaccard - t(D_jaccard))) > 1e-12) {
    stop("Distance-matrix symmetry check failed.")
  }
  write_output(metadata, "Community_metadata.csv")
  write_output(data.frame(Sample_ID = rownames(D_bray), D_bray, check.names = FALSE),
               "Bray_Curtis_distance_matrix.csv")
  write_output(data.frame(Sample_ID = rownames(D_jaccard), D_jaccard, check.names = FALSE),
               "Jaccard_distance_matrix.csv")

  # 2. Descriptive Table S5b --------------------------------------------------
  s5b <- do.call(rbind, lapply(reaches, function(reach) {
    wet <- which(metadata$Reach == reach & metadata$Campaign == "Wet_2023")
    dry <- which(metadata$Reach == reach & metadata$Campaign == "Dry_2023")
    data.frame(Reach_code = reach, Reach = site_name(reach),
               Wet_group = metadata$Sample_ID[wet], Dry_group = metadata$Sample_ID[dry],
               Bray_Curtis = D_bray[wet, dry], Jaccard = D_jaccard[wet, dry],
               Dataset = "Isotope-analyzed fish species counts", stringsAsFactors = FALSE)
  }))
  write_output(s5b, "Table_S5b_Descriptive_Dissimilarity.csv")

  # 3. Enumerate the EXACT restricted allocation sets ------------------------
  # The identity is included ONCE when computing the exact probability.
  # For an optional vegan check it is removed, because vegan adds the observed
  # allocation when calculating the permutation probability.
  permute_vector <- function(x) {
    if (length(x) == 1L) return(matrix(x, nrow = 1L))
    do.call(rbind, lapply(seq_along(x), function(i) {
      remainder <- permute_vector(x[-i])
      cbind(x[i], remainder)
    }))
  }
  identity <- seq_len(n)
  flips <- expand.grid(rep(list(c(FALSE, TRUE)), length(reaches)), KEEP.OUT.ATTRS = FALSE)
  campaign_indices <- t(vapply(seq_len(nrow(flips)), function(k) {
    index <- identity
    for (r in seq_along(reaches)) {
      pair <- which(metadata$Reach == reaches[r])
      if (as.logical(flips[k, r])) index[pair] <- rev(pair)
    }
    as.integer(index)
  }, integer(n)))
  wet_indices <- which(metadata$Campaign == "Wet_2023")
  dry_indices <- which(metadata$Campaign == "Dry_2023")
  wet_permutations <- permute_vector(wet_indices)
  dry_permutations <- permute_vector(dry_indices)
  reach_indices <- matrix(NA_integer_, nrow(wet_permutations) * nrow(dry_permutations), n)
  k <- 0L
  for (a in seq_len(nrow(wet_permutations))) {
    for (b in seq_len(nrow(dry_permutations))) {
      k <- k + 1L
      index <- identity
      index[wet_indices] <- wet_permutations[a, ]
      index[dry_indices] <- dry_permutations[b, ]
      reach_indices[k, ] <- index
    }
  }
  is_identity <- function(p) apply(p, 1L, function(row) all(row == identity))
  stopifnot(nrow(campaign_indices) == 8L, nrow(reach_indices) == 36L,
            sum(is_identity(campaign_indices)) == 1L, sum(is_identity(reach_indices)) == 1L,
            nrow(unique(as.data.frame(campaign_indices))) == 8L,
            nrow(unique(as.data.frame(reach_indices))) == 36L)
  permutation_audit <- rbind(
    data.frame(Test = "Campaign", Allocation = seq_len(8L),
               Observed = is_identity(campaign_indices), campaign_indices),
    data.frame(Test = "Reach", Allocation = seq_len(36L),
               Observed = is_identity(reach_indices), reach_indices)
  )
  names(permutation_audit)[4L:9L] <- paste0("Response_row_", seq_len(n))
  write_output(permutation_audit, "Exact_permutation_allocations.csv")

  # Standard distance-based ONE-FACTOR PERMANOVA sums of squares.
  # No covariate adjustment and no Reach * Campaign interaction is fitted.
  permanova_stat <- function(D, labels) {
    levels <- unique(labels)
    ss_total <- sum(D^2) / (2 * nrow(D))
    ss_residual <- sum(vapply(levels, function(level) {
      ix <- which(labels == level)
      sum(D[ix, ix, drop = FALSE]^2) / (2 * length(ix))
    }, numeric(1L)))
    ss_effect <- ss_total - ss_residual
    df_effect <- length(levels) - 1L
    df_residual <- nrow(D) - length(levels)
    if (ss_total <= 0 || ss_residual <= .Machine$double.eps * ss_total) {
      stop("PERMANOVA has zero total or residual variation; review the data/design.")
    }
    c(SS_effect = ss_effect, SS_residual = ss_residual, SS_total = ss_total,
      df_effect = df_effect, df_residual = df_residual,
      R_squared = ss_effect / ss_total,
      Pseudo_F = (ss_effect / df_effect) / (ss_residual / df_residual))
  }
  tests <- list(Campaign = list(labels = metadata$Campaign, index = campaign_indices),
                Reach = list(labels = metadata$Reach, index = reach_indices))
  distances <- list(Bray_Curtis = D_bray, Jaccard = D_jaccard)
  s5c_rows <- null_rows <- vegan_rows <- list()
  null_index <- 0L
  row_index <- 0L
  vegan_available <- isTRUE(check_vegan) && requireNamespace("vegan", quietly = TRUE)
  for (metric in names(distances)) {
    D <- distances[[metric]]
    for (effect in names(tests)) {
      row_index <- row_index + 1L
      labels <- tests[[effect]]$labels
      permutations <- tests[[effect]]$index
      observed <- permanova_stat(D, labels)
      null_f <- apply(permutations, 1L, function(index) {
        permanova_stat(D[index, index, drop = FALSE], labels)["Pseudo_F"]
      })
      p_exact <- mean(null_f >= observed["Pseudo_F"] - 1e-12)
      s5c_rows[[row_index]] <- data.frame(
        Distance = metric, Effect = effect, Model = paste("Community ~", effect),
        df_effect = unname(observed["df_effect"]),
        df_residual = unname(observed["df_residual"]),
        SS_effect = unname(observed["SS_effect"]),
        SS_residual = unname(observed["SS_residual"]),
        SS_total = unname(observed["SS_total"]),
        R_squared = unname(observed["R_squared"]),
        Pseudo_F = unname(observed["Pseudo_F"]), Exact_P = unname(p_exact),
        Allocations_including_observed = nrow(permutations),
        Nonidentity_label_permutations = nrow(permutations) - 1L,
        Restriction = if (effect == "Campaign") "Campaign swaps within reaches" else "Reach permutations within campaigns",
        Interpretation = "Exploratory conditional one-factor test", stringsAsFactors = FALSE
      )
      null_index <- null_index + 1L
      null_rows[[null_index]] <- data.frame(
        Distance = metric, Effect = effect, Allocation = seq_along(null_f),
        Observed_allocation = is_identity(permutations), Pseudo_F = as.numeric(null_f)
      )
      if (vegan_available) {
        distance_object <- stats::as.dist(D)
        design <- data.frame(group = factor(labels))
        p_nonidentity <- permutations[!is_identity(permutations), , drop = FALSE]
        fit <- vegan::adonis2(distance_object ~ group, data = design,
                              permutations = p_nonidentity, by = NULL, parallel = 1L)
        difference <- max(abs(c(fit$F[1L] - observed["Pseudo_F"],
                                fit$R2[1L] - observed["R_squared"],
                                fit$`Pr(>F)`[1L] - p_exact)))
        vegan_rows[[row_index]] <- data.frame(
          Distance = metric, Effect = effect, Status = if (difference < 1e-8) "Passed" else "Failed",
          Vegan_pseudo_F = fit$F[1L], Vegan_R_squared = fit$R2[1L],
          Vegan_P = fit$`Pr(>F)`[1L], Maximum_difference = unname(difference),
          stringsAsFactors = FALSE
        )
        if (difference >= 1e-8) stop("Base-R and vegan results differ; inspect before reporting.")
      }
    }
  }
  s5c <- do.call(rbind, s5c_rows)
  rownames(s5c) <- NULL
  write_output(s5c, "Table_S5c_PERMANOVA.csv")
  write_output(do.call(rbind, null_rows), "PERMANOVA_null_distributions.csv")
  if (vegan_available) {
    vegan_check <- do.call(rbind, vegan_rows)
  } else {
    vegan_check <- data.frame(Status = "Not run", Reason = if (check_vegan)
      "Optional vegan package is not installed; base-R analysis is complete without it."
      else "Optional vegan cross-check disabled.", stringsAsFactors = FALSE)
  }
  write_output(vegan_check, "Vegan_crosscheck.csv")

  # Reference check applies ONLY to the original, unchanged community file.
  community_md5 <- unname(tools::md5sum(community_path))
  reference_check <- data.frame(
    Distance = c("Bray_Curtis", "Bray_Curtis", "Jaccard", "Jaccard"),
    Effect = rep(c("Campaign", "Reach"), 2L),
    Reference_R_squared = c(0.25078541583994834, 0.49991704647924756,
                            0.26186451490795026, 0.505628102596766),
    Reference_pseudo_F = c(1.3389243676888922, 1.499502361437227,
                           1.4190593472161508, 1.534153049311633),
    Reference_P = rep(c(0.25, 1/6), 2L), stringsAsFactors = FALSE
  )
  reference_check$Actual_R_squared <- s5c$R_squared
  reference_check$Actual_pseudo_F <- s5c$Pseudo_F
  reference_check$Actual_P <- s5c$Exact_P
  reference_check$Original_input_unchanged <- community_md5 == "98a57525ca482ba3b215def2e12f0211"
  reference_check$Matches_reference <- abs(reference_check$Actual_R_squared - reference_check$Reference_R_squared) < 1e-10 &
    abs(reference_check$Actual_pseudo_F - reference_check$Reference_pseudo_F) < 1e-10 &
    abs(reference_check$Actual_P - reference_check$Reference_P) < 1e-10
  if (reference_check$Original_input_unchanged[1L] && !all(reference_check$Matches_reference)) {
    stop("Unchanged community input did not reproduce reference results.")
  }
  write_output(reference_check, "Reference_results_check.csv")

  # 4. Validate reconstructed CPUE record groups ------------------------------
  events <- read_input(event_path)
  required_columns(events, c("Event_ID", "Reach_code", "Gear", "Campaign", "Catch_n", "Effort_net_h"),
                   "CPUE event input")
  if (!nrow(events) || anyNA(events$Event_ID) || anyDuplicated(events$Event_ID)) {
    stop("Event_ID must be present and unique for each reconstructed group.")
  }
  if (anyNA(events$Reach_code) || anyNA(events$Gear) || anyNA(events$Campaign) ||
      any(!nzchar(events$Event_ID)) || any(!events$Reach_code %in% reaches) ||
      any(!events$Campaign %in% campaigns) || any(!events$Gear %in% c("Gill net", "Trap net"))) {
    stop("CPUE identifiers must match the study reaches, two campaigns and Gill net/Trap net gears.")
  }
  events$Catch_n <- numeric_column(events$Catch_n, "Catch_n")
  events$Effort_net_h <- numeric_column(events$Effort_net_h, "Effort_net_h")
  if (any(events$Catch_n < 0) || any(abs(events$Catch_n - round(events$Catch_n)) > 1e-8) ||
      any(events$Effort_net_h <= 0)) stop("Catch must be a nonnegative integer; effort must be positive.")
  effort_discrepancies <- NA_integer_
  if (all(c("Number_nets", "Deployment_duration_h") %in% names(events))) {
    nets <- numeric_column(events$Number_nets, "Number_nets")
    hours <- numeric_column(events$Deployment_duration_h, "Deployment_duration_h")
    if (any(nets <= 0) || any(hours <= 0) || any(abs(nets - round(nets)) > 1e-8)) {
      stop("Recorded net numbers must be positive integers and durations must be positive.")
    }
    events$Effort_recalculated_net_h <- nets * hours
    events$Effort_difference_net_h <- events$Effort_net_h - events$Effort_recalculated_net_h
    effort_discrepancies <- sum(abs(events$Effort_difference_net_h) > 1e-8)
    if (effort_discrepancies > 0L) stop("Effort differs from recorded net numbers times duration; reconcile inputs.")
  }
  outside_windows <- NA_integer_
  outside_catch <- NA_real_
  if ("Date" %in% names(events)) {
    date_strings <- trimws(as.character(events$Date))
    events$Source_date_parsed <- suppressWarnings(as.Date(date_strings, format = "%Y-%m-%d"))
    na_idx <- is.na(events$Source_date_parsed)
    if (any(na_idx)) {
      events$Source_date_parsed[na_idx] <- suppressWarnings(as.Date(date_strings[na_idx], format = "%m/%d/%Y"))
    }
    na_idx <- is.na(events$Source_date_parsed)
    if (any(na_idx)) {
      events$Source_date_parsed[na_idx] <- suppressWarnings(as.Date(date_strings[na_idx], format = "%Y/%m/%d"))
    }
    if (anyNA(events$Source_date_parsed)) stop("Date could not be parsed; expected YYYY-MM-DD (or M/D/YYYY).")
    events$Date <- format(events$Source_date_parsed, "%Y-%m-%d")
    window_ix <- match(events$Campaign, stated_windows$Campaign)
    if (anyNA(window_ix)) stop("No stated date window found for an event campaign.")
    events$Outside_stated_campaign_window <- events$Source_date_parsed < stated_windows$Start_date[window_ix] |
      events$Source_date_parsed > stated_windows$End_date[window_ix]
    outside_windows <- sum(events$Outside_stated_campaign_window)
    outside_catch <- sum(events$Catch_n[events$Outside_stated_campaign_window])
  }
  write_output(events, "CPUE_record_audit.csv")

  # Current data explicitly withhold unverified CPUE. A blank analysis column
  # must never be bypassed by recomputing Catch_n / Effort_net_h silently.
  cpue_check <- hr_cpue_ready(events)
  if (!isTRUE(cpue_check$ready) && !isTRUE(run_conditional_cpue)) {
    write_output(events, "CPUE_record_audit.csv")
    write_output(data.frame(Status = "BLOCKED", Reason = cpue_check$reasons,
                           stringsAsFactors = FALSE), "CPUE_VALIDATION_REQUIRED.csv")
    s5b_word <- data.frame(Reach = s5b$Reach,
      Bray_Curtis = format_number(s5b$Bray_Curtis),
      Jaccard = format_number(s5b$Jaccard), stringsAsFactors = FALSE)
    s5c_word <- data.frame(Distance = s5c$Distance, Effect = s5c$Effect,
      Effect_df = s5c$df_effect, Residual_df = s5c$df_residual,
      R_squared = format_number(s5c$R_squared), Pseudo_F = format_number(s5c$Pseudo_F),
      Exact_P = format_number(s5c$Exact_P),
      Allocations_including_observed = s5c$Allocations_including_observed,
      stringsAsFactors = FALSE)
    write_output(s5b_word, "Table_S5b_Word_ready.csv")
    write_output(s5c_word, "Table_S5c_Word_ready.csv")
    header_map <- data.frame(Column_position = 8L:11L,
      Current_heading = c("p_bar", "delta_p", "z wet", "z dry"),
      Correct_heading = c("z wet", "z dry", "p_bar", "delta-p"),
      Action = "Rename headings only; preserve body values", stringsAsFactors = FALSE)
    write_output(header_map, "Table_S3e_Header_Correction.csv")
    write_output(data.frame(Item = c("Community_samples", "Isotope_fish_records", "Recorded_event_groups", "Recorded_event_catch", "CPUE_status"),
      Value = c(nrow(counts), sum(counts), nrow(events), sum(events$Catch_n), "BLOCKED"),
      stringsAsFactors = FALSE), "Input_audit_summary.csv")
    writeLines(c("Community results from the current run; CPUE was withheld.",
      capture.output(print(s5c_word, row.names = FALSE)),
      "CPUE source verification still required:", cpue_check$reasons),
      file.path(output_dir, "Manuscript_replacement_text.txt"))
    writeLines(c("Community analyses completed. CPUE estimates were not calculated.",
      cpue_check$reasons,
      "Verify deployment identities, complete catches, zero-catch coverage and dates against field logs."),
      file.path(output_dir, "CPUE_STATUS.txt"))
    capture.output(utils::sessionInfo(), file = file.path(output_dir, "SessionInfo.txt"))
    result <- list(S3e_header_map = header_map, S5b = s5b, S5c = s5c, CPUE_status = "BLOCKED",
      CPUE_reasons = cpue_check$reasons, Record_audit = events)
    if (save_objects) saveRDS(result, file.path(output_dir, "Reviewer4_analysis_objects.rds"))
    message("Reviewer 4 community analyses completed; CPUE remains withheld. See CPUE_STATUS.txt.")
    return(invisible(result))
  }

  # 5. Conditional event-group bootstrap for S9 ------------------------------
  # R's RNG and seed are explicit. Every replicate resamples complete groups,
  # retaining catch and effort together, then recomputes their pooled ratio.
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(seed)
  probs <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)
  quantile_pair <- function(x) as.numeric(stats::quantile(x, probs = probs, names = FALSE, type = 7))
  s9a_rows <- s9b_rows <- bootstrap_draws <- list()
  a_index <- b_index <- 0L
  for (reach in reaches) {
    gears <- sort(unique(events$Gear[events$Reach_code == reach]))
    if (!setequal(gears, c("Gill net", "Trap net"))) stop("Each study reach must have both gear categories.")
    for (gear in gears) {
      point <- setNames(rep(NA_real_, 2L), campaigns)
      draws <- setNames(vector("list", 2L), campaigns)
      group_n <- setNames(integer(2L), campaigns)
      for (campaign in campaigns) {
        cell <- events[events$Reach_code == reach & events$Gear == gear & events$Campaign == campaign, , drop = FALSE]
        if (!nrow(cell)) stop("Missing CPUE reach-gear-campaign cell.")
        group_n[campaign] <- nrow(cell)
        point[campaign] <- sum(cell$Catch_n) / sum(cell$Effort_net_h)
        interval <- c(NA_real_, NA_real_)
        if (nrow(cell) >= 2L) {
          draws[[campaign]] <- replicate(n_boot, {
            index <- sample.int(nrow(cell), size = nrow(cell), replace = TRUE)
            sum(cell$Catch_n[index]) / sum(cell$Effort_net_h[index])
          })
          interval <- quantile_pair(draws[[campaign]])
        }
        a_index <- a_index + 1L
        s9a_rows[[a_index]] <- data.frame(
          Reach_code = reach, Reach = site_name(reach), Gear = gear, Campaign = campaign,
          N_reconstructed_event_groups = nrow(cell), Recorded_catch_n = sum(cell$Catch_n),
          N_zero_catch_groups = sum(cell$Catch_n == 0), Inferred_net_hours = sum(cell$Effort_net_h),
          CPUE = unname(point[campaign]), Conditional_CI_low = interval[1L], Conditional_CI_high = interval[2L],
          Confidence_level = ci_level, Bootstrap_replicates = if (nrow(cell) >= 2L) n_boot else 0L,
          Interval_status = if (nrow(cell) >= 2L) "Conditional record-group interval" else "Not estimable: fewer than two groups",
          stringsAsFactors = FALSE
        )
      }
      change_point <- if (point["Wet_2023"] > 0) 100 * (point["Dry_2023"] / point["Wet_2023"] - 1) else NA_real_
      interval <- c(NA_real_, NA_real_)
      undefined <- NA_integer_
      change_draws <- numeric(0L)
      status <- "Not estimable: fewer than two groups in a campaign"
      if (all(group_n >= 2L)) {
        undefined <- sum(draws[["Wet_2023"]] == 0)
        if (point["Wet_2023"] > 0 && undefined == 0L) {
          change_draws <- 100 * (draws[["Dry_2023"]] / draws[["Wet_2023"]] - 1)
          interval <- quantile_pair(change_draws)
          status <- "Conditional record-group interval"
        } else {
          status <- "Not estimable: zero wet denominator in observation or resampling"
        }
      }
      b_index <- b_index + 1L
      s9b_rows[[b_index]] <- data.frame(
        Reach_code = reach, Reach = site_name(reach), Gear = gear,
        Wet_CPUE = unname(point["Wet_2023"]), Dry_CPUE = unname(point["Dry_2023"]),
        Percent_change = unname(change_point), Conditional_CI_low_percent = interval[1L],
        Conditional_CI_high_percent = interval[2L], Undefined_wet_denominator_replicates = undefined,
        Confidence_level = ci_level, Bootstrap_replicates = if (all(group_n >= 2L)) n_boot else 0L,
        Interval_status = status, stringsAsFactors = FALSE
      )
      bootstrap_draws[[paste(reach, gear, sep = "|")]] <- list(CPUE = draws, Percent_change = change_draws)
    }
  }
  s9a <- do.call(rbind, s9a_rows)
  s9b <- do.call(rbind, s9b_rows)
  rownames(s9a) <- rownames(s9b) <- NULL
  stopifnot(sum(s9a$N_reconstructed_event_groups) == nrow(events),
            abs(sum(s9a$Recorded_catch_n) - sum(events$Catch_n)) < 1e-8,
            abs(sum(s9a$Inferred_net_hours) - sum(events$Effort_net_h)) < 1e-8)
  write_output(s9a, "Table_S9a_Conditional_CPUE.csv")
  write_output(s9b, "Table_S9b_Conditional_Change.csv")

  # 6. Formatted copies for easy insertion into Word --------------------------
  s5b_word <- data.frame(Reach = s5b$Reach, Bray_Curtis = format_number(s5b$Bray_Curtis),
                        Jaccard = format_number(s5b$Jaccard), stringsAsFactors = FALSE)
  s5c_word <- data.frame(Distance = s5c$Distance, Effect = s5c$Effect,
                        Effect_df = s5c$df_effect, Residual_df = s5c$df_residual,
                        R_squared = format_number(s5c$R_squared), Pseudo_F = format_number(s5c$Pseudo_F),
                        Exact_P = format_number(s5c$Exact_P),
                        Allocations_including_observed = s5c$Allocations_including_observed,
                        stringsAsFactors = FALSE)
  s9a_word <- data.frame(Reach = s9a$Reach, Gear = s9a$Gear,
                        Campaign = ifelse(s9a$Campaign == "Wet_2023", "Wet", "Dry"),
                        Recorded_event_groups_n = s9a$N_reconstructed_event_groups,
                        Catch_n = s9a$Recorded_catch_n, Inferred_net_hours = s9a$Inferred_net_hours,
                        CPUE = format_number(s9a$CPUE),
                        Conditional_interval = interval_text(s9a$Conditional_CI_low, s9a$Conditional_CI_high),
                        stringsAsFactors = FALSE)
  s9b_word <- data.frame(Reach = s9b$Reach, Gear = s9b$Gear,
                        Change_percent = ifelse(is.na(s9b$Percent_change), "Not estimable", sprintf("%+.1f", s9b$Percent_change)),
                        Conditional_interval_percent = interval_text(s9b$Conditional_CI_low_percent,
                                                                      s9b$Conditional_CI_high_percent, 1L),
                        stringsAsFactors = FALSE)
  write_output(s5b_word, "Table_S5b_Word_ready.csv")
  write_output(s5c_word, "Table_S5c_Word_ready.csv")
  write_output(s9a_word, "Table_S9a_Word_ready.csv")
  write_output(s9b_word, "Table_S9b_Word_ready.csv")
  header_map <- data.frame(Column_position = 8L:11L,
                           Current_heading = c("p_bar", "delta_p", "z wet", "z dry"),
                           Correct_heading = c("z wet", "z dry", "p_bar", "delta_p"),
                           Action = "Rename heading in all 12 S3e panels; preserve body values",
                           stringsAsFactors = FALSE)
  write_output(header_map, "Table_S3e_Header_Correction.csv")

  # 7. Manuscript text generated FROM THIS R RUN -----------------------------
  get_test <- function(metric, effect) s5c[s5c$Distance == metric & s5c$Effect == effect, , drop = FALSE]
  bc_campaign <- get_test("Bray_Curtis", "Campaign")
  bc_reach <- get_test("Bray_Curtis", "Reach")
  ja_campaign <- get_test("Jaccard", "Campaign")
  ja_reach <- get_test("Jaccard", "Reach")
  ci_percent <- 100 * ci_level
  methods_permanova <- paste(
    "Campaign and reach effects were assessed in separate one-factor PERMANOVA models.",
    "The campaign test enumerated all eight campaign-label allocations obtained by swapping wet and dry labels within each reach.",
    "The reach test enumerated all 36 reach-label allocations within the two campaigns.",
    "Exact probabilities included the observed allocation. These six-group tests were interpreted as exploratory and conditional on the exchangeability assumptions of each restriction scheme; no reach-by-campaign interaction was estimated."
  )
  results_permanova <- sprintf(
    paste0("Exploratory PERMANOVA of Bray-Curtis dissimilarities attributed %.1f%% of variation to campaign ",
           "(pseudo-F = %.3f, exact P %s) and %.1f%% to reach (pseudo-F = %.3f, exact P %s) in separate one-factor models. ",
           "Corresponding Jaccard estimates were %.1f%% and %.1f%%, with exact P %s and P %s, respectively (Table S5c)."),
    100 * bc_campaign$R_squared, bc_campaign$Pseudo_F, format_p(bc_campaign$Exact_P),
    100 * bc_reach$R_squared, bc_reach$Pseudo_F, format_p(bc_reach$Exact_P),
    100 * ja_campaign$R_squared, 100 * ja_reach$R_squared,
    format_p(ja_campaign$Exact_P), format_p(ja_reach$Exact_P)
  )
  methods_cpue <- sprintf(
    paste0("Fishing effort was reconstructed from reported net numbers multiplied by fishing duration and expressed as net-hours. ",
           "Gill nets and trap nets were analyzed separately. Pooled CPUE was total recorded catch divided by total inferred net-hours. ",
           "Wet-to-dry change was 100 x (CPUEdry/CPUEwet - 1). Conditional uncertainty was summarized using %d resamples of complete reconstructed event groups ",
           "within each reach-by-campaign-by-gear stratum, retaining catch and effort together and recalculating their pooled ratio. ",
           "R's Mersenne-Twister generator and Rejection sampling method were used with seed %d. Percentile %.0f%% intervals used quantile type 7. ",
           "Intervals were conditional on the recorded event groups and reported effort and did not account for uncertain deployment identity, unrecorded zero-catch deployments or effort-recording errors."),
    n_boot, seed, ci_percent
  )
  gear_summaries <- vapply(seq_len(nrow(s9b)), function(i) {
    r <- s9b[i, ]
    sprintf("%s %s: %+.1f%% (conditional %.0f%% interval: %s to %s%%)",
            r$Reach, r$Gear, r$Percent_change, ci_percent,
            format_number(r$Conditional_CI_low_percent, 1L), format_number(r$Conditional_CI_high_percent, 1L))
  }, character(1L))
  supported <- which(!is.na(s9b$Conditional_CI_low_percent) &
                      (s9b$Conditional_CI_low_percent > 0 | s9b$Conditional_CI_high_percent < 0))
  supported_sentence <- if (length(supported)) {
    paste0("Conditional intervals excluded zero for ",
           paste(paste(s9b$Reach[supported], s9b$Gear[supported]), collapse = "; "), ".")
  } else "No estimable conditional interval excluded zero change."
  results_cpue <- paste(
    "Recorded CPUE showed reach- and gear-specific wet-to-dry changes (Table S9):",
    paste0(paste(gear_summaries, collapse = "; "), "."), supported_sentence,
    "These record-based summaries are conditional on the reconstructed event groups and reported effort; they do not establish population-abundance changes or a repeatable seasonal response across years."
  )
  zero_note <- if (sum(events$Catch_n == 0) == 0L) {
    "No zero-catch groups occur in the supplied records; completeness of zero-catch coverage is unverified."
  } else {
    sprintf("The supplied records include %d zero-catch groups; completeness of deployment coverage remains unverified.", sum(events$Catch_n == 0))
  }
  provenance_note <- paste0(
    "Counts refer to ", nrow(events), " reconstructed catch-record groups containing ", sum(events$Catch_n),
    " recorded fish, not independently verified physical deployments. ", zero_note,
    if (!is.na(outside_windows)) sprintf(" %d groups containing %g fish have source dates outside the manuscript's stated isotope/environmental windows. These records were retained for audit; dates and physical deployments require field-log verification.", outside_windows, outside_catch) else ""
  )
  software_note <- paste0("The supplementary community PERMANOVA reconstructions and conditional CPUE resampling were implemented in R ",
                          getRversion(), ". ",
                          if (vegan_available) paste0("PERMANOVA calculations were also checked against vegan::adonis2, vegan version ",
                                                     utils::packageVersion("vegan"), ".")
                          else "These calculations used base R and required no additional packages.")
  text <- c(
    "REVIEWER COMMENT 4 - WORD REPLACEMENTS FROM THE CURRENT R RUN", "",
    "S3e header correction", "Rename columns 8-11 to: z wet, z dry, p-bar, delta-p. Apply to all 12 panels.",
    "This is a heading correction; this script does not recalculate all species isotope decomposition terms.", "",
    "Community PERMANOVA methods", methods_permanova, "",
    "Community PERMANOVA results", results_permanova, "",
    "CPUE methods", methods_cpue, "", "CPUE results", results_cpue, "",
    "S9 provenance note", provenance_note, "",
    "Software wording", software_note, "",
    "Availability wording", "Analysis-ready community counts, recorded catch-event data and executable R code for Tables S5b, S5c and S9 are supplied as supplementary files.", "",
    "Scope", "Retain the original brms, SIBER and environmental source files for full-manuscript reproducibility. Updating Word headings and citations does not require running Python or R.",
    "Use confidence intervals from this R run consistently; do not copy the earlier Python intervals into these R tables."
  )
  writeLines(enc2utf8(text), file.path(output_dir, "Manuscript_replacement_text.txt"), useBytes = TRUE)

  # 8. Analysis settings, source audit and reusable R objects ------------------
  audit <- data.frame(
    Item = c("Community_samples", "Isotope_fish_records", "Species_columns",
             "Reconstructed_CPUE_event_groups", "Recorded_catch_n", "Zero_catch_groups",
             "Effort_formula_discrepancies", "Groups_outside_stated_windows", "Catches_outside_stated_windows"),
    Value = c(nrow(counts), sum(counts), ncol(counts), nrow(events), sum(events$Catch_n),
              sum(events$Catch_n == 0), effort_discrepancies, outside_windows, outside_catch),
    stringsAsFactors = FALSE
  )
  write_output(audit, "Input_audit_summary.csv")
  settings <- list(
    Created_at = format(Sys.time(), tz = "UTC", usetz = TRUE), R_version = as.character(getRversion()),
    Input_directory = input_dir, Output_directory = output_dir,
    Community_filename = community_filename, Community_md5 = community_md5,
    Event_filename = event_filename, Event_md5 = unname(tools::md5sum(event_path)),
    Bootstrap_replicates = n_boot, Bootstrap_seed = seed, Confidence_level = ci_level,
    RNG = RNGkind(), Bootstrap_unit = "Complete reconstructed catch-record group",
    Quantile_type = 7L, Group_order = "Reach and gear sorted; Wet_2023 then Dry_2023; original source-row order within cell",
    Permutation_allocations_including_observed = c(Campaign = 8L, Reach = 36L),
    Models = c("Community ~ Campaign", "Community ~ Reach"),
    Vegan_crosscheck_performed = vegan_available,
    CPUE_scope = provenance_note,
    Inference_limit = "One aggregate community sample per reach-campaign; CPUE intervals conditional on compiled record groups"
  )
  dput(settings, file = file.path(output_dir, "Analysis_settings.R"))
  capture.output(utils::sessionInfo(), file = file.path(output_dir, "SessionInfo.txt"))
  result <- list(S3e_header_map = header_map, S5b = s5b, S5c = s5c, S9a = s9a, S9b = s9b,
                 Community_metadata = metadata, Community_counts = counts,
                 Bray_Curtis = D_bray, Jaccard = D_jaccard,
                 Record_audit = events, Input_audit = audit,
                 Permutation_allocations = permutation_audit, Settings = settings)
  if (save_objects) saveRDS(list(results = result, conditional_bootstrap_draws = bootstrap_draws),
                             file = file.path(output_dir, "Reviewer4_analysis_objects.rds"))

  cat("\n--- Table S5c: exact community PERMANOVA ---\n")
  print(s5c_word, row.names = FALSE)
  cat("\n--- Table S9a: conditional record-based CPUE ---\n")
  print(s9a_word, row.names = FALSE)
  cat("\n--- Table S9b: conditional percentage changes ---\n")
  print(s9b_word, row.names = FALSE)
  cat("\n", provenance_note, "\n", sep = "")
  cat("\nFinished. Output folder:\n", output_dir, "\n", sep = "")
  cat("Use the Word_ready CSV files and Manuscript_replacement_text.txt to revise the paper.\n")
  cat("Use intervals from this R run consistently; Monte Carlo bounds can differ from NumPy's.\n")
  invisible(result)
}

# ----------------------------- RUN EVERYTHING ------------------------------
reviewer4_results <- run_reviewer4_analysis(
  input_dir = INPUT_DIR,
  output_dir = OUTPUT_DIR,
  community_filename = COMMUNITY_FILENAME,
  event_filename = EVENT_FILENAME,
  n_boot = BOOTSTRAP_REPLICATES,
  seed = BOOTSTRAP_SEED,
  ci_level = CONFIDENCE_LEVEL,
  check_vegan = CHECK_WITH_VEGAN_IF_INSTALLED,
  save_objects = SAVE_ANALYSIS_OBJECTS,
  run_conditional_cpue = RUN_CONDITIONAL_CPUE,
  stated_windows = STATED_CAMPAIGN_WINDOWS
)

# =============================================================================
# LOWER HAN RIVER FISH COMMUNITY COMPOSITION
# NMDS, PERMANOVA, JACCARD, AND BRAY-CURTIS ANALYSIS
# CORRECTED FINAL READY-TO-USE R SCRIPT
#
# THIS IS NOT A BAYESIAN ANALYSIS.
#
# PURPOSE
#   Uses the 609 individual fish records to construct six reach-campaign
#   species-composition samples and then:
#     1. calculates Bray-Curtis dissimilarity from relative abundance;
#     2. calculates Jaccard dissimilarity from presence-absence;
#     3. performs two-dimensional NMDS for both dissimilarities;
#     4. performs blocked exploratory PERMANOVA tests for Campaign and Reach;
#     5. checks multivariate dispersion with PERMDISP; and
#     6. summarizes all pairwise and within-reach wet-dry dissimilarities.
#
# CRITICAL DESIGN LIMITATION
#   There is one aggregated community sample for each Reach x Campaign cell.
#   Therefore:
#     - the Reach x Campaign interaction cannot be tested;
#     - Campaign is tested using permutations blocked within Reach;
#     - Reach is tested using permutations blocked within Campaign;
#     - only 8 campaign permutations and 36 reach permutations are admissible;
#     - PERMANOVA results are exploratory and have coarse p-value resolution.
#
# INPUT WORKBOOK
#   Han_River_Fish_Community_Composition_R_Dataset.xlsx
#
# PRIMARY R INPUT SHEETS
#   Individual_Data
#   Sample_Metadata
#   Pairwise_Design
#
# REQUIRED R PACKAGES
#   vegan, permute, readxl, dplyr, ggplot2
#
# HOW TO RUN
#   1. Select Session -> Restart R.
#   2. Open this script.
#   3. Click Source once.
#   4. Select the Excel workbook if requested.
# =============================================================================


# =============================================================================
# 1. USER SETTINGS
# =============================================================================

SEED <- 20260731L

NMDS_DIMENSIONS <- 2L
NMDS_TRY <- 50L
NMDS_TRYMAX <- 999L

PERMUTATIONS_REQUESTED <- 9999L

set.seed(SEED)
graphics.off()


# =============================================================================
# 2. PACKAGE CHECK
# =============================================================================

required_packages <- c(
  "vegan",
  "permute",
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
  library(vegan)
  library(permute)
  library(readxl)
  library(dplyr)
  library(ggplot2)
})


# =============================================================================
# 3. MAIN ANALYSIS WRAPPER
# =============================================================================

run_fish_community_analysis <- function() {

  # ---------------------------------------------------------------------------
  # 3.1 Locate the input workbook
  # ---------------------------------------------------------------------------

  workbook_pattern <- paste0(
    "^Han_River_Fish_Community_Composition_",
    "R_Dataset.*[.]xlsx$"
  )

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
      "\nSelect the fish-community Excel workbook."
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
  # 3.2 Create result folders and status files
  # ---------------------------------------------------------------------------

  RESULTS_DIR <- file.path(
    dirname(INPUT_FILE),
    "Fish_community_composition_results"
  )

  TABLE_DIR <- file.path(
    RESULTS_DIR,
    "Tables"
  )

  FIGURE_DIR <- file.path(
    RESULTS_DIR,
    "Figures"
  )

  OBJECT_DIR <- file.path(
    RESULTS_DIR,
    "R_Objects"
  )

  LOG_DIR <- file.path(
    RESULTS_DIR,
    "Logs"
  )

  for (directory_path in c(
    RESULTS_DIR,
    TABLE_DIR,
    FIGURE_DIR,
    OBJECT_DIR,
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

  limitation_file <- file.path(
    LOG_DIR,
    "PERMANOVA_DESIGN_LIMITATION.txt"
  )

  unlink(
    c(
      completed_file,
      error_file
    )
  )

  writeLines(
    c(
      "Fish-community composition analysis is running.",
      paste("Started:", Sys.time()),
      paste("Input workbook:", INPUT_FILE)
    ),
    running_file
  )

  writeLines(
    c(
      "PERMANOVA DESIGN LIMITATION",
      "",
      "Only one aggregated community sample is available for each reach-campaign cell.",
      "The Reach x Campaign interaction is not estimable.",
      "Campaign is tested with permutations restricted within reach.",
      "Reach is tested with permutations restricted within campaign.",
      "The tests have only 8 and 36 admissible permutations, respectively.",
      "Interpret p-values as exploratory and report R2 effect sizes and dissimilarities."
    ),
    limitation_file
  )


  # ---------------------------------------------------------------------------
  # 3.3 Read and validate workbook sheets
  # ---------------------------------------------------------------------------

  available_sheets <- readxl::excel_sheets(
    INPUT_FILE
  )

  required_sheets <- c(
    "Individual_Data",
    "Sample_Metadata",
    "Pairwise_Design"
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
    sheet = "Individual_Data"
  )

  sample_metadata <- readxl::read_excel(
    INPUT_FILE,
    sheet = "Sample_Metadata"
  )

  pairwise_design <- readxl::read_excel(
    INPUT_FILE,
    sheet = "Pairwise_Design"
  )

  fish_data <- as.data.frame(
    fish_data,
    stringsAsFactors = FALSE
  )

  sample_metadata <- as.data.frame(
    sample_metadata,
    stringsAsFactors = FALSE
  )

  pairwise_design <- as.data.frame(
    pairwise_design,
    stringsAsFactors = FALSE
  )

  required_fish_columns <- c(
    "Fish_id",
    "Reach_code",
    "Campaign",
    "Species_code",
    "Feeding_guild",
    "Sample_ID"
  )

  required_metadata_columns <- c(
    "Sample_ID",
    "Reach_code",
    "Reach_name",
    "Campaign",
    "Campaign_label",
    "n_fish",
    "species_richness"
  )

  required_pair_columns <- c(
    "Pair_ID",
    "Sample_1",
    "Sample_2",
    "Comparison_type",
    "Manuscript_priority",
    "Same_reach",
    "Same_campaign",
    "Bray_Curtis_relative_abundance",
    "Jaccard_presence_absence"
  )

  missing_fish_columns <- setdiff(
    required_fish_columns,
    names(fish_data)
  )

  missing_metadata_columns <- setdiff(
    required_metadata_columns,
    names(sample_metadata)
  )

  missing_pair_columns <- setdiff(
    required_pair_columns,
    names(pairwise_design)
  )

  if (length(missing_fish_columns) > 0L) {
    stop(
      paste0(
        "Individual_Data is missing column(s): ",
        paste(missing_fish_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (length(missing_metadata_columns) > 0L) {
    stop(
      paste0(
        "Sample_Metadata is missing column(s): ",
        paste(missing_metadata_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (length(missing_pair_columns) > 0L) {
    stop(
      paste0(
        "Pairwise_Design is missing column(s): ",
        paste(missing_pair_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (anyNA(
    fish_data[
      ,
      required_fish_columns
    ]
  )) {
    stop(
      "Missing values were found in required individual-level variables.",
      call. = FALSE
    )
  }

  if (anyNA(
    sample_metadata[
      ,
      required_metadata_columns
    ]
  )) {
    stop(
      "Missing values were found in Sample_Metadata.",
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

  expected_sample_ids <- c(
    "HC__Wet_2023",
    "HC__Dry_2023",
    "QJ__Wet_2023",
    "QJ__Dry_2023",
    "SY__Wet_2023",
    "SY__Dry_2023"
  )

  if (nrow(
    fish_data
  ) != 609L) {
    stop(
      paste0(
        "Expected 609 individual fish, but found ",
        nrow(fish_data),
        "."
      ),
      call. = FALSE
    )
  }

  observed_sample_ids <- unique(
    as.character(
      fish_data$Sample_ID
    )
  )

  if (!setequal(
    observed_sample_ids,
    expected_sample_ids
  )) {
    stop(
      paste0(
        "Unexpected Sample_ID values. Observed: ",
        paste(
          observed_sample_ids,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }

  sample_metadata <- sample_metadata[
    match(
      expected_sample_ids,
      sample_metadata$Sample_ID
    ),
    ,
    drop = FALSE
  ]

  if (anyNA(
    sample_metadata$Sample_ID
  )) {
    stop(
      "Sample_Metadata did not contain all expected sample identifiers.",
      call. = FALSE
    )
  }

  fish_data$Sample_ID <- factor(
    as.character(
      fish_data$Sample_ID
    ),
    levels = expected_sample_ids
  )

  fish_data$Species_code <- factor(
    as.character(
      fish_data$Species_code
    ),
    levels = sort(
      unique(
        as.character(
          fish_data$Species_code
        )
      )
    )
  )

  sample_metadata$Sample_ID <- factor(
    as.character(
      sample_metadata$Sample_ID
    ),
    levels = expected_sample_ids
  )

  sample_metadata$Reach_code <- factor(
    as.character(
      sample_metadata$Reach_code
    ),
    levels = expected_reaches
  )

  sample_metadata$Campaign <- factor(
    as.character(
      sample_metadata$Campaign
    ),
    levels = expected_campaigns
  )

  if (
    anyNA(
      sample_metadata$Reach_code
    ) ||
    anyNA(
      sample_metadata$Campaign
    )
  ) {
    stop(
      "Unexpected reach or campaign values were found.",
      call. = FALSE
    )
  }

  if (nlevels(
    fish_data$Species_code
  ) != 65L) {
    stop(
      paste0(
        "Expected 65 species, but found ",
        nlevels(
          fish_data$Species_code
        ),
        "."
      ),
      call. = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # 3.4 Construct community matrices
  # ---------------------------------------------------------------------------

  # Build the contingency table first, then explicitly convert it to a
  # plain numeric matrix. Using as.matrix(table(...)) can retain the "table"
  # class in some R versions. vegan::vegdist() then calls isSymmetric(), for
  # which no "table" method exists, producing:
  #   no applicable method for 'isSymmetric' applied to an object of class "table"
  count_table <- table(
    factor(
      fish_data$Sample_ID,
      levels = expected_sample_ids
    ),
    factor(
      fish_data$Species_code,
      levels = levels(
        fish_data$Species_code
      )
    )
  )

  count_matrix <- matrix(
    data = as.numeric(
      count_table
    ),
    nrow = length(
      expected_sample_ids
    ),
    ncol = length(
      levels(
        fish_data$Species_code
      )
    ),
    byrow = FALSE,
    dimnames = list(
      expected_sample_ids,
      levels(
        fish_data$Species_code
      )
    )
  )

  storage.mode(
    count_matrix
  ) <- "double"

  relative_abundance_working <- vegan::decostand(
    count_matrix,
    method = "total",
    MARGIN = 1
  )

  # Rebuild as a plain matrix to remove all package- or source-specific
  # classes and attributes before distance calculations.
  relative_abundance_matrix <- matrix(
    data = as.numeric(
      relative_abundance_working
    ),
    nrow = nrow(
      count_matrix
    ),
    ncol = ncol(
      count_matrix
    ),
    byrow = FALSE,
    dimnames = dimnames(
      count_matrix
    )
  )

  presence_absence_matrix <- matrix(
    data = as.numeric(
      count_matrix > 0
    ),
    nrow = nrow(
      count_matrix
    ),
    ncol = ncol(
      count_matrix
    ),
    byrow = FALSE,
    dimnames = dimnames(
      count_matrix
    )
  )

  storage.mode(
    relative_abundance_matrix
  ) <- "double"

  storage.mode(
    presence_absence_matrix
  ) <- "double"

  if (
    inherits(
      count_matrix,
      "table"
    ) ||
    inherits(
      relative_abundance_matrix,
      "table"
    ) ||
    inherits(
      presence_absence_matrix,
      "table"
    )
  ) {
    stop(
      "Internal matrix conversion failed: a community matrix still has class 'table'.",
      call. = FALSE
    )
  }

  if (
    !is.matrix(
      count_matrix
    ) ||
    !is.matrix(
      relative_abundance_matrix
    ) ||
    !is.matrix(
      presence_absence_matrix
    )
  ) {
    stop(
      "Internal matrix conversion failed: community objects are not plain matrices.",
      call. = FALSE
    )
  }

  calculated_n <- rowSums(
    count_matrix
  )

  calculated_richness <- rowSums(
    presence_absence_matrix
  )

  expected_n <- sample_metadata$n_fish[
    match(
      rownames(
        count_matrix
      ),
      as.character(
        sample_metadata$Sample_ID
      )
    )
  ]

  expected_richness <- sample_metadata$species_richness[
    match(
      rownames(
        count_matrix
      ),
      as.character(
        sample_metadata$Sample_ID
      )
    )
  ]

  if (!isTRUE(
    all.equal(
      as.numeric(
        calculated_n
      ),
      as.numeric(
        expected_n
      ),
      tolerance = 0
    )
  )) {
    stop(
      "Count-matrix row totals did not match Sample_Metadata.",
      call. = FALSE
    )
  }

  if (!isTRUE(
    all.equal(
      as.numeric(
        calculated_richness
      ),
      as.numeric(
        expected_richness
      ),
      tolerance = 0
    )
  )) {
    stop(
      "Presence-absence richness did not match Sample_Metadata.",
      call. = FALSE
    )
  }

  if (max(
    abs(
      rowSums(
        relative_abundance_matrix
      ) -
      1
    )
  ) > 1e-12) {
    stop(
      "Relative-abundance rows did not sum to one.",
      call. = FALSE
    )
  }

  write.csv(
    data.frame(
      Sample_ID = rownames(
        count_matrix
      ),
      count_matrix,
      check.names = FALSE
    ),
    file.path(
      TABLE_DIR,
      "Table_S1_community_count_matrix.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    data.frame(
      Sample_ID = rownames(
        relative_abundance_matrix
      ),
      relative_abundance_matrix,
      check.names = FALSE
    ),
    file.path(
      TABLE_DIR,
      "Table_S2_relative_abundance_matrix.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    data.frame(
      Sample_ID = rownames(
        presence_absence_matrix
      ),
      presence_absence_matrix,
      check.names = FALSE
    ),
    file.path(
      TABLE_DIR,
      "Table_S3_presence_absence_matrix.csv"
    ),
    row.names = FALSE
  )

  sample_summary <- sample_metadata |>
    dplyr::mutate(
      Calculated_n_fish = as.numeric(
        calculated_n[
          as.character(
            Sample_ID
          )
        ]
      ),

      Calculated_species_richness = as.numeric(
        calculated_richness[
          as.character(
            Sample_ID
          )
        ]
      )
    )

  write.csv(
    sample_summary,
    file.path(
      TABLE_DIR,
      "Table_2_sample_metadata_and_richness.csv"
    ),
    row.names = FALSE
  )

  cat(
    "\nCommunity matrices created.",
    "\nSamples: ",
    nrow(
      count_matrix
    ),
    "\nSpecies: ",
    ncol(
      count_matrix
    ),
    "\nCount-matrix class: ",
    paste(
      class(
        count_matrix
      ),
      collapse = ", "
    ),
    "\nRelative-abundance-matrix class: ",
    paste(
      class(
        relative_abundance_matrix
      ),
      collapse = ", "
    ),
    "\nPresence-absence-matrix class: ",
    paste(
      class(
        presence_absence_matrix
      ),
      collapse = ", "
    ),
    "\n",
    sep = ""
  )


  # ---------------------------------------------------------------------------
  # 3.5 Calculate Bray-Curtis and Jaccard dissimilarities
  # ---------------------------------------------------------------------------

  # vegdist() now receives plain double matrices, not table objects.
  bray_distance <- vegan::vegdist(
    x = relative_abundance_matrix,
    method = "bray",
    binary = FALSE,
    na.rm = FALSE
  )

  jaccard_distance <- vegan::vegdist(
    x = presence_absence_matrix,
    method = "jaccard",
    binary = TRUE,
    na.rm = FALSE
  )

  if (
    !inherits(
      bray_distance,
      "dist"
    ) ||
    !inherits(
      jaccard_distance,
      "dist"
    )
  ) {
    stop(
      "Distance calculation failed to return valid 'dist' objects.",
      call. = FALSE
    )
  }

  cat(
    "\nDissimilarity matrices calculated successfully.",
    "\nBray-Curtis object class: ",
    paste(
      class(
        bray_distance
      ),
      collapse = ", "
    ),
    "\nJaccard object class: ",
    paste(
      class(
        jaccard_distance
      ),
      collapse = ", "
    ),
    "\n",
    sep = ""
  )

  distance_to_long <- function(
    distance_object,
    value_name
  ) {

    distance_matrix <- as.matrix(
      distance_object
    )

    sample_pairs <- utils::combn(
      rownames(
        distance_matrix
      ),
      2
    )

    output <- data.frame(
      Sample_1 = sample_pairs[1, ],
      Sample_2 = sample_pairs[2, ],
      stringsAsFactors = FALSE
    )

    output[[value_name]] <- vapply(
      seq_len(
        ncol(
          sample_pairs
        )
      ),
      function(pair_index) {

        distance_matrix[
          sample_pairs[1, pair_index],
          sample_pairs[2, pair_index]
        ]
      },
      FUN.VALUE = numeric(1)
    )

    output
  }

  bray_long <- distance_to_long(
    bray_distance,
    "Bray_Curtis"
  )

  jaccard_long <- distance_to_long(
    jaccard_distance,
    "Jaccard"
  )

  canonical_pair_key <- function(
    sample_1,
    sample_2
  ) {

    vapply(
      seq_along(
        sample_1
      ),
      function(index) {

        paste(
          sort(
            c(
              sample_1[index],
              sample_2[index]
            )
          ),
          collapse = "||"
        )
      },
      FUN.VALUE = character(1)
    )
  }

  bray_long$Pair_key <- canonical_pair_key(
    bray_long$Sample_1,
    bray_long$Sample_2
  )

  jaccard_long$Pair_key <- canonical_pair_key(
    jaccard_long$Sample_1,
    jaccard_long$Sample_2
  )

  pairwise_design$Pair_key <- canonical_pair_key(
    as.character(
      pairwise_design$Sample_1
    ),
    as.character(
      pairwise_design$Sample_2
    )
  )

  calculated_pairwise <- merge(
    bray_long[
      ,
      c(
        "Pair_key",
        "Bray_Curtis"
      )
    ],
    jaccard_long[
      ,
      c(
        "Pair_key",
        "Jaccard"
      )
    ],
    by = "Pair_key",
    all = TRUE,
    sort = FALSE
  )

  pairwise_results <- merge(
    pairwise_design,
    calculated_pairwise,
    by = "Pair_key",
    all.x = TRUE,
    sort = FALSE
  )

  pairwise_results <- pairwise_results[
    match(
      pairwise_design$Pair_ID,
      pairwise_results$Pair_ID
    ),
    ,
    drop = FALSE
  ]

  pairwise_results$Bray_validation_error <- (
    pairwise_results$Bray_Curtis -
    pairwise_results$Bray_Curtis_relative_abundance
  )

  pairwise_results$Jaccard_validation_error <- (
    pairwise_results$Jaccard -
    pairwise_results$Jaccard_presence_absence
  )

  maximum_bray_error <- max(
    abs(
      pairwise_results$Bray_validation_error
    ),
    na.rm = TRUE
  )

  maximum_jaccard_error <- max(
    abs(
      pairwise_results$Jaccard_validation_error
    ),
    na.rm = TRUE
  )

  if (
    maximum_bray_error > 1e-10 ||
    maximum_jaccard_error > 1e-10
  ) {
    stop(
      paste0(
        "R-calculated dissimilarities did not match the workbook validation values. ",
        "Maximum Bray-Curtis error = ",
        maximum_bray_error,
        "; maximum Jaccard error = ",
        maximum_jaccard_error,
        "."
      ),
      call. = FALSE
    )
  }

  write.csv(
    pairwise_results,
    file.path(
      TABLE_DIR,
      "Table_6_pairwise_dissimilarities.csv"
    ),
    row.names = FALSE
  )

  within_reach_results <- pairwise_results |>
    dplyr::filter(
      as.logical(
        Same_reach
      )
    ) |>
    dplyr::select(
      Pair_ID,
      Sample_1,
      Sample_2,
      Comparison_type,
      Bray_Curtis,
      Jaccard
    )

  write.csv(
    within_reach_results,
    file.path(
      TABLE_DIR,
      "Table_7_within_reach_wet_dry_dissimilarity.csv"
    ),
    row.names = FALSE
  )

  distance_validation <- data.frame(
    Metric = c(
      "Bray-Curtis",
      "Jaccard"
    ),
    Maximum_absolute_error = c(
      maximum_bray_error,
      maximum_jaccard_error
    ),
    Criterion = "< 1e-10",
    Validation_status = c(
      ifelse(
        maximum_bray_error < 1e-10,
        "Passed",
        "Failed"
      ),
      ifelse(
        maximum_jaccard_error < 1e-10,
        "Passed",
        "Failed"
      )
    )
  )

  write.csv(
    distance_validation,
    file.path(
      TABLE_DIR,
      "Table_S4_distance_validation.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.6 Fit two-dimensional NMDS ordinations
  # ---------------------------------------------------------------------------

  run_nmds <- function(
    community_matrix,
    distance_method,
    use_binary,
    seed_offset
  ) {

    set.seed(
      SEED + seed_offset
    )

    vegan::metaMDS(
      community_matrix,
      distance = distance_method,
      binary = use_binary,
      k = NMDS_DIMENSIONS,
      try = NMDS_TRY,
      trymax = NMDS_TRYMAX,
      autotransform = FALSE,
      noshare = FALSE,
      wascores = TRUE,
      trace = FALSE
    )
  }

  nmds_bray <- run_nmds(
    community_matrix = relative_abundance_matrix,
    distance_method = "bray",
    use_binary = FALSE,
    seed_offset = 101L
  )

  nmds_jaccard <- run_nmds(
    community_matrix = presence_absence_matrix,
    distance_method = "jaccard",
    use_binary = TRUE,
    seed_offset = 202L
  )

  saveRDS(
    nmds_bray,
    file.path(
      OBJECT_DIR,
      "NMDS_Bray_Curtis.rds"
    )
  )

  saveRDS(
    nmds_jaccard,
    file.path(
      OBJECT_DIR,
      "NMDS_Jaccard.rds"
    )
  )

  nmds_diagnostics <- data.frame(
    Metric = c(
      "Bray-Curtis relative abundance",
      "Jaccard presence-absence"
    ),

    Dimensions = NMDS_DIMENSIONS,

    Stress = c(
      nmds_bray$stress,
      nmds_jaccard$stress
    ),

    Converged = c(
      isTRUE(
        nmds_bray$converged
      ),
      isTRUE(
        nmds_jaccard$converged
      )
    ),

    Number_of_tries = c(
      nmds_bray$tries,
      nmds_jaccard$tries
    ),

    Stress_interpretation = c(
      ifelse(
        nmds_bray$stress < 0.10,
        "Good two-dimensional representation",
        ifelse(
          nmds_bray$stress < 0.20,
          "Usable with caution",
          "Poor two-dimensional representation"
        )
      ),
      ifelse(
        nmds_jaccard$stress < 0.10,
        "Good two-dimensional representation",
        ifelse(
          nmds_jaccard$stress < 0.20,
          "Usable with caution",
          "Poor two-dimensional representation"
        )
      )
    ),

    Design_note = paste(
      "Only six samples; ordination is descriptive."
    ),

    stringsAsFactors = FALSE
  )

  write.csv(
    nmds_diagnostics,
    file.path(
      TABLE_DIR,
      "Table_3_NMDS_diagnostics.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.7 Blocked exploratory PERMANOVA
  #
  # Campaign:
  #   permutations restricted within reach
  #   2^3 = 8 admissible arrangements including the observed arrangement
  #
  # Reach:
  #   permutations restricted within campaign
  #   (3!)^2 = 36 admissible arrangements including the observed arrangement
  # ---------------------------------------------------------------------------

  extract_adonis_result <- function(
    adonis_result,
    metric,
    tested_factor,
    blocking_factor,
    admissible_permutations
  ) {

    result_table <- as.data.frame(
      adonis_result
    )

    result_row <- result_table[
      1,
      ,
      drop = FALSE
    ]

    p_column <- grep(
      "^Pr",
      names(
        result_row
      ),
      value = TRUE
    )[1]

    data.frame(
      Distance_metric = metric,
      Tested_factor = tested_factor,
      Blocking_factor = blocking_factor,
      Df = as.numeric(
        result_row$Df
      ),
      Sum_of_squares = as.numeric(
        result_row$SumOfSqs
      ),
      R2 = as.numeric(
        result_row$R2
      ),
      Pseudo_F = as.numeric(
        result_row$F
      ),
      Permutation_p = as.numeric(
        result_row[[
          p_column
        ]]
      ),
      Permutations_requested = PERMUTATIONS_REQUESTED,
      Admissible_permutations_including_identity = admissible_permutations,
      Minimum_possible_p_approximately = 1 / admissible_permutations,
      Inference = "Exploratory: one aggregated sample per reach-campaign cell",
      stringsAsFactors = FALSE
    )
  }

  run_permanova_set <- function(
    distance_object,
    metric_label
  ) {

    set.seed(
      SEED + 303L
    )

    campaign_result <- vegan::adonis2(
      distance_object ~ Campaign,
      data = sample_metadata,
      permutations = PERMUTATIONS_REQUESTED,
      strata = sample_metadata$Reach_code,
      by = NULL,
      parallel = 1
    )

    set.seed(
      SEED + 404L
    )

    reach_result <- vegan::adonis2(
      distance_object ~ Reach_code,
      data = sample_metadata,
      permutations = PERMUTATIONS_REQUESTED,
      strata = sample_metadata$Campaign,
      by = NULL,
      parallel = 1
    )

    dplyr::bind_rows(
      extract_adonis_result(
        adonis_result = campaign_result,
        metric = metric_label,
        tested_factor = "Campaign",
        blocking_factor = "Reach_code",
        admissible_permutations = 8L
      ),

      extract_adonis_result(
        adonis_result = reach_result,
        metric = metric_label,
        tested_factor = "Reach_code",
        blocking_factor = "Campaign",
        admissible_permutations = 36L
      )
    )
  }

  permanova_results <- dplyr::bind_rows(
    run_permanova_set(
      distance_object = bray_distance,
      metric_label = "Bray-Curtis relative abundance"
    ),

    run_permanova_set(
      distance_object = jaccard_distance,
      metric_label = "Jaccard presence-absence"
    )
  )

  write.csv(
    permanova_results,
    file.path(
      TABLE_DIR,
      "Table_4_PERMANOVA_results.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.8 PERMDISP diagnostics
  # ---------------------------------------------------------------------------

  run_permdisp <- function(
    distance_object,
    group_factor,
    metric_label,
    factor_label,
    seed_offset
  ) {

    set.seed(
      SEED + seed_offset
    )

    dispersion_object <- vegan::betadisper(
      d = distance_object,
      group = group_factor,
      type = "median",
      bias.adjust = TRUE
    )

    permutation_result <- vegan::permutest(
      dispersion_object,
      permutations = PERMUTATIONS_REQUESTED,
      pairwise = FALSE
    )

    permutation_table <- as.data.frame(
      permutation_result$tab
    )

    p_column <- grep(
      "^Pr",
      names(
        permutation_table
      ),
      value = TRUE
    )[1]

    factor_row <- permutation_table[
      1,
      ,
      drop = FALSE
    ]

    group_distances <- tapply(
      dispersion_object$distances,
      group_factor,
      mean
    )

    test_row <- data.frame(
      Distance_metric = metric_label,
      Grouping_factor = factor_label,
      Df = as.numeric(
        factor_row$Df
      ),
      Sum_of_squares = as.numeric(
        factor_row$`Sum Sq`
      ),
      Mean_square = as.numeric(
        factor_row$`Mean Sq`
      ),
      F_value = as.numeric(
        factor_row$F
      ),
      Permutation_p = as.numeric(
        factor_row[[
          p_column
        ]]
      ),
      Permutations_requested = PERMUTATIONS_REQUESTED,
      Inference_note = paste(
        "Low power because only six community samples are available."
      ),
      stringsAsFactors = FALSE
    )

    group_rows <- data.frame(
      Distance_metric = metric_label,
      Grouping_factor = factor_label,
      Group_level = names(
        group_distances
      ),
      Mean_distance_to_group_median = as.numeric(
        group_distances
      ),
      stringsAsFactors = FALSE
    )

    list(
      test = test_row,
      groups = group_rows
    )
  }

  permdisp_objects <- list(
    run_permdisp(
      distance_object = bray_distance,
      group_factor = sample_metadata$Campaign,
      metric_label = "Bray-Curtis relative abundance",
      factor_label = "Campaign",
      seed_offset = 501L
    ),

    run_permdisp(
      distance_object = bray_distance,
      group_factor = sample_metadata$Reach_code,
      metric_label = "Bray-Curtis relative abundance",
      factor_label = "Reach_code",
      seed_offset = 502L
    ),

    run_permdisp(
      distance_object = jaccard_distance,
      group_factor = sample_metadata$Campaign,
      metric_label = "Jaccard presence-absence",
      factor_label = "Campaign",
      seed_offset = 503L
    ),

    run_permdisp(
      distance_object = jaccard_distance,
      group_factor = sample_metadata$Reach_code,
      metric_label = "Jaccard presence-absence",
      factor_label = "Reach_code",
      seed_offset = 504L
    )
  )

  permdisp_results <- dplyr::bind_rows(
    lapply(
      permdisp_objects,
      function(object) {
        object$test
      }
    )
  )

  permdisp_group_distances <- dplyr::bind_rows(
    lapply(
      permdisp_objects,
      function(object) {
        object$groups
      }
    )
  )

  write.csv(
    permdisp_results,
    file.path(
      TABLE_DIR,
      "Table_5_PERMDISP_results.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    permdisp_group_distances,
    file.path(
      TABLE_DIR,
      "Table_S5_group_dispersion_distances.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.9 Design and inference table
  # ---------------------------------------------------------------------------

  design_limitations <- data.frame(
    Analysis_component = c(
      "Community sampling unit",
      "Reach x Campaign interaction",
      "Campaign PERMANOVA",
      "Reach PERMANOVA",
      "NMDS",
      "Abundance interpretation"
    ),

    Implementation = c(
      "One aggregated species-composition sample per reach-campaign group",
      "Not tested",
      "Campaign with permutations blocked within reach",
      "Reach_code with permutations blocked within campaign",
      "Two-dimensional descriptive ordination",
      "Bray-Curtis uses relative abundance of analysed fish"
    ),

    Limitation_or_interpretation = c(
      "Six total community samples",
      "No within-cell replication; interaction is not estimable",
      "Only 8 admissible permutations; coarse p-value resolution",
      "Only 36 admissible permutations; coarse p-value resolution",
      "Stress and sample size must be reported",
      "Represents true abundance only if sampling effort and catchability were comparable"
    ),

    stringsAsFactors = FALSE
  )

  write.csv(
    design_limitations,
    file.path(
      TABLE_DIR,
      "Table_1_design_and_inference_limitations.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.10 Prepare ordination scores
  # ---------------------------------------------------------------------------

  extract_site_scores <- function(
    nmds_object,
    metric_label
  ) {

    scores_matrix <- vegan::scores(
      nmds_object,
      display = "sites",
      choices = c(
        1,
        2
      )
    )

    scores_data <- data.frame(
      Sample_ID = rownames(
        scores_matrix
      ),
      NMDS1 = scores_matrix[
        ,
        1
      ],
      NMDS2 = scores_matrix[
        ,
        2
      ],
      stringsAsFactors = FALSE
    )

    scores_data <- merge(
      scores_data,
      sample_metadata[
        ,
        c(
          "Sample_ID",
          "Reach_code",
          "Reach_name",
          "Campaign",
          "Campaign_label"
        ),
        drop = FALSE
      ],
      by = "Sample_ID",
      all.x = TRUE,
      sort = FALSE
    )

    scores_data$Distance_metric <- metric_label

    scores_data[
      match(
        expected_sample_ids,
        scores_data$Sample_ID
      ),
      ,
      drop = FALSE
    ]
  }

  nmds_bray_scores <- extract_site_scores(
    nmds_bray,
    "Bray-Curtis relative abundance"
  )

  nmds_jaccard_scores <- extract_site_scores(
    nmds_jaccard,
    "Jaccard presence-absence"
  )

  nmds_scores <- dplyr::bind_rows(
    nmds_bray_scores,
    nmds_jaccard_scores
  )

  write.csv(
    nmds_scores,
    file.path(
      TABLE_DIR,
      "Table_S6_NMDS_site_scores.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 3.11 Publication figures
  # ---------------------------------------------------------------------------

  campaign_colours <- c(
    "Wet campaign" = "#0072B2",
    "Dry campaign" = "#D55E00"
  )

  reach_shapes <- c(
    "Hanchuan" = 21,
    "Qianjiang" = 22,
    "Shayang" = 24
  )

  short_sample_labels <- c(
    "HC__Wet_2023" = "HC-Wet",
    "HC__Dry_2023" = "HC-Dry",
    "QJ__Wet_2023" = "QJ-Wet",
    "QJ__Dry_2023" = "QJ-Dry",
    "SY__Wet_2023" = "SY-Wet",
    "SY__Dry_2023" = "SY-Dry"
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

        legend.title = ggplot2::element_text(
          face = "bold"
        ),

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


  create_nmds_plot <- function(
    scores_data,
    stress_value,
    filename,
    distance_label
  ) {

    scores_data$Reach_name <- factor(
      scores_data$Reach_name,
      levels = c(
        "Hanchuan",
        "Qianjiang",
        "Shayang"
      )
    )

    scores_data$Campaign_label <- factor(
      scores_data$Campaign_label,
      levels = c(
        "Wet campaign",
        "Dry campaign"
      )
    )

    scores_data$Short_label <- short_sample_labels[
      scores_data$Sample_ID
    ]

    arrow_data <- merge(
      scores_data[
        scores_data$Campaign == "Wet_2023",
        c(
          "Reach_code",
          "Reach_name",
          "NMDS1",
          "NMDS2"
        )
      ],
      scores_data[
        scores_data$Campaign == "Dry_2023",
        c(
          "Reach_code",
          "NMDS1",
          "NMDS2"
        )
      ],
      by = "Reach_code",
      suffixes = c(
        "_wet",
        "_dry"
      )
    )

    ggplot2::ggplot(
      scores_data,
      ggplot2::aes(
        x = NMDS1,
        y = NMDS2
      )
    ) +
      ggplot2::geom_segment(
        data = arrow_data,
        mapping = ggplot2::aes(
          x = NMDS1_wet,
          y = NMDS2_wet,
          xend = NMDS1_dry,
          yend = NMDS2_dry,
          linetype = Reach_name
        ),
        inherit.aes = FALSE,
        linewidth = 0.7,
        colour = "grey35",
        arrow = grid::arrow(
          length = grid::unit(
            0.17,
            "cm"
          ),
          type = "closed"
        )
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          fill = Campaign_label,
          shape = Reach_name
        ),
        size = 4,
        stroke = 0.7,
        colour = "black"
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          label = Short_label
        ),
        nudge_y = 0.035,
        size = 3.1,
        check_overlap = TRUE
      ) +
      ggplot2::scale_fill_manual(
        values = campaign_colours
      ) +
      ggplot2::scale_shape_manual(
        values = reach_shapes
      ) +
      ggplot2::coord_equal() +
      ggplot2::labs(
        x = "NMDS1",
        y = "NMDS2",
        fill = "Campaign",
        shape = "Reach",
        linetype = "Reach",
        caption = paste0(
          distance_label,
          "; two-dimensional stress = ",
          format(
            stress_value,
            digits = 3
          ),
          ". Arrows indicate wet-to-dry displacement within each reach."
        )
      ) +
      theme_publication(
        base_size = 10
      )
  }


  figure_1 <- create_nmds_plot(
    scores_data = nmds_bray_scores,
    stress_value = nmds_bray$stress,
    filename = "Figure_1_Bray_Curtis_NMDS",
    distance_label = "Bray-Curtis dissimilarity based on species relative abundance"
  )

  save_figure(
    figure_1,
    "Figure_1_Bray_Curtis_NMDS",
    width = 6.9,
    height = 5.4
  )


  figure_2 <- create_nmds_plot(
    scores_data = nmds_jaccard_scores,
    stress_value = nmds_jaccard$stress,
    filename = "Figure_2_Jaccard_NMDS",
    distance_label = "Jaccard dissimilarity based on species presence-absence"
  )

  save_figure(
    figure_2,
    "Figure_2_Jaccard_NMDS",
    width = 6.9,
    height = 5.4
  )


  # Heatmap helper
  distance_matrix_to_long <- function(
    distance_object,
    metric_label
  ) {

    matrix_value <- as.matrix(
      distance_object
    )

    output <- expand.grid(
      Row_sample = rownames(
        matrix_value
      ),
      Column_sample = colnames(
        matrix_value
      ),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )

    output$Dissimilarity <- mapply(
      function(
        row_sample,
        column_sample
      ) {

        matrix_value[
          row_sample,
          column_sample
        ]
      },
      output$Row_sample,
      output$Column_sample
    )

    output$Row_label <- short_sample_labels[
      output$Row_sample
    ]

    output$Column_label <- short_sample_labels[
      output$Column_sample
    ]

    output$Metric <- metric_label

    output
  }


  create_heatmap <- function(
    heatmap_data,
    filename,
    caption_text
  ) {

    heatmap_data$Row_label <- factor(
      heatmap_data$Row_label,
      levels = rev(
        short_sample_labels[
          expected_sample_ids
        ]
      )
    )

    heatmap_data$Column_label <- factor(
      heatmap_data$Column_label,
      levels = short_sample_labels[
        expected_sample_ids
      ]
    )

    heatmap_plot <- ggplot2::ggplot(
      heatmap_data,
      ggplot2::aes(
        x = Column_label,
        y = Row_label,
        fill = Dissimilarity
      )
    ) +
      ggplot2::geom_tile(
        colour = "white",
        linewidth = 0.5
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          label = sprintf(
            "%.2f",
            Dissimilarity
          )
        ),
        size = 3
      ) +
      ggplot2::scale_fill_gradient(
        low = "white",
        high = "#2166AC",
        limits = c(
          0,
          1
        ),
        name = "Dissimilarity"
      ) +
      ggplot2::coord_fixed() +
      ggplot2::labs(
        x = NULL,
        y = NULL,
        caption = caption_text
      ) +
      theme_publication(
        base_size = 9.5
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 40,
          hjust = 1
        ),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank()
      )

    save_figure(
      heatmap_plot,
      filename,
      width = 6.6,
      height = 5.8
    )
  }


  bray_heatmap_data <- distance_matrix_to_long(
    bray_distance,
    "Bray-Curtis"
  )

  create_heatmap(
    bray_heatmap_data,
    "Figure_3_Bray_Curtis_dissimilarity_heatmap",
    "Bray-Curtis dissimilarity based on species relative abundance."
  )


  jaccard_heatmap_data <- distance_matrix_to_long(
    jaccard_distance,
    "Jaccard"
  )

  create_heatmap(
    jaccard_heatmap_data,
    "Figure_4_Jaccard_dissimilarity_heatmap",
    "Jaccard dissimilarity based on species presence-absence."
  )


  # Within-reach wet-dry dissimilarity figure
  within_reach_plot_data <- dplyr::bind_rows(
    data.frame(
      Pair_ID = within_reach_results$Pair_ID,
      Sample_1 = within_reach_results$Sample_1,
      Sample_2 = within_reach_results$Sample_2,
      Metric = "Bray-Curtis",
      Dissimilarity = within_reach_results$Bray_Curtis,
      stringsAsFactors = FALSE
    ),

    data.frame(
      Pair_ID = within_reach_results$Pair_ID,
      Sample_1 = within_reach_results$Sample_1,
      Sample_2 = within_reach_results$Sample_2,
      Metric = "Jaccard",
      Dissimilarity = within_reach_results$Jaccard,
      stringsAsFactors = FALSE
    )
  )

  within_reach_plot_data$Reach_code <- sub(
    "__.*$",
    "",
    within_reach_plot_data$Sample_1
  )

  within_reach_plot_data$Reach_name <- factor(
    dplyr::recode(
      within_reach_plot_data$Reach_code,
      HC = "Hanchuan",
      QJ = "Qianjiang",
      SY = "Shayang"
    ),
    levels = c(
      "Hanchuan",
      "Qianjiang",
      "Shayang"
    )
  )

  figure_5 <- ggplot2::ggplot(
    within_reach_plot_data,
    ggplot2::aes(
      x = Reach_name,
      y = Dissimilarity,
      colour = Metric,
      shape = Metric,
      group = Metric
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.6
    ) +
    ggplot2::geom_point(
      size = 3
    ) +
    ggplot2::scale_y_continuous(
      limits = c(
        0,
        1
      )
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "Bray-Curtis" = "#0072B2",
        "Jaccard" = "#D55E00"
      )
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "Bray-Curtis" = 16,
        "Jaccard" = 17
      )
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Wet-dry community dissimilarity",
      colour = "Metric",
      shape = "Metric",
      caption = "Values are descriptive because each reach-campaign combination is represented by one aggregated sample."
    ) +
    theme_publication(
      base_size = 10
    )

  save_figure(
    figure_5,
    "Figure_5_within_reach_wet_dry_dissimilarity",
    width = 6.5,
    height = 4.5
  )


  # ---------------------------------------------------------------------------
  # 3.12 Optional combined Excel results workbook
  # ---------------------------------------------------------------------------

  if (
    requireNamespace(
      "writexl",
      quietly = TRUE
    )
  ) {

    writexl::write_xlsx(
      list(
        Design_Limitations = design_limitations,
        Sample_Summary = sample_summary,
        NMDS_Diagnostics = nmds_diagnostics,
        PERMANOVA = permanova_results,
        PERMDISP = permdisp_results,
        Pairwise_Dissimilarity = pairwise_results,
        Wet_Dry_Dissimilarity = within_reach_results,
        NMDS_Site_Scores = nmds_scores
      ),
      path = file.path(
        RESULTS_DIR,
        "Fish_community_composition_results.xlsx"
      )
    )
  }


  # ---------------------------------------------------------------------------
  # 3.13 Save software information and completion status
  # ---------------------------------------------------------------------------

  software_versions <- data.frame(
    Software = c(
      "R",
      "vegan",
      "permute",
      "readxl",
      "dplyr",
      "ggplot2"
    ),

    Version = c(
      R.version.string,
      as.character(
        utils::packageVersion(
          "vegan"
        )
      ),
      as.character(
        utils::packageVersion(
          "permute"
        )
      ),
      as.character(
        utils::packageVersion(
          "readxl"
        )
      ),
      as.character(
        utils::packageVersion(
          "dplyr"
        )
      ),
      as.character(
        utils::packageVersion(
          "ggplot2"
        )
      )
    ),

    stringsAsFactors = FALSE
  )

  write.csv(
    software_versions,
    file.path(
      TABLE_DIR,
      "Table_S7_software_versions.csv"
    ),
    row.names = FALSE
  )

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
      "Fish-community composition analysis completed successfully.",
      paste("Completed:", Sys.time()),
      paste("Input workbook:", INPUT_FILE),
      paste("Results folder:", RESULTS_DIR),
      paste("Bray-Curtis NMDS stress:", nmds_bray$stress),
      paste("Jaccard NMDS stress:", nmds_jaccard$stress),
      paste("Maximum Bray-Curtis validation error:", maximum_bray_error),
      paste("Maximum Jaccard validation error:", maximum_jaccard_error),
      "PERMANOVA inference remains exploratory because the design has six aggregated samples."
    ),
    completed_file
  )

  cat(
    "\n============================================================\n",
    "FISH-COMMUNITY COMPOSITION ANALYSIS COMPLETED\n",
    "============================================================\n",
    "\nThis is a non-Bayesian multivariate community analysis.\n",
    "\nResults folder:\n",
    RESULTS_DIR,
    "\n\nMain tables:\n",
    "Table_1_design_and_inference_limitations.csv\n",
    "Table_2_sample_metadata_and_richness.csv\n",
    "Table_3_NMDS_diagnostics.csv\n",
    "Table_4_PERMANOVA_results.csv\n",
    "Table_5_PERMDISP_results.csv\n",
    "Table_6_pairwise_dissimilarities.csv\n",
    "Table_7_within_reach_wet_dry_dissimilarity.csv\n",
    "\nMain figures:\n",
    "Figure_1_Bray_Curtis_NMDS\n",
    "Figure_2_Jaccard_NMDS\n",
    "Figure_3_Bray_Curtis_dissimilarity_heatmap\n",
    "Figure_4_Jaccard_dissimilarity_heatmap\n",
    "Figure_5_within_reach_wet_dry_dissimilarity\n",
    "\nImportant: Reach x Campaign interaction was not tested because it is not estimable.\n",
    "============================================================\n",
    sep = ""
  )

  invisible(
    list(
      count_matrix = count_matrix,
      relative_abundance_matrix = relative_abundance_matrix,
      presence_absence_matrix = presence_absence_matrix,
      bray_distance = bray_distance,
      jaccard_distance = jaccard_distance,
      nmds_bray = nmds_bray,
      nmds_jaccard = nmds_jaccard,
      permanova_results = permanova_results,
      permdisp_results = permdisp_results,
      pairwise_results = pairwise_results,
      results_directory = RESULTS_DIR
    )
  )
}


# =============================================================================
# 4. RUN THE ANALYSIS AND SAVE A CLEAR ERROR LOG
# =============================================================================

analysis_result <- tryCatch(
  run_fish_community_analysis(),

  error = function(e) {

    error_message <- paste0(
      "Fish-community composition analysis stopped with an error.\n",
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
      "FISH_COMMUNITY_ANALYSIS_ERROR.txt"
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

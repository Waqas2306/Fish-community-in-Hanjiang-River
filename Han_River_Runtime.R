# Shared paths, source validation and refreshed SIBER dependencies.
# All generated workbooks are written under Results; supplied inputs are retained.

hr_files <- c(
  master = "Han_River 21.11.2026.xlsx",
  bayesian = "Bayesian Multivariate analysis.xlsx",
  community = "Han_River_Fish_Community_Composition_R_Dataset(1).xlsx",
  siber = "Han_River_SIBER_Bayesian_Isotopic_Niche_R_Dataset(1).xlsx",
  shared = "Han_River_Shared_Species_Wet_Dry_R_Dataset(1).xlsx",
  tp = "Han_River_Trophic_Position_R_Dataset(1).xlsx",
  decomposition = "Han_River_Within_Species_Turnover_Decomposition_R_Dataset(1).xlsx",
  morphometric = "Morphometric_Sensitivity_Analysis_Data(2).xlsx",
  isotope = "Isotope_Reorganization_Analysis(2).xlsx",
  linkage = "Han_River_Community_Turnover_Trophic_Reorganization_R_Dataset_Reviewer3.xlsx",
  figure8 = "Han_River_Figure8_Analysis_Ready_Reviewer3.xlsx",
  community_csv = "Reviewer4_Community_Counts (1).csv",
  cpue = "Reviewer4_CPUE_Event_Records (1).csv"
)

hr_root <- function() {
  p <- getOption("hanriver.project", Sys.getenv("HAN_RIVER_PROJECT", unset = ""))
  if (!nzchar(p)) p <- getwd()
  normalizePath(p, winslash = "/", mustWork = TRUE)
}

hr_raw_input <- function(module) {
  if (!module %in% names(hr_files)) stop("Unknown analysis module: ", module)
  p <- file.path(hr_root(), unname(hr_files[[module]]))
  if (!file.exists(p)) stop("Required input is missing: ", p,
                           "\nExtract the complete package before running START_HERE.R.")
  normalizePath(p, winslash = "/", mustWork = TRUE)
}

# The attached Figure 8 workbook may retain the download suffix (1).
if (!file.exists(file.path(hr_root(), hr_files[["figure8"]])) &&
    file.exists(file.path(hr_root(), "Han_River_Figure8_Analysis_Ready_Reviewer3(1).xlsx"))) {
  hr_files[["figure8"]] <- "Han_River_Figure8_Analysis_Ready_Reviewer3(1).xlsx"
}

hr_signature <- function() {
  paths <- unique(c(file.path(hr_root(), unname(hr_files)),
                    list.files(hr_root(), pattern = "[.]R$", full.names = TRUE)))
  paths <- sort(paths[file.exists(paths)])
  h <- tools::md5sum(paths)
  versions <- vapply(c("brms", "rstan", "SIBER", "rjags", "readxl"), function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "not-installed"
  }, character(1))
  signature_text <- c(paste(basename(paths), unname(h), sep = "="),
                      paste(names(versions), versions, sep = "="),
                      paste0("R=", getRversion()))
  f <- tempfile(fileext = ".txt")
  on.exit(unlink(f), add = TRUE)
  writeLines(signature_text, f, useBytes = TRUE)
  substr(unname(tools::md5sum(f)), 1L, 16L)
}

hr_results <- function(module) {
  p <- file.path(hr_root(), "Results", hr_signature(), module)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(p)) stop("Cannot create results folder: ", p)
  p
}

hr_read <- function(path, sheet) {
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Run SETUP_HAN_RIVER.R first: readxl is missing.")
  as.data.frame(hr_read_excel(path, sheet = sheet), stringsAsFactors = FALSE)
}

hr_master <- function() hr_read(hr_raw_input("master"), "Fish_Isotopes")

hr_compare_fish <- function(x, master = hr_master(), require_all = TRUE) {
  id <- if ("Fish_id" %in% names(x)) "Fish_id" else "Fish_isotope_record_id"
  if (!id %in% names(x)) stop("Fish identifier column is missing.")
  if (anyNA(x[[id]]) || anyDuplicated(x[[id]])) stop("Missing or duplicated fish IDs in input.")
  ix <- match(as.character(x[[id]]), as.character(master$Fish_isotope_record_id))
  if (anyNA(ix)) stop("Input contains fish IDs absent from the supplied master.")
  if (require_all && (!setequal(x[[id]], master$Fish_isotope_record_id))) stop("Input fish population differs from master.")
  map <- c(Reach_code = "Reach_code", Campaign = "Campaign", Species_code = "Species_code",
           Scientific_name = "Scientific_name", Feeding_guild = "Feeding_guild",
           d13C_permil = "d13C_permil", d15N_permil = "d15N_permil")
  for (nm in intersect(names(map), names(x))) {
    a <- x[[nm]]; b <- master[[map[[nm]]]][ix]
    ok <- if (is.numeric(b)) is.finite(as.numeric(a)) & abs(as.numeric(a) - b) < 1e-9 else as.character(a) == as.character(b)
    if (anyNA(ok) || !all(ok)) stop("Input does not match master: ", nm)
  }
  invisible(TRUE)
}

hr_validate_workbook <- function(path, module) {
  sheets <- c(bayesian = "Model_Data", community = "Individual_Data", siber = "Individual_Data",
              shared = "Shared_Species_Long", tp = "Fish_Data", decomposition = "Individual_Data",
              morphometric = "Fish_Isotopes", isotope = "Fish_Isotopes", figure8 = "Fish_Isotopes")
  if (module %in% names(sheets)) hr_compare_fish(hr_read(path, sheets[[module]]), require_all = module != "shared")
  invisible(path)
}

hr_morph_counts <- function() {
  iso <- hr_master()
  morph <- hr_read(hr_raw_input("master"), "Fish_Morphometrics")
  morph <- morph[morph$Principal_isotope_analysis == "Yes" & morph$Reach_code %in% c("HC", "QJ", "SY"), , drop = FALSE]
  key <- function(d) paste(d$Reach_code, d$Scientific_name, sep = "|")
  ki <- key(iso); km <- key(morph)
  shared <- intersect(unique(ki[iso$Campaign == "Wet_2023"]), unique(ki[iso$Campaign == "Dry_2023"]))
  complete <- matched <- logical(length(shared))
  for (i in seq_along(shared)) {
    checks <- lapply(c("Wet_2023", "Dry_2023"), function(camp) {
      mm <- morph[km == shared[i] & morph$Campaign == camp, , drop = FALSE]
      c(complete = any(is.finite(mm$Total_length_mm)) && any(is.finite(mm$Body_weight_g)),
        matched = sum(ki == shared[i] & iso$Campaign == camp) == nrow(mm))
    })
    complete[i] <- all(vapply(checks, function(z) z[["complete"]], logical(1)))
    matched[i] <- complete[i] && all(vapply(checks, function(z) z[["matched"]], logical(1)))
  }
  c(length(shared), sum(complete), sum(matched))
}

hr_validate_project <- function() {
  invisible(lapply(names(hr_files), hr_raw_input))
  master <- hr_master()
  if (anyDuplicated(master$Fish_isotope_record_id) || anyNA(master[c("d13C_permil", "d15N_permil")])) stop("Master isotope records are invalid.")
  modules <- c("bayesian", "community", "siber", "shared", "tp", "decomposition", "morphometric", "isotope", "figure8")
  for (module in modules) hr_validate_workbook(hr_raw_input(module), module)
  sib <- hr_read(hr_raw_input("siber"), "SIBER_Input")
  sib_ids <- hr_read(hr_raw_input("siber"), "Individual_Data")
  if (nrow(sib) != nrow(sib_ids) ||
      any(abs(as.numeric(sib$iso1) - as.numeric(sib_ids$d13C_permil)) > 1e-9) ||
      any(abs(as.numeric(sib$iso2) - as.numeric(sib_ids$d15N_permil)) > 1e-9) ||
      !identical(as.character(sib$group), as.character(sib_ids$group)) ||
      !identical(as.character(sib$community), as.character(sib_ids$community))) {
    stop("SIBER_Input is not synchronized with its ID-bearing Individual_Data sheet.")
  }
  tp <- hr_read(hr_raw_input("tp"), "Baseline_Inputs")
  primary <- toupper(as.character(tp$Primary_scenario)) %in% c("TRUE", "1", "YES")
  tp <- tp[primary, , drop = FALSE]
  master_base <- hr_read(hr_raw_input("master"), "TP_Baseline_SiteSeason")
  ix <- match(paste(tp$Reach_code, tp$Campaign, sep = "|"), master_base$Fish_group_key)
  if (nrow(tp) != 6L || anyNA(ix) ||
      any(tp$Scenario_ID != "External_PrimaryConsumer") ||
      any(tp$n_baseline_species != 1) ||
      any(abs(tp$Baseline_mean_d15N - master_base$Mean_d15N_permil[ix]) > 1e-9) ||
      any(abs(tp$Baseline_SD_d15N - master_base$SE_d15N_permil[ix]) > 1e-9)) {
    stop("Primary measured-snail TP baselines do not match the master.")
  }
  counts <- utils::read.csv(hr_raw_input("community_csv"), check.names = FALSE)
  for (i in seq_len(nrow(counts))) {
    k <- strsplit(as.character(counts$Sample_ID[i]), "__", fixed = TRUE)[[1L]]
    z <- table(master$Species_code[master$Reach_code == k[1L] & master$Campaign == k[2L]])
    for (sp in setdiff(names(counts), "Sample_ID")) {
      expected <- if (sp %in% names(z)) as.numeric(z[[sp]]) else 0
      if (counts[[sp]][i] != expected) stop("Community CSV does not match master: ", counts$Sample_ID[i], "/", sp)
    }
  }
  sizes <- hr_morph_counts()
  result <- data.frame(Check = c("Master isotope records", "Shared reach-species pairs", "Morphometric-complete pairs", "Count-matched pairs"),
                       Value = c(nrow(master), sizes), stringsAsFactors = FALSE)
  print(result, row.names = FALSE)
  utils::write.csv(result, file.path(hr_root(), "PREFLIGHT_DATA_CHECK.csv"), row.names = FALSE)
  invisible(result)
}

hr_source_module <- function(filename) {
  p <- file.path(hr_root(), filename)
  if (!file.exists(p)) stop("Analysis script is missing: ", filename)
  old <- options()
  on.exit(options(old), add = TRUE)
  source(p, local = new.env(parent = globalenv()), encoding = "UTF-8")
  invisible(TRUE)
}

hr_ensure_main_model <- function() {
  p <- hr_results("main_model")
  required <- file.path(p, c("Han_River_Bayesian_MV_final.rds", "model_data_used.csv"))
  if (!all(file.exists(required))) {
    message("The current project has no fitted multivariate model. Starting the main-model script.")
    hr_source_module("Han_River_Bayesian_MV_FINAL_CLEAN_v2.R")
  }
  if (!all(file.exists(required))) stop("The main model did not produce its required files.")
  used <- utils::read.csv(required[2L], stringsAsFactors = FALSE)
  hr_compare_fish(used)
  invisible(p)
}

hr_refresh_dependencies <- function() {
  base <- hr_results("siber")
  tabdir <- file.path(base, "Tables")
  filenames <- c(sea = "Table_3_Bayesian_SEA_summary.csv", seapair = "Table_4_pairwise_SEA_comparisons.csv",
                 cent = "Table_5_posterior_centroid_locations.csv", centpair = "Table_6_pairwise_centroid_differences.csv",
                 overlap = "Table_8_Bayesian_95pct_ellipse_overlap.csv")
  paths <- file.path(tabdir, filenames)
  if (!all(file.exists(paths))) {
    message("Fresh SIBER outputs are required. Starting SIBER; this may take some time.")
    hr_source_module("Han_River_SIBER_Bayesian_Isotopic_Niche_Analysis(1).R")
  }
  if (!all(file.exists(paths))) stop("SIBER did not produce all posterior summary tables.")
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Run SETUP_HAN_RIVER.R: openxlsx is missing.")
  tables <- lapply(paths, utils::read.csv, stringsAsFactors = FALSE, check.names = FALSE)
  names(tables) <- names(filenames)
  assert_columns <- function(d, cols, label) {
    missing <- setdiff(cols, names(d))
    if (length(missing)) stop(label, " is missing: ", paste(missing, collapse = ", "))
  }
  assert_columns(tables$sea, c("SIBER_code", "Posterior_median_SEA", "Lower_95", "Upper_95"), "SIBER SEA summary")
  assert_columns(tables$cent, c("SIBER_code", "Posterior_median_d13C", "Posterior_median_d15N"), "SIBER centroids")
  assert_columns(tables$centpair, c("Pair_ID", "Median_centroid_distance", "Median_delta_d13C_2_minus_1", "Median_delta_d15N_2_minus_1"), "SIBER centroid pairs")
  assert_columns(tables$seapair, c("Pair_ID", "Median_SEA_difference_2_minus_1", "Median_SEA_ratio_2_over_1"), "SIBER SEA pairs")
  assert_columns(tables$overlap, c("Pair_ID", "Median_Jaccard_percent"), "SIBER overlaps")
  get_rows <- function(d, ids, field) {
    if (anyDuplicated(d[[field]])) stop("Duplicate posterior summary identifiers: ", field)
    ix <- match(as.character(ids), as.character(d[[field]]))
    if (anyNA(ix)) stop("Posterior summaries do not cover the requested groups/pairs.")
    d[ix, , drop = FALSE]
  }
  linkage_source <- hr_raw_input("linkage")
  meta <- hr_read(linkage_source, "Group_Metadata")
  pairs <- hr_read(linkage_source, "Pairwise_Linkage")
  seasonal <- hr_read(linkage_source, "Seasonal_Linkage")
  cc <- get_rows(tables$cent, meta$SIBER_code, "SIBER_code")
  ss <- get_rows(tables$sea, meta$SIBER_code, "SIBER_code")
  for (nm in c("Posterior_median_d13C", "Lower_95_d13C", "Upper_95_d13C", "Posterior_median_d15N", "Lower_95_d15N", "Upper_95_d15N")) meta[[nm]] <- cc[[nm]]
  meta$Posterior_median_SEA <- ss$Posterior_median_SEA
  meta$Lower_95_SEA <- ss$Lower_95; meta$Upper_95_SEA <- ss$Upper_95
  cp <- get_rows(tables$centpair, pairs$Pair_ID, "Pair_ID")
  sp <- get_rows(tables$seapair, pairs$Pair_ID, "Pair_ID")
  ov <- get_rows(tables$overlap, pairs$Pair_ID, "Pair_ID")
  pairs$Median_centroid_distance <- cp$Median_centroid_distance
  pairs$Lower_95_centroid_distance <- cp$Lower_95_distance
  pairs$Upper_95_centroid_distance <- cp$Upper_95_distance
  pairs$Absolute_delta_d13C <- abs(cp$Median_delta_d13C_2_minus_1)
  pairs$Absolute_delta_d15N <- abs(cp$Median_delta_d15N_2_minus_1)
  pairs$Median_SEA_difference_2_minus_1 <- sp$Median_SEA_difference_2_minus_1
  pairs$Absolute_SEA_difference <- abs(sp$Median_SEA_difference_2_minus_1)
  pairs$Median_SEA_ratio_2_over_1 <- sp$Median_SEA_ratio_2_over_1
  pairs$Absolute_log_SEA_ratio <- abs(log(sp$Median_SEA_ratio_2_over_1))
  pairs$Probability_SEA_2_greater_1 <- sp$Probability_SEA_2_greater_1
  pairs$Median_ellipse_Jaccard_overlap_percent <- ov$Median_Jaccard_percent
  pairs$Niche_nonoverlap_proportion <- 1 - ov$Median_Jaccard_percent / 100
  pairs$Lower_95_ellipse_Jaccard_overlap_percent <- ov$Lower_95_Jaccard_percent
  pairs$Upper_95_ellipse_Jaccard_overlap_percent <- ov$Upper_95_Jaccard_percent
  primary <- c("Median_centroid_distance", "Absolute_delta_d13C", "Absolute_delta_d15N", "Absolute_SEA_difference", "Absolute_log_SEA_ratio", "Niche_nonoverlap_proportion")
  if (nrow(meta) != 6L || nrow(pairs) != 15L || anyNA(pairs[primary]) || !all(vapply(pairs[primary], function(x) all(is.finite(x)), logical(1)))) stop("Refreshed posterior inputs are incomplete.")
  pp <- get_rows(pairs, seasonal$Pair_ID, "Pair_ID")
  smap <- c(Bray_Curtis_turnover = "Bray_Curtis", Jaccard_turnover = "Jaccard", Median_centroid_distance = "Median_centroid_distance",
            Lower_95_centroid_distance = "Lower_95_centroid_distance", Upper_95_centroid_distance = "Upper_95_centroid_distance",
            Absolute_delta_d13C = "Absolute_delta_d13C", Absolute_delta_d15N = "Absolute_delta_d15N",
            Median_SEA_change = "Median_SEA_difference_2_minus_1", Absolute_SEA_change = "Absolute_SEA_difference",
            Median_SEA_ratio = "Median_SEA_ratio_2_over_1", Absolute_log_SEA_ratio = "Absolute_log_SEA_ratio",
            Median_ellipse_Jaccard_overlap_percent = "Median_ellipse_Jaccard_overlap_percent", Niche_nonoverlap_proportion = "Niche_nonoverlap_proportion")
  for (nm in names(smap)) seasonal[[nm]] <- pp[[smap[[nm]]]]
  seasonal$Inference_status <- "Fresh SIBER summaries; descriptive seasonal comparison; three reaches"
  destdir <- file.path(hr_results("prepared_inputs"))
  stamp_refresh <- function(wb, source) {
    notes <- as.data.frame(hr_read_excel(source, sheet = "README", col_names = FALSE))
    labels <- ifelse(is.na(notes[[1L]]), "", as.character(notes[[1L]]))
    rows <- which(grepl("REFIT_REQUIRED|Partial result refresh required|Posterior refit required", labels, ignore.case = TRUE))
    for (row in rows) openxlsx::writeData(wb, "README",
      matrix(c("SIBER_REFRESHED", "Fresh SIBER posterior inputs generated for the current data/code signature. Dependent tests are recalculated by the analysis scripts; inspect their output tables and sampler diagnostics."), nrow = 1L),
      startRow = row, colNames = FALSE)
    validation <- as.data.frame(hr_read_excel(source, sheet = "Validation_Summary", col_names = FALSE))
    rows <- which(as.character(validation[[1L]]) == "Cached posterior-result status")
    for (row in rows) openxlsx::writeData(wb, "Validation_Summary", "SIBER_INPUTS_REFRESHED",
      startRow = row, startCol = 2L, colNames = FALSE)
  }
  linkage_dest <- file.path(destdir, hr_files[["linkage"]])
  wb <- hr_load_analysis_workbook(linkage_source)
  openxlsx::writeData(wb, "Group_Metadata", meta)
  openxlsx::writeData(wb, "Pairwise_Linkage", pairs)
  openxlsx::writeData(wb, "Seasonal_Linkage", seasonal)
  matrix_map <- c(Centroid_Distance_Matrix = "Median_centroid_distance", SEA_Difference_Matrix = "Absolute_SEA_difference",
                  Log_SEA_Ratio_Matrix = "Absolute_log_SEA_ratio", Niche_Nonoverlap_Matrix = "Niche_nonoverlap_proportion")
  for (sn in names(matrix_map)) {
    z <- matrix(0, 6L, 6L, dimnames = list(meta$Sample_ID, meta$Sample_ID))
    for (i in seq_len(nrow(pairs))) {
      a <- as.character(pairs$Sample_1[i]); b <- as.character(pairs$Sample_2[i]); z[a,b] <- z[b,a] <- pairs[[matrix_map[[sn]]]][i]
    }
    zz <- data.frame(Sample_ID = rownames(z), z, check.names = FALSE)
    openxlsx::writeData(wb, sn, zz)
  }
  stamp_refresh(wb, linkage_source)
  openxlsx::saveWorkbook(wb, linkage_dest, overwrite = TRUE)
  figure_source <- hr_raw_input("figure8")
  figure_dest <- file.path(destdir, hr_files[["figure8"]])
  fmeta <- hr_read(figure_source, "Biological_Group_Metadata")
  fpairs <- hr_read(figure_source, "Biological_Pairwise")
  fm <- get_rows(meta, fmeta$Sample_ID, "Sample_ID"); fp <- get_rows(pairs, fpairs$Pair_ID, "Pair_ID")
  for (nm in intersect(names(fmeta), names(fm))) fmeta[[nm]] <- fm[[nm]]
  for (nm in intersect(names(fpairs), names(fp))) fpairs[[nm]] <- fp[[nm]]
  wb <- hr_load_analysis_workbook(figure_source)
  openxlsx::writeData(wb, "Biological_Group_Metadata", fmeta)
  openxlsx::writeData(wb, "Biological_Pairwise", fpairs)
  stamp_refresh(wb, figure_source)
  openxlsx::saveWorkbook(wb, figure_dest, overwrite = TRUE)
  utils::write.csv(data.frame(File = basename(paths), MD5 = unname(tools::md5sum(paths))),
                   file.path(destdir, "SIBER_REFRESH_SOURCES.csv"), row.names = FALSE)
  invisible(c(linkage = linkage_dest, figure8 = figure_dest))
}

hr_input <- function(module) {
  if (module %in% c("linkage", "figure8")) return(unname(hr_refresh_dependencies()[[module]]))
  hr_raw_input(module)
}

hr_cpue_ready <- function(events) {
  truth <- function(x) !is.na(x) & toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES")
  needed <- c("Physical_deployment_identity_verified", "Complete_zero_catch_coverage_verified", "Date", "Reach_code", "Campaign", "Catch_n")
  reasons <- character()
  if (!all(needed %in% names(events))) return(list(ready = FALSE, reasons = "Verification fields or source dates are missing."))
  if (!all(truth(events$Physical_deployment_identity_verified))) reasons <- c(reasons, "Physical deployment identities are not all verified.")
  if (!all(truth(events$Complete_zero_catch_coverage_verified))) reasons <- c(reasons, "Complete zero-catch coverage is not verified.")
  dates <- suppressWarnings(as.Date(events$Date))
  wet <- events$Campaign == "Wet_2023"
  start <- as.Date(ifelse(wet, "2023-06-20", "2023-11-14")); end <- as.Date(ifelse(wet, "2023-06-27", "2023-11-23"))
  if (anyNA(dates) || any(dates < start | dates > end, na.rm = TRUE)) reasons <- c(reasons, "Event dates do not all match the stated campaign windows.")
  master <- hr_read(hr_raw_input("master"), "Fish_Morphometrics")
  for (reach in c("HC", "QJ", "SY")) for (camp in c("Wet_2023", "Dry_2023")) {
    n <- sum(master$Reach_code == reach & master$Campaign == camp)
    got <- sum(as.numeric(events$Catch_n[events$Reach_code == reach & events$Campaign == camp]))
    if (!is.finite(got) || got != n) reasons <- c(reasons, paste(reach, camp, "event catch differs from master."))
  }
  list(ready = !length(reasons), reasons = reasons)
}

# Figure typography: point sizes at the exported physical dimensions.
# These are readability choices for this manuscript, not journal-mandated sizes.
hr_figure_settings <- function() {
  defaults <- list(base = 14, axis_title = 16, axis_text = 14,
                   facet = 14, legend = 13, tag = 18, annotation_mm = 4.25,
                   family = "sans", raster_dpi = 600)
  utils::modifyList(defaults, getOption("hanriver.figure_style", list()))
}

hr_typography <- function() {
  z <- hr_figure_settings()
  ggplot2::theme(
    text = ggplot2::element_text(family = z$family, size = z$base, colour = "black"),
    axis.title = ggplot2::element_text(size = z$axis_title),
    axis.title.x = ggplot2::element_text(size = z$axis_title,
                                       margin = ggplot2::margin(t = 9)),
    axis.title.y = ggplot2::element_text(size = z$axis_title,
                                       margin = ggplot2::margin(r = 9)),
    axis.text = ggplot2::element_text(size = z$axis_text, colour = "black"),
    axis.text.x = ggplot2::element_text(size = z$axis_text),
    axis.text.y = ggplot2::element_text(size = z$axis_text),
    strip.text = ggplot2::element_text(size = z$facet, face = "bold"),
    strip.text.x = ggplot2::element_text(size = z$facet),
    strip.text.y = ggplot2::element_text(size = z$facet),
    legend.title = ggplot2::element_text(size = z$legend, face = "bold"),
    legend.text = ggplot2::element_text(size = z$legend),
    plot.title = ggplot2::element_text(size = 18, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = z$legend),
    plot.caption = ggplot2::element_text(size = z$legend),
    plot.tag = ggplot2::element_text(size = z$tag, face = "bold"),
    panel.spacing = grid::unit(12, "pt"),
    legend.key.size = grid::unit(14, "pt"),
    legend.spacing.x = grid::unit(5, "pt"),
    plot.margin = ggplot2::margin(12, 18, 12, 12)
  )
}

hr_style_plot <- function(plot) {
  # Adding theme properties preserves data, scales, rotations and statistical layers.
  if (inherits(plot, "patchwork")) return(plot & hr_typography())
  if (inherits(plot, "ggplot")) return(plot + hr_typography())
  plot
}

hr_ggsave <- function(filename, plot = ggplot2::last_plot(), device = NULL,
                      path = NULL, scale = 1, width = NA, height = NA,
                      units = "in", dpi = 600, limitsize = TRUE,
                      bg = "white", companions = FALSE, ...) {
  plot <- hr_style_plot(plot)
  ext <- tolower(tools::file_ext(filename))
  dpi <- if (is.numeric(dpi)) max(dpi, hr_figure_settings()$raster_dpi) else
    hr_figure_settings()$raster_dpi
  extra <- list(...)
  # Cairo embeds available fonts. The standard PDF fallback remains portable.
  if (ext == "pdf") {
    if (capabilities("cairo")) {
      device <- grDevices::cairo_pdf
      extra$useDingbats <- NULL
    } else {
      device <- grDevices::pdf
      extra$useDingbats <- FALSE
    }
  }
  args <- c(list(filename = filename, plot = plot, device = device, path = path,
                 scale = scale, width = width, height = height, units = units,
                 dpi = dpi, limitsize = limitsize, bg = bg), extra)
  do.call(ggplot2::ggsave, args)
  if (isTRUE(companions) && ext == "png") {
    stem <- sub("[.]png$", "", filename, ignore.case = TRUE)
    hr_ggsave(paste0(stem, ".tiff"), plot = plot, device = "tiff", path = path,
              scale = scale, width = width, height = height, units = units,
              dpi = dpi, limitsize = limitsize, bg = bg, compression = "lzw")
    hr_ggsave(paste0(stem, ".pdf"), plot = plot, path = path,
              scale = scale, width = width, height = height, units = units,
              dpi = dpi, limitsize = limitsize, bg = bg)
  }
  invisible(filename)
}

# Reporting labels only. Internal reach codes, fish IDs, dates and measurements
# remain unchanged so source validation, joins and model coefficients stay valid.
hr_site_text <- function(x) {
  map <- c(Hanchuan = "S1", Qianjiang = "S2", Shayang = "S3")
  for (old in names(map)) x <- gsub(old, unname(map[[old]]), x, ignore.case = TRUE)
  for (site in unname(map)) {
    x <- gsub(paste0(site, " (", site, ")"), site, x, fixed = TRUE)
  }
  x
}

hr_normalize_site_names <- function(x) {
  if (is.data.frame(x)) {
    for (j in seq_along(x)) {
      if (is.factor(x[[j]])) {
        levels(x[[j]]) <- hr_site_text(levels(x[[j]]))
      } else if (is.character(x[[j]])) {
        x[[j]] <- hr_site_text(x[[j]])
      }
    }
    return(x)
  }
  if (is.character(x)) return(hr_site_text(x))
  x
}

hr_read_excel <- function(...) {
  hr_normalize_site_names(readxl::read_excel(...))
}


# Compatibility reader for analysis workbooks exported with prefixed SpreadsheetML.
# Some openxlsx releases look for <sheets>/<sheet> and fail on <x:sheets>/<x:sheet>.
# Only temporary XML markup is normalized. Source files and cell values are untouched.
hr_load_analysis_workbook <- function(path) {
  expected_sheets <- readxl::excel_sheets(path)
  attempt <- tryCatch(openxlsx::loadWorkbook(path), error = function(e) e)
  if (!inherits(attempt, "error") && identical(names(attempt), expected_sheets)) {
    count <- attr(attempt$sharedStrings, "uniqueCount")
    if (length(count) != 1L || is.na(count)) {
      attr(attempt$sharedStrings, "uniqueCount") <- length(attempt$sharedStrings)
    }
    return(attempt)
  }
  original_message <- if (inherits(attempt, "error")) conditionMessage(attempt) else
    "The workbook loader did not preserve all worksheet names."
  message("Normalizing a temporary workbook copy for openxlsx compatibility: ", basename(path))
  folder <- tempfile("hanriver_workbook_xml_")
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  utils::unzip(path, exdir = folder)
  xml_files <- list.files(folder, pattern = "[.]xml$", recursive = TRUE, full.names = TRUE)
  main_namespace <- "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
  binding_pattern <- paste0('xmlns:[A-Za-z_][A-Za-z0-9_.-]*="', main_namespace, '"')
  changed <- 0L
  for (file in xml_files) {
    text <- paste(readLines(file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    bindings <- unique(regmatches(text, gregexpr(binding_pattern, text, perl = TRUE))[[1L]])
    if (!length(bindings)) next
    for (binding in bindings) {
      prefix <- sub('^xmlns:([^=]+)=.*$', '\\1', binding)
      text <- gsub(paste0("<", prefix, ":"), "<", text, fixed = TRUE)
      text <- gsub(paste0("</", prefix, ":"), "</", text, fixed = TRUE)
      default_binding <- paste0('xmlns="', main_namespace, '"')
      text <- gsub(binding, if (grepl(default_binding, text, fixed = TRUE)) "" else default_binding,
                   text, fixed = TRUE)
    }
    writeLines(text, file, useBytes = TRUE)
    changed <- changed + 1L
  }
  if (!changed) stop("Workbook compatibility repair could not normalize the XML. Original error: ",
                    original_message, call. = FALSE)
  result <- tryCatch(openxlsx::loadWorkbook(folder, isUnzipped = TRUE), error = function(e) {
    stop("Workbook compatibility repair failed for ", basename(path), ": ", conditionMessage(e),
         "\nOriginal error: ", original_message, call. = FALSE)
  })
  if (!identical(names(result), expected_sheets)) stop("Workbook repair changed worksheet names.")
  # uniqueCount is optional in SpreadsheetML but required by openxlsx's writer.
  count <- attr(result$sharedStrings, "uniqueCount")
  if (length(count) != 1L || is.na(count)) {
    attr(result$sharedStrings, "uniqueCount") <- length(result$sharedStrings)
  }
  result
}

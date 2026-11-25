# Plot-only refresh: run this file from the existing complete analysis folder.
# Reuses completed TP/Mantel tables; does not rerun SIBER or fit any model.
# Original Excel, worksheet and figure filenames are retained.
TABLE_DIR <- ""  # Optional: exact Tables folder from your completed Figure 8 run.
WORKBOOK <- ""   # Optional: exact Figure 8 Excel path.

# Figure layout only: no model fitting or changes to numerical results.
f8_colours <- c(HC = "#0072B2", QJ = "#D55E00", SY = "#009E73")
f8_shapes <- c(HC = 21, QJ = 22, SY = 24)
f8_sites <- c(HC = "S1", QJ = "S2", SY = "S3")
f8_variables <- c("Water_temperature_C", "DO_mg_L", "Flow_velocity_m_s",
                  "Water_TN_mg_L", "Water_TP_mg_L", "Chlorophyll_a_ug_L")
f8_variable_labels <- c("(a) Water temperature\n(°C)", "(b) Dissolved oxygen\n(mg L⁻¹)",
  "(c) Flow velocity\n(m s⁻¹)", "(d) Total nitrogen\n(mg L⁻¹)",
  "(e) Total phosphorus\n(mg L⁻¹)", "(f) Chlorophyll-a\n(µg L⁻¹)")
f8_theme <- function() {
  ggplot2::theme_classic(base_size = 11, base_family = "sans") + ggplot2::theme(
    text = ggplot2::element_text(colour = "#202020"),
    axis.title = ggplot2::element_text(size = 12.5, face = "bold"),
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 7)),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 7)),
    axis.text = ggplot2::element_text(size = 11, colour = "#202020"),
    panel.border = ggplot2::element_rect(fill = NA, colour = "#383838", linewidth = 0.35),
    axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_line(linewidth = 0.3),
    panel.grid.major.y = ggplot2::element_line(colour = "#E5E5E5", linewidth = 0.25),
    panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "#F2F2F2", colour = "#555555", linewidth = 0.3),
    strip.text = ggplot2::element_text(size = 10.5, face = "bold", lineheight = 1.05,
                                     margin = ggplot2::margin(5, 3, 5, 3)),
    panel.spacing = grid::unit(8, "pt"), legend.position = "top",
    legend.title = ggplot2::element_text(size = 11, face = "bold"),
    legend.text = ggplot2::element_text(size = 11),
    legend.key.size = grid::unit(12, "pt"), legend.spacing.x = grid::unit(4, "pt"),
    plot.title = ggplot2::element_text(size = 13, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 11), plot.caption = ggplot2::element_text(size = 10),
    plot.tag = ggplot2::element_text(size = 14, face = "bold"),
    plot.margin = ggplot2::margin(7, 9, 7, 7))
}
f8_campaign <- function(x) {
  x <- as.character(x); x[x == "Wet_2023"] <- "Wet"; x[x == "Dry_2023"] <- "Dry"
  factor(x, levels = c("Wet", "Dry"))
}
f8_site_factor <- function(x) {
  x <- as.character(x); x[x == "S1"] <- "HC"; x[x == "S2"] <- "QJ"; x[x == "S3"] <- "SY"
  factor(x, levels = names(f8_sites))
}
f8_environment_figure <- function(summary, raw) {
  prepare <- function(d) {
    d <- d[d$Variable %in% f8_variables, , drop = FALSE]
    d$Campaign <- f8_campaign(d$Campaign); d$Reach_code <- f8_site_factor(d$Reach_code)
    d$Panel <- factor(d$Variable, levels = f8_variables, labels = f8_variable_labels); d
  }
  summary <- prepare(summary); raw <- prepare(raw)
  ggplot2::ggplot(summary, ggplot2::aes(Campaign, Mean, fill = Reach_code,
      shape = Reach_code, group = Reach_code)) +
    ggplot2::geom_point(data = raw, ggplot2::aes(y = Value), alpha = 0.28,
      size = 1.7, stroke = 0.3, colour = "black",
      position = ggplot2::position_jitterdodge(jitter.width = 0.055, dodge.width = 0.58, seed = 20260731L)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean - SE, ymax = Mean + SE, colour = Reach_code),
      position = ggplot2::position_dodge(width = 0.58), width = 0.10, linewidth = 0.5, show.legend = FALSE) +
    ggplot2::geom_point(size = 2.9, stroke = 0.55, colour = "black",
      position = ggplot2::position_dodge(width = 0.58)) +
    ggplot2::facet_wrap(~Panel, scales = "free_y", ncol = 3) +
    ggplot2::scale_fill_manual(values = f8_colours, breaks = names(f8_sites), labels = f8_sites, name = "Site") +
    ggplot2::scale_shape_manual(values = f8_shapes, breaks = names(f8_sites), labels = f8_sites, name = "Site") +
    ggplot2::scale_colour_manual(values = f8_colours, guide = "none") +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.4)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.07, 0.15))) +
    ggplot2::labs(x = "Campaign", y = "Measured value") + f8_theme() +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1,
      override.aes = list(shape = unname(f8_shapes), size = 3.0, colour = "black", alpha = 1)), shape = "none")
}
f8_tp_figure <- function(d) {
  d$Campaign <- f8_campaign(d$Campaign); d$Reach_code <- f8_site_factor(d$Reach_code)
  ggplot2::ggplot(d, ggplot2::aes(Campaign, Point_TP, colour = Reach_code,
      fill = Reach_code, shape = Reach_code, group = Reach_code)) +
    ggplot2::geom_hline(yintercept = 2, linetype = "dashed", colour = "#888888", linewidth = 0.35) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower_95_TP, ymax = Upper_95_TP), width = 0.08, linewidth = 0.5) +
    ggplot2::geom_point(size = 3.0, colour = "black", stroke = 0.55) +
    ggplot2::scale_colour_manual(values = f8_colours, breaks = names(f8_sites), labels = f8_sites, name = "Site") +
    ggplot2::scale_fill_manual(values = f8_colours, breaks = names(f8_sites), labels = f8_sites, name = "Site") +
    ggplot2::scale_shape_manual(values = f8_shapes, breaks = names(f8_sites), labels = f8_sites, name = "Site") +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.35)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.12))) +
    ggplot2::labs(x = "Campaign", y = "Fish trophic position\n(gastropod reference)") + f8_theme() +
    ggplot2::guides(colour = "none", fill = ggplot2::guide_legend(nrow = 1,
      override.aes = list(shape = unname(f8_shapes), colour = "black", size = 3.0)), shape = "none")
}
f8_mantel_figure <- function(d, show_p = FALSE) {
  biological <- c(Bray_Curtis = "Bray–Curtis\ncommunity", Jaccard = "Jaccard\ncommunity",
    Centroid_Distance = "Isotope-centroid\ndistance", Absolute_SEA_Difference = "|SEA difference|",
    Absolute_Log_SEA_Ratio = "|log SEA ratio|", Niche_Nonoverlap = "Niche non-overlap",
    Revised_TP_Difference = "Gastropod-referenced\nTP difference")
  environmental <- c(Overall_Environment = "Overall\nenvironment", Hydrology_Habitat = "Hydrology–\nhabitat",
    Thermal_Oxygen = "Thermal–\noxygen", Nutrient_Productivity = "Nutrient–\nproductivity")
  d$X <- factor(environmental[d$Environmental_metric], levels = unname(environmental))
  d$Y <- factor(biological[d$Biological_metric], levels = rev(unname(biological)))
  if (anyNA(d$X) || anyNA(d$Y)) stop("Unrecognized metric names in the Mantel results.")
  d$Label <- if (show_p) sprintf("%.2f\nP = %.3f", d$Mantel_r, d$Exact_one_sided_P) else sprintf("%.2f", d$Mantel_r)
  d$Text_colour <- ifelse(abs(d$Mantel_r) >= 0.58, "white", "black")
  ggplot2::ggplot(d, ggplot2::aes(X, Y, fill = Mantel_r)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = Label, colour = Text_colour), size = 3.5, lineheight = 1.0, show.legend = FALSE) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
      labels = c("−1", "−0.5", "0", "0.5", "1"), name = "Spearman Mantel ρ") +
    ggplot2::coord_cartesian(expand = FALSE) + ggplot2::labs(x = NULL, y = NULL) + f8_theme() +
    ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 10.5, lineheight = 1.05),
      axis.text.y = ggplot2::element_text(size = 10.5, lineheight = 1.05),
      legend.position = "bottom", legend.title = ggplot2::element_text(size = 10.5)) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(direction = "horizontal", title.position = "top",
      title.hjust = 0.5, barwidth = grid::unit(55, "mm"), barheight = grid::unit(3.5, "mm")))
}
f8_ggsave <- function(...) {
  # Bypass the global large-font export override; keep the figure's own layout.
  ggplot2::ggsave(...)
}
f8_export <- function(plot, folder, stem, width, height) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  for (extension in c("png", "tiff", "pdf")) {
    args <- list(filename = file.path(folder, paste0(stem, ".", extension)), plot = plot,
                 width = width, height = height, units = "in", dpi = 600, bg = "white")
    if (extension == "tiff") args$compression <- "lzw"
    if (extension == "pdf") {
      args$device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
      if (!capabilities("cairo")) args$useDingbats <- FALSE
    }
    do.call(ggplot2::ggsave, args)
  }
}
f8_combine <- function(a, b) {
  patchwork::wrap_plots(a + ggplot2::labs(tag = "(a)"), b + ggplot2::labs(tag = "(b)"),
                        ncol = 1, heights = c(1, 1.8), guides = "keep")
}

.frames <- vapply(sys.frames(), function(x) if (is.character(x$ofile) && length(x$ofile)==1L) x$ofile else "", "")
.frames <- .frames[nzchar(.frames)]
.rawargs <- commandArgs(FALSE)
.args <- sub("^--file=", "", grep("^--file=", .rawargs, value=TRUE))
if(!length(.args)) {
  .fileflag <- which(.rawargs %in% c("-f", "--file"))
  if(length(.fileflag) && .fileflag[1L]<length(.rawargs)) .args <- .rawargs[.fileflag[1L]+1L]
}
ROOT <- if(length(.frames)) dirname(normalizePath(tail(.frames,1L))) else if(length(.args)) dirname(normalizePath(.args[1L])) else getwd()
if (!length(.frames) && !length(.args) && requireNamespace("rstudioapi", quietly=TRUE)) {
  .editor <- rstudioapi::getActiveDocumentContext()$path
  if(nzchar(.editor)) ROOT <- dirname(.editor)
}
packages <- c("readxl", "ggplot2", "patchwork", "dplyr", "tidyr")
missing <- packages[!vapply(packages, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]
if(length(missing)) utils::install.packages(missing, repos="https://cloud.r-project.org")
OUTPUT <- file.path(ROOT, "Figures_layout_corrected")
dir.create(OUTPUT, recursive=TRUE, showWarnings=FALSE)
if(!nzchar(WORKBOOK)) {
  .candidates <- file.path(ROOT, c("Han_River_Figure8_Analysis_Ready_Reviewer3.xlsx", "Han_River_Figure8_Analysis_Ready_Reviewer3(1).xlsx"))
  .candidates <- .candidates[file.exists(.candidates)]
  if(length(.candidates)) WORKBOOK <- .candidates[1L] else if(interactive()) WORKBOOK <- file.choose()
}
if(!nzchar(WORKBOOK) || !file.exists(WORKBOOK)) stop("Set WORKBOOK to your Figure 8 Excel input, then Source again.")
environment <- as.data.frame(readxl::read_excel(WORKBOOK, sheet="Environmental_Transects"))
needed <- c("Reach_code", "Campaign", f8_variables)
if(!all(needed %in% names(environment))) stop("Environmental_Transects has missing required columns.")
if(nrow(environment)!=18L) stop("Expected the 18 transect rows used by the existing analysis.")
raw <- tidyr::pivot_longer(environment, cols=tidyr::all_of(f8_variables), names_to="Variable", values_to="Value")
summary <- dplyr::summarise(dplyr::group_by(raw, Reach_code, Campaign, Variable),
  Mean=mean(Value), SE=stats::sd(Value)/sqrt(dplyr::n()), .groups="drop")
f8_export(f8_environment_figure(summary, raw), OUTPUT,
  "Figure_3_selected_environmental_variables", 7.48, 6.1)
cat("Figure 3 saved in:", OUTPUT, "\n")
required_tables <- c("Table_8_revised_trophic_position_group_summary.csv",
  "Table_10_exact_environmental_biological_Mantel_results.csv", "Table_S8_domain_specific_Mantel_results.csv")
if(!nzchar(TABLE_DIR)) {
  .search_roots <- unique(c(ROOT, dirname(normalizePath(WORKBOOK))))
  .tables <- unique(unlist(lapply(.search_roots, function(r) list.files(r,
    pattern="^Table_8_revised_trophic_position_group_summary[.]csv$", recursive=TRUE, full.names=TRUE))))
  .folders <- unique(dirname(.tables))
  .folders <- .folders[vapply(.folders, function(f) all(file.exists(file.path(f,required_tables))), logical(1))]
  if(length(.folders)) TABLE_DIR <- .folders[which.max(vapply(.folders,
    function(f) as.numeric(file.info(file.path(f,required_tables[1L]))$mtime), numeric(1)))]
}
if(!nzchar(TABLE_DIR) || !all(file.exists(file.path(TABLE_DIR, required_tables)))) {
  message("Figure 3 is complete. To regenerate Figure 8, set TABLE_DIR at the top to the Tables folder from your completed Figure 8 analysis, then Source again.")
} else {
  tables <- lapply(file.path(TABLE_DIR,required_tables), utils::read.csv, check.names=FALSE, stringsAsFactors=FALSE)
  tp <- tables[[1L]]
  tp_required <- c("Reach_code","Campaign","Point_TP","Lower_95_TP","Upper_95_TP","Fish_n","Fish_mean_d15N","Baseline_mean_d15N")
  if(!all(tp_required %in% names(tp))) stop("The selected TP table does not have the current analysis columns. Rerun the full corrected Figure 8 analysis first.")
  if(nrow(tp)!=6L || sum(tp$Fish_n)!=609L) stop("The TP table does not cover all six groups and 609 fish.")
  fish <- as.data.frame(readxl::read_excel(WORKBOOK,sheet="Fish_Isotopes"))
  fish_key <- paste(fish$Reach_code,fish$Campaign,sep="|"); tp_key <- paste(tp$Reach_code,tp$Campaign,sep="|")
  for(i in seq_len(nrow(tp))) {
    values <- fish$d15N_permil[fish_key==tp_key[i]]
    if(length(values)!=tp$Fish_n[i] || abs(mean(values)-tp$Fish_mean_d15N[i])>1e-8)
      stop("The saved TP table does not match the current Excel fish records. Rerun the full analysis before plotting.")
  }
  baselines <- as.data.frame(readxl::read_excel(WORKBOOK,sheet="TP_Baselines"))
  bidx <- match(tp_key, baselines$Fish_group_key)
  if(anyNA(bidx) || any(abs(tp$Baseline_mean_d15N-baselines$Mean_d15N_permil[bidx])>1e-8))
    stop("The saved TP table does not match the current baseline inputs. Rerun the full analysis before plotting.")
  if(any(abs(tp$Point_TP-(2+(tp$Fish_mean_d15N-tp$Baseline_mean_d15N)/3.4))>1e-8))
    stop("The TP table does not use the current fixed-TDF primary calculation.")
  # Tables are read together from one completed run; uncertainty intervals are retained.
  m <- dplyr::bind_rows(tables[[2L]],tables[[3L]])
  needed_m <- c("Environmental_metric","Biological_metric","Mantel_r")
  if(!all(needed_m %in% names(m)) || nrow(m)!=28L ||
     anyDuplicated(paste(m$Environmental_metric,m$Biological_metric)) || any(!is.finite(m$Mantel_r)))
    stop("Expected 28 distinct finite Mantel correlations from one completed run.")
  a <- f8_tp_figure(tp); b <- f8_mantel_figure(m)
  f8_export(a,OUTPUT,"Figure_5_revised_trophic_position",7.48,3.55)
  f8_export(b,OUTPUT,"Figure_6_environmental_biological_Mantel_heatmap",7.48,5.25)
  f8_export(f8_combine(a,b),OUTPUT,"Figure_8_gastropod_TP_environmental_concordance",7.48,8.6)
  utils::write.csv(data.frame(Input=normalizePath(c(WORKBOOK,file.path(TABLE_DIR,required_tables)),winslash="/"),
    MD5=unname(tools::md5sum(c(WORKBOOK,file.path(TABLE_DIR,required_tables))))),
    file.path(OUTPUT,"PLOT_SOURCE_FILES.csv"),row.names=FALSE)
  cat("Figures 3, 5, 6 and 8 regenerated from existing results in:",OUTPUT,"\n")
}

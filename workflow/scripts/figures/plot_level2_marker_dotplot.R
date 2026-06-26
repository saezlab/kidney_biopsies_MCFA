if (exists("snakemake")) {
  rds_path <- snakemake@input[["rds"]]
  markers_csv <- snakemake@input[["markers_csv"]]
  annotation_csv <- snakemake@input[["annotation_csv"]]
  plot_path <- snakemake@output[["plot"]]
  cluster_column <- snakemake@params[["cluster_column"]]
  level1_column <- snakemake@params[["level1_column"]]
  level2_column <- snakemake@params[["level2_column"]]
  top_n <- as.integer(snakemake@params[["top_n"]])
} else {
  rds_path <- "results/seurat/subset_high_quality/sct_harmony.rds"
  markers_csv <- "results/seurat/subset_high_quality/markers/findallmarkers_level2.csv"
  annotation_csv <- "config/RNA_level1_annotations.csv"
  plot_path <- "plots/figures/Figure1/level2_marker_dotplot.pdf"
  cluster_column <- "SCT_snn_res.0.6"
  level1_column <- "level1"
  level2_column <- "level2"
  top_n <- 3L
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(Seurat)
  library(tibble)
})

level1_panel_order <- c("PT", "TAL", "DCT", "CNT_PC", "IC", "EC", "POD", "PEC", "Stromal", "Immune")
level2_plot_order <- c(
  "PT", "aPT", "dPT", "TAL", "DCT", "CNT", "PC",
  "IC-A", "IC-B", "EC", "POD", "PEC", "MES", "VSM-P", "FIB", "Immune"
)

clean_label <- function(x) {
  x <- trimws(as.character(x))
  x <- dplyr::na_if(x, "")
  dplyr::na_if(x, "NA")
}

select_top_markers <- function(markers_tbl, top_n) {
  markers_sig <- markers_tbl %>%
    dplyr::filter(!is.na(.data$p_val_adj), .data$p_val_adj < 0.05)

  if (nrow(markers_sig) == 0L) {
    return(markers_sig[0, , drop = FALSE])
  }

  markers_sig %>%
    dplyr::arrange(.data$cluster, dplyr::desc(.data$avg_log2FC), .data$p_val_adj, .data$gene) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.data$cluster, .data$gene, .keep_all = TRUE) %>%
    dplyr::filter(!duplicated(.data$gene))
}

prepare_marker_panel <- function(marker_panel, feature_universe) {
  marker_panel_present <- lapply(marker_panel, function(features) {
    intersect(unique(as.character(features)), feature_universe)
  })
  marker_panel_present <- marker_panel_present[lengths(marker_panel_present) > 0L]

  if (length(marker_panel_present) == 0L) {
    stop("No marker genes were present in the object.", call. = FALSE)
  }

  marker_panel_present
}

flatten_features <- function(marker_panel) {
  unique(unlist(marker_panel, use.names = FALSE))
}

build_marker_panel_from_markers <- function(marker_tbl, annotation_df, level1_levels, level2_levels) {
  marker_tbl <- marker_tbl %>%
    dplyr::left_join(
      annotation_df %>% dplyr::distinct(level2, level1),
      by = c("cluster" = "level2")
    ) %>%
    dplyr::filter(!is.na(level1)) %>%
    dplyr::mutate(
      level1 = factor(.data$level1, levels = level1_levels),
      cluster = factor(.data$cluster, levels = level2_levels)
    ) %>%
    dplyr::arrange(.data$level1, .data$cluster, dplyr::desc(.data$avg_log2FC), .data$p_val_adj, .data$gene)

  if (nrow(marker_tbl) == 0L) {
    return(list())
  }

  marker_split <- split(as.character(marker_tbl$gene), as.character(marker_tbl$level1))
  marker_split <- marker_split[level1_levels[level1_levels %in% names(marker_split)]]
  marker_split[lengths(marker_split) > 0L]
}

save_dotplot <- function(obj, marker_panel, group_col, out_path) {
  p_dot <- Seurat::DotPlot(
    object = obj,
    features = marker_panel,
    assay = "SCT",
    group.by = group_col,
    dot.min = 0,
    dot.scale = 5,
    scale = TRUE,
    cols = c("lightgrey", "darkblue")
  ) +
    Seurat::RotatedAxis() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 8),
      axis.text.y = ggplot2::element_text(size = 9),
      axis.title.x = ggplot2::element_text(size = 10),
      axis.title.y = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 10)
    )

  plot_width <- min(30, max(12, 0.35 * length(flatten_features(marker_panel))))
  plot_height <- max(8, 0.35 * length(levels(obj[[group_col, drop = TRUE]])))
  ggplot2::ggsave(
    filename = out_path,
    plot = p_dot,
    width = plot_width,
    height = plot_height,
    units = "in"
  )
}

required_marker_cols <- c("cluster", "gene", "avg_log2FC", "p_val_adj")
required_annotation_cols <- c(cluster_column, level1_column, level2_column)

obj <- readRDS(rds_path)
markers_tbl <- readr::read_csv(markers_csv, show_col_types = FALSE)
annotation_df <- utils::read.csv(
  annotation_csv,
  fill = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

missing_marker_cols <- setdiff(required_marker_cols, colnames(markers_tbl))
if (length(missing_marker_cols) > 0) {
  stop(
    "Markers CSV is missing required columns: ",
    paste(missing_marker_cols, collapse = ", "),
    call. = FALSE
  )
}

if (!cluster_column %in% colnames(obj@meta.data)) {
  stop(
    "Missing expected cluster metadata column: ",
    cluster_column,
    call. = FALSE
  )
}

missing_annotation_cols <- setdiff(required_annotation_cols, colnames(annotation_df))
if (length(missing_annotation_cols) > 0) {
  stop(
    "Annotation file is missing required columns: ",
    paste(missing_annotation_cols, collapse = ", "),
    call. = FALSE
  )
}

annotation_df <- annotation_df %>%
  dplyr::select(
    cluster = dplyr::all_of(cluster_column),
    level1 = dplyr::all_of(level1_column),
    level2 = dplyr::all_of(level2_column)
  ) %>%
  dplyr::mutate(
    cluster = clean_label(cluster),
    level1 = clean_label(level1),
    level2 = clean_label(level2)
  ) %>%
  dplyr::filter(!is.na(cluster), !is.na(level1), !is.na(level2)) %>%
  dplyr::distinct(cluster, .keep_all = TRUE)

observed_level1 <- unique(annotation_df$level1)
observed_level2 <- unique(annotation_df$level2)
if (length(observed_level2) == 0L) {
  stop("No complete level2 annotations were found.", call. = FALSE)
}

missing_level1_order <- setdiff(observed_level1, level1_panel_order)
if (length(missing_level1_order) > 0L) {
  stop(
    "Missing level1 panel order entries for: ",
    paste(missing_level1_order, collapse = ", "),
    call. = FALSE
  )
}

missing_level2_order <- setdiff(observed_level2, level2_plot_order)
if (length(missing_level2_order) > 0L) {
  stop(
    "Missing level2 plot order entries for: ",
    paste(missing_level2_order, collapse = ", "),
    call. = FALSE
  )
}

level1_levels <- level1_panel_order[level1_panel_order %in% observed_level1]
level2_levels <- level2_plot_order[level2_plot_order %in% observed_level2]
level2_axis_levels <- rev(level2_levels)

cell_groups <- tibble::tibble(
  cell = colnames(obj),
  cluster = clean_label(obj[[cluster_column, drop = TRUE]])
) %>%
  dplyr::left_join(annotation_df, by = "cluster") %>%
  dplyr::filter(!is.na(level1), !is.na(level2))

if (nrow(cell_groups) == 0L) {
  stop("No cells had complete level1/level2 annotations after cluster mapping.", call. = FALSE)
}

obj <- obj[, cell_groups$cell]
level2_lookup <- stats::setNames(cell_groups$level2, cell_groups$cell)
level1_lookup <- stats::setNames(cell_groups$level1, cell_groups$cell)
obj[[level2_column]] <- factor(
  unname(level2_lookup[colnames(obj)]),
  levels = level2_axis_levels
)
obj[[level1_column]] <- factor(
  unname(level1_lookup[colnames(obj)]),
  levels = level1_levels
)

top_markers <- select_top_markers(markers_tbl, top_n = top_n)
if (nrow(top_markers) == 0L) {
  stop("No significant level2 markers were available to plot.", call. = FALSE)
}

marker_panel <- build_marker_panel_from_markers(top_markers, annotation_df, level1_levels, level2_levels)
marker_panel <- prepare_marker_panel(marker_panel, rownames(obj))

if (!exists("snakemake")) {
  dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
}

save_dotplot(
  obj = obj,
  marker_panel = marker_panel,
  group_col = level2_column,
  out_path = plot_path
)

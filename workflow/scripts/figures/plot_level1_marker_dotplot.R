if (exists("snakemake")) {
  rds_path <- snakemake@input[["rds"]]
  markers_csv <- snakemake@input[["markers_csv"]]
  annotation_csv <- snakemake@input[["annotation_csv"]]
  plot_path <- snakemake@output[["plot"]]
  cluster_column <- snakemake@params[["cluster_column"]]
  group_column <- snakemake@params[["group_column"]]
  top_n <- as.integer(snakemake@params[["top_n"]])
} else {
  rds_path <- "results/seurat/subset_high_quality/sct_harmony.rds"
  markers_csv <- "results/seurat/subset_high_quality/markers/findallmarkers_level1.csv"
  annotation_csv <- "config/RNA_level1_annotations.csv"
  plot_path <- "plots/figures/Figure1/level1_marker_dotplot.pdf"
  cluster_column <- "SCT_snn_res.0.6"
  group_column <- "level1"
  top_n <- 3L
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(Seurat)
  library(tibble)
})

level1_plot_order <- c("PT", "TAL", "DCT", "CNT_PC", "IC", "EC", "POD", "PEC", "Stromal", "Immune")

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
    dplyr::arrange(dplyr::desc(.data$avg_log2FC), .data$p_val_adj, .data$gene, .data$cluster) %>%
    dplyr::filter(!duplicated(.data$gene))
}

build_marker_panel_from_markers <- function(marker_tbl, group_levels) {
  if (nrow(marker_tbl) == 0L) {
    return(list())
  }

  marker_tbl <- marker_tbl %>%
    dplyr::mutate(cluster = factor(.data$cluster, levels = group_levels)) %>%
    dplyr::arrange(.data$cluster)

  marker_split <- split(as.character(marker_tbl$gene), as.character(marker_tbl$cluster))
  marker_split <- marker_split[group_levels[group_levels %in% names(marker_split)]]
  marker_split[lengths(marker_split) > 0L]
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
      axis.text.x = ggplot2::element_text(size = 9),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.title.x = ggplot2::element_text(size = 10),
      axis.title.y = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 10)
    )

  plot_width <- min(24, max(10, 0.45 * length(flatten_features(marker_panel))))
  ggplot2::ggsave(
    filename = out_path,
    plot = p_dot,
    width = plot_width,
    height = 8,
    units = "in"
  )
}

required_marker_cols <- c("cluster", "gene", "avg_log2FC", "p_val_adj")
required_annotation_cols <- c(cluster_column, group_column)

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
    group = dplyr::all_of(group_column)
  ) %>%
  dplyr::mutate(
    cluster = clean_label(cluster),
    group = clean_label(group)
  ) %>%
  dplyr::filter(!is.na(cluster), !is.na(group)) %>%
  dplyr::distinct(cluster, .keep_all = TRUE)

observed_groups <- unique(annotation_df$group)
if (length(observed_groups) == 0L) {
  stop("No complete level1 annotations were found.", call. = FALSE)
}

missing_order_groups <- setdiff(observed_groups, level1_plot_order)
if (length(missing_order_groups) > 0L) {
  stop(
    "Missing level1 plot order entries for: ",
    paste(missing_order_groups, collapse = ", "),
    call. = FALSE
  )
}

group_panel_levels <- level1_plot_order[level1_plot_order %in% observed_groups]
group_axis_levels <- rev(group_panel_levels)

cell_groups <- tibble::tibble(
  cell = colnames(obj),
  cluster = clean_label(obj[[cluster_column, drop = TRUE]])
) %>%
  dplyr::left_join(annotation_df, by = "cluster") %>%
  dplyr::filter(!is.na(group))

if (nrow(cell_groups) == 0L) {
  stop("No cells had complete level1 annotations after cluster mapping.", call. = FALSE)
}

obj <- obj[, cell_groups$cell]
cell_group_lookup <- stats::setNames(cell_groups$group, cell_groups$cell)
obj[[group_column]] <- factor(
  unname(cell_group_lookup[colnames(obj)]),
  levels = group_axis_levels
)

top_markers <- select_top_markers(markers_tbl, top_n = top_n)
if (nrow(top_markers) == 0L) {
  stop("No significant level1 markers were available to plot.", call. = FALSE)
}

marker_panel <- build_marker_panel_from_markers(top_markers, group_panel_levels)
marker_panel <- prepare_marker_panel(marker_panel, rownames(obj))

if (!exists("snakemake")) {
  dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
}

save_dotplot(
  obj = obj,
  marker_panel = marker_panel,
  group_col = group_column,
  out_path = plot_path
)

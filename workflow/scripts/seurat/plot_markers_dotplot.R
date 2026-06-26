if (exists("snakemake")) {
  rds_path <- snakemake@input[["rds"]]
  markers_csv <- snakemake@input[["csv"]]
  plot_path <- snakemake@output[["plot"]]
  resolution <- snakemake@params[["resolution"]]
} else {
  rds_path <- "results/seurat/subset_high_quality/sct_harmony.rds"
  markers_csv <- "results/seurat/subset_high_quality/markers/findallmarkers_res_0.6.csv"
  plot_path <- "plots/seurat/subset_high_quality/dotplot_res_0.6.pdf"
  resolution <- "0.6"
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

select_top_markers <- function(markers_tbl, top_n = 3L) {
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

build_marker_panel_from_markers <- function(marker_tbl) {
  if (nrow(marker_tbl) == 0L) {
    return(list())
  }

  marker_split <- split(as.character(marker_tbl$gene), as.character(marker_tbl$cluster))
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
cluster_col <- paste0("SCT_snn_res.", resolution)

obj <- readRDS(rds_path)
markers_tbl <- readr::read_csv(markers_csv, show_col_types = FALSE)

missing_marker_cols <- setdiff(required_marker_cols, colnames(markers_tbl))
if (length(missing_marker_cols) > 0) {
  stop(
    "Markers CSV is missing required columns: ",
    paste(missing_marker_cols, collapse = ", "),
    call. = FALSE
  )
}

if (!cluster_col %in% colnames(obj@meta.data)) {
  stop(
    "Missing expected cluster metadata column: ",
    cluster_col,
    call. = FALSE
  )
}

top_markers <- select_top_markers(markers_tbl, top_n = 3L)
if (nrow(top_markers) == 0L) {
  stop(
    "No significant markers were available to plot for resolution ",
    resolution,
    ".",
    call. = FALSE
  )
}

marker_panel <- build_marker_panel_from_markers(top_markers)
marker_panel <- prepare_marker_panel(marker_panel, rownames(obj))

if (!exists("snakemake")) {
  dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
}

save_dotplot(
  obj = obj,
  marker_panel = marker_panel,
  group_col = cluster_col,
  out_path = plot_path
)

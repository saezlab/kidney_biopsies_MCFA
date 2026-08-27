if (exists("snakemake")) {
  metadata_csv <- snakemake@input[["metadata"]]
  annotation_csv <- snakemake@input[["annotation_csv"]]
  pdf_path <- snakemake@output[["pdf"]]
  png_path <- snakemake@output[["png"]]
  cluster_column <- snakemake@params[["cluster_column"]]
  level2_column <- snakemake@params[["level2_column"]]
} else {
  metadata_csv <- "results/seurat/subset_high_quality/metadata.csv"
  annotation_csv <- "config/RNA_level1_annotations.csv"
  pdf_path <- "plots/figures/Figure1Ext/level2_qc_metrics.pdf"
  png_path <- "plots/figures/Figure1Ext/level2_qc_metrics.png"
  cluster_column <- "SCT_snn_res.0.6"
  level2_column <- "level2"
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
})

qc_metrics <- c(
  nCount_RNA = "RNA counts",
  nFeature_RNA = "RNA features",
  pct_counts_mt = "Mitochondrial counts (%)",
  scDblFinder.score = "scDblFinder score"
)
log10_metrics <- c("nCount_RNA", "nFeature_RNA", "scDblFinder.score")
level2_plot_order <- c(
  "PT", "aPT", "dPT", "TAL", "DCT", "CNT", "PC",
  "IC-A", "IC-B", "EC", "POD", "PEC", "MES", "VSM-P", "FIB", "Immune"
)

clean_label <- function(x) {
  x <- trimws(as.character(x))
  x <- dplyr::na_if(x, "")
  dplyr::na_if(x, "NA")
}

metadata_df <- readr::read_csv(metadata_csv, show_col_types = FALSE)
annotation_df <- utils::read.csv(
  annotation_csv,
  fill = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_metadata_cols <- c(cluster_column, names(qc_metrics))
missing_metadata_cols <- setdiff(required_metadata_cols, colnames(metadata_df))
if (length(missing_metadata_cols) > 0L) {
  stop(
    "Metadata file is missing required columns: ",
    paste(missing_metadata_cols, collapse = ", "),
    call. = FALSE
  )
}

required_annotation_cols <- c(cluster_column, level2_column)
missing_annotation_cols <- setdiff(required_annotation_cols, colnames(annotation_df))
if (length(missing_annotation_cols) > 0L) {
  stop(
    "Annotation file is missing required columns: ",
    paste(missing_annotation_cols, collapse = ", "),
    call. = FALSE
  )
}

annotation_df <- annotation_df %>%
  dplyr::select(
    cluster = dplyr::all_of(cluster_column),
    level2 = dplyr::all_of(level2_column)
  ) %>%
  dplyr::mutate(
    cluster = clean_label(.data$cluster),
    level2 = clean_label(.data$level2)
  ) %>%
  dplyr::filter(!is.na(.data$cluster), !is.na(.data$level2)) %>%
  dplyr::distinct(.data$cluster, .keep_all = TRUE)

plot_df <- metadata_df %>%
  dplyr::select(cluster = dplyr::all_of(cluster_column), dplyr::all_of(names(qc_metrics))) %>%
  dplyr::mutate(cluster = clean_label(.data$cluster)) %>%
  dplyr::filter(!is.na(.data$cluster)) %>%
  dplyr::left_join(annotation_df, by = "cluster") %>%
  dplyr::filter(!is.na(.data$level2))

if (nrow(plot_df) == 0L) {
  stop("No cells had complete level2 annotations after cluster mapping.", call. = FALSE)
}

observed_level2 <- unique(plot_df$level2)
missing_level2_order <- setdiff(observed_level2, level2_plot_order)
if (length(missing_level2_order) > 0L) {
  stop(
    "Missing level2 plot order entries for: ",
    paste(missing_level2_order, collapse = ", "),
    call. = FALSE
  )
}

level2_levels <- level2_plot_order[level2_plot_order %in% observed_level2]
plot_df <- plot_df %>%
  dplyr::mutate(level2 = factor(.data$level2, levels = level2_levels))

base_theme <- theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.margin = margin(6, 6, 6, 6)
  )

metric_fills <- c(
  nCount_RNA = "#4C956C",
  nFeature_RNA = "#7CAFC2",
  pct_counts_mt = "#F4A259",
  scDblFinder.score = "#B56576"
)

make_metric_plot <- function(metric) {
  use_log10 <- metric %in% log10_metrics
  axis_label <- unname(qc_metrics[[metric]])
  if (use_log10) {
    axis_label <- paste0(axis_label, " (log10 scale)")
  }

  p <- ggplot(plot_df, aes(x = .data$level2, y = .data[[metric]])) +
    geom_boxplot(
      fill = metric_fills[[metric]],
      color = "#333333",
      width = 0.65,
      outlier.shape = 16,
      outlier.size = 0.25,
      outlier.alpha = 0.25,
      linewidth = 0.25
    ) +
    labs(
      title = axis_label,
      y = axis_label
    ) +
    base_theme

  if (use_log10) {
    p <- p + scale_y_log10()
  }

  p
}

plots <- lapply(names(qc_metrics), make_metric_plot)

combined_plot <- (plots[[1]] | plots[[2]]) /
  (plots[[3]] | plots[[4]]) +
  patchwork::plot_annotation(
    title = "QC metrics by level2",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

if (!exists("snakemake")) {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
}

ggsave(pdf_path, combined_plot, width = 12, height = 10, units = "in")
ggsave(png_path, combined_plot, width = 12, height = 10, units = "in", dpi = 300)

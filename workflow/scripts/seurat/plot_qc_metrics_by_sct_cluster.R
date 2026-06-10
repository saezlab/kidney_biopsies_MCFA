if (exists("snakemake")) {
  metadata_csv <- snakemake@input[["metadata"]]
  plot_path <- snakemake@output[["plot"]]
  cluster_column <- snakemake@params[["cluster_column"]]
} else {
  metadata_csv <- "results/seurat/metadata.csv"
  plot_path <- "plots/seurat/annotation/qc_metrics_by_sct_res_0_2.pdf"
  cluster_column <- "SCT_snn_res.0.2"
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
required_cols <- c(cluster_column, names(qc_metrics))

metadata_df <- readr::read_csv(metadata_csv, show_col_types = FALSE)

missing_cols <- setdiff(required_cols, colnames(metadata_df))
if (length(missing_cols) > 0) {
  stop(
    "Metadata file is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

plot_df <- metadata_df %>%
  dplyr::select(cluster = dplyr::all_of(cluster_column), dplyr::all_of(names(qc_metrics))) %>%
  dplyr::filter(!is.na(cluster), cluster != "") %>%
  dplyr::mutate(cluster = as.character(cluster))

if (nrow(plot_df) == 0) {
  stop("No cells with non-missing cluster labels were found in the metadata file.", call. = FALSE)
}

cluster_levels <- unique(plot_df$cluster)
cluster_levels <- cluster_levels[order(suppressWarnings(as.numeric(cluster_levels)), cluster_levels)]
plot_df <- plot_df %>%
  dplyr::mutate(cluster = factor(cluster, levels = cluster_levels))

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

  p <- ggplot(plot_df, aes(x = cluster, y = .data[[metric]])) +
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
    title = paste("QC metrics by", cluster_column),
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(plot_path, combined_plot, width = 12, height = 10)

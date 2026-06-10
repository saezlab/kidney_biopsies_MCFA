if (exists("snakemake")) {
  metadata_csv <- snakemake@input[["metadata"]]
  annotation_csv <- snakemake@input[["annotation_csv"]]
  plot_path <- snakemake@output[["plot"]]
  cluster_column <- snakemake@params[["cluster_column"]]
} else {
  metadata_csv <- "results/seurat/metadata.csv"
  annotation_csv <- "config/biopsydata_withkpmp_final_annotations(in).csv"
  plot_path <- "plots/seurat/annotation/main_class_proportions_by_sct_res_0_2.pdf"
  cluster_column <- "SCT_snn_res.0.2"
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

required_metadata_cols <- c("cellid", cluster_column)
required_annotation_cols <- c("barcode_combined_id", "main_class")
unannotated_label <- "Unannotated"
unannotated_color <- "#BDBDBD"
palette_21 <- c(
  "#e6194B", "#4363d8", "#ffe119", "#3cb44b", "#f58231", "#911eb4",
  "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#469990", "#dcbeff",
  "#9A6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1",
  "#000075", "#a9a9a9", "#000000"
)

label_sort_key <- function(labels) {
  prefix <- suppressWarnings(as.numeric(sub("_.*", "", labels)))
  prefix[is.na(prefix)] <- Inf
  prefix
}

metadata_df <- readr::read_csv(metadata_csv, show_col_types = FALSE)
annotation_df <- readr::read_csv(annotation_csv, show_col_types = FALSE)

missing_metadata_cols <- setdiff(required_metadata_cols, colnames(metadata_df))
if (length(missing_metadata_cols) > 0) {
  stop(
    "Metadata file is missing required columns: ",
    paste(missing_metadata_cols, collapse = ", "),
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
  dplyr::select(barcode_combined_id, main_class) %>%
  dplyr::filter(!is.na(barcode_combined_id), barcode_combined_id != "") %>%
  dplyr::mutate(
    main_class = trimws(main_class),
    main_class = dplyr::na_if(main_class, ""),
    main_class = dplyr::na_if(main_class, "NA")
  ) %>%
  dplyr::distinct(barcode_combined_id, .keep_all = TRUE)

class_levels <- annotation_df %>%
  dplyr::filter(!is.na(main_class)) %>%
  dplyr::distinct(main_class) %>%
  dplyr::pull(main_class)
class_levels <- class_levels[order(label_sort_key(class_levels), class_levels)]

if (length(class_levels) != length(palette_21)) {
  stop(
    "Expected 21 main_class levels to match palette_21, found ",
    length(class_levels),
    ".",
    call. = FALSE
  )
}

plot_df <- metadata_df %>%
  dplyr::select(cellid, cluster = dplyr::all_of(cluster_column)) %>%
  dplyr::left_join(annotation_df, by = c("cellid" = "barcode_combined_id")) %>%
  dplyr::mutate(
    cluster = as.character(cluster),
    main_class = dplyr::coalesce(main_class, unannotated_label),
    main_class = factor(main_class, levels = c(class_levels, unannotated_label))
  )

if (nrow(plot_df) == 0) {
  stop("No cells were found in the metadata file.", call. = FALSE)
}

cluster_levels <- unique(plot_df$cluster)
cluster_levels <- cluster_levels[order(suppressWarnings(as.numeric(cluster_levels)), cluster_levels)]

proportion_df <- plot_df %>%
  dplyr::filter(!is.na(cluster), cluster != "") %>%
  dplyr::count(cluster, main_class, name = "n") %>%
  tidyr::complete(
    cluster = cluster_levels,
    main_class = factor(c(class_levels, unannotated_label), levels = c(class_levels, unannotated_label)),
    fill = list(n = 0)
  ) %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(proportion = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cluster = factor(cluster, levels = cluster_levels))

class_colors <- stats::setNames(c(palette_21, unannotated_color), c(class_levels, unannotated_label))

p <- ggplot(proportion_df, aes(x = cluster, y = proportion, fill = main_class)) +
  geom_col(width = 0.85, color = "white", linewidth = 0.1) +
  scale_fill_manual(values = class_colors, drop = FALSE) +
  scale_y_continuous(
    labels = function(x) paste0(round(x * 100), "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = paste("main_class proportions by", cluster_column),
    x = cluster_column,
    y = "Proportion of nuclei",
    fill = "main_class"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    legend.key.height = grid::unit(0.45, "cm"),
    legend.text = element_text(size = 8),
    plot.title = element_text(face = "bold")
  ) +
  guides(fill = guide_legend(ncol = 1))

grDevices::pdf(plot_path, width = 12, height = 7)
print(p)
grDevices::dev.off()

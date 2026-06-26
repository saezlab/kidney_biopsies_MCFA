if (exists("snakemake")) {
  umap_csv <- snakemake@input[["umap"]]
  metadata_csv <- snakemake@input[["metadata"]]
  annotation_csv <- snakemake@input[["annotation_csv"]]
  plot_path <- snakemake@output[["plot"]]
  cluster_column <- snakemake@params[["cluster_column"]]
} else {
  umap_csv <- "results/seurat/subset_high_quality/umap_embeddings.csv"
  metadata_csv <- "results/seurat/subset_high_quality/metadata.csv"
  annotation_csv <- "config/RNA_level1_annotations.csv"
  plot_path <- "plots/figures/Figure1/rna_umap_level1_level2.pdf"
  cluster_column <- "SCT_snn_res.0.6"
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(png)
  library(readr)
  library(ragg)
  library(tibble)
})

required_umap_cols <- c("cellid", "umap_1", "umap_2")
required_metadata_cols <- c("cellid", cluster_column)
required_annotation_cols <- c(cluster_column, "level1", "level2")
point_colors <- c(
  "PT||PT" = "#C05A3A",
  "PT||aPT" = "#D66F45",
  "PT||dPT" = "#8F4B2F",
  "TAL||TAL" = "#C2B35A",
  "DCT||DCT" = "#3F8D3A",
  "CNT_PC||CNT" = "#6BD868",
  "CNT_PC||PC" = "#4FAE78",
  "IC||IC-A" = "#3D8F93",
  "IC||IC-B" = "#8AD5CC",
  "Stromal||FIB" = "#7C55A5",
  "Stromal||VSM-P" = "#5E3E88",
  "Stromal||MES" = "#B38AD5",
  "EC||EC" = "#5E8FBA",
  "Immune||Immune" = "#C83FA8",
  "POD||POD" = "#D8908C",
  "PEC||PEC" = "#E8B1AA"
)

clean_label <- function(x) {
  x <- trimws(as.character(x))
  x <- dplyr::na_if(x, "")
  dplyr::na_if(x, "NA")
}

render_point_raster <- function(plot_data, point_colors, x_limits, y_limits) {
  tmp_png <- tempfile(fileext = ".png")
  on.exit(unlink(tmp_png), add = TRUE)

  ragg::agg_png(
    filename = tmp_png,
    width = 1800,
    height = 1500,
    units = "px",
    background = "transparent",
    res = 220
  )

  point_plot <- ggplot(plot_data, aes(x = umap_1, y = umap_2, color = color_key)) +
    geom_point(size = 0.025, alpha = 0.25) +
    scale_color_manual(values = point_colors, drop = FALSE) +
    coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.position = "none"
    )

  print(point_plot)
  grDevices::dev.off()
  png::readPNG(tmp_png)
}

add_text_layer <- function(plot, data, fontface, size, color) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    plot +
      ggrepel::geom_text_repel(
        data = data,
        aes(x = umap_1, y = umap_2, label = label),
        inherit.aes = FALSE,
        fontface = fontface,
        size = size,
        color = color,
        min.segment.length = 0,
        box.padding = 0.25,
        point.padding = 0.1,
        segment.color = "grey60",
        segment.size = 0.2,
        max.overlaps = Inf
      )
  } else {
    plot +
      geom_text(
        data = data,
        aes(x = umap_1, y = umap_2, label = label),
        inherit.aes = FALSE,
        fontface = fontface,
        size = size,
        color = color,
        check_overlap = TRUE
      )
  }
}

umap_df <- readr::read_csv(umap_csv, show_col_types = FALSE)
metadata_df <- readr::read_csv(metadata_csv, show_col_types = FALSE)
annotation_df <- utils::read.csv(
  annotation_csv,
  fill = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

missing_umap_cols <- setdiff(required_umap_cols, colnames(umap_df))
if (length(missing_umap_cols) > 0) {
  stop("UMAP file is missing required columns: ", paste(missing_umap_cols, collapse = ", "), call. = FALSE)
}

missing_metadata_cols <- setdiff(required_metadata_cols, colnames(metadata_df))
if (length(missing_metadata_cols) > 0) {
  stop("Metadata file is missing required columns: ", paste(missing_metadata_cols, collapse = ", "), call. = FALSE)
}

missing_annotation_cols <- setdiff(required_annotation_cols, colnames(annotation_df))
if (length(missing_annotation_cols) > 0) {
  stop("Annotation file is missing required columns: ", paste(missing_annotation_cols, collapse = ", "), call. = FALSE)
}

annotation_levels <- annotation_df %>%
  dplyr::select(
    cluster = dplyr::all_of(cluster_column),
    level1,
    level2
  ) %>%
  dplyr::mutate(
    cluster = clean_label(cluster),
    level1 = clean_label(level1),
    level2 = clean_label(level2)
  ) %>%
  dplyr::filter(!is.na(cluster), !is.na(level1), !is.na(level2)) %>%
  dplyr::distinct(cluster, .keep_all = TRUE)

if (nrow(annotation_levels) == 0L) {
  stop("No complete level1/level2 annotations were found.", call. = FALSE)
}

level1_levels <- unique(annotation_levels$level1)
observed_color_keys <- unique(paste(annotation_levels$level1, annotation_levels$level2, sep = "||"))
missing_color_keys <- setdiff(observed_color_keys, names(point_colors))
if (length(missing_color_keys) > 0) {
  stop(
    "Missing colors for level1/level2 annotations: ",
    paste(missing_color_keys, collapse = ", "),
    call. = FALSE
  )
}

plot_df <- umap_df %>%
  dplyr::inner_join(
    metadata_df %>%
      dplyr::select(cellid, cluster = dplyr::all_of(cluster_column)) %>%
      dplyr::mutate(cluster = clean_label(cluster)),
    by = "cellid"
  ) %>%
  dplyr::left_join(annotation_levels, by = "cluster") %>%
  dplyr::filter(!is.na(level1), !is.na(level2)) %>%
  dplyr::mutate(
    level1 = factor(level1, levels = level1_levels),
    color_key = paste(as.character(level1), level2, sep = "||"),
    color_key = factor(color_key, levels = names(point_colors))
  )

if (nrow(plot_df) == 0L) {
  stop("No cells remained after joining annotations and removing incomplete level1/level2 labels.", call. = FALSE)
}

x_range <- range(plot_df$umap_1, na.rm = TRUE)
y_range <- range(plot_df$umap_2, na.rm = TRUE)
x_padding <- diff(x_range) * 0.04
y_padding <- diff(y_range) * 0.04
if (!is.finite(x_padding) || x_padding == 0) {
  x_padding <- 1
}
if (!is.finite(y_padding) || y_padding == 0) {
  y_padding <- 1
}
x_limits <- c(x_range[1] - x_padding, x_range[2] + x_padding)
y_limits <- c(y_range[1] - y_padding, y_range[2] + y_padding)

level1_labels <- plot_df %>%
  dplyr::group_by(level1) %>%
  dplyr::summarise(
    umap_1 = stats::median(umap_1, na.rm = TRUE),
    umap_2 = stats::median(umap_2, na.rm = TRUE),
    label = as.character(dplyr::first(level1)),
    .groups = "drop"
  )

level2_labels <- plot_df %>%
  dplyr::group_by(level1) %>%
  dplyr::filter(dplyr::n_distinct(level2) > 1L) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(level1, level2) %>%
  dplyr::summarise(
    umap_1 = stats::median(umap_1, na.rm = TRUE),
    umap_2 = stats::median(umap_2, na.rm = TRUE),
    label = dplyr::first(level2),
    .groups = "drop"
  )

point_raster <- render_point_raster(
  plot_data = plot_df,
  point_colors = point_colors,
  x_limits = x_limits,
  y_limits = y_limits
)

if (!exists("snakemake")) {
  dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
}

grDevices::pdf(plot_path, width = 7.2, height = 7)

p <- ggplot() +
  annotation_raster(
    raster = as.raster(point_raster),
    xmin = x_limits[1],
    xmax = x_limits[2],
    ymin = y_limits[1],
    ymax = y_limits[2]
  ) +
  coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme_classic(base_size = 11) +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),
    plot.background = element_blank()
  )

p <- add_text_layer(p, level2_labels, fontface = "italic", size = 3.1, color = "grey25")
p <- add_text_layer(p, level1_labels, fontface = "bold", size = 4.2, color = "black")

print(p)
grDevices::dev.off()

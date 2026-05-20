if (exists("snakemake")) {
  umap_csv <- snakemake@input[["umap"]]
  metadata_csv <- snakemake@input[["metadata"]]
  plot_path <- snakemake@output[["plot"]]
  resolution <- snakemake@params[["resolution"]]
} else {
  umap_csv <- "results/seurat/umap_embeddings.csv"
  metadata_csv <- "results/seurat/metadata.csv"
  plot_path <- "plots/seurat/umap_res_0.2.pdf"
  resolution <- "0.2"
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(png)
  library(readr)
  library(ragg)
})

required_umap_cols <- c("cellid", "umap_1", "umap_2")
cluster_col <- paste0("SCT_snn_res.", resolution)

umap_df <- readr::read_csv(umap_csv, show_col_types = FALSE)
metadata_df <- readr::read_csv(metadata_csv, show_col_types = FALSE)

missing_umap_cols <- setdiff(required_umap_cols, colnames(umap_df))
if (length(missing_umap_cols) > 0) {
  stop(
    "UMAP embeddings file is missing required columns: ",
    paste(missing_umap_cols, collapse = ", "),
    call. = FALSE
  )
}

if (!"cellid" %in% colnames(metadata_df)) {
  stop("Metadata file is missing required column: cellid", call. = FALSE)
}

if (!cluster_col %in% colnames(metadata_df)) {
  stop(
    "Metadata file is missing required cluster column: ",
    cluster_col,
    call. = FALSE
  )
}

plot_df <- umap_df %>%
  dplyr::inner_join(
    metadata_df %>%
      dplyr::select(cellid, dplyr::all_of(cluster_col)),
    by = "cellid"
  ) %>%
  dplyr::transmute(
    umap_1 = .data$umap_1,
    umap_2 = .data$umap_2,
    cluster = as.character(.data[[cluster_col]])
  )

if (nrow(plot_df) == 0) {
  stop("Joining UMAP embeddings with metadata by cellid produced zero rows.", call. = FALSE)
}

cluster_levels <- unique(plot_df$cluster)
cluster_values <- suppressWarnings(as.numeric(cluster_levels))
if (all(!is.na(cluster_values))) {
  cluster_levels <- cluster_levels[order(cluster_values)]
} else {
  cluster_levels <- sort(cluster_levels)
}
plot_df <- plot_df %>%
  dplyr::mutate(cluster = factor(cluster, levels = cluster_levels))

base_theme <- theme(
  panel.background = element_blank(),
  plot.background = element_blank(),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  legend.position = "right",
  legend.title = element_text(size = 11),
  legend.text = element_text(size = 9),
  legend.key.height = grid::unit(0.45, "cm")
)

point_raster_width <- 1600
point_raster_height <- 1400
palette_21 <- c(
  "#e6194B", "#4363d8", "#ffe119", "#3cb44b", "#f58231", "#911eb4",
  "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#469990", "#dcbeff",
  "#9A6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1",
  "#000075", "#a9a9a9", "#000000"
)

render_point_raster <- function(plot_data, cluster_colors, x_limits, y_limits) {
  tmp_png <- tempfile(fileext = ".png")
  on.exit(unlink(tmp_png), add = TRUE)

  ragg::agg_png(
    filename = tmp_png,
    width = point_raster_width,
    height = point_raster_height,
    units = "px",
    background = "transparent",
    res = 200
  )

  point_plot <- ggplot(plot_data, aes(x = umap_1, y = umap_2, color = cluster)) +
    geom_point(size = 0.2, alpha = 0.8) +
    scale_color_manual(values = cluster_colors, drop = FALSE) +
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

cluster_palette <- rep_len(palette_21, length(cluster_levels))
cluster_colors <- stats::setNames(cluster_palette, cluster_levels)

legend_df <- plot_df %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    umap_1 = dplyr::first(umap_1),
    umap_2 = dplyr::first(umap_2),
    .groups = "drop"
  )

point_raster <- render_point_raster(
  plot_data = plot_df,
  cluster_colors = cluster_colors,
  x_limits = x_limits,
  y_limits = y_limits
)

grDevices::pdf(plot_path, width = 7.2, height = 7)

p <- ggplot() +
  annotation_raster(
    raster = as.raster(point_raster),
    xmin = x_limits[1],
    xmax = x_limits[2],
    ymin = y_limits[1],
    ymax = y_limits[2]
  ) +
  geom_point(
    data = legend_df,
    aes(x = umap_1, y = umap_2, color = cluster),
    alpha = 0,
    size = 0,
    show.legend = TRUE
  ) +
  scale_color_manual(values = cluster_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE) +
  labs(
    title = paste0("UMAP by ", cluster_col),
    x = "UMAP1",
    y = "UMAP2",
    color = "cluster"
  ) +
  theme_classic() +
  base_theme

print(p)
grDevices::dev.off()

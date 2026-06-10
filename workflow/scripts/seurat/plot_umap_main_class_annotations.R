if (exists("snakemake")) {
  umap_csv <- snakemake@input[["umap"]]
  annotation_csv <- snakemake@input[["annotation_csv"]]
  plot_path <- snakemake@output[["plot"]]
} else {
  umap_csv <- "results/seurat/umap_embeddings.csv"
  annotation_csv <- "config/biopsydata_withkpmp_final_annotations(in).csv"
  plot_path <- "plots/seurat/annotation/umap_main_class_annotations.pdf"
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
required_annotation_cols <- c("barcode_combined_id", "main_class")

umap_df <- readr::read_csv(umap_csv, show_col_types = FALSE)
annotation_df <- readr::read_csv(annotation_csv, show_col_types = FALSE)

missing_umap_cols <- setdiff(required_umap_cols, colnames(umap_df))
if (length(missing_umap_cols) > 0) {
  stop(
    "UMAP embeddings file is missing required columns: ",
    paste(missing_umap_cols, collapse = ", "),
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

plot_df <- umap_df %>%
  dplyr::left_join(annotation_df, by = c("cellid" = "barcode_combined_id"))

annotated_df <- plot_df %>%
  dplyr::filter(!is.na(main_class))

background_df <- plot_df %>%
  dplyr::filter(is.na(main_class))

if (nrow(annotated_df) == 0) {
  stop("No overlapping annotated cells were found between the UMAP and annotation files.", call. = FALSE)
}

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
background_color <- "#CFCFCF"
palette_21 <- c(
  "#e6194B", "#4363d8", "#ffe119", "#3cb44b", "#f58231", "#911eb4",
  "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#469990", "#dcbeff",
  "#9A6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1",
  "#000075", "#a9a9a9", "#000000"
)

render_point_raster <- function(plot_data, class_colors, x_limits, y_limits) {
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

  point_plot <- ggplot() +
    geom_point(
      data = background_df,
      aes(x = umap_1, y = umap_2),
      color = background_color,
      size = 0.15,
      alpha = 0.35
    ) +
    geom_point(
      data = plot_data,
      aes(x = umap_1, y = umap_2, color = main_class),
      size = 0.2,
      alpha = 0.8
    ) +
    scale_color_manual(values = class_colors, drop = FALSE) +
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

class_levels <- unique(annotated_df$main_class)
class_levels <- class_levels[order(as.numeric(sub("_.*", "", class_levels)))]
annotated_df <- annotated_df %>%
  dplyr::mutate(main_class = factor(main_class, levels = class_levels))

if (length(class_levels) != length(palette_21)) {
  stop(
    "Expected 21 main_class levels to match palette_21, found ",
    length(class_levels),
    ".",
    call. = FALSE
  )
}

class_colors <- stats::setNames(palette_21, class_levels)

legend_df <- annotated_df %>%
  dplyr::group_by(main_class) %>%
  dplyr::summarise(
    umap_1 = dplyr::first(umap_1),
    umap_2 = dplyr::first(umap_2),
    .groups = "drop"
  )

point_raster <- render_point_raster(
  plot_data = annotated_df,
  class_colors = class_colors,
  x_limits = x_limits,
  y_limits = y_limits
)

grDevices::pdf(plot_path, width = 12, height = 7)

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
    aes(x = umap_1, y = umap_2, color = main_class),
    alpha = 0,
    size = 0,
    show.legend = TRUE
  ) +
  scale_color_manual(values = class_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE) +
  labs(
    title = "UMAP by main_class",
    x = "UMAP1",
    y = "UMAP2",
    color = "main_class"
  ) +
  theme_classic() +
  base_theme

print(p)
grDevices::dev.off()

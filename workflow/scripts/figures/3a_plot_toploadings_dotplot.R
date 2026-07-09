if (exists("snakemake")) {
  loadings_path <- snakemake@input[["loadings"]]
  centering_means_path <- snakemake@input[["centering_means"]]
  pdf_path <- snakemake@output[["pdf"]]
  png_path <- snakemake@output[["png"]]
  factor_name <- snakemake@params[["factor"]]
  top_n <- as.integer(snakemake@params[["top_n"]])
  centering_group <- snakemake@params[["centering_group"]]
  view_order <- unlist(snakemake@params[["views"]])
} else {
  loadings_path <- "results/MCFA/loadings~hvg~5.csv"
  centering_means_path <- "results/MCFA/centering_means~hvg~5.csv"
  pdf_path <- "plots/figures/Figure3/3a_toploadings_dotplot.pdf"
  png_path <- "plots/figures/Figure3/3a_toploadings_dotplot.png"
  factor_name <- "Factor1"
  top_n <- 10L
  centering_group <- "group1"
  view_order <- c("PT", "DCT", "PEC", "TAL", "CNT_PC", "Stromal")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

strip_view_prefix <- function(variable) {
  sub("^[^:]+:", "", as.character(variable))
}

read_table_with_optional_index <- function(path) {
  tbl <- read.csv(path, check.names = FALSE)
  if (names(tbl)[1] %in% c("", "X")) {
    tbl <- tbl[, -1, drop = FALSE]
  }
  tbl
}

require_columns <- function(tbl, required_columns, file_label) {
  missing_columns <- setdiff(required_columns, names(tbl))
  if (length(missing_columns) > 0L) {
    stop(
      file_label,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

loadings <- read_table_with_optional_index(loadings_path)
centering_means <- read_table_with_optional_index(centering_means_path)

require_columns(loadings, c("view", "variable", "Factor", "Loading"), "Loadings file")
require_columns(
  centering_means,
  c("view", "variable", "group", "centering_mean"),
  "Centering means file"
)

factor_loadings <- loadings %>%
  dplyr::filter(.data$Factor == factor_name, .data$view %in% view_order) %>%
  dplyr::mutate(
    view = factor(.data$view, levels = view_order),
    gene = as.character(.data$variable)
  )

if (nrow(factor_loadings) == 0L) {
  stop("No loadings found for ", factor_name, " in selected views.", call. = FALSE)
}

missing_loading_views <- setdiff(view_order, unique(as.character(factor_loadings$view)))
if (length(missing_loading_views) > 0L) {
  stop(
    "Missing selected views in loadings file: ",
    paste(missing_loading_views, collapse = ", "),
    call. = FALSE
  )
}

selected_genes <- factor_loadings %>%
  dplyr::mutate(abs_loading = abs(.data$Loading)) %>%
  dplyr::group_by(.data$view) %>%
  dplyr::slice_max(order_by = .data$abs_loading, n = top_n, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(.data$view, dplyr::desc(.data$Loading), .data$gene) %>%
  dplyr::mutate(
    source_view = as.character(.data$view),
    gene_id = paste(.data$source_view, .data$gene, sep = ":"),
    source_loading = .data$Loading
  ) %>%
  dplyr::select(dplyr::all_of(c("source_view", "gene", "gene_id", "source_loading")))

if (nrow(selected_genes) != length(view_order) * top_n) {
  stop(
    "Expected ",
    length(view_order) * top_n,
    " selected genes, but found ",
    nrow(selected_genes),
    ".",
    call. = FALSE
  )
}

centering_means <- centering_means %>%
  dplyr::filter(.data$group == centering_group, .data$view %in% view_order) %>%
  dplyr::mutate(gene = strip_view_prefix(.data$variable)) %>%
  dplyr::select(dplyr::all_of(c("view", "gene", "centering_mean")))

if (nrow(centering_means) == 0L) {
  stop("No centering means found for group ", centering_group, ".", call. = FALSE)
}

plot_data <- selected_genes %>%
  dplyr::left_join(
    factor_loadings %>%
      dplyr::select(dplyr::all_of(c("view", "gene", "Loading"))),
    by = "gene",
    relationship = "many-to-many"
  ) %>%
  dplyr::filter(.data$view %in% view_order) %>%
  dplyr::left_join(
    centering_means,
    by = c("view", "gene")
  ) %>%
  dplyr::mutate(
    source_view = factor(.data$source_view, levels = view_order),
    view = factor(.data$view, levels = view_order),
    gene_id = factor(.data$gene_id, levels = rev(unique(selected_genes$gene_id)))
  )

missing_centering <- plot_data %>%
  dplyr::filter(is.na(.data$centering_mean)) %>%
  dplyr::distinct(.data$view, .data$gene)

if (nrow(missing_centering) > 0L) {
  stop(
    "Missing centering means for ",
    nrow(missing_centering),
    " selected view/gene combinations.",
    call. = FALSE
  )
}

max_abs_loading <- max(abs(plot_data$Loading), na.rm = TRUE)
gene_labels <- stats::setNames(selected_genes$gene, selected_genes$gene_id)

p_dot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = .data$view, y = .data$gene_id)
) +
  ggplot2::geom_point(
    ggplot2::aes(size = .data$centering_mean, color = .data$Loading),
    alpha = 0.9
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$source_view),
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::scale_y_discrete(labels = gene_labels) +
  ggplot2::scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-max_abs_loading, max_abs_loading),
    name = "Loading"
  ) +
  ggplot2::scale_size_continuous(
    range = c(1, 6),
    name = "Mean log1p"
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    axis.text.y = ggplot2::element_text(size = 7),
    panel.spacing.y = grid::unit(5, "mm"),
    strip.background = ggplot2::element_blank(),
    legend.position = "right"
  )

if (!exists("snakemake")) {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
}

plot_width <- 4
plot_height <- max(10, 0.18 * nrow(selected_genes))

ggplot2::ggsave(
  filename = pdf_path,
  plot = p_dot,
  width = plot_width,
  height = plot_height,
  units = "in"
)

ggplot2::ggsave(
  filename = png_path,
  plot = p_dot,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 300
)

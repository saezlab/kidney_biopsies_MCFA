library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

if (exists("snakemake")) {
  standard_loadings_path <- snakemake@input[["standard_loadings"]]
  exclude_tin_loadings_path <- snakemake@input[["exclude_tin_loadings"]]
  exclude_controls_loadings_path <- snakemake@input[["exclude_controls_loadings"]]
  pdf_path <- snakemake@output[["pdf"]]
  png_path <- snakemake@output[["png"]]
  correlations_csv <- snakemake@output[["correlations_csv"]]
} else {
  standard_loadings_path <- "results/MCFA/loadings~hvg~5.csv"
  exclude_tin_loadings_path <- "results/MCFA/loadings~hvg~5~exclude-TIN.csv"
  exclude_controls_loadings_path <- "results/MCFA/loadings~hvg~5~exclude-controls.csv"
  pdf_path <- "plots/figures/Figure2Ext/ext2_loadings_sensitivity_correlations.pdf"
  png_path <- "plots/figures/Figure2Ext/ext2_loadings_sensitivity_correlations.png"
  correlations_csv <- "results/MCFA/loadings_sensitivity_correlations~hvg~5.csv"
}

read_loadings <- function(path, factor_name, label, loading_column) {
  loadings <- readr::read_csv(path, show_col_types = FALSE)
  required_columns <- c("view", "variable", "Factor", "Loading")
  missing_columns <- setdiff(required_columns, colnames(loadings))
  if (length(missing_columns) > 0) {
    stop(
      "Missing required column(s) in ",
      path,
      ": ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  loadings %>%
    dplyr::filter(.data$Factor == factor_name) %>%
    dplyr::transmute(
      view = as.character(.data$view),
      variable = as.character(.data$variable),
      model = label,
      factor = factor_name,
      !!loading_column := as.numeric(.data$Loading)
    )
}

standard_loadings <- read_loadings(
  standard_loadings_path,
  factor_name = "Factor1",
  label = "Standard hvg~5",
  loading_column = "standard_loading"
)

make_comparison <- function(sensitivity_path, sensitivity_factor, comparison_label) {
  sensitivity_loadings <- read_loadings(
    sensitivity_path,
    factor_name = sensitivity_factor,
    label = comparison_label,
    loading_column = "sensitivity_loading"
  )

  standard_loadings %>%
    dplyr::inner_join(
      sensitivity_loadings,
      by = c("view", "variable"),
      suffix = c("_standard", "_sensitivity")
    ) %>%
    dplyr::mutate(
      comparison = comparison_label,
      standard_factor = "Factor1",
      sensitivity_factor = sensitivity_factor
    )
}

plot_data <- dplyr::bind_rows(
  make_comparison(
    exclude_tin_loadings_path,
    sensitivity_factor = "Factor2",
    comparison_label = "Exclude TIN: Factor2"
  ),
  make_comparison(
    exclude_controls_loadings_path,
    sensitivity_factor = "Factor1",
    comparison_label = "Exclude controls: Factor1"
  )
)

if (nrow(plot_data) == 0L) {
  stop("No overlapping view/variable loadings were found for sensitivity comparisons.", call. = FALSE)
}

correlations <- plot_data %>%
  dplyr::group_by(.data$comparison, .data$view, .data$standard_factor, .data$sensitivity_factor) %>%
  dplyr::summarise(
    n = sum(stats::complete.cases(.data$standard_loading, .data$sensitivity_loading)),
    pearson_r = ifelse(
      n > 1,
      stats::cor(.data$standard_loading, .data$sensitivity_loading, use = "complete.obs", method = "pearson"),
      NA_real_
    ),
    p_value = ifelse(
      n > 2,
      stats::cor.test(.data$standard_loading, .data$sensitivity_loading, method = "pearson")$p.value,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(label = paste0("r = ", sprintf("%.2f", .data$pearson_r), "\nn = ", .data$n))

readr::write_csv(correlations, correlations_csv)

view_order <- c("PT", "TAL", "DCT", "CNT_PC", "IC", "PEC", "POD", "EC", "Immune", "Stromal")
plot_data <- plot_data %>%
  dplyr::mutate(
    view = factor(.data$view, levels = intersect(view_order, unique(.data$view))),
    comparison = factor(.data$comparison, levels = c("Exclude TIN: Factor2", "Exclude controls: Factor1"))
  )
correlations <- correlations %>%
  dplyr::mutate(
    view = factor(.data$view, levels = levels(plot_data$view)),
    comparison = factor(.data$comparison, levels = levels(plot_data$comparison))
  )

p <- ggplot(plot_data, aes(x = .data$standard_loading, y = .data$sensitivity_loading)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.3) +
  geom_point(alpha = 0.28, size = 0.35, color = "#2C7FB8") +
  geom_text(
    data = correlations,
    aes(label = .data$label),
    x = -Inf,
    y = Inf,
    hjust = -0.08,
    vjust = 1.15,
    inherit.aes = FALSE,
    size = 2.7
  ) +
  facet_grid(comparison ~ view) +
  labs(
    x = "Standard hvg~5 Factor1 loading",
    y = "Sensitivity model loading",
    title = "Sensitivity model loadings compared with standard Factor1"
  ) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "grey92", color = "grey75"),
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 8),
    axis.text = element_text(size = 6),
    panel.spacing = unit(0.45, "lines")
  )

ggsave(pdf_path, p, width = 17, height = 6.5, units = "in")
ggsave(png_path, p, width = 17, height = 6.5, units = "in", dpi = 300)

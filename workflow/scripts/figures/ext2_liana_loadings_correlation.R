library(dplyr)
library(ggplot2)
library(readr)

if (exists("snakemake")) {
  standard_loadings_path <- snakemake@input[["standard_loadings"]]
  liana_loadings_path <- snakemake@input[["liana_loadings"]]
  pdf_path <- snakemake@output[["pdf"]]
  png_path <- snakemake@output[["png"]]
  correlations_csv <- snakemake@output[["correlations_csv"]]
} else {
  standard_loadings_path <- "results/MCFA/loadings~hvg~5.csv"
  liana_loadings_path <- "results/MCFA/loadings~liana_hvg~5.csv"
  pdf_path <- "plots/figures/Figure2Ext/ext2_liana_loadings_correlation.pdf"
  png_path <- "plots/figures/Figure2Ext/ext2_liana_loadings_correlation.png"
  correlations_csv <- "results/MCFA/loadings_liana_correlations.csv"
}

read_loadings <- function(path, factor_name, loading_column) {
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
    dplyr::filter(.data$Factor == factor_name, !stringr::str_detect(.data$view, "&")) %>%
    dplyr::transmute(
      view = as.character(.data$view),
      variable = as.character(.data$variable),
      !!loading_column := as.numeric(.data$Loading)
    )
}

standard_loadings <- read_loadings(
  standard_loadings_path,
  factor_name = "Factor1",
  loading_column = "standard_loading"
)

liana_loadings <- read_loadings(
  liana_loadings_path,
  factor_name = "Factor1",
  loading_column = "liana_loading"
)

plot_data <- standard_loadings %>%
  dplyr::inner_join(liana_loadings, by = c("view", "variable"))

if (nrow(plot_data) == 0L) {
  stop("No overlapping standard transcriptome view loadings were found.", call. = FALSE)
}

correlations <- plot_data %>%
  dplyr::group_by(.data$view) %>%
  dplyr::summarise(
    n = sum(stats::complete.cases(.data$standard_loading, .data$liana_loading)),
    pearson_r = ifelse(
      n > 1,
      stats::cor(.data$standard_loading, .data$liana_loading, use = "complete.obs", method = "pearson"),
      NA_real_
    ),
    p_value = ifelse(
      n > 2,
      stats::cor.test(.data$standard_loading, .data$liana_loading, method = "pearson")$p.value,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(label = paste0("r = ", sprintf("%.2f", .data$pearson_r), "\nn = ", .data$n))

readr::write_csv(correlations, correlations_csv)

view_order <- c("PT", "TAL", "DCT", "CNT_PC", "IC", "PEC", "POD", "EC", "Immune", "Stromal")
plot_data <- plot_data %>%
  dplyr::mutate(view = factor(.data$view, levels = intersect(view_order, unique(.data$view))))
correlations <- correlations %>%
  dplyr::mutate(view = factor(.data$view, levels = levels(plot_data$view)))

p <- ggplot(plot_data, aes(x = .data$standard_loading, y = .data$liana_loading)) +
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
  facet_wrap(~ view, nrow = 2) +
  labs(
    x = "Standard hvg~5 Factor1 loading",
    y = "LIANA+ hvg~5 Factor1 loading",
    title = "LIANA+ standard-view loadings compared with standard Factor1"
  ) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "grey92", color = "grey75"),
    strip.text = element_text(size = 8),
    axis.text = element_text(size = 6),
    panel.spacing = unit(0.45, "lines")
  )

ggsave(pdf_path, p, width = 12, height = 6.5, units = "in")
ggsave(png_path, p, width = 12, height = 6.5, units = "in", dpi = 300)

library(dplyr)
library(ggplot2)
library(readr)

if (exists("snakemake")) {
  r2_by_sample_path <- snakemake@input[["r2_by_sample"]]
  pdf_path <- snakemake@output[["pdf"]]
  png_path <- snakemake@output[["png"]]
  csv_path <- snakemake@output[["csv"]]
} else {
  r2_by_sample_path <- "results/MCFA/robustness_r2_by_sample~hvg.csv"
  pdf_path <- "plots/figures/Figure2Ext/ext2a_model_explanatory_power.pdf"
  png_path <- "plots/figures/Figure2Ext/ext2a_model_explanatory_power.png"
  csv_path <- "results/MCFA/model_explanatory_power~hvg.csv"
}

r2_sample <- readr::read_csv(r2_by_sample_path, show_col_types = FALSE)

required_columns <- c("View", "R2", "Factor", "n_factors")
missing_columns <- setdiff(required_columns, colnames(r2_sample))
if (length(missing_columns) > 0) {
  stop(
    "Missing required column(s) in ",
    r2_by_sample_path,
    ": ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

to_plot <- r2_sample %>%
  dplyr::filter(.data$Factor != "All") %>%
  dplyr::filter(!stringr::str_detect(.data$View, "&")) %>%
  dplyr::mutate(
    R2 = pmax(0, as.numeric(.data$R2)),
    n_factors = as.integer(.data$n_factors)
  ) %>%
  dplyr::group_by(.data$View, .data$Factor, .data$n_factors) %>%
  dplyr::summarise(
    quantile_R2 = stats::quantile(.data$R2, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(.data$View, .data$n_factors) %>%
  dplyr::summarise(
    average_95_quantile = mean(.data$quantile_R2, na.rm = TRUE),
    .groups = "drop"
  )

if (nrow(to_plot) == 0L) {
  stop("No transcriptome view R2 values were available to plot.", call. = FALSE)
}

view_order <- c("PT", "TAL", "DCT", "CNT_PC", "IC", "PEC", "POD", "EC", "Immune", "Stromal")
to_plot <- to_plot %>%
  dplyr::mutate(View = factor(.data$View, levels = intersect(view_order, unique(.data$View))))

readr::write_csv(to_plot, csv_path)

p <- ggplot(to_plot, aes(x = .data$n_factors, y = .data$average_95_quantile, color = .data$View)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.4) +
  scale_x_continuous(breaks = sort(unique(to_plot$n_factors))) +
  labs(
    x = "Number of factors in model",
    y = "Model explanatory power*",
    color = "View",
    caption = "*Calculated as the 0.95 quantile of sample-level R2,\nmean averaged over all factors in model."
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25)
  )

ggsave(pdf_path, p, width = 6.5, height = 4.5, units = "in")
ggsave(png_path, p, width = 6.5, height = 4.5, units = "in", dpi = 300)

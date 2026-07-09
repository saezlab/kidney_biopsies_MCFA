if (exists("snakemake")) {
  metadata_csv <- snakemake@input[["metadata"]]
  patient_metadata_csv <- snakemake@input[["patient_metadata"]]
  plot_out <- snakemake@output[["plot"]]
} else {
  metadata_csv <- "results/seurat/metadata.csv"
  patient_metadata_csv <- "data/metadata/patient_metadata_EDA.csv"
  plot_out <- "plots/seurat/exploratory/cluster_5_8_proportions_by_age_egfr.pdf"
}

library(tidyverse)

# Read metadata
metadata_df <- read_csv(metadata_csv, show_col_types = FALSE)

# Read patient metadata
patient_metadata_df <- read_csv(patient_metadata_csv, show_col_types = FALSE)

# Get the sample column (assuming it's "sample" or "orig.ident")
sample_col <- "sample"
if (!sample_col %in% colnames(metadata_df)) {
  sample_col <- "orig.ident"
}

# Calculate proportions of all clusters at resolution SCT_snn_res.0.2
proportions_df <- metadata_df %>%
  dplyr::select(all_of(sample_col), SCT_snn_res.0.2) %>%
  dplyr::rename(sample = all_of(sample_col)) %>%
  dplyr::group_by(sample) %>%
  dplyr::count(SCT_snn_res.0.2, name = "count") %>%
  dplyr::mutate(total = sum(count)) %>%
  dplyr::mutate(proportion = count / total) %>%
  dplyr::ungroup() %>%
  dplyr::select(sample, cluster = SCT_snn_res.0.2, proportion)

# Apply CLR transformation
# CLR(p_i) = log(p_i / geometric_mean(all p))
clr_df <- proportions_df %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(
    log_proportion = log(proportion),
    geometric_mean = exp(mean(log_proportion, na.rm = TRUE)),
    clr_value = log(proportion / geometric_mean)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(sample, cluster, clr_value)

# Extract clusters 5 and 8
clr_5_8 <- clr_df %>%
  dplyr::filter(cluster %in% c(5, 8)) %>%
  tidyr::pivot_wider(
    names_from = cluster,
    values_from = clr_value,
    names_prefix = "clr_cluster_"
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("clr_cluster_"),
    names_to = "cluster",
    names_prefix = "clr_cluster_",
    values_to = "clr_proportion"
  ) %>%
  dplyr::mutate(cluster = as.numeric(cluster))

# Merge with patient metadata
# Match by sample_ID (K1, K2, etc.)
patient_metadata_df <- patient_metadata_df %>%
  dplyr::rename_all(tolower)

plot_data <- clr_5_8 %>%
  dplyr::left_join(
    patient_metadata_df %>%
      dplyr::select(sample_id, date_of_biopsy, egfr_at_time_of_biopsy),
    by = c("sample" = "sample_id")
  ) %>%
  dplyr::select(sample, cluster, clr_proportion, date_of_biopsy, egfr_at_time_of_biopsy)

# Calculate biopsy age (from 2021-12-31)
cutoff_date <- as.Date("2021-12-31")
plot_data <- plot_data %>%
  dplyr::mutate(
    date_of_biopsy = as.Date(date_of_biopsy),
    biopsy_age_years = as.numeric(cutoff_date - date_of_biopsy) / 365.25
  )

# Remove rows with missing age or eGFR
plot_data_clean <- plot_data %>%
  dplyr::filter(!is.na(date_of_biopsy) & !is.na(egfr_at_time_of_biopsy))

# Create 2x2 subplot figure
library(ggplot2)
library(patchwork)

format_lm_stats <- function(data, x_col) {
  model_data <- data %>%
    dplyr::select(clr_proportion, x = all_of(x_col)) %>%
    tidyr::drop_na()

  if (nrow(model_data) < 3 || dplyr::n_distinct(model_data$x) < 2) {
    return("R2 = NA\np = NA")
  }

  fit <- lm(clr_proportion ~ x, data = model_data)
  fit_summary <- summary(fit)
  r_squared <- fit_summary$r.squared
  p_value <- fit_summary$coefficients["x", "Pr(>|t|)"]

  paste0(
    "R2 = ", formatC(r_squared, format = "f", digits = 3),
    "\np = ", format.pval(p_value, digits = 3, eps = 0.001)
  )
}

cluster_5_data <- plot_data_clean %>%
  dplyr::filter(cluster == 5)

cluster_8_data <- plot_data_clean %>%
  dplyr::filter(cluster == 8)

# Cluster 5 vs Age
p1 <- cluster_5_data %>%
  ggplot(aes(x = biopsy_age_years, y = clr_proportion)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, color = "blue", fill = "lightblue", alpha = 0.3) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = format_lm_stats(cluster_5_data, "biopsy_age_years"),
    hjust = 1.1, vjust = 1.2,
    size = 3
  ) +
  labs(
    title = "Cluster 5: CLR-Proportion vs. Biopsy Age",
    x = "Biopsy Age (years from 2021-12-31)",
    y = "CLR-Transformed Proportion"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, hjust = 0.5))

# Cluster 5 vs eGFR
p2 <- cluster_5_data %>%
  ggplot(aes(x = egfr_at_time_of_biopsy, y = clr_proportion)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, color = "blue", fill = "lightblue", alpha = 0.3) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = format_lm_stats(cluster_5_data, "egfr_at_time_of_biopsy"),
    hjust = 1.1, vjust = 1.2,
    size = 3
  ) +
  labs(
    title = "Cluster 5: CLR-Proportion vs. eGFR",
    x = "eGFR (mL/min/1.73m²)",
    y = "CLR-Transformed Proportion"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, hjust = 0.5))

# Cluster 8 vs Age
p3 <- cluster_8_data %>%
  ggplot(aes(x = biopsy_age_years, y = clr_proportion)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, color = "red", fill = "lightcoral", alpha = 0.3) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = format_lm_stats(cluster_8_data, "biopsy_age_years"),
    hjust = 1.1, vjust = 1.2,
    size = 3
  ) +
  labs(
    title = "Cluster 8: CLR-Proportion vs. Biopsy Age",
    x = "Biopsy Age (years from 2021-12-31)",
    y = "CLR-Transformed Proportion"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, hjust = 0.5))

# Cluster 8 vs eGFR
p4 <- cluster_8_data %>%
  ggplot(aes(x = egfr_at_time_of_biopsy, y = clr_proportion)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, color = "red", fill = "lightcoral", alpha = 0.3) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = format_lm_stats(cluster_8_data, "egfr_at_time_of_biopsy"),
    hjust = 1.1, vjust = 1.2,
    size = 3
  ) +
  labs(
    title = "Cluster 8: CLR-Proportion vs. eGFR",
    x = "eGFR (mL/min/1.73m²)",
    y = "CLR-Transformed Proportion"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, hjust = 0.5))

# Combine plots into 2x2 grid
combined_plot <- (p1 + p2) / (p3 + p4) +
  plot_layout(ncol = 2, nrow = 2)

# Save plot
ggsave(plot_out, combined_plot, width = 14, height = 10, dpi = 300)

cat("Plot saved to:", plot_out, "\n")
cat("Number of samples plotted:", n_distinct(plot_data_clean$sample), "\n")

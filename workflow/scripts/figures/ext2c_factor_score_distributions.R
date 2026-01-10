library(magrittr)
library(ggplot2)
library(ggpubr)
library(patchwork)

if (exists("snakemake")){
  associations <- read.csv(snakemake@input$associations)
  plt1 <- snakemake@output$plt1
  plt1png <- snakemake@output$plt1png
} else {
  factors_obs <- read.csv('results/MCFA/factors_obs~all~5.csv', row.names = 1)
}

colnames(factors_obs) <- stringr::str_replace(colnames(factors_obs), "X[:digit:]+_", "")

factors_obs <- factors_obs %>%
  tibble::rownames_to_column("sample") %>%
  mutate(
    disease = ifelse(stringr::str_detect(sample, "KB"), "control_KPMP", disease),
    Group = disease
  ) %>%
  # pivot longer to have "Factor" and "Value" columns
  tidyr::pivot_longer(
    cols = starts_with("Factor"),
    names_to = "Factor",
    values_to = "Value"
  )

# Single ggplot with facets
ggplot(factors_obs, aes(
  y = forcats::fct_rev(factor(
    Group,
    levels = c("DN", "FSGS", "IgA", "MN", "TIN", "control", "control_KPMP"),
    labels = c("DN", "FSGS", "IgAN", "MN", "TIN", "control", "control_KPMP")
  )),
  x = Value,
  color = Group
)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.7) +
  scale_color_brewer(palette = "Set2") +
  facet_wrap(~ Factor, scales = "free_x", nrow=1) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "Condition", x = "Value")

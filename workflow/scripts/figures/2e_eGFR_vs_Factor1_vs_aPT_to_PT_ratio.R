library(magrittr)
library(dplyr)
library(ggplot2)
library(patchwork)
data <- read.csv('results/composition/prop~main_class~atlas.csv')
factors_obs <- read.csv('results/MCFA/factors_obs~hvg~5.csv')

colnames(data) <- stringr::str_replace(colnames(data), "X[:digit:]+_", "")
df <- data %>%
  dplyr::left_join(factors_obs, by = 'Sample') %>%
  dplyr::transmute(Sample,
                   `Log ratio aPT to PT` = log(aPT_aTAL_aDCT/PT),
                   Group = factor(disease,
                                  levels = c("DN", "FSGS", "IgA", "MN", "TIN", "control", "control_KPMP"),
                                  labels = c("DN", "FSGS", "IgAN", "MN", "TIN", "control", "control_KPMP")),
                   Factor1,
                   `Baseline eGFR` = egfr_at_time_of_biopsy) %>%
  tidyr::drop_na(Factor1) # Restrict to samples included in model

p1 <- df %>%
  ggplot(aes(y = `Baseline eGFR`, x = `Log ratio aPT to PT`)) +
  geom_point(aes(color = Group)) + theme_classic() + scale_colour_brewer(palette = "Set2") +
  geom_smooth(method = "lm", se = FALSE, color = "black", alpha = 0.8) +
  ggpubr::stat_cor(method = "pearson", label.y = 160)


p2 <- df %>%
  ggplot(aes(y = `Baseline eGFR`, x = Factor1)) +
  geom_point(aes(color = Group)) + theme_classic() + scale_colour_brewer(palette = "Set2")  +
  geom_smooth(method = "lm", se = FALSE, color = "black", alpha = 0.8) +
  ggpubr::stat_cor(method = "pearson", label.y = 160)

p3 <- df %>%
  ggplot(aes(y = `Log ratio aPT to PT`, x = Factor1)) +
  geom_point(aes(color = Group)) + theme_classic() + scale_colour_brewer(palette = "Set2")  +
  geom_smooth(method = "lm", se = FALSE, color = "black", alpha = 0.8) +
  ggpubr::stat_cor(method = "pearson", label.y = 6)

design <-
"
21
43
"
patchwork::wrap_plots(p1, p2, guide_area(), p3, design = design) +
  plot_layout(guides = "collect")

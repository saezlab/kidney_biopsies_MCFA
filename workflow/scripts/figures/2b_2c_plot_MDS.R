library(magrittr)
library(ggplot2)
library(ggpubr)
library(patchwork)

factors_obs <- read.csv('results/MCFA/factors_obs~all~5.csv', row.names = 1)

colnames(factors_obs) <- stringr::str_replace(colnames(factors_obs), "X[:digit:]+_", "")

factors_obs <- factors_obs %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::mutate(disease = ifelse(stringr::str_detect(sample, "KB"), "control_KPMP", disease)) %>%
  dplyr::rename(Condition = disease)


mat <- factors_obs %>%
  dplyr::select(sample, Factor1, Factor2, Factor3, Factor5) %>% # Factor 4 outlier
  tibble::column_to_rownames('sample') %>%
  as.matrix()
data_matrix <- scale(mat)
dist_matrix <- dist(data_matrix, method = "euclidean")
mds_result <- cmdscale(dist_matrix, k = 2)  # k=2 for 2D representation
mds_df <- as.data.frame(mds_result) %>%
  setNames(c("Dim1", "Dim2")) %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::left_join(factors_obs)

p1 <- ggplot(mds_df, aes(x = Dim1, y = Dim2)) +
  geom_point(aes(color = factor(Condition, levels = c("DN", "FSGS", "IgA", "MN", "TIN", "control", "control_KPMP"),
                                labels = c("DN", "FSGS", "IgAN", "MN", "TIN", "control", "control_KPMP")))) +
  theme_classic() + scale_colour_brewer(palette = "Set2") + labs(color = "Condition") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# No eGFR data available for KPMP reference samples
p2 <- ggplot(mds_df, aes(x = Dim1, y = Dim2)) +
  geom_point(aes(color = egfr_at_time_of_biopsy)) +
  theme_classic() + scale_color_viridis_c() + labs(color = "Baseline eGFR\n(ml/min/1.73m2)") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

p1 + p2

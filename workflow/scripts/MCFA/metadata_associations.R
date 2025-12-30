library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(stats)

df <- read.csv("results/MCFA/factors_obs~adata~LIANA_all~all~5.csv")
df$disease[stringr::str_detect(df$Sample, "KB")] <- "control_KPMP"

# Make empty strings explicitly NAs, and one-hot encode "disease" (Group) column
df <- df %>%
  mutate(across(where(is.character), ~na_if(.x, ""))) %>%
  mutate(disease = ifelse(disease == "IgA", "IgAN", disease)) %>%
  mutate(value = TRUE)  %>% spread(disease, value,  fill = FALSE )

colnames(df) <- stringr::str_remove(colnames(df), "X[:digit:]+_")

# Define which columns are factor scores and which are metadata
factor_columns <- c("Factor1", "Factor2", "Factor3", "Factor5") # Factor 4 - outlier
metadata_columns <- c("control", "control_KPMP", "DN", "FSGS", "TIN", "IgAN", "MN",
                      "age_at_biopsy", "sex", "IFTA", "prop_gs",
                      "ace_inhibition", "angiotensin_ii_inhibition",
                      "statin", "corticosteroids",  "dialysis",
                      "transplant", "egfr_at_time_of_biopsy", "proteinuria_around_biopsy",
                      "jm_slope", "classification")# "total_slope", "immediate_slope", "chronic_slope" removed because we don't use these versions

# Separate continuous and categorical metadata
continuous_metadata <- metadata_columns[sapply(df[metadata_columns], is.numeric)]
categorical_metadata <- metadata_columns[sapply(df[metadata_columns], is.factor) | sapply(df[metadata_columns], is.character) | sapply(df[metadata_columns], is.logical)]
clr_metadata <- colnames(df)[stringr::str_detect(colnames(df), "_clr")]

# Pearson correlation for continuous metadata
corr_results <- data.frame(Factor=character(), Metadata=character(), Correlation=numeric(), P_value=numeric())

for (factor in factor_columns) {
  for (meta in continuous_metadata) {
    if (sum(!is.na(df[[factor]])) > 1 & sum(!is.na(df[[meta]])) > 1) {
      test <- cor.test(df[[factor]], df[[meta]], use="complete.obs", method="pearson")
      corr_results <- rbind(corr_results, data.frame(Factor=factor, Metadata=meta, Correlation=test$estimate, P_value=test$p.value))
    }
  }
}
# Apply FDR correction for Pearson correlation p-values
corr_results$Adjusted_P_value <- p.adjust(corr_results$P_value, method="fdr")
corr_results$annotation <- ifelse(corr_results$Adjusted_P_value < 0.05, "∗", "")

# Pearson correlation for continuous metadata - CLR
clr_results <- data.frame(Factor=character(), Metadata=character(), Correlation=numeric(), P_value=numeric())

for (factor in factor_columns) {
  for (meta in clr_metadata) {
    if (sum(!is.na(df[[factor]])) > 1 & sum(!is.na(df[[meta]])) > 1) {
      test <- cor.test(df[[factor]], df[[meta]], use="complete.obs", method="pearson")
      clr_results <- rbind(clr_results, data.frame(Factor=factor, Metadata=meta, Correlation=test$estimate, P_value=test$p.value))
    }
  }
}
# Apply FDR correction for Pearson correlation p-values
clr_results$Adjusted_P_value <- p.adjust(clr_results$P_value, method="fdr")
clr_results$annotation <- ifelse(clr_results$Adjusted_P_value < 0.05, "∗", "")

# ANOVA test for categorical metadata
anova_results <- data.frame(Factor=character(), Metadata=character(), F_statistic=numeric(), P_value=numeric(), Degrees_freedom=numeric())

for (factor in factor_columns) {
  for (meta in categorical_metadata) {
    print(meta)
    print(factor)
    if (length(unique(df[[meta]])) > 1) { # Ensure multiple groups exist
      model <- aov(df[[factor]] ~ as.factor(df[[meta]]), data=df, na.action = na.exclude)
      test <- summary(model)
      p_value <- test[[1]][["Pr(>F)"]][1]
      f_statistic <- test[[1]][["F value"]][1]
      degrees_freedom <- test[[1]][['Df']][1]
      anova_results <- rbind(anova_results, data.frame(Factor=factor, Metadata=meta, F_statistic=f_statistic, P_value=p_value, Degrees_freedom = degrees_freedom))
    }
  }
}
# Apply FDR correction for Chi-square test p-values
anova_results$Adjusted_P_value <- p.adjust(anova_results$P_value, method="fdr")
anova_results$annotation <- ifelse(anova_results$Adjusted_P_value < 0.05, "∗", "")

write.csv(corr_results, "results/MCFA/metadata_pearson_correlations.csv")
write.csv(anova_results, "results/MCFA/metadata_anova.csv")
write.csv(clr_results, "results/MCFA/clr_pearson_correlations.csv")

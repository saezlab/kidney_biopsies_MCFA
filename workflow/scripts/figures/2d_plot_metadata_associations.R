library(ggplot2)
library(dplyr)
library(magrittr)

corr_results <- read.csv("results/MCFA/metadata_pearson_correlations.csv")
anova_results <- read.csv("results/MCFA/metadata_anova.csv")
clr_results <- read.csv("results/MCFA/clr_pearson_correlations.csv")

# Plot Pearson Correlation Heatmap
p1 <- ggplot(corr_results, aes(y=forcats::fct_rev(factor(Metadata,
                                        levels = c("age_at_biopsy",
                                                   "egfr_at_time_of_biopsy",
                                                   "IFTA",
                                                   "prop_gs",
                                                   "proteinuria_around_biopsy",
                                                   "jm_slope"),
                                        labels = c("Age",
                                                   "Baseline eGFR",
                                                   "IFTA",
                                                   "Glomerulosclerosis",
                                                   "PCR",
                                                   "eGFR slope"))), x=Factor, fill=Correlation)) +
  geom_tile() +
  geom_text(aes(label=annotation), color="black", size = 5, family = "Courier") +
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, limits=c(-0.8, 0.8)) +
  labs(x="Metadata", y="Metadata (continuous)", fill="Pearson Correlation") +
  theme_minimal() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Plot anova Test Results Heatmap
p2 <- ggplot(anova_results, aes(y=forcats::fct_rev(factor(Metadata,
                                         levels = c("sex", 
                                                    "DN", "FSGS", "IgAN", "MN", "TIN", "control", "control_KPMP",
                                                    "dialysis", "transplant",
                                                    "ace_inhibition", "angiotensin_ii_inhibition", "corticosteroids", "statin",
                                                    "classification"),
                                         labels = c("Sex",
                                                    "DN", "FSGS", "IgAN", "MN", "TIN", "control", "control_KPMP",
                                                    "Dialysis", "Transplant",
                                                    "ACE inhibition", "Angiotensin II inhibition", "Corticosteroids", "Statin",
                                                    "Disease progression"))),
                                x=Factor, 
                                fill=F_statistic)) +
  geom_tile() +
  geom_text(aes(label=annotation), color="black", size = 5, family = "Courier") +
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0) +
  labs(x="Factor", y="Metadata (categorical)", fill="F statistic") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p3 <- clr_results %>%
  dplyr::mutate(Metadata = stringr::str_remove(Metadata, "_clr")) %>%
  dplyr::mutate(Metadata = ifelse(Metadata == "aPT",
                       "aPT2",
                       ifelse(Metadata == "aPT_aTAL_aDCT",
                              "aPT", Metadata))) %>%
  ggplot(aes(y=forcats::fct_rev(Metadata), x=Factor, fill=Correlation)) +
  geom_tile() +
  geom_text(aes(label=annotation), color="black", size = 5, family = "Courier") +
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, limits=c(-0.8, 0.8)) +
  labs(x="Factor", y="Cell type abundance", fill="Pearson Correlation", caption="* Adjusted p-value < 0.05") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

design <- "
13
23
23
"

patchwork::wrap_plots(p1, p2, p3, design=design) + plot_layout(guides = "collect")

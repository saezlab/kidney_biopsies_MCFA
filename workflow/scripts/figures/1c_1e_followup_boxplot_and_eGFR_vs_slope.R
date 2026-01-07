library(magrittr)
library(dplyr)
library(ggplot2)
library(patchwork)

if (exists("snakemake")){
  data <- read.csv(snakemake@input$data)
  plt <- snakemake@output$plt
  gg <- snakemake@output$gg
} else{
  data <- read.csv('data/data_overview/follow_up~boxplot.csv')
  df <- read.csv("data/metadata/patient_metadata_EDA.csv")
}

(
p1 <- ggplot(data, aes(y = forcats::fct_rev(disease), x = days_from_biopsy/365.25, group = disease)) +
  geom_boxplot(outlier.shape = NA, width=0.5) +
  geom_jitter(aes(group = interaction(disease, patient_ID)), size = 0.5, height = 0.2) +
  labs(y = element_blank(), x = "Follow-up (y)") +
  theme_classic() +
  scale_colour_brewer(palette = "Set2") 
)


(
  p2 <- df %>%
  dplyr::mutate(disease = ifelse(disease == "IgA", "IgAN", disease)) %>%
  dplyr::filter(disease %in% c("MN", "TIN", "DN", "IgAN", "FSGS")) %>%
  ggplot(aes(x = egfr_at_time_of_biopsy, y = jm_slope)) + 
  geom_point(aes(color = disease)) + theme_classic() +
  scale_color_brewer(palette = "Set2") +
  geom_hline(yintercept = 0, alpha = 0.3) +
  labs(x = "eGFR at biopsy (ml/min/1.73m2)",
       y = "eGFR slope (ml/min/1.73m2/y)",
       color = "Group")
)

design <- "
11222"
patchwork::wrap_plots(p1, p2, design=design)

#ggsave(plot = p1, filename = plt, device = "pdf", width = 20, height = 14, units = "cm")
#saveRDS(p1, file = gg)

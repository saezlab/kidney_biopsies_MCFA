library(magrittr)
library(dplyr)
library(ggplot2)

if (exists("snakemake")){
  data <- read.csv(snakemake@input$data)
  plt <- snakemake@output$plt
  gg <- snakemake@output$gg
} else{
  data <- read.csv('results/composition/clr~main_class~atlas.csv')
  patient_metadata <- read.csv('data/metadata/patient_metadata_EDA.csv')
} 

metadata <- patient_metadata %>%
  dplyr::select(sample_ID, disease, ) %>% 
  dplyr::mutate(disease = ifelse(disease == "IgA", "IgAN", disease)) %>%
  dplyr::distinct()
metadata_ext_controls <- data.frame(sample_ID = unique(data$Sample[stringr::str_detect(data$Sample, "KB")]),
                                    disease = "control_KPMP")
metadata <- dplyr::bind_rows(metadata, metadata_ext_controls)

colnames(data) <- stringr::str_replace(colnames(data), "X[:digit:]+_", "")

to_plot1 <- data %>%
  tidyr::pivot_longer(!Sample, names_to = 'group', values_to = 'prop') %>%
  dplyr::left_join(metadata, by = c("Sample" = "sample_ID")) %>%
  dplyr::filter(!disease == "control_AKI") %>%
  dplyr::mutate(group = ifelse(group == "aPT", "aPT1",
                               ifelse(group == "aPT_aTAL_aDCT", "aPT", group)))

(
  p2 <- to_plot1 %>%
    dplyr::filter(group %in% c("aPT", "PT", "Immune")) %>%
    ggplot(aes(x = disease, y = prop, color = disease)) +
    geom_hline(yintercept = 0, colour = "lightgray") +
    facet_wrap(.~factor(group, levels = c("aPT", "PT", "Immune")), nrow=1) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5) +
    theme_classic() + scale_colour_brewer(palette = "Set2") + 
    labs(x = element_blank(), y = 'CLR-transformed proportions', subtitle = "Cell type abundance by condition") +
    theme(legend.position="none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggpubr::stat_compare_means(method = "t.test",
                               comparisons = list(c("control_KPMP", "control"),
                                                  c("control_KPMP", "DN"),
                                                  c("control_KPMP", "FSGS"),
                                                  c("control_KPMP", "IgAN"),
                                                  c("control_KPMP", "MN"),
                                                  c("control_KPMP", "TIN")), size = 3) +
    ylim(c(-2.5, 9.5))
)

metadata <- patient_metadata %>%
  dplyr::select(sample_ID, disease, egfr_at_time_of_biopsy, IFTA, proteinuria_around_biopsy, classification) %>% 
  dplyr::mutate(disease = ifelse(disease == "IgA", "IgAN", disease)) %>%
  dplyr::distinct()

to_plot <- data %>%
  tidyr::pivot_longer(!Sample, names_to = 'group', values_to = 'prop') %>%
  dplyr::left_join(metadata, by = c("Sample" = "sample_ID")) %>%
  dplyr::filter(!disease == "control_AKI") %>%
  dplyr::mutate(eGFR = ifelse(egfr_at_time_of_biopsy < 30,
                              "< 30",
                              ifelse(egfr_at_time_of_biopsy < 60,
                                     "30-60",
                                     ifelse(egfr_at_time_of_biopsy < 90,
                                            "60-90",
                                            ifelse(egfr_at_time_of_biopsy >= 90,
                                                   ">= 90",
                                                   NA)
                                     )
                              )
  )
  ) %>%
  dplyr::mutate(eGFR = factor(eGFR, levels = c(">= 90", "60-90", "30-60", "< 30"))) %>%
  dplyr::mutate(IFTA_levels = ifelse(IFTA < 0.1,
                                     "< 10%",
                                     ifelse(IFTA < 0.5,
                                            "10-50%",
                                            ifelse(IFTA >= 0.5,
                                                   ">= 50%",
                                                   NA)))) %>%
  dplyr::mutate(IFTA_levels = factor(IFTA_levels, levels = c("< 10%", "10-50%", ">= 50%"))) %>%
  dplyr::mutate(PCR = ifelse(proteinuria_around_biopsy < 1,
                             "< 1 g/g",
                             ifelse(proteinuria_around_biopsy < 3,
                                    "1-3 g/g",
                                    ifelse(proteinuria_around_biopsy >= 3,
                                           ">= 3 g/g",
                                           NA)))) %>%
  dplyr::mutate(PCR = factor(PCR, levels = c("< 1 g/g", "1-3 g/g", ">= 3 g/g"))) %>%
  dplyr::mutate(group = ifelse(group == "aPT", "aPT1",
                               ifelse(group == "aPT_aTAL_aDCT", "aPT", group)))

(
  p3 <- to_plot %>%
    dplyr::filter(group %in% c("aPT", "PT", "Immune")) %>%
    tidyr::drop_na(eGFR) %>%
    ggplot(aes(x = eGFR,
                            y = prop, color = eGFR)) +
    geom_hline(yintercept = 0, colour = "lightgray") +
    facet_wrap(.~factor(group, levels = c("aPT", "PT", "Immune")), nrow=1) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5) +
    theme_classic() + scale_colour_brewer(palette = "Set1") + 
    labs(x = element_blank(), y = 'CLR-transformed proportions', subtitle = "Cell type abundance by eGFR at biopsy (ml/min/1.73m2)") +
    theme(legend.position="none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggpubr::stat_compare_means(method = "t.test",
                               comparisons = list(c(">= 90", "60-90"),
                                                  c(">= 90", "30-60"),
                                                  c(">= 90", "< 30")), size = 3)  +
    ylim(c(-2.5, 9.5))
)

(
  p4 <- to_plot %>%
    dplyr::filter(group %in% c("aPT", "PT", "Immune")) %>%
    tidyr::drop_na(IFTA_levels) %>%
    ggplot(aes(x = IFTA_levels,
                            y = prop, color = IFTA_levels)) +
    geom_hline(yintercept = 0, colour = "lightgray") +
    facet_wrap(.~factor(group, levels = c("aPT", "PT", "Immune")), nrow=1) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5) +
    theme_classic() + scale_colour_brewer(palette = "Set1") + 
    labs(x = element_blank(), y = 'CLR-transformed proportions', subtitle = "Cell type abundance by IF/TA") +
    theme(legend.position="none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggpubr::stat_compare_means(method = "t.test",
                               comparisons = list(c("< 10%", "10-50%"),
                                                  c("< 10%", ">= 50%")), size = 3)  +
    ylim(c(-2.5, 9.5))
)

library(patchwork)
p2+p3+p4

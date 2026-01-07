library(magrittr)
library(dplyr)
library(ggplot2)
library(patchwork)

if (exists("snakemake")){
  data <- read.csv(snakemake@input$data)
  plt <- snakemake@output$plt
  gg <- snakemake@output$gg
} else{
  data <- read.csv('data/data_overview/age~boxplot.csv')
}

(
  p1 <- ggplot(data, aes(y = forcats::fct_rev(factor(disease,
                                    levels = c("DN", "FSGS", "IgAN",
                                               "MN", "TIN", "AKI (control)",
                                               "control"))),
                         x = age_at_biopsy,
                         group = disease)) +
    geom_boxplot(outlier.shape = NA, width=0.5) +
    geom_jitter(aes(group = interaction(disease, patient_ID)), size = 0.5, height = 0.2) +
    labs(y = element_blank(), x = "Age (y)") +
    theme_classic() +
    scale_colour_brewer(palette = "Set2") 
)

  
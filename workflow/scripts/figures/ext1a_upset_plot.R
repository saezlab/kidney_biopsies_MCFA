library(magrittr)
library(dplyr)
library(ggplot2)
library(UpSetR)
if (exists("snakemake")){
  data <- read.csv(snakemake@input$data)
  plt <- snakemake@output$plot
  gg <- snakemake@output$gg
} else{
  data <- read.csv('results/data_overview/modality~heatmap.csv')
}

# Stacked bar plot
df <- data %>%
  dplyr::group_by(modality, disease) %>%
  dplyr::summarise(Count = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(disease, modality, fill = list(Count = 0)) %>%
  dplyr::mutate(modality = factor(modality, levels = c("RNA", "ATAC", "Xenium", "RNA, ATAC", "RNA, Xenium", "ATAC, Xenium", "RNA, ATAC, Xenium")))


df2 <- df %>%
  dplyr::mutate(modality = stringr::str_replace_all(modality, ", ", "&"))
upset(fromExpression(df2$Count %>% setNames(df2$modality)), order.by = "freq", mainbar.y.max = 80)

# saveRDS(p1, file = gg)
# ggsave(plot = p1, filename = plt, device = "pdf", width = 13, height = 6.2, units = "cm")

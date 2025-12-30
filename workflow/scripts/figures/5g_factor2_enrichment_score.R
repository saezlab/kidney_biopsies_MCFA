library(magrittr)
library(ggplot2)

df <- read.csv('results/RelA_KO/factor2_enrichment.csv', row.names = 1) %>%
  tibble::rownames_to_column("replicate") %>%
  dplyr::mutate(condition = ifelse(stringr::str_detect(replicate, "_AA_"), "NT", "KO"))

df %>% ggplot(aes(x = condition, y = act)) + geom_boxplot() + geom_jitter(height = 0) +
  ggpubr::stat_compare_means(comparisons = list(c("KO", "NT")), method = "t.test") + theme_classic() + theme(axis.title.x = element_blank()) + 
                                                        ylab("Factor 2 PT Enrichment Score") + ylim(c(-11,11)) + scale_x_discrete(limits = c("NT", "KO"))



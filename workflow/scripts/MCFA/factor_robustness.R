library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(Hmisc)

if (exists("snakemake")){
  scores <- read.csv(snakemake@input$scores)[,-c(1)]
  loadings <- read.csv(snakemake@input$loadings)[,-c(1)]
  factor_robustness <- snakemake@output$factor_robustness
} else {
  scores <- read.csv('results/MCFA/robustness_scores~all.csv')[,-c(1)]
  loadings <- read.csv('results/MCFA/robustness_loadings~all.csv')[,-c(1)]
  factor_robustness <- 'results/MCFA/factor_robustness.csv'
}

#####
# Correlation between factor scores
#####

pivoted_scores <- scores %>%
  #dplyr::filter(n_factors >=20 & n_factors <=40) %>%
  tidyr::pivot_wider(names_from = c(n_factors, Factor),
                     values_from = Score,
                     values_fill = NULL)

mat_scores <- pivoted_scores %>%
  dplyr::select(-c(Sample)) %>%
  as.matrix()

corr_scores <- Hmisc::rcorr(mat_scores, type="pearson")$r

#####
# Correlation between factor loadings
#####


pivoted_loadings <- loadings %>%
  tidyr::pivot_wider(names_from = c(n_factors, Factor),
                     values_from = Loading,
                     values_fill = NULL)

mat_loadings <- pivoted_loadings %>%
  dplyr::select(-c(view, variable)) %>%
  as.matrix()

corr_loadings <- Hmisc::rcorr(mat_loadings, type="pearson")$r


####
# Product of score correlation and loading correlation defines factor similarity
####
mat <- corr_scores*corr_loadings

write.csv(mat, factor_robustness)
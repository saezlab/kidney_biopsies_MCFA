library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(Hmisc)

mat <- read.csv('results/MCFA/factor_robustness.csv')
r2_sample <- read.csv("results/MCFA/robustness_r2_by_sample~all.csv") %>%
  dplyr::rename('Sample' = 'Group')


# Define a factor as equivalent if the correlation between the scores 
# multiplied by the correlation between the loadings is creater than 0.9^2
add_stars_make_diag <- function(j, i, x, y, w, h, fill){
  if(i <= j) {
    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
    if (mat[i, j] >= 0.9^2) {
      grid.text("*", x, y)
    }
  }
}

add_stars <- function(j, i, x, y, w, h, fill){
  if (mat[i, j] >= 0.9^2) {
    grid.text("*", x, y)
  }
}

f1 <- circlize::colorRamp2(seq(-1, 1, length = 3),
                           c("blue", "#EEEEEE", "red"))

colnames(mat) <- NULL

mat_subset <- mat[c('5_Factor1', '5_Factor2', '5_Factor3', '5_Factor4', '5_Factor5'),]
add_stars_make_diag2 <- function(j, i, x, y, w, h, fill){
  if(i <= j) {
    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
    if (mat_subset[i, j] >= 0.9^2) {
      grid.text("*", x, y)
    }
  }
}

Heatmap(mat_subset,
        rect_gp = gpar(type = "none"), 
        cell_fun = add_stars_make_diag2,
        row_names_gp = gpar(fontsize = 8),
        cluster_rows = F,
        cluster_row_slices = F,
        cluster_column_slices = F,
        row_title = NULL,
        column_split = rep(c(3, 4, 5, 10, 15, 20, 25, 30), times = c(3, 4, 5, 10, 15, 20, 25, 30)),
        row_names_side = "right",
        cluster_columns = F,
        col = f1,
        heatmap_legend_param = list(title = "r"))

# Take the 95th quantile sample-level R2 as my summary statistic for factor
# explanatory power: a useful factor should have explanatory power for more thqn 5% of samples.
# And define the model explanatory power as the
# average 95th quantile sample-level R2 across all factors.
to_plot <- r2_sample %>% dplyr::filter(!Factor == "All") %>%
  dplyr::filter(!stringr::str_detect(View, "&")) %>%
  dplyr::mutate(R2 = pmax(0, R2)) %>% # set minimum R2 to 0
  dplyr::group_by(View, Factor, n_factors) %>%
  dplyr::summarise(max_R2 = max(R2, na.rm = TRUE),
                   quantile_R2 = quantile(R2, 0.95, na.rm=TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(View, n_factors) %>%
  dplyr::summarise(average_95_quantile = mean(quantile_R2, na.rm = TRUE))

ggplot(to_plot, aes(x = n_factors, y = average_95_quantile, color = factor(View))) +
  geom_line() +
  labs(x = "Number of factors in model",
       y = "Model explanatory power*",
       color = "View",
       caption = "*Calculated as the 0.95 quantile of sample-level R2,\nmean averaged over all factors in model.") +  theme_minimal()


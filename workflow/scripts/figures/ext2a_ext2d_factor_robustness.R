library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(Hmisc)

scores <- read.csv('results/MCFA/robustness_scores~adata~LIANA_all~all.csv')[,-c(1)]
loadings <- read.csv('results/MCFA/robustness_loadings~adata~LIANA_all~all.csv')[,-c(1)]
r2_sample <- read.csv("results/MCFA/robustness_r2_by_sample~adata~LIANA_all~all.csv") %>%
  dplyr::rename('Sample' = 'Group')
r2_total <- read.csv("results/MCFA/robustness_r2_total~adata~LIANA_all~all.csv")

# Find best explained view by model and factor across all samples
bev <- r2_total %>% dplyr::group_by(n_factors, Factor) %>%
  dplyr::slice_max(order_by = R2, n = 1) %>%
  dplyr::select(n_factors, Factor, View)
# Drop samples who have is.na(R2) for the best explained view(s)
drop_samples <- r2_sample %>%
  dplyr::semi_join(bev, by = c("n_factors", "Factor", "View")) %>%
  dplyr::filter(is.na(R2)) %>%
  dplyr::select(n_factors, Factor, Sample, View) %>%
  dplyr::distinct()

scores_filtered <- scores %>%
  dplyr::anti_join(drop_samples, by = c("n_factors", "Factor", "Sample"))

#####
# Factor robustness
#####
pivoted_scores_all <- scores_filtered %>%
  tidyr::pivot_wider(names_from = c(n_factors,
                                    Factor),
                     values_from = Score,
                     values_fill = NULL)

pivoted_loadings_all <- loadings %>%
  tidyr::pivot_wider(names_from = c(n_factors, 
                                    Factor),
                     values_from = Loading,
                     values_fill = NULL)

mat_scores_all <- pivoted_scores_all %>%
  dplyr::select(-c(Sample)) %>%
  as.matrix()

corr_scores_all <- Hmisc::rcorr(mat_scores_all,
                                type="pearson")$r

mat_loadings_all <- pivoted_loadings_all %>%
  dplyr::select(-c(view, variable)) %>%
  as.matrix()

corr_loadings_all <- Hmisc::rcorr(mat_loadings_all,
                                  type="pearson")$r


all_correlations <- corr_scores_all*corr_loadings_all

df <- all_correlations %>%
  data.frame %>%
  rownames_to_column("factor1") %>%
  dplyr::mutate(factor1 = stringr::str_c("X", factor1)) %>%
  tidyr::separate(col = factor1, into = c("model1", "factor1")) %>%
  pivot_longer(cols = where(is.numeric),
               names_to = c("model2", "factor2"),
               names_sep = "_",
               values_to = "correlation")

df2 <- df %>%
  dplyr::distinct() %>%
  dplyr::filter(!(model1 == model2 & factor1 == factor2)) %>%
  dplyr::group_by(factor1, model1, model2) %>%
  dplyr::summarise(max_corr = max(correlation)) %>%
  ungroup() %>%
  dplyr::mutate(model1 = as.integer(str_extract(model1, "[:digit:]+")),
                model2 = as.integer(str_extract(model2, "[:digit:]+")))

# Take the mean of the max correlation amongst all models which are bigger:
df3 <- df2 %>%
  dplyr::group_by(factor1, model1) %>%
  dplyr::filter(model2 > model1) %>%
  dplyr::summarise(robustness = mean(max_corr))

# Take the mean across the max correlation of all other models:  
factor_robustness <- df2 %>%
  dplyr::group_by(factor1, model1) %>%
  dplyr::summarise(robustness = mean(max_corr)) %>%
  dplyr::rename(factor = factor1, model = model1)

factor_mean_r2 <- r2_sample %>%
  dplyr::rename(factor = Factor, model = n_factors) %>%
  dplyr::filter(factor != "All") %>%
  dplyr::group_by(factor, model) %>%
  dplyr::summarise(mean_R2 = mean(R2, na.rm = T))

mat <- all_correlations
HM <- Heatmap(mat,
              row_km = 15,
              row_km_repeats = 500,
              column_names_gp = grid::gpar(fontsize = 8),
              row_names_gp = grid::gpar(fontsize = 8))#, column_km = 20, column_km_repeats = 50)
HM <- draw(HM)

HM2 <- Heatmap(mat,
               row_order = unlist(row_order(HM)),
               column_order = unlist(row_order(HM)),
               column_names_gp = grid::gpar(fontsize = 4),
               row_names_gp = grid::gpar(fontsize = 4))#, column_km = 20, column_km_repeats = 50)
HM2 <- draw(HM2)
clustering <- row_order(HM)
factor_clusters <- lapply(1:length(clustering),
                          function(i) row.names(mat[clustering[[i]],])) %>%
  setNames(1:length(clustering)) %>%
  map2_df(seq_along(.), ~ tibble(factor = .x, cluster = .y)) %>%
  tidyr::separate(col = factor, into = c("model", "factor")) %>%
  dplyr::mutate(model = as.integer(model))

df <- purrr::reduce(list(factor_clusters, factor_mean_r2, factor_robustness),
                    dplyr::inner_join,
                    by = c("model", "factor")) 

p1 <- df %>%
  ggplot(aes(x = mean_R2, y = robustness, colour = factor(model))) +
  geom_point() + theme_minimal() + scale_color_discrete("Model")

p2 <- df %>% ggplot(aes(x = mean_R2, y = robustness, colour = factor(cluster))) +
  geom_point() + theme_minimal() + scale_color_discrete("Cluster")

p1 + p2
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
  #dplyr::filter(n_factors >=20 & n_factors <=40) %>%
  tidyr::pivot_wider(names_from = c(n_factors, Factor),
                     values_from = Loading,
                     values_fill = NULL)

mat_loadings <- pivoted_loadings %>%
  dplyr::select(-c(view, variable)) %>%
  as.matrix()

corr_loadings <- Hmisc::rcorr(mat_loadings, type="pearson")$r


####
# Combination of scores and loadings
####
mat <- corr_scores*corr_loadings

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

Heatmap(mat,
        rect_gp = gpar(type = "none"), 
        cell_fun = add_stars_make_diag,
        row_names_gp = gpar(fontsize = 8),
        #column_names_gp = gpar(fontsize = 8),
        cluster_rows = F,
        cluster_row_slices = F,
        #cluster_columns = F,
        cluster_column_slices = F,
        row_split = rep(c(3, 4, 5, 10, 15, 20, 25, 30), times = c(3, 4, 5, 10, 15, 20, 25, 30)),
        row_title = NULL,
        column_split = rep(c(3, 4, 5, 10, 15, 20, 25, 30), times = c(3, 4, 5, 10, 15, 20, 25, 30)),
        row_names_side = "right",
        cluster_columns = F,
        col = f1,
        heatmap_legend_param = list(title = "r"))

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
        #column_names_gp = gpar(fontsize = 8),
        cluster_rows = F,
        cluster_row_slices = F,
        #cluster_columns = F,
        cluster_column_slices = F,
        #row_split = rep(c(3, 4, 5, 10, 15, 20, 25, 30), times = c(3, 4, 5, 10, 15, 20, 25, 30)),
        row_title = NULL,
        column_split = rep(c(3, 4, 5, 10, 15, 20, 25, 30), times = c(3, 4, 5, 10, 15, 20, 25, 30)),
        row_names_side = "right",
        cluster_columns = F,
        col = f1,
        heatmap_legend_param = list(title = "r"))

###
# Median all-factors R2 over all samples, per view
###
r2_all <- r2_sample %>% dplyr::filter(Factor == "All")
to_plot <- r2_all %>%
  dplyr::filter(!stringr::str_detect(View, "&")) %>%
  dplyr::group_by(View, n_factors) %>%
  dplyr::summarise(median_R2 = median(R2, na.rm=TRUE))

ggplot(to_plot, aes(x = n_factors, y = median_R2, color = factor(View))) +
  geom_line() +
  labs(x = "Number of factors in model",
       y = "Median all-factors R2 across samples",
       color = "View") +  theme_minimal()

# Median keeps increasing, but is this just because more and more factors are explaining individual outliers? I need to penalise this somehow.
# For example, max(R2)/median(R2) per factor should not be too high. Or median(R2)/max(R2) should not be too low.
to_plot <- r2_sample %>% dplyr::filter(!Factor == "All") %>%
  dplyr::filter(!stringr::str_detect(View, "&")) %>%
  dplyr::mutate(R2 = pmax(0, R2)) %>% # set minimum R2 to 0
  dplyr::group_by(View, Factor, n_factors) %>%
  dplyr::summarise(max_R2 = max(R2, na.rm = TRUE),
                   median_R2 = median(R2, na.rm=TRUE)) %>%
  dplyr::mutate(explanatory_score = median_R2/max_R2) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(View, n_factors) %>%
  dplyr::summarise(average_explanatory_score = sum(explanatory_score)/n_factors)

ggplot(to_plot, aes(x = n_factors, y = average_explanatory_score, color = factor(View))) +
  geom_line() +
  labs(x = "Number of factors in model",
       y = "Outlier penalty",
       color = "View") +  theme_minimal()

# Probably I want to find out the number of factors for which the explanatory score is bad. Let's define a bad explanatory score
# as where the median sample-level R2 for a factor view is less than 10% of the maximum sample-level R2.
# For example, imagine Factor 4 has an R2 of 50 for one sample in the view CNT, but over all samples, the median R2 is less than 5.
# It's not a perfect metric, because what happens if a Factor explains four samples really well and the rest of them not at all?
# Maybe instead of the median, I need something like the 90th quantile.
# So, Factor with limited explanatory power in terms of number of samples well explained
# is one where the 90th quantile sample-level R2 is less than 10% of the maximum sample-level R2.

to_plot <- r2_sample %>% dplyr::filter(!Factor == "All") %>%
  dplyr::filter(!stringr::str_detect(View, "&")) %>%
  dplyr::mutate(R2 = pmax(0, R2)) %>% # set minimum R2 to 0
  dplyr::group_by(View, Factor, n_factors) %>%
  dplyr::summarise(max_R2 = max(R2, na.rm = TRUE),
                   ninetieth_quantile = quantile(R2, 0.9, na.rm=TRUE)) %>%
  dplyr::mutate(explanatory_score = ninetieth_quantile/max_R2) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(View, n_factors) %>%
  dplyr::summarise(average_explanatory_score = sum(explanatory_score)/n_factors)


# Even simpler would be just to take the 95th quantile sample-level R2 as my summary statistic.
# I want to plot the average 95th quantile sample-level R2 across all factors.
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

####
# R2
####

r2_all <- r2_sample %>% dplyr::filter(Factor == "All")

n_sample_views <- r2_all %>%
  tidyr::drop_na() %>%
  dplyr::select(Sample, View) %>%
  dplyr::distinct() %>% nrow()

thresholds <- c(0, 10, 25, 50)

to_plot <- lapply(thresholds, function(threshold){
  df <- r2_all %>% 
    dplyr::filter(R2 >= threshold) %>%
    dplyr::group_by(n_factors) %>%
    dplyr::summarize(cumul_prop = n_distinct(paste(View,
                                                   Sample,
                                                   sep = "_")) / n_sample_views)
  df$threshold <- threshold
  return(df)
})

to_plot <- dplyr::bind_rows(to_plot)

ggplot(to_plot, aes(x = n_factors, y = cumul_prop, color = factor(threshold))) +
  geom_line() +
  labs(x = "Number of factors in model",
       y = "% views reaching R2 threshold",
       colour = "R2 threshold") +
  ylim(0,1) +  theme_minimal()

df <- r2_all %>%
  dplyr::group_by(View, Factor, n_factors) %>%
  dplyr::summarise(median_R2 = median(R2, na.rm = T),
                   mean_R2 = mean(R2, na.rm = T))

library(ggrepel)

p1 <- r2_total %>% dplyr::filter(Factor == 'All') %>%
  dplyr::filter(stringr::str_detect(View, "&", negate = TRUE)) %>%
  dplyr::mutate(label = ifelse(n_factors == max(n_factors),
                               as.character(View),
                               NA_character_)) %>%
  ggplot(aes(x = n_factors,
             y = R2,
             color = View)) +
  geom_line() + theme_minimal() +
  geom_label_repel(aes(label = label),
                   size = 2.5,
                   direction = "y",
                   box.padding = 0.1,
                   nudge_x = 0.5,
                   hjust = 0,
                   na.rm = TRUE,
                   segment.alpha = 0.5,
                   segment.linetype = 3,
                   max.overlaps = 10,
                   label.padding = 0.1,
                   label.size = NA) + 
  scale_color_discrete(guide = "none") +
  xlim(0,40) +
  ylim(0,100) +
  labs(x = "Number of factors in model",
       y = expression("Total" ~ R^2),
       subtitle = "Cell type transcriptome views")

liana <- r2_total %>%
  dplyr::filter(Factor == 'All') %>%
  dplyr::filter(stringr::str_detect(View, "&", negate=FALSE)) %>%
  tidyr::separate_wider_delim(cols = View, delim = "&", names = c("Sender", "Receiver"))

receiver <- liana %>%
  dplyr::group_by(Receiver, Factor, n_factors) %>%
  dplyr::summarise(mean_r2 = mean(R2, na.rm = TRUE))

sender <- liana %>%
  dplyr::group_by(Sender, Factor, n_factors) %>%
  dplyr::summarise(mean_r2 = mean(R2, na.rm = TRUE))

p2 <- receiver %>%
  dplyr::mutate(label = ifelse(n_factors == max(n_factors),
                               as.character(Receiver),
                               NA_character_)) %>%
  ggplot(aes(x = n_factors,
             y = mean_r2,
             color = Receiver)) +
  geom_line() + theme_minimal() +
  geom_label_repel(aes(label = label),
                   size = 2.5,
                   direction = "y",
                   box.padding = 0.1,
                   nudge_x = 0.5,
                   hjust = 0,
                   na.rm = TRUE,
                   segment.alpha = 0.5,
                   segment.linetype = 3,
                   max.overlaps = 10,
                   label.padding = 0.1,
                   label.size = NA) + 
  scale_color_discrete(guide = "none") +
  xlim(0,40) +
  ylim(0,100) +
  labs(x = "Number of factors in model",
       y = expression("Mean total" ~ R^2),
       subtitle = "L-R views: Summarised by receiver cell type")

p3 <- sender %>%
  dplyr::mutate(label = ifelse(n_factors == max(n_factors),
                               as.character(Sender),
                               NA_character_)) %>%
  ggplot(aes(x = n_factors,
             y = mean_r2,
             color = Sender)) +
  geom_line() + theme_minimal() +
  geom_label_repel(aes(label = label),
                   size = 2.5,
                   direction = "y",
                   box.padding = 0.1,
                   nudge_x = 0.5,
                   hjust = 0,
                   na.rm = TRUE,
                   segment.alpha = 0.5,
                   segment.linetype = 3,
                   max.overlaps = 10,
                   label.padding = 0.1,
                   label.size = NA) + 
  scale_color_discrete(guide = "none") +
  xlim(0,40) +
  ylim(0,100) +
  labs(x = "Number of factors in model",
       y = expression("Mean total" ~ R^2),
       subtitle = "L-R views: Summarised by sender cell type")
design <- "
1122
1122
3344
3344
"
patchwork::wrap_plots(p1, plot_spacer(), p3, p2, design = design)

df %>% dplyr::filter(Factor == 'All') %>%
  dplyr::mutate(label = ifelse(n_factors == max(n_factors),
                               as.character(View),
                               NA_character_)) %>%
  ggplot(aes(x = n_factors,
             y = median_R2,
             color = View)) +
  geom_line() + theme_minimal() +
  geom_label_repel(aes(label = label),
                   size = 3.5,
                   direction = "y",
                   box.padding = 0.1,
                   nudge_x = 0.5,
                   hjust = 0,
                   na.rm = TRUE,
                   segment.alpha = 0.5,
                   segment.linetype = 3,
                   label.size = NA) + 
  scale_color_discrete(guide = "none") +
  xlim(0,25) +
  labs(x = "Number of factors in model",
       y = expression("Median"~ R^2 ~ "across samples"))

df %>% dplyr::filter(Factor == 'All') %>%
  dplyr::mutate(label = ifelse(n_factors == max(n_factors),
                               as.character(View),
                               NA_character_)) %>%
  ggplot(aes(x = n_factors,
             y = mean_R2,
             color = View)) +
  geom_line() + theme_minimal() +
  geom_label_repel(aes(label = label),
                   size = 3.5,
                   direction = "y",
                   box.padding = 0.1,
                   nudge_x = 0.5,
                   hjust = 0,
                   na.rm = TRUE,
                   segment.alpha = 0.5,
                   segment.linetype = 3,
                   label.size = NA) + 
  scale_color_discrete(guide = "none") +
  xlim(0,25) +
  labs(x = "Number of factors in model",
       y = expression("Mean"~ R^2 ~ "across samples"))


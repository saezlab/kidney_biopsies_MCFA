library(ggplot2)
library(scales)
library(dplyr)
library(patchwork)
library(tidyr)

loadings <- read.csv("results/MCFA/loadings~adata~LIANA_all~all~5.csv",
                     row.names = 1)

loadings_f2 <- loadings %>% dplyr::filter(Factor == "Factor2")
to_plot_classic <- loadings_f2 %>%
  dplyr::filter(!stringr::str_detect(view, "&"))
to_plot_liana <- loadings_f2 %>%
  dplyr::filter(stringr::str_detect(view, "&"))

n = 15
views_of_interest <- c("CNT_PC", "DCT", "dPT_dTAL_dDCT", "TAL", "PT", "Fib_VSMC", "aPT_aTAL_aDCT")
top_variables <- lapply(views_of_interest, function(vof){
  to_plot_classic %>%
    dplyr::filter(view == vof) %>%
    dplyr::slice_max(abs(Loading), n = n) %>%
    dplyr::select(variable) %>% unlist()
}) %>% unlist()


top_loadings <- to_plot_classic %>%
  dplyr::filter(variable %in% top_variables) %>%
  dplyr::filter(view %in% views_of_interest) %>%
  dplyr::arrange(desc(Loading)) %>%
  dplyr::mutate(variable = forcats::fct_inorder(variable)) %>%
  dplyr::mutate(variable = forcats::fct_rev(variable))


library(ComplexHeatmap)
library(circlize)

# Filter data
top_loadings_filtered <- top_loadings %>%
  filter(view %in% views_of_interest) %>%
  pivot_wider(names_from = view, values_from = Loading)  # Reshape data

# Convert to matrix
loading_matrix <- top_loadings_filtered %>%
  tibble::column_to_rownames('variable') %>%
  dplyr::select(-Factor) %>% as.matrix

# Perform hierarchical clustering on rows and define row splits
row_dend <- hclust(dist(loading_matrix))
row_clusters <- cutree(row_dend, k = 9)  # Adjust k to change number of row groups

# Create heatmap with row splits
ht <- Heatmap(
  loading_matrix,
  name = "Loading",
  col = colorRamp2(c(min(loading_matrix, na.rm = TRUE), 0, max(loading_matrix, na.rm = TRUE)), 
                   c("blue", "white", "red")),
  cluster_rows = TRUE,  
  cluster_columns = TRUE,  
  show_column_names = TRUE,
  show_row_names = TRUE,
  #row_split = row_clusters,  # Split rows into clusters
  column_names_rot = 45,
  heatmap_legend_param = list(title = "Loading"),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10)
)

draw(ht, padding = unit(c(5, 10, 5, 5), "mm"))
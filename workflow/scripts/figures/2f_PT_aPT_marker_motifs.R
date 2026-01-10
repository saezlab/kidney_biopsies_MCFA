library(magrittr)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(tidyr)
library(tibble)

disruptedMotifs <- read.csv("data/ATAC/top_disrupted_tfs.csv")
seEnrich <- readRDS("data/ATAC/enrich_motifs.rds")

mat <- assays(seEnrich)[["mlog10Padj"]]
mat_unfiltered <- assays(seEnrich)[["mlog10Padj"]]
write.csv(mat_unfiltered, "results/ATAC/motif_enrichment_mlog10Padj.csv")

padj <- mat_unfiltered %>%
  dplyr::select(aPT, PT) %>%
  tibble::rownames_to_column("motif") %>%
  tidyr::pivot_longer(cols = where(is.numeric),
                      names_to = "cluster_label",
                      values_to = "mlog10Padj") %>%
  dplyr::mutate(padj = 10^(-mlog10Padj))

enrichment <- assays(seEnrich)[["Enrichment"]] %>%
  dplyr::select(aPT, PT) %>%
  tibble::rownames_to_column("motif") %>%
  tidyr::pivot_longer(cols = where(is.numeric),
                      names_to = "cluster_label",
                      values_to = "ES") %>%
  dplyr::left_join(padj) %>%
  dplyr::mutate(annotation = ifelse(padj < 0.05, "∗", "")) %>%
  dplyr::mutate(geneSymbol = stringr::str_extract(motif, "^[:alnum:]+(?=_)")) %>%
  dplyr::left_join(disruptedMotifs)

motif_to_plot <-  enrichment %>% dplyr::group_by(cluster_label) %>%
  dplyr::slice_max(order_by = ES, n = 20) %>%
  dplyr::pull(motif) %>% unique()

motif_to_plot2 <- enrichment %>%
  dplyr::group_by(motif) %>%
  dplyr::slice_max(order_by = ES, n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::slice_max(order_by = ES, n = 30) %>%
  dplyr::pull(motif) %>% unique()

df <- enrichment %>%
  dplyr::filter(motif %in% motif_to_plot2)
# Step 1: Convert data to a matrix for the heatmap (rows = cluster_label, columns = motif)
mat <-  df %>%
  dplyr::select(cluster_label, motif, ES) %>%
  tidyr::pivot_wider(names_from = motif, values_from = ES) %>%
  tibble::column_to_rownames("cluster_label") %>%
  as.matrix() %>% t()

# Step 2: Prepare a corresponding annotation matrix for text (e.g., "*" marks)
text_mat <- df %>%
  dplyr::select(cluster_label, motif, annotation) %>%
  tidyr::pivot_wider(names_from = motif, values_from = annotation) %>%
  tibble::column_to_rownames("cluster_label") %>%
  as.matrix() %>% t()

# Define color scale for ES
col_fun <- colorRamp2(
  seq(0, max(df$ES, na.rm = TRUE), length.out = 100),
  viridis(100, option = "E")
)

# Step 5: Create the heatmap
ht <- Heatmap(
  mat,
  name = "ES",
  col = col_fun,
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(text_mat[i, j], x, y, gp = gpar(col = "white", fontsize = 20, fontfamily = "Courier"))
  },
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_title = NULL,
  row_title = NULL,
  show_row_dend = FALSE,
  column_order = c("PT", "aPT"),
  column_names_rot = 45,
  row_names_side = "left"
  #row_names_rot = 90
)

# Step 6: Draw it
draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     padding = unit(c(2, 2, 10, 10), "pt"))
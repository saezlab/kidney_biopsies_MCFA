library(magrittr)
library(ggplot2)
library(dplyr)

annotations <- read.csv('data/annotations/label_transfer_RNA_ATAC_V3_clean.csv', row.names = 1) %>%
  tibble::rownames_to_column("cell")
coordinates <- read.csv('data/ATAC/coordinates.csv', row.names = 1)

to_plot <- annotations %>%
  dplyr::left_join(coordinates, by = "cell") %>%
  dplyr::mutate(cluster_label = ifelse(cluster_label=="Pod", "POD", cluster_label))

cluster_vs_level1 <- to_plot %>%
  dplyr::select(cluster, level1) %>%
  group_by(cluster, level1) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(percentage = count/sum(count)*100)
cluster_vs_level1 %>%
  ggplot(aes(y = forcats::fct_rev(cluster), x = level1, fill = percentage)) +
  geom_tile() + scale_fill_viridis_c() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "snATAC-seq cluster", x = "snATAC-seq manual annotation") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45))

level1_vs_label_transfer <- to_plot %>%
  dplyr::select(level1, clean_label) %>%
  group_by(level1, clean_label) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(level1) %>%
  mutate(percentage = count / sum(count) * 100)

level1_vs_label_transfer %>%
  ggplot(aes(x = level1, y = clean_label, fill = percentage)) +
  geom_tile() + scale_fill_viridis_c() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "snATAC-seq manual annotation", y = "Label transfer annotation", fill = "% matching\nprelim.\nannotation") 

cluster_vs_label_transfer <- to_plot %>%
  dplyr::select(cluster, clean_label) %>%
  group_by(cluster, clean_label) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::complete(cluster, clean_label, fill = list(count = 0)) %>%
  group_by(cluster) %>%
  mutate(proportion = count / sum(count))

cluster_vs_label_transfer %>%
  ggplot(aes(x = forcats::fct_rev(cluster), y = clean_label, fill = proportion)) +
  geom_tile() + scale_fill_viridis_c() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "snATAC-seq cluster", y = "Label transfer from snRNA-seq", fill = "Proportion") + coord_flip()

final_vs_label_transfer <- to_plot %>%
  dplyr::select(cluster_label, clean_label) %>%
  group_by(cluster_label, clean_label) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::complete(cluster_label, clean_label, fill = list(count = 0)) %>%
  group_by(cluster_label) %>%
  mutate(proportion = count / sum(count))

final_vs_label_transfer %>%
  ggplot(aes(x = cluster_label, y = forcats::fct_rev(clean_label), fill = proportion)) +
  geom_tile() + scale_fill_viridis_c() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "top") +
  labs(x = "snATAC-seq annotation", y = "Label transfer from snRNA-seq", fill = "Proportion") 


label.df_2 <- to_plot %>% 
  group_by(cluster_label) %>% 
  summarize(X = ifelse(median(X) < 0, median(X)-1.2, median(X) + 1.2), Y = ifelse(median(Y) < 0, median(Y) -1.2, median(Y) + 1.2))

distinct_colors <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999", # gray
  "#A65628", # brown
  "#FF69B4", # hot pink
  "#87CEEB", # light blue
  "#228B22", # forest green
  "#FFD700", # gold
  "#6A5ACD", # slate blue
  "#DC143C", # crimson
  "#00CED1", # dark turquoise
  "#8B0000", # dark red
  "#B8860B",  # dark goldenrod
  "#9400D3", # dark violet
  "#40E0D0"  # turquoise
)

to_plot %>%
  ggplot(aes(x = X, y = Y, color = cluster_label)) + geom_point(size = 0.2) + theme_classic() +
  geom_label(data = label.df_2, aes(label = cluster_label)) + scale_color_manual(values = distinct_colors) +
  theme(legend.position = "none", axis.ticks = element_blank(), axis.text = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2")



library(magrittr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(grid)
library(RColorBrewer)

if (exists("snakemake")){
  data <- read.csv(snakemake@input$data)
  plt <- snakemake@output$plot
  gg <- snakemake@output$gg
} else{
  data <- read.csv('results/data_overview/modality_with_duplicates~heatmap.csv')
}

### Version with Xenium information
df <- data
# Build modality matrix
mat <- t(df[, c("rna","atac", "xenium")])
colnames(mat) <- df$patient_ID
rownames(mat) <- c("RNA","ATAC", "Xenium")

# Define ordering key
order_key <- with(df,
                  ifelse(rna & atac & xenium, 1,
                         ifelse(rna & atac & !xenium, 2,
                                ifelse(!rna & atac & xenium, 3, 
                                       ifelse(rna & !atac & xenium, 4,
                                              ifelse(rna & !atac & !xenium, 5,
                                                     ifelse(!rna & atac & !xenium, 6, 7)
                                                     )
                                              )
                                       )
                                )
                         )
                  )

# Sort columns by (disease, then order_key)
disease_levels <- c("DN", "FSGS", "IgAN", "MN", "TIN",
                   "AKI (control)", "control", "control_KPMP",
                   "control_preKTX")

df$disease <- factor(df$disease, levels = disease_levels)
ord <- order(df$disease, order_key)

mat_ord <- mat[, ord, drop = FALSE]
disease_ord <- df$disease[ord]

# Define colors
mod_cols <- c(
  `0` = "#f0f0f0",
  `1` = "#1b9e77",
  `2` = "#01665e"
)


set2_9 <- c(
  "#66C2A5", # teal
  "#FC8D62", # salmon
  "#8DA0CB", # periwinkle
  "#E78AC3", # pink
  "#A6D854", # light green
  "#FFD92F", # yellow
  "#E5C494", # beige
  "#B3B3B3", # gray
  "#B39EB5"  # lavender (added)
)
disease_colors <- setNames(set2_9, disease_levels)

# Compute counts per disease
disease_counts <- table(disease_ord)

# Use formatted titles like "A (n=3)" or just "n=3"
custom_titles <- disease_counts#paste0("n=", disease_counts)

# Heatmap
ha <- HeatmapAnnotation(
  Condition = disease_ord,
  col = list(Condition = disease_colors)
)

ht <- Heatmap(
  mat_ord[c('RNA', 'ATAC', 'Xenium'),],
  name = "Samples",
  col = mod_cols,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  top_annotation = ha,
  #bottom_annotation = hb,
  column_split = disease_ord,   # Conditions by disease
  rect_gp = gpar(col = "white"),
  row_title=NULL,
  column_title=custom_titles, #"Donor",
  show_column_names = FALSE,
  column_title_side = "bottom"
)
draw(ht, padding = unit(c(15, 2, 15, 2), "mm"))

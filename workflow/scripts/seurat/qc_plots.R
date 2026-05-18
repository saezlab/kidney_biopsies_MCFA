if (exists("snakemake")) {
  rds_path <- snakemake@input[["rds"]]
  plot_path <- snakemake@output[["plot"]]
} else {
  rds_path <- "results/seurat/sce~K1.rds"
  plot_path <- "plots/seurat/qc~K1.pdf"
}

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(ggplot2)
  library(patchwork)
})

sce <- readRDS(rds_path)
counts <- SummarizedExperiment::assay(sce, "counts")
gene_names <- rownames(sce)

if (is.null(gene_names)) {
  stop("SCE rownames are required to compute QC feature classes.", call. = FALSE)
}

is_mito <- grepl("^MT-", gene_names)
is_ribo <- grepl("^RPS|^RPL", gene_names)

if (!any(is_mito)) {
  warning("No mitochondrial genes matched '^MT-'; pct_counts_mt will be zero.")
}
if (!any(is_ribo)) {
  warning("No ribosomal genes matched '^RPS|^RPL'; pct_counts_ribo will be zero.")
}

total_counts <- Matrix::colSums(counts)
n_genes_by_counts <- Matrix::colSums(counts > 0)

mt_counts <- if (any(is_mito)) Matrix::colSums(counts[is_mito, , drop = FALSE]) else rep(0, ncol(counts))
ribo_counts <- if (any(is_ribo)) Matrix::colSums(counts[is_ribo, , drop = FALSE]) else rep(0, ncol(counts))

pct_counts_mt <- ifelse(total_counts > 0, mt_counts / total_counts * 100, 0)
pct_counts_ribo <- ifelse(total_counts > 0, ribo_counts / total_counts * 100, 0)

qc_df <- data.frame(
  total_counts = as.numeric(total_counts),
  n_genes_by_counts = as.numeric(n_genes_by_counts),
  pct_counts_mt = as.numeric(pct_counts_mt),
  pct_counts_ribo = as.numeric(pct_counts_ribo)
)

p1 <- ggplot(qc_df, aes(x = total_counts)) +
  geom_histogram(bins = 100, color = "white", fill = "#4C956C") +
  theme_classic() +
  labs(x = "Total counts", y = "Cells", title = "Library size distribution")

p2 <- ggplot(qc_df, aes(x = "", y = pct_counts_mt)) +
  geom_violin(fill = "#F4A259", color = "#A44A3F", trim = FALSE) +
  theme_classic() +
  labs(x = NULL, y = "Percent mitochondrial counts", title = "Mitochondrial content") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p3 <- ggplot(qc_df, aes(x = total_counts, y = n_genes_by_counts, color = pct_counts_mt)) +
  geom_point(alpha = 0.7, size = 0.6) +
  scale_color_viridis_c(option = "C") +
  theme_classic() +
  labs(
    x = "Total counts",
    y = "Genes by counts",
    color = "% MT",
    title = "Counts vs detected genes"
  )

p4 <- ggplot(qc_df, aes(x = "", y = pct_counts_ribo)) +
  geom_violin(fill = "#7CAFC2", color = "#2A4D69", trim = FALSE) +
  theme_classic() +
  labs(x = NULL, y = "Percent ribosomal counts", title = "Ribosomal content") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

combined_plot <- (p1 | p2) / (p3 | p4)

dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
ggsave(plot_path, combined_plot, width = 12, height = 10)

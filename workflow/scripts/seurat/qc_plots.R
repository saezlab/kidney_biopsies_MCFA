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

col_data <- as.data.frame(SummarizedExperiment::colData(sce))
if (!"scDblFinder.score" %in% colnames(col_data) || !"scDblFinder.class" %in% colnames(col_data)) {
  stop(
    "scDblFinder metadata is required for QC plotting. Expected scDblFinder.score and scDblFinder.class in colData.",
    call. = FALSE
  )
}

doublet_score <- as.numeric(col_data[["scDblFinder.score"]])
doublet_class <- as.character(col_data[["scDblFinder.class"]])
doublet_threshold <- min(doublet_score[doublet_class == "doublet"])

mt_over_5 <- pct_counts_mt > 5
n_mt_over_5 <- sum(mt_over_5)
pct_mt_over_5 <- 100 * n_mt_over_5 / length(mt_over_5)

doublets_called <- doublet_class == "doublet"
n_doublets <- sum(doublets_called)
pct_doublets <- 100 * n_doublets / length(doublets_called)

qc_df <- data.frame(
  total_counts = as.numeric(total_counts),
  n_genes_by_counts = as.numeric(n_genes_by_counts),
  pct_counts_mt = as.numeric(pct_counts_mt),
  pct_counts_ribo = as.numeric(pct_counts_ribo),
  doublet_score = doublet_score,
  doublet_class = doublet_class
)

p1 <- ggplot(qc_df, aes(x = total_counts)) +
  geom_histogram(bins = 100, color = "white", fill = "#4C956C") +
  theme_classic() +
  labs(x = "Total counts", y = "Cells", title = "Library size distribution")

p2 <- ggplot(qc_df, aes(x = pct_counts_mt)) +
  geom_histogram(bins = 100, color = "white", fill = "#F4A259") +
  geom_vline(xintercept = 5, linetype = "dashed", color = "#A44A3F", linewidth = 0.8) +
  theme_classic() +
  labs(x = "Percent mitochondrial counts", y = "Cells", title = "Mitochondrial content") +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf(">5%% MT: %d cells (%.1f%%)", n_mt_over_5, pct_mt_over_5),
    hjust = 1.05,
    vjust = 1.5,
    size = 3.5
  )

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

p4 <- ggplot(qc_df, aes(x = pct_counts_ribo)) +
  geom_histogram(bins = 100, color = "white", fill = "#7CAFC2") +
  theme_classic() +
  labs(x = "Percent ribosomal counts", y = "Cells", title = "Ribosomal content")

p5 <- ggplot(qc_df, aes(x = total_counts, y = n_genes_by_counts, color = doublet_score)) +
  geom_point(alpha = 0.7, size = 0.6) +
  scale_color_viridis_c(option = "D") +
  theme_classic() +
  labs(
    x = "Total counts",
    y = "Genes by counts",
    color = "Doublet score",
    title = "Counts vs detected genes by doublet score"
  )

p6 <- ggplot(qc_df, aes(x = doublet_score)) +
  geom_histogram(bins = 100, color = "white", fill = "#B56576") +
  geom_vline(xintercept = doublet_threshold, linetype = "dashed", color = "#6D597A", linewidth = 0.8) +
  theme_classic() +
  labs(x = "Doublet score", y = "Cells", title = "scDblFinder score distribution") +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf("Doublets: %d cells (%.1f%%)", n_doublets, pct_doublets),
    hjust = 1.05,
    vjust = 1.5,
    size = 3.5
  ) +
  annotate(
    "text",
    x = doublet_threshold,
    y = Inf,
    label = sprintf("threshold = %.3f", doublet_threshold),
    hjust = -0.05,
    vjust = 1.5,
    angle = 90,
    size = 3
  )

combined_plot <- (p1 | p2) / (p3 | p4) / (p5 | p6)

dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
ggsave(plot_path, combined_plot, width = 12, height = 15)

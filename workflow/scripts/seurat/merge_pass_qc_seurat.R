if (exists("snakemake")) {
  rds_paths <- snakemake@input[["rds"]]
  out_path <- snakemake@output[["rds"]]
} else {
  rds_paths <- c("results/seurat/sce~K1.rds")
  out_path <- "results/seurat/merged_seurat.rds"
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
})

extract_sample <- function(path) {
  sub("^sce~(.*)\\.rds$", "\\1", basename(path))
}

make_seurat <- function(path) {
  sample <- extract_sample(path)
  sce <- readRDS(path)

  counts <- SummarizedExperiment::assay(sce, "counts")
  meta <- as.data.frame(SummarizedExperiment::colData(sce))

  if (!"scDblFinder.class" %in% colnames(meta)) {
    stop("Missing scDblFinder.class in ", sample, call. = FALSE)
  }

  mt <- grepl("^MT-", rownames(sce))
  total_counts <- Matrix::colSums(counts)
  mt_counts <- Matrix::colSums(counts[mt, , drop = FALSE])
  pct_mt <- mt_counts / total_counts * 100

  keep <- meta$scDblFinder.class == "singlet" & pct_mt < 5
  if (!any(keep)) {
    stop("No cells remained after filtering for ", sample, call. = FALSE)
  }

  seu <- Seurat::CreateSeuratObject(
    counts = counts[, keep, drop = FALSE],
    meta.data = meta[keep, , drop = FALSE],
    project = sample
  )

  seu$sample <- sample
  seu$pct_counts_mt <- pct_mt[keep]

  Seurat::RenameCells(seu, add.cell.id = sample)
}

objs <- lapply(rds_paths, make_seurat)
names(objs) <- vapply(rds_paths, extract_sample, character(1))

merged <- if (length(objs) == 1) {
  objs[[1]] # Dummy merge for checking script outside of Snakemake
} else {
  merge( # S3 method, not exported from Seurat, but for which Seurat provides a method for its own object class
    x = objs[[1]],
    y = objs[-1],
    add.cell.ids = NULL,
    merge.data = FALSE
  )
}

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(merged, out_path)

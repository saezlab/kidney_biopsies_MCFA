if (exists("snakemake")) {
  rds_path <- snakemake@input[["rds"]]
  mtx_path <- snakemake@output[["mtx"]]
  barcodes_path <- snakemake@output[["barcodes"]]
  features_path <- snakemake@output[["features"]]
  variable_features_path <- snakemake@output[["variable_features"]]
} else {
  rds_path <- "results/seurat/sct_harmony.rds"
  mtx_path <- "results/seurat/raw_counts.mtx.gz"
  barcodes_path <- "results/seurat/raw_counts_barcodes.tsv.gz"
  features_path <- "results/seurat/raw_counts_features.tsv.gz"
  variable_features_path <- "results/seurat/sct_variable_features.csv"
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(readr)
  library(tibble)
})

write_tsv_gz <- function(df, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x = df, file = con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

print("Reading Seurat object...")
obj <- readRDS(rds_path)

if (!"SCT" %in% names(obj@assays)) {
  stop("Seurat object does not contain an SCT assay.", call. = FALSE)
}

top_sct_genes <- VariableFeatures(obj[["SCT"]])

write.csv(top_sct_genes, variable_features_path)

if (!"RNA" %in% names(obj@assays)) {
  stop("Seurat object does not contain an RNA assay.", call. = FALSE)
}

print("Joining RNA layers...")
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

print("Extracting RNA counts matrix...")
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")

if (nrow(counts) == 0 || ncol(counts) == 0) {
  stop("RNA counts matrix is empty.", call. = FALSE)
}
if (!inherits(counts, "dgCMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

print(counts)

feature_ids <- rownames(counts)
cell_ids <- colnames(counts)

print(head(feature_ids))
print(head(cell_ids))

if (is.null(feature_ids) || anyNA(feature_ids) || any(feature_ids == "")) {
  stop("RNA counts matrix must have non-empty feature names.", call. = FALSE)
}

if (is.null(cell_ids) || anyNA(cell_ids) || any(cell_ids == "")) {
  stop("RNA counts matrix must have non-empty cell names.", call. = FALSE)
}

print("Writing RNA counts matrix...")
Matrix::writeMM(obj = counts, file = mtx_path)
write_tsv_gz(data.frame(feature = feature_ids, stringsAsFactors = FALSE), features_path)
write_tsv_gz(data.frame(cellid = cell_ids, stringsAsFactors = FALSE), barcodes_path)

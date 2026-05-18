if (exists("snakemake")) {
  matrix_path <- snakemake@input[["matrix"]]
  features_path <- snakemake@input[["features"]]
  barcodes_path <- snakemake@input[["barcodes"]]
  out_path <- snakemake@output[["rds"]]
} else {
  matrix_path <- "data/RNA/SoupX/K1/matrix.mtx.gz"
  features_path <- "data/RNA/SoupX/K1/features.tsv.gz"
  barcodes_path <- "data/RNA/SoupX/K1/barcodes.tsv.gz"
  out_path <- "results/seurat/sce~K1.rds"
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
})

required_paths <- c(matrix_path, features_path, barcodes_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop(
    "Missing SoupX input files: ",
    paste(missing_paths, collapse = ", "),
    call. = FALSE
  )
}

data_dir <- dirname(matrix_path)
counts <- Read10X(data.dir = data_dir)
seurat_obj <- CreateSeuratObject(
  counts = counts,
  min.cells = 3,
  min.features = 200
)
sce <- as.SingleCellExperiment(seurat_obj)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(sce, out_path)

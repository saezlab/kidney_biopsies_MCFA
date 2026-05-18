if (exists("snakemake")) {
  rds_path <- snakemake@input[["rds"]]
  flag_path <- snakemake@output[["flag"]]
} else {
  rds_path <- "results/seurat/sce~K1.rds"
  flag_path <- "results/seurat/flags/scDblFinder~K1.txt"
}

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scDblFinder)
})

set.seed(1234)
sce <- readRDS(rds_path)
sce <- scDblFinder(sce)
saveRDS(sce, rds_path)

dir.create(dirname(flag_path), recursive = TRUE, showWarnings = FALSE)
writeLines("scDblFinder complete; seed=1234", con = flag_path)

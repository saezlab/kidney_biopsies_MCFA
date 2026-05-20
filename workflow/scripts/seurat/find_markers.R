if (exists("snakemake")) {
  rds_path <- snakemake@input[["rds"]]
  markers_csv <- snakemake@output[["csv"]]
  resolution <- snakemake@params[["resolution"]]
} else {
  rds_path <- "results/seurat/sct_harmony.rds"
  markers_csv <- "results/seurat/markers/findallmarkers_res_0.2.csv"
  resolution <- "0.2"
}

suppressPackageStartupMessages({
  library(Seurat)
  library(readr)
  library(tibble)
})

cluster_col <- paste0("SCT_snn_res.", resolution)
obj <- readRDS(rds_path)

if (!cluster_col %in% colnames(obj@meta.data)) {
  stop(
    "Missing expected cluster metadata column: ",
    cluster_col,
    call. = FALSE
  )
}

res_label <- sub("^SCT_snn_res\\.", "", cluster_col)
obj <- Seurat::PrepSCTFindMarkers(obj)
Seurat::Idents(obj) <- obj[[cluster_col, drop = TRUE]]
resolution_markers <- Seurat::FindAllMarkers(
  object = obj,
  assay = "SCT",
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1,
  verbose = TRUE
)

if (is.null(resolution_markers) || ncol(resolution_markers) == 0) {
  resolution_markers <- tibble::tibble()
}

readr::write_csv(resolution_markers, markers_csv)

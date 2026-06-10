if (exists("snakemake")) {
  rds_path <- snakemake@input[["rds"]]
  annotation_rds <- snakemake@output[["rds"]]
  metadata_csv <- snakemake@output[["metadata_csv"]]
  pca_embeddings_csv <- snakemake@output[["pca_embeddings_csv"]]
  harmony_embeddings_csv <- snakemake@output[["harmony_embeddings_csv"]]
  umap_pre_embeddings_csv <- snakemake@output[["umap_pre_embeddings_csv"]]
  umap_embeddings_csv <- snakemake@output[["umap_embeddings_csv"]]
  pca_sd_csv <- snakemake@output[["pca_sd_csv"]]
  clustree_plot <- snakemake@output[["clustree_plot"]]
  session_info_csv <- snakemake@output[["session_info_csv"]]
} else {
  rds_path <- "results/seurat/merged_seurat.rds"
  annotation_rds <- "results/seurat/sct_harmony.rds"
  metadata_csv <- "results/seurat/metadata.csv"
  pca_embeddings_csv <- "results/seurat/pca_embeddings.csv"
  harmony_embeddings_csv <- "results/seurat/harmony_embeddings.csv"
  umap_pre_embeddings_csv <- "results/seurat/umap_pre_harmony_embeddings.csv"
  umap_embeddings_csv <- "results/seurat/umap_embeddings.csv"
  pca_sd_csv <- "results/seurat/pca_stdev.csv"
  clustree_plot <- "plots/seurat/annotation/clustree_sct_harmony.pdf"
  session_info_csv <- "results/seurat/session_info.csv"
}

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(clustree)
  library(readr)
  library(tibble)
  library(dplyr)
})
#Error in getGlobalsAndPackages(expr, envir = envir, globals = globals) : 
# The total size of the 19 globals exported for future expression (‘FUN()’) is 760.16 MiB.
# This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. 
# This limit is set to protect against transfering too large objects to parallel workers by mistake,
# which may not be intended and could be costly. 
# See help("future.globals.maxSize", package = "future") 
# for how to adjust or remove the default threshold via an R option 
# The three largest globals are ‘FUN’ (740.11 MiB of class ‘function’), 
# ‘umi_bin’ (19.22 MiB of class ‘numeric’) and ‘data_step1’ (771.95 KiB of class ‘list’)
options(future.globals.maxSize = Inf)

log_progress <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), paste0(..., collapse = ""))
}

normalize_cluster_columns <- function(meta, resolution_labels) {
  cluster_cols <- grep("^SCT_snn_res\\.", colnames(meta), value = TRUE)

  for (resolution_label in resolution_labels) {
    expected_col <- paste0("SCT_snn_res.", resolution_label)
    if (expected_col %in% colnames(meta)) {
      next
    }

    target_resolution <- as.numeric(resolution_label)
    numeric_matches <- vapply(
      cluster_cols,
      function(colname) {
        suppressWarnings(as.numeric(sub("^SCT_snn_res\\.", "", colname)))
      },
      numeric(1)
    )

    matched_cols <- cluster_cols[!is.na(numeric_matches) & numeric_matches == target_resolution]
    if (length(matched_cols) == 1) {
      colnames(meta)[colnames(meta) == matched_cols] <- expected_col
      cluster_cols[cluster_cols == matched_cols] <- expected_col
    }
  }

  meta
}

dims_use <- 1:30
resolution_grid <- seq(0.1, 2.0, by = 0.1)
resolution_labels <- formatC(
  resolution_grid,
  format = "f",
  digits = 1
)

log_progress("Reading merged Seurat object: ", rds_path)
obj <- readRDS(rds_path)
log_progress("Loaded object with ", ncol(obj), " cells and ", nrow(obj), " features.")

if ("SCT" %in% names(obj@assays)) {
  log_progress("Using existing SCT assay.")
  DefaultAssay(obj) <- "SCT"
} else {
  log_progress("Running SCTransform.")
  DefaultAssay(obj) <- "RNA"
  obj <- Seurat::SCTransform(
    object = obj,
    assay = "RNA",
    new.assay.name = "SCT",
    vst.flavor = "v2",
    verbose = FALSE
  )
  DefaultAssay(obj) <- "SCT"
  log_progress("Finished SCTransform.")
}

log_progress("Running PCA on SCT assay.")
obj <- Seurat::RunPCA(
  object = obj,
  assay = "SCT",
  npcs = max(dims_use),
  verbose = FALSE,
  reduction.name = "pca_sct"
)
log_progress("Finished PCA.")

log_progress("Running pre-Harmony UMAP.")
obj <- Seurat::RunUMAP(
  object = obj,
  reduction = "pca_sct",
  dims = dims_use,
  reduction.name = "umap_pre_harmony",
  verbose = FALSE
)
log_progress("Finished pre-Harmony UMAP.")

log_progress("Running Harmony.")
obj <- harmony::RunHarmony(
  object = obj,
  group.by.vars = "sample",
  reduction = "pca_sct",
  assay.use = "SCT",
  reduction.save = "harmony",
  max_iter = 20L,
  early_stop = FALSE
)
log_progress("Finished Harmony.")

log_progress("Running Harmony UMAP.")
obj <- Seurat::RunUMAP(
  object = obj,
  reduction = "harmony",
  dims = dims_use,
  reduction.name = "umap",
  verbose = FALSE
)
log_progress("Finished Harmony UMAP.")

log_progress("Finding neighbors.")
obj <- Seurat::FindNeighbors(
  object = obj,
  reduction = "harmony",
  dims = dims_use,
  verbose = FALSE
)
log_progress("Finished finding neighbors.")

log_progress("Finding clusters for resolutions ", paste(resolution_labels, collapse = ", "), ".")
obj <- Seurat::FindClusters(
  object = obj,
  resolution = resolution_grid,
  verbose = FALSE
)
log_progress("Finished clustering.")

obj@meta.data <- normalize_cluster_columns(obj@meta.data, resolution_labels)

cluster_cols <- paste0("SCT_snn_res.", resolution_labels)
missing_cluster_cols <- setdiff(cluster_cols, colnames(obj@meta.data))
if (length(missing_cluster_cols) > 0) {
  stop(
    "Missing expected cluster metadata columns: ",
    paste(missing_cluster_cols, collapse = ", "),
    call. = FALSE
  )
}

log_progress("Preparing metadata and embedding exports.")
metadata_out <- obj@meta.data |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "cellid")
pca_out <- Seurat::Embeddings(obj, reduction = "pca_sct") |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "cellid")
harmony_out <- Seurat::Embeddings(obj, reduction = "harmony") |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "cellid")
umap_pre_out <- Seurat::Embeddings(obj, reduction = "umap_pre_harmony") |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "cellid")
umap_out <- Seurat::Embeddings(obj, reduction = "umap") |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "cellid")
pca_sd_out <- tibble::tibble(
  pc = seq_along(obj[["pca_sct"]]@stdev),
  stdev = as.numeric(obj[["pca_sct"]]@stdev)
)

log_progress("Writing clustree plot: ", clustree_plot)
grDevices::pdf(clustree_plot, width = 10, height = 8)
print(clustree::clustree(obj@meta.data, prefix = "SCT_snn_res."))
grDevices::dev.off()

log_progress("Writing metadata and embeddings.")
readr::write_csv(metadata_out, metadata_csv)
readr::write_csv(pca_out, pca_embeddings_csv)
readr::write_csv(harmony_out, harmony_embeddings_csv)
readr::write_csv(umap_pre_out, umap_pre_embeddings_csv)
readr::write_csv(umap_out, umap_embeddings_csv)
readr::write_csv(pca_sd_out, pca_sd_csv)
log_progress("Saving SCT/Harmony Seurat object: ", annotation_rds)
saveRDS(obj, file = annotation_rds)
log_progress("Writing session info: ", session_info_csv)
utils::write.csv(
  utils::capture.output(utils::sessionInfo()),
  file = session_info_csv,
  row.names = FALSE
)
log_progress("Done.")

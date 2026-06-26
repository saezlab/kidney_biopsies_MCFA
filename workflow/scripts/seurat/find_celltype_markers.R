if (exists("snakemake")) {
  rds_path <- snakemake@input[["rds"]]
  annotation_csv <- snakemake@input[["annotation_csv"]]
  markers_csv <- snakemake@output[["csv"]]
  cluster_column <- snakemake@params[["cluster_column"]]
  group_column <- snakemake@params[["group_column"]]
} else {
  rds_path <- "results/seurat/subset_high_quality/sct_harmony.rds"
  annotation_csv <- "config/RNA_level1_annotations.csv"
  markers_csv <- "results/seurat/subset_high_quality/markers/findallmarkers_level1.csv"
  cluster_column <- "SCT_snn_res.0.6"
  group_column <- "level1"
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(Seurat)
  library(tibble)
})

obj <- readRDS(rds_path)
annotation_df <- readr::read_csv(annotation_csv, show_col_types = FALSE)

if (!cluster_column %in% colnames(obj@meta.data)) {
  stop(
    "Missing expected cluster metadata column: ",
    cluster_column,
    call. = FALSE
  )
}

missing_annotation_cols <- setdiff(c(cluster_column, group_column), colnames(annotation_df))
if (length(missing_annotation_cols) > 0) {
  stop(
    "Annotation file is missing required columns: ",
    paste(missing_annotation_cols, collapse = ", "),
    call. = FALSE
  )
}

annotation_df <- annotation_df %>%
  dplyr::select(
    cluster = dplyr::all_of(cluster_column),
    group = dplyr::all_of(group_column)
  ) %>%
  dplyr::mutate(
    cluster = as.character(cluster),
    group = trimws(as.character(group)),
    group = dplyr::na_if(group, ""),
    group = dplyr::na_if(group, "NA")
  ) %>%
  dplyr::filter(!is.na(cluster), cluster != "") %>%
  dplyr::distinct(cluster, .keep_all = TRUE)

cell_groups <- tibble::tibble(
  cell = colnames(obj),
  cluster = as.character(obj[[cluster_column, drop = TRUE]])
) %>%
  dplyr::left_join(annotation_df, by = "cluster") %>%
  dplyr::filter(!is.na(group))

if (nrow(cell_groups) == 0L) {
  stop(
    "No cells had non-empty ",
    group_column,
    " annotations after mapping from ",
    cluster_column,
    ".",
    call. = FALSE
  )
}

obj <- obj[, cell_groups$cell]
cell_group_lookup <- stats::setNames(cell_groups$group, cell_groups$cell)
obj[[group_column]] <- unname(cell_group_lookup[colnames(obj)])

obj <- Seurat::PrepSCTFindMarkers(obj)
Seurat::Idents(obj) <- obj[[group_column, drop = TRUE]]
celltype_markers <- Seurat::FindAllMarkers(
  object = obj,
  assay = "SCT",
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1,
  verbose = TRUE
)

if (is.null(celltype_markers) || ncol(celltype_markers) == 0) {
  celltype_markers <- tibble::tibble()
}

readr::write_csv(celltype_markers, markers_csv)

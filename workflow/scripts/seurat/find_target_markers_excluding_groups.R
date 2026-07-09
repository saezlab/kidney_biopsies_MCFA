if (exists("snakemake")) {
  rds_path <- snakemake@input[["rds"]]
  metadata_csv <- snakemake@input[["metadata"]]
  markers_csv <- snakemake@output[["csv"]]
  group_column <- snakemake@params[["group_column"]]
  target_group <- snakemake@params[["target_group"]]
  exclude_groups <- unlist(snakemake@params[["exclude_groups"]])
} else {
  rds_path <- "results/seurat/subset_high_quality/sct_harmony.rds"
  metadata_csv <- "results/seurat/subset_high_quality/metadata_annotations.csv"
  markers_csv <- "results/seurat/subset_high_quality/markers/findmarkers_level2_dPT_excluding_PT_aPT.csv"
  group_column <- "level2"
  target_group <- "dPT"
  exclude_groups <- c("PT", "aPT")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(Seurat)
  library(tibble)
})

clean_label <- function(x) {
  x <- trimws(as.character(x))
  x <- dplyr::na_if(x, "")
  dplyr::na_if(x, "NA")
}

validate_columns <- function(tbl, required_cols, tbl_name) {
  missing_cols <- setdiff(required_cols, colnames(tbl))
  if (length(missing_cols) > 0L) {
    stop(
      tbl_name,
      " is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
}

obj <- readRDS(rds_path)
metadata <- readr::read_csv(metadata_csv, show_col_types = FALSE)

validate_columns(metadata, c("cellid", group_column), "Metadata CSV")

cell_groups <- metadata %>%
  dplyr::transmute(
    cell = as.character(.data$cellid),
    group = clean_label(.data[[group_column]])
  ) %>%
  dplyr::filter(!is.na(.data$cell), .data$cell != "", !is.na(.data$group)) %>%
  dplyr::filter(!.data$group %in% exclude_groups) %>%
  dplyr::distinct(.data$cell, .keep_all = TRUE) %>%
  dplyr::filter(.data$cell %in% colnames(obj))

if (nrow(cell_groups) == 0L) {
  stop("No annotated cells remained after applying exclusions and matching to the Seurat object.", call. = FALSE)
}

if (!target_group %in% cell_groups$group) {
  stop(
    "Target group '",
    target_group,
    "' was not found after applying exclusions.",
    call. = FALSE
  )
}

comparison_groups <- setdiff(unique(cell_groups$group), target_group)
if (length(comparison_groups) == 0L) {
  stop(
    "No comparison cells remained after excluding groups and selecting target group '",
    target_group,
    "'.",
    call. = FALSE
  )
}

obj <- obj[, cell_groups$cell]
cell_group_lookup <- stats::setNames(cell_groups$group, cell_groups$cell)
obj[[group_column]] <- unname(cell_group_lookup[colnames(obj)])

obj <- Seurat::PrepSCTFindMarkers(obj)
Seurat::Idents(obj) <- obj[[group_column, drop = TRUE]]

target_markers <- Seurat::FindMarkers(
  object = obj,
  ident.1 = target_group,
  ident.2 = NULL,
  assay = "SCT",
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1,
  recorrect_umi = FALSE,
  verbose = TRUE
)

if (is.null(target_markers) || ncol(target_markers) == 0L) {
  target_markers <- tibble::tibble()
} else {
  target_markers <- target_markers %>%
    tibble::rownames_to_column("gene")
}

readr::write_csv(target_markers, markers_csv)

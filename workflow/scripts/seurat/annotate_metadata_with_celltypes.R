if (exists("snakemake")) {
  metadata_csv <- snakemake@input[["metadata"]]
  annotation_csv <- snakemake@input[["annotation_csv"]]
  output_csv <- snakemake@output[["csv"]]
  cluster_column <- snakemake@params[["cluster_column"]]
} else {
  metadata_csv <- "results/seurat/subset_high_quality/metadata.csv"
  annotation_csv <- "config/RNA_level1_annotations.csv"
  output_csv <- "results/seurat/subset_high_quality/metadata_annotations.csv"
  cluster_column <- "SCT_snn_res.0.6"
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
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

metadata <- readr::read_csv(metadata_csv, show_col_types = FALSE)
annotations <- readr::read_csv(annotation_csv, show_col_types = FALSE)

validate_columns(metadata, c("cellid", cluster_column), "Metadata CSV")
validate_columns(annotations, c(cluster_column, "level1", "level2"), "Annotation CSV")

annotations_clean <- annotations %>%
  dplyr::transmute(
    cluster = as.character(.data[[cluster_column]]),
    level1 = clean_label(.data$level1),
    level2 = clean_label(.data$level2)
  )

metadata_annotated <- metadata %>%
  dplyr::mutate(cluster = as.character(.data[[cluster_column]])) %>%
  dplyr::left_join(annotations_clean, by = "cluster") %>%
  dplyr::filter(!is.na(.data$level1), !is.na(.data$level2)) %>%
  dplyr::select(-"cluster")

readr::write_csv(metadata_annotated, output_csv)

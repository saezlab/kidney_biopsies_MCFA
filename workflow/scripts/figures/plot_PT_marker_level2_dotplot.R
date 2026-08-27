if (exists("snakemake")) {
  rds_path <- snakemake@input[["rds"]]
  annotation_csv <- snakemake@input[["annotation_csv"]]
  pdf_path <- snakemake@output[["pdf"]]
  png_path <- snakemake@output[["png"]]
  genes <- unlist(snakemake@params[["genes"]])
  cluster_column <- snakemake@params[["cluster_column"]]
  level1_column <- snakemake@params[["level1_column"]]
  level2_column <- snakemake@params[["level2_column"]]
} else {
  rds_path <- "results/seurat/subset_high_quality/sct_harmony.rds"
  annotation_csv <- "config/RNA_level1_annotations.csv"
  pdf_path <- "plots/figures/Figure1Ext/PT_marker_level2_dotplot.pdf"
  png_path <- "plots/figures/Figure1Ext/PT_marker_level2_dotplot.png"
  genes <- c(
    "LRP2",
    "CUBN",
    "SLC5A12",
    "SLC22A6",
    "SLC22A7",
    "TPM1",
    "ITGB8",
    "PROM1",
    "SPP1",
    "CD24",
    "CIRBP",
    "SLC12A1",
    "SLC8A1"
  )
  cluster_column <- "SCT_snn_res.0.6"
  level1_column <- "level1"
  level2_column <- "level2"
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(Seurat)
  library(tibble)
})

level1_panel_order <- c("PT", "TAL", "DCT", "CNT_PC", "IC", "EC", "POD", "PEC", "Stromal", "Immune")
level2_plot_order <- c(
  "PT", "aPT", "dPT", "TAL", "DCT", "CNT", "PC",
  "IC-A", "IC-B", "EC", "POD", "PEC", "MES", "VSM-P", "FIB", "Immune"
)

clean_label <- function(x) {
  x <- trimws(as.character(x))
  x <- dplyr::na_if(x, "")
  dplyr::na_if(x, "NA")
}

genes <- unique(as.character(genes))
genes <- genes[!is.na(genes) & nzchar(genes)]
if (length(genes) == 0L) {
  stop("No genes were provided for the PT marker level2 dotplot.", call. = FALSE)
}

obj <- readRDS(rds_path)
annotation_df <- utils::read.csv(
  annotation_csv,
  fill = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_annotation_cols <- c(cluster_column, level1_column, level2_column)
missing_annotation_cols <- setdiff(required_annotation_cols, colnames(annotation_df))
if (length(missing_annotation_cols) > 0L) {
  stop(
    "Annotation file is missing required columns: ",
    paste(missing_annotation_cols, collapse = ", "),
    call. = FALSE
  )
}

if (!cluster_column %in% colnames(obj@meta.data)) {
  stop("Missing expected cluster metadata column: ", cluster_column, call. = FALSE)
}
if (!"SCT" %in% Seurat::Assays(obj)) {
  stop("The Seurat object is missing the required 'SCT' assay.", call. = FALSE)
}

annotation_df <- annotation_df %>%
  dplyr::select(
    cluster = dplyr::all_of(cluster_column),
    level1 = dplyr::all_of(level1_column),
    level2 = dplyr::all_of(level2_column)
  ) %>%
  dplyr::mutate(
    cluster = clean_label(.data$cluster),
    level1 = clean_label(.data$level1),
    level2 = clean_label(.data$level2)
  ) %>%
  dplyr::filter(!is.na(.data$cluster), !is.na(.data$level1), !is.na(.data$level2)) %>%
  dplyr::distinct(.data$cluster, .keep_all = TRUE)

observed_level1 <- unique(annotation_df$level1)
observed_level2 <- unique(annotation_df$level2)
if (length(observed_level2) == 0L) {
  stop("No complete level2 annotations were found.", call. = FALSE)
}

missing_level1_order <- setdiff(observed_level1, level1_panel_order)
if (length(missing_level1_order) > 0L) {
  stop(
    "Missing level1 panel order entries for: ",
    paste(missing_level1_order, collapse = ", "),
    call. = FALSE
  )
}

missing_level2_order <- setdiff(observed_level2, level2_plot_order)
if (length(missing_level2_order) > 0L) {
  stop(
    "Missing level2 plot order entries for: ",
    paste(missing_level2_order, collapse = ", "),
    call. = FALSE
  )
}

level2_levels <- level2_plot_order[level2_plot_order %in% observed_level2]
level2_axis_levels <- rev(level2_levels)

cell_groups <- tibble::tibble(
  cell = colnames(obj),
  cluster = clean_label(obj[[cluster_column, drop = TRUE]])
) %>%
  dplyr::left_join(annotation_df, by = "cluster") %>%
  dplyr::filter(!is.na(.data$level1), !is.na(.data$level2))

if (nrow(cell_groups) == 0L) {
  stop("No cells had complete level1/level2 annotations after cluster mapping.", call. = FALSE)
}

missing_genes <- setdiff(genes, rownames(obj[["SCT"]]))
if (length(missing_genes) > 0L) {
  stop(
    "Requested gene(s) were not found in the Seurat object: ",
    paste(missing_genes, collapse = ", "),
    call. = FALSE
  )
}

obj <- obj[, cell_groups$cell]
level2_lookup <- stats::setNames(cell_groups$level2, cell_groups$cell)
obj[[level2_column]] <- factor(
  unname(level2_lookup[colnames(obj)]),
  levels = level2_axis_levels
)

p_dot <- Seurat::DotPlot(
  object = obj,
  features = genes,
  assay = "SCT",
  group.by = level2_column,
  dot.min = 0,
  dot.scale = 5,
  scale = TRUE,
  cols = c("lightgrey", "darkblue")
) +
  Seurat::RotatedAxis() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 9),
    axis.text.y = ggplot2::element_text(size = 9),
    axis.title.x = ggplot2::element_text(size = 10),
    axis.title.y = ggplot2::element_text(size = 10),
    strip.text = ggplot2::element_text(size = 10)
  )

plot_width <- max(10, 0.55 * length(genes))
plot_height <- max(8, 0.35 * length(level2_axis_levels))

if (!exists("snakemake")) {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
}

ggplot2::ggsave(pdf_path, p_dot, width = plot_width, height = plot_height, units = "in")
ggplot2::ggsave(png_path, p_dot, width = plot_width, height = plot_height, units = "in", dpi = 300)

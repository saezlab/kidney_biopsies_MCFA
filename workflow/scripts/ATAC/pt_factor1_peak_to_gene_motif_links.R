if (exists("snakemake")) {
  peak_to_gene_path <- snakemake@input[["peak_to_gene"]]
  motif_peak_matches_path <- snakemake@input[["motif_peak_matches"]]
  loadings_path <- snakemake@input[["loadings"]]
  out_csv <- snakemake@output[["csv"]]
  summary_csv <- snakemake@output[["summary"]]
} else {
  peak_to_gene_path <- "results/ATAC/circe_peak_to_gene.csv"
  motif_peak_matches_path <- "results/ATAC/motif_peak_matches.csv"
  loadings_path <- "results/MCFA/loadings~hvg~5.csv"
  out_csv <- "results/ATAC/pt_factor1_peak_to_gene_motif_links.csv"
  summary_csv <- "results/ATAC/pt_factor1_peak_to_gene_motif_links_summary.csv"
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

read_csv_with_optional_index <- function(path) {
  tbl <- readr::read_csv(path, show_col_types = FALSE, name_repair = "minimal")
  first_name <- names(tbl)[1]
  if (first_name %in% c("", "X", "...1")) {
    tbl <- dplyr::select(tbl, -1)
  }
  tbl
}

require_columns <- function(tbl, required_columns, file_label) {
  missing_columns <- setdiff(required_columns, names(tbl))
  if (length(missing_columns) > 0L) {
    stop(
      file_label,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

peak_to_gene <- read_csv_with_optional_index(peak_to_gene_path)
motif_peak_matches <- read_csv_with_optional_index(motif_peak_matches_path)
loadings <- read_csv_with_optional_index(loadings_path)

require_columns(
  peak_to_gene,
  c("peak_id", "gene_short_name", "coaccess"),
  "Peak-to-gene table"
)
require_columns(
  motif_peak_matches,
  c(
    "direction", "open_group", "motif", "mlog10Padj", "Enrichment", "CompareProportion",
    "BackgroundProporition", "seqnames", "start", "end", "Log2FC",
    "FDR", "Pval", "Mean", "MeanDiff", "AUC"
  ),
  "Motif peak match table"
)
require_columns(
  loadings,
  c("view", "variable", "Factor", "Loading"),
  "MCFA loadings table"
)

peak_to_gene <- peak_to_gene %>%
  dplyr::transmute(
    peak_id = as.character(.data$peak_id),
    gene_short_name = as.character(.data$gene_short_name),
    coaccess = as.numeric(.data$coaccess)
  )

motif_peak_matches <- motif_peak_matches %>%
  dplyr::mutate(
    peak_id = paste(.data$seqnames, .data$start, .data$end, sep = "_"),
    Log2FC = as.numeric(.data$Log2FC),
    FDR = as.numeric(.data$FDR)
  )

pt_factor1_loadings <- loadings %>%
  dplyr::filter(.data$view == "PT", .data$Factor == "Factor1") %>%
  dplyr::transmute(
    gene_short_name = as.character(.data$variable),
    pt_factor1_loading = as.numeric(.data$Loading)
  )

if (nrow(pt_factor1_loadings) == 0L) {
  stop("No PT Factor1 loadings found in ", loadings_path, ".", call. = FALSE)
}

joined <- peak_to_gene %>%
  dplyr::inner_join(
    motif_peak_matches,
    by = "peak_id",
    relationship = "many-to-many"
  ) %>%
  dplyr::inner_join(
    pt_factor1_loadings,
    by = "gene_short_name",
    relationship = "many-to-many"
  ) %>%
  dplyr::select(
    "direction",
    "open_group",
    "motif",
    "mlog10Padj",
    "Enrichment",
    "CompareProportion",
    "BackgroundProporition",
    "peak_id",
    "seqnames",
    "start",
    "end",
    peak_log2fc = "Log2FC",
    peak_fdr = "FDR",
    peak_pval = "Pval",
    peak_mean = "Mean",
    peak_mean_diff = "MeanDiff",
    peak_auc = "AUC",
    "gene_short_name",
    "coaccess",
    "pt_factor1_loading"
  ) %>%
  dplyr::arrange(.data$motif, .data$peak_id, .data$gene_short_name)

summary <- tibble::tibble(
  motif_match_rows = nrow(motif_peak_matches),
  motif_matches_fdr_le_0_1_log2fc_lt_0 = sum(
    motif_peak_matches$FDR <= 0.1 & motif_peak_matches$Log2FC < 0,
    na.rm = TRUE
  ),
  motif_matches_fdr_le_0_1_log2fc_gt_0 = sum(
    motif_peak_matches$FDR <= 0.1 & motif_peak_matches$Log2FC > 0,
    na.rm = TRUE
  ),
  joined_rows = nrow(joined),
  distinct_motifs = dplyr::n_distinct(joined$motif),
  distinct_peaks = dplyr::n_distinct(joined$peak_id),
  distinct_genes = dplyr::n_distinct(joined$gene_short_name)
)

readr::write_csv(joined, out_csv)
readr::write_csv(summary, summary_csv)

cat("Joined rows: ", nrow(joined), "\n", sep = "")
cat(
  "Negative Log2FC motif matches at FDR <= 0.1: ",
  summary$motif_matches_fdr_le_0_1_log2fc_lt_0,
  "\n",
  sep = ""
)

if (exists("snakemake")) {
  project_path <- snakemake@input[["project_path"]]
  marker_peaks_rds <- snakemake@input[["marker_peaks"]]
  out_csv <- snakemake@output[["csv"]]
  threads <- snakemake@threads
} else {
  project_path <- "data/ATAC/proj_peaks2"
  marker_peaks_rds <- "results/ATAC/PT_aPT_differential_accessibility_marker_peaks.rds"
  out_csv <- "results/ATAC/PT_aPT_differential_accessibility_marker_peaks.csv"
  threads <- 1L
}

suppressPackageStartupMessages({
  library(ArchR)
  library(SummarizedExperiment)
  library(readr)
})

threads <- max(1L, as.integer(threads))

cat("project_path: ", project_path, "\n", sep = "")
cat("marker_peaks_rds: ", marker_peaks_rds, "\n", sep = "")
cat("out_csv: ", out_csv, "\n", sep = "")
cat("threads: ", threads, "\n", sep = "")

addArchRGenome("hg38")
addArchRThreads(threads = threads)

proj <- loadArchRProject(path = project_path)
marker_peaks <- readRDS(marker_peaks_rds)
peaks <- getPeakSet(proj)

if (length(peaks) != nrow(marker_peaks)) {
  stop(
    "ArchR peak set length (",
    length(peaks),
    ") does not match marker peak rows (",
    nrow(marker_peaks),
    ").",
    call. = FALSE
  )
}

required_assays <- c("Log2FC", "FDR", "Pval")
missing_assays <- setdiff(required_assays, assayNames(marker_peaks))
if (length(missing_assays) > 0L) {
  stop(
    "Marker peaks object is missing required assay(s): ",
    paste(missing_assays, collapse = ", "),
    call. = FALSE
  )
}

peak_table <- data.frame(
  peak_idx = seq_along(peaks),
  seqnames = as.character(seqnames(peaks)),
  start = start(peaks),
  end = end(peaks),
  Log2FC = as.numeric(assay(marker_peaks, "Log2FC")[, 1]),
  FDR = as.numeric(assay(marker_peaks, "FDR")[, 1]),
  Pval = as.numeric(assay(marker_peaks, "Pval")[, 1]),
  stringsAsFactors = FALSE
)
peak_table$peak_id <- paste(
  peak_table$seqnames,
  peak_table$start,
  peak_table$end,
  sep = "_"
)

if (anyDuplicated(peak_table$peak_id) > 0L) {
  stop("Exported peak IDs are not unique.", call. = FALSE)
}

readr::write_csv(peak_table, out_csv, na = "NA")

cat("Exported peak rows: ", nrow(peak_table), "\n", sep = "")

if (exists("snakemake")) {
  project_path <- snakemake@input[["project_path"]]
  marker_peaks_rds <- snakemake@input[["marker_peaks"]]
  apt_open_enrich_motifs_rds <- snakemake@input[["apt_open_enrich_motifs"]]
  pt_open_enrich_motifs_rds <- snakemake@input[["pt_open_enrich_motifs"]]
  out_rds <- snakemake@output[["rds"]]
  out_csv <- snakemake@output[["csv"]]
  threads <- snakemake@threads
} else {
  project_path <- "data/ATAC/proj_peaks2"
  marker_peaks_rds <- "results/ATAC/PT_aPT_differential_accessibility_marker_peaks.rds"
  apt_open_enrich_motifs_rds <- "results/ATAC/aPT_open_motif_enrichment.rds"
  pt_open_enrich_motifs_rds <- "results/ATAC/PT_open_motif_enrichment.rds"
  out_rds <- "results/ATAC/motif_peak_matches.rds"
  out_csv <- "results/ATAC/motif_peak_matches.csv"
  threads <- 1L
}

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(SummarizedExperiment)
  library(Matrix)
})

threads <- max(1L, as.integer(threads))

cat("project_path: ", project_path, "\n", sep = "")
cat("marker_peaks_rds: ", marker_peaks_rds, "\n", sep = "")
cat("apt_open_enrich_motifs_rds: ", apt_open_enrich_motifs_rds, "\n", sep = "")
cat("pt_open_enrich_motifs_rds: ", pt_open_enrich_motifs_rds, "\n", sep = "")
cat("out_rds: ", out_rds, "\n", sep = "")
cat("out_csv: ", out_csv, "\n", sep = "")
cat("threads: ", threads, "\n", sep = "")

addArchRGenome("hg38")
addArchRThreads(threads = threads)

proj <- loadArchRProject(path = project_path)

peaks <- getPeakSet(proj)

cat("Peak set head:\n")
print(head(peaks))
cat("Peak set length: ", length(peaks), "\n", sep = "")

motif_matches <- getMatches(ArchRProj = proj, name = "Motif")

cat("Motif match assays:\n")
print(assayNames(motif_matches))

motif_match_matrix <- assay(motif_matches)
cat("Motif match matrix class:\n")
print(class(motif_match_matrix))
print(dim(motif_match_matrix))

marker_peaks <- readRDS(marker_peaks_rds)
apt_open_enrich_motifs <- readRDS(apt_open_enrich_motifs_rds)
pt_open_enrich_motifs <- readRDS(pt_open_enrich_motifs_rds)

cat("Marker peaks object:\n")
print(class(marker_peaks))
print(dim(marker_peaks))
print(assayNames(marker_peaks))
cat("Marker peaks rowData head:\n")
print(head(rowData(marker_peaks)))

cat("aPT-open enriched motifs object:\n")
print(class(apt_open_enrich_motifs))
print(dim(apt_open_enrich_motifs))
print(assayNames(apt_open_enrich_motifs))
cat("aPT-open enriched motif rownames head:\n")
print(head(rownames(apt_open_enrich_motifs)))

cat("PT-open enriched motifs object:\n")
print(class(pt_open_enrich_motifs))
print(dim(pt_open_enrich_motifs))
print(assayNames(pt_open_enrich_motifs))
cat("PT-open enriched motif rownames head:\n")
print(head(rownames(pt_open_enrich_motifs)))

motif_padj_max <- 0.05
motif_mlog10padj_min <- -log10(motif_padj_max)
motif_enrichment_min <- 1
peak_fdr_max <- 0.1
peak_log2fc_min <- 0.5

peak_df <- data.frame(
  peak_idx = seq_along(peaks),
  seqnames = as.character(seqnames(peaks)),
  start = start(peaks),
  end = end(peaks),
  stringsAsFactors = FALSE
)

cat("DA peak table head:\n")
print(head(peak_df[da_peak_idx, , drop = FALSE]))

for (assay_name in intersect(
  c("Log2FC", "FDR", "Pval", "Mean", "MeanDiff", "AUC"),
  assayNames(marker_peaks)
)) {
  peak_df[[assay_name]] <- as.vector(assay(marker_peaks, assay_name)[, 1])
}

motif_df <- data.frame(
  motif = character(),
  stringsAsFactors = FALSE
)

build_motif_table <- function(enrich_motifs) {
  motif_table <- data.frame(
    motif = rownames(enrich_motifs),
    stringsAsFactors = FALSE
  )

  for (assay_name in intersect(
    c("mlog10Padj", "Enrichment", "CompareProportion", "BackgroundProporition"),
    assayNames(enrich_motifs)
  )) {
    motif_table[[assay_name]] <- as.vector(assay(enrich_motifs, assay_name)[, 1])
  }

  motif_table
}

match_direction <- function(enrich_motifs, direction, open_group, peak_idx) {
  motif_table <- build_motif_table(enrich_motifs)
  selected_motifs <- motif_table$motif[
    motif_table$mlog10Padj >= motif_mlog10padj_min &
      motif_table$Enrichment >= motif_enrichment_min
  ]

  cat(direction, " selected motifs: ", length(selected_motifs), "\n", sep = "")
  print(head(selected_motifs, 20))
  cat(direction, " DA peak indices: ", length(peak_idx), "\n", sep = "")
  print(head(peak_idx, 20))

  if (length(selected_motifs) == 0L || length(peak_idx) == 0L) {
    return(NULL)
  }

  missing_motifs <- setdiff(selected_motifs, colnames(motif_match_matrix))
  if (length(missing_motifs) > 0L) {
    stop(
      "Selected enriched motif(s) are missing from motif match matrix: ",
      paste(missing_motifs, collapse = ", "),
      call. = FALSE
    )
  }

  match_tables <- lapply(selected_motifs, function(motif) {
    motif_hit_idx <- Matrix::which(
      motif_match_matrix[, motif, drop = FALSE] != 0,
      arr.ind = TRUE
    )

    motif_peak_idx <- motif_hit_idx[, "row"]
    matched_idx <- intersect(peak_idx, motif_peak_idx)
    if (length(matched_idx) == 0L) {
      return(NULL)
    }

    motif_stats <- motif_table[match(motif, motif_table$motif), , drop = FALSE]
    cbind(
      data.frame(
        direction = direction,
        open_group = open_group,
        stringsAsFactors = FALSE
      ),
      motif_stats[rep(1L, length(matched_idx)), , drop = FALSE],
      peak_df[matched_idx, , drop = FALSE],
      row.names = NULL
    )
  })

  do.call(rbind, match_tables)
}

apt_open_peak_idx <- which(
  as.vector(assay(marker_peaks, "FDR")[, 1]) <= peak_fdr_max &
    as.vector(assay(marker_peaks, "Log2FC")[, 1]) >= peak_log2fc_min
)

pt_open_peak_idx <- which(
  as.vector(assay(marker_peaks, "FDR")[, 1]) <= peak_fdr_max &
    as.vector(assay(marker_peaks, "Log2FC")[, 1]) <= -peak_log2fc_min
)

matches_df <- do.call(
  rbind,
  list(
    match_direction(
      apt_open_enrich_motifs,
      direction = "aPT_open",
      open_group = "pt_C6",
      peak_idx = apt_open_peak_idx
    ),
    match_direction(
      pt_open_enrich_motifs,
      direction = "PT_open",
      open_group = "pt_C3",
      peak_idx = pt_open_peak_idx
    )
  )
)

cat("Final motif-peak match table:\n")
print(dim(matches_df))
print(head(matches_df))

if (is.null(matches_df)) {
  matches_df <- data.frame()
}

output <- list(
  matches = matches_df,
  apt_open_peak_idx = apt_open_peak_idx,
  pt_open_peak_idx = pt_open_peak_idx,
  thresholds = list(
    motif_padj_max = motif_padj_max,
    motif_mlog10padj_min = motif_mlog10padj_min,
    motif_enrichment_min = motif_enrichment_min,
    peak_fdr_max = peak_fdr_max,
    peak_log2fc_min = peak_log2fc_min
  )
)

saveRDS(output, file = out_rds)
write.csv(matches_df, file = out_csv, row.names = FALSE)

cat("aPT-open DA peaks: ", length(apt_open_peak_idx), "\n", sep = "")
cat("PT-open DA peaks: ", length(pt_open_peak_idx), "\n", sep = "")
cat("Motif-peak matches: ", nrow(matches_df), "\n", sep = "")

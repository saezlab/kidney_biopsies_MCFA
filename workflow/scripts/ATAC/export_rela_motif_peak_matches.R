if (exists("snakemake")) {
  project_path <- snakemake@input[["project_path"]]
  out_csv <- snakemake@output[["csv"]]
  motif_prefix <- snakemake@params[["motif_prefix"]]
  threads <- snakemake@threads
} else {
  project_path <- "data/ATAC/proj_peaks2"
  out_csv <- "results/ATAC/rela_motif_peak_matches.csv"
  motif_prefix <- "RELA"
  threads <- 1L
}

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
  library(readr)
})

threads <- max(1L, as.integer(threads))

cat("project_path: ", project_path, "\n", sep = "")
cat("out_csv: ", out_csv, "\n", sep = "")
cat("motif_prefix: ", motif_prefix, "\n", sep = "")
cat("threads: ", threads, "\n", sep = "")

addArchRGenome("hg38")
addArchRThreads(threads = threads)

proj <- loadArchRProject(path = project_path)
peaks <- getPeakSet(proj)
motif_matches <- getMatches(ArchRProj = proj, name = "Motif")
motif_match_matrix <- assay(motif_matches)

if (nrow(motif_match_matrix) != length(peaks)) {
  stop(
    "Motif match matrix rows (",
    nrow(motif_match_matrix),
    ") do not match ArchR peak set length (",
    length(peaks),
    ").",
    call. = FALSE
  )
}

motif_pattern <- paste0("^", motif_prefix, "_")
rela_motifs <- grep(motif_pattern, colnames(motif_match_matrix), value = TRUE)

if (length(rela_motifs) == 0L) {
  stop(
    "No motif columns matched prefix pattern: ",
    motif_pattern,
    call. = FALSE
  )
}

motif_hit_idx <- Matrix::which(
  motif_match_matrix[, rela_motifs, drop = FALSE] != 0,
  arr.ind = TRUE
)

if (nrow(motif_hit_idx) == 0L) {
  rela_peak_table <- data.frame(
    peak_id = character(),
    motif = character(),
    seqnames = character(),
    start = integer(),
    end = integer(),
    peak_idx = integer(),
    stringsAsFactors = FALSE
  )
} else {
  peak_idx <- motif_hit_idx[, "row"]
  motif_idx <- motif_hit_idx[, "col"]
  rela_peak_table <- data.frame(
    peak_id = paste(
      as.character(seqnames(peaks))[peak_idx],
      start(peaks)[peak_idx],
      end(peaks)[peak_idx],
      sep = "_"
    ),
    motif = rela_motifs[motif_idx],
    seqnames = as.character(seqnames(peaks))[peak_idx],
    start = start(peaks)[peak_idx],
    end = end(peaks)[peak_idx],
    peak_idx = peak_idx,
    stringsAsFactors = FALSE
  )
}

rela_peak_table <- rela_peak_table[order(rela_peak_table$peak_idx, rela_peak_table$motif), ]

readr::write_csv(rela_peak_table, out_csv)

cat("RELA motif columns: ", paste(rela_motifs, collapse = ", "), "\n", sep = "")
cat("RELA motif peak rows: ", nrow(rela_peak_table), "\n", sep = "")
cat("Distinct RELA motif peaks: ", length(unique(rela_peak_table$peak_id)), "\n", sep = "")

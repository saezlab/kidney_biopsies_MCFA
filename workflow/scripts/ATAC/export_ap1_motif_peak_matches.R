if (exists("snakemake")) {
  project_path <- snakemake@input[["project_path"]]
  out_csv <- snakemake@output[["csv"]]
  motif_prefixes <- unlist(snakemake@params[["motif_prefixes"]])
  threads <- snakemake@threads
} else {
  project_path <- "data/ATAC/proj_peaks2"
  out_csv <- "results/ATAC/ap1_motif_peak_matches.csv"
  motif_prefixes <- c("JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2", "BATF")
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
cat("motif_prefixes: ", paste(motif_prefixes, collapse = ", "), "\n", sep = "")
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

motif_pattern <- paste0("^(", paste(motif_prefixes, collapse = "|"), ")_")
ap1_motifs <- grep(motif_pattern, colnames(motif_match_matrix), value = TRUE)

if (length(ap1_motifs) == 0L) {
  stop(
    "No motif columns matched prefix pattern: ",
    motif_pattern,
    call. = FALSE
  )
}

motif_hit_idx <- Matrix::which(
  motif_match_matrix[, ap1_motifs, drop = FALSE] != 0,
  arr.ind = TRUE
)

if (nrow(motif_hit_idx) == 0L) {
  ap1_peak_table <- data.frame(
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
  ap1_peak_table <- data.frame(
    peak_id = paste(
      as.character(seqnames(peaks))[peak_idx],
      start(peaks)[peak_idx],
      end(peaks)[peak_idx],
      sep = "_"
    ),
    motif = ap1_motifs[motif_idx],
    seqnames = as.character(seqnames(peaks))[peak_idx],
    start = start(peaks)[peak_idx],
    end = end(peaks)[peak_idx],
    peak_idx = peak_idx,
    stringsAsFactors = FALSE
  )
}

ap1_peak_table <- ap1_peak_table[order(ap1_peak_table$peak_idx, ap1_peak_table$motif), ]

readr::write_csv(ap1_peak_table, out_csv)

cat("AP-1 motif columns: ", paste(ap1_motifs, collapse = ", "), "\n", sep = "")
cat("AP-1 motif peak rows: ", nrow(ap1_peak_table), "\n", sep = "")
cat("Distinct AP-1 motif peaks: ", length(unique(ap1_peak_table$peak_id)), "\n", sep = "")

# (C) Sergio Alejandro Gómez
# Helpers of the accessibility scripts (02, 03, 07, 08): the 10 kidney populations, ArchR peak names and the
# population x peak count matrices built from the cluster x library pseudobulk.

suppressPackageStartupMessages({ library(Matrix); library(GenomicRanges) })

POPULATIONS <- c("POD", "EC", "Immune", "IC", "TAL", "CNT_PC", "DCT", "PT", "PEC", "Stromal")
DISPLAY_ORDER <- c("PT", "PEC", "TAL", "DCT", "CNT_PC", "IC", "POD", "EC", "Stromal", "Immune")
DISPLAY_LABEL <- c(PT = "PT", PEC = "PEC", TAL = "TAL", DCT = "DCT", CNT_PC = "CNT/PC", IC = "IC",
                   POD = "Podocyte", EC = "Endothelial", Stromal = "Stromal", Immune = "Immune")
STANDARD_CHROMOSOMES <- c(paste0("chr", 1:22), "chrX", "chrY")
PEAK_UNIVERSE_LABEL <- "ArchR_reproducible_peakset_proj_peaks2"


COHORT_EXCLUDED_GROUPS <- c("BKVAN", "KIN", "control_AKI")

parse_archr_peaks <- function(peak_names) {
  m <- regmatches(peak_names, regexec("^(.+)_([0-9]+)_([0-9]+)$", peak_names))
  stopifnot(all(lengths(m) == 4))
  data.frame(chr = vapply(m, `[`, "", 2), start = as.integer(vapply(m, `[`, "", 3)),
             end = as.integer(vapply(m, `[`, "", 4)), peak = peak_names, stringsAsFactors = FALSE)
}

peak_granges <- function(peaks) {
  pk <- peaks[peaks$is_standard_chr, ]
  gr <- GRanges(pk$chr, IRanges(pk$start, pk$end)); names(gr) <- pk$peak
  gr
}

to_cpm <- function(m, scale = 1e6) sweep(m, 2, colSums(m) / scale, "/")

population_counts <- function(pb, groups, exclude_disease_groups = COHORT_EXCLUDED_GROUPS, exclude_clusters = character(0)) {
  stopifnot(identical(colnames(pb), groups$group))
  keep <- !(groups$diseaseNames %in% exclude_disease_groups) & !(groups$cluster %in% exclude_clusters)
  labs <- sort(unique(groups$level1[keep]))
  ind <- sparseMatrix(i = seq_len(sum(keep)), j = match(groups$level1[keep], labs), x = 1, dims = c(sum(keep), length(labs)))
  cnt <- as.matrix(pb[, keep, drop = FALSE] %*% ind); dimnames(cnt) <- list(rownames(pb), labs)
  storage.mode(cnt) <- "double"
  list(counts = cnt, groups_kept = groups[keep, ])
}

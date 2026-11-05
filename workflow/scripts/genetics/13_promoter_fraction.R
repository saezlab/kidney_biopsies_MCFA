# (C) Sergio Alejandro Gómez
# Aim: Describe how promoter-proximal the peaks of the QTL comparison are.
# Inputs: peak universe, specificity matrix, Ensembl 108 gene annotation.
# Outputs: one row per population with the promoter-proximal percentages of specific and other peaks.
# Usage: Rscript workflow/scripts/genetics/13_promoter_fraction.R --peak_universe=results/genetics/with_ptC5/objects/peak_universe.rds --specificity=results/genetics/with_ptC5/objects/specificity_matrix.rds --gene_annotation=results/genetics/qtl_resources/gene_annotation_hg38_ensembl108.rds --out=results/genetics/with_ptC5/tables/promoter_fraction_by_population.tsv

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "qtl.R"))

args <- parse_args()
promoter_window_bp <- as.integer(arg(args, "promoter_window_bp", "2000"))

peaks <- readRDS(arg(args, "peak_universe"))
spec <- readRDS(arg(args, "specificity"))
stopifnot(identical(rownames(spec), peaks$peak), setequal(colnames(spec), POPULATIONS))

# Promoter-proximal: the peak overlaps TSS +/- window of any protein-coding transcript.
pk <- parse_archr_peaks(peaks$peak)
peaks_gr <- GRanges(pk$chr, IRanges(pk$start, pk$end))
tss_pc <- protein_coding_tss(readRDS(arg(args, "gene_annotation")))
win <- GRanges(seqnames(tss_pc), IRanges(pmax(1L, start(tss_pc) - promoter_window_bp), start(tss_pc) + promoter_window_bp))
prom <- countOverlaps(peaks_gr, win, ignore.strand = TRUE) > 0

specific_to_none <- rowSums(spec) == 0
out <- do.call(rbind, lapply(POPULATIONS, function(p) {
  s <- spec[, p]
  data.frame(population = p, n_specific = sum(s), n_other = sum(!s),
             pct_promoter_specific = 100 * mean(prom[s]), pct_promoter_other = 100 * mean(prom[!s]),
             ratio_other_vs_specific = mean(prom[!s]) / mean(prom[s]),
             pct_other_specific_to_no_population = 100 * mean(specific_to_none[!s]), stringsAsFactors = FALSE)
}))
out <- out[order(-out$n_specific), ]
log_msg(sprintf("%d of %d peaks promoter-proximal (%.1f%%); peaks specific to no population %.1f%%",
                sum(prom), length(prom), 100 * mean(prom), 100 * mean(specific_to_none)))
write_tsv(out, arg(args, "out"))

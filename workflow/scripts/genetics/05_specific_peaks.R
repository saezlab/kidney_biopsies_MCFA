# (C) Sergio Alejandro Gómez
# Aim: Define cell-type-specific peaks. A peak is specific to population c when CPM_c > 0 and CPM_c >= 2 x the mean CPM of
#   the other nine populations (CPM without pseudocount).
# Inputs: eGFR enrichment table (04), peak universe and CPM matrix (03).
# Outputs: specific-peak summary table; peak sets of the enriched populations and the full specificity matrix (.rds).
# Usage: Rscript workflow/scripts/genetics/05_specific_peaks.R --egfr_enrichment=results/genetics/with_ptC5/tables/egfr_enrichment.tsv --peak_universe=results/genetics/with_ptC5/objects/peak_universe.rds --cpm=results/genetics/with_ptC5/objects/accessibility_cpm.rds --include_ptC5=true --out_summary=results/genetics/with_ptC5/tables/specific_peaks_summary.tsv --out_specific_peaks=results/genetics/with_ptC5/objects/specific_peaks.rds --out_specificity=results/genetics/with_ptC5/objects/specificity_matrix.rds

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "accessibility.R"))

args <- parse_args()
version <- if (arg_true(args, "include_ptC5", NULL)) "with_ptC5" else "without_ptC5"
fold <- as.numeric(arg(args, "specificity_fold", "2"))
fdr <- as.numeric(arg(args, "fdr", "0.05"))

egfr <- read_tsv(arg(args, "egfr_enrichment"))
enriched <- egfr$cell_type[egfr$z > 0 & egfr$FDR < fdr]
log_msg("populations enriched for eGFR: ", paste(enriched, collapse = ", "))

peaks <- readRDS(arg(args, "peak_universe")); pk <- peaks[peaks$is_standard_chr, ]
cpm <- readRDS(arg(args, "cpm"))[pk$peak, ]   # stored (alphabetical) column order
stopifnot(setequal(colnames(cpm), POPULATIONS))
universe_text <- sprintf("%s ArchR peaks (%s)", format(nrow(cpm), big.mark = ","), version)
rule <- sprintf("CPM_c > 0 AND CPM_c >= %g x mean(CPM other 9 populations); CPM without pseudocount; %s", fold, universe_text)

# A peak open only in c (mean of the others 0) is specific to c; the inequality is evaluated directly to avoid 0/0.
spec <- matrix(FALSE, nrow(cpm), ncol(cpm), dimnames = dimnames(cpm))
for (p in colnames(cpm)) {
  mean_other <- rowMeans(cpm[, setdiff(colnames(cpm), p), drop = FALSE])
  spec[, p] <- (cpm[, p] > 0) & (cpm[, p] >= fold * mean_other)
}

summ <- do.call(rbind, lapply(enriched, function(p) {
  s <- spec[, p]; mean_other <- rowMeans(cpm[, setdiff(colnames(cpm), p), drop = FALSE])
  data.frame(cell_type = p, gchromvar_z = egfr$z[egfr$cell_type == p], gchromvar_FDR = egfr$FDR[egfr$cell_type == p],
             n_specific_peaks_cpm = sum(s), pct_of_universe = round(100 * sum(s) / nrow(cpm), 2),
             n_specific_exclusive_among_enriched = sum(s & rowSums(spec[, enriched, drop = FALSE]) == 1),
             n_specific_with_zero_mean_others = sum(s & mean_other == 0), median_cpm_in_specific = median(cpm[s, p]), stringsAsFactors = FALSE)
}))
summ$rule <- rule
summ$n_union_all_enriched <- sum(rowSums(spec[, enriched, drop = FALSE]) > 0)
summ$ptc5_version <- version
write_tsv(summ, arg(args, "out_summary"))

peak_sets <- setNames(lapply(enriched, function(p) pk$peak[spec[, p]]), enriched)
save_rds(list(enriched_populations = enriched, peak_sets = peak_sets, rule = rule, fold = fold, ptc5_version = version),
         arg(args, "out_specific_peaks"))
save_rds(spec, arg(args, "out_specificity"))

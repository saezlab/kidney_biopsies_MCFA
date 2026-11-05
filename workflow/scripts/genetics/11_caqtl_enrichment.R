# (C) Sergio Alejandro Gómez
# Aim: Test whether cell-type-specific peaks are enriched for kidney chromatin-accessibility QTLs. 
# Inputs: annotated peak table from the eQTL enrichment, bASA/RASQUAL target peaks and unique caQTL variants (hg38),
#   eGFR enrichment table, eQTL enrichment table.
# Outputs: enrichment table, model coefficients, peak counts, and the annotated peak table with caQTL status.
# Usage: Rscript workflow/scripts/genetics/11_caqtl_enrichment.R --peak_annotation=results/genetics/with_ptC5/objects/peak_annotation_eqtl.rds --caqtl_peaks=results/genetics/qtl_resources/basa_caqtl_peaks_hg38.rds --caqtl_variants=results/genetics/qtl_resources/basa_unique_caqtl_variants_hg38.rds --egfr_enrichment=results/genetics/with_ptC5/tables/egfr_enrichment.tsv --eqtl_enrichment=results/genetics/with_ptC5/tables/eqtl_enrichment.tsv --out_enrichment=results/genetics/with_ptC5/tables/caqtl_enrichment.tsv --out_coefficients=results/genetics/with_ptC5/tables/caqtl_enrichment_coefficients.tsv --out_counts=results/genetics/with_ptC5/tables/caqtl_peak_counts.tsv --out_peak_annotation=results/genetics/with_ptC5/objects/peak_annotation_caqtl.rds

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "qtl.R"))

args <- parse_args()
cluster_bin_bp <- as.numeric(arg(args, "cluster_bin_bp", "1e6"))
fdr <- as.numeric(arg(args, "fdr", "0.05"))

ann <- as.data.table(readRDS(arg(args, "peak_annotation")))
stopifnot(all(c("mean_accessibility", "gc_content", "distance_to_tss", "eqtl_positive") %in% names(ann)))
peaks_gr <- GRanges(ann$chr, IRanges(ann$start, ann$end)); names(peaks_gr) <- ann$peak
pgr <- readRDS(arg(args, "caqtl_peaks"))
vgr <- readRDS(arg(args, "caqtl_variants"))

# ---- caQTL status of each peak ----
ovP <- findOverlaps(peaks_gr, pgr)
pa <- data.table(peak_idx = queryHits(ovP), target = names(pgr)[subjectHits(ovP)], n_snps = pgr$n_caqtl_snps_hg19[subjectHits(ovP)], minp = pgr$min_pvalue_peak[subjectHits(ovP)])
pa <- pa[, .(n_caqtl_target_peaks = uniqueN(target), caqtl_target_peaks = paste(target, collapse = ";"), n_caqtl_snps_regulating = sum(n_snps), min_caqtl_pvalue = min(minp)), by = peak_idx]
ann[, `:=`(n_caqtl_target_peaks = 0L, caqtl_target_peaks = NA_character_, n_caqtl_snps_regulating = 0L, min_caqtl_pvalue = NA_real_)]
ann[pa$peak_idx, `:=`(n_caqtl_target_peaks = pa$n_caqtl_target_peaks, caqtl_target_peaks = pa$caqtl_target_peaks, n_caqtl_snps_regulating = pa$n_caqtl_snps_regulating, min_caqtl_pvalue = pa$min_caqtl_pvalue)]
ann[, caqtl_peak_positive := n_caqtl_target_peaks >= 1]
ovV <- findOverlaps(peaks_gr, vgr)
ann[, n_caqtl_variants_in_peak := as.integer(table(factor(queryHits(ovV), levels = seq_len(nrow(ann)))))]
n_pos <- sum(ann$caqtl_peak_positive); n_target_in <- length(unique(subjectHits(ovP)))
log_msg(sprintf("caQTL target peaks overlapping a peak: %d of %d; caQTL-positive peaks %d (%.2f%%)", n_target_in, length(pgr), n_pos, 100 * n_pos / nrow(ann)))
save_rds(as.data.frame(ann), arg(args, "out_peak_annotation"))

counts <- data.table(metric = c("caqtl_target_peaks_unique_hg38", "caqtl_target_peaks_overlapping_atlas_peak", "atlas_peaks_standard", "caqtl_positive_peaks", "pct_caqtl_positive_peaks",
                                "caqtl_variants_unique_hg38", "caqtl_variants_inside_atlas_peaks", "peaks_both_caqtl_positive_and_eqtl_positive"),
                     value = c(length(pgr), n_target_in, nrow(ann), n_pos, round(100 * n_pos / nrow(ann), 3), length(vgr), length(unique(subjectHits(ovV))), sum(ann$caqtl_peak_positive & ann$eqtl_positive)))
write_tsv(counts, arg(args, "out_counts"))

# ---- Adjusted logistic model; universe = peaks on chromosomes that carry caQTL target peaks ----
caqtl_chr <- unique(as.character(seqnames(pgr)))
uni <- ann[chr %in% caqtl_chr]
log_msg(sprintf("model universe: %d peaks on %d chromosomes with caQTL target peaks (%d excluded)", nrow(uni), length(caqtl_chr), nrow(ann) - nrow(uni)))
fit <- fit_peak_enrichment(uni, "caqtl_peak_positive", cluster_bin_bp = cluster_bin_bp, fdr = fdr)
enr <- add_egfr_columns(fit$enrichment, arg(args, "egfr_enrichment"))
eq <- read_tsv(arg(args, "eqtl_enrichment"))
enr$gtex_eqtl_adjusted_OR <- eq$adjusted_OR[match(enr$cell_type, eq$cell_type)]; enr$gtex_eqtl_FDR <- eq$FDR_BH[match(enr$cell_type, eq$cell_type)]
enr$universe <- sprintf("%d atlas peaks on %s", nrow(uni), paste(caqtl_chr, collapse = ","))
stopifnot(all(enr$converged), all(enr$all_coefficients_finite))
write_tsv(enr, arg(args, "out_enrichment"))
write_tsv(fit$coefficients, arg(args, "out_coefficients"))
log_msg("adjusted OR: ", paste(sprintf("%s %.2f [%.2f-%.2f] FDR %.2g", enr$cell_type, enr$adjusted_OR, enr$CI_low, enr$CI_high, enr$FDR_BH), collapse = "; "))

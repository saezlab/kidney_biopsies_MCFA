# (C) Sergio Alejandro Gómez
# Aim: Test whether cell-type-specific peaks are enriched for kidney-cortex eQTLs.
# Inputs: peak universe, CPM matrix, specificity matrix, eGFR enrichment table, GTEx kidney-cortex unique eVariants,
#   Ensembl 108 gene annotation.
# Outputs: enrichment table, model coefficients, peak counts, and the annotated peak table (covariates + eQTL status),
#   which the caQTL enrichment reuses.
# Usage: Rscript workflow/scripts/genetics/10_eqtl_enrichment.R --peak_universe=results/genetics/with_ptC5/objects/peak_universe.rds --cpm=results/genetics/with_ptC5/objects/accessibility_cpm.rds --specificity=results/genetics/with_ptC5/objects/specificity_matrix.rds --egfr_enrichment=results/genetics/with_ptC5/tables/egfr_enrichment.tsv --gtex_evariants=results/genetics/qtl_resources/gtex_kidney_unique_evariants_hg38.rds --gene_annotation=results/genetics/qtl_resources/gene_annotation_hg38_ensembl108.rds --out_enrichment=results/genetics/with_ptC5/tables/eqtl_enrichment.tsv --out_coefficients=results/genetics/with_ptC5/tables/eqtl_enrichment_coefficients.tsv --out_counts=results/genetics/with_ptC5/tables/eqtl_peak_counts.tsv --out_peak_annotation=results/genetics/with_ptC5/objects/peak_annotation_eqtl.rds

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "qtl.R"))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))

args <- parse_args()
cluster_bin_bp <- as.numeric(arg(args, "cluster_bin_bp", "1e6"))
fdr <- as.numeric(arg(args, "fdr", "0.05"))

peaks <- readRDS(arg(args, "peak_universe")); pk <- as.data.table(peaks[peaks$is_standard_chr, ])
peaks_gr <- GRanges(pk$chr, IRanges(pk$start, pk$end)); names(peaks_gr) <- pk$peak
spec <- readRDS(arg(args, "specificity"))
stopifnot(identical(rownames(spec), pk$peak), setequal(colnames(spec), POPULATIONS))
cpm <- readRDS(arg(args, "cpm"))[pk$peak, POPULATIONS]
ev <- readRDS(arg(args, "gtex_evariants"))

# ---- eQTL status of each peak ----
ov <- findOverlaps(ev, peaks_gr)
evd <- as.data.table(mcols(ev))
hit <- data.table(peak_idx = subjectHits(ov), variant_id = names(ev)[queryHits(ov)], egenes = evd$egenes[queryHits(ov)])
pa <- hit[, .(n_unique_evariants = uniqueN(variant_id), n_unique_egenes = uniqueN(unlist(strsplit(egenes, ";"))),
              egenes = paste(sort(unique(unlist(strsplit(egenes, ";")))), collapse = ";")), by = peak_idx]
ann <- data.table(peak = pk$peak, chr = pk$chr, start = pk$start, end = pk$end, peak_width = pk$width, n_unique_evariants = 0L, n_unique_egenes = 0L, egenes = NA_character_)
ann[pa$peak_idx, `:=`(n_unique_evariants = pa$n_unique_evariants, n_unique_egenes = pa$n_unique_egenes, egenes = pa$egenes)]
ann[, eqtl_positive := n_unique_evariants >= 1]
n_in <- length(unique(queryHits(ov))); n_pos <- sum(ann$eqtl_positive)
log_msg(sprintf("eVariants inside peaks: %d of %d; eQTL-positive peaks %d (%.2f%%)", n_in, length(ev), n_pos, 100 * n_pos / nrow(ann)))

# ---- Peak covariates: mean accessibility, GC content, distance to the nearest protein-coding TSS ----
ann[, mean_accessibility := rowMeans(cpm)]
gc <- numeric(nrow(pk))
for (ii in split(seq_len(nrow(pk)), ceiling(seq_len(nrow(pk)) / 40000))) {
  s <- getSeq(BSgenome.Hsapiens.UCSC.hg38, peaks_gr[ii]); gc[ii] <- as.numeric(Biostrings::letterFrequency(s, "GC", as.prob = TRUE)); rm(s)
}
ann[, gc_content := gc]
tss_pc <- protein_coding_tss(readRDS(arg(args, "gene_annotation")))
d2t <- distanceToNearest(peaks_gr, tss_pc, ignore.strand = TRUE)
dist <- rep(NA_integer_, nrow(pk)); dist[queryHits(d2t)] <- mcols(d2t)$distance
ann[, distance_to_tss := dist]
stopifnot(!anyNA(ann$distance_to_tss), !anyNA(ann$gc_content))
for (p in POPULATIONS) ann[, (paste0("specific_", p)) := spec[, p]]
save_rds(as.data.frame(ann), arg(args, "out_peak_annotation"))

counts <- data.table(metric = c("unique_evariants", "evariants_in_atlas_peaks", "atlas_peaks_standard", "eqtl_positive_peaks", "pct_eqtl_positive_peaks", "unique_egenes_in_peaks"),
                     value = c(length(ev), n_in, nrow(ann), n_pos, round(100 * n_pos / nrow(ann), 3), length(unique(unlist(strsplit(ann$egenes[ann$eqtl_positive], ";"))))))
write_tsv(counts, arg(args, "out_counts"))

# ---- Adjusted logistic model; universe = peaks on chromosomes that carry eVariants ----
ev_chr <- unique(as.character(seqnames(ev)))
uni <- ann[chr %in% ev_chr]
log_msg(sprintf("model universe: %d peaks on %d chromosomes with eVariants (%d excluded)", nrow(uni), length(ev_chr), nrow(ann) - nrow(uni)))
fit <- fit_peak_enrichment(uni, "eqtl_positive", cluster_bin_bp = cluster_bin_bp, fdr = fdr)
enr <- add_egfr_columns(fit$enrichment, arg(args, "egfr_enrichment"))
enr$universe <- sprintf("%d atlas peaks on %s", nrow(uni), paste(ev_chr, collapse = ","))
stopifnot(all(enr$converged), all(enr$all_coefficients_finite))
write_tsv(enr, arg(args, "out_enrichment"))
write_tsv(fit$coefficients, arg(args, "out_coefficients"))
log_msg("adjusted OR: ", paste(sprintf("%s %.2f [%.2f-%.2f] FDR %.2g", enr$cell_type, enr$adjusted_OR, enr$CI_low, enr$CI_high, enr$FDR_BH), collapse = "; "))

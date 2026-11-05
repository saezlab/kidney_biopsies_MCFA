# (C) Sergio Alejandro Gómez
# Aim: Test which kidney populations are enriched for CKD genetic risk, with the same gchromVAR machinery as for eGFR.
# Inputs: peak universe and CPM + 1 matrix (03), fine-mapped CKD variants (hg38 GRanges) and the CKD locus QC table.
# Outputs: CKD enrichment table and per-seed summary, with and without the UMOD locus.
# Usage: Rscript workflow/scripts/genetics/18_ckd_enrichment.R --peak_universe=results/genetics/with_ptC5/objects/peak_universe.rds --cpm_plus1=results/genetics/with_ptC5/objects/accessibility_cpm_plus1.rds --ckd_variants=results/genetics/ckd_finemapping/ckd_finemapped_variants_hg38.rds --ckd_locus_qc=results/genetics/ckd_finemapping/ckd_locus_qc.tsv --out_enrichment=results/genetics/with_ptC5/tables/ckd_enrichment.tsv --out_by_seed=results/genetics/with_ptC5/tables/ckd_enrichment_by_seed.tsv --out_without_umod=results/genetics/with_ptC5/tables/ckd_enrichment_without_umod.tsv --out_without_umod_by_seed=results/genetics/with_ptC5/tables/ckd_enrichment_without_umod_by_seed.tsv

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "accessibility.R"))
source(file.path(script_dir, "R", "gchromvar.R"))

args <- parse_args()
seeds <- as.integer(strsplit(arg(args, "seeds", "2024,1,2,3,4,5"), ",")[[1]])
umod_locus <- arg(args, "umod_locus", "CKD_L22_chr16_19350163_21421556")

peaks_gr <- peak_granges(readRDS(arg(args, "peak_universe")))
cpm1 <- readRDS(arg(args, "cpm_plus1"))[names(peaks_gr), , drop = FALSE]   # stored (alphabetical) column order
qc <- read_tsv(arg(args, "ckd_locus_qc"))
variants <- readRDS(arg(args, "ckd_variants"))
eligible <- qc$locus[qc$final_class %in% c("PASS", "PASS_WITH_WARNING")]
stopifnot(setequal(unique(variants$locus), eligible), umod_locus %in% eligible)

ckd_enrichment <- function(vgr, dropped = character(0)) {
  g <- run_gchromvar(cpm1, peaks_gr, vgr, "Kidney_CKD", seeds)
  res <- g$result
  res$population <- res$cell_type
  res$display_label <- unname(DISPLAY_LABEL[res$cell_type])
  res$direction <- ifelse(res$z > 0, "enriched", "depleted")
  res$n_variants <- length(vgr)
  cls <- qc$final_class[!(qc$locus %in% dropped)]
  res$included_loci <- length(setdiff(eligible, dropped)); res$pass_loci <- sum(cls == "PASS"); res$warning_loci <- sum(cls == "PASS_WITH_WARNING")
  res$blocked_loci <- sum(qc$final_class == "BLOCKED_BY_QC"); res$failed_loci <- sum(qc$final_class == "FAILED")
  res <- res[order(match(res$cell_type, DISPLAY_ORDER)), c("population", "display_label", "z", "deviation", "p", "FDR", "rank", "direction", "significance",
             "seed", "seed_mean_z", "seed_sd_z", "n_peaks_with_weight", "weight_sum", "n_variants", "included_loci", "pass_loci", "warning_loci",
             "blocked_loci", "failed_loci")]
  rownames(res) <- NULL
  list(result = res, by_seed = g$by_seed)
}

all_loci <- ckd_enrichment(variants)
write_tsv(all_loci$result, arg(args, "out_enrichment"))
write_tsv(all_loci$by_seed, arg(args, "out_by_seed"))

no_umod <- ckd_enrichment(variants[variants$locus != umod_locus], dropped = umod_locus)
res <- no_umod$result
res$dropped_locus <- umod_locus
res$dropped_locus_class <- qc$final_class[qc$locus == umod_locus]
res$n_variants_removed <- sum(variants$locus == umod_locus)
res$sum_PIP_removed <- sum(variants$pip[variants$locus == umod_locus])
write_tsv(res, arg(args, "out_without_umod"))
write_tsv(no_umod$by_seed, arg(args, "out_without_umod_by_seed"))
log_msg("CKD: ", paste(sprintf("%s %.2f", all_loci$result$population, all_loci$result$z), collapse = " "))
log_msg("CKD without UMOD: ", paste(sprintf("%s %.2f", res$population, res$z), collapse = " "))

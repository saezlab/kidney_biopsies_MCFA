# (C) Sergio Alejandro Gómez
# Aim: Test which kidney populations are enriched for eGFR genetic risk.
# Inputs: peak universe and CPM + 1 matrix (03_accessibility.R), fine-mapped eGFR variants (1-bp, hg38).
# Outputs: enrichment table (seed 2024, with mean and SD of z over the six seeds), per-seed summary.
# Usage: Rscript workflow/scripts/genetics/04_egfr_enrichment.R --peak_universe=results/genetics/with_ptC5/objects/peak_universe.rds --cpm_plus1=results/genetics/with_ptC5/objects/accessibility_cpm_plus1.rds --variants=results/genetics/egfr_finemapping/egfr_finemapped_variants_hg38.tsv --include_ptC5=true --out_enrichment=results/genetics/with_ptC5/tables/egfr_enrichment.tsv --out_by_seed=results/genetics/with_ptC5/tables/egfr_enrichment_by_seed.tsv

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "accessibility.R"))
source(file.path(script_dir, "R", "gchromvar.R"))

args <- parse_args()
version <- if (arg_true(args, "include_ptC5", NULL)) "with_ptC5" else "without_ptC5"
seeds <- as.integer(strsplit(arg(args, "seeds", "2024,1,2,3,4,5"), ",")[[1]])

peaks_gr <- peak_granges(readRDS(arg(args, "peak_universe")))
cpm1 <- readRDS(arg(args, "cpm_plus1"))[names(peaks_gr), POPULATIONS]
variants <- read_egfr_variants(arg(args, "variants"))
log_msg(sprintf("%d fine-mapped variants, %d with PIP > 0.1", length(variants), sum(variants$pip > 0.1)))

g <- run_gchromvar(cpm1, peaks_gr, variants, "Kidney_eGFR_CARMA", seeds)
res <- g$result
res$variant_set <- "canonical_1bp"
res$ptc5_version <- version
res$peak_universe <- PEAK_UNIVERSE_LABEL
res <- res[, c("cell_type", "z", "deviation", "p", "FDR", "rank", "significance", "variant_set", "seed", "n_peaks_with_weight",
               "weight_sum", "seed_mean_z", "seed_sd_z", "ptc5_version", "peak_universe")]
write_tsv(res, arg(args, "out_enrichment"))
write_tsv(g$by_seed, arg(args, "out_by_seed"))
log_msg("enriched (z > 0, FDR < 0.05): ", paste(res$cell_type[res$significance == "enriched_FDR<0.05"], collapse = ", "))

# (C) Sergio Alejandro Gómez
# Aim: Prioritise candidate eGFR variants: fine-mapped variants with PIP > 0.1 that lie in an ArchR peak specific to at
#   least one population enriched for eGFR.
# Outputs: candidate table, funnel table, candidate GRanges (.rds) and a BED of the candidates (input of the motif analysis).
# Usage: Rscript workflow/scripts/genetics/06_candidate_variants.R --variants=results/genetics/egfr_finemapping/egfr_finemapped_variants_hg38.tsv --peak_universe=results/genetics/with_ptC5/objects/peak_universe.rds --cpm=results/genetics/with_ptC5/objects/accessibility_cpm.rds --specific_peaks=results/genetics/with_ptC5/objects/specific_peaks.rds --specificity=results/genetics/with_ptC5/objects/specificity_matrix.rds --include_ptC5=true --out_candidates=results/genetics/with_ptC5/tables/candidate_variants.tsv --out_funnel=results/genetics/with_ptC5/tables/candidate_funnel.tsv --out_candidates_rds=results/genetics/with_ptC5/objects/candidate_variants.rds --out_bed=results/genetics/with_ptC5/objects/candidate_variants.bed

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "accessibility.R"))
source(file.path(script_dir, "R", "gchromvar.R"))

args <- parse_args()
version <- if (arg_true(args, "include_ptC5", NULL)) "with_ptC5" else "without_ptC5"
pip_threshold <- as.numeric(arg(args, "pip_threshold", "0.1"))

variants <- read_egfr_variants(arg(args, "variants"))
peaks_gr <- peak_granges(readRDS(arg(args, "peak_universe")))
cpm <- readRDS(arg(args, "cpm"))[names(peaks_gr), ]   # stored (alphabetical) column order
sp <- readRDS(arg(args, "specific_peaks"))
spec_all <- readRDS(arg(args, "specificity"))
enriched <- sp$enriched_populations
stopifnot(identical(sp$ptc5_version, version), identical(rownames(spec_all), names(peaks_gr)))
universe_text <- sprintf("%s ArchR peaks (%s)", format(length(peaks_gr), big.mark = ","), version)

# ---- funnel ----
f1 <- variants[variants$pip > pip_threshold]
ov_any <- findOverlaps(f1, peaks_gr)
f2 <- f1[unique(queryHits(ov_any))]
ov_spec <- findOverlaps(f2, peaks_gr[Reduce(union, sp$peak_sets)])
f3 <- f2[unique(queryHits(ov_spec))]
n_loci <- function(g) length(unique(g$locus))
funnel <- data.frame(
  stage = c("0_all_CARMA_variants_1bp_hg38", "1_PIP_gt_0.1", "2_in_any_ArchR_peak", "3_in_specific_peak_of_enriched_population"),
  rule = c("canonical 1-bp CARMA representation", sprintf("PIP > %g (strict)", pip_threshold),
           sprintf("overlap with >= 1 of the %s", universe_text),
           sprintf("overlap with >= 1 peak specific to a positively enriched population {%s}", paste(enriched, collapse = ","))),
  n_variants = c(length(variants), length(f1), length(f2), length(f3)),
  n_loci = c(n_loci(variants), n_loci(f1), n_loci(f2), n_loci(f3)),
  n_peaks_involved = c(NA, NA, length(unique(subjectHits(ov_any))), length(unique(subjectHits(ov_spec)))),
  pip_median = c(median(variants$pip), median(f1$pip), median(f2$pip), median(f3$pip)),
  n_pip_gt_0.5 = c(sum(variants$pip > .5), sum(f1$pip > .5), sum(f2$pip > .5), sum(f3$pip > .5)), stringsAsFactors = FALSE)
funnel$ptc5_version <- version
write_tsv(funnel, arg(args, "out_funnel"))
log_msg("funnel: ", paste(funnel$n_variants, collapse = " -> "))

# ---- candidate table: peak, populations to which the peak is specific, accessibility of the peak ----
ov_c <- findOverlaps(f3, peaks_gr)
peak_of <- tapply(names(peaks_gr)[subjectHits(ov_c)], queryHits(ov_c), function(x) paste(x, collapse = ";"))
peak_of <- unname(peak_of[as.character(seq_along(f3))])
first_peak <- vapply(strsplit(peak_of, ";"), `[`, "", 1)
spec_enriched <- vapply(first_peak, function(p) paste(enriched[spec_all[p, enriched]], collapse = ";"), "")
spec_any <- vapply(first_peak, function(p) paste(colnames(spec_all)[spec_all[p, ]], collapse = ";"), "")
cand <- data.frame(rsID = f3$rsid, chr = as.character(seqnames(f3)), hg38_position = start(f3), PIP = f3$pip, locus = f3$locus,
                   credible_set = f3$credible_set, REF = f3$REF, ALT = f3$ALT, ref_check = f3$ref_check, effect_allele_GWAS = f3$allele1_effect,
                   gwas_z = f3$gwas_z, maf = f3$maf, overlapping_peak = peak_of,
                   specific_cell_types = unname(spec_enriched), n_specific_cell_types = lengths(strsplit(unname(spec_enriched), ";")),
                   specific_cell_types_any_population = unname(spec_any),
                   gchromVAR_enriched_cell_types = paste(enriched, collapse = ";"), stringsAsFactors = FALSE)
cpm_cols <- cpm[first_peak, , drop = FALSE]; colnames(cpm_cols) <- paste0("peak_cpm_", colnames(cpm_cols))
cand <- cbind(cand, round(cpm_cols, 3))
cand <- cand[order(cand$chr, cand$hg38_position), ]; rownames(cand) <- NULL
stopifnot(!anyDuplicated(cand$rsID), all(cand$n_specific_cell_types >= 1), all(grepl("matches_genome", cand$ref_check)))
cand$ptc5_version <- version
cand$peak_universe <- PEAK_UNIVERSE_LABEL
write_tsv(cand, arg(args, "out_candidates"))

f3$overlapping_peak <- cand$overlapping_peak[match(f3$rsid, cand$rsID)]
f3$specific_cell_types <- cand$specific_cell_types[match(f3$rsid, cand$rsID)]
save_rds(f3, arg(args, "out_candidates_rds"))

# BED (0-based start) with the hg38-verified alleles in the name, as the motif analysis reads it
bed <- data.frame(chr = cand$chr, start = cand$hg38_position - 1L, end = cand$hg38_position,
                  name = paste(cand$chr, cand$hg38_position, cand$REF, cand$ALT, sep = ":"), score = cand$PIP, strand = "*")
dir.create(dirname(arg(args, "out_bed")), recursive = TRUE, showWarnings = FALSE)
write.table(bed, arg(args, "out_bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
log_msg(nrow(cand), " candidate variants (", version, ")")

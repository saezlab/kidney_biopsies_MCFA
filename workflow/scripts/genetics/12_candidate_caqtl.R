# (C) Sergio Alejandro Gómez
# Aim: Describe the kidney caQTL evidence of each eGFR candidate variant.
# Inputs: candidate variants, candidate-gene pairs from 09_target_genes.R, motif disruption per variant,
#   bASA/RASQUAL unique caQTL variants, target peaks and associations (hg38).
# Outputs: one row per candidate (summary), caQTL matches, evidence layers per candidate-gene pair, and the pairs
#   with an exact caQTL variant plus a gene-linking source.
# Usage: Rscript workflow/scripts/genetics/12_candidate_caqtl.R --candidates=results/genetics/with_ptC5/tables/candidate_variants.tsv --target_genes=results/genetics/with_ptC5/tables/candidate_target_gene_long.tsv --motif_disruption=results/genetics/with_ptC5/tables/motif_disruption_by_variant.tsv --caqtl_variants=results/genetics/qtl_resources/basa_unique_caqtl_variants_hg38.rds --caqtl_peaks=results/genetics/qtl_resources/basa_caqtl_peaks_hg38.rds --caqtl_associations=results/genetics/qtl_resources/basa_caqtl_associations_hg38.tsv.gz --out_summary=results/genetics/with_ptC5/tables/candidate_caqtl_summary.tsv --out_exact_matches=results/genetics/with_ptC5/tables/candidate_caqtl_exact_matches.tsv --out_multilayer=results/genetics/with_ptC5/tables/multilayer_candidate_gene_evidence.tsv --out_multilayer_high=results/genetics/with_ptC5/tables/multilayer_high_support_examples.tsv

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "qtl.R"))

args <- parse_args()

cand <- as.data.table(read_tsv(arg(args, "candidates")))
stopifnot(nrow(cand) > 0)
pk <- parse_archr_peaks(cand$overlapping_peak); cand[, `:=`(atlas_peak_chr = pk$chr, atlas_peak_start = pk$start, atlas_peak_end = pk$end)]
cgr <- GRanges(cand$chr, IRanges(cand$hg38_position, cand$hg38_position)); names(cgr) <- cand$rsID
agr <- GRanges(cand$atlas_peak_chr, IRanges(cand$atlas_peak_start, cand$atlas_peak_end)); names(agr) <- cand$rsID
vgr <- readRDS(arg(args, "caqtl_variants"))
pgr <- readRDS(arg(args, "caqtl_peaks"))
assoc <- fread(arg(args, "caqtl_associations"), showProgress = FALSE)

# ---- A. exact caQTL variant ----
vm <- as.data.table(mcols(vgr)); vm[, `:=`(chr_hg38 = as.character(seqnames(vgr)), pos_hg38 = start(vgr))]
pos_hit <- match(paste(cand$chr, cand$hg38_position), paste(vm$chr_hg38, vm$pos_hg38))
cand[, caqtl_position_match := !is.na(pos_hit)]
cand[, caqtl_rsid := vm$rs_ID[pos_hit]]; cand[, caqtl_REF_hg38 := vm$REF_hg38[pos_hit]]; cand[, caqtl_ALT_hg38 := vm$ALT_hg38[pos_hit]]
cand[, exact_variant_caqtl := caqtl_position_match & REF == caqtl_REF_hg38 & ALT == caqtl_ALT_hg38]
cand[, caqtl_allele_relation := fifelse(!caqtl_position_match, "NO_POSITION_MATCH", fifelse(exact_variant_caqtl, "EXACT_REF_ALT",
                                 fifelse(REF == caqtl_ALT_hg38 & ALT == caqtl_REF_hg38, "SAME_ALLELES_SWAPPED_ORIENTATION", "ALLELE_MISMATCH")))]
log_msg(sprintf("A. exact caQTL variants: %d of %d", sum(cand$exact_variant_caqtl), nrow(cand)))

ex <- merge(cand[exact_variant_caqtl == TRUE, .(rsID, PIP, locus, specific_cell_types, chr, hg38_position, REF, ALT, overlapping_peak, atlas_peak_chr, atlas_peak_start, atlas_peak_end)],
            assoc[, .(rsID = rs_ID, Feature_ID, peak_hg38, peak_chr_hg38, peak_start_hg38, peak_end_hg38, peak_usable_hg38, usable_hg38, snp_inside_target_peak_hg38,
                      Distance_hg19 = Distance, Effect_size_pi = Effect_size, ChiSquare_statistic, Pvalue, Log10_BH_Qvalue = Log10_Benjamini_Hochberg_Qvalue, Validated, Allele_frequency, MAF)],
            by = "rsID", all.x = TRUE)
ex[, own_atlas_peak_overlaps_this_target_peak := peak_usable_hg38 & atlas_peak_chr == peak_chr_hg38 & atlas_peak_start <= peak_end_hg38 & atlas_peak_end >= peak_start_hg38]
ex[is.na(own_atlas_peak_overlaps_this_target_peak), own_atlas_peak_overlaps_this_target_peak := FALSE]
ex[, candidate_inside_this_target_peak := peak_usable_hg38 & hg38_position >= peak_start_hg38 & hg38_position <= peak_end_hg38]
ex[is.na(candidate_inside_this_target_peak), candidate_inside_this_target_peak := FALSE]
ex[, effect_direction_note := "RASQUAL effect size pi is relative to the resource reference allele; orientation vs the eGFR effect allele is not asserted"]
setorder(ex, -PIP, rsID, Pvalue)
write_tsv(ex, arg(args, "out_exact_matches"))

# ---- B. inside a caQTL target peak; C. own peak is the target of its own caQTL ----
ovB <- findOverlaps(cgr, pgr)
pm <- data.table(rsID = names(cgr)[queryHits(ovB)], Feature_ID = names(pgr)[subjectHits(ovB)], peak_hg38 = pgr$peak_hg38[subjectHits(ovB)])
cand[, candidate_inside_caqtl_target_peak := rsID %in% pm$rsID]
cand[, n_caqtl_target_peaks_containing_candidate := as.integer(table(factor(pm$rsID, levels = rsID)))]
cand[, caqtl_target_peaks_containing_candidate := vapply(rsID, function(r) paste(pm$peak_hg38[pm$rsID == r], collapse = ";"), "")]
own <- ex[own_atlas_peak_overlaps_this_target_peak == TRUE, .(n_own_target_peaks = .N, own_target_peaks = paste(peak_hg38, collapse = ";"), own_target_min_pvalue = min(Pvalue),
                                                              own_target_effect_size_pi = paste(sprintf("%.3f", Effect_size_pi), collapse = ";")), by = rsID]
cand[, own_peak_is_caqtl_target := rsID %in% own$rsID]
cand <- merge(cand, own, by = "rsID", all.x = TRUE, sort = FALSE); cand[is.na(n_own_target_peaks), n_own_target_peaks := 0L]
cand[, direct_kidney_caqtl_support := exact_variant_caqtl & own_peak_is_caqtl_target]
ovC <- findOverlaps(agr, pgr)
cand[, atlas_peak_overlaps_any_caqtl_target_peak := rsID %in% names(agr)[queryHits(ovC)]]
ovV <- findOverlaps(agr, vgr)
cand[, n_caqtl_variants_in_atlas_peak := as.integer(table(factor(names(agr)[queryHits(ovV)], levels = rsID)))]
log_msg(sprintf("B. inside a target peak: %d; C. own peak is the caQTL target: %d; direct support: %d",
                sum(cand$candidate_inside_caqtl_target_peak), sum(cand$own_peak_is_caqtl_target), sum(cand$direct_kidney_caqtl_support)))
exs <- ex[, .(n_caqtl_target_peaks = uniqueN(Feature_ID), caqtl_target_peaks_hg38 = paste(peak_hg38[!is.na(peak_hg38)], collapse = ";"), caqtl_min_pvalue = min(Pvalue),
              caqtl_effect_size_pi = paste(sprintf("%.3f", Effect_size_pi), collapse = ";"), caqtl_validated = paste(unique(Validated), collapse = ";")), by = rsID]
cand <- merge(cand, exs, by = "rsID", all.x = TRUE, sort = FALSE); cand[is.na(n_caqtl_target_peaks), n_caqtl_target_peaks := 0L]
cand[, caqtl_class := fifelse(direct_kidney_caqtl_support, "DIRECT_KIDNEY_CAQTL_SUPPORT", fifelse(exact_variant_caqtl, "EXACT_CAQTL_VARIANT_target_peak_elsewhere",
                       fifelse(candidate_inside_caqtl_target_peak, "INSIDE_CAQTL_TARGET_PEAK_not_the_caqtl_variant", "NO_CAQTL_OVERLAP")))]

# Target-gene fields from 09_target_genes.R
long <- as.data.table(read_tsv(arg(args, "target_genes")))
pl <- function(x) paste(sort(unique(x)), collapse = ";")
tg <- long[, .(open4gene_direct_genes = pl(gene[open4gene_variant_in_linked_peak]), n_open4gene_direct_genes = sum(open4gene_variant_in_linked_peak),
               gtex_egenes = pl(gene[gtex_exact_kidney_eqtl]), n_gtex_egenes = sum(gtex_exact_kidney_eqtl),
               convergent_genes = pl(gene[n_functional_sources >= 2]), n_convergent_genes = sum(n_functional_sources >= 2),
               functional_target_genes = pl(gene[n_functional_sources >= 1]), n_functional_target_genes = sum(n_functional_sources >= 1),
               nearest_protein_coding_gene = pl(gene[nearest_gene])), by = rsID]
cand <- merge(cand, tg, by = "rsID", all.x = TRUE, sort = FALSE)
summ <- cand[, .(rsID, PIP, locus, credible_set, chr, hg38_position, REF, ALT, atlas_peak = overlapping_peak, candidate_cell_types = specific_cell_types, n_candidate_cell_types = n_specific_cell_types,
                 exact_variant_caqtl, caqtl_allele_relation, caqtl_rsid, n_caqtl_target_peaks, caqtl_target_peaks_hg38, caqtl_min_pvalue, caqtl_effect_size_pi, caqtl_validated,
                 candidate_inside_caqtl_target_peak, n_caqtl_target_peaks_containing_candidate, caqtl_target_peaks_containing_candidate,
                 own_peak_is_caqtl_target, n_own_target_peaks, own_target_peaks, own_target_min_pvalue, own_target_effect_size_pi, direct_kidney_caqtl_support, caqtl_class,
                 atlas_peak_overlaps_any_caqtl_target_peak, n_caqtl_variants_in_atlas_peak,
                 open4gene_direct_genes, n_open4gene_direct_genes, gtex_egenes, n_gtex_egenes, convergent_genes, n_convergent_genes, functional_target_genes, n_functional_target_genes, nearest_protein_coding_gene)]
setorder(summ, -direct_kidney_caqtl_support, -exact_variant_caqtl, -candidate_inside_caqtl_target_peak, -PIP, rsID)
write_tsv(summ, arg(args, "out_summary"))

# ---- Evidence layers per candidate-gene pair ----
mb <- as.data.table(read_tsv(arg(args, "motif_disruption")))
ml <- merge(long, summ[, .(rsID, exact_variant_caqtl, caqtl_rsid, caqtl_target_peaks_hg38, n_caqtl_target_peaks, caqtl_min_pvalue, caqtl_effect_size_pi, caqtl_validated,
                           candidate_inside_caqtl_target_peak, caqtl_target_peaks_containing_candidate, own_peak_is_caqtl_target, own_target_peaks, direct_kidney_caqtl_support, caqtl_class,
                           atlas_peak_overlaps_any_caqtl_target_peak, n_caqtl_variants_in_atlas_peak)], by = "rsID", all.x = TRUE, sort = FALSE)
ml <- merge(ml, mb[, .(rsID = rsid, motifbreakr_settingA_n_strong_motifs = n_strong_motifs, motifbreakr_settingA_n_distinct_tfs = n_distinct_tfs)], by = "rsID", all.x = TRUE, sort = FALSE)
ml[, L_egfr_carma_pip_gt_0.1 := PIP > 0.1]; ml[, L_cell_specific_open_chromatin := TRUE]
ml[, L_open4gene_peak_gene_link := open4gene_variant_in_linked_peak %in% TRUE]; ml[, L_gtex_exact_kidney_eqtl := gtex_exact_kidney_eqtl %in% TRUE]
ml[, L_kidney_caqtl_exact_variant := exact_variant_caqtl %in% TRUE]
ml[, L_kidney_caqtl_target_peak_context := (own_peak_is_caqtl_target %in% TRUE) | (candidate_inside_caqtl_target_peak %in% TRUE)]
ml[, L_nearest_gene_proximity := nearest_gene %in% TRUE]
ml[, L_motifbreakr_settingA_strong_disruption := !is.na(motifbreakr_settingA_n_strong_motifs) & motifbreakr_settingA_n_strong_motifs >= 1]
ml[, n_counted_evidence_layers := as.integer(L_open4gene_peak_gene_link) + as.integer(L_gtex_exact_kidney_eqtl) + as.integer(L_kidney_caqtl_exact_variant) + as.integer(L_kidney_caqtl_target_peak_context)]
stopifnot(all(as.integer(ml$L_open4gene_peak_gene_link) + as.integer(ml$L_gtex_exact_kidney_eqtl) == ml$n_functional_sources))
ml[, MULTILAYER_REGULATORY_SUPPORT := L_kidney_caqtl_exact_variant & (L_open4gene_peak_gene_link | L_gtex_exact_kidney_eqtl)]
ml[, evidence_layers_present := apply(as.matrix(.SD), 1, function(r) paste(sub("^L_", "", names(r))[r], collapse = ";")), .SDcols = grep("^L_", names(ml), value = TRUE)]
ml[, target_gene_confidence := pair_category]
setorder(ml, -n_counted_evidence_layers, -PIP, rsID, gene)
keep <- c("rsID", "gene", "target_gene_confidence", "pair_category", "n_functional_sources", "functional_sources", "n_counted_evidence_layers", "MULTILAYER_REGULATORY_SUPPORT", "evidence_layers_present",
          grep("^L_", names(ml), value = TRUE), "PIP", "locus", "credible_set", "chr", "hg38_position", "REF", "ALT", "atlas_peak", "candidate_cell_types", "n_candidate_cell_types",
          "open4gene_source_celltypes", "open4gene_best_min_fdr", "celltype_concordant_open4gene", "gtex_variant_id", "gtex_slope", "gtex_gene_qval", "expression_direction_of_eGFR_raising_allele",
          "caqtl_rsid", "caqtl_target_peaks_hg38", "n_caqtl_target_peaks", "caqtl_min_pvalue", "caqtl_effect_size_pi", "caqtl_validated", "caqtl_class", "own_target_peaks",
          "caqtl_target_peaks_containing_candidate", "atlas_peak_overlaps_any_caqtl_target_peak", "n_caqtl_variants_in_atlas_peak",
          "nearest_gene", "distance_to_nearest_protein_coding_tss", "promoter_proximal_le_2kb", "gene_biotype_ensembl108", "motifbreakr_settingA_n_strong_motifs", "motifbreakr_settingA_n_distinct_tfs")
ml_out <- ml[, ..keep]
write_tsv(ml_out, arg(args, "out_multilayer"))
write_tsv(ml_out[MULTILAYER_REGULATORY_SUPPORT == TRUE], arg(args, "out_multilayer_high"))
log_msg(sprintf("evidence layers: %d pairs; %d pairs / %d candidates with an exact caQTL variant plus a gene-linking source",
                nrow(ml_out), sum(ml_out$MULTILAYER_REGULATORY_SUPPORT), uniqueN(ml_out$rsID[ml_out$MULTILAYER_REGULATORY_SUPPORT])))

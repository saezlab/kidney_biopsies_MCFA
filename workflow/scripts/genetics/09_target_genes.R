# (C) Sergio Alejandro Gómez
# Aim: Assign putative target genes to the eGFR candidate variants.
# Inputs: candidate variants (table and GRanges), Open4Gene links lifted to hg38, GTEx kidney-cortex unique eVariants and
#   significant pairs, Ensembl 108 gene annotation.
# Outputs: candidate-gene pairs (long), one row per candidate (summary), evidence counts, convergent pairs, GTEx matches,
#   Open4Gene matches, and Open4Gene links that overlap the candidate's peak but not the variant (context only).
# Usage: Rscript workflow/scripts/genetics/09_target_genes.R --candidates=results/genetics/with_ptC5/tables/candidate_variants.tsv --candidates_rds=results/genetics/with_ptC5/objects/candidate_variants.rds --open4gene=results/genetics/qtl_resources/open4gene_hg38_unique.rds --gtex_evariants=results/genetics/qtl_resources/gtex_kidney_unique_evariants_hg38.rds --gtex_pairs=results/genetics/qtl_resources/gtex_kidney_significant_eqtl_pairs.rds --gene_annotation=results/genetics/qtl_resources/gene_annotation_hg38_ensembl108.rds --out_long=results/genetics/with_ptC5/tables/candidate_target_gene_long.tsv --out_summary=results/genetics/with_ptC5/tables/candidate_target_gene_summary.tsv --out_counts=results/genetics/with_ptC5/tables/candidate_target_gene_evidence_counts.tsv --out_convergent=results/genetics/with_ptC5/tables/convergent_candidate_gene_links.tsv --out_gtex_matches=results/genetics/with_ptC5/tables/candidate_gtex_exact_eqtl_matches.tsv --out_open4gene_matches=results/genetics/with_ptC5/tables/candidate_open4gene_matches.tsv --out_open4gene_peak_context=results/genetics/with_ptC5/tables/candidate_open4gene_atlas_peak_overlap.tsv

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "qtl.R"))

args <- parse_args()
promoter_window_bp <- as.integer(arg(args, "promoter_window_bp", "2000"))

cand <- as.data.table(read_tsv(arg(args, "candidates")))
cgr_in <- readRDS(arg(args, "candidates_rds"))
stopifnot(nrow(cand) > 0, all(cand$rsID %in% cgr_in$rsid))
n_cand <- nrow(cand)
log_msg(n_cand, " candidate variants")

# Direction of the ALT allele on eGFR (+1: ALT associated with higher eGFR), used to orient the GTEx effects.
m <- match(cand$rsID, cgr_in$rsid)
cand[, `:=`(effect_allele_hg38 = cgr_in$effect_allele_hg38_oriented[m], effect_allele_is_ALT = cgr_in$effect_allele_is_ALT[m])]
cand[, egfr_direction_of_ALT := fifelse(effect_allele_is_ALT, sign(gwas_z), -sign(gwas_z))]
cgr <- GRanges(cand$chr, IRanges(cand$hg38_position, cand$hg38_position)); names(cgr) <- cand$rsID
cand_ct <- strsplit(cand$specific_cell_types, ";")

# ---- Open4Gene: the variant lies inside a linked peak ----
o4 <- as.data.table(readRDS(arg(args, "open4gene")))
o4gr <- GRanges(o4$chr_hg38, IRanges(o4$start_hg38, o4$end_hg38))
ov_direct <- findOverlaps(cgr, o4gr)

open4gene_rows <- function(ov) {
  i <- queryHits(ov); j <- subjectHits(ov)
  d <- data.table(rsID = cand$rsID[i], PIP = cand$PIP[i], locus = cand$locus[i], chr = cand$chr[i], hg38_position = cand$hg38_position[i],
                  atlas_peak = cand$overlapping_peak[i], candidate_cell_types = cand$specific_cell_types[i],
                  open4gene_peak_hg38 = o4$peak_hg38[j], open4gene_peak_hg19 = o4$peak_hg19[j], open4gene_row_id = o4$row_id[j],
                  open4gene_source_celltype = o4$Celltype[j], crosswalk_tier = o4$crosswalk_tier[j],
                  source_level1_unambiguous = o4$level1_unambiguous[j], source_level1_lineage = o4$level1_lineage[j],
                  target_gene = o4$gene[j], zinb_zero_beta = o4$zinb.res.zero.beta[j], zinb_zero_fdr = o4$zinb.res.zero.p.All.fdr[j],
                  zinb_count_beta = o4$zinb.res.count.beta[j], zinb_count_fdr = o4$zinb.res.count.p.All.fdr[j],
                  spearman_rho = o4$spearman.rho[j], spearman_p = o4$spearman.p[j], min_fdr = o4$min_fdr[j], dist_peak_to_gene_hg19 = o4$DisPeak2Gene[j])
  # Does the Open4Gene source cell type match a population the candidate's peak is specific to?
  d[, celltype_concordance := ifelse(is.na(source_level1_unambiguous), NA, mapply(function(l, k) l %in% cand_ct[[k]], source_level1_unambiguous, i))]
  d[, celltype_concordance_lineage := ifelse(is.na(source_level1_lineage), NA, mapply(function(l, k) l %in% cand_ct[[k]], source_level1_lineage, i))]
  d
}
per_gene <- function(d) d[, .(n_links = .N, source_celltypes = paste(sort(unique(open4gene_source_celltype)), collapse = ";"),
                               n_celltype_resolved_links = sum(open4gene_source_celltype != "AllCell"),
                               best_min_fdr = min(min_fdr, na.rm = TRUE), max_abs_spearman_rho = max(abs(spearman_rho), na.rm = TRUE),
                               celltype_concordance = if (all(is.na(celltype_concordance))) NA else any(celltype_concordance, na.rm = TRUE),
                               celltype_concordance_lineage = if (all(is.na(celltype_concordance_lineage))) NA else any(celltype_concordance_lineage, na.rm = TRUE),
                               concordant_source_celltypes = paste(sort(unique(open4gene_source_celltype[!is.na(celltype_concordance) & celltype_concordance])), collapse = ";")),
                           by = .(rsID, gene = target_gene)]
direct <- open4gene_rows(ov_direct)
write_tsv(direct, arg(args, "out_open4gene_matches"))
direct_p <- per_gene(direct)

pkd <- parse_archr_peaks(cand$overlapping_peak); apk <- GRanges(pkd$chr, IRanges(pkd$start, pkd$end))
ov_context <- findOverlaps(apk, o4gr)
ov_context <- ov_context[!paste(queryHits(ov_context), subjectHits(ov_context)) %in% paste(queryHits(ov_direct), subjectHits(ov_direct))]
context <- open4gene_rows(ov_context)
write_tsv(context, arg(args, "out_open4gene_peak_context"))
context_p <- per_gene(context)

# ---- GTEx v10 kidney cortex: exact eVariant match (chr, position, REF, ALT) ----
ev <- readRDS(arg(args, "gtex_evariants")); evd <- as.data.table(mcols(ev))
pairs <- as.data.table(readRDS(arg(args, "gtex_pairs")))
evd[, key_exact := paste(chr, pos, ref, alt)]; evd[, key_pos := paste(chr, pos)]
cand[, key_exact := paste(chr, hg38_position, REF, ALT)]; cand[, key_swap := paste(chr, hg38_position, ALT, REF)]; cand[, key_pos := paste(chr, hg38_position)]
cand[, gtex_allele_match := fifelse(key_exact %in% evd$key_exact, "EXACT_REF_ALT", fifelse(key_swap %in% evd$key_exact, "REF_ALT_SWAPPED",
                             fifelse(key_pos %in% evd$key_pos, "POSITION_ONLY_different_alleles", "NO_MATCH")))]
# Candidate and GTEx REF alleles are both checked against the genome, so a swapped match would point to an error.
if (any(cand$gtex_allele_match == "REF_ALT_SWAPPED")) log_msg("warning: REF/ALT-swapped GTEx matches found")
gm <- merge(cand[gtex_allele_match %in% c("EXACT_REF_ALT", "REF_ALT_SWAPPED"),
                 .(rsID, PIP, locus, chr, hg38_position, REF, ALT, atlas_peak = overlapping_peak, candidate_cell_types = specific_cell_types, gtex_allele_match,
                   key_join = fifelse(gtex_allele_match == "EXACT_REF_ALT", key_exact, key_swap), effect_allele_hg38, egfr_direction_of_ALT, gwas_z)],
            evd[, .(key_join = key_exact, gtex_variant_id = variant_id, gtex_variant_class = variant_class, n_egenes, gtex_lead_rsid = lead_rsid)], by = "key_join")
gp <- merge(gm, pairs[, .(gtex_variant_id = variant_id, eGene = gene_name, eGene_id = gene_id, eGene_biotype = biotype, gtex_slope = slope, gtex_slope_se = slope_se,
                          gtex_pval_nominal = pval_nominal, gtex_gene_qval = gene_qval, gtex_tss_distance = tss_distance, gtex_af = af, alt_allele_direction, is_lead_variant_for_gene)],
            by = "gtex_variant_id", allow.cartesian = TRUE)
gp[, rsid_corroborated := ifelse(is.na(gtex_lead_rsid), NA, mapply(function(a, b) a %in% strsplit(b, ";")[[1]], rsID, gtex_lead_rsid))]
gp[, expression_direction_of_eGFR_raising_allele := fifelse(egfr_direction_of_ALT * sign(gtex_slope) > 0, "higher_expression",
                                                    fifelse(egfr_direction_of_ALT * sign(gtex_slope) < 0, "lower_expression", "undetermined"))]
setorder(gp, -PIP, rsID, eGene)
write_tsv(gp, arg(args, "out_gtex_matches"))
gtex_p <- gp[, .(n_gtex_pairs = .N, gtex_variant_id = gtex_variant_id[1], gtex_slope = gtex_slope[1], gtex_pval_nominal = gtex_pval_nominal[1], gtex_gene_qval = gtex_gene_qval[1],
                 gtex_allele_match = gtex_allele_match[1], expression_direction_of_eGFR_raising_allele = expression_direction_of_eGFR_raising_allele[1],
                 rsid_corroborated = rsid_corroborated[1]), by = .(rsID, gene = eGene)]
log_msg(sprintf("GTEx: %d candidates are significant kidney-cortex eVariants (%d candidate-eGene pairs)", length(unique(gp$rsID)), nrow(gtex_p)))

# ---- Proximity: nearest protein-coding transcript TSS (ties kept) ----
ann <- readRDS(arg(args, "gene_annotation"))
tss_pc <- protein_coding_tss(ann)
hits <- nearest(cgr, tss_pc, select = "all", ignore.strand = TRUE)
dd <- distance(cgr[queryHits(hits)], tss_pc[subjectHits(hits)], ignore.strand = TRUE)
near_gene <- tapply(tss_pc$gene_name[subjectHits(hits)], queryHits(hits), function(x) paste(sort(unique(x)), collapse = ";"))
near_dist <- tapply(dd, queryHits(hits), min)
near <- data.table(rsID = cand$rsID, nearest_protein_coding_gene = NA_character_, distance_to_nearest_protein_coding_tss = NA_integer_)
near$nearest_protein_coding_gene[as.integer(names(near_gene))] <- near_gene
near$distance_to_nearest_protein_coding_tss[as.integer(names(near_dist))] <- as.integer(near_dist)
near[, promoter_proximal_le_2kb := distance_to_nearest_protein_coding_tss <= promoter_window_bp]
ovg <- findOverlaps(cgr, ann$genes)
near[, intragenic := seq_along(cgr) %in% queryHits(ovg)]
is_pc <- ann$genes$gene_biotype[subjectHits(ovg)] == "protein_coding"
intra_pc <- tapply(ann$genes$gene_name[subjectHits(ovg)][is_pc], queryHits(ovg)[is_pc], function(x) paste(sort(unique(x)), collapse = ";"))
near[, intragenic_protein_coding_genes := NA_character_]; near$intragenic_protein_coding_genes[as.integer(names(intra_pc))] <- intra_pc
stopifnot(!anyNA(near$nearest_protein_coding_gene))
near_p <- near[, .(gene = unlist(strsplit(nearest_protein_coding_gene, ";")), distance_to_nearest_protein_coding_tss, promoter_proximal_le_2kb), by = rsID]

# ---- Candidate-gene pairs: boolean flags and a count of gene-linking sources (no weighted score) ----
keys <- unique(rbind(direct_p[, .(rsID, gene)], context_p[, .(rsID, gene)], gtex_p[, .(rsID, gene)], near_p[, .(rsID, gene)]))
long <- merge(keys, cand[, .(rsID, chr, hg38_position, REF, ALT, PIP, locus, credible_set, atlas_peak = overlapping_peak,
                             candidate_cell_types = specific_cell_types, n_candidate_cell_types = n_specific_cell_types)], by = "rsID")
long <- merge(long, direct_p[, .(rsID, gene, open4gene_variant_in_linked_peak = TRUE, open4gene_n_links = n_links, open4gene_source_celltypes = source_celltypes,
                                 open4gene_n_celltype_resolved_links = n_celltype_resolved_links, open4gene_best_min_fdr = best_min_fdr,
                                 open4gene_max_abs_spearman_rho = max_abs_spearman_rho, celltype_concordant_open4gene = celltype_concordance,
                                 celltype_concordant_open4gene_lineage = celltype_concordance_lineage, open4gene_concordant_source_celltypes = concordant_source_celltypes)],
              by = c("rsID", "gene"), all.x = TRUE)
long <- merge(long, context_p[, .(rsID, gene, open4gene_atlas_peak_overlap = TRUE, open4gene_atlas_overlap_source_celltypes = source_celltypes)], by = c("rsID", "gene"), all.x = TRUE)
long <- merge(long, gtex_p[, .(rsID, gene, gtex_exact_kidney_eqtl = TRUE, gtex_variant_id, gtex_slope, gtex_pval_nominal, gtex_gene_qval, gtex_allele_match,
                               expression_direction_of_eGFR_raising_allele, gtex_rsid_corroborated = rsid_corroborated)], by = c("rsID", "gene"), all.x = TRUE)
long <- merge(long, near_p[, .(rsID, gene, nearest_gene = TRUE, distance_to_nearest_protein_coding_tss, promoter_proximal_le_2kb)], by = c("rsID", "gene"), all.x = TRUE)
for (v in c("open4gene_variant_in_linked_peak", "open4gene_atlas_peak_overlap", "gtex_exact_kidney_eqtl", "nearest_gene")) long[is.na(get(v)), (v) := FALSE]
long[, n_functional_sources := as.integer(open4gene_variant_in_linked_peak) + as.integer(gtex_exact_kidney_eqtl)]
long[, functional_sources := sub(";$", "", paste0(ifelse(open4gene_variant_in_linked_peak, "Open4Gene_direct;", ""), ifelse(gtex_exact_kidney_eqtl, "GTEx_kidney_exact_eQTL;", "")))]
long[, pair_category := fifelse(n_functional_sources >= 2, "CONVERGENT_FUNCTIONAL", fifelse(n_functional_sources == 1, "FUNCTIONAL_SINGLE_SOURCE",
                         fifelse(nearest_gene, "PROXIMITY_ONLY", "SECONDARY_ONLY_atlas_peak_overlap")))]
gsym <- ann$genes$gene_name
long[, gene_id_ensembl108 := ann$genes$gene_id[match(gene, gsym)]]; long[, gene_biotype_ensembl108 := ann$genes$gene_biotype[match(gene, gsym)]]
setorder(long, -n_functional_sources, -PIP, rsID, gene)
setcolorder(long, c("rsID", "gene", "pair_category", "n_functional_sources", "functional_sources", "open4gene_variant_in_linked_peak", "gtex_exact_kidney_eqtl", "nearest_gene",
                    "open4gene_atlas_peak_overlap", "celltype_concordant_open4gene", "celltype_concordant_open4gene_lineage", "PIP", "locus", "credible_set", "chr", "hg38_position", "REF", "ALT", "atlas_peak", "candidate_cell_types"))
stopifnot(length(unique(long$rsID)) == n_cand, !any(duplicated(long[, .(rsID, gene)])))
write_tsv(long, arg(args, "out_long"))
conv <- long[pair_category == "CONVERGENT_FUNCTIONAL"]
write_tsv(conv, arg(args, "out_convergent"))

# ---- One row per candidate ----
cat_rank <- c(CONVERGENT_FUNCTIONAL = 4L, FUNCTIONAL_SINGLE_SOURCE = 3L, PROXIMITY_ONLY = 2L, SECONDARY_ONLY_atlas_peak_overlap = 1L)
pl <- function(x) paste(sort(unique(x)), collapse = ";")
sm <- long[, .(n_open4gene_direct_genes = sum(open4gene_variant_in_linked_peak), open4gene_direct_genes = pl(gene[open4gene_variant_in_linked_peak]),
               n_open4gene_atlas_overlap_genes = sum(open4gene_atlas_peak_overlap), open4gene_atlas_overlap_genes = pl(gene[open4gene_atlas_peak_overlap]),
               n_gtex_egenes = sum(gtex_exact_kidney_eqtl), gtex_egenes = pl(gene[gtex_exact_kidney_eqtl]),
               n_functional_target_genes = sum(n_functional_sources >= 1), functional_target_genes = pl(gene[n_functional_sources >= 1]),
               n_convergent_genes = sum(n_functional_sources >= 2), convergent_genes = pl(gene[n_functional_sources >= 2]),
               max_n_functional_sources = max(n_functional_sources),
               celltype_concordant_open4gene_any = if (all(is.na(celltype_concordant_open4gene[open4gene_variant_in_linked_peak]))) NA else any(celltype_concordant_open4gene[open4gene_variant_in_linked_peak], na.rm = TRUE),
               candidate_category = pair_category[which.max(cat_rank[pair_category])]), by = rsID]
summ <- merge(cand[, .(rsID, chr, hg38_position, REF, ALT, PIP, locus, credible_set, atlas_peak = overlapping_peak, candidate_cell_types = specific_cell_types,
                       n_candidate_cell_types = n_specific_cell_types, gtex_allele_match)], sm, by = "rsID")
summ <- merge(summ, near[, .(rsID, nearest_protein_coding_gene, distance_to_nearest_protein_coding_tss, promoter_proximal_le_2kb, intragenic, intragenic_protein_coding_genes)], by = "rsID")
setorder(summ, -max_n_functional_sources, -PIP, rsID)
stopifnot(nrow(summ) == n_cand)
write_tsv(summ, arg(args, "out_summary"))

# ---- Evidence counts ----
ec <- data.table(statistic = c("candidates", "candidates_with_ge1_functional_target_gene", "candidates_with_open4gene_direct_support", "candidates_with_exact_gtex_kidney_eqtl",
                               "candidates_with_convergent_ge2_source_gene", "candidates_proximity_only", "candidates_with_gtex_position_only_allele_mismatch",
                               "candidate_gene_pairs_total", "pairs_convergent_functional", "pairs_functional_single_source", "pairs_proximity_only", "pairs_secondary_only_atlas_peak_overlap",
                               "unique_functional_target_genes", "unique_convergent_genes", "unique_open4gene_direct_genes", "unique_gtex_egenes", "candidates_promoter_proximal_le_2kb", "candidates_intragenic"),
                 n = c(nrow(summ), sum(summ$n_functional_target_genes >= 1), sum(summ$n_open4gene_direct_genes >= 1), sum(summ$n_gtex_egenes >= 1),
                       sum(summ$n_convergent_genes >= 1), sum(summ$candidate_category == "PROXIMITY_ONLY"), sum(summ$gtex_allele_match == "POSITION_ONLY_different_alleles"),
                       nrow(long), sum(long$pair_category == "CONVERGENT_FUNCTIONAL"), sum(long$pair_category == "FUNCTIONAL_SINGLE_SOURCE"), sum(long$pair_category == "PROXIMITY_ONLY"), sum(long$pair_category == "SECONDARY_ONLY_atlas_peak_overlap"),
                       length(unique(long$gene[long$n_functional_sources >= 1])), length(unique(conv$gene)), length(unique(long$gene[long$open4gene_variant_in_linked_peak])),
                       length(unique(long$gene[long$gtex_exact_kidney_eqtl])), sum(near$promoter_proximal_le_2kb), sum(near$intragenic)))
ec[, percent_of_candidates := ifelse(grepl("^candidates", statistic), round(100 * n / n_cand, 1), NA)]
write_tsv(ec, arg(args, "out_counts"))
log_msg(sprintf("target genes: %d/%d candidates with >= 1 gene-linking source; %d convergent pairs", ec$n[2], ec$n[1], ec$n[9]))

# (C) Sergio Alejandro Gómez
# Aim:     Compare the predicted TF motif disruption of the candidates attributed to PT with that of the candidates attributed
#          to TAL. 
# Inputs:  tables/candidate_variants.tsv; tables/motif_disruption_by_variant.tsv and tables/motif_disruption_by_tf.tsv
#          (14_motif_disruption.R).
# Outputs: tables/motif_pt_tal_group_counts.tsv, tables/motif_pt_tal_by_tf.tsv, tables/motif_pt_tal_top30_tfs.tsv,
#          tables/motif_pt_tal_similarity.tsv, tables/motif_pt_tal_topk_overlap.tsv, tables/motif_tf_load_by_group.tsv
# Usage: Rscript workflow/scripts/genetics/15_motif_pt_vs_tal.R --candidates=results/genetics/with_ptC5/tables/candidate_variants.tsv --by_variant=results/genetics/with_ptC5/tables/motif_disruption_by_variant.tsv --by_tf=results/genetics/with_ptC5/tables/motif_disruption_by_tf.tsv --out_group_counts=results/genetics/with_ptC5/tables/motif_pt_tal_group_counts.tsv --out_by_tf=results/genetics/with_ptC5/tables/motif_pt_tal_by_tf.tsv --out_top30=results/genetics/with_ptC5/tables/motif_pt_tal_top30_tfs.tsv --out_similarity=results/genetics/with_ptC5/tables/motif_pt_tal_similarity.tsv --out_topk=results/genetics/with_ptC5/tables/motif_pt_tal_topk_overlap.tsv --out_tf_load=results/genetics/with_ptC5/tables/motif_tf_load_by_group.tsv

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "motif.R"))

args <- parse_args()
top_n      <- as.integer(arg(args, "top_n_tfs", "30"))
top_k      <- as.integer(strsplit(arg(args, "top_k", "10,20,30"), ",", fixed = TRUE)[[1]])
diff_band  <- as.numeric(arg(args, "diff_band", "0.10"))       # |PT - TAL fraction| below this: SIMILAR
similar_rho        <- as.numeric(arg(args, "similar_rho", "0.80"))
partly_similar_rho <- as.numeric(arg(args, "partly_similar_rho", "0.50"))

read_chr <- function(path) fread(path, sep = "\t", header = TRUE, colClasses = "character", na.strings = NULL)
cand <- read_chr(arg(args, "candidates"))
by_variant <- read_chr(arg(args, "by_variant"))
by_tf <- read_chr(arg(args, "by_tf"))
n_cand <- nrow(cand)

# ---- candidate x TF membership of the strong effects ------------------------------------------------------------
tok <- strsplit(by_variant$tfs, ";", fixed = TRUE)
long <- data.table(SNP_id = rep(by_variant$SNP_id, lengths(tok)), geneSymbol = unlist(tok))
setorder(long, SNP_id, geneSymbol)
universe <- sort(unique(long$geneSymbol))
stopifnot(setequal(universe, by_tf$geneSymbol))
pooled <- by_tf[, .(geneSymbol, n_motif_matrices = as.integer(n_motif_matrices), pooled_n_distinct_snps = as.integer(n_distinct_disrupted_snps),
                    pooled_denominator = n_cand, pooled_fraction = as.integer(n_distinct_disrupted_snps) / n_cand,
                    pooled_rank_canonical = as.integer(rank_by_distinct_snps))]

# ---- PT and TAL groups (from the populations each candidate is attributed to) ------------------------------------
sets <- strsplit(cand$specific_cell_types, ";", fixed = TRUE)
has <- function(p) vapply(sets, function(x) p %in% x, TRUE)
g <- data.table(SNP_id = paste0(cand$chr, ":", cand$hg38_position, ":", cand$REF, ":", cand$ALT), locus = cand$locus,
                is_PT = has("PT"), is_TAL = has("TAL"), is_DCT = has("DCT"), is_IC = has("IC"), is_CNT_PC = has("CNT_PC"),
                PT_ONLY = cand$specific_cell_types == "PT", TAL_ONLY = cand$specific_cell_types == "TAL", multi_population = lengths(sets) > 1L)
g[, `:=`(PT_AND_TAL = is_PT & is_TAL, PT_NOT_TAL = is_PT & !is_TAL, TAL_NOT_PT = !is_PT & is_TAL, NEITHER_PT_NOR_TAL = !is_PT & !is_TAL)]
stopifnot(all(long$SNP_id %in% g$SNP_id))

group_def <- list(
  ALL_CANDIDATES     = list(rep(TRUE, n_cand), "all canonical candidates (PIP > 0.1, in an enriched-population-specific peak)"),
  PT_ATTRIBUTED      = list(g$is_PT, "PT is among the enriched specific-cell-type attributions"),
  TAL_ATTRIBUTED     = list(g$is_TAL, "TAL is among the enriched specific-cell-type attributions"),
  DCT_ATTRIBUTED     = list(g$is_DCT, "DCT is among the attributions"),
  IC_ATTRIBUTED      = list(g$is_IC, "IC is among the attributions"),
  CNT_PC_ATTRIBUTED  = list(g$is_CNT_PC, "CNT_PC is among the attributions"),
  MULTI_POPULATION   = list(g$multi_population, "attributed to > 1 enriched population"),
  PT_ONLY            = list(g$PT_ONLY, "attribution is exactly PT and no other population"),
  TAL_ONLY           = list(g$TAL_ONLY, "attribution is exactly TAL and no other population"),
  PT_AND_TAL         = list(g$PT_AND_TAL, "attributed to both PT and TAL (shared; contributes to both primary strata)"),
  PT_NOT_TAL         = list(g$PT_NOT_TAL, "attributed to PT but not to TAL (other populations allowed)"),
  TAL_NOT_PT         = list(g$TAL_NOT_PT, "attributed to TAL but not to PT (other populations allowed)"),
  NEITHER_PT_NOR_TAL = list(g$NEITHER_PT_NOR_TAL, "attributed to neither PT nor TAL"))
group_counts <- rbindlist(lapply(names(group_def), function(n) {
  fl <- group_def[[n]][[1]]
  data.table(group = n, definition = group_def[[n]][[2]], n_candidates = sum(fl), n_loci = uniqueN(g$locus[fl]),
             fraction_of_all_candidates = sum(fl) / n_cand)
}))
stopifnot(group_counts[group == "PT_ATTRIBUTED", n_candidates] > 0, group_counts[group == "TAL_ATTRIBUTED", n_candidates] > 0)

# ---- fraction of disrupted candidates per TF, in each group -------------------------------------------------------
strata <- list(PT_ATTRIBUTED = g[is_PT == TRUE, SNP_id], TAL_ATTRIBUTED = g[is_TAL == TRUE, SNP_id],
               PT_ONLY = g[PT_ONLY == TRUE, SNP_id], TAL_ONLY = g[TAL_ONLY == TRUE, SNP_id],
               PT_AND_TAL = g[PT_AND_TAL == TRUE, SNP_id], PT_NOT_TAL = g[PT_NOT_TAL == TRUE, SNP_id], TAL_NOT_PT = g[TAL_NOT_PT == TRUE, SNP_id],
               ALL_CANDIDATES = g$SNP_id)
prev <- lapply(c("PT_ATTRIBUTED", "TAL_ATTRIBUTED", "PT_AND_TAL", "PT_NOT_TAL", "TAL_NOT_PT"),
               function(s) setorder(stratum_prevalence(long, strata[[s]], universe, s), geneSymbol))
names(prev) <- c("pt", "tal", "sh", "pnt", "tnp")
pt <- prev$pt; tal <- prev$tal; sh <- prev$sh; pnt <- prev$pnt; tnp <- prev$tnp
stopifnot(all(pt$n_distinct_disrupted_snps == pnt$n_distinct_disrupted_snps + sh$n_distinct_disrupted_snps),
          all(tal$n_distinct_disrupted_snps == tnp$n_distinct_disrupted_snps + sh$n_distinct_disrupted_snps))

by_tf_pt_tal <- merge(data.table(geneSymbol = pt$geneSymbol), pooled, by = "geneSymbol")
setorder(by_tf_pt_tal, geneSymbol)
by_tf_pt_tal[, `:=`(PT_n_candidates = pt$n_candidates_group, PT_n_distinct_disrupted_snps = pt$n_distinct_disrupted_snps, PT_fraction = pt$fraction_disrupted,
                    PT_rank_by_count = pt$rank_by_count, PT_rank_by_fraction = pt$rank_by_fraction,
                    TAL_n_candidates = tal$n_candidates_group, TAL_n_distinct_disrupted_snps = tal$n_distinct_disrupted_snps, TAL_fraction = tal$fraction_disrupted,
                    TAL_rank_by_count = tal$rank_by_count, TAL_rank_by_fraction = tal$rank_by_fraction,
                    SHARED_PT_AND_TAL_n_candidates = sh$n_candidates_group, SHARED_PT_AND_TAL_n_disrupted = sh$n_distinct_disrupted_snps,
                    PT_NOT_TAL_n_candidates = pnt$n_candidates_group, PT_NOT_TAL_n_disrupted = pnt$n_distinct_disrupted_snps,
                    TAL_NOT_PT_n_candidates = tnp$n_candidates_group, TAL_NOT_PT_n_disrupted = tnp$n_distinct_disrupted_snps)]
by_tf_pt_tal[, diff_PT_minus_TAL_fraction := PT_fraction - TAL_fraction]
by_tf_pt_tal[, abs_diff_fraction := abs(diff_PT_minus_TAL_fraction)]
by_tf_pt_tal[, rank_diff_PT_minus_TAL := PT_rank_by_count - TAL_rank_by_count]
by_tf_pt_tal[, present_in_PT := PT_n_distinct_disrupted_snps > 0L][, present_in_TAL := TAL_n_distinct_disrupted_snps > 0L]
by_tf_pt_tal[, presence_class := fifelse(present_in_PT & present_in_TAL, "OBSERVED_IN_BOTH", fifelse(present_in_PT, "OBSERVED_ONLY_IN_PT",
                                  fifelse(present_in_TAL, "OBSERVED_ONLY_IN_TAL", "ABSENT_IN_BOTH")))]
by_tf_pt_tal[, diff_band_class := fifelse(abs(diff_PT_minus_TAL_fraction) < diff_band, "SIMILAR",
                                  fifelse(diff_PT_minus_TAL_fraction >= diff_band, "PT_RELATIVELY_MORE_FREQUENT", "TAL_RELATIVELY_MORE_FREQUENT"))]
by_tf_pt_tal[, diff_percentile_among_all_tfs := 100 * (frank(diff_PT_minus_TAL_fraction, ties.method = "average") - 0.5) / .N]
by_tf_pt_tal[, comparison_type := "NON_EXCLUSIVE_OVERLAPPING_STRATA_DESCRIPTIVE_ONLY"]
setorder(by_tf_pt_tal, pooled_rank_canonical)
setcolorder(by_tf_pt_tal, c("geneSymbol", "n_motif_matrices", "pooled_rank_canonical", "pooled_n_distinct_snps", "pooled_denominator", "pooled_fraction"))

# ---- similarity of the two profiles, over all TFs ------------------------------------------------------------------
sim <- similarity_metrics(pt, tal, "PRIMARY_PT_ATTRIBUTED_vs_TAL_ATTRIBUTED", "candidate_variant", "PT_ATTRIBUTED", "TAL_ATTRIBUTED",
                          n_shared_units = length(strata$PT_AND_TAL), diff_band = diff_band, similar_rho = similar_rho,
                          partly_similar_rho = partly_similar_rho, top_k = top_k)
log_msg(sprintf("PT (n = %d) vs TAL (n = %d), %d shared: Spearman rho of TF fractions = %.4f over %d TFs (%s)",
                length(strata$PT_ATTRIBUTED), length(strata$TAL_ATTRIBUTED), length(strata$PT_AND_TAL),
                sim$summary$spearman_rho_fraction, nrow(pt), sim$summary$similarity_band))

# ---- the top TFs of all candidates, with their PT and TAL fractions -------------------------------------------------
top30 <- by_tf_pt_tal[pooled_rank_canonical <= top_n, .(pooled_rank_canonical, geneSymbol, n_motif_matrices, pooled_n_distinct_snps, pooled_denominator,
                                                       pooled_fraction, PT_n_candidates, PT_n_distinct_disrupted_snps, PT_fraction, PT_rank_by_count,
                                                       TAL_n_candidates, TAL_n_distinct_disrupted_snps, TAL_fraction, TAL_rank_by_count,
                                                       SHARED_PT_AND_TAL_n_disrupted, diff_PT_minus_TAL_fraction, diff_band_class)]
setorder(top30, pooled_rank_canonical)

# ---- number of distinct TFs disrupted per candidate, by group --------------------------------------------------------
load_snp <- long[, .(n_distinct_tfs = uniqueN(geneSymbol)), by = SNP_id]
tf_load <- rbindlist(lapply(names(strata), function(s) {
  x <- load_snp[SNP_id %in% strata[[s]], n_distinct_tfs]
  data.table(stratum = s, n_candidates = length(x), median_n_distinct_tfs = as.numeric(median(x)), mean_n_distinct_tfs = mean(x),
             min_n_distinct_tfs = min(x), max_n_distinct_tfs = max(x), mean_fraction_of_universe = mean(x) / length(universe))
}))

# ---- write -------------------------------------------------------------------------------------------------------------
write_exact <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.table(df), path, sep = "\t", quote = FALSE, na = "NA", eol = "\n", scipen = 999)
  log_msg("wrote ", path, " (", nrow(df), " rows)")
}
write_exact(group_counts, arg(args, "out_group_counts"))
write_exact(by_tf_pt_tal, arg(args, "out_by_tf"))
write_exact(top30, arg(args, "out_top30"))
write_exact(sim$summary, arg(args, "out_similarity"))
write_exact(sim$topk, arg(args, "out_topk"))
write_exact(tf_load, arg(args, "out_tf_load"))

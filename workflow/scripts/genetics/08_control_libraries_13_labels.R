# (C) Sergio Alejandro Gómez
# Aim: Sensitivity analysis repeating the eGFR enrichment on the 13 cell labels, which keep the proximal-tubule states apart (PT, adaptive aPT and
#   PT_aPT, degenerative dPT), in all 43 libraries and after removing the libraries of the control, AKI-control, KIN and
#   BKVAN groups (32 CKD biopsy libraries).
# Inputs: pseudobulk of all nuclei (full_cohort, 02), cluster -> 13-label table, fine-mapped eGFR variants,
#   eGFR seed summary of the same version (04; reference for background-seed noise).
# Outputs: per-label comparison, nuclei per label in both runs, libraries kept and removed, summary metrics.
# Usage: Rscript workflow/scripts/genetics/08_control_libraries_13_labels.R --pseudobulk=results/genetics/full_cohort/objects/pseudobulk.rds --labels13=resources/genetics/atac_13_labels.csv --variants=results/genetics/egfr_finemapping/egfr_finemapped_variants_hg38.tsv --egfr_by_seed=results/genetics/with_ptC5/tables/egfr_enrichment_by_seed.tsv --include_ptC5=true --out_comparison=results/genetics/with_ptC5/tables/control_libraries_13_labels.tsv --out_nuclei=results/genetics/with_ptC5/tables/control_libraries_13_labels_nuclei.tsv --out_libraries=results/genetics/with_ptC5/tables/control_libraries_13_labels_libraries.tsv --out_summary=results/genetics/with_ptC5/tables/control_libraries_13_labels_summary.tsv

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "accessibility.R"))
source(file.path(script_dir, "R", "gchromvar.R"))

args <- parse_args()
version <- if (arg_true(args, "include_ptC5", NULL)) "with_ptC5" else "without_ptC5"
seed <- as.integer(arg(args, "seed", "2024"))
exclude_clusters <- if (version == "with_ptC5") character(0) else "pt_C5"

EXCLUDED_GROUPS <- c("control", "control_AKI", "KIN", "BKVAN")
EXCLUSION_SET <- "control_controlAKI_KIN_BKVAN"   # label of this set of disease groups in the tables
LABELS_13 <- c("aPT", "CNT_PC", "DCT", "dPT", "EC", "Fib_VSMC", "IC", "Immune", "PEC", "Pod", "PT", "PT_aPT", "TAL")
STATE_GROUP <- c(PT = "homeostatic", TAL = "homeostatic", DCT = "homeostatic", CNT_PC = "homeostatic", IC = "homeostatic",
                 PEC = "homeostatic", EC = "homeostatic", Immune = "homeostatic", Fib_VSMC = "homeostatic", Pod = "homeostatic",
                 aPT = "adaptive", PT_aPT = "adaptive", dPT = "degenerative")
STATE_NOTE <- c(PT_aPT = "mixed PT / adaptive-PT cluster (pt_C4); 84 % of its nuclei are control_AKI",
                aPT = "adaptive PT (pt_C6)", dPT = "degenerative PT (pt_C2)")

full <- readRDS(arg(args, "pseudobulk"))
lab <- read.csv(arg(args, "labels13"), stringsAsFactors = FALSE, strip.white = TRUE)
label_of <- setNames(lab$label13, lab$cluster)
variants <- read_egfr_variants(arg(args, "variants"))
egfr_by_seed <- read_tsv(arg(args, "egfr_by_seed"))
peaks_gr <- peak_granges(full$peaks)

# ---- 13-label counts of the two runs ----
g <- full$groups
stopifnot(all(g$cluster %in% names(label_of)), setequal(unique(label_of), LABELS_13))
use <- !(g$cluster %in% exclude_clusters)
meta <- data.frame(sample = g$sample[use], label = unname(label_of[g$cluster[use]]), diseaseNames = g$diseaseNames[use],
                   n_cells = as.integer(g$n_nuclei[use]), stringsAsFactors = FALSE)
counts <- full$counts[, use, drop = FALSE]
keep <- !(meta$diseaseNames %in% EXCLUDED_GROUPS)
label_counts <- function(k) {
  ind <- sparseMatrix(i = seq_len(sum(k)), j = match(meta$label[k], LABELS_13), x = 1, dims = c(sum(k), length(LABELS_13)))
  out <- as.matrix(counts[, k, drop = FALSE] %*% ind); dimnames(out) <- list(rownames(counts), LABELS_13)
  storage.mode(out) <- "double"
  out
}
cnt_b <- label_counts(rep(TRUE, nrow(meta)))
cnt_s <- label_counts(keep)

cells_b <- tapply(meta$n_cells, meta$label, sum); cells_s <- tapply(meta$n_cells[keep], meta$label[keep], sum)
samp_b <- tapply(meta$sample, meta$label, function(x) length(unique(x)))
samp_s <- tapply(meta$sample[keep], meta$label[keep], function(x) length(unique(x)))
pcc <- data.frame(population = LABELS_13, state_group = unname(STATE_GROUP[LABELS_13]),
                  n_nuclei_baseline = as.integer(cells_b[LABELS_13]), n_nuclei_sensitivity = as.integer(cells_s[LABELS_13]),
                  n_nuclei_removed = as.integer(cells_b[LABELS_13] - cells_s[LABELS_13]),
                  frac_nuclei_retained = as.numeric(cells_s[LABELS_13] / cells_b[LABELS_13]),
                  n_samples_with_nuclei_baseline = as.integer(samp_b[LABELS_13]), n_samples_with_nuclei_sensitivity = as.integer(samp_s[LABELS_13]),
                  total_counts_baseline = colSums(cnt_b)[LABELS_13], total_counts_sensitivity = colSums(cnt_s)[LABELS_13],
                  zero_nuclei_after_exclusion = as.integer(cells_s[LABELS_13]) == 0, stringsAsFactors = FALSE)
pcc <- pcc[order(-pcc$n_nuclei_baseline), ]; rownames(pcc) <- NULL
write_tsv(pcc, arg(args, "out_nuclei"))

excluded <- sort(unique(meta$sample[!keep])); retained <- sort(unique(meta$sample[keep]))
libs <- data.frame(item = c("exclusion_set", "excluded_disease_labels", "excluded_samples", "n_excluded_samples", "retained_samples", "n_retained_samples",
                            "n_nuclei_baseline", "n_nuclei_retained", "n_pseudobulk_groups_baseline", "n_pseudobulk_groups_retained", "ptc5_version"),
                   value = c(EXCLUSION_SET, paste(EXCLUDED_GROUPS, collapse = ";"), paste(excluded, collapse = ";"), length(excluded),
                             paste(retained, collapse = ";"), length(retained), sum(meta$n_cells), sum(meta$n_cells[keep]), nrow(meta), sum(keep), version),
                   stringsAsFactors = FALSE)
write_tsv(libs, arg(args, "out_libraries"))
log_msg(sprintf("all libraries: %d / %d nuclei; after removal: %d / %d nuclei", length(unique(meta$sample)), sum(meta$n_cells), length(retained), sum(meta$n_cells[keep])))

# ---- gchromVAR of both runs ----
b <- run_gchromvar((to_cpm(cnt_b) + 1)[names(peaks_gr), ], peaks_gr, variants, "Kidney_eGFR_CARMA", seed)$result
s <- run_gchromvar((to_cpm(cnt_s) + 1)[names(peaks_gr), ], peaks_gr, variants, "Kidney_eGFR_CARMA", seed)$result

# ---- per-label comparison ----
call_of <- function(x) ifelse(grepl("^enriched", x), "enriched", ifelse(grepl("^depleted", x), "depleted", "ns"))
m <- match(b$cell_type, s$cell_type); k <- match(b$cell_type, pcc$population)
cmp <- data.frame(population = b$cell_type, cell_group = unname(STATE_GROUP[b$cell_type]),
                  baseline_Z = b$z, baseline_p = b$p, baseline_FDR = b$FDR, baseline_rank = b$rank,
                  sensitivity_Z = s$z[m], sensitivity_p = s$p[m], sensitivity_FDR = s$FDR[m], sensitivity_rank = s$rank[m], stringsAsFactors = FALSE)
cmp$delta_Z <- cmp$sensitivity_Z - cmp$baseline_Z; cmp$delta_rank <- cmp$sensitivity_rank - cmp$baseline_rank
cmp$baseline_significant <- call_of(b$significance); cmp$sensitivity_significant <- call_of(s$significance[m])
cmp$fdr_call_changed <- cmp$baseline_significant != cmp$sensitivity_significant
cmp$sign_reversal_at_fdr <- (cmp$baseline_significant == "enriched" & cmp$sensitivity_significant == "depleted") |
  (cmp$baseline_significant == "depleted" & cmp$sensitivity_significant == "enriched")
cmp$n_nuclei_baseline <- pcc$n_nuclei_baseline[k]; cmp$n_nuclei_sensitivity <- pcc$n_nuclei_sensitivity[k]; cmp$frac_nuclei_retained <- pcc$frac_nuclei_retained[k]
cmp$n_samples_baseline <- pcc$n_samples_with_nuclei_baseline[k]; cmp$n_samples_sensitivity <- pcc$n_samples_with_nuclei_sensitivity[k]
cmp$baseline_deviation <- b$deviation; cmp$sensitivity_deviation <- s$deviation[m]
cmp$cell_group_note <- ifelse(cmp$population %in% names(STATE_NOTE), unname(STATE_NOTE[cmp$population]), "cluster-level cell type, no injury-state split")
cmp <- cmp[order(cmp$baseline_rank), ]; rownames(cmp) <- NULL
cmp$ptc5_version <- version
write_tsv(cmp, arg(args, "out_comparison"))

# ---- summary metrics ----
labs_where <- function(x) paste(cmp$population[x], collapse = ";")
z_change <- function(p) sprintf("%.2f -> %.2f", cmp$baseline_Z[cmp$population == p], cmp$sensitivity_Z[cmp$population == p])
leading <- cmp$population[cmp$baseline_significant == "enriched" & cmp$cell_group == "homeostatic"]
lead <- cmp[cmp$population %in% leading, ]
lost <- cmp$population[cmp$baseline_significant == "enriched" & cmp$sensitivity_significant != "enriched"]
metrics <- list(
  exclusion_set = EXCLUSION_SET, excluded_labels = paste(EXCLUDED_GROUPS, collapse = ";"),
  n_samples_baseline = length(unique(meta$sample)), n_samples_sensitivity = length(retained),
  n_nuclei_baseline = sum(pcc$n_nuclei_baseline), n_nuclei_sensitivity = sum(pcc$n_nuclei_sensitivity),
  n_populations = length(LABELS_13), n_populations_zero_nuclei_after_exclusion = sum(pcc$zero_nuclei_after_exclusion),
  enriched_baseline = sum(cmp$baseline_significant == "enriched"), enriched_sensitivity = sum(cmp$sensitivity_significant == "enriched"),
  depleted_baseline = sum(cmp$baseline_significant == "depleted"), depleted_sensitivity = sum(cmp$sensitivity_significant == "depleted"),
  ns_baseline = sum(cmp$baseline_significant == "ns"), ns_sensitivity = sum(cmp$sensitivity_significant == "ns"),
  leading_set_definition = "homeostatic populations enriched at baseline (z > 0 & FDR < 0.05)",
  leading_set_baseline = paste(leading, collapse = ";"), leading_set_lost_enrichment = paste(intersect(lost, leading), collapse = ";"),
  enrichment_lost_any = paste(lost, collapse = ";"),
  enrichment_gained_any = labs_where(cmp$baseline_significant != "enriched" & cmp$sensitivity_significant == "enriched"),
  depletion_lost_any = labs_where(cmp$baseline_significant == "depleted" & cmp$sensitivity_significant != "depleted"),
  depletion_gained_any = labs_where(cmp$baseline_significant != "depleted" & cmp$sensitivity_significant == "depleted"),
  n_fdr_call_changes = sum(cmp$fdr_call_changed), sign_reversals_at_fdr = sum(cmp$sign_reversal_at_fdr),
  top2_baseline = labs_where(cmp$baseline_rank <= 2), top2_sensitivity = labs_where(cmp$sensitivity_rank <= 2),
  top3_baseline = labs_where(cmp$baseline_rank <= 3), top3_sensitivity = labs_where(cmp$sensitivity_rank <= 3),
  rank1_baseline = cmp$population[cmp$baseline_rank == 1], rank1_sensitivity = cmp$population[cmp$sensitivity_rank == 1],
  TAL_z_baseline_to_sensitivity = z_change("TAL"), PT_z_baseline_to_sensitivity = z_change("PT"), DCT_z_baseline_to_sensitivity = z_change("DCT"),
  aPT_z_baseline_to_sensitivity = z_change("aPT"), PT_aPT_z_baseline_to_sensitivity = z_change("PT_aPT"), dPT_z_baseline_to_sensitivity = z_change("dPT"),
  mean_abs_delta_z = sprintf("%.3f", mean(abs(cmp$delta_Z))), max_abs_delta_z = sprintf("%.3f", max(abs(cmp$delta_Z))),
  population_max_abs_delta_z = cmp$population[which.max(abs(cmp$delta_Z))],
  largest_increase_z = sprintf("%s %+.2f", cmp$population[which.max(cmp$delta_Z)], max(cmp$delta_Z)),
  largest_decrease_z = sprintf("%s %+.2f", cmp$population[which.min(cmp$delta_Z)], min(cmp$delta_Z)),
  spearman_rank_correlation_z = sprintf("%.3f", cor(cmp$baseline_Z, cmp$sensitivity_Z, method = "spearman")),
  pearson_correlation_z = sprintf("%.3f", cor(cmp$baseline_Z, cmp$sensitivity_Z)),
  sum_z_homeostatic_enriched_baseline = sprintf("%.2f", sum(lead$baseline_Z)), sum_z_homeostatic_enriched_sensitivity = sprintf("%.2f", sum(lead$sensitivity_Z)),
  mean_z_adaptive_baseline = sprintf("%.2f", mean(cmp$baseline_Z[cmp$cell_group == "adaptive"])),
  mean_z_adaptive_sensitivity = sprintf("%.2f", mean(cmp$sensitivity_Z[cmp$cell_group == "adaptive"])),
  # mean SD of the eGFR z over the six background seeds: the size of z changes expected from the background draw alone
  seed_noise_reference_mean_sd_z = sprintf("%.2f", mean(egfr_by_seed$z_sd_6seeds)),
  ptc5_version = version, peak_universe = PEAK_UNIVERSE_LABEL, variant_representation = "canonical_1bp")
summ <- data.frame(metric = names(metrics), value = vapply(metrics, as.character, ""), stringsAsFactors = FALSE)
write_tsv(summ, arg(args, "out_summary"))
log_msg(sprintf("13 labels (%s): Spearman rho %s; calls changed: %s", version, metrics$spearman_rank_correlation_z, metrics$n_fdr_call_changes))

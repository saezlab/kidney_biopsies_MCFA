# (C) Sergio Alejandro Gómez
# Aim: Check that the eGFR enrichment of the 10 populations does not depend on the five control libraries.
# Inputs: pseudobulk object (02), peak universe (03), fine-mapped eGFR variants.
# Outputs: per-population comparison of the two runs, every seed run of both, one-row summary.
# Usage: Rscript workflow/scripts/genetics/07_control_libraries_10_populations.R --pseudobulk=results/genetics/full_cohort/objects/pseudobulk.rds --peak_universe=results/genetics/with_ptC5/objects/peak_universe.rds --variants=results/genetics/egfr_finemapping/egfr_finemapped_variants_hg38.tsv --include_ptC5=true --out_comparison=results/genetics/with_ptC5/tables/control_libraries_10_populations.tsv --out_by_seed=results/genetics/with_ptC5/tables/control_libraries_10_populations_by_seed.tsv --out_summary=results/genetics/with_ptC5/tables/control_libraries_10_populations_summary.tsv

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "accessibility.R"))
source(file.path(script_dir, "R", "gchromvar.R"))

args <- parse_args()
version <- if (arg_true(args, "include_ptC5", NULL)) "with_ptC5" else "without_ptC5"
seeds <- as.integer(strsplit(arg(args, "seeds", "2024,1,2,3,4,5"), ",")[[1]])
exclude_clusters <- if (version == "with_ptC5") character(0) else "pt_C5"
ARMS <- list(baseline = COHORT_EXCLUDED_GROUPS, exclusion = c(COHORT_EXCLUDED_GROUPS, "control"))

pbo <- readRDS(arg(args, "pseudobulk"))
peaks_gr <- peak_granges(readRDS(arg(args, "peak_universe")))
variants <- read_egfr_variants(arg(args, "variants"))

pc <- lapply(ARMS, function(ex) population_counts(pbo$counts, pbo$groups, ex, exclude_clusters))
arm_cells <- function(g) {
  nuc <- tapply(g$n_nuclei, g$level1, sum); lib <- tapply(g$sample, g$level1, function(s) length(unique(s)))
  data.frame(population = POPULATIONS, n_nuclei = as.integer(nuc[POPULATIONS]), n_libraries = as.integer(lib[POPULATIONS]), stringsAsFactors = FALSE)
}
cells <- lapply(pc, function(x) arm_cells(x$groups_kept))
n_lib <- sapply(pc, function(x) length(unique(x$groups_kept$sample)))
n_nuc <- sapply(pc, function(x) sum(x$groups_kept$n_nuclei))
removed <- setdiff(unique(pc$baseline$groups_kept$sample), unique(pc$exclusion$groups_kept$sample))
log_msg(sprintf("baseline %d libraries / %d nuclei; without control libraries %d / %d (removed: %s)",
                n_lib[["baseline"]], n_nuc[["baseline"]], n_lib[["exclusion"]], n_nuc[["exclusion"]], paste(sort(removed), collapse = ",")))

run_arm <- function(arm) {
  cpm1 <- (to_cpm(pc[[arm]]$counts) + 1)[names(peaks_gr), POPULATIONS]
  run_gchromvar(cpm1, peaks_gr, variants, "Kidney_eGFR_CARMA", seeds)
}
g <- list(baseline = run_arm("baseline"), exclusion = run_arm("exclusion"))

# ---- per-population comparison ----
call_of <- c(`enriched_FDR<0.05` = "enriched", `depleted_FDR<0.05` = "depleted", ns = "ns")
arm_table <- function(arm) {
  r <- g[[arm]]$result[match(POPULATIONS, g[[arm]]$result$cell_type), ]
  s <- g[[arm]]$by_seed[match(POPULATIONS, g[[arm]]$by_seed$cell_type), ]
  cc <- cells[[arm]]
  out <- data.frame(population = POPULATIONS, n_nuclei = cc$n_nuclei, n_libraries = cc$n_libraries, z = r$z, p = r$p, FDR = r$FDR,
                    call = unname(call_of[r$significance]), rank = r$rank, z_mean_6seeds = s$z_mean_6seeds, z_sd_6seeds = s$z_sd_6seeds,
                    n_seeds_FDR_lt_0.05 = s$n_seeds_FDR_lt_0.05, rank_range_6seeds = s$rank_range, stringsAsFactors = FALSE)
  names(out)[-1] <- paste0(names(out)[-1], "_", arm)
  out
}
res <- merge(arm_table("baseline"), arm_table("exclusion"), by = "population", sort = FALSE)
res <- res[order(-res$z_baseline), ]
res$display_label <- DISPLAY_LABEL[res$population]
res$delta_z <- res$z_exclusion - res$z_baseline
res$delta_rank <- res$rank_exclusion - res$rank_baseline
res$call_change <- ifelse(res$call_baseline == res$call_exclusion, "none", paste(res$call_baseline, "->", res$call_exclusion))
res$sign_reversal_at_FDR <- res$call_baseline != "ns" & res$call_exclusion != "ns" & sign(res$z_baseline) != sign(res$z_exclusion)
res$z_sign_change <- sign(res$z_baseline) != sign(res$z_exclusion)
res$ptc5_version <- version
res$seed_primary <- seeds[1]; res$seeds_summary <- paste(seeds, collapse = ",")
res$n_peaks_with_weight <- g$baseline$result$n_peaks_with_weight[1]; res$weight_sum <- g$baseline$result$weight_sum[1]
res <- res[, c("population", "display_label", setdiff(names(res), c("population", "display_label")))]
rownames(res) <- NULL
write_tsv(res, arg(args, "out_comparison"))

runs <- rbind(transform(g$baseline$all_runs, arm = "baseline"), transform(g$exclusion$all_runs, arm = "exclusion"))
runs$ptc5_version <- version
write_tsv(runs, arg(args, "out_by_seed"))

# ---- summary ----
chg <- res[res$call_change != "none", ]
in_order <- function(calls, k) paste(res$population[calls == k][order(res$rank_exclusion[calls == k])], collapse = ",")
summ <- data.frame(
  ptc5_version = version,
  libraries_baseline = n_lib[["baseline"]], nuclei_baseline = n_nuc[["baseline"]],
  libraries_exclusion = n_lib[["exclusion"]], nuclei_exclusion = n_nuc[["exclusion"]],
  libraries_removed = paste(sort(removed), collapse = ","),
  spearman_rho_z_seed2024 = cor(res$z_baseline, res$z_exclusion, method = "spearman"),
  spearman_rho_z_mean_6seeds = cor(res$z_mean_6seeds_baseline, res$z_mean_6seeds_exclusion, method = "spearman"),
  n_call_changes = nrow(chg),
  call_changes = if (nrow(chg)) paste(sprintf("%s %s (z %.2f -> %.2f; FDR %.3g -> %.3g)", chg$population, chg$call_change, chg$z_baseline,
                                              chg$z_exclusion, chg$FDR_baseline, chg$FDR_exclusion), collapse = "; ") else "none",
  n_sign_reversal_at_FDR = sum(res$sign_reversal_at_FDR), n_z_sign_change = sum(res$z_sign_change),
  order_baseline = paste(res$population[order(res$rank_baseline)], collapse = " > "),
  order_exclusion = paste(res$population[order(res$rank_exclusion)], collapse = " > "),
  enriched_exclusion = in_order(res$call_exclusion, "enriched"), depleted_exclusion = in_order(res$call_exclusion, "depleted"),
  stringsAsFactors = FALSE)
write_tsv(summ, arg(args, "out_summary"))
log_msg(sprintf("Spearman rho (seed 2024) = %.3f; call changes: %s", summ$spearman_rho_z_seed2024, summ$call_changes))

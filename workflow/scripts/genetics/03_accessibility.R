# (C) Sergio Alejandro Gómez
# Aim: Pseudobulk accessibility of the 10 kidney populations in the analysis cohort (37 libraries; BKVAN, KIN and
#   AKI-control libraries excluded).
# Inputs: pseudobulk object (02_pseudobulk.R).
# Outputs: peak universe, raw counts, CPM and CPM + 1 matrices, nuclei metadata (.rds); cohort composition table.
# Usage: Rscript workflow/scripts/genetics/03_accessibility.R --pseudobulk=results/genetics/full_cohort/objects/pseudobulk.rds --include_ptC5=true --out_peak_universe=results/genetics/with_ptC5/objects/peak_universe.rds --out_counts=results/genetics/with_ptC5/objects/accessibility_counts.rds --out_cpm=results/genetics/with_ptC5/objects/accessibility_cpm.rds --out_cpm_plus1=results/genetics/with_ptC5/objects/accessibility_cpm_plus1.rds --out_cell_metadata=results/genetics/with_ptC5/objects/cell_metadata.rds --out_composition=results/genetics/with_ptC5/tables/cohort_composition.tsv

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "accessibility.R"))

args <- parse_args()
version <- if (arg_true(args, "include_ptC5", NULL)) "with_ptC5" else "without_ptC5"
exclude_groups <- strsplit(arg(args, "exclude_disease_groups", paste(COHORT_EXCLUDED_GROUPS, collapse = ",")), ",")[[1]]
exclude_clusters <- if (version == "with_ptC5") character(0) else "pt_C5"
pseudocount <- as.numeric(arg(args, "pseudocount", "1"))

pbo <- readRDS(arg(args, "pseudobulk"))
cells <- pbo$cells
cells$in_cohort <- !(cells$diseaseNames %in% exclude_groups) & !(cells$cluster %in% exclude_clusters)
n_cells <- sum(cells$in_cohort)
log_msg(sprintf("%s cohort: %d libraries, %d nuclei", version, length(unique(cells$sample[cells$in_cohort])), n_cells))
save_rds(cells, arg(args, "out_cell_metadata"))
save_rds(pbo$peaks, arg(args, "out_peak_universe"))

# ---- population x peak matrices ----
pc <- population_counts(pbo$counts, pbo$groups, exclude_groups, exclude_clusters)
cnt <- pc$counts
stopifnot(setequal(colnames(cnt), POPULATIONS), sum(pc$groups_kept$n_nuclei) == n_cells)
cpm <- to_cpm(cnt)
cpm1 <- cpm + pseudocount
log_msg(sprintf("peaks with zero count in every population: %d (kept; the peak set is fixed)", sum(rowSums(cnt) == 0)))
save_rds(cnt, arg(args, "out_counts"))
save_rds(cpm, arg(args, "out_cpm"))          # cell-type specificity (no pseudocount)
save_rds(cpm1, arg(args, "out_cpm_plus1"))   # gchromVAR counts assay

# ---- cohort composition ----
comp <- do.call(rbind, lapply(DISPLAY_ORDER, function(l) {
  idx <- cells$in_cohort & cells$level1 == l
  data.frame(population = l, display_label = DISPLAY_LABEL[[l]], n_nuclei = sum(idx), pct_nuclei = round(100 * sum(idx) / n_cells, 2),
             n_libraries_with_nuclei = length(unique(cells$sample[idx])), stringsAsFactors = FALSE)
}))
comp$ptc5_version <- version
comp$peak_universe <- PEAK_UNIVERSE_LABEL
comp$n_peaks_zero_in_cohort <- sum(rowSums(cnt) == 0)
write_tsv(comp, arg(args, "out_composition"))

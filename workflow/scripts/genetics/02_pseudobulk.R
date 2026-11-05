# (C) Sergio Alejandro Gómez
# Aim: Load the cluster x library pseudobulk of the ArchR PeakMatrix exported by 01_export_peak_matrix.py and attach the
#   peak table, the nuclei and the 10-population annotation of each cluster.
# Inputs: export directory of 01_export_peak_matrix.py, cluster -> population table (resources/genetics/atac_population_labels.csv).
# Outputs: pseudobulk object (.rds): list(counts = dgCMatrix peaks x cluster|library, groups, peaks, cells, provenance).
# Usage: Rscript workflow/scripts/genetics/02_pseudobulk.R --export_dir=results/genetics/full_cohort/peak_matrix_export --population_labels=resources/genetics/atac_population_labels.csv --out_pseudobulk=results/genetics/full_cohort/objects/pseudobulk.rds

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "accessibility.R"))

args <- parse_args()
d <- arg(args, "export_dir")

ann <- read.csv(arg(args, "population_labels"), stringsAsFactors = FALSE, strip.white = TRUE)
population_of <- setNames(ann$level1, ann$cluster)

meta <- jsonlite::fromJSON(file.path(d, "export_summary.json"))
peaks <- parse_archr_peaks(read.delim(file.path(d, "pseudobulk_peaks.tsv"), colClasses = "character")$peak)
groups <- read.delim(file.path(d, "pseudobulk_groups.tsv"), colClasses = c("character", "character", "character", "character", "integer"))
cells <- read.delim(gzfile(file.path(d, "cells.tsv.gz")), colClasses = "character")
np <- nrow(peaks); ng <- nrow(groups)
stopifnot(np == meta$n_peaks, ng == meta$n_groups)

# the binary file is the peaks x groups integer matrix, column-major
con <- file(file.path(d, "pseudobulk_counts_int32_le.bin"), "rb")
v <- readBin(con, what = "integer", n = np * ng, size = 4, endian = "little"); close(con)
stopifnot(length(v) == np * ng)
m <- matrix(v, nrow = np, ncol = ng, dimnames = list(peaks$peak, groups$group)); rm(v)
pb <- list(counts = as(m, "CsparseMatrix"), groups = groups[, c("group", "cluster", "sample", "n_nuclei")]); rm(m); invisible(gc())
stopifnot(isTRUE(all.equal(sum(pb$counts@x), as.numeric(meta$total_count))),
          identical(as.integer(table(factor(paste(cells$cluster, cells$sample, sep = "|"), levels = groups$group))), groups$n_nuclei))

# ---- peaks, nuclei and groups ----
peaks$is_standard_chr <- peaks$chr %in% STANDARD_CHROMOSOMES
peaks$width <- peaks$end - peaks$start + 1L
stopifnot(!anyDuplicated(peaks$peak), identical(rownames(pb$counts), peaks$peak), all(cells$cluster %in% names(population_of)))
library_disease <- unique(cells[, c("sample", "diseaseNames")])
stopifnot(!anyDuplicated(library_disease$sample))   # one disease group per library
cells$level1 <- unname(population_of[cells$cluster])
pb$groups$diseaseNames <- setNames(library_disease$diseaseNames, library_disease$sample)[pb$groups$sample]
pb$groups$level1 <- unname(population_of[pb$groups$cluster])
log_msg(sprintf("%d peaks x %d cluster|library groups; %d nuclei in %d libraries", np, ng, nrow(cells), length(unique(cells$sample))))

provenance <- list(input_format = "h5ad_export", input = normalizePath(d, winslash = "/"),
                   n_nuclei_object = nrow(cells), n_peaks = nrow(peaks), nnz = meta$nnz, total_count = meta$total_count,
                   value_spectrum = paste(names(unlist(meta$value_spectrum)), as.numeric(unlist(meta$value_spectrum)), sep = ":", collapse = ";"),
                   created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), R = R.version.string, Matrix = as.character(packageVersion("Matrix")))
save_rds(list(counts = pb$counts, groups = pb$groups, peaks = peaks, cells = cells, provenance = provenance), arg(args, "out_pseudobulk"))

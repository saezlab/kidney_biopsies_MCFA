# (C) Sergio Alejandro Gómez
# Helpers of the gchromVAR scripts (04, 07, 08, 18): reading the fine-mapped eGFR variants, PIP weights per peak and the
# weighted gchromVAR enrichment over cell populations. The same machinery serves eGFR and CKD.

suppressPackageStartupMessages(library(GenomicRanges))

read_egfr_variants <- function(path) {
  df <- read.delim(path, colClasses = "character", check.names = FALSE, quote = "", comment.char = "")
  for (v in intersect(c("pos_hg38", "pip", "gwas_z", "maf", "gwas_n", "pos_hg19"), names(df))) df[[v]] <- as.numeric(df[[v]])
  for (v in intersect(c("strand_flipped_hg19_to_hg38", "effect_allele_is_ALT"), names(df))) df[[v]] <- as.logical(df[[v]])
  gr <- GRanges(df$chr, IRanges(as.integer(df$pos_hg38), as.integer(df$pos_hg38)))
  mcols(gr) <- df[, setdiff(names(df), c("chr", "pos_hg38")), drop = FALSE]
  stopifnot(all(width(gr) == 1), all(gr$pip >= 0 & gr$pip <= 1))
  gr
}

peak_weights <- function(peaks_gr, variants, pip = variants$pip) {
  ov <- findOverlaps(peaks_gr, variants)
  w <- numeric(length(peaks_gr))
  s <- tapply(pip[subjectHits(ov)], queryHits(ov), sum)
  w[as.integer(names(s))] <- s
  w
}

gchromvar_table <- function(z, dev, cell_types, run_id, fdr = 0.05) {
  p <- 2 * pnorm(-abs(z))
  df <- data.frame(run_id = run_id, cell_type = cell_types, z = as.numeric(z), deviation = as.numeric(dev),
                   p = p, FDR = p.adjust(p, method = "BH"), stringsAsFactors = FALSE)
  df <- df[order(-df$z), ]; df$rank <- seq_len(nrow(df))
  df$significance <- ifelse(df$FDR < fdr & df$z > 0, "enriched_FDR<0.05", ifelse(df$FDR < fdr & df$z < 0, "depleted_FDR<0.05", "ns"))
  rownames(df) <- NULL
  df
}

run_gchromvar <- function(cpm1, peaks_gr, variants, trait, seeds = c(2024L, 1:5), n_background = 50L, fdr = 0.05) {
  suppressPackageStartupMessages({
    library(chromVAR); library(gchromVAR); library(SummarizedExperiment)
    library(BSgenome.Hsapiens.UCSC.hg38); library(BiocParallel)
  })
  register(SerialParam())
  stopifnot(identical(rownames(cpm1), names(peaks_gr)), all(width(variants) == 1))
  w <- matrix(peak_weights(peaks_gr, variants), ncol = 1, dimnames = list(NULL, trait))
  log_msg(sprintf("%s: %d variants, %d peaks carry weight (sum %.3f)", trait, length(variants), sum(w > 0), sum(w)))
  stopifnot(sum(w > 0) > 0)
  SE <- SummarizedExperiment(assays = list(counts = as.matrix(cpm1)), rowRanges = peaks_gr,
                             colData = DataFrame(cell_type = colnames(cpm1), row.names = colnames(cpm1)))
  SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg38)
  SE_w <- SummarizedExperiment(assays = list(weights = w), rowRanges = peaks_gr, colData = DataFrame(trait = trait))
  runs <- lapply(seeds, function(s) {
    set.seed(s)
    bg <- getBackgroundPeaks(SE, niterations = n_background)
    dev <- computeWeightedDeviations(object = SE, weights = SE_w, background_peaks = bg)
    tab <- gchromvar_table(assays(dev)[["z"]][1, ], assays(dev)[["deviations"]][1, ], colnames(cpm1), sprintf("%s_seed%d", trait, s), fdr)
    tab$seed <- s; tab$n_peaks_with_weight <- sum(w > 0); tab$weight_sum <- sum(w)
    log_msg(sprintf("seed %d: %s", s, paste(sprintf("%s=%.2f", tab$cell_type, tab$z), collapse = " ")))
    tab
  })
  all_runs <- do.call(rbind, runs)
  by_seed <- do.call(rbind, lapply(split(all_runs, all_runs$cell_type), function(d) data.frame(
    cell_type = d$cell_type[1], z_seed2024 = d$z[d$seed == seeds[1]], z_mean_6seeds = mean(d$z), z_sd_6seeds = sd(d$z),
    z_min = min(d$z), z_max = max(d$z), n_seeds_FDR_lt_0.05 = sum(d$FDR < fdr),
    rank_range = paste(range(d$rank), collapse = "-"), stringsAsFactors = FALSE)))
  by_seed <- by_seed[order(-by_seed$z_seed2024), ]; rownames(by_seed) <- NULL
  res <- runs[[1]]
  res$seed_mean_z <- by_seed$z_mean_6seeds[match(res$cell_type, by_seed$cell_type)]
  res$seed_sd_z <- by_seed$z_sd_6seeds[match(res$cell_type, by_seed$cell_type)]
  list(result = res, by_seed = by_seed, all_runs = all_runs)
}

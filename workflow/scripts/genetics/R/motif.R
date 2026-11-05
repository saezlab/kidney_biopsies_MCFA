# (C) Sergio Alejandro Gómez
# Helpers of the motif disruption scripts (14_motif_disruption.R, 15_motif_pt_vs_tal.R).

# ---- motifbreakR --------------------------------------------------------------------------------------------

score_chunk <- function(f, motifs, method, bkg, threshold) {
  suppressPackageStartupMessages({ library(motifbreakR); library(MotifDb); library(BSgenome.Hsapiens.UCSC.hg38) })
  s <- motifbreakR::snps.from.file(file = f, search.genome = BSgenome.Hsapiens.UCSC.hg38, format = "bed")
  motifbreakR::motifbreakR(snpList = s, filterp = TRUE, pwmList = motifs, show.neutral = FALSE, method = method,
                           bkg = bkg, threshold = threshold, BPPARAM = BiocParallel::SerialParam())
}

run_motifbreakr <- function(bed_file, chunk_dir, n_workers, seed, method, threshold,
                            bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)) {
  motifs <- query(MotifDb, "hsapiens", "jaspar2022")
  log_msg("MotifDb ", as.character(packageVersion("MotifDb")), ": ", length(motifs), " human JASPAR 2022 matrices; motifbreakR ",
          as.character(packageVersion("motifbreakR")))
  snps <- snps.from.file(file = bed_file, search.genome = BSgenome.Hsapiens.UCSC.hg38, format = "bed")

  bed <- read.table(bed_file, sep = "\t", stringsAsFactors = FALSE)
  n_chunks <- min(nrow(bed), n_workers * 4L)
  chunk_id <- ceiling(seq_len(nrow(bed)) * n_chunks / nrow(bed))
  dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)
  chunk_files <- vapply(seq_len(n_chunks), function(k) {
    f <- file.path(chunk_dir, sprintf("chunk_%02d.bed", k))
    write.table(bed[chunk_id == k, ], f, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    f
  }, "")

  t0 <- Sys.time()
  bpp <- SnowParam(workers = n_workers, type = "SOCK", RNGseed = seed, exportglobals = FALSE)
  res_list <- bplapply(chunk_files, score_chunk, motifs = motifs, method = method, bkg = bkg, threshold = threshold, BPPARAM = bpp)
  bpstop(bpp)
  unlink(chunk_dir, recursive = TRUE)

  res <- do.call(c, unname(res_list))
  attributes(res)$genome.package <- "BSgenome.Hsapiens.UCSC.hg38"
  attributes(res)$motifs <- motifs[mcols(motifs)$providerId %in% unique(res$providerId)]
  log_msg(sprintf("motifbreakR finished in %.1f min: %d variant x matrix rows, %d variants, %d strong",
                  as.numeric(difftime(Sys.time(), t0, units = "mins")), length(res), length(unique(res$SNP_id)), sum(res$effect == "strong")))
  list(res = res, snps = snps, motifs = motifs, n_chunks = n_chunks)
}

summarise_strong_effects <- function(res, snps, cand) {
  df <- as.data.frame(mcols(res)); df$SNP_id <- as.character(df$SNP_id); df$chr <- as.character(seqnames(res)); df$pos <- start(res)
  df$motifPos <- vapply(df$motifPos, function(x) paste(x, collapse = ":"), "")
  for (cc in colnames(df)) if (is.list(df[[cc]])) df[[cc]] <- vapply(df[[cc]], function(x) paste(x, collapse = ";"), "")
  strong <- df[df$effect == "strong", c("chr", "pos", "SNP_id", "rsid", "PIP", "REF", "ALT", "geneSymbol", "providerId", "dataSource", "effect",
                                        "alleleDiff", "alleleEffectSize", "pctRef", "pctAlt", "scoreRef", "scoreAlt", "motifPos", "specific_cell_types")]
  strong <- strong[order(strong$geneSymbol, -abs(strong$alleleDiff)), ]

  tf <- strong |> group_by(geneSymbol) |>
    summarise(n_distinct_disrupted_snps = n_distinct(SNP_id), n_strong_snp_motif_events = n(), n_motif_matrices = n_distinct(providerId),
              motif_matrices = paste(sort(unique(providerId)), collapse = ";"),
              n_snps_stronger_ref = n_distinct(SNP_id[alleleDiff < 0]), n_snps_stronger_alt = n_distinct(SNP_id[alleleDiff > 0]),
              mean_abs_alleleDiff = mean(abs(alleleDiff)), .groups = "drop") |>
    arrange(desc(n_distinct_disrupted_snps), desc(n_strong_snp_motif_events), geneSymbol) |> mutate(rank_by_distinct_snps = row_number()) |>
    arrange(desc(n_strong_snp_motif_events), geneSymbol) |> mutate(rank_by_events = row_number()) |>
    arrange(rank_by_distinct_snps) |> as.data.frame()
  tf$event_inflation_ratio <- round(tf$n_strong_snp_motif_events / tf$n_distinct_disrupted_snps, 3)
  tf$rsids_disrupting <- vapply(tf$geneSymbol, function(g) paste(sort(unique(strong$rsid[strong$geneSymbol == g])), collapse = ";"), "")

  by_variant <- strong |> group_by(SNP_id) |>
    summarise(n_strong_motifs = n(), n_distinct_tfs = n_distinct(geneSymbol), tfs = paste(sort(unique(geneSymbol)), collapse = ";"), .groups = "drop") |>
    arrange(desc(n_distinct_tfs)) |> as.data.frame()
  by_variant$rsid <- snps$rsid[match(by_variant$SNP_id, names(snps))]
  by_variant$PIP <- cand$PIP[match(by_variant$rsid, cand$rsID)]
  by_variant$specific_cell_types <- cand$specific_cell_types[match(by_variant$rsid, cand$rsID)]
  list(strong = strong, tf = tf, by_variant = by_variant)
}

# ---- promoter-proximal / distal peaks -----------------------------------------------------------------------

peak_class <- function(peak_names, gene_annotation, window) {
  pk <- unique(peak_names)
  m <- regmatches(pk, regexec("^(.+)_([0-9]+)_([0-9]+)$", pk))
  stopifnot(all(lengths(m) == 4))
  peaks_gr <- GRanges(vapply(m, `[`, "", 2), IRanges(as.integer(vapply(m, `[`, "", 3)), as.integer(vapply(m, `[`, "", 4))))
  tss <- gene_annotation$transcript_tss[gene_annotation$transcript_tss$gene_biotype == "protein_coding"]
  d2t <- distanceToNearest(peaks_gr, tss, ignore.strand = TRUE)
  dist <- rep(NA_integer_, length(pk)); dist[queryHits(d2t)] <- mcols(d2t)$distance
  win <- GRanges(seqnames(tss), IRanges(pmax(1L, start(tss) - window), start(tss) + window))
  n_win <- countOverlaps(peaks_gr, win, ignore.strand = TRUE)
  out <- data.frame(peak = pk, peak_class = ifelse(n_win > 0, "promoter_proximal", "distal"), distance_to_tss = dist, stringsAsFactors = FALSE)
  out[match(peak_names, out$peak), , drop = FALSE]
}

# ---- TF prevalence in groups of candidates (script 15) --------------------------------------------------------

rank_by_count <- function(count, label) {
  o <- order(-count, label, method = "radix")
  r <- integer(length(count)); r[o] <- seq_along(count); r
}

boundary_tie_n <- function(count, rank, k) { ck <- count[rank == k]; if (length(ck) == 0) NA_integer_ else sum(count == ck) }

stratum_prevalence <- function(long, snp_ids, universe, stratum) {
  sub <- long[SNP_id %in% snp_ids]
  cnt <- sub[, .(n_distinct_disrupted_snps = uniqueN(SNP_id)), by = geneSymbol]
  out <- data.table(stratum = stratum, geneSymbol = universe)
  out <- merge(out, cnt, by = "geneSymbol", all.x = TRUE, sort = FALSE)
  out[is.na(n_distinct_disrupted_snps), n_distinct_disrupted_snps := 0L]
  out[, n_distinct_disrupted_snps := as.integer(n_distinct_disrupted_snps)]
  out[, n_candidates_group := length(snp_ids)]
  out[, fraction_disrupted := n_distinct_disrupted_snps / n_candidates_group]
  out[, rank_by_count := rank_by_count(n_distinct_disrupted_snps, geneSymbol)]
  out[, rank_by_fraction := rank_by_count]   # the same within a group (fixed denominator)
  out[, present_in_stratum := n_distinct_disrupted_snps > 0L]
  setcolorder(out, c("stratum", "geneSymbol", "n_candidates_group", "n_distinct_disrupted_snps", "fraction_disrupted", "rank_by_count", "rank_by_fraction", "present_in_stratum"))
  setorder(out, rank_by_count)
  out
}

similarity_metrics <- function(A, B, comparison, unit, label_A, label_B, n_shared_units, diff_band, similar_rho, partly_similar_rho, top_k) {
  stopifnot(identical(A$geneSymbol, B$geneSymbol))
  fa <- A$fraction_disrupted; fb <- B$fraction_disrupted
  rho_frac <- suppressWarnings(stats::cor(fa, fb, method = "spearman"))
  rho_rank <- suppressWarnings(stats::cor(A$rank_by_count, B$rank_by_count, method = "spearman"))
  r_pearson <- suppressWarnings(stats::cor(fa, fb, method = "pearson"))
  tau <- suppressWarnings(stats::cor(fa, fb, method = "kendall"))
  d <- fa - fb
  band <- if (rho_frac >= similar_rho) "GLOBALLY_SIMILAR" else if (rho_frac >= partly_similar_rho) "PARTIALLY_SIMILAR" else "DIFFERENT"
  summ <- data.table(comparison = comparison, unit = unit, group_A = label_A, group_B = label_B,
                     n_units_A = A$n_candidates_group[1], n_units_B = B$n_candidates_group[1], n_shared_units = n_shared_units,
                     n_tf_universe = nrow(A), n_tf_present_both = sum(A$present_in_stratum & B$present_in_stratum),
                     n_tf_only_A = sum(A$present_in_stratum & !B$present_in_stratum), n_tf_only_B = sum(!A$present_in_stratum & B$present_in_stratum),
                     n_tf_absent_both = sum(!A$present_in_stratum & !B$present_in_stratum),
                     spearman_rho_fraction = rho_frac, spearman_rho_rank = rho_rank, pearson_r_fraction = r_pearson, kendall_tau_fraction = tau,
                     mean_fraction_A = mean(fa), mean_fraction_B = mean(fb), median_fraction_A = median(fa), median_fraction_B = median(fb),
                     mean_diff_A_minus_B = mean(d), median_diff_A_minus_B = median(d), mean_abs_diff = mean(abs(d)), max_abs_diff = max(abs(d)),
                     n_tf_abs_diff_ge_band = sum(abs(d) >= diff_band), n_tf_A_higher_ge_band = sum(d >= diff_band), n_tf_B_higher_ge_band = sum(d <= -diff_band),
                     similarity_band = band, rank_rule = "desc(n_distinct_disrupted_snps), geneSymbol",
                     spearman_note = "Spearman rho on fractions uses average ranks for ties (stats::cor, method = 'spearman'); rank-rank rho uses the deterministic row_number ranks")
  topk <- rbindlist(lapply(top_k, function(k) {
    sa <- A$geneSymbol[A$rank_by_count <= k]; sa <- sa[order(A$rank_by_count[match(sa, A$geneSymbol)])]
    sb <- B$geneSymbol[B$rank_by_count <= k]; sb <- sb[order(B$rank_by_count[match(sb, B$geneSymbol)])]
    ov <- intersect(sa, sb)
    data.table(comparison = comparison, unit = unit, k = k, group_A = label_A, group_B = label_B,
               overlap_n = length(ov), union_n = length(union(sa, sb)), jaccard = length(ov) / length(union(sa, sb)),
               boundary_tie_n_A = boundary_tie_n(A$n_distinct_disrupted_snps, A$rank_by_count, k),
               boundary_tie_n_B = boundary_tie_n(B$n_distinct_disrupted_snps, B$rank_by_count, k),
               top_k_A_in_rank_order = paste(sa, collapse = ";"), top_k_B_in_rank_order = paste(sb, collapse = ";"),
               overlap_labels = paste(sort(ov), collapse = ";"))
  }))
  list(summary = summ, topk = topk)
}

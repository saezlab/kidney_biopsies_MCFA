# (C) Sergio Alejandro Gómez
# Aim: Define the CKD fine-mapping loci from the CKDGen summary statistics and prepare their inputs.
# Inputs: --ckd_gwas  CKDGen CKD summary statistics (European ancestry; Wuttke et al., Nat Genet 2019; GRCh37)
#         --ld_panel  PLINK prefix of the 1000 Genomes Phase 3 European reference panel (503 individuals; GRCh37)
# Outputs: --out_loci       one row per locus: window, lead variant, variant counts, LD-reference coverage
#          --out_variants   harmonised variants of every locus (effect allele = CKDGen Allele1, z = beta / SE)
#          --out_genotypes  list of reference genotype matrices, one per locus (individuals x variants, effect-allele counts)
# Usage: Rscript workflow/scripts/genetics/16_ckd_loci_and_ld.R --ckd_gwas=data/genetics/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt --ld_panel=data/genetics/g1000_eur --out_loci=results/genetics/ckd_finemapping/ckd_loci.tsv --out_variants=results/genetics/ckd_finemapping/ckd_locus_variants.tsv --out_genotypes=results/genetics/ckd_finemapping/ckd_locus_genotypes.rds

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "plink.R"))
suppressPackageStartupMessages(library(GenomicRanges))

args <- parse_args()
gws_p          <- as.numeric(arg(args, "gws_p", "5e-8"))      # genome-wide significance
window_bp      <- as.numeric(arg(args, "window_bp", "1e6"))   # locus = significant variants +/- 1 Mb, overlapping windows merged
min_n_fraction <- as.numeric(arg(args, "min_n_fraction", "0.8"))   # keep variants with N >= 80% of the largest N in the locus
min_maf        <- as.numeric(arg(args, "min_maf", "0.01"))    # common variants: minor allele frequency >= 1%
max_missing    <- as.numeric(arg(args, "max_missing", "0.05"))
panel          <- arg(args, "ld_panel")

# ---- summary statistics -------------------------------------------------------------------------------------------
log_msg("reading ", arg(args, "ckd_gwas"))
gwas <- fread(arg(args, "ckd_gwas"), showProgress = FALSE)
setnames(gwas, c("Chr", "Pos_b37", "RSID", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value", "n_total_sum"),
               c("chr", "pos", "rsid", "A1", "A2", "freq1", "beta", "se", "p", "n"))
gwas[, `:=`(chr = as.character(chr), pos = as.integer(pos), A1 = toupper(A1), A2 = toupper(A2), freq1 = as.numeric(freq1),
            beta = as.numeric(beta), se = as.numeric(se), p = as.numeric(p), n = as.numeric(n))]
log_msg(nrow(gwas), " variants")

# ---- loci ----------------------------------------------------------------------------------------------------------
gws <- gwas[p < gws_p][order(p)]
merged <- reduce(GRanges(gws$chr, IRanges(pmax(1, gws$pos - window_bp), gws$pos + window_bp)))
merged <- merged[order(as.integer(as.character(seqnames(merged))), start(merged))]
loci <- data.frame(chr = as.character(seqnames(merged)), fm_start = start(merged), fm_end = end(merged), stringsAsFactors = FALSE)
loci$locus <- sprintf("CKD_L%02d_chr%s_%d_%d", seq_len(nrow(loci)), loci$chr, loci$fm_start, loci$fm_end)

# lead variant = smallest P in the locus
hits <- findOverlaps(merged, GRanges(gws$chr, IRanges(gws$pos, width = 1)))
lead <- sapply(seq_len(nrow(loci)), function(i) { j <- subjectHits(hits)[queryHits(hits) == i]; j[which.min(gws$p[j])] })
loci$lead_variant <- gws$rsid[lead]; loci$lead_pos <- gws$pos[lead]; loci$lead_p <- gws$p[lead]
loci$n_gws_variants <- tabulate(queryHits(hits), nrow(loci))
log_msg(nrow(gws), " genome-wide significant variants in ", nrow(loci), " loci")

# ---- reference panel -------------------------------------------------------------------------------------------------
fam <- read_fam(panel)
bim <- read_bim(panel)[chr %in% loci$chr]
log_msg("reference panel: ", nrow(fam), " individuals")

harmonise_locus <- function(g, b) {
  acgt <- c("A", "C", "G", "T")
  dup_pos <- duplicated(g$pos) | duplicated(g$pos, fromLast = TRUE)   # any two GWAS rows at one position (incl. indels)
  g[, valid := is.finite(p) & p > 0 & p <= 1 & is.finite(beta) & is.finite(se) & se > 0 & is.finite(n) &
               A1 %in% acgt & A2 %in% acgt & A1 != A2 & !dup_pos]
  g[, enough_n := valid & n >= min_n_fraction * max(n[valid])]
  g[, palindromic := paste0(A1, A2) %in% c("AT", "TA", "CG", "GC")]

  b <- b[!(duplicated(pos) | duplicated(pos, fromLast = TRUE))]   # positions with more than one panel SNV are not used
  j <- match(g$pos, b$pos)
  pa1 <- b$A1[j]; pa2 <- b$A2[j]; comp <- function(a) chartr("ACGT", "TGCA", a)
  same    <- (g$A1 == pa1 & g$A2 == pa2) | (comp(g$A1) == pa1 & comp(g$A2) == pa2)
  swapped <- (g$A1 == pa2 & g$A2 == pa1) | (comp(g$A1) == pa2 & comp(g$A2) == pa1)
  g[, aligned := !is.na(j) & (same | swapped)]
  g[, effect_is_panel_a2 := aligned & swapped]                    # the panel dosage counts A2: flip when the effect allele is A1
  g[, panel_idx := b$idx[j]]

  # LD-reference coverage: fraction of the common, N-eligible GWAS SNVs that are present in the panel
  maf <- pmin(g$freq1, 1 - g$freq1); common <- is.finite(maf) & maf >= min_maf
  coverage <- round(sum(g$enough_n & common & g$aligned) / max(1, sum(g$enough_n & common)), 4)
  v <- g[enough_n & aligned & !palindromic]
  setorder(v, pos, A1, A2)
  list(variants = v, coverage = coverage)
}

# ---- per locus: variants and reference genotypes -----------------------------------------------------------------
variants <- list(); genotypes <- list(); locus_stats <- list()
for (i in seq_len(nrow(loci))) {
  id <- loci$locus[i]
  g <- gwas[chr == loci$chr[i] & pos >= loci$fm_start[i] & pos <= loci$fm_end[i]]
  b <- bim[chr == loci$chr[i] & pos >= loci$fm_start[i] & pos <= loci$fm_end[i]]
  h <- harmonise_locus(g, b); v <- h$variants

  # reference genotypes, oriented so that they count the GWAS effect allele
  G <- read_bed_dosage(panel, v$panel_idx, nrow(fam))
  flip <- which(!v$effect_is_panel_a2)
  G[, flip] <- 2L - G[, flip, drop = FALSE]
  af <- colMeans(G, na.rm = TRUE) / 2
  keep <- is.finite(pmin(af, 1 - af)) & pmin(af, 1 - af) >= min_maf & colMeans(is.na(G)) <= max_missing
  v <- v[keep]; G <- G[, keep, drop = FALSE]; af <- af[keep]

  v[, `:=`(locus = id, variant_id = paste(chr, pos, A1, A2, sep = ":"), z = beta / se, is_gws = p < gws_p, reference_af = af)]
  # a GWAS effect-allele frequency that matches the panel frequency of the other allele points to an allele-coding error
  v[, af_matches_other_allele := abs(freq1 - af) > 0.5 & abs(freq1 - (1 - af)) < 0.1]
  colnames(G) <- v$variant_id; rownames(G) <- fam$IID

  ok <- is.finite(v$freq1) & is.finite(af)
  locus_stats[[id]] <- data.frame(n_variants = nrow(v), n_gws_retained = sum(v$is_gws), ld_coverage = h$coverage,
                                  af_r = if (sum(ok) > 2) cor(v$freq1[ok], af[ok]) else NA, n_locus = median(v$n))
  variants[[id]] <- v[, .(locus, variant_id, chr, pos, rsid, A1, A2, freq1, beta, se, p, n, z, is_gws, reference_af, af_matches_other_allele)]
  genotypes[[id]] <- G
  log_msg(sprintf("%s: %d GWAS variants -> %d retained (%d genome-wide significant); LD-reference coverage %.4f",
                  id, nrow(g), nrow(v), sum(v$is_gws), h$coverage))
}

loci <- cbind(loci[, c("locus", "chr", "fm_start", "fm_end", "lead_variant", "lead_pos", "lead_p", "n_gws_variants")],
              do.call(rbind, locus_stats))
write_tsv(loci, arg(args, "out_loci"))
write_tsv(rbindlist(variants), arg(args, "out_variants"))
save_rds(genotypes, arg(args, "out_genotypes"))

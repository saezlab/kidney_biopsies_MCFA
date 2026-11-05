# (C) Sergio Alejandro Gómez
# Aim:     Predict which transcription-factor (TF) motifs the eGFR candidate variants disrupt.
# Inputs:  tables/candidate_variants.tsv and objects/candidate_variants.bed (candidate variants); the Ensembl 108 gene
#          annotation object (transcript TSSs), used to label each candidate's peak promoter-proximal or distal.
# Outputs: tables/motif_disruption_events.tsv (strong variant x matrix effects), tables/motif_disruption_by_tf.tsv,
#          tables/motif_disruption_by_variant.tsv, tables/motif_disruption_top30_tfs.tsv, tables/motifbreakr_settings.tsv,
#          tables/candidate_peak_class_counts.tsv, objects/motifbreakr_results.rds (full motifbreakR result).
# Usage: Rscript workflow/scripts/genetics/14_motif_disruption.R --candidates=results/genetics/with_ptC5/tables/candidate_variants.tsv --bed=results/genetics/with_ptC5/objects/candidate_variants.bed --gene_annotation=results/genetics/qtl_resources/gene_annotation_hg38_ensembl108.rds --out_events=results/genetics/with_ptC5/tables/motif_disruption_events.tsv --out_by_tf=results/genetics/with_ptC5/tables/motif_disruption_by_tf.tsv --out_by_variant=results/genetics/with_ptC5/tables/motif_disruption_by_variant.tsv --out_top30=results/genetics/with_ptC5/tables/motif_disruption_top30_tfs.tsv --out_settings=results/genetics/with_ptC5/tables/motifbreakr_settings.tsv --out_peak_class=results/genetics/with_ptC5/tables/candidate_peak_class_counts.tsv --out_results=results/genetics/with_ptC5/objects/motifbreakr_results.rds

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "motif.R"))
suppressPackageStartupMessages({
  library(GenomicRanges); library(dplyr)
  library(motifbreakR); library(MotifDb); library(BSgenome.Hsapiens.UCSC.hg38); library(BiocParallel)
})

args <- parse_args()
n_workers <- as.integer(arg(args, "n_workers", "6"))        # SOCK workers; the BED is cut into n_workers x 4 chunks
seed      <- as.integer(arg(args, "seed", "2024"))          # RNG seed of the worker cluster
method    <- arg(args, "method", "ic")                      # motifbreakR scoring: information content
threshold <- as.numeric(arg(args, "threshold", "0.85"))     # p-value threshold (filterp = TRUE)
top_n     <- as.integer(arg(args, "top_n_tfs", "30"))
window    <- as.integer(arg(args, "promoter_window_bp", "2000"))
bkg       <- c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)      # uniform background

# ---- load ------------------------------------------------------------------------------------------------------
cand <- read_tsv(arg(args, "candidates"), colClasses = c(REF = "character", ALT = "character"))
bed_file <- arg(args, "bed")
gene_annotation <- readRDS(arg(args, "gene_annotation"))
out_results <- arg(args, "out_results")
log_msg(nrow(cand), " candidate variants")

# ---- motifbreakR -------------------------------------------------------------------------------------------------
mb <- run_motifbreakr(bed_file, file.path(dirname(out_results), "motifbreakr_chunks"), n_workers, seed, method, threshold, bkg)
res <- mb$res; snps <- mb$snps
stopifnot(length(snps) == nrow(cand))   # every REF allele matches hg38

key_snp <- paste(as.character(seqnames(snps)), start(snps)); key_cand <- paste(cand$chr, cand$hg38_position)
snps$rsid <- cand$rsID[match(key_snp, key_cand)]
stopifnot(!anyNA(snps$rsid), all(as.character(snps$REF) == cand$REF[match(key_snp, key_cand)]))
res$rsid <- snps$rsid[match(as.character(res$SNP_id), names(snps))]
res$PIP <- cand$PIP[match(res$rsid, cand$rsID)]
res$specific_cell_types <- cand$specific_cell_types[match(res$rsid, cand$rsID)]
metadata(res)$setting <- list(filterp = TRUE, threshold = threshold, method = method, bkg = bkg, motif_query = c("hsapiens", "jaspar2022"),
                              n_motifs = length(mb$motifs), motifbreakR = as.character(packageVersion("motifbreakR")),
                              MotifDb = as.character(packageVersion("MotifDb")))
save_rds(res, out_results)

sm <- summarise_strong_effects(res, snps, cand)

# ---- promoter-proximal / distal label of each candidate's peak (descriptive) ----------------------------------------
cand$peak <- vapply(strsplit(cand$overlapping_peak, ";", fixed = TRUE), `[`, "", 1)
pc <- peak_class(cand$peak, gene_annotation, window)
cand$peak_class <- pc$peak_class
cand$distance_to_nearest_protein_coding_tss <- pc$distance_to_tss
stopifnot(!anyNA(cand$peak_class))

add_peak_class <- function(tab) {
  tab$overlapping_peak <- cand$peak[match(tab$rsid, cand$rsID)]
  i <- match(tab$rsid, cand$rsID)
  tab$peak_class <- cand$peak_class[i]
  tab$distance_to_nearest_protein_coding_tss <- cand$distance_to_nearest_protein_coding_tss[i]
  tab
}
strong <- add_peak_class(sm$strong)
by_variant <- add_peak_class(sm$by_variant)
tf <- sm$tf
for (cl in c("promoter_proximal", "distal"))
  tf[[paste0("n_distinct_disrupted_snps_", cl)]] <- vapply(tf$geneSymbol, function(g)
    length(unique(strong$SNP_id[strong$geneSymbol == g & strong$peak_class == cl])), 0L)
top30 <- head(tf[order(tf$rank_by_distinct_snps), ], top_n)[, c("rank_by_distinct_snps", "geneSymbol", "n_distinct_disrupted_snps",
                                                                 "n_motif_matrices", "n_strong_snp_motif_events", "motif_matrices")]

# candidates and peaks by class, for all candidates and for each attributed population (non-exclusive)
attribution <- strsplit(cand$specific_cell_types, ";", fixed = TRUE)
groups <- c(list(ALL_CANDIDATES = rep(TRUE, nrow(cand))),
            lapply(setNames(nm = sort(unique(unlist(attribution)))), function(p) vapply(attribution, function(x) p %in% x, TRUE)))
peak_counts <- do.call(rbind, lapply(names(groups), function(g) do.call(rbind, lapply(c("promoter_proximal", "distal"), function(cl) {
  k <- groups[[g]] & cand$peak_class == cl
  data.frame(group = g, peak_class = cl, n_candidates = sum(k), n_peaks = length(unique(cand$peak[k])), stringsAsFactors = FALSE)
}))))

settings <- data.frame(
  item = c("motif_source", "query", "n_motifs", "MotifDb_version", "motifbreakR_version", "BSgenome_version", "method", "background",
           "filterp", "show.neutral", "threshold", "effect_class_used", "n_workers", "n_chunks", "rng_seed", "n_candidates_input",
           "n_snps_loaded", "n_rows_total", "n_rows_strong", "n_snps_with_strong", "n_tfs_with_strong"),
  value = c("MotifDb (Bioconductor) human JASPAR 2022 subset", "hsapiens,jaspar2022", length(mb$motifs), as.character(packageVersion("MotifDb")),
            as.character(packageVersion("motifbreakR")), as.character(packageVersion("BSgenome.Hsapiens.UCSC.hg38")), method,
            paste(names(bkg), bkg, sep = "=", collapse = ","), "TRUE", "FALSE", threshold, "strong", n_workers, mb$n_chunks, seed,
            nrow(cand), length(snps), length(res), nrow(strong), length(unique(strong$SNP_id)), length(unique(strong$geneSymbol))),
  stringsAsFactors = FALSE)

# ---- write -------------------------------------------------------------------------------------------------------
write_tsv(strong, arg(args, "out_events"))
write_tsv(tf, arg(args, "out_by_tf"))
write_tsv(by_variant, arg(args, "out_by_variant"))
write_tsv(top30, arg(args, "out_top30"))
write_tsv(settings, arg(args, "out_settings"))
write_tsv(peak_counts, arg(args, "out_peak_class"))
log_msg(length(unique(strong$SNP_id)), " of ", nrow(cand), " candidates carry a strong effect; top TF ", tf$geneSymbol[1],
        " (", tf$n_distinct_disrupted_snps[1], " candidates)")

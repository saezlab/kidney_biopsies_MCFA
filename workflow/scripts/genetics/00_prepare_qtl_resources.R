# (C) Sergio Alejandro Gómez
# Aim: Prepare the public kidney gene-regulation resources used by scripts 09-13.
# Inputs: Open4Gene_62K_Cells_Peak_to_Gene_Significant_Associations.txt.gz, the GTEx significant pairs as gzip TSV
#   (00_gtex_parquet_to_tsv.py), Kidney_Cortex.v10.eGenes.txt.gz, Homo_sapiens.GRCh38.108.gtf,
#   bASA_RASQUAL_82_Kidney_Significant_Associations.txt.gz, UCSC hg19ToHg38.over.chain.
# Outputs: open4gene_hg38_unique.rds, gtex_kidney_unique_evariants_hg38.rds, gtex_kidney_significant_eqtl_pairs.rds,
#   gene_annotation_hg38_ensembl108.rds, basa_unique_caqtl_variants_hg38.rds, basa_caqtl_peaks_hg38.rds,
#   basa_caqtl_associations_hg38.tsv.gz
# Usage: Rscript workflow/scripts/genetics/00_prepare_qtl_resources.R --open4gene=data/genetics/Open4Gene_62K_Cells_Peak_to_Gene_Significant_Associations.txt.gz --gtex_pairs=results/genetics/qtl_resources/gtex_kidney_signif_pairs.tsv.gz --gtex_egenes=data/genetics/Kidney_Cortex.v10.eGenes.txt.gz --gtf=data/genetics/Homo_sapiens.GRCh38.108.gtf --basa=data/genetics/bASA_RASQUAL_82_Kidney_Significant_Associations.txt.gz --chain=data/genetics/hg19ToHg38.over.chain --out_open4gene=results/genetics/qtl_resources/open4gene_hg38_unique.rds --out_gtex_evariants=results/genetics/qtl_resources/gtex_kidney_unique_evariants_hg38.rds --out_gtex_pairs=results/genetics/qtl_resources/gtex_kidney_significant_eqtl_pairs.rds --out_gene_annotation=results/genetics/qtl_resources/gene_annotation_hg38_ensembl108.rds --out_caqtl_variants=results/genetics/qtl_resources/basa_unique_caqtl_variants_hg38.rds --out_caqtl_peaks=results/genetics/qtl_resources/basa_caqtl_peaks_hg38.rds --out_caqtl_associations=results/genetics/qtl_resources/basa_caqtl_associations_hg38.tsv.gz

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
suppressPackageStartupMessages({ library(GenomicRanges); library(rtracklayer); library(BSgenome.Hsapiens.UCSC.hg38) })

args <- parse_args()
min_mapped_fraction <- as.numeric(arg(args, "liftover_min_mapped_fraction", "0.95"))
max_span_ratio <- as.numeric(arg(args, "liftover_max_span_ratio", "2"))
bs <- BSgenome.Hsapiens.UCSC.hg38
std <- c(paste0("chr", 1:22), "chrX", "chrY")
chain <- import.chain(arg(args, "chain"))

# ---- liftOver helpers ----
# Whole-interval liftOver. df: id, chr, start, end (1-based closed). One row per interval with its category.
liftover_intervals <- function(df) {
  gr <- GRanges(df$chr, IRanges(df$start, df$end)); names(gr) <- df$id
  lifted <- liftOver(gr, chain); np <- lengths(lifted)
  ul <- unlist(lifted, use.names = FALSE)
  pieces <- data.table(i = rep(seq_along(lifted), np), chr = as.character(seqnames(ul)), strand = as.character(strand(ul)), start = start(ul), end = end(ul), w = width(ul))
  agg <- pieces[, .(n_pieces = .N, n_chr = uniqueN(chr), n_strand = uniqueN(strand), chr_hg38 = if (uniqueN(chr) == 1) chr[1] else paste(unique(chr), collapse = ";"),
                    start_hg38 = if (uniqueN(chr) == 1) min(start) else NA_integer_, end_hg38 = if (uniqueN(chr) == 1) max(end) else NA_integer_, mapped_bases = sum(w)), by = i]
  per <- data.table(i = seq_along(lifted), n_pieces = 0L, n_chr = 0L, n_strand = 0L, chr_hg38 = NA_character_, start_hg38 = NA_integer_, end_hg38 = NA_integer_, mapped_bases = 0L)
  per[agg$i, `:=`(n_pieces = agg$n_pieces, n_chr = agg$n_chr, n_strand = agg$n_strand, chr_hg38 = agg$chr_hg38, start_hg38 = agg$start_hg38, end_hg38 = agg$end_hg38, mapped_bases = agg$mapped_bases)]
  per[, i := NULL]
  per <- cbind(as.data.table(df), width_hg19 = width(gr), per)
  per[, width_hg38 := ifelse(is.na(start_hg38), NA_integer_, end_hg38 - start_hg38 + 1L)]
  per[, mapped_fraction := mapped_bases / width_hg19]
  per[, span_ratio := width_hg38 / width_hg19]
  per[, chromosome_changed := !is.na(chr_hg38) & n_chr == 1 & chr_hg38 != chr]
  per[, category := fifelse(n_pieces == 0, "UNMAPPED",
                    fifelse(n_chr > 1 | n_strand > 1, "MULTI_CHROMOSOME_OR_STRAND",
                    fifelse(chromosome_changed, "CHROMOSOME_CHANGED",
                    fifelse(!(chr_hg38 %in% std), "NON_STANDARD_TARGET",
                    fifelse(span_ratio > max_span_ratio, "SPAN_TOO_LARGE",
                    fifelse(mapped_fraction < min_mapped_fraction, "PARTIAL_MAPPING", "UNIQUE"))))))]
  per[, id_hg38 := ifelse(category == "UNIQUE", paste(chr_hg38, start_hg38, end_hg38, sep = "-"), NA_character_)]
  per
}

# 1-bp variant liftOver. df: chr, pos, ref, alt. Alleles are complemented when the base lifts to the minus strand;
# a variant is accepted when the hg38 base equals the oriented REF, or the oriented ALT (then REF and ALT are swapped).
complement_allele <- function(a) chartr("ACGTacgt", "TGCAtgca", a)
liftover_variants_1bp <- function(df) {
  lo <- liftOver(GRanges(df$chr, IRanges(df$pos, df$pos)), chain); n <- lengths(lo)
  out <- data.frame(df, n_pieces = n, chr_hg38 = NA_character_, pos_hg38 = NA_integer_, strand_hg38 = NA_character_, stringsAsFactors = FALSE)
  one <- which(n == 1)
  u <- unlist(lo[one])
  out$chr_hg38[one] <- as.character(seqnames(u)); out$pos_hg38[one] <- start(u); out$strand_hg38[one] <- as.character(strand(u))
  out$width_hg38 <- NA_integer_; out$width_hg38[one] <- width(u)
  out$strand_flipped <- !is.na(out$strand_hg38) & out$strand_hg38 == "-"
  out$ref_hg38_oriented <- ifelse(out$strand_flipped, complement_allele(out$ref), out$ref)
  out$alt_hg38_oriented <- ifelse(out$strand_flipped, complement_allele(out$alt), out$alt)
  out$genome_base_hg38 <- NA_character_
  ok <- !is.na(out$chr_hg38) & out$chr_hg38 %in% std & out$width_hg38 == 1
  if (any(ok)) {
    sl <- seqlengths(bs)[out$chr_hg38[ok]]
    inb <- out$pos_hg38[ok] >= 1 & out$pos_hg38[ok] <= sl
    idx <- which(ok)[inb]
    out$genome_base_hg38[idx] <- as.character(Biostrings::getSeq(bs, GRanges(out$chr_hg38[idx], IRanges(out$pos_hg38[idx], out$pos_hg38[idx]))))
  }
  out$category <- ifelse(n == 0, "UNMAPPED",
                  ifelse(n > 1, "MULTI_MAPPED",
                  ifelse(!(out$chr_hg38 %in% std), "NON_STANDARD_TARGET",
                  ifelse(out$chr_hg38 != out$chr, "CHROMOSOME_CHANGED",
                  ifelse(is.na(out$width_hg38) | out$width_hg38 != 1, "ABNORMAL_WIDTH",
                  ifelse(is.na(out$genome_base_hg38), "OUT_OF_RANGE",
                  ifelse(out$genome_base_hg38 == out$ref_hg38_oriented, "UNIQUE_REF_MATCH",
                  ifelse(out$genome_base_hg38 == out$alt_hg38_oriented, "UNIQUE_SWAPPED_ALT_MATCHES_GENOME", "UNIQUE_ALLELE_MISMATCH"))))))))
  sw <- out$category == "UNIQUE_SWAPPED_ALT_MATCHES_GENOME"
  out$REF_hg38 <- ifelse(out$category == "UNIQUE_REF_MATCH", out$ref_hg38_oriented, ifelse(sw, out$alt_hg38_oriented, NA_character_))
  out$ALT_hg38 <- ifelse(out$category == "UNIQUE_REF_MATCH", out$alt_hg38_oriented, ifelse(sw, out$ref_hg38_oriented, NA_character_))
  out$ref_alt_swapped_vs_source <- sw
  out$accepted <- out$category %in% c("UNIQUE_REF_MATCH", "UNIQUE_SWAPPED_ALT_MATCHES_GENOME")
  out
}

# ---- 1. Open4Gene peak-to-gene links ----
og <- fread(arg(args, "open4gene"), sep = "\t", showProgress = FALSE); og[, row_id := seq_len(.N)]
pu <- unique(og[, .(peak, Chr, PeakStart, PeakEnd)])
stopifnot(!any(duplicated(pu$peak)))
per <- liftover_intervals(data.frame(id = pu$peak, chr = pu$Chr, start = pu$PeakStart, end = pu$PeakEnd, stringsAsFactors = FALSE))
log_msg("Open4Gene peak liftOver: ", paste(per[, .N, by = category][order(-N), paste(category, N, sep = "=")], collapse = "; "))
setnames(og, c("Chr", "PeakStart", "PeakEnd", "seqnames", "start", "end", "PeakCenter"),
         c("chr_hg19", "start_hg19", "end_hg19", "gene_chr_hg19_ensembl_style", "gene_start_hg19", "gene_end_hg19", "peak_center_hg19"))
setnames(og, "peak", "peak_hg19")
lk <- merge(og, per[, .(peak_hg19 = id, peak_hg38 = id_hg38, chr_hg38, start_hg38, end_hg38, width_hg19, width_hg38, n_pieces, mapped_fraction, span_ratio, category)],
            by = "peak_hg19", all.x = TRUE)
setorder(lk, row_id)

crosswalk <- data.frame(
  open4gene_celltype = c("PT", "DCT", "CNT", "PC", "IC_A", "IC_B", "Podo", "PEC", "Endo_Peritubular", "Endo_GC", "Endo_Lympho", "Mac", "T_Lympho", "B_Lympho", "Fibro",
                         "Injured_PT", "Mes", "LOH", "Macula_Densa", "AllCell"),
  level1 = c("PT", "DCT", "CNT_PC", "CNT_PC", "IC", "IC", "POD", "PEC", "EC", "EC", "EC", "Immune", "Immune", "Immune", "Stromal", "PT", "Stromal", "TAL", "TAL", NA),
  tier = c(rep("UNAMBIGUOUS", 15), rep("LINEAGE", 4), "NOT_MAPPED"), stringsAsFactors = FALSE)
cw <- match(lk$Celltype, crosswalk$open4gene_celltype)
lk[, level1_unambiguous := ifelse(crosswalk$tier[cw] == "UNAMBIGUOUS", crosswalk$level1[cw], NA_character_)]
lk[, level1_lineage := ifelse(crosswalk$tier[cw] %in% c("UNAMBIGUOUS", "LINEAGE"), crosswalk$level1[cw], NA_character_)]
lk[, crosswalk_tier := crosswalk$tier[cw]]
lk[, min_fdr := pmin(zinb.res.zero.p.All.fdr, zinb.res.count.p.All.fdr, na.rm = TRUE)]
stopifnot(!anyNA(lk$crosswalk_tier))
save_rds(as.data.frame(lk[category == "UNIQUE"]), arg(args, "out_open4gene"))

# ---- 2. GTEx v10 kidney cortex: significant pairs and unique eVariants ----
pairs <- fread(arg(args, "gtex_pairs"), sep = "\t", showProgress = FALSE)
stopifnot(!anyNA(pairs), !any(duplicated(pairs[, .(gene_id, variant_id)])))
vp <- tstrsplit(pairs$variant_id, "_", fixed = TRUE)      # variant_id = chr_pos_ref_alt_b38
pairs[, `:=`(chr = vp[[1]], pos = as.integer(vp[[2]]), ref = vp[[3]], alt = vp[[4]], build = vp[[5]])]
stopifnot(all(pairs$build == "b38"), all(pairs$chr %in% std), all(grepl("^[ACGT]+$", pairs$ref) & grepl("^[ACGT]+$", pairs$alt)))
pairs[, variant_class := fifelse(nchar(ref) == 1 & nchar(alt) == 1, "SNV", fifelse(nchar(ref) == nchar(alt), "MNV", "INDEL"))]
eg <- fread(arg(args, "gtex_egenes"), sep = "\t", showProgress = FALSE)
egm <- eg[, .(gene_id, gene_name, biotype, gene_chr, gene_start, gene_end, gene_strand = strand, gene_qval = qval, gene_pval_nominal_threshold = pval_nominal_threshold,
              lead_variant_id = variant_id, lead_rsid = rs_id_dbSNP155_GRCh38p13, lead_slope = slope, lead_pval_nominal = pval_nominal)]
pairs <- merge(pairs, egm, by = "gene_id", all.x = TRUE, sort = FALSE)
stopifnot(all(pairs$gene_qval <= 0.05), all(pairs$pval_nominal <= pairs$pval_nominal_threshold))
pairs[, alt_allele_direction := fifelse(slope > 0, "ALT_increases_expression", fifelse(slope < 0, "ALT_decreases_expression", "no_effect"))]
pairs[, is_lead_variant_for_gene := variant_id == lead_variant_id]

uv <- unique(pairs[, .(variant_id, chr, pos, ref, alt, variant_class)])
gseq <- as.character(getSeq(bs, GRanges(uv$chr, IRanges(uv$pos, uv$pos + nchar(uv$ref) - 1L))))
uv[, ref_matches_hg38 := gseq == ref]
stopifnot(all(uv$ref_matches_hg38))
agg <- pairs[order(gene_name), .(n_egenes = uniqueN(gene_id), egenes = paste(unique(gene_name), collapse = ";"), egene_ids = paste(unique(gene_id), collapse = ";"),
                                 n_protein_coding_egenes = uniqueN(gene_id[biotype == "protein_coding"]),
                                 min_pval_nominal = min(pval_nominal), max_abs_slope = max(abs(slope)), min_gene_pval_beta = min(pval_beta), min_gene_qval = min(gene_qval),
                                 n_genes_as_lead_variant = sum(is_lead_variant_for_gene), lead_rsid = paste(unique(lead_rsid[is_lead_variant_for_gene]), collapse = ";"),
                                 min_abs_tss_distance = min(abs(tss_distance)), af_gtex = af[1]), by = variant_id]
uv <- merge(uv, agg, by = "variant_id", sort = FALSE)
uv[lead_rsid == "", lead_rsid := NA_character_]
setorder(uv, chr, pos)
ev_gr <- GRanges(uv$chr, IRanges(uv$pos, uv$pos), seqinfo = Seqinfo(std, seqlengths = seqlengths(bs)[std], genome = "hg38"))
mcols(ev_gr) <- as.data.frame(uv)
names(ev_gr) <- uv$variant_id
save_rds(ev_gr, arg(args, "out_gtex_evariants"))
setcolorder(pairs, c("variant_id", "chr", "pos", "ref", "alt", "variant_class", "gene_id", "gene_name", "biotype"))
save_rds(as.data.frame(pairs), arg(args, "out_gtex_pairs"))
log_msg(sprintf("GTEx: %d significant pairs, %d unique eVariants, %d eGenes", nrow(pairs), nrow(uv), uniqueN(pairs$gene_id)))

# ---- 3. Ensembl 108 genes and transcript TSSs (UCSC chromosome names; scaffolds dropped) ----
gtf <- arg(args, "gtf")
con <- file(gtf, "r"); hdr <- character(0)
repeat { l <- readLines(con, n = 1); if (!length(l) || !startsWith(l, "#")) break; hdr <- c(hdr, l) }
close(con)
stopifnot(any(grepl("genome-version GRCh38", hdr)))
d <- fread(gtf, sep = "\t", header = FALSE, quote = "", skip = length(hdr), select = c(1, 3, 4, 5, 7, 9),
           col.names = c("seqname", "feature", "start", "end", "strand", "attr"), showProgress = FALSE)
n_lines <- nrow(d)
d <- d[feature %in% c("gene", "transcript")]
getattr <- function(x, key) { m <- regmatches(x, regexpr(paste0(key, ' "[^"]*"'), x)); out <- rep(NA_character_, length(x)); hit <- grepl(paste0(key, ' "'), x); out[hit] <- sub(paste0('^', key, ' "'), "", sub('"$', "", m)); out }
d[, gene_id := getattr(attr, "gene_id")]; d[, gene_name := getattr(attr, "gene_name")]; d[, gene_biotype := getattr(attr, "gene_biotype")]
d[, transcript_id := getattr(attr, "transcript_id")]; d[, transcript_biotype := getattr(attr, "transcript_biotype")]
d[, ensembl_canonical := grepl('tag "Ensembl_canonical"', attr, fixed = TRUE)]
d[, attr := NULL]
d[, chr := ifelse(seqname %in% c(as.character(1:22), "X", "Y"), paste0("chr", seqname), ifelse(seqname == "MT", "chrM", NA_character_))]
d <- d[!is.na(chr)]
d[is.na(gene_name) | gene_name == "", gene_name := gene_id]
d[, tss := ifelse(strand == "+", start, end)]
genes <- d[feature == "gene"]; tx <- d[feature == "transcript"]
stopifnot(!any(duplicated(genes$gene_id)), !any(duplicated(tx$transcript_id)), all(tx$gene_id %in% genes$gene_id))
genes_gr <- GRanges(genes$chr, IRanges(genes$start, genes$end), strand = genes$strand, gene_id = genes$gene_id, gene_name = genes$gene_name, gene_biotype = genes$gene_biotype, tss = genes$tss)
names(genes_gr) <- genes$gene_id
tss_gr <- GRanges(tx$chr, IRanges(tx$tss, tx$tss), strand = tx$strand, transcript_id = tx$transcript_id, gene_id = tx$gene_id, gene_name = tx$gene_name, gene_biotype = tx$gene_biotype,
                  transcript_biotype = tx$transcript_biotype, ensembl_canonical = tx$ensembl_canonical)

save_rds(list(genes = genes_gr, transcript_tss = tss_gr, source = gtf, source_sha256 = digest::digest(file = gtf, algo = "sha256"), header = hdr,
              build = "GRCh38 (Ensembl 108; chromosome names converted to UCSC style)", n_gtf_lines = n_lines), arg(args, "out_gene_annotation"))
log_msg(sprintf("Ensembl 108: %d genes (%d protein-coding), %d transcript TSSs", nrow(genes), sum(genes$gene_biotype == "protein_coding"), nrow(tx)))

# ---- 4. bASA/RASQUAL kidney caQTLs: variants, target peaks and associations ----
ca <- fread(arg(args, "basa"), sep = "\t", showProgress = FALSE); ca[, row_id := .I]
# 4a. variants: 1-bp liftOver; two hg19 variants landing on one hg38 base are both excluded
vsum <- ca[, .(chr_hg19 = Chromosome[1], pos_hg19 = SNP_position[1], ref_hg19 = Ref_allele[1], alt_hg19 = Alt_allele[1],
               allele_frequency_as_provided = Allele_frequency[1], maf_as_provided = MAF[1],
               n_target_peaks = uniqueN(Feature_ID), min_pvalue = min(Pvalue), min_log10_bh_q = min(Log10_Benjamini_Hochberg_Qvalue), max_abs_effect_size_pi_minus_0.5 = max(abs(Effect_size - 0.5)),
               n_associations_snp_inside_target_peak = sum(Distance == 0), validated_any_matrixeqtl_or_stratas = any(Validated != "None"), n_nonconverged_associations = sum(Convergence_status_0success != 0)), by = rs_ID]
lv <- as.data.table(liftover_variants_1bp(data.frame(rs_ID = vsum$rs_ID, chr = vsum$chr_hg19, pos = vsum$pos_hg19, ref = vsum$ref_hg19, alt = vsum$alt_hg19, stringsAsFactors = FALSE)))
lv[, hg38_key := ifelse(accepted, paste(chr_hg38, pos_hg38), NA_character_)]
coll <- lv[accepted == TRUE, .N, by = hg38_key][N > 1, hg38_key]
lv[hg38_key %in% coll, `:=`(category = "HG38_POSITION_COLLISION", accepted = FALSE)]
log_msg("caQTL variant liftOver: ", paste(lv[, .N, by = category][order(-N), paste(category, N, sep = "=")], collapse = "; "))
acc <- lv[accepted == TRUE]
stopifnot(all(acc$REF_hg38 == acc$genome_base_hg38))
var_tab <- merge(acc[, .(rs_ID, chr_hg38, pos_hg38, REF_hg38, ALT_hg38, ref_alt_swapped_vs_source, strand_flipped, chr_hg19 = chr, pos_hg19 = pos, ref_hg19 = ref, alt_hg19 = alt, liftover_category = category)],
                 vsum[, -c("chr_hg19", "pos_hg19", "ref_hg19", "alt_hg19")], by = "rs_ID")
setorder(var_tab, chr_hg38, pos_hg38)
vgr <- GRanges(var_tab$chr_hg38, IRanges(var_tab$pos_hg38, var_tab$pos_hg38)); names(vgr) <- var_tab$rs_ID
mcols(vgr) <- as.data.frame(var_tab[, -c("chr_hg38", "pos_hg38")])
save_rds(vgr, arg(args, "out_caqtl_variants"))

# 4b. target peaks: whole-interval liftOver
pw <- unique(ca[, .(Feature_ID, Peak_chr, Peak_start, Peak_end)])
stopifnot(!any(duplicated(pw$Feature_ID)))
psum <- ca[, .(n_caqtl_snps_hg19 = uniqueN(rs_ID), n_caqtl_snps_inside_peak = uniqueN(rs_ID[Distance == 0]), min_pvalue_peak = min(Pvalue)), by = Feature_ID]
per <- liftover_intervals(data.frame(id = pw$Feature_ID, chr = pw$Peak_chr, start = pw$Peak_start, end = pw$Peak_end, stringsAsFactors = FALSE))
setnames(per, c("id", "chr", "start", "end"), c("Feature_ID", "chr_hg19", "start_hg19", "end_hg19"))
per <- merge(per, psum, by = "Feature_ID", sort = FALSE)
log_msg("caQTL target peak liftOver: ", paste(per[, .N, by = category][order(-N), paste(category, N, sep = "=")], collapse = "; "))
uniq <- per[category == "UNIQUE"]
pgr <- GRanges(uniq$chr_hg38, IRanges(uniq$start_hg38, uniq$end_hg38)); names(pgr) <- uniq$Feature_ID
mcols(pgr) <- as.data.frame(uniq[, .(Feature_ID, peak_hg38 = id_hg38, chr_hg19, start_hg19, end_hg19, width_hg19, width_hg38, n_pieces, mapped_fraction, span_ratio, n_caqtl_snps_hg19, n_caqtl_snps_inside_peak, min_pvalue_peak)])
save_rds(pgr, arg(args, "out_caqtl_peaks"))

# 4c. every association row in hg38, with flags telling whether its variant and its peak could be lifted
assoc <- merge(ca, var_tab[, .(rs_ID, chr_hg38, pos_hg38, REF_hg38, ALT_hg38, ref_alt_swapped_vs_source, strand_flipped, variant_liftover_category = liftover_category)], by = "rs_ID", all.x = TRUE, sort = FALSE)
assoc <- merge(assoc, per[, .(Feature_ID, peak_hg38 = id_hg38, peak_chr_hg38 = chr_hg38, peak_start_hg38 = start_hg38, peak_end_hg38 = end_hg38, peak_liftover_category = category)], by = "Feature_ID", all.x = TRUE, sort = FALSE)
setorder(assoc, row_id)
assoc[is.na(variant_liftover_category), variant_liftover_category := lv$category[match(rs_ID, lv$rs_ID)]]
assoc[, variant_usable_hg38 := !is.na(pos_hg38)]; assoc[, peak_usable_hg38 := !is.na(peak_hg38)]
assoc[, usable_hg38 := variant_usable_hg38 & peak_usable_hg38]
assoc[, hg38_status := fifelse(usable_hg38, "BOTH_LIFTED", fifelse(!variant_usable_hg38 & !peak_usable_hg38, "VARIANT_AND_PEAK_EXCLUDED", fifelse(!variant_usable_hg38, "VARIANT_EXCLUDED", "PEAK_EXCLUDED")))]
assoc[usable_hg38 == TRUE, snp_inside_target_peak_hg38 := pos_hg38 >= peak_start_hg38 & pos_hg38 <= peak_end_hg38]
stopifnot(nrow(assoc) == nrow(ca))
setnames(assoc, c("Chromosome", "SNP_position", "Ref_allele", "Alt_allele", "Peak_chr", "Peak_start", "Peak_end"),
         c("snp_chr_hg19", "snp_pos_hg19", "ref_hg19", "alt_hg19", "peak_chr_hg19", "peak_start_hg19", "peak_end_hg19"))
out_assoc <- arg(args, "out_caqtl_associations")
dir.create(dirname(out_assoc), recursive = TRUE, showWarnings = FALSE)
con <- gzfile(out_assoc, "w"); write.table(assoc, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"); close(con)
log_msg(sprintf("wrote %s (%d rows; %d usable in hg38)", out_assoc, nrow(assoc), sum(assoc$usable_hg38)))

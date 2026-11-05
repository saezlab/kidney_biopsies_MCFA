# (C) Sergio Alejandro Gómez
# Aim: Prepare the fine-mapped eGFR variants for the enrichment analyses.
# Inputs: CARMA table (semicolon-separated, decimal commas), UCSC hg19ToHg38 chain, BSgenome.Hsapiens.UCSC.hg38.
# Outputs: variant table (one row per variant, hg38 1-bp position, oriented alleles).
# Usage: Rscript workflow/scripts/genetics/00_prepare_egfr_variants.R --carma=data/genetics/finemapping_carma.csv --chain=data/genetics/hg19ToHg38.over.chain --out_variants=results/genetics/egfr_finemapping/egfr_finemapped_variants_hg38.tsv

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
suppressPackageStartupMessages({ library(GenomicRanges); library(rtracklayer); library(BSgenome.Hsapiens.UCSC.hg38) })

args <- parse_args()
raw <- read.csv(arg(args, "carma"), sep = ";", stringsAsFactors = FALSE, check.names = FALSE)
colnames(raw) <- c("locus", "chr_hg19", "rsid", "pos_hg19", "allele1_effect", "allele2_other", "maf", "gwas_z", "gwas_n", "pip", "credible_set")
for (v in c("maf", "gwas_z", "pip")) raw[[v]] <- as.numeric(gsub(",", ".", raw[[v]], fixed = TRUE))   # decimal commas
stopifnot(!anyDuplicated(raw$rsid), all(raw$allele1_effect %in% c("A", "C", "G", "T")), all(raw$allele2_other %in% c("A", "C", "G", "T")),
          all(raw$pip >= 0 & raw$pip <= 1))
log_msg(sprintf("%d variants, %d with PIP > 0.1", nrow(raw), sum(raw$pip > 0.1)))
chain <- import.chain(arg(args, "chain"))

# 1-bp point liftOver (a 2-bp window would shift the base through minus-strand chain segments)
lo <- liftOver(GRanges(paste0("chr", raw$chr_hg19), IRanges(raw$pos_hg19, raw$pos_hg19), rsid = raw$rsid), chain)
stopifnot(all(lengths(lo) == 1))
# The hg38 positions distributed with the fine-mapping were lifted as 2-bp windows ending at the variant; their end is kept
# as a column, to show where it differs by one base from the point liftOver.
lo2 <- liftOver(GRanges(paste0("chr", raw$chr_hg19), IRanges(raw$pos_hg19 - 1, raw$pos_hg19)), chain)
stopifnot(all(lengths(lo2) == 1))
pos_doc <- end(unlist(lo2))
gr38 <- unlist(lo)
stopifnot(identical(gr38$rsid, raw$rsid))
chr <- as.character(seqnames(gr38)); pos <- start(gr38)

# REF/ALT against the hg38 base, tried in this order: allele2 is the reference; allele1 is (swapped); the complement of
# allele2 or allele1 is (variant lifted through a minus-strand chain segment)
base <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, GRanges(chr, IRanges(pos, pos))))
comp <- c(A = "T", C = "G", G = "C", T = "A")
allele1 <- raw$allele1_effect; allele2 <- raw$allele2_other
ref_check <- ifelse(base == allele2, "REF=allele2_matches_genome",
             ifelse(base == allele1, "REF=allele1_matches_genome_SWAPPED",
             ifelse(base == comp[allele2], "REF=revcomp(allele2)_matches_genome_STRAND_FLIPPED",
             ifelse(base == comp[allele1], "REF=revcomp(allele1)_matches_genome_STRAND_FLIPPED_SWAPPED", "NEITHER_ALLELE_MATCHES_GENOME"))))
ref <- ifelse(ref_check == "REF=allele2_matches_genome", allele2, ifelse(ref_check == "REF=allele1_matches_genome_SWAPPED", allele1,
       ifelse(ref_check == "REF=revcomp(allele2)_matches_genome_STRAND_FLIPPED", unname(comp[allele2]),
       ifelse(ref_check == "REF=revcomp(allele1)_matches_genome_STRAND_FLIPPED_SWAPPED", unname(comp[allele1]), NA))))
alt <- ifelse(ref_check == "REF=allele2_matches_genome", allele1, ifelse(ref_check == "REF=allele1_matches_genome_SWAPPED", allele2,
       ifelse(ref_check == "REF=revcomp(allele2)_matches_genome_STRAND_FLIPPED", unname(comp[allele1]),
       ifelse(ref_check == "REF=revcomp(allele1)_matches_genome_STRAND_FLIPPED_SWAPPED", unname(comp[allele2]), NA))))
flipped <- grepl("STRAND_FLIPPED", ref_check)
effect_hg38 <- ifelse(flipped, unname(comp[allele1]), allele1)
log_msg("REF orientation: ", paste(names(table(ref_check)), table(ref_check), sep = "=", collapse = "; "))

out <- data.frame(chr = chr, pos_hg38 = pos, rsid = raw$rsid, locus = raw$locus, credible_set = raw$credible_set, pip = raw$pip,
                  allele1_effect = allele1, allele2_other = allele2, REF = ref, ALT = alt, maf = raw$maf, gwas_z = raw$gwas_z, gwas_n = raw$gwas_n,
                  chr_hg19 = paste0("chr", raw$chr_hg19), pos_hg19 = raw$pos_hg19, pos_hg38_documented_end = pos_doc,
                  position_source = ifelse(pos == pos_doc, "1bp_liftover=documented_end", "1bp_liftover(documented_end_off_by_one)"),
                  ref_check = ref_check, hg38_reference_base = base, strand_flipped_hg19_to_hg38 = flipped,
                  effect_allele_hg38_oriented = effect_hg38, effect_allele_is_ALT = effect_hg38 == alt, stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(paste(out$chr, out$pos_hg38)), !anyNA(out$REF))
write_tsv(out, arg(args, "out_variants"))

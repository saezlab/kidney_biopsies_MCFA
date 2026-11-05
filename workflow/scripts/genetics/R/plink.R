# (C) Sergio Alejandro Gómez
# Minimal reader for PLINK 1 binary genotype files (.bed / .bim / .fam), used for the 1000 Genomes reference panel.

suppressPackageStartupMessages(library(data.table))

read_fam <- function(prefix) {
  fam <- read.table(paste0(prefix, ".fam"), header = FALSE, stringsAsFactors = FALSE,
                    colClasses = c("character", "character", "character", "character", "integer", "character"))
  names(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")
  fam
}

read_bim <- function(prefix) {
  bim <- fread(paste0(prefix, ".bim"), header = FALSE, showProgress = FALSE,
               colClasses = list(character = c(1, 2, 5, 6), integer = c(3, 4)))
  setnames(bim, c("chr", "snp", "cm", "pos", "A1", "A2"))
  bim[, idx := .I]   # row of the variant in the .bim = its block in the .bed
  bim
}

# Genotypes of the variants at .bim rows `idx`, as copies of the A2 allele (NA if missing); individuals x variants.
read_bed_dosage <- function(prefix, idx, n_samples) {
  bed <- paste0(prefix, ".bed")
  stopifnot(identical(readBin(bed, "raw", 3), as.raw(c(0x6c, 0x1b, 0x01))))
  bytes_per_variant <- ceiling(n_samples / 4)
  code <- c(0L, NA, 1L, 2L)                                                     # two-bit codes 00, 01, 10, 11
  lut <- sapply(0:3, function(k) code[bitwAnd(bitwShiftR(0:255, 2 * k), 3L) + 1L])   # byte value -> its four genotypes

  # the variants of one locus lie in one stretch of the file: read it at once, then decode the requested variants
  first <- min(idx)
  con <- file(bed, "rb"); on.exit(close(con))
  seek(con, 3 + bytes_per_variant * (first - 1))
  block <- readBin(con, "raw", bytes_per_variant * (max(idx) - first + 1))
  bytes <- matrix(as.integer(block), nrow = bytes_per_variant)[, idx - first + 1, drop = FALSE]

  out <- matrix(NA_integer_, nrow = n_samples, ncol = length(idx))
  for (j in seq_along(idx)) out[, j] <- as.vector(t(lut[bytes[, j] + 1L, ]))[seq_len(n_samples)]
  out
}

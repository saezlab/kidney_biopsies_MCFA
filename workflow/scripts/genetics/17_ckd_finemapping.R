# (C) Sergio Alejandro Gómez
# Aim: Fine-map the CKD loci with SuSiE-RSS.
# Inputs: --loci, --variants, --genotypes  outputs of 16_ckd_loci_and_ld.R
#         --chain       UCSC hg19ToHg38.over.chain
#         --susier_lib  optional: R library holding the development version of susieR (0.16.6, commit ef213fee) when it
#                       is not in the default library (workflow/envs/genetics_postinstall.R installs it there)
# Outputs: --out_variants      GRanges (hg38) of the variants of the included loci, with PIPs and credible sets
#          --out_locus_qc      one row per locus: QC class, LD-reference diagnostics and fine-mapping summary
#          --out_pip_by_locus  PIP summary of every included locus
# Usage: Rscript workflow/scripts/genetics/17_ckd_finemapping.R --loci=results/genetics/ckd_finemapping/ckd_loci.tsv --variants=results/genetics/ckd_finemapping/ckd_locus_variants.tsv --genotypes=results/genetics/ckd_finemapping/ckd_locus_genotypes.rds --chain=data/genetics/hg19ToHg38.over.chain --threads=4 --out_variants=results/genetics/ckd_finemapping/ckd_finemapped_variants_hg38.rds --out_locus_qc=results/genetics/ckd_finemapping/ckd_locus_qc.tsv --out_pip_by_locus=results/genetics/ckd_finemapping/ckd_pip_by_locus.tsv

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(script_dir, "R", "utils.R"))
args <- parse_args()
susier_lib <- arg(args, "susier_lib", "")
if (nzchar(susier_lib)) .libPaths(c(susier_lib, .libPaths()))
suppressPackageStartupMessages({ library(GenomicRanges); library(rtracklayer); library(BSgenome.Hsapiens.UCSC.hg38) })
stopifnot(packageVersion("susieR") == "0.16.6")
log_msg("susieR ", packageVersion("susieR"), " from ", dirname(find.package("susieR")))

par <- list(
  seed = as.integer(arg(args, "seed", "2024")),
  L = as.integer(arg(args, "L", "10")),                       # at most 10 causal signals per locus
  coverage = as.numeric(arg(args, "coverage", "0.95")),       # 95% credible sets
  min_abs_corr = as.numeric(arg(args, "min_abs_corr", "0.5")),  # credible sets need purity (min |r| between members) >= 0.5
  max_iter = as.integer(arg(args, "max_iter", "200")),
  kriging_logLR = as.numeric(arg(args, "kriging_logLR", "2")),
  kriging_abs_z = as.numeric(arg(args, "kriging_abs_z", "2")),
  min_variants = as.integer(arg(args, "min_variants", "50")),
  min_ld_coverage = as.numeric(arg(args, "min_ld_coverage", "0.8")),
  max_r_over_B = as.numeric(arg(args, "max_r_over_B", "0.2")))
threads <- as.integer(arg(args, "threads", "1"))

loci <- read_tsv(arg(args, "loci"))
variants <- fread(arg(args, "variants"), showProgress = FALSE, colClasses = list(character = c("A1", "A2", "rsid", "variant_id")))
genotypes <- readRDS(arg(args, "genotypes"))

# ---- one locus ------------------------------------------------------------------------------------------------------
# Quality control. A locus is excluded (BLOCKED_BY_QC) when
#   - fewer than 50 variants remain, the panel covers < 80% of its common variants, or no genome-wide-significant
#     variant remains (then it is not fitted);
#   - SuSiE-RSS did not converge;
#   - r_over_B > 0.20: the effective rank of the reference LD matrix (r) is too large relative to the 503 reference
#     individuals (B) for the panel to be adequate;
#   - R_sensitivity_flag: the evidence for a credible set depends strongly on the finite-reference correction
#     (Bayes-factor attenuation >= 30);
#   - a genome-wide-significant variant had to be removed as an allele-coding error (see below).
# A locus is kept with a warning (PASS_WITH_WARNING) when
#   - WARNING_KRIGING_UNEXPLAINED: kriging_rss finds a variant whose z-score disagrees with the z-scores of the variants
#     in LD with it (logLR > 2 and |z| > 2) and there is no sign of an allele-coding error; the variant is kept;
#   - WARNING_R_RELIABILITY: R_reliability_flag, the package's warning of a possible mismatch between the summary
#     statistics and the reference LD (raised when the residual-artifact score Q_art exceeds 0.10 or by
#     R_sensitivity_flag), while no exclusion criterion of the fit applies;
#   - WARNING_SER_M2_DISCORDANCE: the top variant of the single-effect model (SER, Maller et al. 2012) lies outside
#     every credible set;
#   - WARNING_OTHER: any other warning from susie_rss.
# Otherwise the locus passes (PASS). PASS and PASS_WITH_WARNING loci enter the CKD enrichment.
# Columns prefixed M2 describe this SuSiE-RSS model (finite reference + LD mismatch).
fit_locus <- function(info, v, X, par) {
  input_failures <- function(v) c(
    "final variants after harmonisation >= 50" = nrow(v) < par$min_variants,
    "LD-reference representation among N80-eligible common SNVs >= 0.80" = info$ld_coverage < par$min_ld_coverage,
    "at least one retained original genome-wide-significant variant" = sum(v$is_gws) == 0)
  n_locus <- median(v$n)   # sample size given to SuSiE-RSS: median per-variant N
  qc <- list(fit_attempted = FALSE, hard_fail_reason = "", warning_types = "", n_variants = nrow(v), n_GWS_retained = sum(v$is_gws),
             n_locus = n_locus, M2_package_warnings = NA_character_)
  fail <- input_failures(v)
  if (any(fail)) {
    qc$hard_fail_reason <- paste(names(fail)[fail], collapse = " | ")
    return(list(qc = qc, final_class = "BLOCKED_BY_QC", pips = NULL))
  }

  # LD matrix of the reference panel; its eigen-decomposition is shared by the susieR diagnostics
  storage.mode(X) <- "double"
  R <- suppressWarnings(cor(X, use = "pairwise.complete.obs")); dimnames(R) <- list(v$variant_id, v$variant_id)
  eig <- eigen(R, symmetric = TRUE); attr(R, "eigen") <- eig
  stopifnot(isSymmetric(unclass(R), tol = 1e-10), all(is.finite(R)), max(abs(diag(R) - 1)) < 1e-8,
            min(eig$values) > -1e-6 * max(eig$values), all(is.finite(v$z)), !anyDuplicated(v$pos))

  # LD-mismatch diagnostics: s (inconsistency between z-scores and reference LD) and kriging outliers
  s <- sapply(c("null-mle", "null-partialmle", "null-pseudomle"),
              function(m) tryCatch(susieR::estimate_s_rss(v$z, R, n = n_locus, method = m), error = function(e) NA_real_))
  stopifnot(!is.na(s[["null-mle"]]))
  kr <- susieR::kriging_rss(v$z, R, n = n_locus, s = s[["null-mle"]])$conditional_dist
  outlier <- kr$logLR > par$kriging_logLR & abs(kr$z) > par$kriging_abs_z
  coding_error <- outlier & v$af_matches_other_allele %in% TRUE
  coding_error_gws <- any(v$is_gws[coding_error])
  # outliers explained by an allele-coding error are removed, never flipped
  v <- v[!coding_error]; X <- X[, !coding_error, drop = FALSE]

  # SuSiE-RSS with the reference genotypes: R_finite models the finite size of the reference panel,
  # R_mismatch = "eb" estimates the mismatch between the summary statistics and the reference LD
  fit_warnings <- character(0)
  set.seed(par$seed)
  fit <- withCallingHandlers(
    susieR::susie_rss(z = v$z, X = X, n = n_locus, L = min(par$L, nrow(v)), coverage = par$coverage, min_abs_corr = par$min_abs_corr,
                      max_iter = par$max_iter, estimate_residual_variance = FALSE, refine = FALSE, R_finite = TRUE, R_mismatch = "eb"),
    warning = function(w) { fit_warnings <<- c(fit_warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  d <- fit$R_finite_diagnostics
  stopifnot(!is.null(d), length(fit$pip) == nrow(v), all(fit$pip >= 0 & fit$pip <= 1))
  first <- function(x, cast = as.numeric) if (is.null(x)) NA else cast(x[1])

  cs <- fit$sets$cs; if (is.null(cs)) cs <- list()
  cs_id <- rep(NA_character_, nrow(v)); for (k in seq_along(cs)) cs_id[cs[[k]]] <- names(cs)[k]
  in_cs <- unique(unlist(cs))

  # single-effect model of the package vs the credible sets: top variant and overlap of its 95% set with the credible sets
  alpha <- d$ser_model$alpha; alpha <- if (is.matrix(alpha)) as.numeric(alpha[1, ]) else as.numeric(alpha)
  stopifnot(length(alpha) == nrow(v), !anyNA(alpha))
  ord <- order(-alpha); k95 <- which(cumsum(alpha[ord]) >= 0.95)[1]; if (is.na(k95)) k95 <- length(alpha)
  ser95 <- ord[seq_len(k95)]; ser_top <- which.max(alpha)
  ser_top_in_cs <- if (length(cs)) ser_top %in% in_cs else NA

  r_over_B <- first(d$r_over_B); sensitivity <- first(d$R_sensitivity_flag, as.logical); reliability <- first(d$R_reliability_flag, as.logical)
  fit_fail <- c(
    "M2 converged" = !isTRUE(fit$converged),
    "r_over_B <= 0.20" = !isTRUE(r_over_B <= par$max_r_over_B),
    "R_sensitivity_flag == FALSE" = !identical(sensitivity, FALSE),
    "no genome-wide-significant variant removed as an allele-coding error" = coding_error_gws)
  fail <- c(input_failures(v), fit_fail)
  other_warnings <- fit_warnings[!grepl("Possible summary statistics and R reference mismatch detected", fit_warnings, fixed = TRUE)]
  warn <- c(WARNING_KRIGING_UNEXPLAINED = any(outlier & !coding_error),
            WARNING_R_RELIABILITY = isTRUE(reliability) && !any(fit_fail),
            WARNING_SER_M2_DISCORDANCE = isFALSE(ser_top_in_cs),
            WARNING_OTHER = length(other_warnings) > 0)
  final_class <- if (any(fail)) "BLOCKED_BY_QC" else if (any(warn)) "PASS_WITH_WARNING" else "PASS"

  qc <- list(fit_attempted = TRUE, hard_fail_reason = paste(names(fail)[fail], collapse = " | "), warning_types = paste(names(warn)[warn], collapse = ";"),
             n_variants = nrow(v), n_GWS_retained = sum(v$is_gws), n_locus = n_locus,
             estimate_s_rss = s[["null-mle"]], estimate_s_rss_partialmle = s[["null-partialmle"]], estimate_s_rss_pseudomle = s[["null-pseudomle"]],
             kriging_candidates = sum(outlier), kriging_unexplained = sum(outlier & !coding_error), demonstrated_coding_errors = sum(coding_error),
             r_over_B = r_over_B, Q_art = first(d$Q_art), artifact_flag = first(d$artifact_flag, as.logical), lambda_bias = first(d$lambda_bias),
             B_corrected = first(d$B_corrected), R_sensitivity_flag = sensitivity, R_reliability_flag = reliability, effective_rank = first(d$effective_rank),
             M2_converged = isTRUE(fit$converged), M2_niter = fit$niter, max_PIP = max(fit$pip), M2_top_variant = v$rsid[which.max(fit$pip)],
             n_CS = length(cs), minimum_CS_purity = if (length(cs)) min(fit$sets$purity$min.abs.corr) else NA,
             SER_top_variant = v$rsid[ser_top], SER_top_in_M2_CS = ser_top_in_cs,
             SER_jaccard = length(intersect(ser95, in_cs)) / length(union(ser95, in_cs)),
             M2_package_warnings = paste(unique(fit_warnings), collapse = " / "))
  pips <- data.frame(locus = info$locus, chr = v$chr, pos_b37 = v$pos, rsID = v$rsid, effect_allele = v$A1, other_allele = v$A2,
                     beta = v$beta, z = v$z, p = v$p, is_gws = v$is_gws, PIP = as.numeric(fit$pip), CS_id = cs_id, stringsAsFactors = FALSE)
  list(qc = qc, final_class = final_class, pips = pips)
}

# ---- fine-map every locus (one locus per worker when --threads > 1) ----------------------------------------------------
jobs <- lapply(seq_len(nrow(loci)), function(i) list(info = loci[i, ], v = variants[locus == loci$locus[i]], X = genotypes[[loci$locus[i]]]))
run_job <- function(job) {
  t0 <- Sys.time()
  res <- fit_locus(job$info, job$v, job$X, par)
  log_msg(sprintf("%s: %s (%.1f min)", job$info$locus, res$final_class, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  res
}
largest_first <- order(-loci$n_variants)   # the largest loci take longest: start them first
if (threads > 1) {
  cl <- parallel::makeCluster(threads, outfile = "")
  parallel::clusterExport(cl, c("susier_lib", "fit_locus", "par", "log_msg"))
  parallel::clusterEvalQ(cl, { if (nzchar(susier_lib)) .libPaths(c(susier_lib, .libPaths())); suppressPackageStartupMessages(library(data.table)) })
  results <- parallel::clusterApplyLB(cl, jobs[largest_first], run_job)
  parallel::stopCluster(cl)
} else {
  results <- lapply(jobs[largest_first], run_job)
}
results <- results[order(largest_first)]   # back to locus order
names(results) <- loci$locus

# ---- locus QC table ---------------------------------------------------------------------------------------------------
by_lead_p <- order(loci$lead_p, loci$chr, loci$fm_start)
pilot <- rep("NOT_PILOT", nrow(loci))
pilot[by_lead_p[c(1, ceiling(nrow(loci) / 2), nrow(loci))]] <- c("PILOT_STRONG", "PILOT_TYPICAL", "PILOT_WEAK")

qc_cols <- c("estimate_s_rss", "estimate_s_rss_partialmle", "estimate_s_rss_pseudomle", "kriging_candidates", "kriging_unexplained",
             "demonstrated_coding_errors", "r_over_B", "Q_art", "artifact_flag", "lambda_bias", "B_corrected", "R_sensitivity_flag",
             "R_reliability_flag", "effective_rank", "M2_converged", "M2_niter", "max_PIP", "M2_top_variant", "n_CS", "minimum_CS_purity",
             "SER_top_variant", "SER_top_in_M2_CS", "SER_jaccard")
locus_qc <- do.call(rbind, lapply(seq_len(nrow(loci)), function(i) {
  r <- results[[i]]; q <- r$qc
  q[setdiff(qc_cols, names(q))] <- NA
  data.frame(locus = loci$locus[i], chr = loci$chr[i], fm_start = loci$fm_start[i], fm_end = loci$fm_end[i], lead_variant = loci$lead_variant[i],
             lead_p = loci$lead_p[i], pilot_status = pilot[i], fit_attempted = q$fit_attempted, final_class = r$final_class,
             hard_fail_reason = q$hard_fail_reason, warning_count = if (nzchar(q$warning_types)) length(strsplit(q$warning_types, ";")[[1]]) else 0L,
             warning_types = q$warning_types, n_variants = q$n_variants, n_GWS_retained = q$n_GWS_retained, n_GWS_original = loci$n_gws_variants[i],
             LD_coverage = loci$ld_coverage[i], AF_r = loci$af_r[i], n_locus = q$n_locus, q[qc_cols], M2_package_warnings = q$M2_package_warnings,
             stringsAsFactors = FALSE)
}))
write_tsv(locus_qc, arg(args, "out_locus_qc"))
log_msg(paste(sprintf("%s %d", names(table(locus_qc$final_class)), table(locus_qc$final_class)), collapse = ", "))

# ---- PIPs of the included loci ----------------------------------------------------------------------------------------
included <- locus_qc$locus[locus_qc$final_class %in% c("PASS", "PASS_WITH_WARNING")]
pips <- rbindlist(lapply(included, function(l) cbind(results[[l]]$pips, locus_class = results[[l]]$final_class,
                                                      warning_types = results[[l]]$qc$warning_types)))
setorder(pips, chr, pos_b37, rsID)
pip_by_locus <- pips[, .(locus_class = locus_class[1], warning_types = warning_types[1], n_variants = .N, sum_PIP = sum(PIP), max_PIP = max(PIP),
                         top_variant = rsID[which.max(PIP)], n_PIP_gt_0.10 = sum(PIP > 0.10), n_PIP_gt_0.50 = sum(PIP > 0.50),
                         n_in_CS = sum(!is.na(CS_id)), n_CS = uniqueN(na.omit(CS_id))), by = locus]
setorder(pip_by_locus, locus)
write_tsv(pip_by_locus, arg(args, "out_pip_by_locus"))

# ---- lift to hg38 -----------------------------------------------------------------------------------------------------
lift_to_hg38 <- function(chr, pos, other, effect, chain, genome) {
  hits <- liftOver(GRanges(chr, IRanges(pos, width = 1)), chain)
  one <- lengths(hits) == 1
  u <- unlist(hits[one])
  chr38 <- rep(NA_character_, length(pos)); pos38 <- rep(NA_integer_, length(pos)); width38 <- pos38; flipped <- rep(FALSE, length(pos))
  chr38[one] <- as.character(seqnames(u)); pos38[one] <- start(u); width38[one] <- width(u); flipped[one] <- as.character(strand(u)) == "-"
  ok <- one & chr38 == chr & width38 == 1
  ok[ok] <- pos38[ok] >= 1 & pos38[ok] <= seqlengths(genome)[chr38[ok]]
  base <- rep(NA_character_, length(pos))
  base[ok] <- as.character(getSeq(genome, GRanges(chr38[ok], IRanges(pos38[ok], width = 1))))
  comp <- function(a) chartr("ACGTacgt", "TGCAtgca", a)
  other38 <- ifelse(flipped, comp(other), other); effect38 <- ifelse(flipped, comp(effect), effect)
  ref_is_other <- ok & base == other38; ref_is_effect <- ok & !ref_is_other & base == effect38
  data.frame(chr = chr38, pos = pos38, accepted = ref_is_other | ref_is_effect,
             REF = ifelse(ref_is_other, other38, effect38), ALT = ifelse(ref_is_other, effect38, other38), effect_allele = effect38,
             stringsAsFactors = FALSE)
}
lift <- lift_to_hg38(paste0("chr", pips$chr), pips$pos_b37, pips$other_allele, pips$effect_allele, import.chain(arg(args, "chain")),
                     BSgenome.Hsapiens.UCSC.hg38)
log_msg(sum(lift$accepted), " of ", nrow(pips), " variants lifted to hg38 (sum PIP ", round(sum(pips$PIP[lift$accepted]), 3), " of ", round(sum(pips$PIP), 3), ")")
a <- pips[lift$accepted]; l <- lift[lift$accepted, ]
gr <- GRanges(l$chr, IRanges(l$pos, l$pos), rsid = a$rsID, locus = a$locus, locus_class = a$locus_class, warning_types = a$warning_types,
              pip = a$PIP, credible_set_id = a$CS_id, is_gws = a$is_gws, REF = l$REF, ALT = l$ALT, effect_allele_hg38 = l$effect_allele,
              beta = a$beta, z = a$z, p = a$p, chr_hg19 = paste0("chr", a$chr), pos_hg19 = a$pos_b37)
save_rds(gr, arg(args, "out_variants"))

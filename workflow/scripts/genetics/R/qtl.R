# (C) Sergio Alejandro Gómez
# Helpers of the kidney QTL scripts (09-13): population names, ArchR peak names, the protein-coding TSS set
# and the peak-level logistic enrichment model shared by the eQTL and caQTL analyses.

suppressPackageStartupMessages({ library(GenomicRanges); library(data.table) })

POPULATIONS <- c("POD", "EC", "Immune", "IC", "TAL", "CNT_PC", "DCT", "PT", "PEC", "Stromal")
DISPLAY_ORDER <- c("PT", "PEC", "TAL", "DCT", "CNT_PC", "IC", "POD", "EC", "Stromal", "Immune")

ENRICHMENT_FORMULA <- "OUTCOME ~ specific + splines::ns(log1p(mean_accessibility), df = 4) + splines::ns(gc_content, df = 4) + factor(chr) + splines::ns(log1p(distance_to_tss), df = 4)"

parse_archr_peaks <- function(peak_names) {
  m <- regmatches(peak_names, regexec("^(.+)_([0-9]+)_([0-9]+)$", peak_names))
  stopifnot(all(lengths(m) == 4))
  data.frame(chr = vapply(m, `[`, "", 2), start = as.integer(vapply(m, `[`, "", 3)),
             end = as.integer(vapply(m, `[`, "", 4)), peak = peak_names, stringsAsFactors = FALSE)
}

protein_coding_tss <- function(gann, biotype = "protein_coding") gann$transcript_tss[gann$transcript_tss$gene_biotype == biotype]

cluster_robust_se <- function(fit, cluster) {
  X <- model.matrix(fit); U <- X * residuals(fit, type = "response")
  Uc <- rowsum(U, cluster); G <- nrow(Uc); n <- nrow(X); k <- ncol(X)
  bread <- vcov(fit); V <- (G / (G - 1)) * ((n - 1) / (n - k)) * bread %*% crossprod(Uc) %*% bread
  sqrt(diag(V))
}

fit_peak_enrichment <- function(uni, outcome, cluster_bin_bp = 1e6, fdr = 0.05, populations = POPULATIONS) {
  form <- as.formula(sub("OUTCOME", outcome, ENRICHMENT_FORMULA, fixed = TRUE))
  bins <- paste(uni$chr, floor(uni$start / cluster_bin_bp))
  vif_ok <- requireNamespace("car", quietly = TRUE)
  fit_pop <- function(p) {
    df <- data.frame(y = as.integer(uni[[outcome]]), specific = uni[[paste0("specific_", p)]],
                     mean_accessibility = uni$mean_accessibility, gc_content = uni$gc_content,
                     chr = uni$chr, distance_to_tss = uni$distance_to_tss)
    names(df)[1] <- outcome
    warns <- character(0)
    fit <- withCallingHandlers(glm(form, family = binomial(), data = df),
                               warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
    cf <- summary(fit)$coefficients; s <- cf["specificTRUE", ]
    est <- s[["Estimate"]]; se <- s[["Std. Error"]]
    big <- rownames(cf)[abs(cf[, "Estimate"]) > 15 | cf[, "Std. Error"] > 10]
    gv <- if (vif_ok) tryCatch({ v <- car::vif(fit); if (is.matrix(v)) v["specific", "GVIF"] else v[["specific"]] }, error = function(e) NA_real_) else NA_real_
    se_cl <- cluster_robust_se(fit, bins)[["specificTRUE"]]
    y <- df[[outcome]]
    res <- data.frame(cell_type = p, outcome = outcome, adjusted_OR = exp(est),
                      CI_low = exp(est - qnorm(0.975) * se), CI_high = exp(est + qnorm(0.975) * se),
                      log_OR = est, SE = se, z = s[["z value"]], p_wald = s[["Pr(>|z|)"]],
                      SE_cluster_robust_1Mb = se_cl, CI_low_cluster_robust = exp(est - qnorm(0.975) * se_cl),
                      CI_high_cluster_robust = exp(est + qnorm(0.975) * se_cl), p_cluster_robust = 2 * pnorm(-abs(est / se_cl)),
                      n_peaks_universe = nrow(df), n_specific = sum(df$specific), n_positive_specific = sum(y[df$specific]),
                      prop_positive_specific = mean(y[df$specific]), n_nonspecific = sum(!df$specific),
                      n_positive_nonspecific = sum(y[!df$specific]), prop_positive_nonspecific = mean(y[!df$specific]),
                      converged = fit$converged, iterations = fit$iter, all_coefficients_finite = all(is.finite(coef(fit))),
                      n_events = sum(y), quasi_separation_terms = paste(big, collapse = ";"),
                      n_warnings = length(warns), warnings = paste(unique(warns), collapse = " | "), GVIF_specific = gv,
                      null_deviance = fit$null.deviance, residual_deviance = fit$deviance, AIC = AIC(fit),
                      n_coefficients = length(coef(fit)), stringsAsFactors = FALSE)
    co <- data.frame(cell_type = p, outcome = outcome, term = rownames(cf), estimate = cf[, "Estimate"],
                     std_error = cf[, "Std. Error"], z = cf[, "z value"], p = cf[, "Pr(>|z|)"], stringsAsFactors = FALSE)
    log_msg(sprintf("  %s: OR %.3f", p, exp(est)))
    list(res = res, coef = co)
  }
  fits <- lapply(populations, fit_pop)
  enr <- do.call(rbind, lapply(fits, `[[`, "res")); coefs <- do.call(rbind, lapply(fits, `[[`, "coef"))
  enr$FDR_BH <- p.adjust(enr$p_wald, method = "BH")
  enr$interpretation <- ifelse(enr$FDR_BH < fdr & enr$adjusted_OR > 1, "enriched (FDR < 0.05)",
                        ifelse(enr$FDR_BH < fdr & enr$adjusted_OR < 1, "depleted (FDR < 0.05)", "not significant"))
  enr$model <- deparse1(form)
  enr <- enr[match(DISPLAY_ORDER, enr$cell_type), ]; rownames(enr) <- NULL
  rownames(coefs) <- NULL
  list(enrichment = enr, coefficients = coefs)
}

add_egfr_columns <- function(enr, egfr_path) {
  g <- read_tsv(egfr_path)
  enr$gchromvar_z <- g$z[match(enr$cell_type, g$cell_type)]
  enr$gchromvar_eGFR_status <- g$significance[match(enr$cell_type, g$cell_type)]
  enr
}

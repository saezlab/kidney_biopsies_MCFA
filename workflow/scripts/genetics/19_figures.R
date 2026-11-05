# (C) Sergio Alejandro Gómez
# Aim: Draw the genetics figures of one pt_C5 version from its result tables.
# Inputs: the tables folder of the version (egfr_enrichment.tsv, candidate_variants.tsv, motif_disruption_top30_tfs.tsv,
#   motif_pt_tal_by_tf.tsv, motif_pt_tal_similarity.tsv, eqtl_enrichment.tsv, caqtl_enrichment.tsv, ckd_enrichment.tsv,
#   control_libraries_10_populations.tsv, control_libraries_10_populations_summary.tsv).
# Outputs: <out_dir>/<figure>.pdf (vector) and .png (600 dpi) for egfr_enrichment, motif_disruption_top30, motif_pt_vs_tal,
#   eqtl_enrichment, caqtl_enrichment, ckd_enrichment and egfr_enrichment_control_libraries.
# Usage: Rscript workflow/scripts/genetics/19_figures.R --tables=results/genetics/with_ptC5/tables --include_ptC5=true --out_dir=results/genetics/with_ptC5/figures

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "R", "utils.R"))
suppressPackageStartupMessages(library(ggplot2))

args <- parse_args()
version <- if (arg_true(args, "include_ptC5", NULL)) "with_ptC5" else "without_ptC5"
tables <- arg(args, "tables")
out_dir <- arg(args, "out_dir")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tab <- function(name) read_tsv(file.path(tables, name))

# ---- Style shared by all figures ----
POP_ORDER <- c("TAL", "PT", "IC", "DCT", "CNT_PC", "PEC", "Stromal", "POD", "EC", "Immune")
POP_LABEL <- c(TAL = "TAL", PT = "PT", IC = "IC", DCT = "DCT", CNT_PC = "CNT/PC", PEC = "PEC", Stromal = "Stromal", POD = "POD",
               EC = "EC", Immune = "Immune")
VERSION_LAB <- c(with_ptC5 = "ArchR peak set, with pt_C5 nuclei", without_ptC5 = "ArchR peak set, without pt_C5 nuclei")[[version]]
CALL_COL <- c(enriched = "#B91D1E", ns = "#CDCCCC", depleted = "#3E4B9A")
CALL_LAB <- c(enriched = "Enriched (FDR < 0.05)", ns = "Not significant", depleted = "Depleted (FDR < 0.05)")
CALL_LAB_SHORT <- c(enriched = "enriched, FDR < 0.05", ns = "not significant", depleted = "depleted, FDR < 0.05")
FONT <- "Arial"
PT7 <- 7                      # base size of the small panels (eGFR, CKD, top-30 TFs, control libraries)
NAVY <- "#1F4E9E"
OUTLINE <- "grey15"
QTL_BREAKS <- c(0.6, 0.8, 1, 1.25, 1.5, 2, 2.5, 3)

theme_house <- function() {
  theme_classic(base_size = 9, base_family = FONT) +
    theme(plot.title = element_text(face = "bold", size = 9, hjust = 0), plot.title.position = "plot",
          plot.subtitle = element_text(colour = "grey45", size = 7.5, hjust = 0, margin = margin(b = 4)),
          plot.tag = element_text(face = "bold", size = 10), axis.text = element_text(colour = "black", size = 8),
          axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 7.5),
          legend.key.size = unit(3, "mm"))
}
theme_small <- function() {
  theme_classic(base_size = PT7, base_family = FONT) +
    theme(text = element_text(colour = "black"), axis.text = element_text(colour = "black", size = PT7 - 0.5),
          axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          plot.title = element_text(size = PT7 - 0.5, hjust = 0.5, margin = margin(b = 1)), plot.title.position = "plot",
          plot.subtitle = element_text(size = PT7 - 1, hjust = 0.5, colour = "grey45", margin = margin(b = 3)))
}
call_of <- function(stat, fdr, null = 0) {
  ifelse(!is.na(fdr) & fdr < 0.05 & stat > null, "enriched", ifelse(!is.na(fdr) & fdr < 0.05 & stat < null, "depleted", "ns"))
}
plus_lab <- function(x) ifelse(x > 0, paste0("+", x), as.character(x))
z_breaks <- function(xl) { b <- seq(-8, 12, 4); b[b >= xl[1] & b <= xl[2]] }

close_device <- function() {
  grid::convertHeight(grid::grobHeight(grid::textGrob("M", gp = grid::gpar(fontsize = 1.5))), "in")
  invisible(grDevices::dev.off())
}

save_fig <- function(p, name, width_mm, height_mm) {
  base <- file.path(out_dir, name)
  w <- width_mm / 25.4; h <- height_mm / 25.4
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = w, height = h, family = FONT); print(p); close_device()
  grDevices::png(paste0(base, ".png"), width = w, height = h, units = "in", res = 600, type = "cairo", family = FONT); print(p); close_device()
  log_msg("wrote ", base, ".pdf/.png")
}

# ---- eGFR and CKD enrichment: horizontal z bars, TAL at the top ----
z_limits <- function(z) c(min(-9, floor(min(z)) - 2), max(10, ceiling(max(z)) + 2))
z_bar_plot <- function(d, title, xl = z_limits(d$z)) {
  d$y <- factor(d$cell_type, levels = rev(POP_ORDER))
  d$call <- factor(d$call, levels = names(CALL_COL))
  off <- diff(xl) * 0.015
  ggplot(d, aes(x = z, y = y)) +
    geom_col(aes(fill = call), colour = OUTLINE, linewidth = 0.25, width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "black") +
    geom_text(aes(x = ifelse(z >= 0, z + off, z - off), label = sprintf("%.2f", z), hjust = ifelse(z >= 0, 0, 1)),
              size = (PT7 - 1) / .pt, family = FONT) +
    scale_fill_manual(values = CALL_COL, labels = CALL_LAB_SHORT, breaks = names(CALL_COL), drop = FALSE, name = NULL) +
    scale_x_continuous(breaks = z_breaks(xl), labels = plus_lab) +
    scale_y_discrete(labels = POP_LABEL) +
    coord_cartesian(xlim = xl, expand = FALSE, clip = "off") +
    labs(x = "Enrichment Z-score (gchromVAR)", y = NULL, title = title, subtitle = "Cell-type specific enrichment of fine-mapped variants") +
    theme_small() +
    # legend inside, lower right, starting just right of the zero line so that it never covers a negative bar
    theme(legend.position = "inside", legend.position.inside = c((0 - xl[1]) / diff(xl) + 0.02, 0.02), legend.justification = c(0, 0),
          legend.margin = margin(0, 0, 0, 0), legend.text = element_text(size = PT7 - 1.5), legend.key.size = unit(2.4, "mm"),
          legend.background = element_blank(), plot.margin = margin(4, 6, 3, 3))
}

egfr <- tab("egfr_enrichment.tsv")
stopifnot(setequal(egfr$cell_type, POP_ORDER), !anyNA(egfr$z))
if (version == "with_ptC5") stopifnot(identical(egfr$cell_type[order(-egfr$z)], POP_ORDER))
egfr <- egfr[match(POP_ORDER, egfr$cell_type), ]
egfr$call <- call_of(egfr$z, egfr$FDR)
save_fig(z_bar_plot(egfr, "eGFR GWAS signal enrichment"), "egfr_enrichment", 54, 70)

# ---- The 30 TFs disrupted by the most candidates ----
n_cand <- nrow(tab("candidate_variants.tsv"))
top <- tab("motif_disruption_top30_tfs.tsv")
top <- top[order(top$rank_by_distinct_snps), ][seq_len(min(30, nrow(top))), ]
top$tf <- factor(top$geneSymbol, levels = rev(top$geneSymbol))        # rank 1 at the top
p <- ggplot(top, aes(x = n_distinct_disrupted_snps, y = tf)) +
  geom_segment(aes(x = 0, xend = n_distinct_disrupted_snps, yend = tf), colour = "grey60", linewidth = 0.3) +
  geom_point(shape = 21, size = 1.3, fill = NAVY, colour = "black", stroke = 0.2) +
  scale_x_continuous(limits = c(0, max(top$n_distinct_disrupted_snps) * 1.08), breaks = scales::pretty_breaks(4), expand = expansion(mult = 0)) +
  labs(x = "Disrupted candidate variants", y = NULL, title = "Top genetically disrupted TFs",
       subtitle = sprintf("Strong motif disruptions by %s candidate variants", format(n_cand, big.mark = ","))) +
  theme_classic(base_size = PT7, base_family = FONT) +
  theme(text = element_text(colour = "black"), axis.text.x = element_text(colour = "black", size = PT7 - 0.5),
        axis.text.y = element_text(colour = "black", size = PT7 - 2), axis.ticks.y = element_blank(),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        plot.title = element_text(size = PT7 - 0.5, hjust = 0, margin = margin(b = 1)), plot.title.position = "plot",
        plot.subtitle = element_text(size = PT7 - 1.5, hjust = 0, colour = "grey45", margin = margin(b = 3)),
        plot.margin = margin(4, 5, 3, 3))
save_fig(p, "motif_disruption_top30", 46, 70)

# ---- PT- versus TAL-attributed candidates: fraction disrupted per TF ----
spread_labels <- function(y, gap, lo, hi) {
  o <- order(y); v <- y[o]
  for (it in 1:200) {
    moved <- FALSE
    for (k in seq_along(v)[-1]) if (v[k] - v[k - 1] < gap) { m <- (v[k] + v[k - 1]) / 2; v[k - 1] <- m - gap / 2; v[k] <- m + gap / 2; moved <- TRUE }
    v <- pmin(pmax(v, lo), hi)
    if (!moved) break
  }
  out <- numeric(length(y)); out[o] <- v; out
}
uncross <- function(x, y, lx, ly) {
  orient <- function(ax, ay, bx, by, cx, cy) sign((bx - ax) * (cy - ay) - (by - ay) * (cx - ax))
  crosses <- function(i, j) {
    o1 <- orient(x[i], y[i], lx, ly[i], x[j], y[j]); o2 <- orient(x[i], y[i], lx, ly[i], lx, ly[j])
    o3 <- orient(x[j], y[j], lx, ly[j], x[i], y[i]); o4 <- orient(x[j], y[j], lx, ly[j], lx, ly[i])
    o1 * o2 < 0 && o3 * o4 < 0
  }
  for (it in 1:500) {                                  # each swap shortens the total line length, so this ends
    swapped <- FALSE
    for (i in seq_along(x)) for (j in seq_along(x)) if (i < j && crosses(i, j)) { t <- ly[i]; ly[i] <- ly[j]; ly[j] <- t; swapped <- TRUE }
    if (!swapped) break
  }
  ly
}
ptt <- tab("motif_pt_tal_by_tf.tsv")
sim <- tab("motif_pt_tal_similarity.tsv")
sim <- sim[sim$group_A == "PT_ATTRIBUTED" & sim$group_B == "TAL_ATTRIBUTED" & sim$unit == "candidate_variant", ]
stopifnot(nrow(sim) == 1)
lab <- ptt[order(ptt$pooled_rank_canonical), ][1:10, ]
lim <- c(0, max(c(ptt$PT_fraction, ptt$TAL_fraction), na.rm = TRUE) * 1.12)
lab$lx <- max(lab$PT_fraction) + 0.04 * diff(lim)
lab$ly <- spread_labels(lab$TAL_fraction, gap = 0.045 * diff(lim), lo = lim[1], hi = lim[2])
lab$ly <- uncross(lab$PT_fraction, lab$TAL_fraction, lab$lx[1], lab$ly)
p <- ggplot(ptt, aes(PT_fraction, TAL_fraction)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = 2, colour = "grey40") +
  geom_point(size = 0.7, colour = "grey55", alpha = 0.7) +
  geom_point(data = lab, size = 1.2, colour = "black") +
  geom_segment(data = lab, aes(x = PT_fraction, y = TAL_fraction, xend = lx - 0.004, yend = ly), linewidth = 0.2, colour = "grey35") +
  geom_text(data = lab, aes(x = lx, y = ly, label = geneSymbol), size = 2.2, family = FONT, hjust = 0, vjust = 0.5) +
  coord_equal(xlim = c(lim[1], lim[2] + 0.16 * diff(lim)), ylim = lim, clip = "off") +
  labs(x = sprintf("Fraction of PT-attributed candidates (n = %d)", sim$n_units_A),
       y = sprintf("Fraction of TAL-attributed candidates (n = %d)", sim$n_units_B),
       title = "Predicted motif disruption in PT- and TAL-attributed candidates",
       subtitle = sprintf("%s;\none point per TF (n = %d); Spearman rho = %.2f; %d candidates in both strata;\nnon-exclusive strata, descriptive only; labelled: the 10 most frequent TFs over all candidates",
                          VERSION_LAB, nrow(ptt), sim$spearman_rho_fraction, sim$n_shared_units)) +
  theme_house()
save_fig(p, "motif_pt_vs_tal", 110, 116)

# ---- Kidney eQTL and caQTL enrichment of cell-type-specific peaks ----
qtl_limits <- function(v) {
  lo <- max(QTL_BREAKS[QTL_BREAKS <= min(v)], min(QTL_BREAKS)); hi <- min(QTL_BREAKS[QTL_BREAKS >= max(v)], max(QTL_BREAKS))
  c(min(lo, min(v)), max(hi, max(v)))
}
qtl_plot <- function(file, layer) {
  d <- tab(file)
  d <- d[match(POP_ORDER, d$cell_type), ]
  d <- d[!is.na(d$cell_type), ]
  stopifnot(!anyNA(d$p_cluster_robust))
  d$call <- call_of(d$adjusted_OR, d$FDR_BH, null = 1)
  d$pop <- factor(POP_LABEL[d$cell_type], levels = POP_LABEL)
  d$not_robust <- d$call != "ns" & d$p_cluster_robust >= 0.05
  log_msg(layer, " calls not robust to clustering: ", if (any(d$not_robust)) paste(d$cell_type[d$not_robust], collapse = ", ") else "none")
  ggplot(d, aes(pop, adjusted_OR, colour = call)) +
    geom_hline(yintercept = 1, linewidth = 0.3, linetype = 2, colour = "grey40") +
    geom_errorbar(aes(ymin = CI_low_cluster_robust, ymax = CI_high_cluster_robust), width = 0.3, linewidth = 0.3) +
    geom_linerange(aes(ymin = CI_low, ymax = CI_high), linewidth = 1.3) +
    geom_point(shape = 16, size = 2.2) +
    geom_point(data = d[d$not_robust, ], shape = 21, size = 2.2, stroke = 0.8, fill = "white", show.legend = FALSE) +
    scale_colour_manual(values = CALL_COL, labels = CALL_LAB, breaks = names(CALL_COL), drop = FALSE) +
    scale_y_log10(breaks = QTL_BREAKS, labels = function(x) formatC(x, format = "g", digits = 3),
                  limits = qtl_limits(c(d$CI_low, d$CI_high, d$CI_low_cluster_robust, d$CI_high_cluster_robust, 1))) +
    labs(x = NULL, y = "Adjusted odds ratio (95% CI)", title = sprintf("Kidney %s enrichment in cell-type-specific peaks", layer),
         subtitle = sprintf("%s;\nlogistic model per population, adjusted for accessibility, GC and TSS distance;\nBH FDR over %d populations; thick bars, Wald 95%% CI; thin bars, cluster-robust\n95%% CI (1-Mb bins); open symbols, calls not reproduced under the cluster-robust model",
                            VERSION_LAB, nrow(d))) +
    theme_house() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
save_fig(qtl_plot("eqtl_enrichment.tsv", "eQTL"), "eqtl_enrichment", 110, 143)
save_fig(qtl_plot("caqtl_enrichment.tsv", "caQTL"), "caqtl_enrichment", 110, 143)

# ---- CKD enrichment ----
ckd <- tab("ckd_enrichment.tsv")
ckd <- ckd[match(POP_ORDER, ckd$population), ]
stopifnot(!anyNA(ckd$population), !anyNA(ckd$z))
ckd <- data.frame(cell_type = ckd$population, z = ckd$z, call = call_of(ckd$z, ckd$FDR))
save_fig(z_bar_plot(ckd, "CKD GWAS signal enrichment", xl = z_limits(c(egfr$z, ckd$z))), "ckd_enrichment", 54, 70)

# ---- eGFR enrichment with and without the five control libraries ----
x <- tab("control_libraries_10_populations.tsv")
s <- tab("control_libraries_10_populations_summary.tsv")
stopifnot(setequal(x$population, POP_ORDER), nrow(x) == 10, !anyNA(c(x$z_baseline, x$z_exclusion)))
x <- x[match(POP_ORDER, x$population), ]
a1 <- sprintf("Analysis cohort, %d libraries (%s nuclei)", s$libraries_baseline, format(s$nuclei_baseline, big.mark = ","))
a2 <- sprintf("Without the 5 control libraries (%d libraries, %s nuclei)", s$libraries_exclusion, format(s$nuclei_exclusion, big.mark = ","))
d <- rbind(data.frame(pop = x$population, z = x$z_baseline, call = x$call_baseline, arm = a1),
           data.frame(pop = x$population, z = x$z_exclusion, call = x$call_exclusion, arm = a2))
d$pop <- factor(d$pop, levels = rev(POP_ORDER))
d$call <- factor(d$call, levels = names(CALL_COL))
d$arm <- factor(d$arm, levels = c(a1, a2))
d$y <- as.numeric(d$pop)
w <- data.frame(y = as.numeric(factor(x$population, levels = rev(POP_ORDER))), z1 = x$z_baseline, z2 = x$z_exclusion)
xl <- c(min(-8, floor(min(d$z)) - 1), max(10, ceiling(max(d$z)) + 1))
p <- ggplot(d, aes(z, y)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "black") +
  geom_segment(data = w, aes(x = z1, xend = z2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4, colour = "grey55") +
  geom_point(aes(fill = call, shape = arm), size = 2.1, stroke = 0.3, colour = OUTLINE) +
  scale_fill_manual(values = CALL_COL, labels = CALL_LAB_SHORT, breaks = names(CALL_COL), drop = FALSE, name = NULL) +
  scale_shape_manual(values = c(21, 24), name = NULL) +
  scale_y_continuous(breaks = seq_along(POP_ORDER), labels = POP_LABEL[rev(POP_ORDER)], expand = expansion(add = 0.5)) +
  scale_x_continuous(breaks = z_breaks(xl), labels = plus_lab) +
  coord_cartesian(xlim = xl, expand = TRUE) +
  guides(fill = guide_legend(order = 1, nrow = 1, override.aes = list(shape = 21, size = 2.1)),
         shape = guide_legend(order = 2, ncol = 1, override.aes = list(fill = "white"))) +
  labs(x = "Enrichment Z-score (gchromVAR)", y = NULL, title = "eGFR GWAS signal enrichment without control libraries",
       subtitle = sprintf("%s; BH FDR over 10 populations;\nbackground seed 2024; Spearman rho = %.2f", VERSION_LAB, s$spearman_rho_z_seed2024)) +
  theme_small() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0.5, "mm"), legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = PT7 - 1.5), legend.key.size = unit(2.6, "mm"), plot.margin = margin(4, 6, 3, 3))
save_fig(p, "egfr_enrichment_control_libraries", 80, 82)

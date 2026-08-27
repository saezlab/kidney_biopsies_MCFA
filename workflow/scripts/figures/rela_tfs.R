if (exists("snakemake")) {
  rela_tfs_path <- snakemake@input[["rela_tfs"]]
  factor_tfs_path <- snakemake@input[["factor_tfs"]]
  pdf_path <- snakemake@output[["pdf"]]
  png_path <- snakemake@output[["png"]]
  selected_view <- snakemake@params[["view"]]
  factor_label <- snakemake@params[["factor_label"]]
} else {
  rela_tfs_path <- "results/RelA_KO/collectri_TFs.csv"
  factor_tfs_path <- "results/MCFA/acts_collectri.csv"
  pdf_path <- "plots/figures/Figure5/rela_tfs.pdf"
  png_path <- "plots/figures/Figure5/rela_tfs.png"
  selected_view <- "PT"
  factor_label <- "Factor 1"
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(readr)
  library(tidyr)
})

rela_tfs <- readr::read_csv(rela_tfs_path, show_col_types = FALSE)
factor_tfs <- readr::read_csv(factor_tfs_path, show_col_types = FALSE) %>%
  dplyr::filter(.data$view == selected_view)

missing_rela_columns <- setdiff(c("TF", "padj", "act"), names(rela_tfs))
missing_factor_columns <- setdiff(c("TF", "padj", "act", "view"), names(factor_tfs))
if (length(missing_rela_columns) > 0L) {
  stop("RelA CollecTRI file is missing required columns: ", paste(missing_rela_columns, collapse = ", "), call. = FALSE)
}
if (length(missing_factor_columns) > 0L) {
  stop("MCFA CollecTRI file is missing required columns: ", paste(missing_factor_columns, collapse = ", "), call. = FALSE)
}
if (nrow(factor_tfs) == 0L) {
  stop("No MCFA CollecTRI TF activities found for view ", selected_view, ".", call. = FALSE)
}

to_plot <- dplyr::full_join(rela_tfs, factor_tfs, by = "TF", suffix = c("_KO", "_F1")) %>%
  dplyr::mutate(
    `padj < 0.05` = dplyr::case_when(
      .data$padj_F1 < 0.05 & .data$padj_KO < 0.05 ~ paste(factor_label, "& RelA KO"),
      .data$padj_F1 < 0.05 ~ factor_label,
      .data$padj_KO < 0.05 ~ "RelA KO",
      TRUE ~ "NS"
    ),
    `padj < 0.05` = factor(
      .data$`padj < 0.05`,
      levels = c(paste(factor_label, "& RelA KO"), factor_label, "RelA KO", "NS")
    ),
    lab = ifelse(.data$`padj < 0.05` %in% c(factor_label, paste(factor_label, "& RelA KO")), .data$TF, "")
  ) %>%
  tidyr::drop_na(act_F1, act_KO)

pal <- c(scales::pal_brewer(palette = "Dark2")(3), "#F0F0F0")

p <- to_plot %>%
  ggplot2::ggplot(ggplot2::aes(x = .data$act_F1, y = .data$act_KO)) +
  ggplot2::geom_vline(xintercept = 0, alpha = 0.3) +
  ggplot2::geom_hline(yintercept = 0, alpha = 0.3) +
  ggplot2::geom_point(ggplot2::aes(color = .data$`padj < 0.05`)) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, color = "darkgray", linewidth = 0.5) +
  ggrepel::geom_label_repel(ggplot2::aes(label = .data$lab), max.overlaps = 10) +
  ggpubr::stat_cor(label.x.npc = "center") +
  ggplot2::scale_colour_discrete(palette = pal) +
  ggplot2::labs(
    x = paste("Estimated TF activities,", factor_label, selected_view),
    y = "Estimated TF activities, RelA KO PT",
    color = "padj < 0.05"
  ) +
  ggplot2::theme_classic(base_size = 9)

if (!exists("snakemake")) {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
}

ggplot2::ggsave(pdf_path, p, width = 4.5, height = 4, units = "in")
ggplot2::ggsave(png_path, p, width = 4.5, height = 4, units = "in", dpi = 300)

if (exists("snakemake")) {
  enrichment_path <- snakemake@input[["enrichment"]]
  pdf_path <- snakemake@output[["pdf"]]
  png_path <- snakemake@output[["png"]]
  y_label <- snakemake@params[["y_label"]]
} else {
  enrichment_path <- "results/RelA_KO/factor1_enrichment.csv"
  pdf_path <- "plots/figures/Figure5/rela_factor_score.pdf"
  png_path <- "plots/figures/Figure5/rela_factor_score.png"
  y_label <- "Factor 1 PT Enrichment Score"
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(readr)
  library(stringr)
  library(tibble)
})

df <- readr::read_csv(enrichment_path, show_col_types = FALSE)
if (names(df)[1] %in% c("", "X", "...1")) {
  names(df)[1] <- "replicate"
}
if (!"replicate" %in% names(df)) {
  df <- tibble::rownames_to_column(as.data.frame(df), "replicate")
}

missing_columns <- setdiff(c("replicate", "act"), names(df))
if (length(missing_columns) > 0L) {
  stop("Factor enrichment file is missing required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

df <- df %>%
  dplyr::mutate(condition = ifelse(stringr::str_detect(.data$replicate, "_AA_"), "NT", "KO"))

p <- df %>%
  ggplot2::ggplot(ggplot2::aes(x = .data$condition, y = .data$act)) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::geom_jitter(height = 0, width = 0.12) +
  ggpubr::stat_compare_means(comparisons = list(c("KO", "NT")), method = "t.test") +
  ggplot2::scale_x_discrete(limits = c("NT", "KO")) +
  ggplot2::coord_cartesian(ylim = c(-11, 11)) +
  ggplot2::labs(x = NULL, y = y_label) +
  ggplot2::theme_classic(base_size = 9)

if (!exists("snakemake")) {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
}

ggplot2::ggsave(pdf_path, p, width = 3, height = 4, units = "in")
ggplot2::ggsave(png_path, p, width = 3, height = 4, units = "in", dpi = 300)

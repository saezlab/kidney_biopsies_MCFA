if (exists("snakemake")) {
  pathways_path <- snakemake@input[["pathways"]]
  pdf_path <- snakemake@output[["pdf"]]
  png_path <- snakemake@output[["png"]]
} else {
  pathways_path <- "results/RelA_KO/progeny_pathways.csv"
  pdf_path <- "plots/figures/Figure5/rela_pathways.pdf"
  png_path <- "plots/figures/Figure5/rela_pathways.png"
}

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(readr)
})

pathways <- readr::read_csv(pathways_path, show_col_types = FALSE)
missing_columns <- setdiff(c("Pathway", "padj", "act"), names(pathways))
if (length(missing_columns) > 0L) {
  stop("PROGENy pathways file is missing required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

p <- pathways %>%
  dplyr::mutate(sig = ifelse(.data$padj < 0.05, "*", "")) %>%
  ggplot2::ggplot(ggplot2::aes(x = .data$act, y = forcats::fct_reorder(.data$Pathway, .data$act))) +
  ggplot2::geom_col(fill = "#214C9B") +
  ggplot2::geom_text(
    ggplot2::aes(label = .data$sig, hjust = ifelse(.data$act < 0, 1.25, -0.25)),
    position = ggplot2::position_dodge(width = 1.5),
    vjust = 0.8,
    size = 15 / .pt
  ) +
  ggplot2::coord_cartesian(xlim = c(-55, 25)) +
  ggplot2::labs(x = "Enrichment score", y = "PROGENy Pathway") +
  ggplot2::theme_classic(base_size = 9)

if (!exists("snakemake")) {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
}

ggplot2::ggsave(pdf_path, p, width = 4, height = 3, units = "in")
ggplot2::ggsave(png_path, p, width = 4, height = 3, units = "in", dpi = 300)

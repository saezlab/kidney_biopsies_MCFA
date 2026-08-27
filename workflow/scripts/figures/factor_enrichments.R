if (exists("snakemake")) {
  collectri_path <- snakemake@input[["collectri"]]
  progeny_path <- snakemake@input[["progeny"]]
  hallmark_path <- snakemake@input[["hallmark"]]
  pdf_path <- snakemake@output[["pdf"]]
  png_path <- snakemake@output[["png"]]
  factor_name <- snakemake@params[["factor"]]
  top_hallmark_sets <- as.integer(snakemake@params[["top_hallmark_sets"]])
  top_collectri_tfs <- as.integer(snakemake@params[["top_collectri_tfs"]])
  view_order <- unlist(snakemake@params[["views"]])
} else {
  collectri_path <- "results/MCFA/acts_collectri.csv"
  progeny_path <- "results/MCFA/acts_progeny.csv"
  hallmark_path <- "results/MCFA/hallmark_Factor1.csv"
  pdf_path <- "plots/figures/Figure5/factor_enrichments.pdf"
  png_path <- "plots/figures/Figure5/factor_enrichments.png"
  factor_name <- "Factor1"
  top_hallmark_sets <- 10L
  top_collectri_tfs <- 20L
  view_order <- c("PT", "DCT", "PEC", "TAL", "CNT_PC", "Stromal")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(stringr)
})

read_csv_drop_index <- function(path) {
  tbl <- readr::read_csv(path, show_col_types = FALSE)
  if (names(tbl)[1] %in% c("", "X", "...1")) {
    tbl <- tbl[, -1, drop = FALSE]
  }
  names(tbl) <- make.names(names(tbl), unique = TRUE)
  tbl
}

require_columns <- function(tbl, required_columns, label) {
  missing_columns <- setdiff(required_columns, names(tbl))
  if (length(missing_columns) > 0L) {
    stop(label, " is missing required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
}

neg_log10 <- function(x) {
  -log10(pmax(as.numeric(x), .Machine$double.xmin))
}

collectri <- read_csv_drop_index(collectri_path)
progeny <- read_csv_drop_index(progeny_path)
hallmark <- read_csv_drop_index(hallmark_path)

require_columns(collectri, c("TF", "padj", "act", "Factor", "view"), "CollecTRI activities")
require_columns(progeny, c("Pathway", "padj", "act", "Factor", "view"), "PROGENy activities")
require_columns(hallmark, c("Term", "FDR.p.value", "view", "direction"), "Hallmark enrichment")

collectri <- collectri %>%
  dplyr::filter(.data$Factor == factor_name, .data$view %in% view_order) %>%
  dplyr::mutate(view = factor(.data$view, levels = view_order))
progeny <- progeny %>%
  dplyr::filter(.data$Factor == factor_name, .data$view %in% view_order) %>%
  dplyr::mutate(view = factor(.data$view, levels = view_order))
hallmark <- hallmark %>%
  dplyr::filter(.data$view %in% view_order) %>%
  dplyr::mutate(
    Term = stringr::str_remove(.data$Term, "^HALLMARK_"),
    view = factor(.data$view, levels = view_order)
  )

if (nrow(collectri) == 0L || nrow(progeny) == 0L || nrow(hallmark) == 0L) {
  stop(
    "No enrichment data were available for ",
    factor_name,
    " in requested views: ",
    paste(view_order, collapse = ", "),
    ".",
    call. = FALSE
  )
}

sets_up <- hallmark %>%
  dplyr::filter(.data$direction == "up") %>%
  dplyr::group_by(.data$Term) %>%
  dplyr::summarise(max_score = max(neg_log10(.data$FDR.p.value), na.rm = TRUE), .groups = "drop") %>%
  dplyr::slice_max(order_by = .data$max_score, n = top_hallmark_sets, with_ties = FALSE) %>%
  dplyr::pull(.data$Term)

sets_down <- hallmark %>%
  dplyr::filter(.data$direction == "down") %>%
  dplyr::group_by(.data$Term) %>%
  dplyr::summarise(max_score = max(neg_log10(.data$FDR.p.value), na.rm = TRUE), .groups = "drop") %>%
  dplyr::slice_max(order_by = .data$max_score, n = top_hallmark_sets, with_ties = FALSE) %>%
  dplyr::pull(.data$Term)

hallmark_fill_limit <- hallmark %>%
  dplyr::filter(
    (.data$direction == "up" & .data$Term %in% sets_up) |
      (.data$direction == "down" & .data$Term %in% sets_down)
  ) %>%
  dplyr::summarise(max_score = max(neg_log10(.data$FDR.p.value), na.rm = TRUE)) %>%
  dplyr::pull(.data$max_score)

p_hallmark_up <- hallmark %>%
  dplyr::filter(.data$direction == "up", .data$Term %in% sets_up) %>%
  dplyr::mutate(sig = ifelse(.data$FDR.p.value < 0.05, "*", "")) %>%
  ggplot2::ggplot(ggplot2::aes(x = .data$view, y = forcats::fct_reorder(.data$Term, neg_log10(.data$FDR.p.value)), fill = neg_log10(.data$FDR.p.value))) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(ggplot2::aes(label = .data$sig), color = "white") +
  ggplot2::scale_fill_viridis_c(option = "rocket", limits = c(0, hallmark_fill_limit)) +
  ggplot2::labs(subtitle = "Top 200 positive loadings", y = "Hallmark Gene Set", fill = "-log10 FDR") +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank())

p_hallmark_down <- hallmark %>%
  dplyr::filter(.data$direction == "down", .data$Term %in% sets_down) %>%
  dplyr::mutate(sig = ifelse(.data$FDR.p.value < 0.05, "*", "")) %>%
  ggplot2::ggplot(ggplot2::aes(x = .data$view, y = forcats::fct_rev(forcats::fct_reorder(.data$Term, neg_log10(.data$FDR.p.value))), fill = neg_log10(.data$FDR.p.value))) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(ggplot2::aes(label = .data$sig), color = "white") +
  ggplot2::scale_fill_viridis_c(option = "mako", limits = c(0, hallmark_fill_limit)) +
  ggplot2::labs(subtitle = "Top 200 negative loadings", x = "View", y = "Hallmark Gene Set", fill = "-log10 FDR") +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1, vjust = 1, angle = 45))

p_hallmark <- patchwork::wrap_plots(p_hallmark_up, p_hallmark_down, design = "1\n2") +
  patchwork::plot_layout(axis_titles = "collect")

max_abs_progeny <- max(abs(progeny$act), na.rm = TRUE)
p_progeny <- progeny %>%
  dplyr::arrange(dplyr::desc(.data$act)) %>%
  dplyr::mutate(
    Pathway = forcats::fct_rev(forcats::fct_inorder(.data$Pathway)),
    sig = ifelse(.data$padj < 0.05, "*", "")
  ) %>%
  ggplot2::ggplot(ggplot2::aes(y = .data$Pathway, x = .data$view, fill = .data$act)) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(ggplot2::aes(label = .data$sig), color = "white") +
  ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-max_abs_progeny, max_abs_progeny)) +
  ggplot2::labs(x = "View", y = "PROGENy Pathway", fill = "Activity") +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))

select_tfs <- collectri %>%
  dplyr::group_by(.data$TF) %>%
  dplyr::summarise(max_abs_score = max(abs(.data$act), na.rm = TRUE), .groups = "drop") %>%
  dplyr::slice_max(n = top_collectri_tfs, order_by = .data$max_abs_score, with_ties = FALSE) %>%
  dplyr::pull(.data$TF)

max_abs_collectri <- max(abs(collectri$act), na.rm = TRUE)
p_collectri <- collectri %>%
  dplyr::filter(.data$TF %in% select_tfs) %>%
  dplyr::arrange(dplyr::desc(.data$act)) %>%
  dplyr::mutate(
    TF = forcats::fct_rev(forcats::fct_inorder(.data$TF)),
    sig = ifelse(.data$padj < 0.05, "*", "")
  ) %>%
  ggplot2::ggplot(ggplot2::aes(y = .data$TF, x = .data$view, fill = .data$act)) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(ggplot2::aes(label = .data$sig), color = "white") +
  ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-max_abs_collectri, max_abs_collectri)) +
  ggplot2::labs(x = "View", y = "CollecTRI TF", fill = "Activity") +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))

p <- patchwork::wrap_plots(p_hallmark, p_progeny, p_collectri, design = "123")

if (!exists("snakemake")) {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
}

ggplot2::ggsave(pdf_path, p, width = 13, height = 7, units = "in")
ggplot2::ggsave(png_path, p, width = 13, height = 7, units = "in", dpi = 300)

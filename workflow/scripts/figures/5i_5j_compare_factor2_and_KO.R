library(magrittr)
library(ggplot2)

degs <- read.table("data/RelA_KO/DESEQ2.EC_AB.versus.EC_AA.txt", header = T)
tfs_ko <- read.csv("results/RelA_KO/collectri_TFs.csv")
tfs_f2 <- read.csv("results/MCFA/acts_collectri.csv") %>%
  dplyr::filter(view == "PT")
pathways_ko <- read.csv("results/RelA_KO/progeny_pathways.csv")
pathways_f2 <- read.csv("results/MCFA/acts_progeny.csv") %>%
  dplyr::filter(view == "PT")

# Pathways

pathways_ko %>%
  dplyr::mutate(sig = ifelse(padj < 0.05, "*", "")) %>%
  ggplot(aes(x = act, y = forcats::fct_reorder(Pathway, act))) +
  geom_col(fill = "#214C9B") + theme_classic() +
  geom_text(aes(label = sig, hjust = ifelse(act < 0, 1.25, -0.25)),
            position = position_dodge(width = 1.5),
            vjust = 0.8,
            size = 15 / .pt) + xlim(c(-55, 25)) +
  labs(x = "Enrichment score", y = "PROGENy Pathway")

to_plot <- dplyr::full_join(pathways_ko, pathways_f2, by = "Pathway", suffix = c("_KO", "_F2")) %>%
  dplyr::mutate("padj < 0.05" = ifelse(padj_F2 < 0.05 & padj_KO < 0.05, "Factor 2 & RelA KO",
                                       ifelse(padj_F2 < 0.05, "Factor 2",
                                              ifelse(padj_KO < 0.05, "RelA KO", "NS")))) %>%
  dplyr::mutate('padj < 0.05' = factor(`padj < 0.05`, levels = c("Factor 2 & RelA KO",
                                                                 "Factor 2",
                                                                 "RelA KO",
                                                                 "NS"))) %>%
  tidyr::drop_na('padj < 0.05') %>%
  dplyr::mutate(lab = ifelse(`padj < 0.05` %in% c("Factor 2", "RelAKO", "Factor 2 & RelA KO"), Pathway, ""))

# TFs
to_plot <- dplyr::full_join(tfs_ko, tfs_f2, by = "TF", suffix = c("_KO", "_F2")) %>%
  dplyr::mutate("padj < 0.05" = ifelse(padj_F2 < 0.05 & padj_KO < 0.05, "Factor 2 & RelA KO",
                                       ifelse(padj_F2 < 0.05, "Factor 2",
                                              ifelse(padj_KO < 0.05, "RelA KO", "NS")))) %>%
  dplyr::mutate('padj < 0.05' = factor(`padj < 0.05`, levels = c("Factor 2 & RelA KO",
                                                                 "Factor 2",
                                                                 "RelA KO",
                                                                 "NS"))) %>%
  tidyr::drop_na('padj < 0.05') %>%
  dplyr::mutate(lab = ifelse(`padj < 0.05` %in% c("Factor 2", "RelAKO", "Factor 2 & RelA KO"), TF, ""))

pal <- c(scales::pal_brewer(palette = "Dark2")(3), "#F0F0F0")

to_plot %>%
  ggplot(aes(x = act_F2, y = act_KO)) +
  geom_vline(xintercept = 0, alpha = 0.3) + 
  geom_hline(yintercept = 0, alpha = 0.3) +
  geom_point(aes(color = `padj < 0.05`)) + theme_classic() + 
  geom_smooth(method = "lm", se = FALSE, color = "darkgray", linewidth = 0.5) +
  ggrepel::geom_label_repel(aes(label = lab), max.overlaps = 10) +
  ggpubr::stat_cor(label.x.npc = "center") + scale_colour_discrete(palette = pal) +
  labs(x = "Estimated TF activities, Factor 2 PT", y = "Estimated TF activities, RelA KO PT")

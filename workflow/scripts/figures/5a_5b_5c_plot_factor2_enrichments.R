library(magrittr)
library(ggplot2)
library(patchwork)

factor2_collectri_TFs <- read.csv("results/MCFA/acts_collectri.csv")  %>%
  dplyr::mutate(view = ifelse(view == "aPT_aTAL_aDCT", "aPT", view))

factor2_progeny_pathways <- read.csv("results/MCFA/acts_progeny.csv")  %>%
  dplyr::mutate(view = ifelse(view == "aPT_aTAL_aDCT", "aPT", view))

factor2_hallmark <- read.csv('results/MCFA/hallmark_Factor2.csv')[,-1]%>%
  dplyr::mutate(view = ifelse(view == "aPT_aTAL_aDCT", "aPT", view))%>%
  dplyr::mutate(Term = stringr::str_remove(Term, "HALLMARK_"))

# Hallmark
sets_up <- factor2_hallmark %>%
  dplyr::filter(direction == "up") %>%
  dplyr::group_by(Term) %>%
  dplyr::summarise(max = max(-log10(FDR.p.value))) %>%
  dplyr::slice_max(order_by = max, n = 10) %>%
  dplyr::select(Term) %>% unlist()

sets_down <- factor2_hallmark %>%
  dplyr::filter(direction == "down") %>%
  dplyr::group_by(Term) %>%
  dplyr::summarise(max = max(-log10(FDR.p.value))) %>%
  dplyr::slice_max(order_by = max, n = 10) %>%
  dplyr::select(Term) %>% unlist()


p1 <- factor2_hallmark %>%
  dplyr::filter(direction == "up", Term %in% sets_up) %>%
  dplyr::mutate(sig = ifelse(FDR.p.value < 0.05, "*", "")) %>%
  ggplot(aes(x = view, y = forcats::fct_reorder(Term, -log10(FDR.p.value)), fill = -log10(FDR.p.value))) + geom_tile() +
  theme_classic() +
  geom_text(aes(label = sig), color = "white") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())  +
  labs(subtitle = "Top 200 positive loadings", y = "Hallmark Gene Set") +
  scale_fill_viridis_c(option = "rocket")#scale_fill_gradient(low = "white", high = "darkred")


p2 <- factor2_hallmark %>%
  dplyr::filter(direction == "down", Term %in% sets_down) %>%
  dplyr::mutate(sig = ifelse(FDR.p.value < 0.05, "*", "")) %>%
  ggplot(aes(x = view, y = forcats::fct_rev(forcats::fct_reorder(Term, -log10(FDR.p.value))), fill = -log10(FDR.p.value))) + geom_tile() +
  theme_classic() +
  geom_text(aes(label = sig), color = "white") +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45)) +
  labs(subtitle = "Top 200 negative loadings", y = "Hallmark Gene Set") +
  xlab("View") +
  scale_fill_viridis_c(option = "mako")# scale_fill_gradient(low = "white", high = "darkblue")

d <- "
1
2"
p_hallmark <- patchwork::wrap_plots(p1, p2, design = d) + plot_layout(axis_titles = "collect")


# PROGENy
max_abs_score <- max(abs(factor2_progeny_pathways$act))
(
  p_progeny <- factor2_progeny_pathways %>%
    dplyr::arrange(desc(act)) %>%
    dplyr::mutate(Pathway = forcats::fct_inorder(Pathway)) %>%
    dplyr::mutate(Pathway = forcats::fct_rev(Pathway)) %>%
    dplyr::mutate(sig = ifelse(padj < 0.05, "*", ""),
                  Score = act, View = view) %>%
    ggplot(aes(y = Pathway, x = View, fill = Score)) +
    geom_tile() + 
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "fill",
      limits = c(-max_abs_score, max_abs_score)
    )  + theme_classic() +
    geom_text(aes(label = sig), color = "white") +
    labs(x = "View", y = "PROGENy Pathway") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = rel(1)), axis.text.y = element_text(size = rel(1)))
    
)

# CollecTRI
select_tfs <- factor2_collectri_TFs %>%
  dplyr::group_by(TF) %>%
  dplyr::summarise(max_abs_score = max(abs(act))) %>%
  dplyr::slice_max(n = 20, order_by = max_abs_score) %>%
  dplyr::select(TF) %>% unlist()
max_abs_score <- max(abs(factor2_collectri_TFs$act))
(
  p_collectri <- factor2_collectri_TFs %>%
    dplyr::filter(TF %in% select_tfs) %>%
    dplyr::arrange(desc(act)) %>%
    dplyr::mutate(TF = forcats::fct_inorder(TF)) %>%
    dplyr::mutate(TF = forcats::fct_rev(TF)) %>%
    dplyr::mutate(sig = ifelse(padj < 0.05, "*", ""),
                  Score = act, View = view) %>%
    ggplot(aes(y = TF, x = View, fill = Score)) +
    geom_tile() + 
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "fill",
      limits = c(-max_abs_score, max_abs_score)
    )  + theme_classic() +
    geom_text(aes(label = sig), color = "white") +
    labs(x = "View", y = "CollecTRI TF") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = rel(1)), axis.text.y = element_text(size = rel(1)))
)

library(patchwork)
design = "
123"
patchwork::wrap_plots(p_hallmark, p_progeny, p_collectri, design = design)

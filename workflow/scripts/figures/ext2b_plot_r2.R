library(ggplot2)
library(magrittr)
library(dplyr)

if (exists("snakemake")){
  r2 <- read.csv(snakemake@input$r2)[,-1] # skip rownames (contains duplicates)
  data <- snakemake@output$data
  plt <- snakemake@output$plt1
  gg <- snakemake@output$gg1
} else {
  r2 <- read.csv('results/MCFA/r2_total~all~5.csv')
}

to_plot <- r2 %>%
  dplyr::filter(Factor == "Factor2") %>%
  dplyr::filter(!stringr::str_detect(View, "&"))
(
  p1 <- ggplot(to_plot, aes(x = forcats::fct_rev(View), y = R2, fill = R2))
  + geom_col()
  + scale_fill_viridis_c(option = "plasma", limits = c(0, 17))  + theme_classic()
  + labs(y = expr(R^2), x = 'View', subtitle = "Transcriptome views")
  + theme(legend.position = "none")
  + coord_flip()
)


to_plot2 <- r2 %>%
  dplyr::filter(Factor == "Factor2") %>%
  dplyr::filter(stringr::str_detect(View, "&")) %>%
  tidyr::separate_wider_delim(View, delim = "&", names = c("Sender", "Receiver"))
(
  p2 <- ggplot(to_plot2, aes(x = Receiver, y = forcats::fct_rev(Sender), fill = R2))
  + geom_tile()
  + scale_fill_viridis_c(option = "plasma", limits = c(0, 17))  + theme_classic()
  + labs(x = "Receiver", y = "Sender", subtitle = "L-R views") 
  + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
)


design <- "
112222
"
p1 + p2 + plot_layout(guides = "collect", design = design) + plot_annotation(title = "Variance explained - Factor 2")

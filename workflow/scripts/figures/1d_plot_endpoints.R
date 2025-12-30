library(magrittr)
library(ggplot2)
endpoints_atlas <- read.csv("results/data_overview/endpoint_atlas~heatmap.csv", row.names = 1)
proteinuria_over_time <- read.csv("results/data_overview/proteinuria_over_time.csv", row.names = 1)

df <- endpoints_atlas %>%
  dplyr::left_join(proteinuria_over_time) %>%
  dplyr::transmute(patient_ID,
                   disease,
                   dialysis_transplant = dialysis | transplant,
                   sustained_egfr_drop_40pc,
                   sustained_doubling_creatinine,
                   proteinuria_increase = proteinuria_over_time == "increased") %>%
  tibble::column_to_rownames("patient_ID") %>%
  dplyr::mutate(across(where(is.logical), ~ ifelse(.x, 1, 0))) %>%
  rowwise() %>%
  mutate(any = max(c_across(where(is.numeric))))

endpoints_long_atlas <- reshape2::melt(df, variable.name = "Endpoint",
                                     value.name = "Reached")
data <- endpoints_long_atlas %>%
  dplyr::group_by(disease, Endpoint, Reached) %>%
  dplyr::summarise(Count = n()) %>%
  dplyr::mutate(Reached = factor(ifelse(Reached == 1, "True", "False"), levels = c("True", "False"))) %>%
  dplyr::ungroup() %>%
  tidyr::complete(disease, Endpoint , Reached, fill = list(Count = 0)) %>%
  dplyr::group_by(disease, Endpoint) %>%
  dplyr::mutate(total = sum(Count)) %>%
  dplyr::filter(Reached == "True") %>%
  dplyr::mutate(percent_reached = Count/total*100) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(disease) %>%
  dplyr::mutate(any = max(Count)) %>%
  dplyr::mutate(disease = stringr::str_c(disease, "\n", any, "/", total, "\n(", sprintf("%.2f", (any/total)*100), "%)")) %>%
  dplyr::filter(Endpoint != "any")


(
  p1 <- ggplot(data, aes(x = forcats::fct_rev(Endpoint), y = forcats::fct_rev(Reached), fill = percent_reached)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "plasma", name = "Number of Patients", limits = c(0, 60)) + 
    labs(x = element_blank(), y = element_blank(), title = "Endpoint reached") +
    theme_minimal() +
    facet_wrap(~ disease, ncol = 5) + coord_flip() +
    geom_text(aes(label = stringr::str_c(Count, " (", sprintf("%.0f%%", percent_reached), ")")), color = "white", size = 3) +
    theme(legend.position = "none", axis.text.x = element_blank()) +
    scale_x_discrete(labels = c("dialysis_transplant" = "Dialysis/transplant",
                                "sustained_doubling_creatinine" = "Sustained doubling of creatinine",
                                "sustained_egfr_drop_40pc" = "Sustained 40% drop in eGFR",
                                "proteinuria_increase" = "Increased proteinuria"))
)



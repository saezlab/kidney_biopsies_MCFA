library(ggplot2)
library(magrittr)
library(patchwork)
metadata <- read.csv('data/data_overview/patient_metadata_including_xenium_and_duplicate_samples.csv', row.names = 1)
rna_annotations <- read.csv('data/annotations/biopsydata_withkpmp_final_annotations(in).csv', row.names = 1)
atac_annotations <- read.csv('data/annotations/label_transfer_RNA_ATAC_V3_clean.csv')

sample_group <- metadata %>%
  dplyr::transmute(sample_ID, Condition = disease) %>%
  dplyr::distinct()

existing_color_map <- c(
  "#F3786F", # 0_PT
  "#17BDCF", # 1_TAL
  "#36A3DC", # 2_CNT_PC
  "#6390CB", # 3_dPT_dTAL_dDCT
  "#8583BE", # 4_aPT
  "#9A7BB8", # 5_DCT
  "#AE74B2" # 6_EC
) # these are too close to one another

# Retrieved from https://sashamaps.net/docs/resources/20-colors/
palette_21 <- c('#e6194B', '#4363d8', '#ffe119',
                '#3cb44b', '#f58231', '#911eb4',
                '#42d4f4', '#f032e6', '#bfef45',
                '#fabed4', '#469990', '#dcbeff',
                '#9A6324', '#fffac8', '#800000',
                '#aaffc3', '#808000', '#ffd8b1',
                '#000075', '#a9a9a9', '#000000')


##########################
# RNA, by condition
##########################

rna_to_plot <- rna_annotations %>%
  dplyr::transmute(sample_ID = stringr::str_extract(barcode_combined_id, "^[:alnum:]*(?=\\_)"),
                   main_class) %>%
  dplyr::left_join(sample_group) %>%
  dplyr::group_by(Condition, main_class) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::group_by(Condition) %>%
  dplyr::mutate(total_by_condition = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(main_class = ifelse(main_class == "4_aPT_aTAL_aDCT",
                                    "4_aPT",
                                    ifelse(main_class == "18_aPT",
                                           "18_aPT2",
                                           main_class)),
                Proportion = n/total_by_condition,
                Condition = factor(Condition,
                                   levels = c("DN", "FSGS", "IgAN",
                                              "MN", "TIN", "AKI (control)",
                                              "control", "control_KPMP")))

p1 <- rna_to_plot %>%
  ggplot(aes(x=Condition,
             y = Proportion,
             fill = forcats::fct_reorder(main_class,
                                         as.numeric(sub("_.*", "", main_class))))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=palette_21) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1)) +
  labs(fill="Cluster annotation (RNA)")  + guides(fill=guide_legend(ncol=2))


##########################
# ATAC, by condition
##########################

palette_atac <- c('#e6194B', '#4363d8', '#ffe119',
                '#3cb44b', '#f58231', '#911eb4',
                '#42d4f4', '#f032e6', '#bfef45',
                '#fabed4', #'#469990', '#dcbeff',
                '#9A6324', '#fffac8',# '#800000',
               # '#aaffc3', 
                '#000075')
order_like_rna <- c("PT", "TAL", "CNT_PC", "dPT", "aPT", "DCT", "EC",
           "Fib_VSMC", "Immune", "IC", "Pod", "PEC", "PT_aPT")

atac_to_plot <- atac_annotations %>%
  dplyr::transmute(sample_ID = stringr::str_extract(X, "^[:alnum:]*(?=#)"),
                   cluster_label) %>%
  dplyr::left_join(sample_group) %>%
  dplyr::mutate(Condition = ifelse(sample_ID %in% c("K116", "K118", "K122"), "Other", Condition)) %>% # K122 BK virus nephropathy, K116 and K118 KIN
  dplyr::group_by(Condition, cluster_label) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::group_by(Condition) %>%
  dplyr::mutate(total_by_condition = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Proportion = n/total_by_condition,
                Condition = factor(Condition,
                                   levels = c("DN", "FSGS", "IgAN",
                                              "MN", "TIN", "AKI (control)", "Other",
                                              "control", "control_KPMP")))
p2 <- atac_to_plot %>%
  ggplot(aes(x=Condition,
             y = Proportion,
             fill = factor(cluster_label, levels = order_like_rna))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=palette_atac) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1)) +
  labs(fill="Cluster annotation (ATAC)") + guides(fill=guide_legend(ncol=2))

patchwork::wrap_plots(p1,p2)
##########################
# RNA, by Sample
##########################
rna_to_plot <- rna_annotations %>%
  dplyr::transmute(sample_ID = stringr::str_extract(barcode_combined_id, "^[:alnum:]*(?=\\_)"),
                   main_class) %>%
  dplyr::left_join(sample_group) %>%
  dplyr::mutate(Condition = factor(Condition,
                                   levels = c("DN", "FSGS", "IgAN",
                                              "MN", "TIN", "AKI (control)",
                                              "control", "control_KPMP"))) %>%
  dplyr::group_by(sample_ID, main_class) %>%
  dplyr::summarise(Condition = unique(Condition), n = n()) %>%
  dplyr::group_by(sample_ID) %>%
  dplyr::mutate(total_by_sample_ID = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(main_class = ifelse(main_class == "4_aPT_aTAL_aDCT",
                                    "4_aPT",
                                    ifelse(main_class == "18_aPT",
                                           "18_aPT2",
                                           main_class)),
                Proportion = n/total_by_sample_ID)

p3 <- rna_to_plot %>%
  ggplot(aes(x=forcats::fct_reorder(sample_ID, as.numeric(Condition)),
             y = n,
             fill = forcats::fct_reorder(main_class,
                                         as.numeric(sub("_.*", "", main_class))))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=palette_21) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1)) +
  labs(fill="Cluster annotation", y = "Nuclei", x = "Sample")

##########################
# ATAC, by sample
##########################

atac_to_plot <- atac_annotations %>%
  dplyr::transmute(sample_ID = stringr::str_extract(X, "^[:alnum:]*(?=#)"),
                   cluster_label) %>%
  dplyr::left_join(sample_group) %>%
  dplyr::mutate(Condition = factor(Condition,
                                   levels = c("DN", "FSGS", "IgAN",
                                              "MN", "TIN", "AKI (control)",
                                              "control"))) %>%
  dplyr::group_by(sample_ID, cluster_label) %>%
  dplyr::summarise(Condition = unique(Condition),n = n()) %>%
  dplyr::group_by(sample_ID) %>%
  dplyr::mutate(total_by_sample_ID = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Proportion = n/total_by_sample_ID)

p4 <- atac_to_plot %>%
  ggplot(aes(x=forcats::fct_reorder(sample_ID, as.numeric(Condition)),
             y = n,
             fill = factor(cluster_label, levels = order_like_rna))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=palette_atac) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1)) +
  labs(fill="Cluster annotation (ATAC)", y = "Nuclei", x = "Sample") +
  guides(fill=guide_legend(ncol=2))


design <- "
1111111111
2222223333"
patchwork::wrap_plots(p3, p4, plot_spacer(), design = design) + patchwork::plot_layout(guides = "collect")


library(magrittr)
library(dplyr)

if (exists("snakemake")){
  patient_metadata <- read.csv(snakemake@input$patient_metadata_EDA)
  clinical_data_extended <- read.csv(snakemake@input$clinical_data_extended)
  rna_annotations <- read.csv(snakemake@input$rna_annotations)
  atac_annotations <- read.csv(snakemake@input$atac_annotations)
  xenium_available <- read.table(snakemake@input$xenium_available, sep = ";", header = TRUE)
  proteinuria_over_time <- read.csv(snakemake@input$proteinuria_over_time)
  follow_up_csv <- snakemake@output$follow_up_csv
  age_csv <- snakemake@output$age_csv
  modality_csv <- snakemake@output$modality_csv
  modality_with_duplicates_csv <- snakemake@output$modality_with_duplicates_csv
  endpoint_atlas_csv <- snakemake@output$endpoint_atlas_csv
} else {
  patient_metadata <- read.csv('data/metadata/patient_metadata_EDA.csv')
  clinical_data_extended <- read.csv('data/metadata/clinical_data_extended.csv')
  rna_annotations <- read.csv('data/annotations/biopsydata_withkpmp_final_annotations(in).csv')
  atac_annotations <- read.csv('data/annotations/label_transfer_RNA_ATAC_V3_clean.csv', header=TRUE, row.names = 1)
  xenium_available <- read.table('data/metadata/spatial_availability.csv', sep = ";", header = TRUE)
  proteinuria_over_time <- read.csv('data/metadata/proteinuria_over_time.csv', row.names = 1)
  follow_up_csv <- 'results/data_overview/follow_up~boxplot.csv'
  age_csv <- 'results/data_overview/age~boxplot.csv'
  modality_csv <- 'results/data_overview/modality~heatmap.csv'
  modality_with_duplicates_csv <- 'results/data_overview/modality_with_duplicates~heatmap.csv'
  endpoint_atlas_csv <- 'results/data_overview/endpoint_atlas~heatmap.csv'
}

# Clarify which libraries passed QC in RNA/ATAC -> if they were annotated, they passed QC
# Reference tissue samples with low eGFR (i.e. AKI (control)) were not included in MCFA analysis, but
# annotated nonetheless. They were previously labelled "control_AKI".
low_eGFR_controls <- patient_metadata %>%
  dplyr::filter(disease == "control_AKI") %>%
  dplyr::select(sample_ID) %>% unlist()

patient_metadata$disease <- ifelse(patient_metadata$disease == "control_AKI",
                                   "AKI (control)",
                                   patient_metadata$disease)

rna_annotated <- rna_annotations %>%
  tidyr::separate_wider_delim(cols=X, delim="_", names=c("sample", "barcode")) %>%
  dplyr::select(sample) %>% unlist() %>% unique()

atac_annotated <- atac_annotations %>%
  tibble::rownames_to_column("cellid") %>%
  tidyr::separate_wider_delim(cols=cellid, delim="#", names=c("sample", "barcode")) %>%
  dplyr::select(sample) %>% unlist() %>% unique()

rna <- patient_metadata %>%
  dplyr::filter(!is.na(ck_number_rna)) %>%
  dplyr::filter(!is.na(disease)) %>%
  dplyr::filter(sample_ID %in% rna_annotated)# Ensure final figures only include annotated samples

atac <- patient_metadata %>%
  dplyr::filter(!is.na(ck_number_atac)) %>%
  dplyr::filter(sample_ID %in% atac_annotated) %>% # Ensure final figures only include annotated samples
  dplyr::filter(!is.na(disease))

kpmp_controls <- data.frame(sample_ID = rna_annotated[stringr::str_detect(rna_annotated, "KB")],
                            disease = "control_KPMP") %>%
  dplyr::mutate(patient_ID = sample_ID) # To avoid doing problems when grouping by patient

# Add Xenium information to metadata
xenium_info <- xenium_available %>%
  dplyr::transmute(sample_ID = sample_ID,
                   xenium_available = Spatial.available == "yes")
patient_metadata <- patient_metadata %>%
  dplyr::full_join(xenium_info) %>%
  dplyr::mutate(disease = ifelse(stringr::str_detect(sample_ID, "Apollot0"), "control_preKTX", disease))

xenium <- patient_metadata %>%
  dplyr::filter(xenium_available) %>%
  dplyr::filter(!is.na(disease)) %>%
  dplyr::mutate(patient_ID = ifelse(stringr::str_detect(sample_ID, "Apollot0"), sample_ID, patient_ID)) # To avoid problems grouping by patient - currently going on the information that there is a one-to-one mapping between sample and patient 25.06.2025

metadata_with_duplicates <- dplyr::bind_rows(rna, atac, kpmp_controls, xenium) %>%
  dplyr::distinct() %>%
  dplyr::filter(!disease %in% c("BKVAN", "KIN")) %>%
  dplyr::mutate(rna = sample_ID %in% rna_annotated,#rna_annotated_loweGFRcontrols_removed,
                atac = sample_ID %in% atac_annotated,
                xenium = sample_ID %in% xenium$sample_ID) %>%
  dplyr::mutate(modality = ifelse((rna & atac & xenium), "RNA, ATAC, Xenium",
                                  ifelse((rna & atac & !xenium), "RNA, ATAC",
                                         ifelse((rna & xenium & !atac), "RNA, Xenium",
                                                ifelse((rna & !xenium & !atac), "RNA",
                                                       ifelse((atac & xenium & !rna), "ATAC, Xenium",
                                                              ifelse((atac & !rna & !xenium), "ATAC",
                                                                     ifelse((xenium & !rna & !atac), "Xenium", NA)))))))) %>%
  dplyr::mutate(disease = ifelse(disease == "IgA", "IgAN", disease))

metadata <- metadata_with_duplicates %>%  dplyr::group_by(patient_ID) %>% # For the purposes of these plots, we want
  dplyr::slice_head(n=1) %>% # to list patients rather than samples.
  dplyr::ungroup() # So remove duplicate entries where we have two samples for one patient.

write.csv(metadata_with_duplicates, "results/data_overview/patient_metadata_including_xenium_and_duplicate_samples.csv")
#######################
# Follow-up boxplot
#######################

# Rename
clinical_data_extended$disease <- ifelse(clinical_data_extended$disease == "control_AKI",
                                   "AKI (control)",
                                   clinical_data_extended$disease)

clinical_data <- clinical_data_extended %>%
  dplyr::filter(patient_ID %in% c(rna$patient_ID, atac$patient_ID, xenium$patient_ID)) %>%
  dplyr::filter(!disease %in% c("control", "AKI (control)", "control_KPMP", "KIN", "BKVAN")) %>%
  dplyr::mutate(disease = ifelse(disease == "IgA", "IgAN", disease))

rna_patientIDs <- metadata %>%
  dplyr::select(patient_ID, sample_ID) %>%
  dplyr::filter(sample_ID %in% rna_annotated) %>% #rna_annotated_loweGFRcontrols_removed) %>%
  dplyr::select(patient_ID) %>% unlist()

atac_patientIDs <- metadata %>%
  dplyr::select(patient_ID, sample_ID) %>%
  dplyr::filter(sample_ID %in% atac_annotated) %>%
  dplyr::select(patient_ID) %>% unlist()

xenium_patientIDs <- metadata %>%
  dplyr::select(patient_ID, sample_ID) %>%
  dplyr::filter(sample_ID %in% xenium$sample_ID) %>%
  dplyr::select(patient_ID) %>% unlist()

follow_up_time <-  clinical_data %>%
  dplyr::group_by(patient_ID) %>%
  dplyr::slice_max(days_from_biopsy) %>%
  dplyr::select(patient_ID, disease, days_from_biopsy) %>%
  dplyr::mutate(rna = patient_ID %in% rna_patientIDs,
                atac = patient_ID %in% atac_patientIDs,
                xenium = patient_ID %in% xenium_patientIDs) %>%
  dplyr::mutate(modality = ifelse((rna & atac & xenium), "RNA, ATAC, Xenium",
                                  ifelse((rna & atac & !xenium), "RNA, ATAC",
                                         ifelse((rna & xenium & !atac), "RNA, Xenium",
                                                ifelse((rna & !xenium & !atac), "RNA",
                                                       ifelse((atac & xenium & !rna), "ATAC, Xenium",
                                                              ifelse((atac & !rna & !xenium), "ATAC",
                                                                     ifelse((xenium & !rna & !atac), "Xenium", NA)))))))) %>%
  dplyr::distinct()



write.csv(follow_up_time, follow_up_csv, row.names = F)

# follow_up_time_PCR <- clinical_data %>%
#   dplyr::filter(type %in% c("Protein-creatinine ratio (sample, recalculated)",
#                             "Protein-creatinine ratio (collection, recalculated)")) %>%
#   dplyr::group_by(patient_ID) %>%
#   dplyr::slice_max(days_from_biopsy) %>%
#   dplyr::select(patient_ID, disease, days_from_biopsy) %>%
#   dplyr::mutate(rna = patient_ID %in% rna_patientIDs,
#                 atac = patient_ID %in% atac_patientIDs,
#                 xenium = patient_ID %in% xenium_patientIDs) %>%
#   dplyr::mutate(modality = ifelse((rna & atac & xenium), "RNA, ATAC, Xenium",
#                                   ifelse((rna & atac & !xenium), "RNA, ATAC",
#                                          ifelse((rna & xenium & !atac), "RNA, Xenium",
#                                                 ifelse((rna & !xenium & !atac), "RNA",
#                                                        ifelse((atac & xenium & !rna), "ATAC, Xenium",
#                                                               ifelse((atac & !rna & !xenium), "ATAC",
#                                                                      ifelse((xenium & !rna & !atac), "Xenium", NA)))))))) %>%
#   dplyr::distinct()
# 
# follow_up_time_eGFR <- clinical_data %>%
#   dplyr::filter(type %in% c("eGFR (recalculated)")) %>%
#   dplyr::group_by(patient_ID) %>%
#   dplyr::slice_max(days_from_biopsy) %>%
#   dplyr::select(patient_ID, disease, days_from_biopsy) %>%
#   dplyr::mutate(rna = patient_ID %in% rna_patientIDs,
#                 atac = patient_ID %in% atac_patientIDs,
#                 xenium = patient_ID %in% xenium_patientIDs) %>%
#   dplyr::mutate(modality = ifelse((rna & atac & xenium), "RNA, ATAC, Xenium",
#                                   ifelse((rna & atac & !xenium), "RNA, ATAC",
#                                          ifelse((rna & xenium & !atac), "RNA, Xenium",
#                                                 ifelse((rna & !xenium & !atac), "RNA",
#                                                        ifelse((atac & xenium & !rna), "ATAC, Xenium",
#                                                               ifelse((atac & !rna & !xenium), "ATAC",
#                                                                      ifelse((xenium & !rna & !atac), "Xenium", NA)))))))) %>%
#   dplyr::distinct()
# 
# df1 <- follow_up_time %>% dplyr::rename(last_follow_up_any = days_from_biopsy)
# df2 <- follow_up_time_PCR %>% dplyr::rename(last_follow_up_PCR = days_from_biopsy)
# df3 <- follow_up_time_eGFR %>% dplyr::rename(last_follow_up_eGFR = days_from_biopsy)
# 
# df <- df1 %>% 
#   dplyr::full_join(df2) %>%
#   dplyr::full_join(df3) %>%
#   tidyr::pivot_longer(cols = starts_with("last_follow_up"),
#                       names_to = "type",
#                       values_to = "days_from_biopsy")
# 
# write.csv(df, 'results/Figure1/follow_up_eGFR_PCR_all.csv')
#######################
# Age at biopsy boxplot
#######################

age <- metadata %>%
  dplyr::select(patient_ID, disease, age_at_biopsy, modality) %>%
  tidyr::drop_na() # Remove NAs (KPMP controls)

write.csv(age, age_csv, row.names = F)

#######################
# Modality heatmap
#######################
modality <- metadata %>%
  dplyr::select(c(patient_ID, disease, modality, rna, atac, xenium)) %>%
  dplyr::distinct()

write.csv(modality, modality_csv, row.names = F)

modality_with_duplicates <- metadata_with_duplicates %>%
  dplyr::select(c(sample_ID, patient_ID, disease, modality, rna, atac, xenium)) %>%
  dplyr::distinct() %>%
  dplyr::group_by(patient_ID) %>%
  dplyr::summarise(disease = unique(disease), rna = sum(rna), atac = sum(atac), xenium = sum(xenium))

write.csv(modality_with_duplicates, modality_with_duplicates_csv, row.names = F)
#######################
# Endpoint heatmap
#######################
# Take mean measurement in the first seven days after biopsy as baseline
eGFR_at_baseline <- clinical_data %>%
  dplyr::filter(days_from_biopsy >= 0 & days_from_biopsy <= 7) %>%
  dplyr::filter(type == "eGFR (recalculated)") %>%
  dplyr::group_by(patient_ID) %>%
  dplyr::slice_min(order_by = days_from_biopsy, n = 3) %>%
  dplyr::summarise(baseline = mean(measurement))

eGFR_drop_40pc_from_baseline <- clinical_data %>%
  dplyr::filter(days_from_biopsy >= 0) %>%
  dplyr::filter(type == "eGFR (recalculated)") %>%
  dplyr::left_join(eGFR_at_baseline) %>%
  dplyr::group_by(patient_ID) %>%
  dplyr::mutate(dropped = ifelse(measurement < 0.6*baseline, TRUE, FALSE))

sustained_drop <- eGFR_drop_40pc_from_baseline %>%
  dplyr::ungroup() %>%
  dplyr::mutate(consecutive_dropped_group = ifelse(dropped,
                                                   data.table::rleid(eGFR_drop_40pc_from_baseline$dropped),
                                                   NA)) %>%
  dplyr::filter(!is.na(consecutive_dropped_group)) %>%
  dplyr::group_by(consecutive_dropped_group) %>%
  dplyr::summarise(patient_ID = unique(patient_ID),
                   baseline = unique(baseline),
                   first_dropped_value = min(measurement[which(days_from_biopsy == min(days_from_biopsy))]),
                   first_dropped_time = min(days_from_biopsy),
                   n_consec_dropped_measurements = n(),
                   days_apart = max(days_from_biopsy) - min(days_from_biopsy)) %>%
  dplyr::filter(n_consec_dropped_measurements >= 2, days_apart >= 30)%>%
  dplyr::group_by(patient_ID) %>%
  dplyr::slice_min(order_by = first_dropped_time, n = 1)

# ---
# Creatinine doubling from baseline
# ---
# Take mean measurement of first seven days as baseline
creatinine_at_baseline <- clinical_data %>%
  dplyr::filter(days_from_biopsy >= 0 & days_from_biopsy <= 7) %>%
  dplyr::filter(type == "Creatinine (serum)") %>%
  dplyr::group_by(patient_ID) %>%
  dplyr::slice_min(order_by = days_from_biopsy, n = 3) %>%
  dplyr::summarise(baseline = mean(measurement))

creatinine_doubling_from_baseline <- clinical_data %>%
  dplyr::filter(days_from_biopsy >= 0) %>%
  dplyr::filter(type == "Creatinine (serum)") %>%
  dplyr::left_join(creatinine_at_baseline) %>%
  dplyr::group_by(patient_ID) %>%
  dplyr::mutate(doubled = ifelse(measurement > 2*baseline, TRUE, FALSE))

# Find cases where TRUE for at least two consecutive measurements at least 30 days apart
sustained_doubling <- creatinine_doubling_from_baseline %>%
  dplyr::ungroup() %>%
  dplyr::mutate(consecutive_doubled_group = ifelse(doubled,
                                                   data.table::rleid(creatinine_doubling_from_baseline$doubled),
                                                   NA)) %>%
  dplyr::filter(!is.na(consecutive_doubled_group)) %>%
  dplyr::group_by(consecutive_doubled_group) %>%
  dplyr::summarise(patient_ID = unique(patient_ID),
                   baseline = unique(baseline),
                   first_doubling_value = min(measurement[which(days_from_biopsy == min(days_from_biopsy))]),
                   first_doubling_time = min(days_from_biopsy),
                   n_consec_doubled_measurements = n(),
                   days_apart = max(days_from_biopsy) - min(days_from_biopsy)) %>%
  dplyr::filter(n_consec_doubled_measurements >= 2, days_apart >= 30)%>%
  dplyr::group_by(patient_ID) %>%
  dplyr::slice_min(order_by = first_doubling_time, n = 1)

endpoints_rna <- metadata %>%
  dplyr::filter(rna) %>%
  dplyr::filter(!disease %in% c("control", "AKI (control)", "control_KPMP")) %>%
  dplyr::select(patient_ID, disease, dialysis, transplant) %>%
  dplyr::mutate(#sustained_egfr_below_10 = patient_ID %in% sustained_below_10$patient_ID,
                sustained_egfr_drop_40pc = patient_ID %in% sustained_drop$patient_ID,
                sustained_doubling_creatinine = patient_ID %in% sustained_doubling$patient_ID) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("patient_ID") %>%
  dplyr::mutate(across(where(is.logical), ~ ifelse(.x, 1, 0))) %>%
  rowwise() %>%
  mutate(any = max(c_across(where(is.numeric))))

endpoints_atlas <- metadata %>%
    dplyr::filter(!disease %in% c("control", "AKI (control)", "control_KPMP", "control_preKTX")) %>%
  dplyr::select(patient_ID, disease, dialysis, transplant, sample_ID) %>%
  dplyr::mutate(#sustained_egfr_below_10 = patient_ID %in% sustained_below_10$patient_ID,
                sustained_egfr_drop_40pc = patient_ID %in% sustained_drop$patient_ID,
                sustained_doubling_creatinine = patient_ID %in% sustained_doubling$patient_ID) %>%
  dplyr::distinct()

write.csv(endpoints_atlas, endpoint_atlas_csv)

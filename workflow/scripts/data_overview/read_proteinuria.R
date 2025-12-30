patient_ids_with_proteinuria <- read.csv("data/Figure1/patient_selection_with_proteinuria.csv", skip = 1) %>%
  janitor::clean_names() %>%
  dplyr::transmute(sample_ID = stringr::str_c("K", id_number), proteinuria_over_time)

write.csv(patient_ids_with_proteinuria, "data/Figure1/proteinuria_over_time.csv")

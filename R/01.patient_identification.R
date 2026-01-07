## Import library
library(tidyverse)
library(janitor)
library(gtsummary)
theme_set(theme_classic())
theme_gtsummary_compact()

###################################################################### I ### Load data
path_clinical <- fs::path("", "Volumes", "Gillis_Research", "Lab_Data", "CHEvolutionTMN")
clinical_data <-
  readxl::read_xlsx(paste0(here::here(), #path_clinical, 
                           "/data/raw data/CHevolution_CohortRaw_10R25000197_v5_20251024.xlsx"), 
                           # "/RawData/CHevolution_CohortRaw_10R25000197_v5_20251024.xlsx"), 
                    sheet = "Base_Data") %>% 
  clean_names()
cancer_data <-
  readxl::read_xlsx(paste0(here::here(), #path_clinical, 
                           "/data/raw data/CHevolution_CohortRaw_10R25000197_v5_20251024.xlsx"), 
                    # "/RawData/CHevolution_CohortRaw_10R25000197_v5_20251024.xlsx"), 
                    sheet = "Registry_General") %>% 
  clean_names()
registry_treatment <-
  readxl::read_xlsx(paste0(here::here(), #path_clinical, 
                           "/data/raw data/CHevolution_CohortRaw_10R25000197_v5_20251024.xlsx"), 
                    # "/RawData/CHevolution_CohortRaw_10R25000197_v5_20251024.xlsx"), 
                    sheet = "Registry_Treatment") %>% 
  clean_names()
emr_medication <-
  readxl::read_xlsx(paste0(here::here(), #path_clinical, 
                           "/data/raw data/CHevolution_CohortRaw_10R25000197_v5_20251024.xlsx"), 
                    # "/RawData/CHevolution_CohortRaw_10R25000197_v5_20251024.xlsx"), 
                    sheet = "Emr_Medication") %>% 
  clean_names()
dna_data <-
  readxl::read_xlsx(paste0(here::here(), #path_clinical, 
                           "/data/raw data/CHevolution_CohortRaw_10R25000197_v5_20251024.xlsx"), 
                    # "/RawData/CHevolution_CohortRaw_10R25000197_v5_20251024.xlsx"), 
                    sheet = "Biobanking") %>% 
  clean_names()
molecular_ngs <-
  readxl::read_xlsx(paste0(here::here(), #path_clinical, 
                           "/data/raw data/CHevolution_CohortRaw_10R25000197_v5_20251024.xlsx"), 
                    # "/RawData/CHevolution_CohortRaw_10R25000197_v5_20251024.xlsx"), 
                    sheet = "Molecular_NGS") %>% 
  clean_names()


###################################################################### II ### Data claening
table(dna_data$sample_type)
# Buffy Coat               Genomic DNA                       MNC           MNC less CD138+ 
#   19                        23                       407                        18 
# MNC Slide            Paraffin Block                    Plasma       Plasma Cells CD138+ 
#   150                        13                       177                        22 
# Serum             Stained Slide           Total RNA >10nt          Total RNA >200nt 
# 69                        14                         1                         1 
# Unprocessed Liquid Tissue           Unstained Slide 
# 64                        14 
table(dna_data$best_anatomic_site)
table(dna_data$best_category)

dna_data <- dna_data %>% 
  mutate(mrn = as.character(mrn)) %>% 
  # Select sample type needed for CH detection
  filter(sample_type %in% c(
    "Buffy Coat", "Genomic DNA","MNC", "PBMC",
    "MNC less CD138+", "Unprocessed Liquid Tissue")
  ) %>% 
  select(mrn, sample_family_id, sample_id,
         specimen_collection_dt = collection_dt,
         best_anatomic_site, sample_type) %>%
  # add same sample/same date on the same row
  arrange(mrn, specimen_collection_dt, best_anatomic_site) %>% 
  # Summarize to have 1 sample/day per row and not 1 row for each aliquot of the same sample 
  group_by(mrn, sample_family_id, specimen_collection_dt) %>% 
  summarise_at(vars(sample_id, best_anatomic_site, sample_type), str_c, collapse = "; ") %>%
  arrange(mrn, specimen_collection_dt, best_anatomic_site) %>% 
  group_by(mrn, specimen_collection_dt) %>% 
  summarise_at(vars(sample_family_id, sample_id, best_anatomic_site, sample_type), str_c, collapse = "; ") %>%
  ungroup() #%>% 
  # group_by(mrn) %>% 
  # mutate(number_of_samples = n()) %>% 
  # ungroup() %>% 
  # filter(number_of_samples > 1) %>% 
  # select(-number_of_samples)

write_rds(dna_data, paste0(here::here(),
                           "/data/processed_data",
                           "/dna_data_", today(),".rds"))

set.seed(123)
clinical_data_1 <- clinical_data %>% # will not be in the sample data
  mutate(mrn = as.character(mrn)) %>% 
  # generate de-identy ids
  mutate(study_name = "ch_evolution_") %>%
  group_by(mrn = factor(mrn, levels = sample(unique(mrn)))) %>%
  mutate(random_id = cur_group_id()) %>%
  ungroup() %>%
  mutate(zero = 6 - nchar(random_id)) %>%
  mutate(add_zero = stringi::stri_dup("0", zero)) %>%
  select(c(study_name, add_zero, random_id, mrn, everything())) %>%
  unite(deidentified_patient_id, study_name:random_id, sep = "", remove = FALSE)

# write_csv(clinical_data_1 %>%
#             select(mrn, deidentified_patient_id, random_id),
#           paste0("data/processed_data",
#                  "/CHevolution_Match_mrn_and_deidentified_ids_for_the_CHevolution_project_before_filtering",
#                  today(), ".csv"))
# write_csv(clinical_data_1 %>%
#             select(mrn, deidentified_patient_id, random_id),
#           paste0(path_clinical, "/ProcessedData",
#                  "/CHevolution_Match_mrn_and_deidentified_ids_for_the_CHevolution_project_before_filtering",
#                  today(), ".csv"))


clinical_data_1 <- clinical_data_1 %>%
  # Clean up var created for coding de-id ids
  select(c(deidentified_patient_id, mrn, everything(), -c(zero, study_name, add_zero))) %>% 
  # Filtered out patients without blood samples of interest
  # filter(str_detect(mrn, paste0(dna_data$mrn, collapse = "|"))) %>% 
  select(#deidentified_patient_id, 
         mrn, 
         cdsc_vital_status, last_contact_or_death_dt, last_visit_dt, 
         # non-myeloid info
         non_myeloid_tumor_id, non_myeloid_tumor_seq_num, 
         # non_myeloid_histology_cd, non_myeloid_dx_dt, # added back by code later
         non_myeloid_treatment_type,
         first_non_myeloid_treatment_start_dt, 
         # myeloid info
         # myeloid_tumor_id, myeloid_tumor_seq_num, 
         first_myeloid_dx_dt, 
         time_treatment_to_myeloid_dx = treatment_to_myeloid_dx#,
         # starts_with("myeloid_histology_desc_"),
         # starts_with("myeloid_histology_cd_")
         ) %>% 
  distinct() %>% 
  mutate(non_myeloid_treatment_type = str_to_lower(non_myeloid_treatment_type)) %>% 
  group_by(across(c(-non_myeloid_treatment_type))) %>% 
  summarise_at(vars(non_myeloid_treatment_type), str_c, collapse = "/") %>%
  ungroup()


cancer_data <- cancer_data %>% 
  mutate(mrn = as.character(mrn)) %>% 
  arrange(mrn, dx_dt) %>% 
  select(-starts_with("patient"))

# Patient exclusion based on histology info
cancer_data %>% 
  filter(tumor_seq_num == "1") %>% 
  select(histology_desc, primary_site_group_desc) %>% 
  tbl_summary(sort = everything() ~ "frequency")

cancer_data1 <- cancer_data %>%
  # Verify tumor_seq_num - 1 is wrongly coded as 0 + remove diagnosis 60 as it is non cancer dx -
  filter(tumor_seq_num != 60) %>% 
  group_by(mrn) %>%
  mutate(tumor_seq_num = row_number(),
         .after = tumor_seq_num) %>%
  ungroup() %>%
  # select(-tumor_seq_num) %>% 
  # Exclude patients who only had myeloid diseases
  mutate(exclude = case_when(
    tumor_seq_num == "1" &
      histology_desc %in% c(
        "REFRACTORY ANEMIA WITH EXCESS BLASTS", 
        "MYELODYSPLASTIC SYNDROME NOS", 
        "CHRONIC MYELOMONOCYTIC LEUKEMIA NOS", 
        "POLYCYTHEMIA VERA")                             ~ "Exclude"
  )) %>% 
  group_by(mrn) %>% 
  fill(exclude, .direction = "updown") %>% 
  ungroup() %>% 
  filter(is.na(exclude)) %>% # None are removed with the updated data 
  select(-exclude)
  
cancer_data1 %>% 
  filter(tumor_seq_num != "1") %>% 
  # filter(primary_site_group_desc != "BLADDER",
  #        primary_site_group_desc != "BREAST",
  #        primary_site_group_desc != "MELANOMA OF SKIN",
  #        primary_site_group_desc != "LUNG/BRONCHUS-SMALL CELL") %>% 
  select(primary_site_group_desc, primary_site_desc, histology_desc) %>% 
  tbl_summary(by = primary_site_group_desc,
              sort = everything() ~ "frequency") %>% 
  modify_spanning_header(everything() ~ "primary_site_group_desc")
  
# For checking CDSC dx coding - all good - not needed anymore
# library(data.table)
# cancer_data_wide <- dcast(setDT(cancer_data1), mrn ~ rowid(mrn),
#                           value.var = c(
#                             "dx_dt",
#                             "primary_site_group_desc",
#                             "primary_site_desc",
#                             "histology_cd",
#                             "histology_desc"
#                           )) %>% 
#   select(mrn, ends_with("_1"), ends_with("_2"),
#          ends_with("_3"), ends_with("_4"),
#          ends_with("_5")) %>% 
#   arrange(primary_site_group_desc_1, primary_site_desc_1) %>% 
#   mutate(tmn_event = case_when(
#     primary_site_desc_2 == "BONE MARROW"                 ~ "diagnosis 2",
#     primary_site_desc_3 == "BONE MARROW"                 ~ "diagnosis 3",
#     primary_site_desc_4 == "BONE MARROW"                 ~ "diagnosis 4",
#     primary_site_desc_5 == "BONE MARROW"                 ~ "diagnosis 5",
#   ), .after = mrn) %>% 
#   mutate(tmn_dx_date = case_when(
#     tmn_event == "diagnosis 2"                   ~ dx_dt_2,
#     tmn_event == "diagnosis 3"                   ~ dx_dt_3,
#     tmn_event == "diagnosis 4"                   ~ dx_dt_4,
#     tmn_event == "diagnosis 5"                   ~ dx_dt_5
#   )) %>% 
#   full_join(clinical_data_1 %>% 
#               select(mrn, 
#                      CDSCcoded_non_myeloid_dx_dt = non_myeloid_dx_dt, 
#                      CDSCcoded_non_myeloid_histology_cd = non_myeloid_histology_cd,
#                      CDSCcoded_first_myeloid_dx_dt = first_myeloid_dx_dt), 
#             ., 
#             by = "mrn") %>% 
#   mutate(identical_primary = case_when(
#     dx_dt_1 == CDSCcoded_non_myeloid_dx_dt    ~ "TRUE",
#     dx_dt_1 != CDSCcoded_non_myeloid_dx_dt    ~ "FALSE"
#   ), .after = CDSCcoded_non_myeloid_histology_cd) %>% 
#   mutate(identical_tmn = case_when(
#     tmn_dx_date == CDSCcoded_first_myeloid_dx_dt    ~ "TRUE",
#     tmn_dx_date != CDSCcoded_first_myeloid_dx_dt    ~ "FALSE"
#   ), .after = CDSCcoded_first_myeloid_dx_dt)
# 
# write_csv(cancer_data_wide, 
#           paste0("data/processed_data",
#                  "/Diagnosis and TMN check from new diagnosis", 
#                  today(), ".csv"))
# write_csv(cancer_data_wide, 
#           paste0(path_clinical, "/ProcessedData",
#                  "/Diagnosis and TMN check from new diagnosis", 
#                  today(), ".csv"))



cancer_data1 <- cancer_data1 %>% 
  # The CDSC files doesn't include all info so need to extract them again
  inner_join(clinical_data_1 %>% 
              select(mrn, 
                     first_myeloid_dx_dt), 
            ., 
            by = "mrn") %>% 
  mutate(primary_cancer_site = case_when(
    tumor_seq_num == "1"                                  ~ primary_site_desc
  )) %>%
  mutate(primary_cancer_site_group = case_when(
    tumor_seq_num == "1"                                  ~ primary_site_group_desc
  )) %>%
  mutate(primary_cancer_histology_code = case_when(
    tumor_seq_num == "1"                                  ~ histology_cd
  )) %>%
  mutate(primary_cancer_histology = case_when(
    tumor_seq_num == "1"                                  ~ histology_desc
  )) %>%
  mutate(primary_cancer_dt = case_when(
    tumor_seq_num == "1"                                  ~ dx_dt
  )) %>%
  mutate(primary_cancer_grade_clin = case_when(
    tumor_seq_num == "1"                                  ~ grade_clinical_desc
  )) %>%
  mutate(primary_cancer_grade_path = case_when(
    tumor_seq_num == "1"                                  ~ grade_pathological_desc
  )) %>%
  mutate(primary_cancer_grade_tnm = case_when(
    tumor_seq_num == "1"                                  ~ stage_tnm_cs_mixed_group_desc
  )) %>%
  # mutate(is_tmn_dx = case_when(
  #   dx_dt == first_myeloid_dx_dt                                    ~ "Yes"
  # )) %>% 
  mutate(tmn_site = case_when(
    dx_dt == first_myeloid_dx_dt                                    ~ primary_site_desc
  )) %>% 
  mutate(tmn_site_group = case_when(
    dx_dt == first_myeloid_dx_dt                                    ~ primary_site_group_desc
  )) %>% 
  mutate(tmn_histology_code = case_when(
    dx_dt == first_myeloid_dx_dt                                    ~ histology_cd
  )) %>% 
  mutate(tmn_histology = case_when(
    dx_dt == first_myeloid_dx_dt                                    ~ histology_desc
  )) %>% 
  group_by(mrn) %>% 
  fill(starts_with("primary_cancer_"), starts_with("tmn_"), .direction = "updown") %>% 
  mutate(smoking_history = last(smoking_history)) %>% 
  ungroup() %>% 
  select(mrn, starts_with("primary_cancer_"), first_myeloid_dx_dt, starts_with("tmn_"),
         birth_dt : smoking_history) %>% 
  distinct()
  
  
  # # Code first TMN event
  # mutate(tmn_event = case_when(
  #   tumor_seq_num != 1 &
  #     primary_site_desc == "BONE MARROW"                 ~ "Yes"
  # ), .after = primary_site_desc) %>% 
  # group_by(mrn, tmn_event) %>%
  # mutate(tmn_seq_num = case_when(
  #   tumor_seq_num != 1 &
  #     primary_site_desc == "BONE MARROW"                ~ row_number()
  # ), .after = primary_site_desc) %>% 
  # # mutate(tmn_seq_num = case_when(
  # #   tumor_seq_num != 1 &
  # #     primary_site_desc == "BONE MARROW" &
  # #     primary_site_group_desc == "HEMERETIC"             ~ row_number()
  # # ), .after = primary_site_desc) %>% 
  # ungroup() %>% 
  # # Extract first TMN event info
  # mutate(tmn_tumor_id = case_when(
  #   tmn_seq_num == 1                                     ~ tumor_id
  # )) %>% 
  # mutate(tmn_dx_date = case_when(
  #   tmn_seq_num == 1                                     ~ dx_dt
  # )) %>% 
  # mutate(tmn_primary_site_desc = case_when(
  #   tmn_seq_num == 1                                     ~ primary_site_desc
  # )) %>% 
  # mutate(tmn_primary_site_group_desc = case_when(
  #   tmn_seq_num == 1                                     ~ primary_site_group_desc
  # )) %>% 
  # mutate(tmn_histology_desc = case_when(
  #   tmn_seq_num == 1                                     ~ histology_desc
  # )) %>% 
  # group_by(mrn) %>% 
  # fill(tmn_tumor_id, tmn_dx_date, 
  #      tmn_primary_site_group_desc,
  #      tmn_primary_site_desc, tmn_histology_desc, 
  #      .direction = "updown") %>% 
  # ungroup() %>% 
  # filter(tumor_seq_num == 1)

# write_csv(cancer_data1 %>% 
#   select(mrn, primary_site_desc, histology_desc, tmn_primary_site_desc, tmn_histology_desc), 
#   "data/processed_data/Patient primary cancer info and subsequent TMN info.csv")
# 
# a <- cancer_data1 %>% 
#   select(mrn, primary_site_desc, histology_desc, tmn_primary_site_desc, tmn_histology_desc)

# NGS ----
# molecular_ngs1 <- molecular_ngs %>% 
#   unite(ngs_type, c(ngs_gene, ngs_alteration), sep = " = ", remove = FALSE) %>% 
#   mutate(ngs_type = case_when(
#     !is.na(ngs_gene_alternative)        ~ ngs_alteration,
#     is.na(ngs_gene_alternative)         ~ ngs_type,
#   )) %>% 
#   select(mrn, ngs_type, ngs_significant) %>% 
#   distinct() %>%
#   pivot_wider(id_cols = mrn, names_from = ngs_type, values_from = ngs_significant#, values_fill = "na"
#               )
molecular_ngs <- molecular_ngs %>% 
  mutate(mrn = as.character(mrn))

molecular_ngs1 <- molecular_ngs %>% 
  select(mrn, ngs_report_date) %>% 
  mutate(ngs_report_date = as.POSIXct(format(ngs_report_date), format = "%Y-%m-%d")) %>%
  # mutate(ngs_report_date = as.character(ngs_report_date)) %>% 
  mutate(ngs_done = "ngs_dates") %>% 
  distinct() %>% 
  pivot_wider(id_cols = mrn, 
              names_from = ngs_done, 
              values_from = ngs_report_date, 
              values_fn = function(x) paste(x, collapse = ", "))




################################################################################# III ### Merge data
tmn_patients <- dna_data %>% 
  mutate(has_dna_sample = "Yes") %>% 
  # Merge with Cancer Char for patients with samples available
  full_join(., cancer_data1, by = c("mrn")) %>% 
  full_join(., molecular_ngs1, by = c("mrn")) %>% 
  # Merge with Demographic for patients with samples available
  inner_join(., clinical_data_1, 
             by = c("mrn", "first_myeloid_dx_dt")) %>% 
  select(#deidentified_patient_id, 
         mrn, has_dna_sample, 
         specimen_collection_dt, sample_id, sample_family_id,
         primary_cancer_treatment_type = non_myeloid_treatment_type,
         primary_cancer_treatment_start_dt = first_non_myeloid_treatment_start_dt,
         primary_cancer_tumor_id = non_myeloid_tumor_id, 
         primary_cancer_tumor_seq_num = non_myeloid_tumor_seq_num,
         starts_with("primary_cancer_"), 
         time_treatment_to_myeloid_dx,
         tmn_dx_dt = first_myeloid_dx_dt,
         starts_with("tmn_"),
         # tmn_tumor_id = myeloid_tumor_id, 
         everything()) %>% 

  # add new variables
  # create age
  mutate(age_at_diagnosis = round(interval(start = birth_dt, end = primary_cancer_dt)/
                                    duration(n = 1, units = "years"), 1)
  ) %>%
  mutate(age_at_sample = round(interval(start = birth_dt, end = specimen_collection_dt)/
                                 duration(n = 1, units = "years"), 1),
         .after = specimen_collection_dt
  ) %>% 
  mutate(time_primary_cancer_to_tmn_years = round(interval(start = primary_cancer_dt, end = tmn_dx_dt)/
                                 duration(n = 1, units = "years"), 1)
  ) %>% 
  # Survival
  # OS
  mutate(os_event = case_when(
    cdsc_vital_status == "ALIVE"            ~ 0,
    cdsc_vital_status == "DEAD"             ~ 1
  ), .before = cdsc_vital_status) %>% 
  mutate(os_time_from_dx_months = interval(start = primary_cancer_dt,
                                           end = last_contact_or_death_dt)/
           duration(n = 1, unit = "months"),
         .after = os_event) %>% 
  mutate(os_time_from_firsttx_months = interval(start = primary_cancer_treatment_start_dt,
                                           end = last_contact_or_death_dt)/
           duration(n = 1, unit = "months"),
         .after = os_event) %>% 
  
  mutate(across(c(where(is_character),
                  -sample_id, -sample_type, 
                  -sample_family_id, 
                  -primary_cancer_treatment_type), ~ str_to_sentence(.)))

write_rds(tmn_patients, 
          paste0("data/processed_data",
                 "/all_cancer_type_patients_with_tmn", 
                 str_remove_all(today(), "-"), ".rds"))


################################################################################# IV ### Find samples
tmn_patients <- 
  read_rds(paste0(here::here(), 
                  "/data/processed_data",
                  "/all_cancer_type_patients_with_tmn20251028.rds"))

tmn_blood <- tmn_patients %>% 
  mutate(blood_vs_tmn_sequence = case_when(
    specimen_collection_dt < tmn_dx_dt                 ~ "Sample before TMN",
    specimen_collection_dt == tmn_dx_dt                ~ "Sample at TMN",
    specimen_collection_dt > tmn_dx_dt                 ~ "Sample after TMN",
    is.na(tmn_dx_dt)                                   ~ "No TMN",
    TRUE                                               ~ NA_character_
  ), .after = sample_family_id) %>% 
  group_by(mrn) %>% 
  # mutate(number_of_sample_bf_tmn = case_when(
  #   blood_vs_tmn_sequence == "Yes"                                ~ n()
  # ), .after = blood_vs_tmn_sequence) %>% 
  mutate(sample_sequence_number = row_number(),
         .after = blood_vs_tmn_sequence) %>% 
  ungroup() %>% 
  arrange(mrn, specimen_collection_dt, blood_vs_tmn_sequence) %>% 
  group_by(mrn) %>%
  mutate(sample_lag_days = as.Date(specimen_collection_dt) - lag(as.Date(specimen_collection_dt)),
         sample_lag_days = as.numeric(str_remove(sample_lag_days, " days")),
         .after = sample_sequence_number
  ) %>% 
  ungroup() %>% 
  mutate(interval_sample_to_tmn_days = interval(start = tmn_dx_dt, end = specimen_collection_dt) /
           duration(n = 1, units = "days"),
         .after = sample_lag_days) %>%
  mutate(interval_primary_cancer_to_samples_days = interval(start = primary_cancer_dt, end = specimen_collection_dt) /
           duration(n = 1, units = "days"),
         .after = sample_lag_days) %>% 
  mutate(has_dna_sample = case_when(
    has_dna_sample == "Yes"             ~ "Yes",
    is.na(has_dna_sample)               ~ "No"
  ))

write_csv(tmn_blood, 
          paste0("data/processed_data",
                 "/CHevolution_all_cancer_type_patients_with_tmn_", 
                 str_remove_all(today(), "-"), ".csv"))
write_csv(tmn_blood, 
          paste0(path_clinical, "/ProcessedData",
                 "/CHevolution_all_cancer_type_patients_with_tmn_", 
                 str_remove_all(today(), "-"), ".csv"))
write_rds(tmn_blood, 
          paste0("data/processed_data",
                 "/CHevolution_all_cancer_type_patients_with_tmn_", 
                 str_remove_all(today(), "-"), ".rds"))



################################################################################# V ### Treatments
emr_medication <- emr_medication %>% 
  mutate(mrn = as.character(mrn))

registry_treatment <- registry_treatment %>% 
  mutate(mrn = as.character(mrn))


################################################################################# VI ### Samples list
tmn_blood <-
  read_rds(paste0(here::here(), 
                  "/data/processed_data",
                  "/CHevolution_all_cancer_type_patients_with_tmn_20251028.rds"))

breast_samples <- tmn_blood %>%
  filter(primary_cancer_site_group == "Breast") #%>% 
  # filter(interval_sample_to_tmn_days <= 10)

write_csv(breast_samples %>% select(mrn, sample_sequence_number, sample_id,
                               sample_family_id, interval_sample_to_tmn_days,
                               tmn_dx_dt,
                               specimen_collection_dt, blood_vs_tmn_sequence,
                               best_anatomic_site, sample_type) %>% 
            mutate(tmn_dx_dt = as.Date(tmn_dx_dt)) %>% 
            mutate(specimen_collection_dt = as.Date(specimen_collection_dt)),
          paste0("data/processed_data",
                 "/CHevolution_breast_patients_sampleslist_", 
                 str_remove_all(today(), "-"), ".csv"))
write_csv(breast_samples %>% select(mrn, sample_sequence_number, sample_id,
                                    sample_family_id, interval_sample_to_tmn_days,
                                    tmn_dx_dt,
                                    specimen_collection_dt, blood_vs_tmn_sequence, 
                                    best_anatomic_site, sample_type) %>% 
            mutate(tmn_dx_dt = as.Date(tmn_dx_dt)) %>% 
            mutate(specimen_collection_dt = as.Date(specimen_collection_dt)),
          paste0(path_clinical, "/ProcessedData",
                 "/CHevolution_breast_patients_sampleslist_", 
                 str_remove_all(today(), "-"), ".csv"))
















## Import library
library(tidyverse)
library(janitor)
# library(gtsummary)
# theme_set(theme_classic())
# theme_gtsummary_compact()

###################################################################### I ### Load data
path_clinical <- fs::path("", "Volumes", "Gillis_Research", "Lab_Data", "CHEvolutionTMN")
clinical_data <-
  readxl::read_xlsx(paste0(here::here(), #path_clinical, 
                           "/data/raw data/CHevolution_CohortRaw_10R25000313_v6_20251110.xlsx"), 
                    # "/RawData/CHevolution_CohortRaw_10R25000313_v6_20251110.xlsx"), 
                    sheet = "Base_Data") %>% 
  clean_names()
cancer_data <-
  readxl::read_xlsx(paste0(here::here(), #path_clinical, 
                           "/data/raw data/CHevolution_CohortRaw_10R25000313_v6_20251110.xlsx"), 
                    # "/RawData/CHevolution_CohortRaw_10R25000313_v6_20251110.xlsx"), 
                    sheet = "Registry_General") %>% 
  clean_names()
registry_treatment <-
  readxl::read_xlsx(paste0(here::here(), #path_clinical, 
                           "/data/raw data/CHevolution_CohortRaw_10R25000313_v6_20251110.xlsx"), 
                    # "/RawData/CHevolution_CohortRaw_10R25000313_v6_20251110.xlsx"), 
                    sheet = "Registry_Treatment") %>% 
  clean_names()
emr_medication <-
  readxl::read_xlsx(paste0(here::here(), #path_clinical, 
                           "/data/raw data/CHevolution_CohortRaw_10R25000313_v6_20251110.xlsx"), 
                    # "/RawData/CHevolution_CohortRaw_10R25000313_v6_20251110.xlsx"), 
                    sheet = "Emr_Medication") %>% 
  clean_names()
dna_data <-
  readxl::read_xlsx(paste0(here::here(), #path_clinical, 
                           "/data/raw data/CHevolution_CohortRaw_10R25000313_v6_20251110.xlsx"), 
                    # "/RawData/CHevolution_CohortRaw_10R25000313_v6_20251110.xlsx"), 
                    sheet = "Biobanking") %>% 
  clean_names()
molecular_ngs <-
  readxl::read_xlsx(paste0(here::here(), #path_clinical, 
                           "/data/raw data/CHevolution_CohortRaw_10R25000313_v6_20251110.xlsx"), 
                    # "/RawData/CHevolution_CohortRaw_10R25000313_v6_20251110.xlsx"), 
                    sheet = "Molecular_NGS") %>% 
  clean_names()

tmn_blood <-
  read_rds(paste0(here::here(), 
                  "/data/processed_data",
                  "/CHevolution_all_cancer_type_patients_with_tmn_20251028.rds"))
# tmn_blood <-
#   read_csv(paste0(path_clinical, 
#                   "/ProcessedData/CHevolution_AllPtswTMN_20251028.csv"))
samples_received <-
  readxl::read_xlsx(paste0(path_clinical, 
                    "/CoreRequests/CHevolution_GillisWest_TC_SamplesReceived_EXP2491DNA_20251110.xlsx"), 
                    skip = 9) %>% 
  clean_names() %>% 
  select(-x1) %>% 
  purrr::keep(~!all(is.na(.)))


###################################################################### II ### Data claening
dna_data <- dna_data %>% 
  mutate(mrn = as.character(mrn)) %>% 
  # Include only sample available for research
  filter(sample_available_for_research_indicator == "TRUE") %>% 
  # Select sample type needed for CH detection
  filter(sample_type %in% c(
    "Buffy Coat", "Genomic DNA","MNC", "PBMC",
    "MNC less CD138+", "Unprocessed Liquid Tissue")
  ) %>% 
  select(mrn, sample_family_id, sample_id_ordered = sample_id,
         specimen_collection_dt = collection_dt,
         best_anatomic_site, cdsc_sample_type = sample_type) %>%
  full_join(., samples_received %>% 
               rename(chevolution_dnacore_studyname = study,
                      sample_id_received = sample), 
             by = c("mrn", "sample_family_id", 
                    "sample_id_ordered" = "dna_master_lv_id")) %>% 
  group_by(mrn) %>% 
  fill(chevolution_dnacore_studyname, lv_alias_dna_id, .direction = "updown") %>% 
  ungroup() %>% 
  # add same sample/same date on the same row with Blood as first listed in string
  arrange(mrn, specimen_collection_dt, best_anatomic_site) %>% 
  # Summarize to have 1 sample/day per row and not 1 row for each aliquot of the same sample 
  group_by(mrn, sample_family_id, specimen_collection_dt) %>% 
  summarise_at(vars(sample_id_ordered, best_anatomic_site, cdsc_sample_type, 
                    chevolution_dnacore_studyname, lv_alias_dna_id), str_c, collapse = "; ") %>%
  arrange(mrn, specimen_collection_dt, best_anatomic_site) %>% 
  group_by(mrn, specimen_collection_dt, chevolution_dnacore_studyname, lv_alias_dna_id) %>% 
  summarise_at(vars(sample_family_id, sample_id_ordered, best_anatomic_site, cdsc_sample_type), str_c, collapse = "; ") %>%
  ungroup()

set.seed(123)
clinical_data_1 <- clinical_data %>% 
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
#                  "/CHevolution_MatchMRN_DeidentifiedIds_for_CHEvolution_project_before_filtering",
#                  str_remove_all(today(), "-"), ".csv"))
# write_csv(clinical_data_1 %>%
#             select(mrn, deidentified_patient_id, random_id),
#           paste0(path_clinical, "/ProcessedData",
#                  "/CHevolution_MatchMRN_DeidentifiedIds_for_CHEvolution_project_before_filtering",
#                  str_remove_all(today(), "-"), ".csv"))


# The clinical has multiple row per patients due to different samples included in the data
# and to multiple myeloid dx
# I am coding the myeloid dx myself
clinical_data_1 <- clinical_data_1 %>%
  # Clean up var created for coding de-id ids
  select(c(deidentified_patient_id, mrn, everything(), -c(zero, study_name, add_zero))) %>% 
  # CDSC made a clinical data including the samples (some that are not useful to us)
  # To get unique row clinical data
  # remove variables related to the samples
  select(deidentified_patient_id, 
    mrn, 
    cdsc_vital_status, last_contact_or_death_dt, last_visit_dt, 
    # non-myeloid info
    non_myeloid_tumor_id, non_myeloid_tumor_seq_num, 
    non_myeloid_treatment_type,
    first_non_myeloid_treatment_start_dt, 
    # myeloid info
    # myeloid_tumor_id, myeloid_tumor_seq_num, 
    first_myeloid_dx_dt, 
    time_treatment_to_myeloid_dx = treatment_to_myeloid_dx
  ) %>% 
  distinct() %>% 
  arrange(mrn) %>% 
  mutate(non_myeloid_treatment_type = str_to_lower(non_myeloid_treatment_type)) %>% 
  group_by(across(c(-non_myeloid_treatment_type))) %>% 
  summarise_at(vars(non_myeloid_treatment_type), str_c, collapse = "/") %>%
  ungroup()

cancer_data <- cancer_data %>% 
  mutate(mrn = as.character(mrn)) %>% 
  arrange(mrn, dx_dt) %>% 
  select(-starts_with("patient"))

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

molecular_ngs <- molecular_ngs %>% 
  mutate(mrn = as.character(mrn))

molecular_ngs1 <- molecular_ngs %>% 
  mutate(ngs_report_date = as.POSIXct(format(ngs_report_date), format = "%Y-%m-%d")) %>%
  mutate(ngs_gene = case_when(
    !is.na(ngs_gene_alternative)                   ~ ngs_alteration,
    !is.na(ngs_gene)                               ~ ngs_gene
  )) %>%
  select(mrn, ngs_report_date, ngs_gene, ngs_alteration, ngs_significance = ngs_significant) %>% 
  arrange(mrn, ngs_report_date) %>% 
  group_by(mrn, ngs_gene) %>% 
  mutate(gene_test_seq_num = row_number()) %>% 
  # mutate(ngs_report_date = as.character(ngs_report_date)) %>% 
  # mutate(ngs_done = "ngs_dates") %>% 
  distinct() %>% 
  pivot_wider(id_cols = mrn, 
              names_from = c(ngs_gene, gene_test_seq_num), 
              values_from = c(ngs_report_date, ngs_significance, ngs_alteration), 
              # values_fn = function(x) paste(x, collapse = ", ")
              names_vary = "slowest",
              names_glue = "{ngs_gene}_{str_c(.value, gene_test_seq_num)}"
              )


################################################################################# III ### Merge data
tmn_patients <- dna_data %>% 
  mutate(has_dna_sample = "Yes") %>% 
  # Merge with Cancer Char for patients with samples available
  full_join(., cancer_data1, by = c("mrn")) %>% 
  full_join(., molecular_ngs1, by = c("mrn")) %>% 
  # Merge with Demographic for patients with samples available
  inner_join(., clinical_data_1, 
             by = c("mrn", "first_myeloid_dx_dt")) %>% 
  mutate(has_dna_sample = case_when(
    is.na(chevolution_dnacore_studyname) &
      has_dna_sample == "Yes"            ~ "Not available",
    has_dna_sample == "Yes"              ~ "Yes",
    TRUE                                 ~ "No"
  )) %>% 
  mutate(across(specimen_collection_dt, ~ case_when(
    has_dna_sample == "Not available"   ~ NA_Date_,
    TRUE                                ~ specimen_collection_dt
  )))

tmn_patients <- tmn_patients %>% 
  select(deidentified_patient_id, 
    mrn, has_dna_sample, 
    specimen_collection_dt, sample_id_ordered, sample_family_id,
    chevolution_dnacore_studyname, lv_alias_dna_id,
    sample_id_ordered, best_anatomic_site, cdsc_sample_type,
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
  
  mutate(across(c(where(is_character), -deidentified_patient_id,
                  -sample_id_ordered, -cdsc_sample_type, 
                  -sample_family_id, 
                  -primary_cancer_treatment_type), ~ str_to_sentence(.)))

write_rds(tmn_patients, 
          paste0("data/processed_data",
                 "/all_cancer_type_patients_with_tmn_", 
                 str_remove_all(today(), "-"), ".rds"))


################################################################################# IV ### Find samples
























tmn_patients <- 
  read_rds(paste0(here::here(), 
                  "/data/processed_data",
                  "/all_cancer_type_patients_with_tmn_20260107.rds"))

tmn_blood <- tmn_patients %>% 
  mutate(blood_vs_tmn_sequence = case_when(
    specimen_collection_dt < tmn_dx_dt                 ~ "Sample before TMN",
    specimen_collection_dt == tmn_dx_dt                ~ "Sample at TMN",
    specimen_collection_dt > tmn_dx_dt                 ~ "Sample after TMN",
    is.na(tmn_dx_dt)                                   ~ "No TMN",
    is.na(specimen_collection_dt)                      ~ "No sample",
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
  mutate(interval_tmn_to_sample_days = interval(start = tmn_dx_dt, end = specimen_collection_dt) /
           duration(n = 1, units = "days"),
         .after = sample_lag_days) %>%
  mutate(interval_primary_cancer_to_sample_days = interval(start = primary_cancer_dt, end = specimen_collection_dt) /
           duration(n = 1, units = "days"),
         .after = sample_lag_days)

write_csv(tmn_blood, 
          paste0("data/processed_data",
                 "/CHevolution_AllPtswTMN_", 
                 str_remove_all(today(), "-"), ".csv"))
write_csv(tmn_blood, 
          paste0(path_clinical, "/ProcessedData",
                 "/CHevolution_AllPtswTMN_", 
                 str_remove_all(today(), "-"), ".csv"))
write_rds(tmn_blood, 
          paste0("data/processed_data",
                 "/CHevolution_AllPtswTMN_", 
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
                  "/CHevolution_AllPtswTMN_20260107.rds"))

breast_samples <- tmn_blood %>%
  filter(primary_cancer_site_group == "Breast") #%>% 
# filter(interval_sample_to_tmn_days <= 10)

write_csv(breast_samples,# %>% select(mrn, sample_sequence_number, sample_id_ordered,
            #                         sample_family_id, interval_sample_to_tmn_days,
            #                         tmn_dx_dt,
            #                         specimen_collection_dt, blood_vs_tmn_sequence,
            #                         best_anatomic_site, sample_type) %>% 
            # mutate(tmn_dx_dt = as.Date(tmn_dx_dt)) %>% 
            # mutate(specimen_collection_dt = as.Date(specimen_collection_dt)),
          paste0("data/processed_data",
                 "/CHevolution_BreastSamples_", 
                 str_remove_all(today(), "-"), ".csv"))
write_rds(breast_samples,# %>% select(mrn, sample_sequence_number, sample_id_ordered,
          #                         sample_family_id, interval_sample_to_tmn_days,
          #                         tmn_dx_dt,
          #                         specimen_collection_dt, blood_vs_tmn_sequence,
          #                         best_anatomic_site, sample_type) %>% 
          # mutate(tmn_dx_dt = as.Date(tmn_dx_dt)) %>% 
          # mutate(specimen_collection_dt = as.Date(specimen_collection_dt)),
          paste0("data/processed_data",
                 "/CHevolution_BreastSamples_", 
                 str_remove_all(today(), "-"), ".rds"))
write_csv(breast_samples,# %>% select(mrn, sample_sequence_number, sample_id_ordered,
          #                         sample_family_id, interval_sample_to_tmn_days,
          #                         tmn_dx_dt,
          #                         specimen_collection_dt, blood_vs_tmn_sequence,
          #                         best_anatomic_site, sample_type) %>% 
          # mutate(tmn_dx_dt = as.Date(tmn_dx_dt)) %>% 
          # mutate(specimen_collection_dt = as.Date(specimen_collection_dt)),
          paste0(path_clinical, "/ProcessedData",
                 "/CHevolution_BreastSamples_", 
                 str_remove_all(today(), "-"), ".csv"))



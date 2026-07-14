## Import library
library(tidyverse)
library(janitor)
# library(gtsummary)


###################################################################### I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research", "Lab_Data", "CHEvolutionTMN")

data <-
  read_rds(paste0(here::here(), 
                  # path, "/ProcessedData", "/CHevolution_BreastSamples_20260120.csv"
                  "/data/processed_data",
                  "/CHevolution_BreastSamples_20260120.rds"
                  ))

treatment_dose <-
  readxl::read_xlsx(paste0(here::here(), 
                           # path,
                           # "/ProcessedData/CHEvolution_Treatment_20260726.xlsx"),
                           "/data/processed_data/CHEvolution_Treatment_20260714.xlsx"),
                    sheet = "w doses kg m2", na = "NA") %>% 
  clean_names()


parent_dir_path <- dirname(path)
drug_class <- 
  read.csv(paste0(#here::here(), 
    # "/data/processed_data",
    parent_dir_path,
    "/SharedResources/BoltonDrugCategories",
    "/CHevolution_Updated_BoltonChemoDosing_20260713.csv"))


###################################################################### II ### Cumulative drug doses----
treatment_dose <- treatment_dose %>% 
  select(-c(contains("bsa_")))


treatment_dose1 <- 
  treatment_dose %>% 
  select(-c(drug2, dose2, unit2, drug3, dose3, unit3, rad_dose : rad_site)) %>%
  filter(class != "Radiation") %>% 
  rename(drug = drug1, dose = dose1, unit = unit1) %>% 
  bind_rows(., treatment_dose %>% 
              select(mrn, breast_dx : wt_kg, drug2, dose2, unit2, class : sample7_collection_dt) %>% 
              filter(!is.na(drug2)) %>% 
              rename(drug = drug2, dose = dose2, unit = unit2)) %>% 
  bind_rows(., treatment_dose %>% 
              select(mrn, breast_dx : wt_kg, drug3, dose3, unit3, class : sample7_collection_dt) %>% 
              filter(!is.na(drug3)) %>% 
              rename(drug = drug3, dose = dose3, unit = unit3)) %>% 
  mutate(drug = case_when(
    drug == "nab-paclitaxel"             ~ "paclitaxel",
    TRUE                                 ~ drug
  )) %>% 
  mutate(weight_adjusted_dose = dose / wt_kg, .after = unit) %>% 
  # Add Bolton categories
  left_join(., drug_class, 
            by = c("drug" = "drug_name")) %>% 
  mutate(last_sample_dt = coalesce(sample7_collection_dt,
                                   sample6_collection_dt,
                                   sample5_collection_dt,
                                   sample4_collection_dt,
                                   sample3_collection_dt,
                                   sample2_collection_dt,
                                   sample1_collection_dt
                                   ), .before = sample1_collection_dt) %>% 
  # Exclude drugs that ends before sample1 or start after last
  # mutate(new_dt = tx_end + months(6), .after = tx_end)
  mutate(drug_before_sample1 = case_when(
    tx_end + months(6) <= sample1_collection_dt        ~ "Yes"
  )) %>% 
  filter(is.na(drug_before_sample1)) %>% 
  mutate(drug_after_last = case_when(
    tx_start > last_sample_dt                          ~ "Yes"
  )) %>% 
  filter(is.na(drug_after_last)) %>% 
  select(-drug_before_sample1, -drug_after_last) %>% 
  
  # Approach 2
  mutate(drug_tertile_sample1_to_last = case_when(
    n() > 2                                             ~ ntile(weight_adjusted_dose, 3),
  )) %>% 
  group_by(mrn, drug) %>% 
  mutate(munber_of_time_received_by_patient = n()) %>% 
  group_by(drug) %>% 
  mutate(munber_of_time_received_in_cohort = n()) %>% 
  ungroup() %>% 
  select(mrn, weight_adjusted_dose, drug, drug_tertile_sample1_to_last, 
         munber_of_time_received_by_patient, munber_of_time_received_in_cohort)
  
  write_csv(treatment_dose1, "First tertile calculated over all line received from sample 1 to last sample.csv")
  
  
  
  # # Approach 1
  # # categories drug vs xx sample
  # # mutate(drug_before_first_sample = case_when(
  # #   tx_start <= sample1_collection_dt                    ~ "Yes"
  # # ), .after = drug) %>% 
  # mutate(drugs_between_sample1_to_last = case_when(
  #   tx_end >= sample1_collection_dt &
  #     tx_start <= last_sample_dt                       ~ "Yes"
  # ), .before = last_sample_dt) %>%
  # # calculate cumulative dose per patient per drug
  # group_by(mrn, drug, drugs_between_sample1_to_last) %>%
  #   
  # mutate(drug_cumulative_dose_sample1_to_last = case_when(
  #   drugs_between_sample1_to_last == "Yes"             ~ sum(weight_adjusted_dose)
  # ), .after = drugs_between_sample1_to_last) %>% 
  # # group_by(mrn, chemotherapy_drug, pre_to_last_sample) %>% 
  # # mutate(cumulative_dose_pre_lastseqsample = case_when(
  # #   pre_to_last_sample == "Yes"                         ~ sum(weight_adjusted_dose)
  # # )) %>% 
  # ungroup() %>% 
  # 
  # # Remove rows that were accumulated
  # select(mrn, drug, narrow_drug_class_cytotoxic_only,
  #        drug_cumulative_dose_sample1_to_last) %>%
  # arrange(mrn, drug,
  #         drug_cumulative_dose_sample1_to_last) %>%
  # distinct(mrn, drug, .keep_all = TRUE) %>%
  # # create tertile categories
  # group_by(drug) %>% 
  # mutate(n = n()) %>%
  # mutate(tertile_drug_dose_sample1_to_last = case_when(
  #   n() > 2                                            ~ ntile(drug_cumulative_dose_sample1_to_last, 3)
  # )) %>% 
  
  # Approach 2
  
  
  
  
  
  
    mutate(tertile_pre_lastseqsample = case_when(
      n() > 2                                             ~ ntile(cumulative_dose_pre_lastseqsample, 3),
      n() == 1 &
        chemotherapy_drug == "abraxane"                   ~ 1,
      n() == 2 &
        chemotherapy_drug %in% 
        c("epirubicin", "fluorouracil")                   ~ 2
    )) %>% 
    # mutate(tertile_cum_drug_cat = case_when(
    #   n() > 1 &
    #     ntile(cumulative_dose, 3) == 1                    ~ "Low",
    #   n() > 1 &
    #     ntile(cumulative_dose, 3) == 2                    ~ "Medium",
    #   n() > 1 &
    #     ntile(cumulative_dose, 3) == 3                    ~ "High",
    # )) %>% 
    
    # calculate score - sum scores for each drug in a specific drug class
    group_by(mrn, narrow_drug_class_cytotoxic_only) %>% 
    # mutate(n = n()) %>%
    mutate(drug_class_score_pre_firstseqsample = case_when(
      !is.na(tertile_pre_firstseqsample)                  ~ sum(tertile_pre_firstseqsample, na.rm = TRUE)
    )) %>% 
    mutate(drug_class_score_pre_lastseqsample = case_when(
      !is.na(tertile_pre_lastseqsample)                   ~ sum(tertile_pre_lastseqsample, na.rm = TRUE)
    )) %>% 
    distinct(mrn, narrow_drug_class_cytotoxic_only, .keep_all = TRUE) %>% 
    select(-chemotherapy_drug) %>% 
    # group_by(narrow_drug_class_cytotoxic_only) %>%
    # mutate(n = n()) %>% 
    # Bolton summed score
    group_by(mrn, narrow_drug_class_cytotoxic_only) %>% 
    mutate(sum_bolton_drug_score_pre_firstseqsample_per_class_patient = case_when(
      !is.na(drug_class_score_pre_firstseqsample)         ~ sum(drug_class_score_pre_firstseqsample, na.rm = TRUE)
    )) %>% 
    mutate(sum_bolton_drug_score_pre_lastseqsample_per_class_patient = case_when(
      !is.na(drug_class_score_pre_lastseqsample)          ~ sum(drug_class_score_pre_lastseqsample, na.rm = TRUE)
    )) %>% 
    # Patient's
    group_by(mrn) %>% 
    mutate(sum_all_drug_score_pre_firstseqsample_per_patient = case_when(
      !is.na(drug_class_score_pre_firstseqsample)         ~ sum(drug_class_score_pre_firstseqsample, na.rm = TRUE)
    )) %>% 
    mutate(sum_all_drug_score_pre_lastseqsample_per_patient = case_when(
      !is.na(drug_class_score_pre_lastseqsample)          ~ sum(drug_class_score_pre_lastseqsample, na.rm = TRUE)
    )) %>% 
    ungroup() %>% 
    # # use cut() instead of ntile() to attribute the ties in values
    mutate(tertile_patient_pre_firstseqsample = cut(sum_all_drug_score_pre_firstseqsample_per_patient,
                                                    breaks = quantile(sum_all_drug_score_pre_firstseqsample_per_patient, 
                                                                      probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                                                    include.lowest = TRUE,
                                                    labels = c("Low", "Medium", "High"))) %>%
    mutate(tertile_patient_pre_lastseqsample = cut(sum_all_drug_score_pre_lastseqsample_per_patient,
                                                   breaks = quantile(sum_all_drug_score_pre_lastseqsample_per_patient, 
                                                                     probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                                                   include.lowest = TRUE,
                                                   labels = c("Low", "Medium", "High"))) %>% 
    full_join(sequenced_patient_data %>% 
                select(mrn, deidentified_patient_id) %>% 
                distinct() %>% 
                mutate(mrn = as.character(mrn)), 
              ., 
              by = "mrn") %>% 
    select(mrn, deidentified_patient_id, narrow_drug_class_cytotoxic_only, 
           sum_bolton_drug_score_pre_firstseqsample_per_class_patient, 
           sum_bolton_drug_score_pre_lastseqsample_per_class_patient,
           sum_all_drug_score_pre_firstseqsample_per_patient,
           tertile_patient_pre_firstseqsample,
           sum_all_drug_score_pre_lastseqsample_per_patient,
           tertile_patient_pre_lastseqsample
    ) %>% 
    filter(!is.na(narrow_drug_class_cytotoxic_only))


radiation_dose1 <- 
  treatment_dose %>% 
  select(-c(drug1, dose1, unit1, drug2, dose2, unit2, drug3, dose3, unit3)) %>%
  filter(!is.na(rad_dose)) %>% 
  mutate(last_sample_dt = coalesce(sample7_collection_dt,
                                   sample6_collection_dt,
                                   sample5_collection_dt,
                                   sample4_collection_dt,
                                   sample3_collection_dt,
                                   sample2_collection_dt,
                                   sample1_collection_dt
  ), .before = sample1_collection_dt) %>% 
  # Exclude drugs that ends before sample1 or start after last
  mutate(rad_before_sample1 = case_when(
    tx_end <= sample1_collection_dt                    ~ "Yes"
  ), .before = sample1_collection_dt) %>% 
  # filter(is.na(rad_before_sample1)) %>% 
  mutate(rad_after_last = case_when(
    tx_start > last_sample_dt                          ~ "Yes"
  ), .before = sample1_collection_dt) %>% 
  filter(is.na(rad_after_last)) %>% 
  select(-rad_before_sample1, -rad_after_last) %>% 
  
  
  
  
  
  # select(mrn, deidentified_patient_id, treatment_type, x2gy_fx_eqd2, radiation_body_site) %>% 
  # # Marked patients who have no data
  # mutate(has_dosing_data = case_when(
  #   !is.na(x2gy_fx_eqd2)                                ~ "Yes"
  # ), .after = x2gy_fx_eqd2) %>%
  # group_by(mrn) %>% 
  # fill(has_dosing_data, .direction = "updown") %>% 
  # mutate(has_dosing_data = replace_na(has_dosing_data, "No")) %>% 
  # mutate(sample_sequence = row_number(), .after = mrn) %>% 
  # ungroup() %>% 
  # 
  # mutate(pre_to_first_sample = case_when(
  #   sample_sequence == 1                                ~ "Yes"
  # )) %>% 
  # mutate(pre_to_last_sample = case_when(
  #   sample_sequence >= 1                                ~ "Yes"
  # )) %>% 
  # mutate(is_breast_radiation = case_when(
  #   str_detect(radiation_body_site, "breast") |
  #     str_detect(radiation_body_site, "CW")             ~ "Yes"
  # )) %>% 
  # # calculate cumulative dose per patient
  # group_by(mrn, pre_to_first_sample) %>% 
  # mutate(cumulative_dose_pre_firstseqsample = case_when(
  #   pre_to_first_sample == "Yes" &
  #     has_dosing_data == "Yes"                          ~ x2gy_fx_eqd2
  # )) %>% 
  # group_by(mrn, pre_to_last_sample) %>% 
  # mutate(cumulative_dose_pre_lastseqsample = case_when(
  #   pre_to_last_sample == "Yes" &
  #     has_dosing_data == "Yes"                          ~ sum(x2gy_fx_eqd2, na.rm = TRUE)
  # )) %>% 
  # ungroup() %>% 
  # 
  # # calculate cumulative dose per patient with breast site radiation only
  # group_by(mrn, is_breast_radiation, pre_to_first_sample) %>% 
  # mutate(cumulative_dose_pre_firstseqsample_breast_site_only = case_when(
  #   pre_to_first_sample == "Yes" &
  #     is_breast_radiation == "Yes" &
  #     has_dosing_data == "Yes"                          ~ sum(x2gy_fx_eqd2, na.rm = TRUE)
  # )) %>% 
  # group_by(mrn, is_breast_radiation, pre_to_last_sample) %>% 
  # mutate(cumulative_dose_pre_lastseqsample_breast_site_only = case_when(
  #   pre_to_last_sample == "Yes" &
  #     is_breast_radiation == "Yes" &
  #     has_dosing_data == "Yes"                          ~ sum(x2gy_fx_eqd2, na.rm = TRUE)
  # )) %>% 
  # ungroup() %>% 
  # select(-c(x2gy_fx_eqd2, sample_sequence, 
  #           pre_to_first_sample, pre_to_last_sample, 
  #           is_breast_radiation)) %>%
  # group_by(mrn) %>% 
  # fill(cumulative_dose_pre_firstseqsample, cumulative_dose_pre_lastseqsample,
  #      cumulative_dose_pre_firstseqsample_breast_site_only,
  #      cumulative_dose_pre_lastseqsample_breast_site_only,
  #      radiation_body_site, .direction = "updown") %>% 
  # ungroup() %>% 
  # distinct(mrn, .keep_all = TRUE) %>% 
  # 
  # # create tertile categories
  # mutate(tertile_score_pre_firstseqsample = cut(cumulative_dose_pre_firstseqsample,
  #                                               breaks = quantile(cumulative_dose_pre_firstseqsample, 
  #                                                                 probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
  #                                               include.lowest = TRUE,
  #                                               labels = c("1", "2", "3")), 
  #        .after = cumulative_dose_pre_firstseqsample) %>% 
  # mutate(tertile_score_pre_lastseqsample = cut(cumulative_dose_pre_lastseqsample,
  #                                              breaks = quantile(cumulative_dose_pre_lastseqsample, 
  #                                                                probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
  #                                              include.lowest = TRUE,
  #                                              labels = c("1", "2", "3")), 
  #        .after = cumulative_dose_pre_lastseqsample) %>% 
  # 
  # mutate(tertile_score_pre_firstseqsample_breast_site_only = cut(cumulative_dose_pre_firstseqsample_breast_site_only,
  #                                                                breaks = quantile(cumulative_dose_pre_firstseqsample_breast_site_only, 
  #                                                                                  probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
  #                                                                include.lowest = TRUE,
  #                                                                labels = c("1", "2", "3")), 
  #        .after = cumulative_dose_pre_firstseqsample_breast_site_only) %>% 
  # mutate(tertile_score_pre_lastseqsample_breast_site_only = cut(cumulative_dose_pre_lastseqsample_breast_site_only,
  #                                                               breaks = quantile(cumulative_dose_pre_lastseqsample_breast_site_only, 
  #                                                                                 probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
  #                                                               include.lowest = TRUE,
  #                                                               labels = c("1", "2", "3")), 
  #        .after = cumulative_dose_pre_lastseqsample_breast_site_only)




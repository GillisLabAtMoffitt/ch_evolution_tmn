## Import library
library(tidyverse)
library(janitor)
# library(gtsummary)


###################################################################### I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research", "Lab_Data", "CHEvolutionTMN")

# data <-
#   read_rds(paste0(here::here(), 
#                   # path, "/ProcessedData", "/CHevolution_BreastSamples_20260120.csv"
#                   "/data/processed_data",
#                   "/CHevolution_BreastSamples_20260120.rds"
#                   ))

treatment_dose <-
  readxl::read_xlsx(paste0(#here::here(), 
                           path,
                           "/ProcessedData/CHEvolution_Treatment_20260721.xlsx"),
                           # "/data/processed_data/CHEvolution_Treatment_20260721.xlsx"),
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
  select(-c(contains("bsa_"))) %>% 
  mutate(mrn = as.character(mrn))


treatment_dose1 <- 
  treatment_dose %>% 
  select(-c(drug2, dose2, unit2, raw_dose2, tertile2,
            drug3, dose3, unit3, raw_dose3, tertile3,
            rad_dose : rad_site)) %>%
  filter(class != "Radiation") %>% 
  rename(drug = drug1, dose = dose1, unit = unit1, raw_dose = raw_dose1, tertile = tertile1) %>% 
  bind_rows(., treatment_dose %>% 
              select(mrn, breast_dx : wt_kg, drug2, dose2, unit2, raw_dose2, tertile2, class : sample7_collection_dt) %>% 
              filter(!is.na(drug2)) %>% 
              rename(drug = drug2, dose = dose2, unit = unit2, raw_dose = raw_dose2, tertile = tertile2) %>% 
              mutate(tertile = as.numeric(tertile))) %>% 
  bind_rows(., treatment_dose %>% 
              select(mrn, breast_dx : wt_kg, drug3, dose3, unit3, raw_dose3, tertile3, class : sample7_collection_dt) %>% 
              filter(!is.na(drug3)) %>% 
              rename(drug = drug3, dose = dose3, unit = unit3, raw_dose = raw_dose3, tertile = tertile3)) %>% 
  mutate(weight_adjusted_dose = dose / wt_kg, .after = unit) %>% 
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
  group_by(drug) %>% 
  mutate(drug_tertile_sample1_to_last = case_when(
    n() > 2                                             ~ ntile(weight_adjusted_dose, 3),
  ), .after = tertile) %>% 
  mutate(drug_tertile_sample1_to_last = coalesce(tertile, drug_tertile_sample1_to_last)) %>% 
  # group_by(mrn, drug) %>% 
  # mutate(munber_of_time_received_by_patient = n()) %>% 
  # group_by(drug) %>% 
  # mutate(munber_of_time_received_in_cohort = n()) %>% 
  ungroup() %>% 
  select(mrn, weight_adjusted_dose, drug, 
         drug_tertile_sample1_to_last#,
         # munber_of_time_received_by_patient, munber_of_time_received_in_cohort
         ) %>% 
  # nab-paclitaxel needs it's own tertile as the administration is different
  # But update name now to add Bolton class
  mutate(drug = case_when(
    drug == "nab-paclitaxel"             ~ "paclitaxel",
    TRUE                                 ~ drug
  )) %>%
  # Add Bolton categories
  left_join(., drug_class, 
            by = c("drug" = "drug_name")) %>% 
  # calculate class patient's score - sum scores for each drug in a specific drug class for a patient
  group_by(mrn, narrow_drug_class_cytotoxic_only) %>% 
  mutate(class_score_sample1_to_last = case_when(
    !is.na(drug_tertile_sample1_to_last)                  ~ sum(drug_tertile_sample1_to_last, na.rm = TRUE)
  )) %>% 
  # select(mrn, drug, drug_tertile_sample1_to_last, class_score_sample1_to_last,
  #        narrow_drug_class_cytotoxic_only,
  #        munber_of_time_received_by_patient,
  #        everything()) %>% 
  # Keep 1 score per class per patient
  distinct(mrn, narrow_drug_class_cytotoxic_only, .keep_all = TRUE) %>% 
  select(-c(drug, weight_adjusted_dose, drug_tertile_sample1_to_last)) 


radiation_dose1 <- 
  treatment_dose %>% 
  select(-c(drug1, dose1, unit1, raw_dose1, tertile1,
            drug2, dose2, unit2, raw_dose2, tertile2,
            drug3, dose3, unit3, raw_dose3, tertile3)) %>%
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
  
  group_by(mrn) %>% 
  mutate(cumulative_dose_sample1_to_last = sum(eqd2, na.rm = TRUE)) %>% 
  ungroup() %>% 
  arrange(mrn, cumulative_dose_sample1_to_last) %>% 
  distinct(mrn, .keep_all = TRUE) %>% 
  select(mrn, cumulative_dose_sample1_to_last) %>% 
  # create tertile categories
  mutate(tertile_score_pre_firstseqsample = cut(cumulative_dose_sample1_to_last,
                                                breaks = quantile(cumulative_dose_sample1_to_last, 
                                                                  probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                                                include.lowest = TRUE,
                                                labels = c("1", "2", "3")))
  

write_csv(treatment_dose1,
          paste0("data/processed_data",
                 "/CHevolution_DrugClassDosingScore_", 
                 str_remove_all(today(), "-"), ".csv"))

write_csv(radiation_dose1,
          paste0("data/processed_data",
                 "/CHevolution_RadiationDosingScore_", 
                 str_remove_all(today(), "-"), ".csv"))

write_csv(treatment_dose1,
          paste0(path, "/ProcessedData",
                 "/CHevolution_DrugClassDosingScore_", 
                 str_remove_all(today(), "-"), ".csv"))

write_csv(radiation_dose1,
          paste0(path, "/ProcessedData",
                 "/CHevolution_RadiationDosingScore_", 
                 str_remove_all(today(), "-"), ".csv"))



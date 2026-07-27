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
    
  # Create sample range filter
  select(mrn, contains("collection_dt"), drug, dose, weight_adjusted_dose, drug_tertile_sample1_to_last, tx_start, tx_end) |> 
  # Then use tx_start or it will exclude 1 treatment that start before sample 2 and end after sample 2 
  mutate(range_before_sample2 = case_when(
    tx_start <= sample2_collection_dt                    ~ "Yes"
  )) |> 
  mutate(range_sample2_sample3 = case_when(
    is.na(range_before_sample2) &
      tx_start <= sample3_collection_dt                  ~ "Yes"
  )) |> 
  mutate(range_sample3_sample4 = case_when(
    is.na(range_before_sample2) &
      is.na(range_sample2_sample3) &
      tx_start <= sample4_collection_dt                  ~ "Yes"
  )) |> 
  mutate(range_sample4_sample5 = case_when(
    is.na(range_before_sample2) &
      is.na(range_sample2_sample3) &
      is.na(range_sample3_sample4) &
      tx_start <= sample5_collection_dt                  ~ "Yes"
  )) |> 
  mutate(range_sample5_sample6 = case_when(
    is.na(range_before_sample2) &
      is.na(range_sample2_sample3) &
      is.na(range_sample3_sample4) &
      is.na(range_sample4_sample5) &
      tx_start <= sample6_collection_dt                  ~ "Yes"
  )) |> 
  mutate(range_sample6_sample7 = case_when(
    is.na(range_before_sample2) &
      is.na(range_sample2_sample3) &
      is.na(range_sample3_sample4) &
      is.na(range_sample4_sample5) &
      is.na(range_sample5_sample6) &
      tx_start <= sample7_collection_dt                  ~ "Yes"
  )) |> 

    
    
    
    
    
  select(mrn, weight_adjusted_dose, drug, 
         drug_tertile_sample1_to_last, contains("range_")#,
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
  group_by(mrn, narrow_drug_class_cytotoxic_only, range_before_sample2) %>% 
  mutate(class_score_before_sample2 = case_when(
    !is.na(drug_tertile_sample1_to_last) &
      range_before_sample2 == "Yes"                       ~ sum(drug_tertile_sample1_to_last, na.rm = TRUE)
  )) |> 
  group_by(mrn, narrow_drug_class_cytotoxic_only, range_sample2_sample3) %>% 
  mutate(class_score_sample2_to_sample3 = case_when(
    !is.na(drug_tertile_sample1_to_last) &
      range_sample2_sample3 == "Yes"                      ~ sum(drug_tertile_sample1_to_last, na.rm = TRUE)
  )) |> 
  group_by(mrn, narrow_drug_class_cytotoxic_only, range_sample3_sample4) %>% 
  mutate(class_score_sample3_to_sample4 = case_when(
    !is.na(drug_tertile_sample1_to_last) &
      range_sample3_sample4 == "Yes"                      ~ sum(drug_tertile_sample1_to_last, na.rm = TRUE)
  )) |> 
  group_by(mrn, narrow_drug_class_cytotoxic_only, range_sample4_sample5) %>% 
  mutate(class_score_sample4_to_sample5 = case_when(
    !is.na(drug_tertile_sample1_to_last) &
      range_sample4_sample5 == "Yes"                      ~ sum(drug_tertile_sample1_to_last, na.rm = TRUE)
  )) |> 
  group_by(mrn, narrow_drug_class_cytotoxic_only, range_sample5_sample6) %>% 
  mutate(class_score_sample5_to_sample6 = case_when(
    !is.na(drug_tertile_sample1_to_last) &
      range_sample5_sample6 == "Yes"                      ~ sum(drug_tertile_sample1_to_last, na.rm = TRUE)
  )) |> 
  group_by(mrn, narrow_drug_class_cytotoxic_only, range_sample6_sample7) %>% 
  mutate(class_score_sample6_to_sample7 = case_when(
    !is.na(drug_tertile_sample1_to_last) &
      range_sample6_sample7 == "Yes"                      ~ sum(drug_tertile_sample1_to_last, na.rm = TRUE)
  )) |> 
  group_by(mrn, narrow_drug_class_cytotoxic_only) |> 
  fill(contains("class_score"), .direction = "updown") |> 
  ungroup() |> 

  # Keep 1 score per class per patient
  distinct(mrn, narrow_drug_class_cytotoxic_only, .keep_all = TRUE) %>% 
  select(-c(drug, weight_adjusted_dose, drug_tertile_sample1_to_last, general_drug_class,
            starts_with("range_"))) 


radiation_dose1 <- 
  treatment_dose %>% 
  filter(!is.na(rad_dose)) %>% 
  select(-c(drug1, dose1, unit1, raw_dose1, tertile1,
            drug2, dose2, unit2, raw_dose2, tertile2,
            drug3, dose3, unit3, raw_dose3, tertile3,
            x0 : tmn_dx, wt_kg : d, class)) %>%
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
  select(-rad_before_sample1, -rad_after_last) |> 
  # Approach 2- tertile for each radiation
  mutate(rad_eqd2_tertile = cut(eqd2,
                                breaks = quantile(eqd2, 
                                                  probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                                include.lowest = TRUE,
                                labels = c("1", "2", "3"))) |> 
  mutate(rad_eqd2_tertile = as.numeric(rad_eqd2_tertile)) |> 
  # Approach 2- tertile for each radiation
  # Create sample range filter
  mutate(range_before_sample2 = case_when(
    tx_start <= sample2_collection_dt                    ~ "Yes"
  )) |> 
  mutate(range_sample2_sample3 = case_when(
    is.na(range_before_sample2) &
      tx_start <= sample3_collection_dt                  ~ "Yes"
  )) |> 
  mutate(range_sample3_sample4 = case_when(
    is.na(range_before_sample2) &
      is.na(range_sample2_sample3) &
      tx_start <= sample4_collection_dt                  ~ "Yes"
  )) |> 
  mutate(range_sample4_sample5 = case_when(
    is.na(range_before_sample2) &
      is.na(range_sample2_sample3) &
      is.na(range_sample3_sample4) &
      tx_start <= sample5_collection_dt                  ~ "Yes"
  )) |> 
  mutate(range_sample5_sample6 = case_when(
    is.na(range_before_sample2) &
      is.na(range_sample2_sample3) &
      is.na(range_sample3_sample4) &
      is.na(range_sample4_sample5) &
      tx_start <= sample6_collection_dt                  ~ "Yes"
  )) |> 
  mutate(range_sample6_sample7 = case_when(
    is.na(range_before_sample2) &
      is.na(range_sample2_sample3) &
      is.na(range_sample3_sample4) &
      is.na(range_sample4_sample5) &
      is.na(range_sample5_sample6) &
      tx_start <= sample7_collection_dt                  ~ "Yes"
  )) |> 
  # calculate patient's score - sum tertile for each rad
  group_by(mrn) %>% 
  mutate(rad_score_sample1_to_last = case_when(
    !is.na(rad_eqd2_tertile)                              ~ sum(rad_eqd2_tertile, na.rm = TRUE)
  )) %>% 
  group_by(mrn, range_before_sample2) %>% 
  mutate(rad_score_before_sample2 = case_when(
    !is.na(rad_eqd2_tertile) &
      range_before_sample2 == "Yes"                       ~ sum(rad_eqd2_tertile, na.rm = TRUE)
  )) |> 
  group_by(mrn, range_sample2_sample3) %>% 
  mutate(rad_score_sample2_to_sample3 = case_when(
    !is.na(rad_eqd2_tertile) &
      range_sample2_sample3 == "Yes"                      ~ sum(rad_eqd2_tertile, na.rm = TRUE)
  )) |> 
  group_by(mrn, range_sample3_sample4) %>% 
  mutate(rad_score_sample3_to_sample4 = case_when(
    !is.na(rad_eqd2_tertile) &
      range_sample3_sample4 == "Yes"                      ~ sum(rad_eqd2_tertile, na.rm = TRUE)
  )) |> 
  group_by(mrn, range_sample4_sample5) %>% 
  mutate(rad_score_sample4_to_sample5 = case_when(
    !is.na(rad_eqd2_tertile) &
      range_sample4_sample5 == "Yes"                      ~ sum(rad_eqd2_tertile, na.rm = TRUE)
  )) |> 
  group_by(mrn, range_sample5_sample6) %>% 
  mutate(rad_score_sample5_to_sample6 = case_when(
    !is.na(rad_eqd2_tertile) &
      range_sample5_sample6 == "Yes"                      ~ sum(rad_eqd2_tertile, na.rm = TRUE)
  )) |> 
  group_by(mrn, range_sample6_sample7) %>% 
  mutate(rad_score_sample6_to_sample7 = case_when(
    !is.na(rad_eqd2_tertile) &
      range_sample6_sample7 == "Yes"                      ~ sum(rad_eqd2_tertile, na.rm = TRUE)
  )) |> 
  # Approach 1 -1 tertile for cumulative dose over all ranges
  group_by(mrn) %>% 
  mutate(cumulative_dose_sample1_to_last = sum(eqd2, na.rm = TRUE)) %>% 
  ungroup() %>% 
  # Approach 1 -1 tertile over cumulative dose 
  arrange(mrn, cumulative_dose_sample1_to_last) |> 
  select(-c(tx_start, tx_end, rad_site : sample7_collection_dt)) |> 
  group_by(mrn, across(contains("range_")), across(contains("rad_score")), cumulative_dose_sample1_to_last) |> 
  summarise_at(vars(eqd2, rad_eqd2_tertile), str_c, collapse = " + ") %>% 
  ungroup() |> 
  
  # distinct(mrn, .keep_all = TRUE) %>% # Used summarize instead
  # select(mrn, cumulative_dose_sample1_to_last) %>% 
  # create tertile categories
  mutate(tertile_cumdose = cut(cumulative_dose_sample1_to_last,
                               breaks = quantile(cumulative_dose_sample1_to_last, 
                                                 probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                               include.lowest = TRUE,
                               labels = c("1", "2", "3"))) |> 
  select(mrn, eqd2, cumulative_dose_sample1_to_last, tertile_cumdose,
         everything(), rad_eqd2_tertile) |> 
  mutate(tertile_summedtertile_of_alldoses = cut(rad_score_before_sample2,
                                                 breaks = quantile(rad_score_before_sample2, 
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



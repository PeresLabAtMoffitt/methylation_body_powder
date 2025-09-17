## Import library
library(tidyverse)
# library(haven)


###################################################################### I ### Load data
# path1 <- fs::path("", "Volumes", "Peres_Research", "AACES", "Concept_risk by tumor immune profiles")

path_clinical <- fs::path("", "Volumes", "Peres_Research", "K99_R00", "Image analysis data")
clinical_ncocs <-
  read_csv(paste0(path_clinical, 
                  "/AACES MCC18207 mIF Data/ncocs_all.csv"))
clinical_aaces <-
  read_csv(paste0(path_clinical, 
                  "/AACES MCC18207 mIF Data/aaces_all.csv"))


load("/Volumes/Peres_Research/Methylomics R01/Data from Lucas/cleaned_07082022.rda")
# load("/Volumes/Peres_Research/Methylomics R01/Data from Lucas/cleaned_05042021.rda")
load(file="/Volumes/Peres_Research/Methylomics R01/Data from Lucas/epitoc_with_ids.RData")


###################################################################### II ### Clean clinical data
patient_id <- phenoclean %>% 
  select(Complete.Barcode, suid = Sample_Name, ocwaaid)

clinical_cleaning <- function(data){
  data <- data %>% 
    mutate(suid = as.character(suid)) %>% 
    mutate(casecon = case_when(
      casecon == 1                                       ~ "Case",
      casecon == 2                                       ~ "Control"
    )) %>% 
    mutate(site = case_when(
      site == "BA" |
        site == "LA" |
        site == "AL"                                     ~ "South Central",
      site == "GA" |
        site == "NC" |
        site == "SC" |
        site == "TN"                                     ~ "South/Mid-Atlantic",
      site == "OH" |
        site == "NJ" |
        site == "IL" |
        site == "MI"                                     ~ "North",
      TRUE                                               ~ site
    ), site = factor(site, levels = c("South/Mid-Atlantic", "South Central",
                                      "North"))
    ) %>% 
    mutate(race = case_when(
      race == 1                                          ~ "White",
      race == 2                                          ~ "Black",
      race == 3                                          ~ "Biracial",
      TRUE                                                ~ NA_character_
    )) %>% 
    mutate(hispanic = case_when(
      hispanic == 1                                      ~ "Hispanic",
      hispanic == 2                                      ~ "Non-hispanic",
      TRUE                                                ~ NA_character_
    )) %>% 
    mutate(refage_cat = case_when(
      refage < 50                      ~ "<50",
      refage >= 50 &
        refage < 60                    ~ "50-59",
      refage >= 60 &
        refage < 70                    ~ "60-69",
      refage >= 70 &
        refage < 80                    ~ "70-79",
      refage >= 80                     ~ "≥80"
    )) %>% 
    mutate(BMI_classification = case_when(
      BMI_recent < 18.5	~ "underweight",
      BMI_recent >=18.5 & BMI_recent <25 ~ "normal",
      BMI_recent >=25.0 & BMI_recent <30 ~ "overweight",
      BMI_recent >=30.0 & BMI_recent <35 ~ "obesity I",
      BMI_recent >=35.0 & BMI_recent <40 ~ "obesity II",
      BMI_recent >= 40.0 ~ "obesity III"
    )) %>% 
    mutate(BMI_recent_grp = case_when(
      BMI_recent < 25                  ~ "<25",
      BMI_recent >= 25 &
        BMI_recent < 30                ~ "25-29",
      BMI_recent >= 30 &
        BMI_recent < 35                ~ "30-35",
      BMI_recent >= 35                 ~ "≥35"
    )) %>% 
    mutate(BMI_YA_grp = case_when(
      BMI_YA < 20                      ~ "<20",
      BMI_YA >= 20 &
        BMI_YA < 25                    ~ "20-24",
      BMI_YA >= 25                     ~ "≥25"
    )) %>% 
    mutate(BMI_recent_grp = case_when(
      BMI_recent <25                       ~ "<25 kg/m2",
      BMI_recent >=25 & BMI_recent <= 30   ~ "25-30 kg/m2",
      BMI_recent > 30.0                    ~ ">30 kg/m2"
    )) %>% 
    mutate(BMI_YA_grp = case_when(
      BMI_YA <25                           ~ "<25 kg/m2",
      BMI_YA >=25 & BMI_YA <= 30           ~ "25-30 kg/m2",
      BMI_YA > 30.0                        ~ ">30 kg/m2"
    )) %>% 
    mutate(birthplace = case_when(
      birthplace == 1                                    ~ "born in United States",
      birthplace == 2                                    ~ "born outside of United States",
      TRUE                                                ~ NA_character_
    )) %>% 
    mutate(education = case_when(
      education == 1                                      ~ "high school graduate/GED or less",
      education == 2                                      ~ "some college",
      education == 3                                      ~ "college graduate",
      education == 4                                      ~ "graduate/professional school",
      TRUE                                                ~ NA_character_
    )) %>% 
    mutate(histology = case_when(
      histology == 1                                     ~ "Serous",
      histology == 2                                     ~ "Endometrioid",
      histology == 3                                     ~ "Clear cell",
      histology == 4                                     ~ "Mucinous",
      histology == 5                                     ~ "Carcinosarcoma",
      histology == 6                                     ~ "Carcinoma, NOS",
      histology == 7                                     ~ "Other specified epithelial ovarian cancer \n(e.g. Malignant Brenner, mixed)",
      histology == 8                                     ~ "Epithelial, NOS",
      histology == 9                                     ~ "Synchronous",
      TRUE                                                ~ NA_character_
    )) %>% 
    mutate(behavior = case_when(
      behavior == 1                                      ~ "Borderline",
      behavior == 2                                      ~ "Invasive",
      TRUE                                                ~ NA_character_
    )) %>% 
    mutate(stage_cat = case_when(
      stage == 1 |
        stage == 2                                         ~ "Early",
      stage == 3                                         ~ "Late",
      TRUE                                                ~ NA_character_
    )) %>% 
    mutate(stage = case_when(
      stage == 1                                         ~ "Localized",
      stage == 2                                         ~ "Regional",
      stage == 3                                         ~ "Distant",
      TRUE                                                ~ NA_character_
    )) %>% 
    mutate(stage = case_when(
      stage == "Localized"                 ~ "Regional",
      TRUE                                 ~ stage
    )) %>% 
    mutate(grade = case_when(
      grade == 1                                         ~ "well differentiated",
      grade == 2                                         ~ "moderately differentiated",
      grade == 3                                         ~ "poorly differentiated",
      grade == 4                                         ~ "undifferentiated",
      TRUE                                               ~ NA_character_
    )) %>% 
    mutate(histotype1 = case_when(
      histotype == 1                                     ~ "high-grade serous",
      histotype == 2                                     ~ "low-grade serous",
      histotype == 5 |
        histotype == 10                                  ~ "mucinous",
      histotype == 3                                     ~ "endometrioid",
      histotype == 4                                     ~ "clear cell",
      histotype %in% (6:13)                              ~ "other epithelial", # Will not take the 10
      TRUE                                               ~ NA_character_
    )) %>% 
    mutate(histotype2 = case_when(
      histotype == 1                                     ~ "high-grade serous",
      histotype %in% (2:13)                              ~ "non-high-grade serous",
      TRUE                                               ~ NA_character_
    )) %>%
    mutate(histotype = case_when(
      histotype == 1                                     ~ "high-grade serous",
      histotype == 2                                     ~ "low-grade serous",
      histotype == 3                                     ~ "endometrioid",
      histotype == 4                                     ~ "clear cell",
      histotype == 5                                     ~ "mucinous",
      histotype == 6                                     ~ "carcinosarcoma",
      histotype == 7                                     ~ "other epithelial ovarian cancer \n(e.g. Malignant Brenner, mixed, carcinoma, NOS)",
      histotype == 9                                     ~ "serous borderline",
      histotype == 10                                    ~ "mucinous borderline",
      histotype == 11                                    ~ "other epithelial borderline",
      histotype == 13                                    ~ "synchronous",
      TRUE                                               ~ NA_character_
    )) %>% 
    mutate_at(c("refage", "pregmos", "agefbirth", "agelbirth", 
                "height", "wt_recent", "wt_YA", "wtgain", "BMI_recent", "BMI_YA",
                "talcgenfreq", "talcgendur", 
                "talcnongenfreq", "talcnongendur",
                "talcagegen", "talcagenongen"), 
              ~ case_when(
                . %in% c("88","98", "99", 
                         "888", "998", "999")            ~ NA_real_,
                TRUE                                     ~ as.numeric(.))
    ) %>% 
    mutate_at(c("talcever", "talcpartner", "talc_occ",
                "neoadj_treat", "adj_treat", "dblkstat_treat"),
              ~ case_when(
                . == 1                                              ~ "Yes",
                . == 2                                              ~ "No",
                TRUE                                                ~ NA_character_
              )) %>% 
    mutate(talcever = factor(talcever, levels = c("No","Yes"))) %>% 
    mutate(talcgen = case_when(
      talcever == "No"                                          ~ "No",
      talcgen == 1                                              ~ "Yes",
      talcgen == 2                                              ~ "No",
      TRUE                                                      ~ NA_character_
    ), talcgen = factor(talcgen, levels = c("No","Yes"))
    ) %>% 
    mutate(talcnongen = case_when(
      talcever == "No"                                          ~ "No",
      talcnongen == 1                                           ~ "Yes",
      talcnongen == 2                                           ~ "No",
      TRUE                                                      ~ NA_character_
    ), talcnongen = factor(talcnongen, levels = c("No","Yes"))
    )
}

clinical_aaces <- clinical_cleaning(clinical_aaces)
clinical_ncocs <- clinical_cleaning(clinical_ncocs)


clinical_data <- bind_rows(clinical_aaces, 
                           clinical_ncocs) %>% 
  inner_join(patient_id, ., 
             by = c("suid", "ocwaaid" = "OCWAAID")) 


###################################################################### III ### Clean methylation data----
epitoc <- epitoc %>% as_tibble(rownames = "Complete.Barcode")

clock_data <- inner_join(epitoc, clinical_data,
                         by = "Complete.Barcode") %>%
  mutate(tert_irS = ntile(irS, 3)) %>%
  mutate(tert_irS = case_when(
    tert_irS == 1                                ~ 'Low',
    tert_irS == 2                                ~ 'Medium',
    tert_irS == 3                                ~ 'High'
  ), tert_irS = factor(tert_irS, levels = c("Low", "Medium", "High"))) %>% 
  select(suid, ocwaaid, Complete.Barcode, irS, tert_irS, contains("talc"), everything())


write_rds(clock_data, paste0(here::here(), 
                             "/data/processed data",
                             "/clock_data.rds"))


# END


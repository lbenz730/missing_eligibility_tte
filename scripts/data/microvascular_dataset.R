library(tidyverse)
library(glue)
library(haven)
library(lubridate)
library(arrow)

source('scripts/util/helpers.R')
options(dplyr.summarise.inform = F)

library(furrr)
n_cores <- 12
plan(future::multisession(workers = n_cores))
options(future.globals.maxSize = 8 * 1024^3)


### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

####################################
### Step 1: Create Cases Dataset ###
####################################

### Read in Data Files
cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/durable_cases.sas7bdat'))
case_heights <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/durable_modal_ht.sas7bdat'))
case_races <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/raw_race_cases.sas7bdat'))
weights <- read_parquet(glue('{data_dir}/all_weights.parquet')) ### Cleaned Weights w/ outliers removed

case_races <- 
  case_races %>% 
  select('subject_id' = durable_studyid, 'race' = race1)

df_cases <- 
  cases %>% 
  select('subject_id' = durable_studyid,
         'birth_date' = BIRTH_DATE,
         'site' = site,
         index_date,
         'gender' = GENDER) %>% 
  left_join(case_races, by = 'subject_id')

### Files for Diabetes Diagnosis
hba1c_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/hba1c_glucose_labs_cases.sas7bdat'))
diabetes_icd9_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/raw_diabetes_dx_cases.sas7bdat'))
diabetes_rx_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/all_diabetes_rx.sas7bdat'))

diabetes_dx <- function(icd9_code) {
  dx <- 
    case_when(substr(icd9_code, 1, 3) %in% c('250', 'E08', 'E09', 'E10', 'E11', 'E13') ~ T,
              substr(icd9_code, 1, 5) %in% c('362.0', '357.2') ~ T,
              substr(icd9_code, 1, 6) == '366.41' ~ T,
              T ~ F)
  
  return(dx)
}

### Diabetes via Labs
### A1c Measurement of at least 6.5 or Fasting Glucose of at least 126 
### at the most recent measurement w/in 2 years of surgery
diabetes_cases_labs <- 
  hba1c_cases %>% 
  rename('subject_id' = durable_studyid) %>% 
  filter(lab_date <= index_date) %>% ### Need this for cases to remove measures after surgery
  filter(index_date <= lab_date %m+% years(2)) %>% 
  group_by(subject_id, test_type) %>%  ### Most recent measurement (of this test type)
  filter(lab_date == max(lab_date)) %>% 
  ungroup() %>% 
  filter( (test_type == 'HGBA1C' & result >= 6.5) | (test_type == 'GLU_F' & result >= 126 & toupper(result_unit) == 'MG/DL') ) %>% 
  mutate('diabetic' = 1) %>% 
  distinct(subject_id, index_date)


### Diabetes via Rx
### Current Rx for diabetes medication
### I am assuming that current means index_date is w/in rxsup (supply days) or rxdate
diabetes_cases_rx <- 
  diabetes_rx_cases %>% 
  select('subject_id' = durable_studyid, rxdate, rxsup, index_date) %>% 
  filter(rxdate <= index_date) %>% 
  filter(index_date <= rxdate + rxsup) %>% 
  distinct(subject_id, index_date)

diabetes_cases <- 
  diabetes_cases_labs %>% 
  bind_rows(diabetes_cases_rx) %>% 
  distinct() %>% 
  filter(index_date >= '2005-01-01', 
         index_date <= '2011-12-31') %>% 
  select(-index_date)

### Get Baseline BMI, A1c, and serum creatine levels
baseline_bmi <- 
  weights %>% 
  filter(subject_id %in% diabetes_cases$subject_id) %>%
  filter(measure_date <= index_date) %>% 
  filter(measure_date %m+% years(2)  >= index_date) %>% 
  group_by(subject_id) %>% 
  filter(measure_date == max(measure_date)) %>% 
  filter(bmi == max(bmi)) %>% 
  ungroup() %>% 
  select(subject_id, bmi)


baseline_a1c <- 
  hba1c_cases %>% 
  rename('subject_id' = durable_studyid) %>% 
  filter(subject_id %in% diabetes_cases$subject_id) %>%
  filter(lab_date <= index_date) %>% 
  filter(lab_date %m+% years(2)  >= index_date) %>% 
  group_by(subject_id, test_type) %>%  ### Most recent measurement (of this test type)
  filter(lab_date == max(lab_date)) %>% 
  ungroup() %>% 
  filter( (test_type == 'HGBA1C' & result >= 6.5) | (test_type == 'GLU_F' & result >= 126 & toupper(result_unit) == 'MG/DL') ) %>% 
  mutate('hgba1c' = case_when(test_type == 'HGBA1C' ~ result,
                              test_type == 'GLU_F' ~ (result + 46.7)/28.7)) %>% 
  group_by(subject_id) %>%  ### Most recent measurement 
  filter(lab_date == max(lab_date)) %>% 
  filter(hgba1c == max(hgba1c)) %>% 
  ungroup() %>% 
  select(subject_id, hgba1c)

kidney_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/kidney_lab_cases_final.sas7bdat'))

ckd_epi <- function(scr, age, sex) {
  alpha <- ifelse(sex == 'M', -0.302, -0.241)
  kappa <- ifelse(sex == 'M', 0.9, 0.7)
  
  eGFR <- 142 * pmin(scr/kappa, 1)^alpha * pmax(scr/kappa, 1)^(-1.2) *  0.9938^age * ifelse(sex == 'M', 1, 1.012)
  return(eGFR)
}

eGFR_cases <- 
  kidney_cases %>% 
  filter(test_Type == 'CREATININE') %>% 
  select('subject_id' = durable_studyid,
         lab_date,
         result_num) %>% 
  left_join(df_cases, by = 'subject_id') %>% 
  mutate('eGFR' = ckd_epi(scr = result_num,
                          age = as.numeric((lab_date - birth_date))/365.25,
                          sex = gender))

baseline_eGFR <- 
  eGFR_cases %>% 
  filter(subject_id %in% diabetes_cases$subject_id) %>%
  filter(lab_date <= index_date) %>% 
  filter(lab_date %m+% years(2)  >= index_date) %>% 
  group_by(subject_id) %>%  ### Most recent measurement 
  filter(lab_date == max(lab_date)) %>% 
  filter(eGFR == max(eGFR)) %>% 
  ungroup() %>% 
  select(subject_id, 'scr' = result_num, eGFR)

### Exclusions/Censoring/Outcomes
###
#### Exclusion:
### (1) Not Enrolled continuously for the year before surgery
### (2) Pregnancy in year before surgery
### (3) Pre-existing outcomes
###
### Censoring
### (1) Incident Cancer
### (2) Disenrollment
### (3) Death
### (4) End of Study
### (5) No BMI/BP measures for 1 month
###
### Outcomes:
###
### (1) Neuropathy 
###       dx %in% c('357.2', '250.60', '250.62')
### (2) Retinopathy
###     (a) ICD-9 Code 362.0x
###     (b) Retinopathy PX code in the below along w/ ICD-9 code of 250.5x
### (3) Nephropathy
###       2 measures of egrf of < 60 seperated by at least 90 days w/ no measures between of 60+

case_enrollment <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/all_enrollment_collapsed.sas7bdat'))
include_enrollment <- 
  case_enrollment %>% 
  filter(enr_end >= index_date, enr_start <= index_date %m-% years(1)) %>% 
  distinct('subject_id' = durable_studyid,'censor_enrollment' = enr_end)

### I'm assuming any code here means pregnant
case_pregnancy <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/pregnant_dx.sas7bdat'))
pregnancy_exclude <- 
  case_pregnancy %>% 
  filter(adate <= index_date, adate >= index_date %m-% years(1)) %>% 
  select('subject_id' = durable_studyid)

pregnancy_censor <- 
  case_pregnancy %>% 
  filter(adate > index_date) %>% 
  group_by('subject_id' = durable_studyid) %>% 
  summarise('censor_pregnancy' = min(adate)) 

### Death
death_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/death_cases.sas7bdat'))
death_cases <- 
  death_cases %>% 
  select('subject_id' = durable_studyid,
         'censor_death' = deathdt)

### Incident Cancer
cancer_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/raw_cancer_dx_cases.sas7bdat'))
cancer_cases <- 
  cancer_cases %>% 
  filter(index_date <= adate) %>% 
  group_by('subject_id' = durable_studyid) %>% 
  summarise('censor_cancer' = min(adate))

### No BP/WT measures for 13 months
blood_pressure_cases <- read_parquet(glue('{ehr_dir}/parquet_files/cases/blood_pressure_cases.parquet'))
df_measures <- 
  weights %>% 
  filter(subject_id %in% cases$durable_studyid) %>% 
  select(subject_id, measure_date, index_date) %>% 
  bind_rows(blood_pressure_cases %>% select('subject_id' = durable_studyid, measure_date, index_date)) %>% 
  filter(measure_date >= index_date) %>% 
  arrange(subject_id, measure_date) %>% 
  mutate('censor_measurement' = measure_date + 30.44 * 13) %>% 
  group_by(subject_id) %>% 
  filter(lead(measure_date, default = last(measure_date) +  30.44 * 13) >= censor_measurement) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(subject_id, censor_measurement)


### Outcomes (and exclude patients w/ outcomes before surgery)
### Neuropathy
neuropathy_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/raw_neuropathy_dx_cases.sas7bdat'))
neuropathy_dx <- 
  neuropathy_cases %>% 
  filter(dx %in% c('357.2', '250.60', '250.62')) %>% 
  group_by('subject_id' = durable_studyid) %>% 
  summarise('neuropathy' = min(adate)) %>% 
  ungroup()

### Retinopathy
### (1) ICD-9 Code 362.0x
retinopathy_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/retinopathy_dx.sas7bdat'), encoding = 'latin1')
retinopathy_dx <- 
  retinopathy_cases %>% 
  filter(grepl('^362\\.0', dx)) %>% 
  group_by('subject_id' = durable_studyid, index_date) %>% 
  summarise('retinopathy' = min(adate)) %>% 
  ungroup()

### (2) Retinopathy PX code in the below along w/ ICD-9 code of 250.5x
retinopathy_px <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/retinopathy_px.sas7bdat'), encoding = 'latin1')
retinopathy_dx2 <- 
  inner_join(
    # Retinopathy PX code in the below
    retinopathy_px %>% 
      filter(px %in% c('67005', '67010', '67015', '67025',
                       '67027', '67028', '67036', '67038',
                       '67039', '67040', '67041', '67042',
                       '67043', as.character(67101:67112),
                       '67113', as.character(67141:67145),
                       as.character(67208:67218),
                       as.character(67220:67225),
                       as.character(67227:67228)
      )) %>% 
      group_by('subject_id' = durable_studyid) %>% 
      summarise('cpt_date' = min(adate)),
    
    #ICD-9 code of 250.5x
    diabetes_icd9_cases %>% 
      filter(grepl('^250\\.5', dx)) %>% 
      group_by('subject_id' = durable_studyid) %>% 
      summarise('icd9_date' = min(adate)),
    
    
    by = 'subject_id'
  ) %>% 
  mutate('retinopathy' = pmax(cpt_date, icd9_date))

retinopathy_dx <- 
  retinopathy_dx %>% 
  bind_rows(retinopathy_dx2) %>% 
  group_by(subject_id) %>% 
  summarise('retinopathy' = min(retinopathy))


### Nephropathy
## 2 measures of egrf of < 60 
## seperated by at least 90 days
## no measures between of 60+
diagnose_nephropathy <- function(egrf, date) {
  n <- length(date)
  if(n > 2) {
    for(t1 in 1:(n-1)) {
      for(t2 in 2:n) {
        if(as.numeric(date[t2] - date[t1]) >= 90 & egrf[t1] < 60 & egrf[t2] < 60) {
          if(all(egrf[t1:t2] < 60)) {
            return(date[t2])
          }
        }
      }
    }
  }
  return(NA_Date_)
}

nephropathy_dx <- 
  eGFR_cases %>% 
  group_by(subject_id) %>% 
  group_split() %>% 
  future_map_dfr(~{
    tibble('subject_id' = .x$subject_id[1],
           'nephropathy' = diagnose_nephropathy(.x$eGFR, .x$lab_date))
  },
  .progress = T)


### Final Dataset
case_population <- 
  df_cases %>% 
  ### Gets All subjects w/ diabetes who underwent surgery between 2005-2011
  inner_join(diabetes_cases, by = 'subject_id') %>% 
  
  ### Exclusion Criteria
  inner_join(include_enrollment, by = 'subject_id') %>% 
  anti_join(pregnancy_exclude, by = 'subject_id') %>% 
  
  ### Outcomes
  left_join(neuropathy_dx, by = 'subject_id') %>% 
  left_join(retinopathy_dx, by = 'subject_id') %>% 
  left_join(nephropathy_dx, by = 'subject_id') %>% 
  
  ### Remove Pre-Existing Conditions
  filter(is.na(neuropathy) | neuropathy > index_date) %>% 
  filter(is.na(retinopathy) | retinopathy > index_date) %>% 
  filter(is.na(nephropathy) | nephropathy > index_date) %>% 
  
  ### Baseline Measurements
  inner_join(baseline_bmi,  by = 'subject_id') %>%
  inner_join(baseline_a1c,  by = 'subject_id') %>%
  inner_join(baseline_eGFR,  by = 'subject_id') %>%
  
  ### Censoring (enrollment above)
  left_join(death_cases, by = 'subject_id') %>% 
  left_join(cancer_cases, by = 'subject_id') %>% 
  left_join(df_measures, by = 'subject_id') %>% 
  mutate('end_of_study' = as.Date('2015-09-30')) 

write_parquet(case_population, glue('{data_dir}/microvascular_cases.parquet'))
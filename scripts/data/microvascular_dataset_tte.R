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

if(!dir.exists(glue('{data_dir}/microvascular_tte'))) {
  dir.create(glue('{data_dir}/microvascular_tte'))
}

###########################
###### Demographics #######
###########################
cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/durable_cases.sas7bdat'))
case_races <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/raw_race_cases.sas7bdat'))
controls <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/raw_race_controls.sas7bdat')) 

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

df_controls <- 
  controls %>% 
  select('subject_id' = control_studyid,
         'birth_date' = birth_date,
         'site' = site_id,
         'gender' = Gender,
         'race' = race1)

df_subjects <- 
  bind_rows(df_cases, df_controls)

write_parquet(df_subjects, glue('{data_dir}/microvascular_tte/subjects.parquet'))


#####################################
##### Diabetes Related Files ########
#####################################
hba1c_cases <- read_parquet(glue('{ehr_dir}/parquet_files/cases/hba1c_glucose_labs_cases.parquet'))
hba1c_controls <- read_parquet(glue('{ehr_dir}/parquet_files/controls/hba1c_glucose_labs_controls.parquet'))
diabetes_icd9_cases <- read_parquet(glue('{ehr_dir}/parquet_files/cases/raw_diabetes_dx_cases.parquet'))
diabetes_icd9_controls <- read_parquet(glue('{ehr_dir}/parquet_files/controls/raw_diabetes_dx_controls.parquet'))
diabetes_rx_cases <- read_parquet(glue('{ehr_dir}/parquet_files/cases/all_diabetes_rx.parquet'))
diabetes_rx_controls <- read_parquet(glue('{ehr_dir}/parquet_files/controls/rx_diabetes_controls.parquet'))

diabetes_labs <- 
  hba1c_controls %>% 
  mutate('index_date' = NA_Date_) %>% 
  rename('subject_id' = control_studyid) %>% 
  bind_rows(hba1c_cases %>% rename('subject_id' = durable_studyid)) %>% 
  filter(test_type %in% c('HGBA1C', 'GLU_F')) %>% 
  select(subject_id, lab_date, test_type, result, result_unit)

write_parquet(diabetes_labs, glue('{data_dir}/microvascular_tte/diabetes_labs.parquet'))

diabetes_rx <- 
  diabetes_rx_controls %>% 
  mutate('index_date' = NA_Date_) %>% 
  rename('subject_id' = control_studyid) %>% 
  bind_rows(diabetes_rx_cases %>% rename('subject_id' = durable_studyid)) %>% 
  select(subject_id, ndc, insulin_flg, rxdate, rxsup) %>% 
  mutate('rx_end_date' = rxdate + rxsup) 

write_parquet(diabetes_rx, glue('{data_dir}/microvascular_tte/diabetes_rx.parquet'))

diabetes_dx <- 
  diabetes_icd9_controls %>% 
  mutate('index_date' = NA_Date_) %>% 
  rename('subject_id' = control_studyid) %>% 
  bind_rows(diabetes_icd9_cases %>% rename('subject_id' = durable_studyid)) %>% 
  select(subject_id, adate, dx)

write_parquet(diabetes_dx, glue('{data_dir}/microvascular_tte/diabetes_dx.parquet'))


######################### 
### Kidney Covariates ###
#########################
ckd_epi <- function(scr, age, sex) {
  scr <- pmax(0.1, scr)
  alpha <- ifelse(sex == 'M', -0.302, -0.241)
  kappa <- ifelse(sex == 'M', 0.9, 0.7)
  
  eGFR <- 142 * pmin(scr/kappa, 1)^alpha * pmax(scr/kappa, 1)^(-1.2) *  0.9938^age * ifelse(sex == 'M', 1, 1.012)
  return(eGFR)
}

kidney_cases <- read_parquet(glue('{ehr_dir}/parquet_files/cases/kidney_lab_cases_final.parquet'))
kidney_controls <- read_parquet(glue('{ehr_dir}/parquet_files/controls/kidney_lab_control_final.parquet'))

df_kidney <- 
  kidney_controls %>% 
  filter(test_Type == 'CREATININE') %>% 
  rename('subject_id' = control_studyid) %>% 
  mutate('index_date' = NA_Date_) %>% 
  select(subject_id, index_date, lab_date, result_num) %>% 
  bind_rows(
    kidney_cases %>% 
      filter(test_Type == 'CREATININE') %>% 
      rename('subject_id' = durable_studyid) %>% 
      select(subject_id, index_date, lab_date, result_num)
  ) %>% 
  left_join(df_subjects %>% select(subject_id, birth_date, gender), by = 'subject_id') %>% 
  mutate('eGFR' = ckd_epi(scr = result_num,
                          age = as.numeric((lab_date - birth_date))/365.25,
                          sex = gender)) %>% 
  select(subject_id, lab_date, 'scr' = result_num, eGFR)

write_parquet(df_kidney, glue('{data_dir}/microvascular_tte/kidney_labs.parquet'))

######################
### Smoking Status ###
######################
smoking_cases <- read_parquet(glue('{ehr_dir}/parquet_files/cases/smoking_vdw.parquet'))
smoking_controls <- read_parquet(glue('{ehr_dir}/parquet_files/controls/smoking_vdw_controls.parquet'))

df_smoking <- 
  smoking_controls %>% 
  rename('subject_id' = control_studyid) %>% 
  mutate('index_date' = NA_Date_) %>% 
  bind_rows(smoking_cases %>% rename('subject_id' = durable_studyid)) %>% 
  mutate('smoking_status' = case_when(tobacco_use == 'Y' ~ 'current',
                                      tobacco_use == 'Q' ~ 'former',
                                      tobacco_use %in% c('N', 'P') ~ 'never',
                                      T ~ 'no_self_report')) %>% 
  select(subject_id, contact_date, smoking_status)

write_parquet(df_smoking, glue('{data_dir}/microvascular_tte/smoking.parquet'))


### Create files we need for Exclusions/Censoring/Outcomes
###
#### Exclusion:
### (1) Not Enrolled continuously for the year before trial
### (2) Pregnancy in year before trial
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
###       2 measures of egfr of < 60 seperated by at least 90 days w/ no measures between of 60+

### Enrollment 
case_enrollment <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/all_enrollment_collapsed.sas7bdat'))
control_enrollment <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/enrollment_collapsed_controls.sas7bdat'))

df_enrollment <- 
  case_enrollment %>% 
  rename('subject_id' = durable_studyid) %>% 
  bind_rows(control_enrollment %>% rename('subject_id' = control_studyid)) %>% 
  mutate('enr_1yr' = enr_start %m+% years(1)) %>% 
  select(subject_id, enr_start, enr_1yr, enr_end) 

write_parquet(df_enrollment, glue('{data_dir}/microvascular_tte/enrollment.parquet'))

### Pregnancy
### I'm assuming any code here means pregnant
case_pregnancy <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/pregnant_dx.sas7bdat'))
control_pregnancy <- read_parquet(glue('{ehr_dir}/parquet_files/controls/raw_pregnancy_dx_controls.parquet'))

df_pregnancy <- 
  control_pregnancy %>% 
  rename('subject_id' = control_studyid) %>% 
  bind_rows(case_pregnancy %>% rename('subject_id' = durable_studyid)) %>% 
  select(subject_id, adate, dx)

write_parquet(df_pregnancy, glue('{data_dir}/microvascular_tte/pregnancy.parquet'))

### Death
case_death <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/death_cases.sas7bdat'))
control_death <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/death_controls.sas7bdat'))

df_death <- 
  control_death %>% 
  rename('subject_id' = control_studyid) %>% 
  bind_rows(case_death %>% rename('subject_id' = durable_studyid)) %>% 
  select(subject_id, 'censor_death' = deathdt)

write_parquet(df_death, glue('{data_dir}/microvascular_tte/death.parquet'))

### Incident Cancer
case_cancer <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/raw_cancer_dx_cases.sas7bdat'))
control_cancer <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/raw_cancer_dx_controls.sas7bdat'))

df_cancer <- 
  control_cancer %>% 
  rename('subject_id' = control_studyid) %>% 
  bind_rows(case_cancer %>% rename('subject_id' = durable_studyid)) %>% 
  group_by(subject_id) %>% 
  summarise('censor_cancer' = min(adate))

write_parquet(df_cancer, glue('{data_dir}/microvascular_tte/cancer.parquet'))


### No BP/WT measures for 13 months
blood_pressure_cases <- read_parquet(glue('{ehr_dir}/parquet_files/cases/blood_pressure_cases.parquet'))
blood_pressure_controls <- read_parquet(glue('{ehr_dir}/parquet_files/controls/blood_pressure_controls.parquet'))
df_bp <- 
  blood_pressure_controls %>% 
  rename('subject_id' = control_studyid) %>% 
  bind_rows(blood_pressure_cases %>% rename('subject_id' = durable_studyid)) %>% 
  select(subject_id, measure_date, systolic, diastolic)
write_parquet(df_bp, glue('{data_dir}/blood_pressure.parquet'))

weights <- read_parquet(glue('{data_dir}/all_weights.parquet'))

df_measures <- 
  weights %>% 
  select(subject_id, measure_date) %>% 
  bind_rows(df_bp %>% select(subject_id, measure_date)) %>% 
  arrange(subject_id, measure_date) %>% 
  mutate('censor_measurement' = measure_date + 30.44 * 13) %>% 
  group_by(subject_id) %>% 
  filter(lead(measure_date, default = last(measure_date) +  30.44 * 13) >= censor_measurement) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(subject_id, censor_measurement)
write_parquet(df_measures, glue('{data_dir}/microvascular_tte/censor_measurements.parquet'))


### Outcomes (and exclude patients w/ outcomes before surgery)
### Neuropathy date of first dx
neuropathy_cases <- read_parquet(glue('{ehr_dir}/parquet_files/cases/raw_neuropathy_dx_cases.parquet'))
neuropathy_controls <- read_parquet(glue('{ehr_dir}/parquet_files/controls/raw_neuropathy_dx_controls.parquet'))

neuropathy_dx <- 
  neuropathy_controls %>% 
  rename('subject_id' = control_studyid) %>% 
  bind_rows(neuropathy_cases %>% rename('subject_id' = durable_studyid)) %>% 
  filter(dx %in% c('357.2', '250.60', '250.62')) %>% 
  group_by(subject_id) %>% 
  summarise('neuropathy' = min(adate)) %>% 
  ungroup()

write_parquet(neuropathy_dx, glue('{data_dir}/microvascular_tte/neuropathy_dx.parquet'))

### Retinopathy
retinopathy_dx_cases <- read_parquet(glue('{ehr_dir}/parquet_files/cases/retinopathy_dx.parquet'))
retinopathy_px_cases <- read_parquet(glue('{ehr_dir}/parquet_files/cases/retinopathy_px.parquet'))
retinopathy_dx_controls <- read_parquet(glue('{ehr_dir}/parquet_files/controls/retinopathy_dx_controls.parquet'))
retinopathy_px_controls <- read_parquet(glue('{ehr_dir}/parquet_files/controls/retinopathy_px_controls.parquet'))

### (1) ICD-9 Code 362.0x
retinopathy_dx <- 
  retinopathy_dx_controls %>% 
  rename('subject_id' = control_studyid) %>% 
  bind_rows(retinopathy_dx_cases %>% rename('subject_id' = durable_studyid)) %>% 
  filter(grepl('^362\\.0', dx)) %>% 
  group_by(subject_id) %>% 
  summarise('retinopathy' = min(adate)) %>% 
  ungroup()

### (2) Retinopathy PX code in the below along w/ ICD-9 code of 250.5x
retinopathy_dx2 <- 
  inner_join(
    # Retinopathy PX code in the below
    retinopathy_px_controls %>% 
      rename('subject_id' = control_studyid) %>% 
      bind_rows(retinopathy_px_cases %>% rename('subject_id' = durable_studyid)) %>% 
      filter(px %in% c('67005', '67010', '67015', '67025',
                       '67027', '67028', '67036', '67038',
                       '67039', '67040', '67041', '67042',
                       '67043', as.character(67101:67112),
                       '67113', as.character(67141:67145),
                       as.character(67208:67218),
                       as.character(67220:67225),
                       as.character(67227:67228)
      )) %>% 
      group_by(subject_id) %>% 
      summarise('cpt_date' = min(adate)),
    
    #ICD-9 code of 250.5x
    diabetes_dx %>% 
      filter(grepl('^250\\.5', dx)) %>% 
      group_by(subject_id) %>% 
      summarise('icd9_date' = min(adate)),
    
    
    by = 'subject_id'
  ) %>% 
  mutate('retinopathy' = pmax(cpt_date, icd9_date))

retinopathy_dx <- 
  retinopathy_dx %>% 
  bind_rows(retinopathy_dx2) %>% 
  group_by(subject_id) %>% 
  summarise('retinopathy' = min(retinopathy))

write_parquet(retinopathy_dx, glue('{data_dir}/microvascular_tte/retinopathy_dx.parquet'))


### Nephropathy
## 2 measures of egfr of < 60 
## seperated by at least 90 days
## no measures between of 60+
diagnose_nephropathy <- function(egfr, date) {
  n <- length(date)
  if(n > 2) {
    for(t1 in 1:(n-1)) {
      for(t2 in 2:n) {
        if(as.numeric(date[t2] - date[t1]) >= 90 & egfr[t1] < 60 & egfr[t2] < 60) {
          if(all(egfr[t1:t2] < 60)) {
            return(date[t2])
          }
        }
      }
    }
  }
  return(NA_Date_)
}

nephropathy_dx <- 
  df_kidney %>% 
  group_by(subject_id) %>% 
  group_split() %>% 
  future_map_dfr(~{
    tibble('subject_id' = .x$subject_id[1],
           'nephropathy' = diagnose_nephropathy(.x$eGFR, .x$lab_date))
  },
  .progress = T) %>% 
  filter(!is.na(nephropathy))

write_parquet(nephropathy_dx, glue('{data_dir}/microvascular_tte/nephropathy_dx.parquet'))

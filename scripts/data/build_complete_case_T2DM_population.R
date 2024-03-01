### Build a Population for the notion of complete case which we can use in simulations
###
### To qualify for the population the subject has to 
### (1) Be continuous enrolled for a full year before 2011 (e.g. all of 2010)
### (2) Have some BMI of at least 35 during the year 2011 (we will take the first one)
### (3) Have some HbA1c of at least 6.5% (or equivalent) during 2010/2011 (we will take first one)
###     This is just for diabetes indicator
### (4) Have not previously had surgery prior to 2011
###
### Subsequently, we will simulate trajectories for both BMI and HbA1c

library(tidyverse)
library(glue)
library(haven)
library(lubridate)
library(data.table)
library(furrr)
library(arrow)

source('scripts/util/helpers.R')

n_cores <- 12
plan(future::multisession(workers = n_cores))
options(future.globals.maxSize = 8 * 1024^3)
options(dplyr.summarise.inform = F)

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Read in ID files
cat('Reading in Data\n')
cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/durable_cases.sas7bdat'))
controls <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/raw_race_controls.sas7bdat')) 

### Raw Case and Control Weight measurements
case_weights <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/raw_wts.sas7bdat'))
control_weights <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/raw_wts_controls.sas7bdat'))

case_weights <- 
  case_weights %>% 
  select(measure_date, wt, 'subject_id' = durable_studyid)

control_weights <- 
  control_weights %>% 
  select(measure_date, wt, 'subject_id' = control_studyid)

### Heights (for BMI Calculation) 
case_heights <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/durable_modal_ht.sas7bdat'))
control_heights <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/durable_modal_ht_controls.sas7bdat'))

case_heights <- 
  case_heights %>% 
  select('subject_id' = durable_studyid, 'height' = ModeHeight)

control_heights <- 
  control_heights %>% 
  select('subject_id' = control_studyid, 'height' = ModeHeight)

### Races for cases to join in (control races := ID file) 
case_races <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/raw_race_cases.sas7bdat'))
case_races <- 
  case_races %>% 
  select('subject_id' = durable_studyid, 'race' = race1)

### Build dataset of possible subjects (cases/controls)
cat('Build dataset of possible subjects\n')
df_cases <- 
  cases %>% 
  select('subject_id' = durable_studyid,
         'birth_date' = BIRTH_DATE,
         'site' = site,
         'surg_date' = index_date,
         'gender' = GENDER) %>% 
  left_join(case_races, by = 'subject_id') %>% 
  left_join(case_heights, by = 'subject_id')

df_controls <- 
  controls %>% 
  select('subject_id' = control_studyid,
         'birth_date' = birth_date,
         'site' = site_id,
         'gender' = Gender,
         'race' = race1) %>% 
  left_join(control_heights, by = 'subject_id')

df_subjects <- 
  df_cases %>% 
  bind_rows(df_controls) %>% 
  mutate('case' = !is.na(surg_date))

################################################################################
### (1) Be continuous enrolled for a full year before 2011 (e.g. all of 2010) 
################################################################################

### Enrollment database
### Want only subjects who were continuously enrolled during all of 2010
case_enrollment <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/all_enrollment_collapsed.sas7bdat'))
control_enrollment <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/enrollment_collapsed_controls.sas7bdat'))

case_enrollment <- 
  case_enrollment %>% 
  select('subject_id' = durable_studyid, contains('enr')) %>% 
  filter(enr_start <= '2010-01-01',
         enr_end >= '2011-01-01') 

control_enrollment <- 
  control_enrollment %>% 
  select('subject_id' = control_studyid, contains('enr')) %>% 
  filter(enr_start <= '2010-01-01',
         enr_end >= '2011-01-01') 

df_enrollment <- bind_rows(case_enrollment, control_enrollment)

### Restrict all our files to the set of these subjects continuously enrolled in 2010
cat('Restrict all our files to the set of these subjects continuously enrolled in 2010\n')
df_subjects <- 
  df_subjects %>% 
  inner_join(df_enrollment, by = 'subject_id')

case_weights <- 
  case_weights %>% 
  inner_join(df_enrollment, by = 'subject_id')

control_weights <- 
  control_weights %>% 
  inner_join(df_enrollment, by = 'subject_id')

################################################################################
### (2) Have some BMI of at least 35 during the year 2011 (we will take the first one)
################################################################################

### Smooth weights to find outliers
### Fit loess of weights on date for each subject 
### Remove outlier weights that exceed n_sd SD from smoothed valued
### Since SE can be quite small when we have a lot of measures we keep measures
### within an absolute weight tolerance
### Returns TRUE if weight is outlier
outlier_weight <- function(wt, date, n_sd, weight_tol) {
  ### Don't even bother smoothing on too few obs
  if(length(wt) <= 5) {
    return(rep(F, length(wt)))
  }
  
  ### finer-grained smoothing if we have more
  s <- ifelse(length(wt) > 10, 0.5, 0.75) 
  
  date <- as.numeric(date)
  smoother <- loess(wt ~ date, span = s)
  wt_smooth <- tryCatch(predict(smoother, date, se = T))
  
  if(class(wt_smooth) != 'try-error') {
    upper_bound <- wt_smooth$fit + pmax(n_sd * wt_smooth$se.fit, weight_tol, na.rm = T)
    lower_bound <- wt_smooth$fit - pmax(n_sd * wt_smooth$se.fit, weight_tol, na.rm = T)
    
    ix <- wt < lower_bound | wt > upper_bound
    
  } else {
    return(rep(F, length(wt)))
  }
  
  return(ix)
}

case_subj <- 
  case_weights %>% 
  select(subject_id, wt, measure_date) %>% 
  group_by(subject_id) %>% 
  group_split()

control_subj <- 
  control_weights %>% 
  select(subject_id, wt, measure_date) %>% 
  group_by(subject_id) %>% 
  group_split()

case_outlier_wts <- 
  future_map_dfr(case_subj,  
                 ~{.x %>% 
                     mutate('outlier_wt' = outlier_weight(wt = .x$wt, 
                                                          date = .x$measure_date,
                                                          n_sd = 5,
                                                          weight_tol = 20))
                 })

control_outlier_wts <- 
  future_map_dfr(control_subj, 
                 ~{.x %>% 
                     mutate('outlier_wt' = outlier_weight(wt = .x$wt, 
                                                          date = .x$measure_date,
                                                          n_sd = 5,
                                                          weight_tol = 20))
                 })

outlier_wts <- 
  bind_rows(case_outlier_wts,
            control_outlier_wts)

### Weights in 2011
weights <- 
  bind_rows(case_weights, control_weights) %>% 
  inner_join(outlier_wts, by = c('measure_date', 'wt', 'subject_id')) %>% 
  left_join(select(df_subjects, subject_id, height, birth_date, surg_date), by = 'subject_id') %>% 
  filter(!outlier_wt) %>%  ### Remove Outlier weights
  filter(year(measure_date) == 2011) %>% 
  filter(measure_date <= surg_date | is.na(surg_date)) %>% 
  mutate('bmi' = (wt*0.45359237) / (height*0.0254)^2,
         'age' = as.numeric(measure_date - birth_date)/365.25)

### Just take first measure in 2011 (if any)
baseline_weights <- 
  weights %>% 
  group_by(subject_id) %>% 
  reframe('baseline_bmi' = first(bmi[age <= 80 & age >= 18]),
          'baseline_age' = first(age[age <= 80 & age >= 18]))

### Join in the baseline weights/age to filter the population further.
df_subjects <- 
  df_subjects %>% 
  inner_join(baseline_weights, by = 'subject_id')

################################################################################
### (3) Have some HbA1c of at least 6.5% (or equivalent) during 2010/2011 (we will take first one)
################################################################################
hba1c_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/hba1c_glucose_labs_cases.sas7bdat'))

diabetes_icd9_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/raw_diabetes_dx_cases.sas7bdat'))
diabetes_icd9_controls <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/raw_diabetes_dx_controls.sas7bdat'))
diabetes_rx_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/all_diabetes_rx.sas7bdat'))
diabetes_rx_controls <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/rx_diabetes_controls.sas7bdat'))

diabetes_dx <- function(icd9_code) {
  dx <- 
    case_when(substr(icd9_code, 1, 3) %in% c('250', 'E08', 'E09', 'E10', 'E11', 'E13') ~ T,
              substr(icd9_code, 1, 5) %in% c('362.0', '357.2') ~ T,
              substr(icd9_code, 1, 6) == '366.41' ~ T,
              T ~ F)
  
  return(dx)
}

### Diabetes via DX
diabetes_cases_dx <- 
  diabetes_icd9_cases %>% 
  filter(year(adate) %in% c(2011)) %>% 
  filter(adate <= index_date) %>% ### Need this for cases to remove measures after surgery
  filter(diabetes_dx(dx)) %>% 
  group_by('subject_id' = durable_studyid) %>% 
  summarise('diabetic' = 1)

diabetes_controls_dx <- 
  diabetes_icd9_controls %>% 
  filter(year(adate) %in% c(2011)) %>% 
  filter(diabetes_dx(dx)) %>% 
  group_by('subject_id' = control_studyid) %>% 
  summarise('diabetic' = 1)

### Diabetes via Labs
diabetes_cases_labs <- 
  hba1c_cases %>% 
  filter(year(lab_date) %in% c(2011)) %>% 
  filter(lab_date <= index_date) %>% ### Need this for cases to remove measures after surgery
  filter( (test_type == 'HGBA1C' & result >= 6.5) | (test_type == 'GLU_F' & result >= 126) ) %>% 
  group_by('subject_id' = durable_studyid) %>% 
  summarise('diabetic' = 1)

diabetes_controls_labs <- 
  hba1c_controls %>% 
  filter(year(lab_date) %in% c(2011)) %>% 
  filter( (test_type == 'HGBA1C' & result >= 6.5) | (test_type == 'GLU_F' & result >= 126) ) %>% 
  group_by('subject_id' = control_studyid) %>% 
  summarise('diabetic' = 1)

### Insulin Usage
insulin_cases <- 
  diabetes_rx_cases %>% 
  filter(rxdate <= index_date) %>% 
  filter(rx_year == 2011) %>% 
  group_by('subject_id' = durable_studyid) %>% 
  summarise('insulin' = sum(any(insulin_flg == 1, na.rm = T)))

insulin_controls <- 
  diabetes_rx_controls %>% 
  filter(rx_year == 2011) %>% 
  group_by('subject_id' = control_studyid) %>% 
  summarise('insulin' = sum(any(insulin_flg == 1, na.rm = T)))

diabetes <- 
  bind_rows(diabetes_cases_dx, diabetes_controls_dx, 
            diabetes_cases_labs, diabetes_controls_labs) %>% 
  distinct() %>% 
  left_join(bind_rows(insulin_cases, insulin_controls), by = 'subject_id') %>% 
  mutate('insulin' = ifelse(is.na(insulin), 0, insulin))  

df_subjects <- 
  df_subjects %>% 
  left_join(diabetes, by = 'subject_id') %>% 
  mutate('insulin' = ifelse(is.na(insulin), 0, insulin),
         'diabetic' = ifelse(is.na(diabetic), 0, diabetic))

### "Baseline HGBA1C" 
df_baseline_a1c_cases <- 
  hba1c_cases %>% 
  filter(year(lab_date) %in% c(2011)) %>% 
  filter(lab_date <= index_date) %>% ### Need this for cases to remove measures after surgery
  filter(test_type != 'GLU_RAN') %>%
  filter(result_unit %in% c('%', 'mg/dl', 'mg/dL', 'MG/DL')) %>% 
  mutate('subject_id' = durable_studyid, 
         'hgba1c' = case_when(test_type == 'HGBA1C' ~ result,
                              test_type == 'GLU_F' ~ (result + 46.7)/28.7)) %>% 
  filter(hgba1c < 20) %>% 
  group_by(subject_id) %>% 
  summarise('baseline_hgba1c' = first(hgba1c)) %>% 
  select(subject_id, baseline_hgba1c)

df_baseline_a1c_controls <- 
  hba1c_controls %>% 
  filter(year(lab_date) %in% c(2011)) %>% 
  filter(test_type != 'GLU_RAN') %>%
  filter(result_unit %in% c('%', 'mg/dl', 'mg/dL', 'MG/DL')) %>% 
  mutate('subject_id' = control_studyid, 
         'hgba1c' = case_when(test_type == 'HGBA1C' ~ result,
                              test_type == 'GLU_F' ~ (result + 46.7)/28.7)) %>% 
  filter(hgba1c < 20) %>% 
  group_by(subject_id) %>% 
  summarise('baseline_hgba1c' = first(hgba1c)) %>% 
  select(subject_id, baseline_hgba1c)

df_subjects <- 
  df_subjects %>% 
  inner_join(bind_rows(df_baseline_a1c_cases, df_baseline_a1c_controls),
             by = 'subject_id')

################################################################################
### (4) Have not previously had surgery prior to 2011
################################################################################
df_subjects <- 
  df_subjects %>% 
  filter(is.na(surg_date) | year(surg_date) < 2011)

##########################################
### Smoking Status (at some point in 2011)
##########################################
smoking_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/smoking_vdw.sas7bdat'))
smoking_controls <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/smoking_vdw_controls.sas7bdat'))

df_smoking <- 
  bind_rows(
    ### Cases
    smoking_cases %>% 
      filter(contact_date <= index_date) %>% 
      filter(year(contact_date) == 2011) %>% 
      mutate('smoking_status' = case_when(tobacco_use == 'Y' ~ 'current',
                                          tobacco_use == 'Q' ~ 'former',
                                          tobacco_use %in% c('N', 'P') ~ 'never',
                                          T ~ 'no_self_report')) %>% 
      select('subject_id' = durable_studyid, contact_date, smoking_status) %>% 
      group_by(subject_id) %>% 
      summarise('smoking_status' = first(smoking_status)),
    
    ### Controls
    smoking_controls %>% 
      filter(year(contact_date) == 2011) %>% 
      mutate('smoking_status' = case_when(tobacco_use == 'Y' ~ 'current',
                                          tobacco_use == 'Q' ~ 'former',
                                          tobacco_use %in% c('N', 'P') ~ 'never',
                                          T ~ 'no_self_report')) %>% 
      select('subject_id' = control_studyid, contact_date, smoking_status) %>% 
      group_by(subject_id) %>% 
      summarise('smoking_status' = first(smoking_status))
  )

df_subjects <- 
  df_subjects %>% 
  left_join(df_smoking, by = 'subject_id') %>% 
  mutate('smoking_status' = ifelse(is.na(smoking_status), 'no_self_report', smoking_status))

### Clean Up Measures
df_subjects <- 
  df_subjects %>% 
  mutate('baseline_bmi' = winsorize(baseline_bmi, q = c(0.001, 0.999)),
         'baseline_hgba1c' = winsorize(baseline_hgba1c, q = c(0.001, 0.999)))

### Comorbidities
elix_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/final_elixhauser_score_cases.sas7bdat'))

### Use distribution of computed case scores to inform random sampling on controls
### Have to due this since scores aren't computed for the controls 
### Only have raw comorbidity diagnosis codes
elix_dist <- 
  elix_cases %>% 
  mutate('elix_score' = case_when(combinedscore < 0 ~ -1,
                                  combinedscore == 0 ~ 0,
                                  combinedscore == 1 ~ 1,
                                  combinedscore >= 2 ~ 2)) %>% 
  group_by(elix_score) %>% 
  count() %>% 
  ungroup() %>% 
  mutate('freq' = n/sum(n))

set.seed(132)

df_subjects  <- 
  df_subjects %>% 
  mutate('elix_score' = sample(x = elix_dist$elix_score, 
                               size = nrow(.), 
                               replace = T, 
                               prob = elix_dist$freq))

### Save file
write_parquet(df_subjects, glue('{data_dir}/t2dm_population.parquet'))

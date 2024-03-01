### Script to clean all weight measurements
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
options(future.globals.maxSize = 80 * 1024^3)
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
         'index_date' = index_date,
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
  mutate('case' = !is.na(index_date))


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
                 },
                 .progress = T)

control_outlier_wts <- 
  future_map_dfr(control_subj, 
                 ~{.x %>% 
                     mutate('outlier_wt' = outlier_weight(wt = .x$wt, 
                                                          date = .x$measure_date,
                                                          n_sd = 5,
                                                          weight_tol = 20))
                 },
                 .progress = T)

outlier_wts <- 
  bind_rows(case_outlier_wts,
            control_outlier_wts)

weights <- 
  bind_rows(case_weights, control_weights) %>% 
  inner_join(outlier_wts, by = c('measure_date', 'wt', 'subject_id')) %>% 
  left_join(select(df_subjects, subject_id, height, birth_date, index_date), by = 'subject_id') %>% 
  filter(!outlier_wt) %>%  ### Remove Outlier weights
  mutate('bmi' = (wt*0.45359237) / (height*0.0254)^2,
         'age' = as.numeric(measure_date - birth_date)/365.25)

write_parquet(weights, glue('{data_dir}/all_weights.parquet'))
write_parquet(outlier_wts, glue('{data_dir}/all_outlier_weight_tags.parquet'))

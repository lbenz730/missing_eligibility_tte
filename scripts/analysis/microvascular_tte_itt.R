library(tidyverse)
library(lubridate)
library(arrow)
library(glue)
library(data.table)

source('scripts/util/helpers.R')
source('scripts/simulations/helpers.R')

options(dplyr.summarise.inform = F)

library(furrr)
n_cores <- 12
plan(future::multisession(workers = n_cores))
options(future.globals.maxSize = 80 * 1024^3)

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Load in all data
### Demographics, labs etc.
df_subjects <- read_parquet(glue('{data_dir}/microvascular_tte/subjects.parquet'))
weights <- read_parquet(glue('{data_dir}/all_weights.parquet')) ### Cleaned Weights w/ outliers removed
df_pregnancy <- read_parquet(glue('{data_dir}/microvascular_tte/pregnancy.parquet'))
df_enrollment <- read_parquet(glue('{data_dir}/microvascular_tte/enrollment.parquet'))
df_kidney <- read_parquet(glue('{data_dir}/microvascular_tte/kidney_labs.parquet'))
diabetes_rx <- read_parquet(glue('{data_dir}/microvascular_tte/diabetes_rx.parquet'))
diabetes_labs <- read_parquet(glue('{data_dir}/microvascular_tte/diabetes_labs.parquet'))
smoking <- read_parquet(glue('{data_dir}/microvascular_tte/smoking.parquet'))

### Drop measures on the same date
weights <-
  weights %>% 
  distinct(subject_id, measure_date, .keep_all = T)

smoking <- 
  smoking %>% 
  distinct(subject_id, contact_date, .keep_all = T)

df_kidney <- 
  df_kidney %>% 
  distinct(subject_id, lab_date, .keep_all = T)

diabetes_labs <- 
  diabetes_labs %>% 
  distinct(subject_id, lab_date, .keep_all = T)

### Censoring
censor_measurements <- read_parquet(glue('{data_dir}/microvascular_tte/censor_measurements.parquet'))
df_cancer <- read_parquet(glue('{data_dir}/microvascular_tte/cancer.parquet'))
df_death <- read_parquet(glue('{data_dir}/microvascular_tte/death.parquet'))

### Outcomes
neuropathy_dx <- read_parquet(glue('{data_dir}/microvascular_tte/neuropathy_dx.parquet'))
nephropathy_dx <- read_parquet(glue('{data_dir}/microvascular_tte/nephropathy_dx.parquet'))
retinopathy_dx <- read_parquet(glue('{data_dir}/microvascular_tte/retinopathy_dx.parquet'))


### Function to build a trial (pre-discrete time)
build_trial <- function(trial_id) {
  options(dplyr.summarise.inform = F)
  
  trial_start <- study_start %m+% months(trial_id - 1)
  trial_end <- study_start %m+% months(trial_id) - 1
  bmi_start <- trial_start %m-% months(bmi_lookback)
  a1c_start <- trial_start %m-% months(a1c_lookback)
  
  
  ### Pregnancies to exclude in the last year
  pregancy_1yr <- 
    df_pregnancy %>% 
    filter(adate <= trial_start,
           adate >= trial_start %m-% years(1)) %>% 
    distinct(subject_id)
  
  
  trial_population <-
    df_subjects %>% 
    
    ### Remove old cases and code surgery in momth
    filter(is.na(index_date) | index_date >= trial_start) %>% 
    mutate('surgery' = as.numeric(!is.na(index_date) & index_date >= trial_start & index_date <= trial_end)) %>% 
    
    ### Outcome dates 
    left_join(neuropathy_dx, by = 'subject_id') %>% 
    left_join(retinopathy_dx, by = 'subject_id') %>% 
    left_join(nephropathy_dx, by = 'subject_id') %>% 
    
    ### Censoring
    left_join(censor_measurements, by = 'subject_id') %>% 
    left_join(df_cancer, by = 'subject_id') %>% 
    left_join(df_death, by = 'subject_id') %>% 
    mutate('study_end' = study_end) %>% 
    anti_join(pregancy_1yr, by = 'subject_id') %>% 
    
    ### Enrollment 
    inner_join(df_enrollment %>% 
                 filter(enr_1yr <= trial_start, enr_end >= trial_start),
               by = 'subject_id') %>% 
    
    ### Times
    mutate(
      'censor_time' = pmin(enr_end, 
                           censor_death, 
                           censor_cancer, 
                           censor_measurement,
                           study_end, 
                           na.rm = T),
      
      'outcome_time' = pmin(neuropathy,
                            nephropathy,
                            retinopathy,
                            na.rm = T),
      
      'non_adhere_time' = case_when(surgery == 0 ~ index_date,
                                    T ~ NA_Date_)
      
    ) %>% 
    
    ### Remove pre-existing outcomes
    filter(outcome_time >= trial_start | is.na(outcome_time),
           censor_time >= trial_start | is.na(censor_time)) %>% 
    
    mutate(
      ### Event times (up to nearest month)
      'time_combined' = ceiling(as.numeric(outcome_time - trial_start)/30.25),
      'time_neuropathy' = ceiling(as.numeric(neuropathy - trial_start)/30.25),
      'time_nephropathy' = ceiling(as.numeric(nephropathy - trial_start)/30.25),
      'time_retinopathy' = ceiling(as.numeric(retinopathy - trial_start)/30.25),
      'time_censor' = ceiling(as.numeric(censor_time - trial_start)/30.25),
      'time_non_adhere' = ceiling(as.numeric(non_adhere_time - trial_start)/30.25),
    )
  
  
  ### Determine Eligibility
  df_bmi <- 
    weights %>% 
    filter(measure_date <= trial_end,
           measure_date >= bmi_start) %>% 
    arrange(measure_date) %>% 
    group_by(subject_id) %>% 
    summarise('baseline_bmi' = last(bmi)) %>% 
    ungroup()
  
  df_a1c <- 
    diabetes_labs %>% 
    filter(lab_date <= trial_end, 
           lab_date >= a1c_start) %>% 
    arrange(lab_date) %>%  
    group_by(subject_id, test_type) %>% 
    summarise('baseline_value' = last(result)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = test_type,
                values_from = baseline_value) %>% 
    rename('baseline_gluF' = GLU_F,
           'baseline_a1c' = HGBA1C) %>% 
    mutate('baseline_a1c' = ifelse(is.na(baseline_a1c), (baseline_gluF + 46.7)/28.7, baseline_a1c),
           'baseline_gluF' = ifelse(is.na(baseline_gluF), baseline_a1c * 28.7 - 46.7, baseline_gluF))
  
  df_rx <- 
    diabetes_rx %>% 
    filter(rxdate <= trial_start, rx_end_date >= trial_start) %>% 
    select(subject_id, insulin_flg) %>% 
    group_by(subject_id) %>% 
    summarise('insulin_flg' = any(insulin_flg == 1, na.rm = T)) %>% 
    mutate('diabetes_rx' = 1) %>% 
    ungroup()
  
  ### Most recent Smoking Status
  df_smoking <- 
    smoking %>% 
    filter(contact_date <= study_start) %>% 
    arrange(desc(contact_date)) %>% 
    group_by(subject_id)  %>% 
    slice(1) %>% 
    ungroup()
  
  ### Most Recent Kidney Labs
  df_eGFR <- 
    df_kidney %>% 
    filter(lab_date <= study_start) %>% 
    arrange(desc(lab_date)) %>% 
    group_by(subject_id) %>% 
    slice(1) %>% 
    ungroup() %>% 
    rename('baseline_scr' = scr,
           'baseline_eGFR' = eGFR)
  
  ### Dataset for trial (before expansion)
  df_trial <- 
    trial_population %>% 
    mutate('trial_id' = trial_id,
           'trial_start' = trial_start) %>% 
    mutate('baseline_age' = as.numeric(trial_start - birth_date)/365.25) %>% 
    left_join(df_bmi, by = 'subject_id') %>% 
    left_join(df_a1c, by = 'subject_id') %>% 
    left_join(df_rx, by = 'subject_id') %>% 
    left_join(df_smoking, by = 'subject_id') %>% 
    left_join(df_eGFR, by = 'subject_id') %>% 
    mutate('diabetes_rx' = ifelse(is.na(diabetes_rx), 0, diabetes_rx),
           'insulin_flg' = ifelse(is.na(insulin_flg), 0, insulin_flg)) %>% 
    mutate('eligible' = 
             case_when(
               ### Eligible
               surgery == 1 & (baseline_gluF >= 126 | baseline_a1c >= 6.5 | diabetes_rx == 1) ~ 1,
               surgery == 0 & baseline_bmi >= 35 & (baseline_gluF >= 126 | baseline_a1c >= 6.5 | diabetes_rx == 1) ~ 1,
               
               ### Not Eligible
               surgery == 1 & !(baseline_gluF >= 126 | baseline_a1c >= 6.5 | diabetes_rx == 1) ~ 0,
               surgery == 0 & !(baseline_bmi >= 35 & (baseline_gluF >= 126 | baseline_a1c >= 6.5 | diabetes_rx == 1)) ~ 0,
               
               ### Missing
               T ~ NA_integer_
             )
    ) %>% 
    mutate('R' = as.numeric(!is.na(eligible)))
  
  
  ### Light Data Cleaning
  df_trial <-
    df_trial %>% 
    mutate('race' = ifelse(race %in% c('MU', 'OT', 'UN'), 'OT', race))
  
  
  return(df_trial)
}

### Combinations to consider
df_files <- 
  crossing('bmi_lookback' = c(1,3,6,12),
           'a1c_lookback' = c(1,3,6,12,18,24))
args <- commandArgs(trailingOnly = T)
file_id <- as.numeric(args[1])

### Encode constants
study_start <- as.Date('2005-01-01')
study_end <- as.Date('2015-09-30')
n_trials <- 84
bmi_lookback <- df_files$bmi_lookback[file_id]
a1c_lookback <- df_files$a1c_lookback[file_id]
set.seed(73097)

### Build the starting dataset for all trials
df_trials <-
  map_dfr(1:n_trials, ~{
    cat('Building Trial:', .x, '\n')
    build_trial(.x)
  })

### "Impute" Smoking Status
df_trials <-
  df_trials %>%
  mutate('smoking_status' = ifelse(is.na(smoking_status), 'no_self_report', smoking_status))

### Impute Baseline SCR and compute eGFR off of it
ckd_epi <- function(scr, age, sex) {
  scr <- pmax(0.1, scr)
  alpha <- ifelse(sex == 'M', -0.302, -0.241)
  kappa <- ifelse(sex == 'M', 0.9, 0.7)

  eGFR <- 142 * pmin(scr/kappa, 1)^alpha * pmax(scr/kappa, 1)^(-1.2) *  0.9938^age * ifelse(sex == 'M', 1, 1.012)
  return(eGFR)
}

scr_model <-
  glm(baseline_scr ~ gender + surgery + race + site + baseline_age + insulin_flg,
      family = Gamma(link = "log"),
      data = df_trials)
shape <- 1/summary(scr_model)$dispersion
scr_imputed <-
  rgamma(n = sum(is.na(df_trials$baseline_scr)),
         shape = shape,
         rate = shape/exp(predict(scr_model, newdata = filter(df_trials, is.na(baseline_scr)))))
ix_na <- is.na(df_trials$baseline_scr)
df_trials$baseline_scr[ix_na] <- scr_imputed
df_trials$baseline_eGFR[ix_na] <-
  ckd_epi(scr = df_trials$baseline_scr[ix_na],
          age = df_trials$baseline_age[ix_na],
          sex = df_trials$gender[ix_na])

### Clear uneeded files
rm(censor_measurements)
rm(df_death)
rm(df_enrollment)
rm(df_pregnancy)
rm(df_subjects)
rm(diabetes_rx)
rm(nephropathy_dx)
rm(neuropathy_dx)
rm(retinopathy_dx)
garbage <- gc()

### Save out file for future usage
if(!dir.exists(glue('{data_dir}/microvascular_tte/pooled_trials'))) {
  dir.create(glue('{data_dir}/microvascular_tte/pooled_trials'))
}
write_parquet(df_trials, glue('{data_dir}/microvascular_tte/pooled_trials/trials_{file_id}.parquet'))
df_trials <- read_parquet(glue('{data_dir}/microvascular_tte/pooled_trials/trials_{file_id}.parquet'))

### Analyze Results
### Weights for Missing Data (IPW R)
cat('Fitting Missingness Model\n')
missingness_model <-
  glm(R ~ site + baseline_age + race + gender + surgery + smoking_status + baseline_eGFR,
      data = select(df_trials, R, site, baseline_age, race, gender, surgery, smoking_status, baseline_eGFR),
      family = 'binomial')

df_analysis <- 
  df_trials %>% 
  mutate('ipw_R' = mean(df_trials$R)/missingness_model$fitted.values) %>%
  mutate('ipw_R' = winsorize(ipw_R, q = c(0, 0.99))) %>%
  filter(eligible == 1)

### Delete large model object
missingness_model <- strip_model(missingness_model)
rm(missingness_model)
rm(df_trials) 
garbage <- gc()

### Impute Missing Baseline Covariates
### BMI
cat('Imputing Baseline Covariates\n')
bmi_model <- 
  glm(baseline_bmi ~ site + baseline_age + race + gender + insulin_flg + surgery + smoking_status + baseline_eGFR,
      family = Gamma(link = "log"),
      data = df_analysis) 

shape <- 1/summary(bmi_model)$dispersion
bmi_imputed <- 
  rgamma(n = nrow(df_analysis),
         shape = shape,
         rate = shape/exp(predict(bmi_model, newdata = df_analysis)))
ix_na <- is.na(df_analysis$baseline_bmi)
df_analysis$baseline_bmi[ix_na] <- bmi_imputed[ix_na]

### A1c (add back in site)
a1c_model <- 
  glm(baseline_a1c ~ site + baseline_age + race + gender + baseline_bmi + insulin_flg + surgery + smoking_status + baseline_eGFR,
      family = Gamma(link = "log"),
      data = df_analysis) 

shape <- 1/summary(a1c_model)$dispersion
a1c_imputed <- 
  rgamma(n = nrow(df_analysis),
         shape = shape,
         rate = shape/exp(predict(a1c_model, newdata = df_analysis)))


ix_na <- is.na(df_analysis$baseline_a1c)
df_analysis$baseline_a1c[ix_na] <- a1c_imputed[ix_na]

### Delete large model objects
rm(bmi_model)
rm(a1c_model)
garbage <- gc()

### Weights for Confounding (IPW A)
cat('Fitting Treatment Model\n')
treatment_model <- 
  glm(surgery ~ site + baseline_age + race + gender + insulin_flg + baseline_bmi + baseline_a1c + smoking_status + baseline_eGFR,
      data = df_analysis,
      family = 'binomial')

p_surgery <- treatment_model$fitted.values
df_analysis <- 
  df_analysis %>% 
  mutate('ipw_A' = ifelse(surgery == 1, mean(surgery)/p_surgery, mean(1-surgery)/(1-p_surgery))) %>% 
  group_by(surgery) %>% 
  mutate('ipw_A' = winsorize(x = ipw_A, q = c(0, 0.99))) %>% 
  ungroup()

### Make time discrete
covariates <- 
  c('gender', 'race', 'site', 'insulin_flg', 'smoking_status', 'age',
    'baseline_bmi', 'baseline_a1c', 'baseline_age', 'baseline_eGFR',
    'bmi_tv', 'a1c_tv', 'eGFR_tv', 'smoking_status_tv')

outcomes <- c('time_combined', 'time_retinopathy', 'time_neuropathy', 'time_nephropathy')
df_results <- NULL
for(outcome in outcomes) {
  cat('Computing models for:', outcome, '\n')
  df_analysis$time_outcome <- df_analysis[[outcome]]
  
  df_times <-   
    df_analysis %>% 
    select(trial_id, subject_id, starts_with('time')) %>% 
    group_by(trial_id, subject_id) %>% 
    reframe('time' = 1:pmin(time_outcome, time_censor, na.rm = T),
            'censor' = as.numeric(!is.na(time_censor) & time == time_censor),
            'outcome' = as.numeric(!is.na(time_outcome) & time == time_outcome)) %>% 
    ungroup()
  
  df_full <- 
    df_analysis %>% 
    select(trial_id, subject_id, surgery, any_of(covariates), ipw_R, ipw_A) %>% 
    inner_join(df_times, by = c('subject_id', 'trial_id')) %>% 
    mutate('actual_time' = study_start %m+% months(trial_id - 1 + time - 1)) %>% 
    mutate('age' = baseline_age + (time - 1)/12)
  
  ### Fit Final Model
  ### Note that because I am fitting a saturated (weighted) model I can use my 
  ### glm_quick formula for a closed form solution
  model_unweighted <- 
    glm_quick(Y = df_full$outcome,
              A = df_full$surgery,
              W = rep(1, nrow(df_full)))
  model_A <-
    glm_quick(Y = df_full$outcome,
              A = df_full$surgery,
              W = df_full$ipw_A)
  model_AR <-
    glm_quick(Y = df_full$outcome,
              A = df_full$surgery,
              W = df_full$ipw_A * df_full$ipw_R)

  df_tmp <- 
    tibble('outcome' = case_when(outcome == 'time_combined' ~ 'Any Microvascular Event',
                                 outcome == 'time_retinopathy' ~ 'Retinopathy',
                                 outcome == 'time_neuropathy' ~ 'Neuropathy',
                                 outcome == 'time_nephropathy' ~ 'Nephropathy'),
           'weights' = c('Unweighted', 'W^A', 'W^A x W^R'),
           'estimate' = c(model_unweighted$beta1,
                          model_A$beta1,
                          model_AR$beta1)) %>% 
    mutate('exp_est' = exp(estimate),
           'effect' = 'Intent-To-Treat')
  
  df_results <- bind_rows(df_results, df_tmp)
}

### Save results
if(!dir.exists('data/microvascular_results')) {
  dir.create('data/microvascular_results')
}

df_results <- 
  df_results %>% 
  mutate('file_id' = file_id,
         'bmi_lookback' = bmi_lookback,
         'a1c_lookback' = a1c_lookback)

write_csv(df_results, glue('data/microvascular_results/file_{file_id}_itt.csv'))
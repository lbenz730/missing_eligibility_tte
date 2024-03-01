library(tidyverse)
library(lubridate)
library(arrow)
library(glue)
library(data.table)

source('scripts/util/helpers.R')
source('scripts/simulations/helpers.R')

options(dplyr.summarise.inform = F)

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'


### Load in all data
### Demographics, labs etc.
weights <- read_parquet(glue('{data_dir}/all_weights.parquet')) ### Cleaned Weights w/ outliers removed
df_kidney <- read_parquet(glue('{data_dir}/microvascular_tte/kidney_labs.parquet'))
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

### Combinations to consider
df_settings <- 
  crossing('bmi_lookback' = c(1,3,6,12),
           'a1c_lookback' = c(1,3,6,12,18,24)) %>% 
  mutate('setting_id' = 1:nrow(.))

args <- commandArgs(trailingOnly = T)
file_id <- as.numeric(args[1])


### Encode constants
n_boot <- 1000
df_files <- 
  crossing('boot_id' = 1:n_boot,
           'outcome_id' = 1:4)
study_start <- as.Date('2005-01-01')
study_end <- as.Date('2015-09-30')
n_trials <- 84
bmi_lookback <- 3
a1c_lookback <- 12
setting_id <- df_settings$setting_id[df_settings$bmi_lookback == bmi_lookback & df_settings$a1c_lookback == a1c_lookback]
boot_id <- df_files$boot_id[file_id]
outcome_id <- df_files$outcome_id[file_id]

### Set Seed
set.seed(73097)
seeds <- sample(1:1000000, n_boot, replace = F)
set.seed(seeds[boot_id])

### Boot-strap 
df_trials <- read_parquet(glue('{data_dir}/microvascular_tte/pooled_trials/trials_{setting_id}.parquet'))
df_trials <- as.data.table(df_trials)
subject_ids <- unique(df_trials$subject_id)
ix_subj <- sample(1:length(subject_ids), size = length(subject_ids), replace = T)
dt_subj <- data.table('subject_id' = subject_ids[ix_subj]) 
df_sample <- df_trials[dt_subj, on = 'subject_id']
df_sample[, 'replicate_id' := 1:.N, by = c('subject_id', 'trial_id')]
rm(df_trials)
garbage <- gc()

### Analyze Results
### Weights for Missing Data (IPW R)
cat('Fitting Missingness Model\n')
missingness_model <-
  glm(R ~ site + baseline_age + race + gender + surgery + smoking_status + baseline_eGFR,
      data = select(df_sample, R, site, baseline_age, race, gender, surgery, smoking_status, baseline_eGFR),
      family = 'binomial')

df_analysis <- 
  df_sample %>% 
  mutate('ipw_R' = mean(df_sample$R)/missingness_model$fitted.values) %>%
  mutate('ipw_R' = winsorize(ipw_R, q = c(0, 0.99))) %>%
  filter(eligible == 1)

### Delete large model object
missingness_model <- strip_model(missingness_model)
rm(missingness_model)
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

### Delete large model objects
rm(bmi_model)
rm(a1c_model)
rm(treatment_model)
garbage <- gc()

### Make time discrete
covariates <- 
  c('gender', 'race', 'site', 'insulin_flg', 'smoking_status', 'age',
    'baseline_bmi', 'baseline_a1c', 'baseline_age', 'baseline_eGFR',
    'bmi_tv', 'a1c_tv', 'eGFR_tv', 'smoking_status_tv')

outcomes <- c('time_combined', 'time_retinopathy', 'time_neuropathy', 'time_nephropathy')
df_results <- NULL
outcome <- outcomes[outcome_id]

cat('Computing models for:', outcome, '\n')
df_analysis$time_outcome <- df_analysis[[outcome]]

df_times <-   
  df_analysis %>% 
  select(trial_id, subject_id, replicate_id, starts_with('time')) %>% 
  group_by(trial_id, subject_id, replicate_id) %>% 
  reframe('time' = 1:pmin(time_outcome, time_censor, time_non_adhere, na.rm = T),
          'censor' = as.numeric(!is.na(time_censor) & time == time_censor),
          'non_adhere' = as.numeric(!is.na(time_non_adhere) & time == time_non_adhere),
          'outcome' = as.numeric(!is.na(time_outcome) & time == time_outcome)) %>% 
  ungroup()

df_full <- 
  df_analysis %>% 
  select(trial_id, subject_id, replicate_id, surgery, any_of(covariates), ipw_R, ipw_A) %>% 
  inner_join(df_times, by = c('subject_id', 'trial_id', 'replicate_id')) %>% 
  mutate('actual_time' = study_start %m+% months(trial_id - 1 + time - 1)) %>% 
  mutate('age' = baseline_age + (time - 1)/12)

## Join in most recent (A1c, BMI, smoking_status, eGRF) measures
df_full <-
  df_full %>% 
  left_join(smoking, 
            join_by(subject_id, closest(actual_time >= contact_date)), 
            suffix = c('', '_tv')) %>% 
  left_join(weights %>% select(subject_id, 'bmi_tv' = bmi, measure_date), 
            join_by(subject_id, closest(actual_time >= measure_date))) %>% 
  left_join(diabetes_labs %>% 
              mutate('a1c_tv' = ifelse(test_type == 'HGBA1C', result, (result + 46.7)/28.7 )) %>% 
              select(subject_id, 'a1c_tv' = result, lab_date), 
            join_by(subject_id, closest(actual_time >= lab_date))) %>% 
  left_join(df_kidney %>% 
              select(subject_id, 'eGFR_tv' = eGFR, lab_date),
            join_by(subject_id, closest(actual_time >= lab_date))) %>% 
  mutate('bmi_tv' = ifelse(time == 1, baseline_bmi, bmi_tv),
         'a1c_tv' = ifelse(time == 1, baseline_a1c, a1c_tv),
         'eGFR_tv' = ifelse(time == 1, baseline_eGFR, eGFR_tv)) %>% 
  mutate('bmi_tv' = ifelse(is.na(bmi_tv), baseline_bmi, bmi_tv),
         'a1c_tv' = ifelse(is.na(a1c_tv), baseline_a1c, a1c_tv),
         'eGFR_tv' = ifelse(is.na(eGFR_tv), baseline_eGFR, eGFR_tv),
         'smoking_status_tv' = ifelse(is.na(smoking_status_tv), smoking_status, smoking_status_tv))


### Non-adherence weights
model_nonadhere <- 
  glm(non_adhere ~ 
        site + race + gender + insulin_flg + 
        baseline_bmi + baseline_a1c + baseline_eGFR + 
        age + bmi_tv + a1c_tv + eGFR_tv + smoking_status_tv,
      family = 'binomial',
      data = 
        df_full %>% 
        filter(surgery == 0, time > 1) %>% 
        select(non_adhere, any_of(covariates))
  )

df_full <- 
  df_full %>% 
  mutate('prob_nonadhere' = 0,
         'prob_nonadhere' = replace(prob_nonadhere, surgery == 0 & time > 1, model_nonadhere$fitted.values),
         'prob_nonadhere_sw' = mean(df_full$non_adhere[df_full$time > 1 & df_full$surgery == 0]),
         'prob_nonadhere_sw' = replace(prob_nonadhere_sw, surgery == 1 | time == 1, 0)) %>% 
  group_by(trial_id, subject_id) %>%
  arrange(time) %>% 
  mutate('cum_prob_adhere' = cumprod(1 - prob_nonadhere)) %>%
  mutate('cum_prob_adhere_sw' = cumprod(1 - prob_nonadhere_sw)) %>%
  ungroup() %>%
  mutate('ipw_N' = winsorize(cum_prob_adhere_sw/cum_prob_adhere, q = c(0, 0.99)))

rm(model_nonadhere)
gc()

### Fit Final Model
### Note that because I am fitting a saturated (weighted) model I can use my 
### glm_quick formula for a closed form solution
model_ANR <-
  glm_quick(Y = df_full$outcome,
            A = df_full$surgery,
            W = df_full$ipw_A * df_full$ipw_N * df_full$ipw_R)

df_tmp <- 
  tibble('outcome' = case_when(outcome == 'time_combined' ~ 'Any Microvascular Event',
                               outcome == 'time_retinopathy' ~ 'Retinopathy',
                               outcome == 'time_neuropathy' ~ 'Neuropathy',
                               outcome == 'time_nephropathy' ~ 'Nephropathy'),
         'weights' = 'W^A x W^N x W^R',
         'estimate' = model_ANR$beta1) %>% 
  mutate('exp_est' = exp(estimate),
         'effect' = 'Per-Protocol')

df_results <- bind_rows(df_results, df_tmp)


### Save results
if(!dir.exists('data/final_pp_bootstraps')) {
  dir.create('data/final_pp_bootstraps')
}

df_results <- 
  df_results %>% 
  mutate('file_id' = file_id,
         'bmi_lookback' = bmi_lookback,
         'a1c_lookback' = a1c_lookback)

write_csv(df_results, glue('data/final_pp_bootstraps/file_{file_id}.csv'))
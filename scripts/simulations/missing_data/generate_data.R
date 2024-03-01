library(tidyverse)
library(arrow)
library(data.table)
source('scripts/simulations/missing_data/weight_trajectories_functions.R')

generate_data <- function(params) {
  ### Unpack list
  n_subjects <- params$n_subjects
  n_trials <- params$n_trials
  n_sims <- params$n_sims
  study_duration <- params$study_duration
  weight_models <- params$weight_models
  a1c_models <- params$a1c_models
  treatment_model <- params$treatment_model
  bs_type_model <- params$bs_type_model
  outcome_model <- params$outcome_model
  elig_criteria <- params$eligibility
  complete_case_elig_model <- params$cc_model
  
  subj_ix <- sample(df_pop$subject_id, n_subjects)
  
  df_subjects <- 
    df_pop %>% 
    filter(subject_id %in% subj_ix)
  
  ### (1) Generate a control weight and a1c trajectory for each subject 
  df_weight_traj_control <- 
    sample_weight_trajectories(df = df_subjects, 
                               coeff = weight_models, 
                               surgery = F) 
  
  df_a1c_traj_control <- 
    sample_a1c_trajectories(df = df_subjects, 
                            coeff = a1c_models, 
                            surgery = F) 
  
  df_weight_query_control <- 
    query_trajectories(df_trajectory = df_weight_traj_control,
                       time_points = 0:study_duration,
                       error = weight_models$control$errors,
                       surgery = F,
                       a1c = F)
  
  df_a1c_query_control <- 
    query_trajectories(df_trajectory = df_a1c_traj_control,
                       time_points = 0:study_duration,
                       error = a1c_models$control$errors,
                       surgery = F,
                       a1c = T,
                       baseline_vals = df_subjects$baseline_hgba1c)
  
  df_query_control <- 
    df_weight_query_control %>% 
    select(subject_id, time, 'ptwc' = observed) %>% 
    inner_join(select(df_a1c_query_control, subject_id, time, 'pct_a1c' = observed), 
               by = c('subject_id', 'time'))
  
  ### (2) Update values of time varying covariates
  df_control <- 
    df_subjects %>% 
    inner_join(df_query_control, by = 'subject_id') %>% 
    mutate('bmi' = baseline_bmi * (1 + ptwc),
           'hgba1c' = baseline_hgba1c * (1 + pct_a1c),
           'age' = baseline_age + time/12) %>% 
    mutate('surgery' = 0)
  
  ### (3) Apply treatment model to determine if/when subject gets treatment
  ### (bmi >= 35) required to get surgery
  surgery_times <- 
    df_control %>% 
    mutate('treatment_prob' = (bmi >= 35) * expit(compute_model(df = df_control, beta = treatment_model)),
           'rygb_prob' = (bmi >= 35) * expit(compute_model(df = df_control, beta = bs_type_model))) %>% 
    mutate('surgery' = rbinom(n = nrow(.), size = 1, p = treatment_prob),
           'bs_type' = rbinom(n = nrow(.), size = 1, p = rygb_prob)) %>% 
    mutate('bs_type' = case_when(surgery == 0 ~ 'control',
                                 bs_type == 1 ~ 'rygb',
                                 bs_type == 0 ~ 'sleeve')) %>% 
    group_by(subject_id) %>% 
    summarise('treatment_time' = dplyr::first(time[surgery == 1]),
              'bs_type' = dplyr::first(bs_type[surgery == 1], default = 'control'))
  
  ### (4) For people who do get surgery, add in their trajectory at the point they get surgery
  surgery_ids <- 
    surgery_times %>% 
    filter(!is.na(treatment_time)) %>% 
    pull(subject_id)
  
  df_surg_subjects <- 
    df_control %>% 
    inner_join(surgery_times, by = 'subject_id') %>% 
    filter(time == treatment_time) %>% 
    ### Get current values of these covariates
    mutate('baseline_bmi' = bmi,
           'baseline_age' = age,
           'baseline_hgba1c' = hgba1c)
  
  df_weight_traj_surgery <- 
    sample_weight_trajectories(df = df_surg_subjects, 
                               coeff = weight_models, 
                               surgery = T)
  df_a1c_traj_surgery <- 
    sample_a1c_trajectories(df = df_surg_subjects, 
                            coeff = a1c_models, 
                            surgery = T)
  
  df_weight_query_surgery <- 
    query_trajectories(df_trajectory = df_weight_traj_surgery,
                       time_points = 0:study_duration,
                       error = weight_models$surgery$errors,
                       surgery = T,
                       a1c = F,
                       spline_knots = weight_models$surgery$spline$knots)
  
  df_a1c_query_surgery <- 
    query_trajectories(df_trajectory = df_a1c_traj_surgery,
                       time_points = 0:study_duration,
                       error = a1c_models$surgery$errors,
                       surgery = T,
                       a1c = T,
                       spline_knots = a1c_models$surgery$spline$knots,
                       baseline_vals = df_surg_subjects$baseline_hgba1c)
  
  df_query_surgery <- 
    df_weight_query_surgery %>% 
    select(subject_id, time, 'ptwc' = observed) %>% 
    inner_join(select(df_a1c_query_surgery, subject_id, time, 'pct_a1c' = observed), 
               by = c('subject_id', 'time'))
  
  df_surgery <- 
    df_subjects %>% 
    select(-baseline_bmi) %>% 
    inner_join(select(df_surg_subjects, subject_id, baseline_bmi), by = 'subject_id') %>% 
    inner_join(df_query_surgery, by = 'subject_id') %>% 
    inner_join(surgery_times, by = 'subject_id') %>% 
    mutate('time' = time + treatment_time) %>% 
    filter(time <= study_duration) %>% 
    mutate('bmi' = baseline_bmi * (1 + ptwc),
           'hgba1c' = baseline_hgba1c * (1 + pct_a1c),
           'age' = baseline_age + time/12) %>% 
    mutate('surgery' = 1) 
  
  
  ### (5) Put together surgery and control weight months
  keep_vars <- 
    c('subject_id', 'time', 'treatment_time', 'bmi', 'age', 'hgba1c')
  
  df_control <- 
    df_control %>% 
    inner_join(surgery_times, by = 'subject_id') %>% 
    filter(is.na(treatment_time) | time < treatment_time) %>% 
    mutate('bs_type' = 'control') %>%
    select(any_of(keep_vars), any_of(clean_names(params, names(.))))
  
  df_surgery <- 
    df_surgery %>% 
    select(any_of(keep_vars), any_of(clean_names(params, names(.))))
  
  df_final <- 
    df_control %>% 
    bind_rows(df_surgery) %>% 
    arrange(subject_id, time) 
  
  ### (6) Outcome
  df_final <- 
    df_final %>% 
    mutate('prob_outcome' =  expit(compute_model(df = df_final, beta = outcome_model))) %>% 
    mutate('outcome' = rbinom(n = nrow(df_final), size = 1, prob = prob_outcome)) %>% 
    group_by(subject_id) %>% 
    filter(cumsum(cumsum(outcome)) <= 1) %>%  ### Filter to first occurrence of outcome 
    ungroup() %>% 
    select(-prob_outcome)
  
  ### (7) Compute eligibility for each month 
  if(params$study_design == '(1) RYGB vs. VSG') {
    df_final <- 
      df_final %>% 
      mutate('eligible' = as.numeric(!is.na(treatment_time) & 
                                       time == treatment_time & 
                                       age >= 18 & 
                                       age <= 80 & 
                                       bmi >= 35 & 
                                       hgba1c >= 5.7) ### Diabetic criteria
      )
  } else if(params$study_design == '(2) Surgery vs. No Surgery') {
    df_final <- 
      df_final %>% 
      mutate('eligible' = as.numeric(age >= 18 & 
                                       age <= 80 & 
                                       bmi >= 35 & 
                                       (is.na(treatment_time) | time <= treatment_time))) ### haven't had surgery before
  } else if(params$study_design == '(3) Surgery vs. No Surgery (Diabetic)') {
    df_final <- 
      df_final %>% 
      mutate('eligible' = as.numeric(age >= 18 & 
                                       age <= 80 & 
                                       bmi >= 35 & 
                                       (is.na(treatment_time) | time <= treatment_time) & ### haven't had surgery before
                                       hgba1c >= 5.7)) ### Diabetic criteria 
  }
  
  return(df_final) 
}

impose_missingness <- function(df_final, params) {
  elig_criteria <- params$eligibility
  complete_case_elig_model <- params$cc_model
  
  ### (8) Impose missingness in the eligibility
  ### Currently doing monotone missingness in eligibility
  if(is.null(complete_case_elig_model)) {
    df_final$R <- 1
  } else if(elig_criteria == 'pre_diabetes' | elig_criteria == 'diabetes') {
    df_final <- 
      df_final %>% 
      mutate('prob_R' = expit(compute_model(df = df_final,
                                            beta = complete_case_elig_model))) %>% 
      mutate('R' = rbinom(n = nrow(df_final), size = 1, prob = prob_R)) %>% 
      mutate('hgba1c' = replace(hgba1c, R == 0, NA)) %>% 
      mutate('bmi' = replace(bmi, R == 0, NA)) %>% 
      mutate('eligible' = replace(eligible, R == 0, NA)) %>% 
      select(-prob_R)
  } else if(elig_criteria == 'bmi_only') {
    df_final <- 
      df_final %>% 
      mutate('prob_R' = expit(compute_model(df = df_final,
                                            beta = complete_case_elig_model))) %>% 
      mutate('R' = rbinom(n = nrow(df_final), size = 1, prob = prob_R)) %>%
      ### Impose missingness on BMI value but not on the eligiblity of surgical patients
      ### At some point we might want to call R = 1 here even though BMI could be missing 
      ### at the time of surgery
      mutate('bmi' = replace(bmi, R == 0, NA)) %>% 
      mutate('eligible' = replace(eligible, R == 0 & surgery != 1, NA)) %>% 
      select(-prob_R)
    
  } 
  return(df_final)
}

build_trial <- function(df, trial_ids, params) {
  
  ### Creat duplicates of dataset x # of trial
  dt <- data.table(map_dfr(trial_ids, ~mutate(df, 'trial_id' = .x)))
  
  ### Filter down to eligible person-trials
  df_elig <- 
    dt %>% 
    filter(time == trial_id - 1) %>%
    filter(eligible == 1) %>% 
    select(subject_id, trial_id)
  
  dt <- dt[df_elig, on = c('subject_id', 'trial_id')]
  dt <- dt[time >= trial_id - 1]
  
  setorder(dt, time)
  dt[, ':=' ('censor' = surgery != first(surgery),
             'follow_up' = 0:(.N -1),
             'surgery' = first(surgery), ### Treatment for this trial 
             'baseline_bmi' = first(bmi),
             'baseline_hgba1c' = first(hgba1c),
             'baseline_age' = first(age)),
     by = c('subject_id', 'trial_id')]
  
  ### Remove all censored observations except the first one which is needed for IPCW
  dt[, 'rm' := lag(censor, default = F), by = c('subject_id', 'trial_id')]
  dt <- dt[rm == F]
  
  df_trial <- 
    as_tibble(dt) %>% 
    mutate('trial_id' = trial_id) %>% 
    select(trial_id, subject_id, time, follow_up, outcome, censor, 
           any_of(clean_names(params, names(.))))
  
  
  return(df_trial)
}

build_trial_bootstrap <- function(df, trial_ids, params) {
  
  ### Creat duplicates of dataset x # of trial
  dt <- data.table(map_dfr(trial_ids, ~mutate(df, 'trial_id' = .x)))
  
  ### Filter down to eligible person-trials
  df_elig <- 
    dt %>% 
    filter(time == trial_id - 1) %>%
    filter(eligible == 1) %>% 
    select(subject_id, trial_id, replicate_id)
  
  dt <- dt[df_elig, on = c('subject_id', 'trial_id', 'replicate_id')]
  dt <- dt[time >= trial_id - 1]
  
  setorder(dt, time)
  dt[, ':=' ('censor' = surgery != first(surgery),
             'follow_up' = 0:(.N -1),
             'surgery' = first(surgery), ### Treatment for this trial 
             'baseline_bmi' = first(bmi),
             'baseline_hgba1c' = first(hgba1c),
             'baseline_age' = first(age)),
     by = c('subject_id', 'trial_id', 'replicate_id')]
  
  ### Remove all censored observations except the first one which is needed for IPCW
  dt[, 'rm' := lag(censor, default = F), by = c('subject_id', 'trial_id', 'replicate_id')]
  dt <- dt[rm == F]
  
  df_trial <- 
    as_tibble(dt) %>% 
    mutate('trial_id' = trial_id) %>% 
    select(trial_id, subject_id, replicate_id, time, follow_up, outcome, censor, 
           any_of(clean_names(params, names(.))))
  
  
  return(df_trial)
}

surgery_only_trial <- function(df, params) {
  ### Include eligible at time of surgery
  elig_ids <- 
    df %>% 
    filter(time == treatment_time) %>% 
    filter(eligible == 1) %>% 
    pull(subject_id)
  
  df_trial <- 
    df %>% 
    filter(subject_id %in% elig_ids) %>% 
    filter(time >= treatment_time) %>%  
    arrange(time) %>% 
    group_by(subject_id) %>% 
    mutate('follow_up' = 0:(n()-1)) %>% 
    mutate('baseline_bmi' = first(bmi),
           'baseline_hgba1c' = first(hgba1c),
           'baseline_age' = first(age)) %>% 
    ungroup() %>% 
    mutate('trial_id' = 1) %>% 
    mutate('bs_type' = as.numeric(bs_type == 'rygb')) %>% 
    select(trial_id, subject_id, time, follow_up, outcome, 
           any_of(clean_names(params, names(.))))
  
  return(df_trial)
  
}

### Function to build expanded dataset for analysis via pooled logistic regression
build_analysis_dataset <- function(df, params) {
  ### Don't need sequential trials for the first case
  if(params$study_design == '(1) RYGB vs. VSG') {
    df_trials <- surgery_only_trial(df, params)
  } else { ### need sequential trials
    df_trials <- build_trial(df, 1:params$n_trials, params)
  }
  
  return(df_trials)
}

build_analysis_dataset_bootstrap <- function(df, params) {
  ### need sequential trials
  df_trials <- build_trial_bootstrap(df, 1:params$n_trials, params)
  
  return(df_trials)
}

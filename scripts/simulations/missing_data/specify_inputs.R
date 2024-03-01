library(tidyverse)
source('scripts/simulations/helpers.R')

### Weight Trajectory Model Specification
### Load in Saved out Surgery spline knots
spline_knots_weight <- read_rds('models/eda/surgery_spline_weight.rds')
spline_knots_a1c <- read_rds('models/eda/surgery_spline_a1c.rds')

### Augment this to be in the 5 groups
spline_knots_weight$coeff <- 
  list('1' =  spline_knots_weight$coeff + c(seq(0.01, 0.08, length.out = 5), 0.1, seq(0.1, 0.05, length.out =  5)),
       '2' = spline_knots_weight$coeff + c(seq(0.005, 0.03, length.out = 5), 0.05, seq(0.05, 0.025, length.out =  5)),
       '3' = spline_knots_weight$coeff,
       '4' = spline_knots_weight$coeff -  seq(0.005, 0.08, length.out = 11),
       '5' = spline_knots_weight$coeff - seq(0.01, 0.15, length.out = 11))

spline_knots_a1c$coeff <- 
  list('1' =  spline_knots_a1c$coeff + c(seq(0.01, 0.08, length.out = 6)),
       '2' = spline_knots_a1c$coeff + c(seq(0.005, 0.045, length.out = 6)),
       '3' = spline_knots_a1c$coeff ,
       '4' = spline_knots_a1c$coeff -  seq(0.005, 0.04, length.out = 6),
       '5' = spline_knots_a1c$coeff - seq(0.01, 0.08, length.out = 6))

weight_models <- 
  list('control' = list('slope' = list('(Intercept)' = -5e-4),
                        'quadratic' = list('prob_nonzero' = list('(Intercept)' = logit(0.25)),
                                           'term' = list('(Intercept)' = 0)),
                        'std_dev' = list('slope' = 5e-4,
                                         'quadratic' = 2e-5),
                        'errors' = list('structure' = 'iid',
                                        'std_dev' = sqrt(1e-5))),
       
       'surgery' = list('group_model' = list('intercepts' = logit(cumsum(c(0.05, 0.2, 0.5, 0.2)))),
                        'std_dev' = 0.025,
                        'spline' = spline_knots_weight,
                        'errors' = list('structure' = 'iid',
                                        'std_dev' = sqrt(1e-5))))

a1c_models <- 
  list('control' = list('slope' = list('(Intercept)' = 0),
                        'std_dev' = 3e-4,
                        'errors' = list('structure' = 'heteroskedastic_by_a1c',
                                        'error_function' = function(x) { 1e-4 * x^2 } )),
       'surgery' = list('group_model' = list('intercepts' = logit(cumsum(c(0.05, 0.2, 0.5, 0.2)))),
                        'std_dev' = 0.025,
                        'spline' = spline_knots_a1c,
                        'errors' = list('structure' = 'heteroskedastic_by_a1c',
                                        'error_function' = function(x) { 1e-4 * x^2  })))


### Simulation Settings
df_settings <- 
  crossing('study_design' = c('(1) RYGB vs. VSG', '(2) Surgery vs. No Surgery', '(3) Surgery vs. No Surgery (Diabetic)'),
           'missingness' = c('(1) M-Bias', '(2) Treatment Heterogeneity', '(3) M-Bias w/ Mediator', '(4) MNAR')) %>% 
  mutate('sim_id' = 1:nrow(.)) %>% 
  select(sim_id, everything())
write_csv(df_settings, 'data/simulations/missing_data/inputs/settings.csv')


### Simulation Parameters: 
###
### sim_id: id for simulation setting 
### n_subjects: # of subjects in population
### n_trials: # of sequential trials 
### n_sims: # of simulations
### study_duration: length of study from start of first emulated trial 
### weight_models: list of coefficients/parameters for weight trajectory models
### a1c_models: list of coefficients/parameters for a1c trajectory models
### treatment_model: model for monthly probability of getting surgery
### outcome_model: model for monthly probability of outcome
### cc_mode: monthly missing data models for eligibility 
### eligibility: eligibility criteria
### truth_model: The true model to fit
### analysis_models: list of models for analysis of data

###############################
### RYGB vs. Sleeve Surgery ###
###############################

### Simulation 1 (M-Bias)
params <- 
  list('sim_id' = 1,
       'study_design' = df_settings$study_design[1],
       'missingness' = df_settings$missingness[1],
       'n_subjects' = 50000,
       'n_trials' = 1,
       'n_sims' = 1000,
       'study_duration' = 60,
       'weight_models' = weight_models,
       'a1c_models' = a1c_models,
       'treatment_model' = list('(Intercept)' = logit(0.03)),
       'bs_type_model' = list('(Intercept)' = logit(0.6),
                              'site' = log(2)),
       'outcome_model' = list('(Intercept)' = logit(0.015),
                              'elix_score' = 1.2,
                              'bs_type[rygb]' = log(0.7)),
       'cc_model' =  list('(Intercept)' = logit(0.08),
                          'elix_score' = 0.8,
                          'site' = 1.1),
       'eligibility' = 'pre_diabetes',
       'truth_model' = list('outcome_model' = outcome ~ bs_type),
       'analysis_models' = list(
         ### Main Effects Model
         'Y ~ A' = list('outcome_model' = outcome ~ bs_type),
         
         ### Main Effects Model w/ IPWR
         'Y ~ A (IPWR: R ~ L^{R,A})' = list('outcome_model' = outcome ~ bs_type,
                                            'cc_model' = R ~ site,
                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R,Y})' = list('outcome_model' = outcome ~ bs_type,
                                            'cc_model' = R ~ elix_score,
                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R})' = list('outcome_model' = outcome ~ bs_type,
                                          'cc_model' = R ~ site + elix_score,
                                          'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R} + A)' = list('outcome_model' = outcome ~ bs_type,
                                              'cc_model' = R ~ site + elix_score + bs_type,
                                              'ipwr_purpose' = 'Estimation'))
  )

write_rds(params, 'data/simulations/missing_data/inputs/sim_params_1.rds')


### Simulation 2 (Treatment Effect Heterogeneity)
params <- 
  list('sim_id' = 2,
       'study_design' = df_settings$study_design[2],
       'missingness' = df_settings$missingness[2],
       'n_subjects' = 50000,
       'n_trials' = 1,
       'n_sims' = 1000,
       'study_duration' = 60,
       'weight_models' = weight_models,
       'a1c_models' = a1c_models,
       'treatment_model' = list('(Intercept)' = logit(0.03)),
       'bs_type_model' = list('(Intercept)' = logit(0.6),
                              'site' = log(2)),
       'outcome_model' = list('(Intercept)' = logit(0.02),
                              'elix_score' = 0.25,
                              'bs_type[rygb]:elix_score' = 0.1,
                              'bs_type[rygb]' = log(0.7)),
       'cc_model' =  list('(Intercept)' = logit(0.1),
                          'elix_score' = 0.5,
                          'site' = 0.25),
       'eligibility' = 'pre_diabetes',
       'truth_model' = list('outcome_model' = outcome ~ bs_type),
       'analysis_models' = list(
         ### Main Effects Model 
         'Y ~ A' = list('outcome_model' = outcome ~ bs_type),
         
         ### Main Effects Models w/ IPWR
         'Y ~ A  (IPWR: R ~ L^{R,A})' = list('outcome_model' = outcome ~ bs_type,
                                             'cc_model' = R ~ site,
                                             'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R,Y})' = list('outcome_model' = outcome ~ bs_type,
                                            'cc_model' = R ~ elix_score,
                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R})' = list('outcome_model' = outcome ~ bs_type,
                                          'cc_model' = R ~ site + elix_score,
                                          'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R} + A)' = list('outcome_model' = outcome ~ bs_type,
                                              'cc_model' = R ~ site + elix_score + bs_type,
                                              'ipwr_purpose' = 'Estimation'))
  )


write_rds(params, 'data/simulations/missing_data/inputs/sim_params_2.rds')

### Simulation 3 (M-Bias w/ Mediator)
weight_models_3 <- 
  list('control' = list('slope' = list('(Intercept)' = -5e-4,
                                       'gender' = -1e-4,
                                       'race' = -2e-4,
                                       'insulin' = 5e-4),
                        'quadratic' = list('prob_nonzero' = list('(Intercept)' = logit(0.25),
                                                                 'insulin' = 0.5),
                                           'term' = list('(Intercept)' = 0,
                                                         'insulin' = 1e-5)),
                        'std_dev' = list('slope' = 5e-4,
                                         'quadratic' = 2e-5),
                        'errors' = list('structure' = 'iid',
                                        'std_dev' = sqrt(1e-5))),
       
       'surgery' = list('group_model' = list('intercepts' = logit(cumsum(c(0.05, 0.2, 0.5, 0.2))),
                                             'elix_score' = -5,
                                             'bs_type[rygb]' = 4,
                                             'insulin' = -3),
                        'std_dev' = 0.025,
                        'spline' = spline_knots_weight,
                        'errors' = list('structure' = 'iid',
                                        'std_dev' = sqrt(1e-5))))

a1c_models_3 <- 
  list('control' = list('slope' = list('(Intercept)' = 0,
                                       'insulin' = -1e-3),
                        'std_dev' = 3e-4,
                        'errors' = list('structure' = 'heteroskedastic_by_a1c',
                                        'error_function' = function(x) { 1e-4 * x^2 } )),
       'surgery' = list('group_model' =  list('intercepts' = logit(cumsum(c(0.2, 0.25, 0.4, 0.1))),
                                              'elix_score' = -5,
                                              'bs_type[rygb]' = 4,
                                              'insulin' = -3),
                        'std_dev' = 0.025,
                        'spline' = spline_knots_a1c,
                        'errors' = list('structure' = 'heteroskedastic_by_a1c',
                                        'error_function' = function(x) { 1e-4 * x^2  })))

params <- 
  list('sim_id' = 3,
       'study_design' = df_settings$study_design[3],
       'missingness' = df_settings$missingness[3],
       'n_subjects' = 50000,
       'n_trials' = 1,
       'n_sims' = 1000,
       'study_duration' = 60,
       'weight_models' = weight_models_3,
       'a1c_models' = a1c_models_3,
       'treatment_model' = list('(Intercept)' = logit(0.03)),
       'bs_type_model' = list('(Intercept)' = logit(0.6),
                              'site' = log(2)),
       'outcome_model' = list('(Intercept)' = logit(0.0001),
                              'bmi' = 0.105,
                              'hgba1c' = 0.185),
       'cc_model' =  list('(Intercept)' = logit(0.1),
                          'elix_score' = 0.6,
                          'site' = -0.55,
                          'insulin' = 0.5),
       'eligibility' = 'pre_diabetes',  
       'truth_model' = list('outcome_model' = outcome ~ bs_type),
       'analysis_models' = list(
         ### Main Effects Model
         'Y ~ A' = list('outcome_model' = outcome ~ bs_type),
         
         ### Main Effects Model w/ IPWR
         'Y ~ A (IPWR: R ~ L^{R,A})' = list('outcome_model' = outcome ~ bs_type,
                                            'cc_model' = R ~ site,
                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R,Y})' = list('outcome_model' = outcome ~ bs_type,
                                            'cc_model' = R ~ elix_score + insulin,
                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R})' = list('outcome_model' = outcome ~ bs_type,
                                          'cc_model' = R ~ site + elix_score + insulin,
                                          'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R} + A)' = list('outcome_model' = outcome ~ bs_type,
                                              'cc_model' = R ~ site + elix_score + insulin + bs_type,
                                              'ipwr_purpose' = 'Estimation'))
  )

write_rds(params, 'data/simulations/missing_data/inputs/sim_params_3.rds')

### Simulation 4 (MNAR)
params <- 
  list('sim_id' = 4,
       'study_design' = df_settings$study_design[4],
       'missingness' = df_settings$missingness[4],
       'n_subjects' = 50000,
       'n_trials' = 1,
       'n_sims' = 1000,
       'study_duration' = 60,
       'weight_models' = weight_models,
       'a1c_models' = a1c_models,
       'treatment_model' = list('(Intercept)' = logit(0.03)),
       'bs_type_model' = list('(Intercept)' = logit(0.00075),
                              'bmi' = 0.2),
       'outcome_model' = list('(Intercept)' = logit(5e-4),
                              'baseline_hgba1c' = 0.75,
                              'bs_type[rygb]' = log(0.7)),
       'cc_model' =  list('(Intercept)' = logit(1e-5),
                          'bmi' = 0.08,
                          'hgba1c' = 1.2),
       'eligibility' = 'pre_diabetes',
       'truth_model' = list('outcome_model' = outcome ~ bs_type),
       'analysis_models' = list(
         ## Main Effects Model
         ### Can't even estimate IPWR models
         'Y ~ A' = list('outcome_model' = outcome ~ bs_type)
       )
  )

write_rds(params, 'data/simulations/missing_data/inputs/sim_params_4.rds')


#################################
### Surgery vs. No Surgery ######
#################################

### Simulation 5 (M-Bias)
params <- 
  list('sim_id' = 5,
       'study_design' = df_settings$study_design[5],
       'missingness' = df_settings$missingness[5],
       'n_subjects' = 10000,
       'n_trials' = 12,
       'n_sims' = 1000,
       'study_duration' = 36,
       'weight_models' = weight_models,
       'a1c_models' = a1c_models,
       'treatment_model' = list('(Intercept)' = logit(0.04),
                                'smoking_status[current]' = -2,
                                'smoking_status[former]' = -0.75),
       'bs_type_model' = list('(Intercept)' = logit(0.6),
                              'site' = log(2)),
       'outcome_model' = list('(Intercept)' = logit(0.015),
                              'elix_score' = 1.2,
                              'surgery' = log(0.7)),
       'cc_model' =  list('(Intercept)' = logit(0.08),
                          'elix_score' = 0.8,
                          'site' = 0.5,
                          'smoking_status[former]' = -0.5,
                          'smoking_status[current]' = -1),
       'eligibility' = 'bmi_only',
       'truth_model' = list('outcome_model' = outcome ~ surgery,
                            'censor_model' = censor ~ smoking_status),
       'analysis_models' = list(
         ### Main Effects Model
         'Y ~ A' = list('outcome_model' = outcome ~ surgery),
         
         ### Main Effects Model w/ IPWR
         'Y ~ A (IPWR: R ~ L^{R,A})' = list('outcome_model' = outcome ~ surgery,
                                            'cc_model' = R ~ smoking_status + site,
                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R,Y})' = list('outcome_model' = outcome ~ surgery,
                                            'cc_model' = R ~ elix_score,
                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R})' = list('outcome_model' = outcome ~ surgery,
                                          'cc_model' = R ~ smoking_status + site + elix_score,
                                          'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R} + A)' = list('outcome_model' = outcome ~ surgery,
                                              'cc_model' = R ~ smoking_status + site + elix_score + surgery,
                                              'ipwr_purpose' = 'Estimation'),
         
         ### Main Effects Model w/ IPWR (Single Side)
         'Y ~ A (IPWR: R ~ L^{R,A} Stratified by A)' = list('outcome_model' = outcome ~ surgery,
                                                            'cc_model' = R ~ smoking_status + site,
                                                            'ipwr_non_surgery' = T,
                                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R,Y} Stratified by A)' = list('outcome_model' = outcome ~ surgery,
                                                            'cc_model' = R ~ elix_score,
                                                            'ipwr_non_surgery' = T,
                                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R} Stratified by A)' = list('outcome_model' = outcome ~ surgery,
                                                          'cc_model' = R ~ smoking_status + site + elix_score,
                                                          'ipwr_non_surgery' = T,
                                                          'ipwr_purpose' = 'Estimation'),
         
         ### Main Effects Model (w/ IPCW)
         'Y ~ A (IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                          'censor_model' = censor ~ smoking_status),
         
         ### Main Effects Model w/ IPWR x IPCW
         'Y ~ A (IPWR: R ~ L^{R,A}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                             'cc_model' = R ~ smoking_status + site,
                                                             'ipwr_purpose' = 'Estimation',
                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R,Y}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                             'cc_model' = R ~ elix_score,
                                                             'ipwr_purpose' = 'Estimation',
                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                           'cc_model' = R ~ smoking_status + site + elix_score,
                                                           'ipwr_purpose' = 'Estimation',
                                                           'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R} + A, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                               'cc_model' = R ~ smoking_status + site + elix_score + surgery,
                                                               'ipwr_purpose' = 'Estimation',
                                                               'censor_model' = censor ~ smoking_status),
         
         ### Main Effects Model w/ IPWR x IPCW (Single Side)
         'Y ~ A (IPWR: R ~ L^{R,A} Stratified by A, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                                             'cc_model' = R ~ smoking_status + site,
                                                                             'ipwr_purpose' = 'Estimation',
                                                                             'ipwr_non_surgery' = T,
                                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R,Y} Stratified by A, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                                             'cc_model' = R ~ elix_score,
                                                                             'ipwr_purpose' = 'Estimation',
                                                                             'ipwr_non_surgery' = T,
                                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R} Stratified by A, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                                           'cc_model' = R ~ smoking_status + site + elix_score,
                                                                           'ipwr_purpose' = 'Estimation',
                                                                           'ipwr_non_surgery' = T,
                                                                           'censor_model' = censor ~ smoking_status))
  )

write_rds(params, 'data/simulations/missing_data/inputs/sim_params_5.rds')


### Simulation 6 (Treatment Effect Heterogeneity)
params <- 
  list('sim_id' = 6,
       'study_design' = df_settings$study_design[6],
       'missingness' = df_settings$missingness[6],
       'n_subjects' = 10000,
       'n_trials' = 12,
       'n_sims' = 1000,
       'study_duration' = 36,
       'weight_models' = weight_models,
       'a1c_models' = a1c_models,
       'treatment_model' = list('(Intercept)' = logit(0.04),
                                'smoking_status[current]' = -2,
                                'smoking_status[former]' = -0.75),
       'bs_type_model' = list('(Intercept)' = logit(0.6),
                              'site' = log(2)),
       'outcome_model' = list('(Intercept)' = logit(0.025),
                              'elix_score' = 0.25,
                              'surgery:elix_score' = 0.1,
                              'surgery' = log(0.7)),
       'cc_model' =  list('(Intercept)' = logit(0.08),
                          'elix_score' = 0.5,
                          'site' = 0.25,
                          'smoking_status[former]' = -0.5,
                          'smoking_status[current]' = -1),
       'eligibility' = 'bmi_only',
       'truth_model' = list('outcome_model' = outcome ~ surgery,
                            'interaction' = F,
                            'censor_model' = censor ~ smoking_status),
       'analysis_models' = list(
         ### Main Effects Model
         'Y ~ A' = list('outcome_model' = outcome ~ surgery),
         
         ### Main Effects Model (w/ IPCW)
         'Y ~ A (IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                          'censor_model' = censor ~ smoking_status),
         ### Main Effects Model w/ IPWR
         'Y ~ A (IPWR: R ~ L^{R,A})' = list('outcome_model' = outcome ~ surgery,
                                            'cc_model' = R ~ smoking_status + site,
                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R,Y})' = list('outcome_model' = outcome ~ surgery,
                                            'cc_model' = R ~ elix_score,
                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R})' = list('outcome_model' = outcome ~ surgery,
                                          'cc_model' = R ~ smoking_status + site + elix_score,
                                          'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R} + A)' = list('outcome_model' = outcome ~ surgery,
                                              'cc_model' = R ~ smoking_status + site + elix_score + surgery,
                                              'ipwr_purpose' = 'Estimation'),
         
         ### Main Effects Model w/ IPWR (Single Side)
         'Y ~ A (IPWR: R ~ L^{R,A} Stratified by A)' = list('outcome_model' = outcome ~ surgery,
                                                            'cc_model' = R ~ smoking_status + site,
                                                            'ipwr_non_surgery' = T,
                                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R,Y} Stratified by A)' = list('outcome_model' = outcome ~ surgery,
                                                            'cc_model' = R ~ elix_score,
                                                            'ipwr_non_surgery' = T,
                                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R} Stratified by A)' = list('outcome_model' = outcome ~ surgery,
                                                          'cc_model' = R ~ smoking_status + site + elix_score,
                                                          'ipwr_non_surgery' = T,
                                                          'ipwr_purpose' = 'Estimation'),
         
         ### Main Effects Model w/ IPWR x IPCW
         'Y ~ A (IPWR: R ~ L^{R,A}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                             'cc_model' = R ~  smoking_status + site,
                                                             'ipwr_purpose' = 'Estimation',
                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R,Y}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                             'cc_model' = R ~ elix_score,
                                                             'ipwr_purpose' = 'Estimation',
                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                           'cc_model' = R ~ smoking_status + site + elix_score,
                                                           'ipwr_purpose' = 'Estimation',
                                                           'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R} + A, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                               'cc_model' = R ~ smoking_status + site + elix_score + surgery,
                                                               'ipwr_purpose' = 'Estimation',
                                                               'censor_model' = censor ~ smoking_status),
         
         ### Main Effects Model w/ IPWR x IPCW (Single Side)
         'Y ~ A (IPWR: R ~ L^{R,A} Stratified by A, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                                             'cc_model' = R ~ smoking_status + site,
                                                                             'ipwr_purpose' = 'Estimation',
                                                                             'ipwr_non_surgery' = T,
                                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R,Y} Stratified by A, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                                             'cc_model' = R ~ elix_score,
                                                                             'ipwr_purpose' = 'Estimation',
                                                                             'ipwr_non_surgery' = T,
                                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R} Stratified by A, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                                           'cc_model' = R ~ smoking_status + site + elix_score,
                                                                           'ipwr_purpose' = 'Estimation',
                                                                           'ipwr_non_surgery' = T,
                                                                           'censor_model' = censor ~ smoking_status))
       
  )

write_rds(params, 'data/simulations/missing_data/inputs/sim_params_6.rds')

### Simulation 7 (M-Bias Through Mediator)
params <- 
  list('sim_id' = 7,
       'study_design' = df_settings$study_design[7],
       'missingness' = df_settings$missingness[7],
       'n_subjects' = 10000,
       'n_trials' = 12,
       'n_sims' = 1000,
       'study_duration' = 36,
       'weight_models' = weight_models_3,
       'a1c_models' = a1c_models_3,
       'treatment_model' = list('(Intercept)' = logit(0.04),
                                'smoking_status[current]' = -2,
                                'smoking_status[former]' = -0.75),
       'bs_type_model' = list('(Intercept)' = logit(0.6),
                              'site' = log(2)),
       'outcome_model' = list('(Intercept)' = logit(0.0001),
                              'bmi' = 0.105,
                              'hgba1c' = 0.185),
       'cc_model' =   list('(Intercept)' = logit(0.25),
                           'elix_score' = 0.6,
                           'site' = -0.55,
                           'insulin' = 0.5),
       'eligibility' = 'bmi_only',  
       'truth_model' = list('outcome_model' = outcome ~ surgery,
                            'censor_model' = censor ~ smoking_status),
       'analysis_models' = list(
         ### Main Effects Model
         'Y ~ A' = list('outcome_model' = outcome ~ surgery),
         
         ### Main Effects Model w/ IPWR
         'Y ~ A (IPWR: R ~ L^{R,A})' = list('outcome_model' = outcome ~ surgery,
                                            'cc_model' = R ~ site,
                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R,Y})' = list('outcome_model' = outcome ~ surgery,
                                            'cc_model' = R ~ elix_score + insulin,
                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R})' = list('outcome_model' = outcome ~ surgery,
                                          'cc_model' = R ~ site + elix_score + insulin,
                                          'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R} + A)' = list('outcome_model' = outcome ~ surgery,
                                              'cc_model' = R ~ site + elix_score + insulin + surgery,
                                              'ipwr_purpose' = 'Estimation'),
         
         ### Main Effects Model w/ IPWR (Single Side)
         'Y ~ A (IPWR: R ~ L^{R,A} Stratified by A)' = list('outcome_model' = outcome ~ surgery,
                                                            'cc_model' = R ~ site,
                                                            'ipwr_non_surgery' = T,
                                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R,Y} Stratified by A)' = list('outcome_model' = outcome ~ surgery,
                                                            'cc_model' = R ~ elix_score + insulin,
                                                            'ipwr_non_surgery' = T,
                                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R} Stratified by A)' = list('outcome_model' = outcome ~ surgery,
                                                          'cc_model' = R ~ site + elix_score + insulin,
                                                          'ipwr_non_surgery' = T,
                                                          'ipwr_purpose' = 'Estimation'),
         
         ### Main Effects Model (w/ IPCW)
         'Y ~ A (IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                          'censor_model' = censor ~ smoking_status),
         
         ### Main Effects Model w/ IPWR x IPCW
         'Y ~ A (IPWR: R ~ L^{R,A}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                             'cc_model' = R ~ site,
                                                             'ipwr_purpose' = 'Estimation',
                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R,Y}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                             'cc_model' = R ~ elix_score + insulin,
                                                             'ipwr_purpose' = 'Estimation',
                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                           'cc_model' = R ~  site + elix_score + insulin,
                                                           'ipwr_purpose' = 'Estimation',
                                                           'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R} + A, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                               'cc_model' = R ~ site + elix_score + insulin + surgery,
                                                               'ipwr_purpose' = 'Estimation',
                                                               'censor_model' = censor ~ smoking_status),
         
         ### Main Effects Model w/ IPWR x IPCW (Single Side)
         'Y ~ A (IPWR: R ~ L^{R,A} Stratified by A, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                                             'cc_model' = R ~ site,
                                                                             'ipwr_purpose' = 'Estimation',
                                                                             'ipwr_non_surgery' = T,
                                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R,Y} Stratified by A, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                                             'cc_model' = R ~ elix_score + insulin,
                                                                             'ipwr_purpose' = 'Estimation',
                                                                             'ipwr_non_surgery' = T,
                                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R} Stratified by A, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                                           'cc_model' = R ~ site + elix_score + insulin,
                                                                           'ipwr_purpose' = 'Estimation',
                                                                           'ipwr_non_surgery' = T,
                                                                           'censor_model' = censor ~ smoking_status))
  )

write_rds(params, 'data/simulations/missing_data/inputs/sim_params_7.rds')


### Simulation 8 (MNAR)
params <- 
  list('sim_id' = 8,
       'study_design' = df_settings$study_design[8],
       'missingness' = df_settings$missingness[8],
       'n_subjects' = 10000,
       'n_trials' = 12,
       'n_sims' = 1000,
       'study_duration' = 36,
       'weight_models' = weight_models,
       'a1c_models' = a1c_models,
       'treatment_model' = list('(Intercept)' = logit(0.04),
                                'smoking_status[current]' = -2,
                                'smoking_status[former]' = -0.75),
       'bs_type_model' = list('(Intercept)' = logit(0.00075),
                              'bmi' = 0.2),
       'outcome_model' = list('(Intercept)' = logit(5e-4),
                              'baseline_hgba1c' = 0.75,
                              'surgery' = log(0.7)),
       'cc_model' =  list('(Intercept)' = logit(1e-5),
                          'bmi' = 0.08,
                          'hgba1c' = 1.2),
       'eligibility' = 'bmi_only',
       'truth_model' = list('outcome_model' = outcome ~ surgery,
                            'censor_model' = censor ~ smoking_status),
       'analysis_models' = list(
         ## Main Effects Model
         ### Can't even estimate IPWR models
         'Y ~ A' = list('outcome_model' = outcome ~ surgery),
         
         ### Main Effects Model (w/ IPCW)
         'Y ~ A (IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                          'censor_model' = censor ~ smoking_status))
  )

write_rds(params, 'data/simulations/missing_data/inputs/sim_params_8.rds')



#################################################
### Surgery vs. No Surgery (w/ A1C Cutoff) ######
#################################################

### Simulation 9 (M-Bias)
params <- 
  list('sim_id' = 9,
       'study_design' = df_settings$study_design[9],
       'missingness' = df_settings$missingness[9],
       'n_subjects' = 10000,
       'n_trials' = 12,
       'n_sims' = 1000,
       'study_duration' = 36,
       'weight_models' = weight_models,
       'a1c_models' = a1c_models,
       'treatment_model' = list('(Intercept)' = logit(0.04),
                                'smoking_status[current]' = -2,
                                'smoking_status[former]' = -0.75),
       'bs_type_model' = list('(Intercept)' = logit(0.6),
                              'site' = log(2)),
       'outcome_model' = list('(Intercept)' = logit(0.015),
                              'elix_score' = 1.2,
                              'surgery' = log(0.7)),
       'cc_model' =  list('(Intercept)' = logit(0.08),
                          'elix_score' = 0.8,
                          'site' = 0.5,
                          'smoking_status[former]' = -0.5,
                          'smoking_status[current]' = -1),
       'eligibility' = 'pre_diabetes',
       'truth_model' = list('outcome_model' = outcome ~ surgery,
                            'censor_model' = censor ~ smoking_status),
       'analysis_models' = list(
         ### Main Effects Model
         'Y ~ A' = list('outcome_model' = outcome ~ surgery),
         
         ### Main Effects Model w/ IPWR
         'Y ~ A (IPWR: R ~ L^{R,A})' = list('outcome_model' = outcome ~ surgery,
                                            'cc_model' = R ~ smoking_status + site,
                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R,Y})' = list('outcome_model' = outcome ~ surgery,
                                            'cc_model' = R ~ elix_score,
                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R})' = list('outcome_model' = outcome ~ surgery,
                                          'cc_model' = R ~ smoking_status + site + elix_score,
                                          'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R} + A)' = list('outcome_model' = outcome ~ surgery,
                                              'cc_model' = R ~ smoking_status + site + elix_score + surgery,
                                              'ipwr_purpose' = 'Estimation'),
         
         ### Main Effects Model (w/ IPCW)
         'Y ~ A (IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                          'censor_model' = censor ~ smoking_status),
         
         ### Main Effects Model w/ IPWR x IPCW
         'Y ~ A (IPWR: R ~ L^{R,A}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                             'cc_model' = R ~ smoking_status + site,
                                                             'ipwr_purpose' = 'Estimation',
                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R,Y}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                             'cc_model' = R ~ elix_score,
                                                             'ipwr_purpose' = 'Estimation',
                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                           'cc_model' = R ~ smoking_status + site + elix_score,
                                                           'ipwr_purpose' = 'Estimation',
                                                           'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R} + A, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                               'cc_model' = R ~ smoking_status + site + elix_score + surgery,
                                                               'ipwr_purpose' = 'Estimation',
                                                               'censor_model' = censor ~ smoking_status))
  )

write_rds(params, 'data/simulations/missing_data/inputs/sim_params_9.rds')


### Simulation 10 (Treatment Effect Heterogeneity)
params <- 
  list('sim_id' = 10,
       'study_design' = df_settings$study_design[10],
       'missingness' = df_settings$missingness[10],
       'n_subjects' = 10000,
       'n_trials' = 12,
       'n_sims' = 1000,
       'study_duration' = 36,
       'weight_models' = weight_models,
       'a1c_models' = a1c_models,
       'treatment_model' = list('(Intercept)' = logit(0.04),
                                'smoking_status[current]' = -2,
                                'smoking_status[former]' = -0.75),
       'bs_type_model' = list('(Intercept)' = logit(0.6),
                              'site' = log(2)),
       'outcome_model' = list('(Intercept)' = logit(0.025),
                              'elix_score' = 0.25,
                              'surgery:elix_score' = 0.1,
                              'surgery' = log(0.7)),
       'cc_model' =  list('(Intercept)' = logit(0.08),
                          'elix_score' = 0.5,
                          'site' = 0.25,
                          'smoking_status[former]' = -0.5,
                          'smoking_status[current]' = -1),
       'eligibility' = 'pre_diabetes',
       'truth_model' = list('outcome_model' = outcome ~ surgery,
                            'censor_model' = censor ~ smoking_status),
       'analysis_models' = list(
         ### Main Effects Model
         'Y ~ A' = list('outcome_model' = outcome ~ surgery),
         
         ### Main Effects Model (w/ IPCW)
         'Y ~ A (IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                          'censor_model' = censor ~ smoking_status),
         
         
         ### Main Effects Model w/ IPWR
         'Y ~ A (IPWR: R ~ L^{R,A})' = list('outcome_model' = outcome ~ surgery,
                                            'cc_model' = R ~ smoking_status + site,
                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R,Y})' = list('outcome_model' = outcome ~ surgery,
                                            'cc_model' = R ~ elix_score,
                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R})' = list('outcome_model' = outcome ~ surgery,
                                          'cc_model' = R ~ smoking_status + site + elix_score,
                                          'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R} + A)' = list('outcome_model' = outcome ~ surgery,
                                              'cc_model' = R ~ smoking_status + site + elix_score + surgery,
                                              'ipwr_purpose' = 'Estimation'),
         
         ### Main Effects Model w/ IPWR x IPCW
         'Y ~ A (IPWR: R ~ L^{R,A}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                             'cc_model' = R ~  smoking_status + site,
                                                             'ipwr_purpose' = 'Estimation',
                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R,Y}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                             'cc_model' = R ~ elix_score,
                                                             'ipwr_purpose' = 'Estimation',
                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                           'cc_model' = R ~ smoking_status + site + elix_score,
                                                           'ipwr_purpose' = 'Estimation',
                                                           'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R} + A, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                               'cc_model' = R ~ smoking_status + site + elix_score + surgery,
                                                               'ipwr_purpose' = 'Estimation',
                                                               'censor_model' = censor ~ smoking_status))
  )

write_rds(params, 'data/simulations/missing_data/inputs/sim_params_10.rds')

### Simulation 11 (M-Bias Through Mediator)
params <- 
  list('sim_id' = 11,
       'study_design' = df_settings$study_design[11],
       'missingness' = df_settings$missingness[11],
       'n_subjects' = 10000,
       'n_trials' = 12,
       'n_sims' = 1000,
       'study_duration' = 36,
       'weight_models' = weight_models_3,
       'a1c_models' = a1c_models_3,
       'treatment_model' = list('(Intercept)' = logit(0.04),
                                'smoking_status[current]' = -2,
                                'smoking_status[former]' = -0.75),
       'bs_type_model' = list('(Intercept)' = logit(0.6),
                              'site' = log(2)),
       'outcome_model' = list('(Intercept)' = logit(0.0001),
                              'bmi' = 0.105,
                              'hgba1c' = 0.185),
       'cc_model' =   list('(Intercept)' = logit(0.25),
                           'elix_score' = 0.6,
                           'site' = -0.55,
                           'insulin' = 0.5),
       'eligibility' = 'pre_diabetes',  
       'truth_model' = list('outcome_model' = outcome ~ surgery,
                            'censor_model' = censor ~ smoking_status),
       'analysis_models' = list(
         ### Main Effects Model
         'Y ~ A' = list('outcome_model' = outcome ~ surgery),
         
         ### Main Effects Model w/ IPWR
         'Y ~ A (IPWR: R ~ L^{R,A})' = list('outcome_model' = outcome ~ surgery,
                                            'cc_model' = R ~ site,
                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R,Y})' = list('outcome_model' = outcome ~ surgery,
                                            'cc_model' = R ~ elix_score + insulin,
                                            'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R})' = list('outcome_model' = outcome ~ surgery,
                                          'cc_model' = R ~ site + elix_score + insulin,
                                          'ipwr_purpose' = 'Estimation'),
         'Y ~ A (IPWR: R ~ L^{R} + A)' = list('outcome_model' = outcome ~ surgery,
                                              'cc_model' = R ~ site + elix_score + insulin + surgery,
                                              'ipwr_purpose' = 'Estimation'),
         
         ### Main Effects Model (w/ IPCW)
         'Y ~ A (IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                          'censor_model' = censor ~ smoking_status),
         
         ### Main Effects Model w/ IPWR x IPCW
         'Y ~ A (IPWR: R ~ L^{R,A}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                             'cc_model' = R ~ site,
                                                             'ipwr_purpose' = 'Estimation',
                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R,Y}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                             'cc_model' = R ~ elix_score + insulin,
                                                             'ipwr_purpose' = 'Estimation',
                                                             'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                           'cc_model' = R ~  site + elix_score + insulin,
                                                           'ipwr_purpose' = 'Estimation',
                                                           'censor_model' = censor ~ smoking_status),
         'Y ~ A (IPWR: R ~ L^{R} + A, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                               'cc_model' = R ~ site + elix_score + insulin + surgery,
                                                               'ipwr_purpose' = 'Estimation',
                                                               'censor_model' = censor ~ smoking_status))
  )

write_rds(params, 'data/simulations/missing_data/inputs/sim_params_11.rds')

### Simulation 12 (MNAR)
params <- 
  list('sim_id' = 12,
       'study_design' = df_settings$study_design[12],
       'missingness' = df_settings$missingness[12],
       'n_subjects' = 10000,
       'n_trials' = 12,
       'n_sims' = 1000,
       'study_duration' = 36,
       'weight_models' = weight_models,
       'a1c_models' = a1c_models,
       'treatment_model' = list('(Intercept)' = logit(0.04),
                                'smoking_status[current]' = -2,
                                'smoking_status[former]' = -0.75),
       'bs_type_model' = list('(Intercept)' = logit(0.00075),
                              'bmi' = 0.2),
       'outcome_model' = list('(Intercept)' = logit(5e-4),
                              'baseline_hgba1c' = 0.75,
                              'surgery' = log(0.7)),
       'cc_model' =  list('(Intercept)' = logit(1e-5),
                          'bmi' = 0.08,
                          'hgba1c' = 1.2),
       'eligibility' = 'pre_diabetes',  
       'truth_model' = list('outcome_model' = outcome ~ surgery,
                            'censor_model' = censor ~ smoking_status),
       'analysis_models' = list(
         ## Main Effects Model
         ### Can't even estimate IPWR models
         'Y ~ A' = list('outcome_model' = outcome ~ surgery),
         
         ### Main Effects Model (w/ IPCW)
         'Y ~ A (IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                          'censor_model' = censor ~ smoking_status))
  )

write_rds(params, 'data/simulations/missing_data/inputs/sim_params_12.rds')


### Parameters to show proof of variance sim

### Simulation 13 (M-Bias)
params <- 
  list('sim_id' = 13,
       'study_design' = df_settings$study_design[9],
       'missingness' = df_settings$missingness[9],
       'n_subjects' = 5000,
       'n_trials' = 12,
       'n_sims' = 1000,
       'study_duration' = 36,
       'weight_models' = weight_models,
       'a1c_models' = a1c_models,
       'treatment_model' = list('(Intercept)' = logit(0.04),
                                'smoking_status[current]' = -2,
                                'smoking_status[former]' = -0.75),
       'bs_type_model' = list('(Intercept)' = logit(0.6),
                              'site' = log(2)),
       'outcome_model' = list('(Intercept)' = logit(0.015),
                              'elix_score' = 1.2,
                              'surgery' = log(0.7)),
       'cc_model' =  list('(Intercept)' = logit(0.08),
                          'elix_score' = 0.8,
                          'site' = 0.5,
                          'smoking_status[former]' = -0.5,
                          'smoking_status[current]' = -1),
       'eligibility' = 'pre_diabetes',
       'truth_model' = list('outcome_model' = outcome ~ surgery,
                            'censor_model' = censor ~ smoking_status),
       'analysis_models' = list(
         
         'Y ~ A (IPWR: R ~ L^{R}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                           'cc_model' = R ~ smoking_status + site + elix_score,
                                                           'ipwr_purpose' = 'Estimation',
                                                           'censor_model' = censor ~ smoking_status)
       )
  )

write_rds(params, 'data/simulations/missing_data/inputs/sim_params_13.rds')


### Simulation 14 (Treatment Effect Heterogeneity)
params <- 
  list('sim_id' = 14,
       'study_design' = df_settings$study_design[10],
       'missingness' = df_settings$missingness[10],
       'n_subjects' = 5000,
       'n_trials' = 12,
       'n_sims' = 1000,
       'study_duration' = 36,
       'weight_models' = weight_models,
       'a1c_models' = a1c_models,
       'treatment_model' = list('(Intercept)' = logit(0.04),
                                'smoking_status[current]' = -2,
                                'smoking_status[former]' = -0.75),
       'bs_type_model' = list('(Intercept)' = logit(0.6),
                              'site' = log(2)),
       'outcome_model' = list('(Intercept)' = logit(0.025),
                              'elix_score' = 0.25,
                              'surgery:elix_score' = 0.1,
                              'surgery' = log(0.7)),
       'cc_model' =  list('(Intercept)' = logit(0.08),
                          'elix_score' = 0.5,
                          'site' = 0.25,
                          'smoking_status[former]' = -0.5,
                          'smoking_status[current]' = -1),
       'eligibility' = 'pre_diabetes',
       'truth_model' = list('outcome_model' = outcome ~ surgery,
                            'censor_model' = censor ~ smoking_status),
       'analysis_models' = list(
         'Y ~ A (IPWR: R ~ L^{R}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                           'cc_model' = R ~ smoking_status + site + elix_score,
                                                           'ipwr_purpose' = 'Estimation',
                                                           'censor_model' = censor ~ smoking_status)
       )
  )

write_rds(params, 'data/simulations/missing_data/inputs/sim_params_14.rds')

### Simulation 15 (M-Bias Through Mediator)
params <- 
  list('sim_id' = 15,
       'study_design' = df_settings$study_design[11],
       'missingness' = df_settings$missingness[11],
       'n_subjects' = 5000,
       'n_trials' = 12,
       'n_sims' = 1000,
       'study_duration' = 36,
       'weight_models' = weight_models_3,
       'a1c_models' = a1c_models_3,
       'treatment_model' = list('(Intercept)' = logit(0.04),
                                'smoking_status[current]' = -2,
                                'smoking_status[former]' = -0.75),
       'bs_type_model' = list('(Intercept)' = logit(0.6),
                              'site' = log(2)),
       'outcome_model' = list('(Intercept)' = logit(0.0001),
                              'bmi' = 0.105,
                              'hgba1c' = 0.185),
       'cc_model' =   list('(Intercept)' = logit(0.25),
                           'elix_score' = 0.6,
                           'site' = -0.55,
                           'insulin' = 0.5),
       'eligibility' = 'pre_diabetes',  
       'truth_model' = list('outcome_model' = outcome ~ surgery,
                            'censor_model' = censor ~ smoking_status),
       'analysis_models' = list(
         'Y ~ A (IPWR: R ~ L^{R}, IPWC: C ~ L^{A})' = list('outcome_model' = outcome ~ surgery,
                                                           'cc_model' = R ~  site + elix_score + insulin,
                                                           'ipwr_purpose' = 'Estimation',
                                                           'censor_model' = censor ~ smoking_status)
       )
  )

write_rds(params, 'data/simulations/missing_data/inputs/sim_params_15.rds')
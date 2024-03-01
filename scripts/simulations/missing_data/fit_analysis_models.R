library(tidyverse)

fit_analysis_models <- function(df_trials, df_traj, params, models, interaction) {
  
  ### Inverse Probability of Censoring Weights
  if(!is.null(models$censor_model)) {
    ### Fit censor model to estimate IPCW weights
    ### Can only be censored after first month and if you switch from surgery
    ### to no surgery (not the other way around)
    df_censor <-
      df_trials %>%
      filter(surgery == 0, follow_up > 0)
    
    if(all(df_censor$smoking_status != 'no_self_report')) {
      df_trials$smoking_status[df_trials$smoking_status == 'no_self_report'] <- 'never' 
    }
    
    censor_model <-
      glm(formula = models$censor_model,
          family = 'binomial',
          data = df_censor)
    
    df_trials <-
      df_trials %>%
      mutate('prob_censor' = predict(censor_model, newdata = df_trials, type = 'response')) %>%
      mutate('prob_censor' = replace(prob_censor, surgery == 1 | follow_up == 0, 0)) %>%
      mutate('prob_censor_sw' = mean(df_censor$censor)) %>% ### Numerator weight stabilization
      mutate('prob_censor_sw' = replace(prob_censor_sw, surgery == 1 | follow_up == 0, 0)) %>%
      group_by(trial_id, subject_id) %>%
      arrange(follow_up) %>%
      mutate('cum_prob_uncensored' = cumprod(1 - prob_censor)) %>%
      mutate('cum_prob_uncensored_sw' = cumprod(1 - prob_censor_sw)) %>%
      ungroup() %>%
      mutate('ipcw' = cum_prob_uncensored_sw/cum_prob_uncensored)
    
  } else {
    df_trials$ipcw <- 1
  }
  
  if('censor' %in% names(df_trials)) {
    df_trials <-
      df_trials %>% 
      filter(!censor) 
  }
  
  ### Inverse probability of Selection Weights
  if(!is.null(models$cc_model) & !replace_null(models$ipwr_non_surgery, F)) {
    cc_model <- 
      glm(models$cc_model, 
          data = df_traj,
          family = 'binomial')
    
    df_ipw_R <- 
      df_traj %>% 
      mutate('prob_R' = predict(cc_model, newdata = ., type = 'response')) %>% 
      mutate('ipw_R' = mean(R)/prob_R) %>% ### Weight Stabilization
      select(subject_id, time, ipw_R)
    
    df_trials <- 
      df_trials %>% 
      inner_join(df_ipw_R, by = c('subject_id', 'time')) %>% 
      arrange(time) %>% 
      group_by(subject_id, trial_id) %>% 
      mutate(ipw_R = first(ipw_R)) %>% 
      ungroup()
    
    
  } else if(!is.null(models$cc_model) & replace_null(models$ipwr_non_surgery, F)) {
    ### Non-Surgery IPWR Only
    cc_model <- 
      glm(models$cc_model, 
          data = df_traj %>% filter(surgery == 0),
          family = 'binomial')
    
    df_ipw_R <- 
      df_traj %>% 
      mutate('prob_R' = predict(cc_model, newdata = ., type = 'response')) %>% 
      mutate('ipw_R' = mean(R[surgery == 0])/prob_R) %>% ### Weight Stabilization
      mutate('ipw_R' = ifelse(surgery == 1, 1, ipw_R)) %>% 
      select(subject_id, time, ipw_R)
    
    df_trials <- 
      df_trials %>% 
      inner_join(df_ipw_R, by = c('subject_id', 'time')) %>% 
      arrange(time) %>% 
      group_by(subject_id, trial_id) %>% 
      mutate(ipw_R = first(ipw_R)) %>% 
      ungroup()
    
  } else {
    df_trials$ipw_R <- 1
  }
  
  
  
  if(!interaction) {
    ### IPW = IPCW x IPWR
    df_trials$ipw <- df_trials$ipcw * df_trials$ipw_R
    
    ### Fit model
    outcome_model <- 
      glm(formula = models$outcome_model,
          family = 'binomial',
          weights = ipw,
          data = df_trials) 
    
    results <- 
      broom::tidy(outcome_model) %>% 
      select(term, estimate) 
  } else {
    ### No longer using this but note for if needed to use g-formula to 
    ### compute a marginal effect
    if(is.null(models$ipwr_purpose) | replace_null(models$ipwr_purpose) != 'Marginalization') {
      
      ### Use IPWR as model fitting weights rather than marginalization weights
      ### IPW = IPCW x IPWR
      df_trials$ipw <- df_trials$ipcw * df_trials$ipw_R
      df_trials$ipw_R <- 1
      
    } else if(models$ipwr_purpose == 'Marginalization') {
      ### IPW = IPCW 
      df_trials$ipw <- df_trials$ipcw 
    }
    
    ### Fit model
    outcome_model <- 
      glm(formula = models$outcome_model,
          family = 'binomial',
          weights = ipw,
          data = df_trials) 
    
    ### If interaction, use gformula
    df0 <- df1 <- df_trials
    df0$bs_type <- 0
    df0$surgery <- 0
    df1$bs_type <- 1
    df1$surgery <- 1
    
    
    df_trials$mu0_hat <-  predict(outcome_model, newdata = df0)  
    df_trials$mu1_hat <-  predict(outcome_model, newdata = df1)  
    results <- 
      tibble('estimate' = weighted.mean(df_trials$mu1_hat - df_trials$mu0_hat, df_trials$ipw_R),
             'term' =  ifelse(params$study_design == '(1) RYGB vs. VSG',
                              'bs_type', 
                              'surgery'))
  }
  
  return(results)
}

fit_analysis_models_bootstrap <- function(df_trials, df_traj, params, models, interaction) {
  
  ### Inverse Probability of Censoring Weights
  if(!is.null(models$censor_model)) {
    ### Fit censor model to estimate IPCW weights
    ### Can only be censored after first month and if you switch from surgery
    ### to no surgery (not the other way around)
    df_censor <-
      df_trials %>%
      filter(surgery == 0, follow_up > 0)
    
    if(all(df_censor$smoking_status != 'no_self_report')) {
      df_trials$smoking_status[df_trials$smoking_status == 'no_self_report'] <- 'never' 
    }
    
    censor_model <-
      glm(formula = models$censor_model,
          family = 'binomial',
          data = df_censor)
    
    df_trials <-
      df_trials %>%
      mutate('prob_censor' = 0) %>% 
      mutate('prob_censor' = replace(prob_censor, surgery == 0 & follow_up > 0, censor_model$fitted.values)) %>%
      mutate('prob_censor_sw' = mean(df_censor$censor)) %>% ### Numerator weight stabilization
      mutate('prob_censor_sw' = replace(prob_censor_sw, surgery == 1 | follow_up == 0, 0)) %>%
      group_by(trial_id, subject_id) %>%
      arrange(follow_up) %>%
      mutate('cum_prob_uncensored' = cumprod(1 - prob_censor)) %>%
      mutate('cum_prob_uncensored_sw' = cumprod(1 - prob_censor_sw)) %>%
      ungroup() %>%
      mutate('ipcw' = cum_prob_uncensored_sw/cum_prob_uncensored)
    
  } else {
    df_trials$ipcw <- 1
  }
  
  if('censor' %in% names(df_trials)) {
    df_trials <-
      df_trials %>% 
      filter(!censor) 
  }
  
  ### Inverse probability of Selection Weights
  if(!is.null(models$cc_model) & !replace_null(models$ipwr_non_surgery, F)) {
    cc_model <- 
      glm(models$cc_model, 
          data = df_traj,
          family = 'binomial')
    
    df_ipw_R <- 
      df_traj %>% 
      mutate('prob_R' = predict(cc_model, newdata = ., type = 'response')) %>% 
      mutate('ipw_R' = mean(R)/prob_R) %>% ### Weight Stabilization
      select(subject_id, time, replicate_id, ipw_R)
    
    df_trials <- 
      df_trials %>% 
      inner_join(df_ipw_R, by = c('subject_id', 'replicate_id', 'time')) %>% 
      arrange(time) %>% 
      group_by(subject_id, replicate_id, trial_id) %>% 
      mutate(ipw_R = first(ipw_R)) %>% 
      ungroup()
    
    
  } else if(!is.null(models$cc_model) & replace_null(models$ipwr_non_surgery, F)) {
    ### Non-Surgery IPWR Only
    cc_model <- 
      glm(models$cc_model, 
          data = df_traj %>% filter(surgery == 0),
          family = 'binomial')
    
    df_ipw_R <- 
      df_traj %>% 
      mutate('prob_R' = predict(cc_model, newdata = ., type = 'response')) %>% 
      mutate('ipw_R' = mean(R[surgery == 0])/prob_R) %>% ### Weight Stabilization
      mutate('ipw_R' = ifelse(surgery == 1, 1, ipw_R)) %>% 
      select(subject_id, time, replicate_id, ipw_R)
    
    df_trials <- 
      df_trials %>% 
      inner_join(df_ipw_R, by = c('subject_id', 'replicate_id', 'time')) %>% 
      arrange(time) %>% 
      group_by(subject_id, replicate_id, trial_id,) %>% 
      mutate(ipw_R = first(ipw_R)) %>% 
      ungroup()
    
  } else {
    df_trials$ipw_R <- 1
  }
  
  
  
  if(!interaction) {
    ### IPW = IPCW x IPWR
    df_trials$ipw <- df_trials$ipcw * df_trials$ipw_R
    
    ### Fit model
    outcome_model <- 
      glm(formula = models$outcome_model,
          family = 'binomial',
          weights = ipw,
          data = df_trials) 
    
    results <- 
      broom::tidy(outcome_model) %>% 
      select(term, estimate) 
  } else {
    ### No longer using this but note for if needed to use g-formula to 
    ### compute a marginal effect
    if(is.null(models$ipwr_purpose) | replace_null(models$ipwr_purpose) != 'Marginalization') {
      
      ### Use IPWR as model fitting weights rather than marginalization weights
      ### IPW = IPCW x IPWR
      df_trials$ipw <- df_trials$ipcw * df_trials$ipw_R
      df_trials$ipw_R <- 1
      
    } else if(models$ipwr_purpose == 'Marginalization') {
      ### IPW = IPCW 
      df_trials$ipw <- df_trials$ipcw 
    }
    
    ### Fit model
    outcome_model <- 
      glm(formula = models$outcome_model,
          family = 'binomial',
          weights = ipw,
          data = df_trials) 
    
    ### If interaction, use gformula
    df0 <- df1 <- df_trials
    df0$bs_type <- 0
    df0$surgery <- 0
    df1$bs_type <- 1
    df1$surgery <- 1
    
    
    df_trials$mu0_hat <-  predict(outcome_model, newdata = df0)  
    df_trials$mu1_hat <-  predict(outcome_model, newdata = df1)  
    results <- 
      tibble('estimate' = weighted.mean(df_trials$mu1_hat - df_trials$mu0_hat, df_trials$ipw_R),
             'term' =  ifelse(params$study_design == '(1) RYGB vs. VSG',
                              'bs_type', 
                              'surgery'))
  }
  
  return(results)
}


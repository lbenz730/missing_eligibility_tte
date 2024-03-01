library(tidyverse)
library(splines)
library(mvtnorm)
source('scripts/simulations/helpers.R')

### List of coefficients that control trajectories
### control:
###   slope: vector of fixed effects for slope
###   quadratic: (BMI ONLY, NOT A1C)
###     prob_nonzero: vector of coefficients for probability of non-zero quadratic term
###     term: vector of fixed effects for quadratic term
###   std_dev:
###     slope: sd for slope random effects
###     quadratic: sd for quadratic random effetcs.
###.  errors: 
###     structure: structure of error vectors
###     std_dev: error sd if iid
###     error_function: function that maps baseline value to error SD if heteroskeastic
### surgery:
###   group_model: list of coefficients for multinomial group model
###   std_dev: SD for spline random effects
###   spline: The spline knots/coefficients
###   errors:
###     structure: structure of error vectors
###     std_dev: error sd if iid
###     error_function: function that maps baseline value to error SD if heteroskeastic


### Function to sample mean weight trajectories
###
### df = data frame of covariate values
### coeff = list of trajectory model coefficients
### surgery = logical if we are to sample from surgery trajectories or 
### control trajectories
sample_weight_trajectories <- function(df, coeff, surgery) {
  ### Trajectories for Controls = 
  ### Linear term + quadratic slope w/ some probability
  if(!surgery) {
    df_trajectory <- 
      tibble('subject_id' = df$subject_id, 
             'slope_fixed' = NA,
             'slope_random' = NA,
             'delta' = NA,
             'quadratic_fixed' = NA,
             'quadratic_random' = NA)
    
    ### (1) Compute slope fixed effects and sample slope random effects
    df_trajectory$slope_fixed <- compute_model(df, beta = coeff$control$slope)
    df_trajectory$slope_random <- 
      rnorm(n = nrow(df), 
            mean = 0, 
            sd = coeff$control$std_dev$slope)
    
    ### (2) Sample Probability of Quadratic term being non-zero
    p_nonzero <- expit(compute_model(df, beta = coeff$control$quadratic$prob_nonzero))
    df_trajectory$delta <- rbinom(n = nrow(df), size = 1, p = p_nonzero)
    
    ### (3) Compute quadratic fixed and sample quadratic random effects
    df_trajectory$quadratic_fixed <- 
      compute_model(df, beta = coeff$control$quadratic$term)
    df_trajectory$quadratic_random <-
      rnorm(n = nrow(df),
            mean = 0,
            sd = coeff$control$std_dev$quadratic)
    
  } else {
    ### Trajectories for surgery
    ### Based off some splines
    
    ### (1) Sample Group for weight trajectory
    probs <- compute_multinomial_model(df, beta = coeff$surgery$group_model)
    probs <- expit(probs)
    group_probs <- cbind(probs, NA)
    for(i in 1:nrow(probs)) {
      group_probs[i, ] <- c(probs[i,1], lead(probs[i,], default = 1) - probs[i,])
    }
    
    group <- glmnet::rmult(group_probs)
    
    
    ### (2) Start all trajectories at the spline corresponding to the group
    M <- matrix(NA, nrow = nrow(df), ncol = length(coeff$surgery$spline$coeff[[1]]))
    colnames(M) <- paste0('knot_', 1:length(coeff$surgery$spline$coeff[[1]]))
    df_trajectory <- as_tibble(M)
    
    for(i in 1:length(coeff$surgery$spline$coeff)) {
      ix <- which(group == i)
      df_trajectory[ix,] <- 
        matrix(coeff$surgery$spline$coeff[[i]], 
               nrow = length(ix), 
               ncol = length(coeff$surgery$spline$coeff[[i]]),
               byrow = T)
    }
    
    df_trajectory$subject_id <- df$subject_id
    df_trajectory$group_id <- group
    df_trajectory <- 
      df_trajectory %>% 
      select(subject_id, group_id, everything())
    
    ### (3) Sample random effects around group trajectory
    df_trajectory$random_effect <- 
      rnorm(n = nrow(df), 
            mean = 0, 
            sd = coeff$surgery$std_dev)
    
  }
  
  return(df_trajectory)
}

### Function to sample mean a1c trajectories
###
### df = data frame of covariate values
### coeff = list of trajectory model coefficients
### surgery = logical if we are to sample from surgery trajectories or 
### control trajectories
sample_a1c_trajectories <- function(df, coeff, surgery) {
  if(!surgery) {
    df_trajectory <- 
      tibble('subject_id' = df$subject_id, 
             'slope_fixed' = NA,
             'slope_random' = NA)
    
    ### (1) Compute slope fixed effects and sample slope random effects
    ### Unlike weight trajectories, we don't have quadratic effects
    df_trajectory$slope_fixed <- compute_model(df, beta = coeff$control$slope)
    df_trajectory$slope_random <- 
      rnorm(n = nrow(df), 
            mean = 0, 
            sd = coeff$control$std_dev)
    
  } else {
    ### Trajectories for surgery
    ### Based off some splines
    
    ### (1) Sample Group for a1c trajectory
    probs <- compute_multinomial_model(df, beta = coeff$surgery$group_model)
    probs <- expit(probs)
    group_probs <- cbind(probs, NA)
    for(i in 1:nrow(probs)) {
      group_probs[i, ] <- c(probs[i,1], lead(probs[i,], default = 1) - probs[i,])
    }
    
    group <- glmnet::rmult(group_probs)
    
    
    ### (2) Start all trajectories at the spline corresponding to the group
    M <- matrix(NA, nrow = nrow(df), ncol = length(coeff$surgery$spline$coeff[[1]]))
    colnames(M) <- paste0('knot_', 1:length(coeff$surgery$spline$coeff[[1]]))
    df_trajectory <- as_tibble(M)
    
    for(i in 1:length(coeff$surgery$spline$coeff)) {
      ix <- which(group == i)
      df_trajectory[ix,] <- 
        matrix(coeff$surgery$spline$coeff[[i]], 
               nrow = length(ix), 
               ncol = length(coeff$surgery$spline$coeff[[i]]),
               byrow = T)
    }
    
    df_trajectory$subject_id <- df$subject_id
    df_trajectory$group_id <- group
    df_trajectory <- 
      df_trajectory %>% 
      select(subject_id, group_id, everything())
    
    ### (3) Sample random effects around group trajectory
    df_trajectory$random_effect <- 
      rnorm(n = nrow(df), 
            mean = 0, 
            sd = coeff$surgery$std_dev)
    
  }
  
  return(df_trajectory)
}

### Function to query a weight/a1c trajectory at any time point
### df_tracjectory: A df of subject specific coefficients for trajectories
### time_points: A vector of time points (in months) to query subject weight at 
### surgery = LOGICAL of whether or not trajectory is a surgery trajectory or not
### error: error variance matrix to draw random errors from
### spline_knots: knots in spline
### spline_time_points: Fit the spline basis function on
### a1c = indicator if we are querying trajectory for a1c (T) or weight (F)
###       the reason for the difference is that control trajectories for weight
###       are linear/quadratic model based while a1c trends are spline based
### baseline_vals = vector of baseline values (if heteroskedastic errors)

query_trajectories <- function(df_trajectory, 
                               time_points, 
                               surgery, 
                               a1c, 
                               error, 
                               spline_time_points = 0:60, 
                               spline_knots = NULL,
                               baseline_vals = NULL) {
  
  
  
  ### Build covariance matrix for errors
  if(error$structure == 'iid') {
    Sigma <- diag(rep(error$std_dev^2, length(time_points))) 
    ### 0 so that the first entry is always 0 at time 0 
    Sigma[time_points == 0, time_points == 0] <- 0
    
    ### Sample Errors
    errors <- 
      rmvnorm(n = nrow(df_trajectory),
              mean = rep(0, length(time_points)),
              sigma = Sigma)
    rownames(errors) <- df_trajectory$subject_id
    colnames(errors) <- paste0('time_', time_points)
    
  } else if(error$structure == 'heteroskedastic_by_a1c') {
    errors <- 
      rnorm(n = length(time_points) * length(baseline_vals),
            mean = 0,
            sd = rep(error$error_function(baseline_vals), each = length(time_points)))
    errors <- 
      matrix(errors, 
             nrow = length(baseline_vals), 
             ncol = length(time_points), 
             byrow = T)
    errors[,1] <- 0 ### 0 so that the first entry is always 0 at time 0 
    rownames(errors) <- df_trajectory$subject_id
    colnames(errors) <- paste0('time_', time_points)
  }
  
  ### Covert to df
  df_errors <- 
    errors %>% 
    as_tibble(rownames = 'subject_id') %>% 
    pivot_longer(cols = contains('time_'),
                 names_to = 'time',
                 values_to = 'error',
                 names_prefix = 'time_') %>% 
    mutate('time' = as.numeric(time),
           'subject_id' = as.numeric(subject_id))
  
  ### Weight Control Trends
  if(!surgery & !a1c) {
    df_query <- 
      crossing('time' = time_points,
               'subject_id' = df_trajectory$subject_id) %>% 
      inner_join(df_trajectory, by = 'subject_id') %>% 
      inner_join(df_errors, by = c('subject_id', 'time')) %>% 
      mutate('trend' = time * (slope_fixed + slope_random) + delta * (time^2) * (quadratic_fixed + quadratic_random)) %>% 
      mutate('observed' = trend + error)
  } else if(!surgery & a1c) {
    df_query <- 
      crossing('time' = time_points,
               'subject_id' = df_trajectory$subject_id) %>% 
      inner_join(df_trajectory, by = 'subject_id') %>% 
      inner_join(df_errors, by = c('subject_id', 'time')) %>% 
      mutate('trend' = time * (slope_fixed + slope_random)) %>% 
      mutate('observed' = trend + error)
  } else {
    ### Everything else (spline based)
    ### Use Matrix Algebra to compute trends quickly all at once
    X <- ns(spline_time_points, knots = spline_knots)
    traj_mat <-
      df_trajectory %>% 
      mutate_at(vars(contains('knot_')), ~{.x + random_effect}) %>% 
      select(contains('knot_')) %>% 
      as.matrix()
    trend <- X %*% t(traj_mat)
    colnames(trend) <- paste0('subject_', df_trajectory$subject_id)
    
    ### Covert back to df
    df_query <- 
      trend %>% 
      as_tibble() %>% 
      mutate('time' = spline_time_points) %>% 
      filter(time %in% time_points) %>% 
      pivot_longer(cols = contains('subject_'),
                   names_to = 'subject_id',
                   values_to = 'trend',
                   names_prefix = 'subject_') %>% 
      mutate('subject_id' = as.numeric(subject_id)) %>% 
      inner_join(df_errors, by = c('subject_id', 'time')) %>% 
      mutate('observed' = trend + error) %>% 
      inner_join(select(df_trajectory, subject_id, any_of('group_id')), by = 'subject_id')
  }
  
  return(df_query)
}



plot_trajectories <- function(df_query, plot_title, plot_subtitle, draw_mean = T, by_group = F, facet = F, weight = T) {
  y_lab <- ifelse(weight, '% Weight Change Since Baseline', '% A1C Change Since Baseline')
  df_query <- 
    df_query %>% 
    select(time, subject_id, any_of('group_id'), trend, observed) %>% 
    pivot_longer(cols = c('trend', 'observed'),
                 names_to = 'trend_type',
                 values_to = 'ptwc') %>% 
    mutate('trend_type' = case_when(trend_type == 'trend' ~ 'Subject Mean Trajectory (Fixed Effects + Random Effects)',
                                    trend_type == 'observed' ~ 'Observed Subject Trajectory (Fixed Effects + Random Effects + Random Error)')) %>% 
    mutate('trend_type' = fct_relevel(trend_type, 
                                      'Subject Mean Trajectory (Fixed Effects + Random Effects)',
                                      'Observed Subject Trajectory (Fixed Effects + Random Effects + Random Error)'))
  
  
  
  if(!by_group) {
    p <- 
      ggplot(df_query, aes(x = time, y = ptwc)) + 
      facet_wrap(~trend_type) + 
      geom_line(aes(group = subject_id), alpha = 0.1) + 
      scale_y_continuous(labels = scales::percent) + 
      scale_x_continuous(breaks = seq(0, max(df_query$time), 12), labels = ~.x/12) + 
      labs(x = 'Time since Baseline (Years)',
           y = y_lab,
           title = plot_title,
           subtitle = plot_subtitle)
    
    
    if(draw_mean) {
      df_trend <- 
        df_query %>% 
        filter(trend_type == 'Subject Mean Trajectory (Fixed Effects + Random Effects)') %>% 
        group_by(time) %>% 
        summarise('mean_trend' = mean(ptwc)) %>% 
        ungroup()
      
      p <- 
        p + 
        geom_line(data = df_trend, aes(y = mean_trend), col = 'red', lwd = 1.2)
    }
    
  } else {
    if(facet) {
      p <- 
        ggplot(df_query, aes(x = time, y = ptwc)) + 
        facet_grid(trend_type~paste('Group:', group_id)) +
        geom_line(aes(group = subject_id), alpha = 0.1) + 
        scale_y_continuous(labels = scales::percent) + 
        scale_x_continuous(breaks = seq(0, max(df_query$time), 12), labels = ~.x/12) + 
        labs(x = 'Time since Baseline (Years)',
             y = y_lab,
             title = plot_title,
             subtitle = plot_subtitle) + 
        theme(legend.position = 'none')
    } else {
      p <- 
        ggplot(df_query, aes(x = time, y = ptwc)) + 
        geom_line(aes(group = subject_id, col = as.factor(group_id)), alpha = 0.1) + 
        scale_y_continuous(labels = scales::percent) + 
        scale_x_continuous(breaks = seq(0, max(df_query$time), 12), labels = ~.x/12) + 
        labs(x = 'Time since Baseline (Years)',
             y = y_lab,
             col = 'Group',
             title = plot_title,
             subtitle = plot_subtitle) 
    }
    
    if(draw_mean) {
      df_trend <- 
        df_query %>% 
        filter(trend_type == 'Subject Mean Trajectory (Fixed Effects + Random Effects)') %>% 
        group_by(time, group_id) %>% 
        summarise('mean_trend' = mean(ptwc)) %>% 
        ungroup()
      
      p <- 
        p + 
        geom_line(data = df_trend, aes(y = mean_trend, col = as.factor(group_id)), lwd = 1.5) 
    }
    
    
    
  }
  
  return(p)
}

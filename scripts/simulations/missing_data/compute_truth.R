library(tidyverse)

### Function to compute the true ATE 
compute_truth <- function(dfs_trajectories, params) {
  target_parameter <- 
    ifelse(params$study_design == '(1) RYGB vs. VSG',
           'bs_type', 
           'surgery')
  
  cat('Prepare Analysis Datasets\n')
  t0 <- Sys.time()
  dfs_trials <- 
    map(1:n_batches, ~{
      cat('Computing batch:', .x, 'of', n_batches, '\n')
      map(dfs_trajectories[ batches[[.x]] ], function(dataset) {
        build_analysis_dataset(df = dataset, params)
      })
    })
  dfs_trials <- flatten(dfs_trials)
  t1 <- Sys.time()
  cat('Completed in:', t1-t0, '\n')
  
  cat('Fit Analysis Models\n')
  t0 <- Sys.time()
  if(params$study_design == '(1) RYGB vs. VSG') {
    df_results <- 
      map_dfr(dfs_trials, 
              ~fit_analysis_models(df_trials = .x, 
                                   df_traj = NULL,
                                   params = params,
                                   models = params$truth_model,
                                   interaction = F), 
              .id = 'dataset_id')
  } else {
    df_results <- 
      map_dfr(1:n_batches, ~{
        cat('Computing batch:', .x, 'of', n_batches, '\n')
        map_dfr(dfs_trials[ batches[[.x]] ], function(dataset) {
          fit_analysis_models(df_trials = dataset,
                              df_traj = NULL,
                              params = params,
                              models = params$truth_model,
                              interaction = F)
        }, 
        .id = 'dataset_id') %>% 
          mutate('dataset_id' = as.character(batches[[.x]][as.numeric(dataset_id)]))
      })
    
  }
  t1 <- Sys.time()
  cat('Completed in:', t1-t0, '\n')
  
  ### True ATE
  df_truth <- 
    df_results %>% 
    summarise('true_ate' = mean(estimate[term == target_parameter]))
  
  return(df_truth)
  
}

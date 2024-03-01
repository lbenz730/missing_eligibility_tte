library(tidyverse)
library(arrow)
library(glue)
library(furrr)

source('scripts/simulations/helpers.R')
source('scripts/simulations/missing_data/weight_trajectories_functions.R')
source('scripts/simulations/missing_data/generate_data.R')
source('scripts/simulations/missing_data/fit_analysis_models.R')
source('scripts/simulations/missing_data/compute_truth.R')

### Set up Parallelization
n_cores <- 12
plan(future::multisession(workers = n_cores))
options(future.globals.maxSize = 100 * 1024^3)
set.seed(7123)

### Simulation Parameters
args <- commandArgs(trailingOnly = T)
sim_id <- as.numeric(args[1])
params <- read_rds(glue('data/simulations/missing_data/inputs/sim_params_{sim_id}.rds'))
batch_size <- ifelse(params$study_design == '(1) RYGB vs. VSG', 100, 25)
batches <- split(1:params$n_sims, ceiling(1:params$n_sims / batch_size))
n_batches <- length(batches)

### Read in data and inputs
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'
df_pop <- 
  read_parquet(glue('{data_dir}/t2dm_population.parquet')) %>% 
  mutate('race' = ifelse(race == 'WH', 1, 0),
         'gender' = ifelse(gender == 'F', 1, 0),
         'site' = ifelse(site == 'SC', 0, 1)) %>% 
  select(subject_id, baseline_bmi, baseline_hgba1c, baseline_age,
         any_of(clean_names(params, names(.))))

cat('Generating Weight Trajectories and Outcomes\n')
t0 <- Sys.time()
dfs_trajectories <- 
  future_map(1:params$n_sims, 
             ~generate_data(params),
             .options = furrr_options(seed = 71283, chunk_size = 10))
t1 <- Sys.time()
cat('Completed in:', t1-t0, '\n')

### Compute the true ate given this configuration
cat('Computing True ATE\n')
df_truth <- compute_truth(dfs_trajectories, params)

### Weight Trajectories and Outcomes (on Absolute Time Scale)
### Note that this still includes inelegible people
### At some point may want to compare w/ methods that use imputation when 
### ascertainment of the eligibility criteria is unavailable 
cat('Imposing Missingness in Eligibility Data\n')
t0 <- Sys.time()
dfs_trajectories <- 
  map(1:n_batches, ~{
    cat('Computing batch:', .x, 'of', n_batches, '\n')
    map(dfs_trajectories[ batches[[.x]] ], function(dataset) {
      impose_missingness(df = dataset, params)
    })
  })
dfs_trajectories <- flatten(dfs_trajectories)
t1 <- Sys.time()
cat('Completed in:', t1-t0, '\n')

### Prepare datasets in target trial format
### Note this step filters to only eligible subjects
### If we want to do some kind of imputation of missing eligibility 
### we will want to add it before this step
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

### Fit the analysis models (via looping over different ones to fit)
cat('Fit Analysis Models\n')
t0 <- Sys.time()
df_results <- NULL
for(i in 1:length(params$analysis_models)) {
  model_name <- names(params$analysis_models)[i]
  cat(model_name, '\n')
  df_tmp <- 
    future_map2_dfr(dfs_trials, 
                    dfs_trajectories,
                   ~fit_analysis_models(df_trials = .x, 
                                        df_traj = .y,
                                        params = params, 
                                        models = params$analysis_models[[i]],
                                        interaction = replace_null(params$analysis_models[[i]]$interaction, replace_string = F)), 
                   .id = 'dataset_id',
                   .options = furrr_options(seed = 901273, chunk_size = 10)) %>% 
    mutate('analysis_method' = model_name,
           'outcome_model' = gsub('\\s+\\(.*', '', model_name),
           'ipwr_model' = ifelse(!grepl('IPWR', model_name), 
                                 '---', 
                                 gsub('^.*\\(IPWR:\\s+', '', gsub('\\)$', '', gsub(',\\sIPWC.*', '', model_name)))),
           'ipwc_model' = ifelse(!grepl('IPWC', model_name), 
                                 '---', 
                                 gsub('^.*IPWC:\\s+', '', gsub('\\)$', '', model_name))),
           'ipwr_purpose' = replace_null(params$analysis_models[[i]]$ipwr_purpose, '---'),
           'ipcw_model' = !is.null(params$analysis_models[[i]]$censor_model))  

  ### Save results
  df_results <- 
    df_results %>% 
    bind_rows(df_tmp)
}
t1 <- Sys.time()
cat('Completed in:', t1-t0, '\n')

### Save Results
keep_vars <-
  c('sim_id', 'study_design', 'missingness',
    'n_trials', 'n_sims', 'n_subjects', 'study_duration')

target_parameter <- 
  ifelse(params$study_design == '(1) RYGB vs. VSG',
         'bs_type', 
         'surgery')

df_results <- 
  bind_cols(as_tibble(params[keep_vars]))  %>% 
  bind_cols(df_results) %>% 
  filter(term == target_parameter) %>% 
  rename('estimated_ate' = estimate) %>% 
  bind_cols(df_truth) 


write_parquet(df_results, glue('data/simulations/missing_data/outputs/sim_results_{sim_id}.parquet'))

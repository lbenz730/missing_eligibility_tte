library(tidyverse)
library(arrow)
library(glue)
library(furrr)
library(data.table)

source('scripts/simulations/helpers.R')
source('scripts/simulations/missing_data/weight_trajectories_functions.R')
source('scripts/simulations/missing_data/generate_data.R')
source('scripts/simulations/missing_data/fit_analysis_models.R')
source('scripts/simulations/missing_data/compute_truth.R')

### Set up Parallelization
n_cores <- 48
plan(future::multisession(workers = n_cores))
options(future.globals.maxSize = 100 * 1024^3)
set.seed(7123)

### Simulation Parameters
args <- commandArgs(trailingOnly = T)
boot_id <- as.numeric(args[1])
sim_id <- as.numeric(args[2])
params <- read_rds(glue('data/simulations/missing_data/inputs/sim_params_{sim_id}.rds'))
batch_size <- ifelse(params$study_design == '(1) RYGB vs. VSG', 100, 25)
batch_size <- ifelse(sim_id > 12, 250, batch_size)
batches <- split(1:params$n_sims, ceiling(1:params$n_sims / batch_size))
n_batches <- length(batches)
n_boot <- 1000

### Read in data and inputs
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'
df_truth <- read_parquet(glue('data/simulations/bootstrap_variance/sim_data/truth_{sim_id}.parquet'))
df_traj <- read_parquet(glue('data/simulations/bootstrap_variance/sim_data/traj_{sim_id}.parquet'))
dfs_trajectories <- 
  df_traj %>% 
  group_by(dataset_id) %>% 
  group_split()

### Weight Trajectories and Outcomes (on Absolute Time Scale)
### Note that this still includes inelegible people
### At some point may want to compare w/ methods that use imputation when 
### ascertainment of the eligibility criteria is unavailable 
dfs_trajectories <- 
  map(1:n_batches, ~{
    cat('Computing batch:', .x, 'of', n_batches, '\n')
    map(dfs_trajectories[ batches[[.x]] ], function(dataset) {
      impose_missingness(df = dataset, params)
    })
  })
dfs_trajectories <- flatten(dfs_trajectories)

### All trajectories combined 
dt_all_traj <- 
  map_dfr(dfs_trajectories, ~.x, .id = 'dataset_id') %>% 
  mutate('dataset_id' = as.numeric(dataset_id)) %>% 
  as.data.table()


subj_ids <- map(dfs_trajectories, ~unique(.x$subject_id))
seeds <- sample(1:100000000, n_boot, replace = F)


b <- boot_id
cat('Boostrap Replicate', b, '\n')
boot_subj_ids <- 
  future_map_dfr(subj_ids, 
                 ~tibble('subject_id' = sample(.x, size = params$n_subjects, replace = T)), 
                 .id = 'dataset_id',
                 .options = furrr_options(seed = seeds[b])) %>% 
  mutate('dataset_id' = as.numeric(dataset_id)) %>% 
  as.data.table()

dt_all <- dt_all_traj[boot_subj_ids, on = c('dataset_id', 'subject_id')]
dt_all[, 'replicate_id' := 1:.N, by = c('dataset_id', 'subject_id', 'time')]


dfs_traj_boot <- 
  dt_all %>% 
  group_by(dataset_id) %>% 
  group_split()

### Prepare Analysis Datasets
dfs_trials <- 
  map(1:n_batches, ~{
    map(dfs_traj_boot[ batches[[.x]] ], function(dataset) {
      build_analysis_dataset_bootstrap(df = dataset, params)
    })
  })
dfs_trials <- flatten(dfs_trials)

### Fit Analysis Model
model_name <- names(params$analysis_models)[1]
df_tmp <- 
  future_map2_dfr(dfs_trials, 
                  dfs_traj_boot,
                  ~fit_analysis_models_bootstrap(df_trials = .x, 
                                                 df_traj = .y,
                                                 params = params, 
                                                 models = params$analysis_models[[1]],
                                                 interaction = replace_null(params$analysis_models[[1]]$interaction, replace_string = F)), 
                  .id = 'dataset_id',
                  .options = furrr_options(seed = 901273)) %>% 
  mutate('dataset_id' = as.numeric(dataset_id)) %>% 
  mutate('bootstrap_id' = b,
         'true_ate' = df_truth$true_ate,
         'analysis_method' = model_name,
         'outcome_model' = gsub('\\s+\\(.*', '', model_name),
         'ipwr_model' = ifelse(!grepl('IPWR', model_name), 
                               '---', 
                               gsub('^.*\\(IPWR:\\s+', '', gsub('\\)$', '', gsub(',\\sIPWC.*', '', model_name)))),
         'ipwc_model' = ifelse(!grepl('IPWC', model_name), 
                               '---', 
                               gsub('^.*IPWC:\\s+', '', gsub('\\)$', '', model_name))),
         'ipwr_purpose' = replace_null(params$analysis_models[[1]]$ipwr_purpose, '---'),
         'ipcw_model' = !is.null(params$analysis_models[[1]]$censor_model))  


### Save Results
if(!dir.exists('data/simulations/bootstrap_variance/')) {
  dir.create('data/simulations/bootstrap_variance/')
  dir.create('data/simulations/bootstrap_variance/boot_level')
}

write_parquet(df_tmp, glue('data/simulations/bootstrap_variance/boot_level/file_{sim_id}_boot_{b}.parquet'))



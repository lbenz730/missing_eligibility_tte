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
sim_id <- as.numeric(args[1])
params <- read_rds(glue('data/simulations/missing_data/inputs/sim_params_{sim_id}.rds'))
batch_size <- ifelse(params$study_design == '(1) RYGB vs. VSG', 100, 25)
batch_size <- ifelse(sim_id > 12, 250, batch_size)
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
#df_truth <- compute_truth(dfs_trajectories, params)

if(!dir.exists('data/simulations/bootstrap_variance/sim_data')) {
  dir.create('data/simulations/bootstrap_variance/sim_data')
}

df_traj <- 
  map_dfr(dfs_trajectories, ~.x, .id = 'dataset_id') %>% 
  mutate('dataset_id' = as.numeric(dataset_id))
          
write_parquet(df_traj, glue('data/simulations/bootstrap_variance/sim_data/traj_{sim_id}.parquet'))
write_parquet(df_truth, glue('data/simulations/bootstrap_variance/sim_data/truth_{sim_id}.parquet'))

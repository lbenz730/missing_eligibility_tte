library(tidyverse)
library(arrow)
library(glue)
library(furrr)

source('scripts/simulations/helpers.R')
source('scripts/simulations/missing_data/weight_trajectories_functions.R')
source('scripts/simulations/missing_data/generate_data.R')

### Set up Parallelization
n_cores <- 12
plan(future::multisession(workers = n_cores))
options(future.globals.maxSize = 100 * 1024^3)

sim_id <- 3
params <- read_rds(glue('data/simulations/missing_data/inputs/sim_params_{sim_id}.rds'))
params$n_sims <- 100
params$n_subjects <- 10000

### Read in data and inputs
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'
df_pop <- 
  read_parquet(glue('{data_dir}/t2dm_population.parquet')) %>% 
  mutate('race' = ifelse(race == 'WH', 1, 0),
         'gender' = ifelse(gender == 'F', 1, 0),
         'site' = ifelse(site == 'SC', 0, 1)) %>% 
  select(subject_id, baseline_bmi, baseline_hgba1c, baseline_age,
         any_of(clean_names(params, names(.))))

dfs_trajectories <- 
  future_map(1:params$n_sims, 
             ~generate_data(params),
             .options = furrr_options(seed = 71283, chunk_size = 10),
             .progress = T)


dfs_trials <- 
  map(dfs_trajectories, function(dataset) {
    build_analysis_dataset(df = dataset, params)
  })

df_final <- 
  future_map_dfr(dfs_trials, ~{
    .x %>% 
      mutate('prob_R' = expit(compute_model(df = ., beta = params$cc_model))) %>% 
      mutate('R' = rbinom(n = nrow(.), size = 1, p = prob_R)) %>% 
      group_by(subject_id) %>% 
      mutate('R' = first(R)) %>% 
      ungroup()},
    .options = furrr_options(seed = 7112, chunk_size = 10),
    .id = 'dataset_id',
    .progress = T)


ggplot(df_final %>% filter(follow_up <= 36), aes(x = follow_up, y = bmi)) +
  facet_wrap(~ifelse(R == 1, 'Eligible Complete Case', 'Eligible but Missing Eligibility')) + 
  geom_smooth(aes(col = as.factor(bs_type))) + 
  labs(x = 'Time since Surgery (Months)',
       y = 'BMI',
       title = 'Comparison of BMI by Complete Case/Surgery Status',
       subtitle = 'Subset of Simulations/Patients\nSim ID: 3',
       col = 'Bariatric Surgery Type')
ggsave('figures/simulations/missing_data/eda/08_sim_3_bmi.png', height = 9/1.2, width = 16/1.2)
       

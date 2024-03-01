library(tidyverse)
library(arrow)
library(glue)
library(furrr)
library(knitr)

n_cores <- 12
plan(future::multisession(workers = n_cores))

df_results <- 
  future_map_dfr(dir('data/simulations/bootstrap_variance/boot_level/', full.names = T), ~{
    read_parquet(.x) %>% 
      mutate('sim_id' = gsub('_boot.*$', '', gsub('data/simulations/bootstrap_variance/boot_level//file_', '', .x))) %>% 
      mutate('sim_id' = as.numeric(sim_id)) %>% 
      mutate('dataset_id' = as.numeric(dataset_id))
  })

df_contain <-
  df_results %>% 
  filter(term == 'surgery') %>% 
  group_by(sim_id, dataset_id) %>% 
  summarise('n_boot' = n(),
            'truth' = mean(true_ate),
            'lower' = quantile(estimate, 0.025, na.rm = T),
            'mean' = mean(estimate, na.rm = T),
            'median' = median(estimate, na.rm = T),
            'upper' = quantile(estimate, 0.975, na.rm = T)) %>% 
  ungroup() %>% 
  mutate('interval_contain' = truth >= lower & truth <= upper)

df_contain %>% 
  group_by(sim_id) %>% 
  summarise('coverage' = mean(interval_contain),
            'n_boot' = min(n_boot),
            'mean' = mean(mean),
            'ate' = mean(truth)) %>% 
  mutate('setting' = case_when(sim_id == 13 ~ 'M-Bias',
                               sim_id == 14 ~ 'Treatment Effect Heterogeneity',
                               sim_id == 15 ~ 'M-Bias w/ Mediator')) %>% 
  select(setting, coverage) %>% 
  kable(align = 'c', format = 'latex')


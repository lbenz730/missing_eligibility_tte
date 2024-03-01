library(tidyverse)
library(glue)
library(arrow)
library(patchwork)
source('scripts/simulations/missing_data/weight_trajectories_functions.R') 
source('scripts/simulations/missing_data/generate_data.R')
set.seed(1312)
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Load in Saved out Surgery spline knots
spline_knots_bmi <- read_rds('models/eda/surgery_spline_weight.rds')
spline_knots_a1c <- read_rds('models/eda/surgery_spline_a1c.rds')
spline_knots_a1c$coeff <- spline_knots_a1c$coeff + c(rep(0, 5), 0.06)

### Augment this to be in the 5 groups
spline_knots_bmi$coeff <- 
  list('1' =  spline_knots_bmi$coeff + c(seq(0.01, 0.08, length.out = 5), 0.1, seq(0.1, 0.05, length.out =  5)),
       '2' = spline_knots_bmi$coeff + c(seq(0.005, 0.03, length.out = 5), 0.05, seq(0.05, 0.025, length.out =  5)),
       '3' = spline_knots_bmi$coeff,
       '4' = spline_knots_bmi$coeff -  seq(0.005, 0.08, length.out = 11),
       '5' = spline_knots_bmi$coeff - seq(0.01, 0.15, length.out = 11))

spline_knots_a1c$coeff <- 
  list('1' =  spline_knots_a1c$coeff + c(seq(0.01, 0.08, length.out = 6)),
       '2' = spline_knots_a1c$coeff + c(seq(0.005, 0.045, length.out = 6)),
       '3' = spline_knots_a1c$coeff ,
       '4' = spline_knots_a1c$coeff -  seq(0.005, 0.04, length.out = 6),
       '5' = spline_knots_a1c$coeff - seq(0.01, 0.08, length.out = 6) )

### Plot what the 5 a1c groups look like
time_points <- 0:60
X_a1c <- ns(time_points, knots = spline_knots_a1c$knots)
X_bmi <- ns(time_points, knots = spline_knots_bmi$knots)
df_groups_a1c <- 
  map_dfr(1:5, ~tibble('time' = time_points, 
                       'group' = .x,
                       'outcome' = 'A1c',
                       'ptwc' = as.vector( X_a1c %*% spline_knots_a1c$coeff[[.x]] )))
df_groups_bmi <- 
  map_dfr(1:5, ~tibble('time' = time_points, 
                       'group' = .x,
                       'outcome' = 'BMI',
                       'ptwc' = as.vector( X_bmi %*% spline_knots_bmi$coeff[[.x]] )))

df_groups <- 
  bind_rows(df_groups_a1c, df_groups_bmi)

p_latent <-
  ggplot(df_groups, aes(x = time, y = ptwc)) + 
  facet_wrap(~outcome) + 
  geom_line(aes(col = as.factor(group))) + 
  scale_x_continuous(breaks = seq(0, 60, 12)) + 
  scale_y_continuous(labels = scales::percent) +
  labs(x = 'Time Since Surgery (Months)',
       y = '% Change\nSince Surgery',
       title = 'Latent Post-Surgical Mean Trajectories',
       color = expression(paste('Latent Group for A1c (', xi[k], ')/BMI (', eta[k], ')')))


### Actual Values
df_pop <- 
  read_parquet(glue('{data_dir}/t2dm_population.parquet')) %>% 
  mutate('race' = ifelse(race == 'WH', 1, 0),
         'gender' = ifelse(gender == 'F', 1, 0),
         'site' = ifelse(site == 'SC', 1, 0))
params <- read_rds('data/simulations/missing_data/inputs/sim_params_5.rds')
params$study_duration <- 60
df_traj <- generate_data(params)

df_control <- 
  df_traj %>% 
  filter(surgery == 0) %>% 
  group_by(subject_id) %>% 
  mutate('BMI' = (bmi - first(bmi))/first(bmi),
         'A1c' = (hgba1c - first(hgba1c))/first(hgba1c)) %>% 
  ungroup() %>% 
  select(subject_id, time, BMI, A1c) %>% 
  pivot_longer(cols = c('BMI', 'A1c'),
               names_to = 'measure',
               values_to = 'value')

df_rel <- 
  df_traj %>% 
  group_by(subject_id) %>% 
  mutate('BMI' = (bmi - first(bmi))/first(bmi),
         'A1c' = (hgba1c - first(hgba1c))/first(hgba1c)) %>% 
  ungroup() %>% 
  select(subject_id, time, BMI, A1c) %>% 
  pivot_longer(cols = c('BMI', 'A1c'),
               names_to = 'measure',
               values_to = 'value')

df_surg <- 
  df_traj %>% 
  filter(surgery == 1) %>%
  group_by(subject_id) %>% 
  mutate('BMI' = (bmi - first(bmi))/first(bmi),
         'A1c' = (hgba1c - first(hgba1c))/first(hgba1c)) %>% 
  ungroup() %>% 
  select(subject_id, time, BMI, A1c, treatment_time) %>% 
  pivot_longer(cols = c('BMI', 'A1c'),
               names_to = 'measure',
               values_to = 'value')

p_control <- 
  ggplot(df_control, aes(x = time, y = value)) +
  facet_wrap(~measure) + 
  geom_line(aes(group = subject_id), alpha = 0.02) + 
  scale_x_continuous(breaks = seq(0, 60, 12)) + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = 'Time Since Start of Simulation (Months)',
       y = '% Change\nSince Start of Simulation',
       title = 'Sample Pre-Surgical Subject Trajectories')

p_surg <- 
  ggplot(df_surg, aes(x = time - treatment_time, y = value)) +
  facet_wrap(~measure) + 
  geom_line(aes(group = subject_id), alpha = 0.02) + 
  scale_x_continuous(breaks = seq(0, 60, 12)) + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = 'Time Since Surgery (Months)',
       y = '% Change\nSince Surgery',
       title = 'Sample Post-Surgical Subject Trajectories')


p_total <- 
  ggplot(df_rel, aes(x = time, y = value)) +
  facet_wrap(~measure) + 
  geom_line(aes(group = subject_id), alpha = 0.02) + 
  scale_x_continuous(breaks = seq(0, 60, 12)) + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = 'Time Since Start of Simulation (Months)',
       y = '% Change\nSince Start of Simulation',
       title = 'Sample Subject Trajectories')


(p_control + p_latent)/(p_surg + p_total) + 
  plot_annotation(tag_levels = 'A', tag_suffix = ')')  & 
  theme(plot.tag = element_text(size = 24))

ggsave('figures/paper1/simulation_map.png', height = 9, width = 16)
ggsave('figures/paper1/simulation_map.pdf', height = 9, width = 16)

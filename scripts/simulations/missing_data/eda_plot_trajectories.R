library(tidyverse)
library(glue)
library(arrow)
source('scripts/simulations/missing_data/weight_trajectories_functions.R') 
set.seed(1312)

### Load in Saved out Surgery spline knots
spline_knots <- read_rds('models/eda/surgery_spline_weight.rds')

### Augment this to be in the 5 groups
spline_knots$coeff <- 
  list('1' =  spline_knots$coeff + c(seq(0.01, 0.08, length.out = 5), 0.1, seq(0.1, 0.05, length.out =  5)),
       '2' = spline_knots$coeff + c(seq(0.005, 0.03, length.out = 5), 0.05, seq(0.05, 0.025, length.out =  5)),
       '3' = spline_knots$coeff,
       '4' = spline_knots$coeff -  seq(0.005, 0.08, length.out = 11),
       '5' = spline_knots$coeff - seq(0.01, 0.15, length.out = 11))

### Plot what the 5 Weight groups look like
time_points <- 0:60
X <- ns(time_points, knots = spline_knots$knots)
df_groups <- 
  map_dfr(1:5, ~tibble('time' = time_points, 
                       'group' = .x,
                       'ptwc' = as.vector( X %*% spline_knots$coeff[[.x]] )))

ggplot(df_groups, aes(x = time, y = ptwc)) + 
  geom_line(aes(col = as.factor(group))) + 
  scale_y_continuous(labels = scales::percent) + 
  scale_x_continuous(breaks = seq(0, max(time_points), 12), labels = ~.x/12) + 
  labs(x = 'Time since Baseline (Years)',
       y = '% Weight Change Since Baseline',
       title = 'Weight Loss Group Trends',
       subtitle = 'Surgery',
       color = 'Trajectory Group')
ggsave('figures/simulations/missing_data/eda/02_group_trend_weight_trajectories.png',
       height = 9/1.2,
       width = 16/1.2)


### List of coefficients that control trajectories
### control:
###   slope: vector of fixed effects for slope
###   quadratic:
###     prob_nonzero: vector of coefficients for probability of non-zero quadratic term
###     term: vector of fixed effects for quadratic term
###   std_dev:
###     slope: sd for slope random effects
###     quadratic: sd for quadratic random effetcs.
###.  errors: 
###     structure: structure of error vectors
###     std_dev: error sd if iid
###
### surgery:
###   group_model: vector of covariates for ordinal logistic regression group model
###   std_dev: sd for random effect 
###   spline: spline object 
###.  errors: 
###     structure: structure of error vectors
###     std_dev: error sd if iid

coeff <- 
  list('control' = list('slope' = list('(Intercept)' = -5e-4),
                        'quadratic' = list('prob_nonzero' = list('(Intercept)' = logit(0.25)),
                                           'term' = list('(Intercept)' = 0)),
                        'std_dev' = list('slope' = 5e-4,
                                         'quadratic' = 2e-5),
                        'errors' = list('structure' = 'iid',
                                        'std_dev' =  sqrt(1e-5))),
       
       'surgery' = list('group_model' = list('intercepts' = logit(cumsum(c(0.05, 0.2, 0.5, 0.2)))),
                        'std_dev' = 0.025,
                        'spline' = spline_knots,
                        'errors' = list('structure' = 'iid',
                                        'std_dev' = sqrt(1e-5)))
  )

### Toy dataset
df <- tibble('subject_id' = 1:1000)

### Control Example
df_traj_control <- sample_weight_trajectories(df, coeff, surgery = F) 
df_query_control <- 
  query_trajectories(df_trajectory = df_traj_control,
                     time_points = 0:60,
                     error = coeff$control$errors,
                     surgery = F,
                     a1c = F)

plot_trajectories(df_query = df_query_control,
                  plot_title = 'Sample Simulated Weight Trajectories',
                  plot_subtitle = 'Control',
                  draw_mean = T)

ggsave('figures/simulations/missing_data/eda/01_sample_weight_traj_control.png',
       height = 9/1.2,
       width = 16/1.2)

### Surgery Example
df_traj_surgery <- sample_weight_trajectories(df, coeff, surgery = T) 
df_query_surgery <- 
  query_trajectories(df_trajectory = df_traj_surgery,
                     time_points = 0:60,
                     surgery = T,
                     spline_knots = spline_knots$knots,
                     error = coeff$surgery$errors,
                     a1c = F)

plot_trajectories(df_query = df_query_surgery,
                  plot_title = 'Sample Simulated Weight Trajectories',
                  plot_subtitle = 'Surgery',
                  draw_mean = T)

ggsave('figures/simulations/missing_data/eda/03_sample_weight_traj_surgery.png',
       height = 9/1.2,
       width = 16/1.2)

plot_trajectories(df_query = df_query_surgery,
                  plot_title = 'Sample Simulated Weight Trajectories',
                  plot_subtitle = 'Surgery',
                  by_group = T,
                  draw_mean = T)

ggsave('figures/simulations/missing_data/eda/04_sample_weight_traj_surgery_by_group.png',
       height = 9/1.2,
       width = 16/1.2)

plot_trajectories(df_query = df_query_surgery,
                  plot_title = 'Sample Simulated Weight Trajectories',
                  plot_subtitle = 'Surgery',
                  by_group = T,
                  draw_mean = T,
                  facet = T) +
  theme(strip.text.y = element_text(size = 7))

ggsave('figures/simulations/missing_data/eda/05_sample_weight_traj_surgery_by_group_faceted.png',
       height = 9,
       width = 16)



### HgbA1c Trajectories
### Load in Saved out Surgery spline knots
spline_knots_a1c <- read_rds('models/eda/surgery_spline_a1c.rds')
spline_knots_a1c$coeff <- spline_knots_a1c$coeff + c(rep(0, 5), 0.06)

### Augment this to be in the 5 groups
spline_knots_a1c$coeff <- 
  list('1' =  spline_knots_a1c$coeff + c(seq(0.01, 0.08, length.out = 6)),
       '2' = spline_knots_a1c$coeff + c(seq(0.005, 0.045, length.out = 6)),
       '3' = spline_knots_a1c$coeff ,
       '4' = spline_knots_a1c$coeff -  seq(0.005, 0.04, length.out = 6),
       '5' = spline_knots_a1c$coeff - seq(0.01, 0.08, length.out = 6) )

### Plot what the 5 a1c groups look like
time_points <- 0:60
X <- ns(time_points, knots = spline_knots_a1c$knots)
df_groups <- 
  map_dfr(1:5, ~tibble('time' = time_points, 
                       'group' = .x,
                       'ptwc' = as.vector( X %*% spline_knots_a1c$coeff[[.x]] )))

ggplot(df_groups, aes(x = time, y = ptwc)) + 
  geom_line(aes(col = as.factor(group))) + 
  scale_y_continuous(labels = scales::percent) + 
  scale_x_continuous(breaks = seq(0, max(time_points), 12), labels = ~.x/12) + 
  labs(x = 'Time since Baseline (Years)',
       y = '% A1C Change Since Baseline',
       title = 'A1C Group Trends',
       subtitle = 'Surgery',
       color = 'Trajectory Group')
ggsave('figures/simulations/missing_data/eda/06_group_trend_a1c_trajectories.png',
       height = 9/1.2,
       width = 16/1.2)


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

### Toy dataset
df <- tibble('subject_id' = 1:1000)
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'
df_pop <- 
  read_parquet(glue('{data_dir}/t2dm_population.parquet')) %>% 
  mutate('race' = ifelse(race == 'WH', 1, 0),
         'gender' = ifelse(gender == 'F', 1, 0),
         'site' = ifelse(site == 'SC', 1, 0))
df$baseline <- sample(df_pop$baseline_hgba1c, 1000)


### Control Example
df_traj_control <- sample_a1c_trajectories(df, coeff = a1c_models, surgery = F) 
df_query_control <- 
  query_trajectories(df_trajectory = df_traj_control,
                     time_points = 0:60,
                     error = a1c_models$control$errors,
                     surgery = F,
                     a1c = T, 
                     baseline_vals = df$baseline)

plot_trajectories(df_query = df_query_control,
                  plot_title = 'Sample Simulated A1C Trajectories',
                  plot_subtitle = 'Control',
                  by_group = F,
                  draw_mean = T,
                  facet = F,
                  weight = F)

ggsave('figures/simulations/missing_data/eda/07_control_a1c_trajectories.png',
       height = 9/1.2,
       width = 16/1.2)


df_traj_surgery <- sample_a1c_trajectories(df, coeff = a1c_models, surgery = T) 
df_query_surgery <- 
  query_trajectories(df_trajectory = df_traj_surgery,
                     time_points = 0:48,
                     error = a1c_models$surgery$errors,
                     surgery = T,
                     a1c = T, 
                     baseline_vals = df$baseline,
                     spline_knots = a1c_models$surgery$spline$knots)

plot_trajectories(df_query = df_query_surgery,
                  plot_title = 'Sample Simulated A1C Trajectories',
                  plot_subtitle = 'Surgery',
                  by_group = F,
                  draw_mean = T,
                  facet = F,
                  weight = F)

ggsave('figures/simulations/missing_data/eda/08_surgery_a1c_trajectories.png',
       height = 9/1.2,
       width = 16/1.2)

plot_trajectories(df_query = df_query_surgery,
                  plot_title = 'Sample Simulated A1C Trajectories',
                  plot_subtitle = 'Surgery',
                  by_group = T,
                  draw_mean = T,
                  facet = T,
                  weight = F)

ggsave('figures/simulations/missing_data/eda/09_surgery_a1c_trajectories_faceted.png',
       height = 9/1.2,
       width = 16/1.2)



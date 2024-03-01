library(tidyverse)
library(arrow)
library(glue)
library(splines)
source('scripts/helpers.R')
  
df_trials <- read_parquet('data/trials_combined.parquet')

### Model for analyzing all data together
### Outcome Pooled logistic regression
### IPW Censoring? (can't have non-adherence for surgery group)

### Knots on quantiles of event time scale
knots <- 
  df_trials %>% 
  filter(pct_20_loss == 1) %>% 
  pull(follow_up_time) %>% 
  quantile(., probs = seq(0, 1, length.out = 7)[2:6])

cat('Outcome Pooled logistic regression\n')
outcome_model <-
  glm(pct_20_loss ~ 
        surgery * splines::ns(month, df = 6, knots) + site + 
        gender + race + site + age + baseline_bmi,
      family = 'binomial',
      y = F,
      model = F,
data = df_trials)

outcome_model <- strip_model(outcome_model)
write_rds(outcome_model, 'models/eda/outcome_pooled_LR.rds')


cat('Outcome Pooled logistic regression (for each trial seperately)\n')
for(i in 1:12) {
  cat('Trial', i, '\n')
  trial_outcome_model <- 
    glm(pct_20_loss ~ 
          surgery * splines::ns(month, df = 6, knots) + 
          site + gender + race + site + age + baseline_bmi,
        family = 'binomial',
        y = F,
        model = F,
        data = df_trials %>% filter(trial_id == i))
  trial_outcome_model <- strip_model(trial_outcome_model)
  write_rds(trial_outcome_model, glue('models/eda/outcome_pooled_LR_trial_{i}.rds'))
}




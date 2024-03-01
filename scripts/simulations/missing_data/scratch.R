
params$n_subjects <- 25000
df_traj <- generate_data(params)
df_trial <- build_analysis_dataset(df_traj, params)


# impose_missingness(x, params)
y <- 
  df_trial %>% 
  mutate('prob_R' = expit(compute_model(df = ., beta = params$cc_model))) %>% 
  mutate('R' = rbinom(n = nrow(.), size = 1, p = prob_R)) %>% 
  group_by(subject_id) %>% 
  mutate('R' = first(R)) %>% 
  ungroup()

### CC Rate
y %>% 
  filter(follow_up == 0) %>%
  pull(R) %>% 
  mean()

### Outcome Rate
mean(y$outcome)

y %>% 
  group_by(R, bs_type) %>% 
  summarise(n(),
            mean(bmi),
            mean(baseline_hgba1c),
            mean(outcome)) %>% 
  group_by(R) %>% 
  summarise('ate' = logit(`mean(outcome)`[bs_type == 1]) - logit(`mean(outcome)`[bs_type == 0]))


ggplot(y, aes(x = follow_up, y = bmi)) + 
  geom_smooth(aes(color = as.factor(paste(R, bs_type))))

ggplot(y, aes(x = follow_up, y = hgba1c)) + 
  geom_smooth(aes(color = as.factor(paste(R, bs_type))))




glm(outcome ~ bs_type, 
    family = 'binomial',
    data = y)

glm(outcome ~ bs_type, 
    family = 'binomial',
    data = y %>% filter(R == 1))
library(tidyverse)
library(arrow)
library(glue)
options(dplyr.summarise.inform = F)

source('scripts/util/helpers.R')

### Directory where EHR data is stored
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

df_files <- 
  crossing('bmi_lookback' = c(1,3,6,12),
           'a1c_lookback' = c(1,3,6,12,18,24))

df_missing <- 
  map_dfr(1:nrow(df_files), ~{
    read_parquet(glue('{data_dir}/microvascular_tte/pooled_trials/trials_{.x}.parquet')) %>% 
      mutate('bmi_lookback' = df_files$bmi_lookback[.x],
             'a1c_lookback' = df_files$a1c_lookback[.x]) %>% 
      group_by(bmi_lookback, a1c_lookback) %>% 
      summarise('n_subj_trials' = n(),
                'n_eligible' = sum(eligible == 1, na.rm = T),
                'n_ineligible' = sum(eligible == 0, na.rm = T),
                'n_missing' = sum(R == 0)) %>% 
      ungroup()
    
  })

df_summary <- 
  df_missing %>% 
  pivot_longer(cols = c('n_eligible',  'n_ineligible',  'n_missing'),
               names_to = 'category',
               values_to = 'n_obs',
               names_prefix = 'n_') %>% 
  mutate('pct' = n_obs/n_subj_trials) %>% 
  mutate('category' = tools::toTitleCase(category)) %>% 
  mutate('category' = ifelse(category == 'Missing', 'Missing Eligibility', category)) %>% 
  mutate('category_exp' = paste0('atop("', category, '","',
                                 case_when(category == 'Eligible' ~  '(E"[mk] == "1,"~"R"[mk] == "1)")',
                                           category == 'Ineligible' ~  '(E"[mk] == "0,"~"R"[mk] == "1)")',
                                           category == 'Missing Eligibility' ~ '(R"[mk] == "0)")')),
         'bmi_lookback_exp' = fct_reorder(paste0('"', paste('BMI Lookback:', bmi_lookback, ifelse(bmi_lookback == 1, 'Month', 'Months')), '"'), bmi_lookback))

ggplot(df_summary, aes(x = a1c_lookback, y = pct)) + 
  facet_grid(category_exp ~ bmi_lookback_exp, labeller = labeller(.cols = label_parsed, .rows = label_parsed, .multi_line = TRUE), scales = 'free_y') + 
  geom_line(aes(color = category)) + 
  geom_point(aes(color = category), size = 4) + 
  scale_x_continuous(breaks = seq(0,24,6)) +
  scale_y_continuous(labels = scales::percent) + 
  labs(x = 'Blood Glucose Labs Lookback Time (Months)',
       y = 'Proportion of Possible Subject-Trials',
       title = expression(paste('Distribution of Eligibility Status (', E[mk], ',', R[mk], ')')),
       subtitle = paste('Among', scales::number(df_summary$n_subj_trials[1], big.mark = ','), 'Subject-Trials')) + 
  theme(legend.position = 'none',
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 32),
        plot.subtitle = element_text(size = 24))

ggsave('figures/paper1/eligibility_distribution.png', height = 9/1.2, width = 16/1.2)
ggsave('figures/paper1/eligibility_distribution.pdf', height = 9/1.2, width = 16/1.2)


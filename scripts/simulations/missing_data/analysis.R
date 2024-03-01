library(tidyverse)
library(glue)
library(arrow)
library(ggridges)
library(gt)
library(gtExtras)

source('scripts/simulations/helpers.R')

df_settings <- read_csv('data/simulations/missing_data/inputs/settings.csv')
files <- dir(glue('data/simulations/missing_data/outputs'), full.names = T)
df_results <- 
  map_dfr(files, read_parquet) %>% 
  mutate('bias' = estimated_ate - true_ate, 
         'pct_bias' = (estimated_ate - true_ate)/true_ate) 

df_truth <- 
  df_results %>% 
  distinct(sim_id, study_design, missingness, true_ate)

df_analysis <- 
  df_results %>% 
  group_by(sim_id, study_design, missingness, outcome_model, ipwr_purpose,  ipwr_model, ipwc_model) %>% 
  summarise('true_ate' = first(true_ate),
            'std_dev' = sd(estimated_ate),
            'm_estimated_ate' = median(estimated_ate),
            'estimated_ate' = mean(estimated_ate),
            'm_bias' = median(bias),
            'bias' = mean(bias),
            'pct_m_bias' = median(pct_bias),
            'pct_bias' = mean(pct_bias)) %>% 
  select(sim_id, study_design, missingness, outcome_model, ipwr_purpose,  ipwr_model, ipwc_model,
         true_ate, std_dev, estimated_ate, bias, pct_bias, everything()) %>% 
  group_by(study_design, missingness) %>% 
  mutate('rel_eff' = std_dev/std_dev[1]) %>% 
  select(-std_dev) %>% 
  ungroup() 


for(id in unique(df_analysis$sim_id)) {
  subtitle <- glue('**Simulation**: {id}<br>
                   **Study Design**: {df_settings$study_design[id]}<br>
                   **Missingness**: {df_settings$missingness[id]}<br>
                   **True ATE**: {sprintf("%0.3f", df_truth$true_ate[df_truth$sim_id == id])}') 
  
  df_tmp <- 
    df_analysis %>% 
    filter(sim_id == id) %>% 
    select(-sim_id) %>% 
    mutate('outcome_model' = gsub('\\}', '</sup>', gsub('\\^\\{', '<sup>', outcome_model)),
           'ipwr_model' = gsub('\\}', '</sup>', gsub('\\^\\{', '<sup>', ipwr_model)),
           'ipwc_model' = gsub('\\}', '</sup>', gsub('\\^\\{', '<sup>', ipwc_model))
    ) %>% 
    select(-study_design, -missingness, -true_ate) %>% 
    mutate_if(is.character, ~gsub('---', '--', .x))
  
  ### Single Table
  if(id %in% c(1,2,3,4,8,9,10, 11,12)) {
    
    gt_results <-
      gt(df_tmp) %>% 
      fmt_markdown(columns = c(contains('model'), contains('ipw'))) %>% 
      cols_align('center') %>% 
      fmt_number(c(contains('ate'), 'bias', 'm_bias', 'rel_eff'), decimals = 3) %>% 
      fmt_percent(c('pct_bias', 'pct_m_bias'), decimals = 1) %>% 
      tab_spanner(columns = c('outcome_model', 'ipwc_model', 'ipwr_model', 'ipwr_purpose'), label = 'Analysis Mechanism') %>% 
      tab_spanner(c('estimated_ate', 'bias', 'pct_bias'), label = 'Mean') %>% 
      tab_spanner(c('m_estimated_ate', 'm_bias', 'pct_m_bias'), label = 'Median') %>% 
      tab_style(style = list(cell_borders(sides = "right", color = "black", weight = px(3))),
                locations = list(cells_body(columns = c('ipwr_purpose')))) %>% 
      tab_options(column_labels.font.size = 16,
                  heading.title.font.size = 30,
                  heading.subtitle.font.size = 20,
                  heading.title.font.weight = 'bold',
                  heading.subtitle.font.weight = 'bold',
                  column_labels.font.weight = 'bold',
                  row_group.font.weight = 'bold',
                  row_group.font.size  = 20) %>% 
      cols_label('outcome_model' = 'Outcome Model',
                 'ipwc_model' = 'IPWC Model',
                 'ipwr_model' = 'IPWR Model',
                 'ipwr_purpose' = 'IPWR Usage',
                 'estimated_ate' = 'Estimated ATE',
                 'm_estimated_ate' = 'Estimated ATE',
                 'bias' = 'Bias',
                 'pct_bias' = '% Bias',
                 'm_bias' = 'Bias',
                 'pct_m_bias' = '% Bias',
                 'rel_eff' = 'Relative Efficiency') %>% 
      tab_header(title = md('**Simulation Summary**'),
                 subtitle = md(subtitle)) %>% 
      tab_source_note(md('**IPWR Usage**:')) %>% 
      tab_source_note(md('**Estimation**: Weights used in estimation of analysis outcome model')) %>% 
      tab_source_note(md('**Marginalization**: Weights used in marginalization of estimated individual 
                       treatment effects over baseline covariate distribution'))
    
    
    if(dir.exists('logs/') ) {
      print(gt_results)
    }
    else {
      id <- ifelse(id < 10, paste0('0', id), id)
      gtsave(gt_results, glue('figures/simulations/missing_data/results/tables/gt_results_{id}.png'))
    }
    
    
  } else {
    if(TRUE) {
      gt_results <- 
        df_tmp %>% 
        group_by('stratified' = grepl('Stratified by A', ipwr_model)) %>% 
        group_split() %>% 
        map(~{
          subtitle_ <- paste0(subtitle, '<br>**IPWR Weight Stratification: ', .x$stratified[1], '**')
          
          .x %>% 
            select(-stratified) %>% 
            gt() %>% 
            fmt_markdown(columns = c(contains('model'), contains('ipw'))) %>% 
            cols_align('center') %>% 
            fmt_number(c(contains('ate'), 'bias', 'm_bias', 'rel_eff'), decimals = 3) %>% 
            fmt_percent(c('pct_bias', 'pct_m_bias'), decimals = 1) %>% 
            tab_spanner(columns = c('outcome_model', 'ipwc_model', 'ipwr_model', 'ipwr_purpose'), label = 'Analysis Mechanism') %>% 
            tab_spanner(c('estimated_ate', 'bias', 'pct_bias'), label = 'Mean') %>% 
            tab_spanner(c('m_estimated_ate', 'm_bias', 'pct_m_bias'), label = 'Median') %>% 
            tab_style(style = list(cell_borders(sides = "right", color = "black", weight = px(3))),
                      locations = list(cells_body(columns = c('ipwr_purpose')))) %>% 
            tab_options(column_labels.font.size = 16,
                        heading.title.font.size = 30,
                        heading.subtitle.font.size = 20,
                        heading.title.font.weight = 'bold',
                        heading.subtitle.font.weight = 'bold',
                        column_labels.font.weight = 'bold',
                        row_group.font.weight = 'bold',
                        row_group.font.size  = 20) %>% 
            cols_label('outcome_model' = 'Outcome Model',
                       'ipwc_model' = 'IPWC Model',
                       'ipwr_model' = 'IPWR Model',
                       'ipwr_purpose' = 'IPWR Usage',
                       'estimated_ate' = 'Estimated ATE',
                       'm_estimated_ate' = 'Estimated ATE',
                       'bias' = 'Bias',
                       'pct_bias' = '% Bias',
                       'm_bias' = 'Bias',
                       'pct_m_bias' = '% Bias',
                       'rel_eff' = 'Relative Efficiency') %>% 
            tab_header(title = md('**Simulation Summary**'),
                       subtitle = md(subtitle)) %>% 
            tab_source_note(md('**IPWR Usage**:')) %>% 
            tab_source_note(md('**Estimation**: Weights used in estimation of analysis outcome model')) %>% 
            tab_source_note(md('**Marginalization**: Weights used in marginalization of estimated individual 
                       treatment effects over baseline covariate distribution'))
          
          
        })
    } else {
      gt_results <- 
        df_tmp %>% 
        group_by('marginalized' = grepl('Margin', ipwr_purpose)) %>% 
        group_split() %>% 
        map(~{
          .x %>% 
            select(-marginalized) %>% 
            gt() %>% 
            fmt_markdown(columns = c(contains('model'), contains('ipw'))) %>% 
            cols_align('center') %>% 
            fmt_number(c(contains('ate'), 'bias', 'm_bias', 'rel_eff'), decimals = 3) %>% 
            fmt_percent(c('pct_bias', 'pct_m_bias'), decimals = 1) %>% 
            tab_spanner(columns = c('outcome_model', 'ipwc_model', 'ipwr_model', 'ipwr_purpose'), label = 'Analysis Mechanism') %>% 
            tab_spanner(c('estimated_ate', 'bias', 'pct_bias'), label = 'Mean') %>% 
            tab_spanner(c('m_estimated_ate', 'm_bias', 'pct_m_bias'), label = 'Median') %>% 
            tab_style(style = list(cell_borders(sides = "right", color = "black", weight = px(3))),
                      locations = list(cells_body(columns = c('ipwr_purpose')))) %>% 
            tab_options(column_labels.font.size = 16,
                        heading.title.font.size = 30,
                        heading.subtitle.font.size = 20,
                        heading.title.font.weight = 'bold',
                        heading.subtitle.font.weight = 'bold',
                        column_labels.font.weight = 'bold',
                        row_group.font.weight = 'bold',
                        row_group.font.size  = 20) %>% 
            cols_label('outcome_model' = 'Outcome Model',
                       'ipwc_model' = 'IPWC Model',
                       'ipwr_model' = 'IPWR Model',
                       'ipwr_purpose' = 'IPWR Usage',
                       'estimated_ate' = 'Estimated ATE',
                       'm_estimated_ate' = 'Estimated ATE',
                       'bias' = 'Bias',
                       'pct_bias' = '% Bias',
                       'm_bias' = 'Bias',
                       'pct_m_bias' = '% Bias',
                       'rel_eff' = 'Relative Efficiency') %>% 
            tab_header(title = md('**Simulation Summary**'),
                       subtitle = md(subtitle)) %>% 
            tab_source_note(md('**IPWR Usage**:')) %>% 
            tab_source_note(md('**Estimation**: Weights used in estimation of analysis outcome model')) %>% 
            tab_source_note(md('**Marginalization**: Weights used in marginalization of estimated individual 
                       treatment effects over baseline covariate distribution'))
        })
    }
    
    ### Save figure
    if(dir.exists('logs/')) {
      print(gt_results)
    } else {
      id <- ifelse(id < 10, paste0('0', id), id)
      gt_two_column_layout(gt_results, 
                           output = 'save',
                           vwidth = 2400,
                           filename = glue('figures/simulations/missing_data/results/tables/gt_results_{id}.png'))
    }
  }
}

df_results <- 
  df_results %>% 
  mutate('analysis_method' = paste(outcome_model, ipwc_model,  ipwr_model, ipwr_purpose, sep = ' | '))

ggplot(df_results, aes(x = analysis_method, y = bias)) + 
  coord_flip() +
  facet_grid(study_design ~ missingness, scales = 'free_y') + 
  geom_hline(data = df_truth, aes(yintercept = 0), lty = 2, lwd = 1) + 
  geom_boxplot(aes(fill = analysis_method, col = analysis_method), alpha = 0.2) + 
  theme(legend.position = 'none',
        axis.text.y = element_text(size = 6)) + 
  labs(x = 'Analysis Method\n(Outcome Model | IPWC Model | IPWR Model | IPW Usage',
       y = 'Bias',
       title = 'Distribution of Average Treatment Effect Estimates',
       subtitle = 'Various Missingness Mechanisms and Study Designs')
ggsave('figures/simulations/missing_data/results/bias_boxplots.png', height = 9/1.2, width = 16/1.2)

ggplot(df_results, aes(x = analysis_method, y = pct_bias)) + 
  facet_grid(study_design ~ missingness, scales = 'free_y') + 
  geom_hline(data = df_truth, aes(yintercept = 0), lty = 2, lwd = 1) + 
  geom_boxplot(aes(fill = analysis_method, col = analysis_method), alpha = 0.2) + 
  coord_flip() + 
  theme(legend.position = 'none',
        axis.text.y = element_text(size = 6)) + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = 'Analysis Method\n(Outcome Model | IPWC Model | IPWR Model | IPW Usage',
       y = 'Relative Bias (%)',
       title = 'Distribution of Average Treatment Effect Estimates',
       subtitle = 'Various Missingness Mechanisms and Study Designs')
ggsave('figures/simulations/missing_data/results/pct_bias_boxplots.png', height = 9/1.2, width = 16/1.2)

ggplot(df_results, aes(x = bias, y = analysis_method)) + 
  facet_grid(study_design ~ missingness, scales = 'free_x') + 
  geom_density_ridges(aes(fill = analysis_method), 
                      scale = 1,
                      quantiles = 0.5, 
                      quantile_lines = T,
                      alpha = 0.2) +
  geom_vline(data = df_truth, aes(xintercept = 0), lty = 2, lwd = 1) +
  theme(legend.position = 'none',
        axis.text.y = element_text(size = 6)) + 
  labs(x = 'Bias',
       y = 'Analysis Method\n(Outcome Model | IPWC Model | IPWR Model | IPW Usage',
       title = 'Distribution of Average Treatment Effect Estimates',
       subtitle = 'Various Missingness Mechanisms and Study Designs')
ggsave('figures/simulations/missing_data/results/bias_ridgeplots.png', height = 9/1.2, width = 16/1.2)

ggplot(df_results, aes(x = pct_bias, y = analysis_method)) + 
  facet_grid(study_design ~ missingness) + 
  geom_density_ridges(aes(fill = analysis_method), 
                      scale = 1,
                      quantiles = 0.5, 
                      quantile_lines = T,
                      alpha = 0.2) +
  geom_vline(data = df_truth, aes(xintercept = 0), lty = 2, lwd = 1) +
  theme(legend.position = 'none',
        axis.text.y = element_text(size = 6)) + 
  scale_x_continuous(labels = scales::percent) + 
  labs(x = 'Relative Bias (%)',
       y = 'Analysis Method\n(Outcome Model | IPWC Model | IPWR Model | IPW Usage',
       title = 'Distribution of Average Treatment Effect Estimates',
       subtitle = 'Various Missingness Mechanisms and Study Designs')
ggsave('figures/simulations/missing_data/results/pct_bias_ridgeplots.png', height = 9/1.2, width = 16/1.2)

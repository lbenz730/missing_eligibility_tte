library(tidyverse)
library(glue)
library(arrow)
library(knitr)
library(kableExtra)


### Read In Results
df_settings <- read_csv('data/simulations/missing_data/inputs/settings.csv')
files <- dir(glue('data/simulations/missing_data/outputs'), full.names = T)
df_results <- 
  map_dfr(files, read_parquet) %>% 
  mutate('ipwn_model' = gsub('C', 'N', ipwc_model)) %>% 
  mutate('bias' = estimated_ate - true_ate, 
         'pct_bias' = (estimated_ate - true_ate)/true_ate) 

df_truth <- 
  df_results %>% 
  distinct(sim_id, study_design, missingness, true_ate)

df_analysis <- 
  df_results %>% 
  group_by(sim_id, missingness, ipwr_model, ipwn_model) %>% 
  summarise('true_ate' = first(true_ate),
            'm_bias' = median(estimated_ate) - true_ate,
            'bias' = mean(bias),
            'pct_m_bias' = first((median(estimated_ate) - true_ate)/true_ate),
            'pct_bias' = mean(pct_bias)) %>% 
  select(sim_id, missingness, ipwr_model, ipwn_model, true_ate, bias, pct_bias, everything()) %>% 
  ungroup()


prepare_results <- function(sim_ids, stratify) {
  df <- 
    df_analysis %>% 
    filter(sim_id %in% sim_ids) %>% 
    mutate('stratified' = ifelse(grepl('Stratified', ipwr_model), 'Stratified by $A$', 'Unstratified'))  %>% 
    arrange(missingness, true_ate, desc(stratified)) %>% 
    select(missingness, true_ate, stratified, ipwr_model, ipwn_model,  bias, pct_bias, m_bias, pct_m_bias) %>% 
    mutate('missingness' = case_when(missingness == '(1) M-Bias' ~ 'M-Bias',
                                     missingness == '(2) Treatment Heterogeneity' ~ '\\shortstack{Treatment\\\\Effect\\\\Heterogeneity}',
                                     missingness == '(3) M-Bias w/ Mediator' ~ '\\shortstack{M-Bias\\\\w/ Mediator}')) %>% 
    
    mutate('pct_bias' = 100 * pct_bias,
           'pct_m_bias' = 100 * pct_m_bias, 
           'ipwr_model' = paste0('$', ipwr_model, '$'),
           'ipwr_model' = gsub(',', '', ipwr_model),
           'ipwn_model' = paste0('$', ipwn_model, '$')) %>% 
    mutate('ipwr_model' = gsub('\\$---\\$', '---', ipwr_model),
           'ipwr_model' = gsub('Stratified by A', '', ipwr_model),
           'ipwn_model' = gsub('\\$---\\$', '---', ipwn_model)) %>% 
    mutate('ipwr_model' = gsub('~', '\\\\sim', gsub('L', '\\\\bm{L}', ipwr_model)),
           'ipwn_model' = gsub('~', '\\\\sim', gsub('L', '\\\\bm{L}', ipwn_model)))
  
  if(!stratify) {
    df <- 
      df %>% 
      select(-stratified)
    
    
    names(df) <-
      c('\\textbf{Missingness}', '$\\bm\\psi_{PP}$', '$\\bm W_{mk}^{R}$', '$\\bm W_{mkt}^{N}$',
        '\\textbf{Bias}', '\\textbf{\\% Bias}', '\\textbf{Bias}', '\\textbf{\\% Bias}')
  } else {
    names(df) <-
      c('\\textbf{Missingness}', '$\\bm\\psi_{PP}$', '\\textbf{Statification}', '$\\bm W_{mk}^{R}$', '$\\bm W_{mkt}^{N}$',
        '\\textbf{Bias}', '\\textbf{\\% Bias}', '\\textbf{Bias}', '\\textbf{\\% Bias}')
  }
  
  
  return(df)
  
}

write_latex <- function(sim_ids, caption, label, stratify) {
  df_tbl <- prepare_results(sim_ids, stratify) 
  if(stratify) {
    dig <- c(0,3,0, 0, 0, 3, 1, 3, 1) 
    col <- c(1:4)
  } else {
    dig <-  c(0,3,0, 0, 3, 1, 3, 1)
    col <- c(1:3)
  }
  latex <- 
    df_tbl %>% 
    kbl(align = 'c', digits = dig , escape = F, format = 'latex') %>%
    column_spec(1, border_left = T)  %>%
    column_spec(ncol(df_tbl), border_right = T) %>% 
    collapse_rows(columns = col, valign = 'middle') %>% 
    add_header_above(c('\\\\textbf{Setting}' = 2, '\\\\textbf{IPW Models}' = 2 + stratify, '\\\\textbf{Mean} $\\\\bm{ \\\\hat \\\\psi}_{PP}$' = 2, '\\\\textbf{Median} $\\\\bm{ \\\\hat{ \\\\psi}_{PP}}$' = 2), escape = F) %>%
    gsub('tabular\\}\\[t]', 'tabular\\}', .) %>% 
    gsub('c\\|c\\|c\\|c\\|c\\|c', 'Sc\\|Sc\\|Sc\\|Sc\\|Sc\\|Sc', .) %>%
    gsub('\\multicolumn\\{2\\}\\{c\\|\\}\\{\\\\textbf\\{Setting\\}\\}', '\\multicolumn\\{2\\}\\{\\|c\\|\\}\\{\\\\textbf\\{Setting\\}\\}', .) %>% 
    gsub('\\multicolumn\\{2\\}\\{c\\}\\{\\\\textbf\\{Median\\} \\$\\\\bm\\{ \\\\hat\\{ \\\\psi\\}_\\{PP\\}\\}\\$\\}', 
         '\\multicolumn\\{2\\}\\{Sc\\|\\}\\{\\\\textbf\\{Median\\} \\$\\\\bm\\{ \\\\hat\\{ \\\\psi\\}_\\{PP\\}\\}\\$\\}', .) %>%
    paste0("\\begin{table}[H]\n\\scriptsize\n\\centering", ., label, '\n\\end{table}')
  
  
  latex <- unlist(strsplit(latex, '\\\\end\\{tabular\\}'))
  latex <- paste0(latex[1], paste0('\\end{tabular}', caption), latex[2])
  
  return(latex)
  
}

ltx_1 <- 
  write_latex(sim_ids = c(9,10,11),
              stratify = F,
              label = '\\label{table:results_1}',
              caption = '\\caption{Simulation results from hypothetical study \\#1. Eligibility criteria: BMI $\\geq$ 35 m/kg$^2$, A1c $\\geq$ 5.7 \\%, no previous initiation of bariatric surgery.}')


write(ltx_1, 'figures/paper1/tables/sim_table_1.tex')


ltx_2 <- 
  write_latex(sim_ids = c(5,6,7),
              stratify = T,
              label = '\\label{table:results_2}',
              caption = '\\caption{Simulation results from hypothetical study \\#2. Eligibility criteria: BMI $\\geq$ 35 m/kg$^2$, no previous initiation of bariatric surgery.}')


write(ltx_2, 'figures/paper1/tables/sim_table_2.tex')







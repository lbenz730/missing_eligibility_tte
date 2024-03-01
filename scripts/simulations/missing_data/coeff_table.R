library(tidyverse)
library(glue)
library(arrow)
library(knitr)
library(kableExtra)
library(xtable)


params <- read_rds('data/simulations/missing_data/inputs/sim_params_5.rds')
params_mediator <- read_rds('data/simulations/missing_data/inputs/sim_params_7.rds')


### Weight/A1c Models for Surgery
weight_surg <- 
  tibble('measure' = 'BMI',
         'coeff_1' = params$weight_models$surgery$spline$coeff$`1`,
         'coeff_2' = params$weight_models$surgery$spline$coeff$`2`,
         'coeff_3' = params$weight_models$surgery$spline$coeff$`3`,
         'coeff_4' = params$weight_models$surgery$spline$coeff$`4`,
         'coeff_5' = params$weight_models$surgery$spline$coeff$`5`)

a1c_surg <- 
  tibble('measure' = 'A1c',
         'coeff_1' = params$a1c_models$surgery$spline$coeff$`1`,
         'coeff_2' = params$a1c_models$surgery$spline$coeff$`2`,
         'coeff_3' = params$a1c_models$surgery$spline$coeff$`3`,
         'coeff_4' = params$a1c_models$surgery$spline$coeff$`4`,
         'coeff_5' = params$a1c_models$surgery$spline$coeff$`5`)

df_surg <- bind_rows(weight_surg, a1c_surg)

names(df_surg) <- c('\\textbf{Measure}', glue('$\\bm\\beta_{1:5}$'))
caption <- 
  paste('\\caption{Spline coefficients $\\bm\\beta_j$ for post-surgical BMI and A1c trajectories.', 
        'A subject\'s set of coefficients are determined by latent post-surgical classes $\\eta_k$ (for BMI) and $\\xi_{k}$ (for A1c).',
        'Knots for BMI splines are at $t = 3, 6, 9, 12, 15, 18, 21, 24, 36$, and $48$ months, while knots for A1c splines are at',
        '$3, 9, 15, 30$ and $48$ months.}')

label <- '\\label{table:spline_coef}'

latex <- 
  df_surg %>% 
  kbl(align = 'c', digits = 3 , escape = F, format = 'latex') %>% 
  column_spec(1, border_left = T)  %>%
  column_spec(ncol(weight_surg), border_right = T) %>% 
  collapse_rows(columns = 1, valign = 'middle') %>% 
  paste0("\\begin{table}[H]\n\\scriptsize\n\\centering", ., '\n\\end{table}') %>% 
  gsub('tabular\\}\\[t]', 'tabular\\}', .) %>% 
  gsub('c\\|c\\|c\\|', 'Sc\\|Sc\\|Sc\\|', .) %>% 
  paste0(., label) 


latex <- unlist(strsplit(latex, '\\\\end\\{tabular\\}'))
latex <- paste0(latex[1], paste0('\\end{tabular}', caption), latex[2])

write(latex, 'figures/paper1/tables/post_surgery_spline.tex')


### Model Coeff
params5 <- read_rds('data/simulations/missing_data/inputs/sim_params_5.rds')
params6 <- read_rds('data/simulations/missing_data/inputs/sim_params_6.rds')
params7 <- read_rds('data/simulations/missing_data/inputs/sim_params_7.rds')
params9 <- read_rds('data/simulations/missing_data/inputs/sim_params_9.rds')
params10 <- read_rds('data/simulations/missing_data/inputs/sim_params_10.rds')
params11 <- read_rds('data/simulations/missing_data/inputs/sim_params_11.rds')


gather_coeff <- function(params) {
  bind_rows(
    tibble('model' = 'Pre-Surgery BMI',
           'model_num' = 1,
           'term_num' = 1,
           'coeff_greek' = '$\\bm \\beta_1$',
           'term' = names(params$weight_models$control$slope),
           'coeff_value' = unlist(params$weight_models$control$slope)),
    
    tibble('model' = 'Pre-Surgery BMI',
           'model_num' = 1,
           'term_num' = 2,
           'coeff_greek' = '$\\bm \\beta_2$',
           'term' = names(params$weight_models$control$quadratic$term),
           'coeff_value' = unlist(params$weight_models$control$quadratic$term)),
    
    tibble('model' = 'Pre-Surgery BMI',
           'model_num' = 1,
           'term_num' = 3,
           'coeff_greek' = '$\\bm \\delta$',
           'term' = names(params$weight_models$control$quadratic$prob_nonzero),
           'coeff_value' = unlist(params$weight_models$control$quadratic$prob_nonzero)),
    
    tibble('model' = 'Pre-Surgery BMI',
           'model_num' = 1,
           'term_num' = 4,
           'coeff_greek' = '$\\sigma^2_{\\text{BMI}}$',
           'term' = '---',
           'coeff_value' = params$weight_models$control$errors$std_dev),
    
    tibble('model' = 'Pre-Surgery BMI',
           'model_num' = 1,
           'term_num' = 5,
           'coeff_greek' = '$\\tau^2_1$',
           'term' = '---',
           'coeff_value' = params$weight_models$control$std_dev$slope),
    
    tibble('model' = 'Pre-Surgery BMI',
           'model_num' = 1,
           'term_num' = 6,
           'coeff_greek' = '$\\tau^2_2$',
           'term' = '---',
           'coeff_value' = params$weight_models$control$std_dev$quadratic),
    
    tibble('model' = 'Pre-Surgery A1c',
           'model_num' = 2,
           'term_num' = 1,
           'coeff_greek' = '$\\bm \\beta_3$',
           'term' = names(params$a1c_models$control$slope),
           'coeff_value' = unlist(params$a1c_models$control$slope)),
    
    tibble('model' = 'Pre-Surgery A1c',
           'model_num' = 2,
           'term_num' = 2,
           'coeff_greek' = '$\\sigma^2_{\\text{A1c}}$',
           'term' = 'I(hgba1c\\string^2)',
           'coeff_value' = 0.001),
    
    tibble('model' = 'Pre-Surgery A1c',
           'model_num' = 2,
           'term_num' = 3,
           'coeff_greek' = '$\\tau^3_2$',
           'term' = '---',
           'coeff_value' = params$a1c_models$control$std_dev),
    
    tibble('model' = 'Treatment',
           'model_num' = 3,
           'term_num' = 1,
           'coeff_greek' = '$\\bm\\alpha$',
           'term' = names(params$treatment_model),
           'coeff_value' = unlist(params$treatment_model)),
    
    tibble('model' = 'Treatment',
           'model_num' = 3,
           'term_num' = 2,
           'coeff_greek' = '$\\bm\\pi$',
           'term' = names(params$bs_type_model),
           'coeff_value' = unlist(params$bs_type_model)),
    
    tibble('model' = 'Post-Surgery BMI', 
           'model_num' = 4,
           'term_num' = 1,
           'coeff_greek' = paste0('$\\lambda_{0', 1:4, '}$'),
           'term' = '(Intercept)', 
           'coeff_value' = params$weight_models$surgery$group_model$intercepts),
    
    tibble('model' = 'Post-Surgery BMI', 
           'model_num' = 4,
           'term_num' = 2,
           'coeff_greek' = '$\\bm\\lambda_1$',
           'term' = names(params$weight_models$surgery$group_model[setdiff(names(params$weight_models$surgery$group_model), 'intercepts')]),
           'coeff_value' = unlist(params$weight_models$surgery$group_model[setdiff(names(params$weight_models$surgery$group_model), 'intercepts')])),
    
    tibble('model' = 'Post-Surgery BMI',
           'model_num' = 4,
           'term_num' = 3,
           'coeff_greek' = '$\\tau^2_4$',
           'term' = '---',
           'coeff_value' = params$weight_models$surgery$std_dev),
    
    tibble('model' = 'Post-Surgery A1c', 
           'model_num' = 5,
           'term_num' = 1,
           'coeff_greek' = paste0('$\\phi_{0', 1:4, '}$'),
           'term' = '(Intercept)', 
           'coeff_value' = params$a1c_models$surgery$group_model$intercepts),
    
    tibble('model' = 'Post-Surgery A1c', 
           'model_num' = 5,
           'term_num' = 2,
           'coeff_greek' = '$\\bm\\phi_1$',
           'term' = names(params$a1c_models$surgery$group_model[setdiff(names(params$a1c_models$surgery$group_model), 'intercepts')]),
           'coeff_value' = unlist(params$a1c_models$surgery$group_model[setdiff(names(params$a1c_models$surgery$group_model), 'intercepts')])),
    
    tibble('model' = 'Post-Surgery A1c',
           'model_num' = 5,
           'term_num' = 3,
           'coeff_greek' = '$\\tau^2_5$',
           'term' = '---',
           'coeff_value' = params$a1c_models$surgery$std_dev),
    
    tibble('model' = 'Outcome',
           'model_num' = 6,
           'term_num' = 1,
           'coeff_greek' = '$\\bm\\omega$',
           'term' = names(params$outcome_model),
           'coeff_value' = unlist(params$outcome)),
    
    tibble('model' = 'Missing Eligibility',
           'model_num' = 7,
           'term_num' = 1,
           'coeff_greek' = '$\\bm\\rho$',
           'term' = names(params$cc_model),
           'coeff_value' = unlist(params$cc_model))
    
  ) %>% 
    mutate('term' = gsub('\\_', '\\\\_', term))
  
  
}

sci_notation <- function(x) {
  case_when(is.na(x) ~ '0',
            x == 0 ~ '0',
            abs(x) > 0.1 ~ sprintf('%0.2f', x),
            T ~ sanitize.numbers(format(x, 
                                        digits = 2,
                                        scientific = T),
                                 type = 'latex',
                                 math.style.exponents = T))
  
}



df_tbl <- 
  bind_rows(
    gather_coeff(params5) %>% mutate('params' = '5'),
    gather_coeff(params6) %>% mutate('params' = '6'),
    gather_coeff(params7) %>% mutate('params' = '7')
  ) %>% 
  mutate('term' = paste0('\\texttt{', term, '}')) %>% 
  mutate('coeff_value' = map_chr(coeff_value, sci_notation)) %>% 
  pivot_wider(names_from = 'params',
              values_from = 'coeff_value',
              names_prefix = 'sim_') %>% 
  arrange(model_num, term_num) %>% 
  select(-model_num, -term_num) %>% 
  mutate_at(vars(contains('sim')), ~ifelse(is.na(.x), '0', .x))

names(df_tbl) <- c('\\textbf{Component Model}', '\\textbf{Coefficient}', '\\textbf{Variable}', 
                   '\\textbf{Fig. \\ref{fig:missing_DAGs_A}}', 
                   '\\textbf{Fig. \\ref{fig:missing_DAGs_B}}', 
                   '\\textbf{Fig. \\ref{fig:missing_DAGs_C}}')

caption <- 
  paste('\\caption{Coefficients for simulation settings, summarized for each missing data mechanism in Figure \\ref{fig:missing_DAGs}.',
        'Component variables comprising $\\bm L_{k0}$ include \\texttt{gender}, \\texttt{race}, \\texttt{insulin} (usage), \\texttt{smoking\\_status}',
        'and \\texttt{elix\\_score} (comorbidities). Component variables comprising $\\bm L_{kt}$ include \\texttt{bmi}, \\texttt{hgba1c}. Treatment variables',
        'are denoted by \\texttt{surgery} ($A_{kt}$) and \\texttt{bs\\_type} ($A\'_{kt}$)}')

label <- '\\label{table:sim_coef}'

latex <- 
  df_tbl %>% 
  kbl(align = 'c', digits = 3 , escape = F, format = 'latex') %>% 
  column_spec(1, border_left = T)  %>%
  column_spec(ncol(weight_surg), border_right = T) %>% 
  collapse_rows(columns = c(1, 2), valign = 'middle') %>% 
  add_header_above(c('\\\\textbf{Model Information}' = 3, '\\\\textbf{Missingness Mechanism}' = 3), escape = F) %>%
  paste0("\\begin{table}[H]\n\\tiny\n\\centering", ., '\n\\end{table}') %>% 
  gsub('\\\\multicolumn\\{3\\}\\{c\\}\\{\\\\textbf\\{Missingness Mechanism\\}\\}',
       '\\\\multicolumn\\{3\\}\\{Sc\\|\\}\\{\\\\textbf\\{Missingness Mechanism\\}\\}', .) %>% 
  gsub('\\\\multicolumn\\{3\\}\\{c\\|\\}\\{\\\\textbf\\{Model Information\\}\\}',
       '\\\\multicolumn\\{3\\}\\{\\|Sc\\|\\}\\{\\\\textbf\\{Model Information\\}\\}', .) %>% 
  gsub('tabular\\}\\[t]', 'tabular\\}', .) %>% 
  gsub('c\\|c\\|c\\|', 'Sc\\|Sc\\|Sc\\|', .) %>% 
  paste0(., label) 


latex <- unlist(strsplit(latex, '\\\\end\\{tabular\\}'))
latex <- paste0(latex[1], paste0('\\end{tabular}', caption), latex[2])
write(latex, 'figures/paper1/tables/sim_coeff.tex')





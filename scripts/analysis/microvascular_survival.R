library(tidyverse)
library(lubridate)
library(arrow)
library(glue)
library(survival)
library(survminer)
source('scripts/util/helpers.R')

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Plot Survival of cases for microvascular outcomes
case_population <- read_parquet(glue('{data_dir}/microvascular_cases.parquet'))

df_events <- 
  case_population %>% 
  mutate('censor_time' = pmin(censor_enrollment, 
                              censor_death, 
                              censor_cancer, 
                              censor_measurement,
                              end_of_study, 
                              na.rm = T),
         'outcome_time' = pmin(neuropathy,
                               nephropathy,
                               retinopathy,
                               na.rm = T)) %>% 
  mutate('time_combined' = as.numeric(pmin(outcome_time, censor_time, na.rm = T) - index_date),
         'time_neuropathy' = as.numeric(pmin(neuropathy, censor_time, na.rm = T) - index_date),
         'time_nephropathy' = as.numeric(pmin(neuropathy, censor_time, na.rm = T) - index_date),
         'time_retinopathy' = as.numeric(pmin(retinopathy, censor_time, na.rm = T) - index_date),
         'event_combined' = as.numeric(!is.na(outcome_time) & outcome_time < censor_time),
         'event_neuropathy' = as.numeric(!is.na(neuropathy) & neuropathy < censor_time),
         'event_nephropathy' = as.numeric(!is.na(nephropathy) & nephropathy < censor_time),
         'event_retinopathy' = as.numeric(!is.na(retinopathy) & retinopathy < censor_time)) 


## https://rpkgs.datanovia.com/survminer/reference/ggsurvplot_combine.html
surv_list <- 
  list('Microvascular Event' = survfit(Surv(time_combined, event_combined) ~ 1, data = df_events),
       'Neuropathy' = survfit(Surv(time_neuropathy, event_neuropathy) ~ 1, data = df_events),
       'Nephropathy' = survfit(Surv(time_nephropathy, event_nephropathy) ~ 1, data = df_events),
       'Retinopathy' = survfit(Surv(time_retinopathy, event_retinopathy) ~ 1, data = df_events))

survplot <- 
  ggsurvplot_combine(surv_list,
                     data = df_events,
                     
                     conf.int = F,
                     risk.table = T,
                     fun = 'event',
                     
                     break.x.by = 365,
                     xlim = c(0, 7)* 365.25,
                     ylim = c(0, 0.25),
                     
                     xlab = 'Time since Index Date (Years)',
                     title = 'Bariatric Surgery Group',
                     ylab = 'Cumulative Incidence',
                     subtitle = 'Cumulative Incidence of Microvascular Events',
                     legend.labs = c(names(surv_list)),
                     legend.title = '',
                     
                     
                     tables.theme = 
                       theme_minimal() + 
                       theme(plot.title = element_text(hjust = 0.5),
                             plot.subtitle = element_text(hjust = 0.5),
                             axis.title = element_text(size = 16),
                             axis.text = element_text(size = 12),
                             panel.background = element_blank(),
                             panel.grid = element_blank()
                             ),
                     
                     
                     
                     ggtheme = 
                       theme_minimal() + 
                       theme(plot.title = element_text(hjust = 0.5, size = 24),
                             plot.subtitle = element_text(hjust = 0.5, size = 20),
                             axis.title = element_text(size = 16),
                             axis.text = element_text(size = 12))
                     )

survplot$plot <- 
  survplot$plot + 
  theme(legend.text = element_text(size = 12)) + 
  scale_x_continuous(breaks = c(0:7)* 365.25, labels = function(x) x/365.25) + 
  scale_y_continuous(limits = c(0, 0.25), labels = scales::percent)

survplot$table <- 
  survplot$table + 
  scale_x_continuous(breaks = c(0:7)* 365.25, labels = function(x) x/365.25) 


ggsave(survplot, filename = 'figures/paper1/bariatric_microvascular.png', height = 9/1.2, width = 16/1.2)

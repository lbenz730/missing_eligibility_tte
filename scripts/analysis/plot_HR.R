library(tidyverse)
library(glue)
library(furrr)
library(patchwork)
library(knitr)
library(kableExtra)
plan(multisession(workers = 12))
source('scripts/util/helpers.R')

df_results <- future_map_dfr(dir('data/microvascular_results/', full.names = T), ~read_csv(.x, col_types = cols()))
df_results <- 
  df_results %>% 
  mutate('weights_exp' = gsub('x', '%*%', weights),
         'bmi_lookback_exp' = fct_reorder(paste('BMI Lookback:', bmi_lookback, ifelse(bmi_lookback == 1, 'Month', 'Months')), bmi_lookback),
         'weights_exp_f' = fct_relevel(weights_exp, 'W^A %*% W^N %*% W^R', after = 4))

### ITT
ggplot(df_results %>% filter(effect == 'Intent-To-Treat'), aes(x = a1c_lookback, y = exp_est)) + 
  geom_blank(data = df_results %>% filter(effect == 'Per-Protocol')) + 
  facet_grid(outcome ~ bmi_lookback_exp, scales = 'free_y') + 
  geom_line(aes(color = weights_exp)) +
  geom_point(aes(color = weights_exp), size = 4) + 
  scale_x_continuous(breaks = seq(0,24,6)) +
  scale_color_manual(values = gg_color_hue(5)[c(1,2,4)], labels = function(x) parse(text = x)) + 
  labs(x = 'Blood Glucose Labs Lookback Time (Months)',
       y = expression(paste('Estimated Discrete Hazard Ratio (', e^hat(psi)[ITT], ')')),
       color = 'Inverse Probability Weights',
       title = 'Effect of Bariatric Surgery on Microvascular Outcomes',
       subtitle = 'Sensitivity to Various Eligibility Lookback Times') + 
  theme(strip.text.y = element_text(size = 11),
        strip.text.x = element_text(size = 12),
        panel.spacing=unit(1,"lines"),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 32),
        plot.subtitle = element_text(size = 24))

ggsave('figures/paper1/bariatric_microvascular_sensitivity_itt.png', height = 12/1.2, width = 16/1.2)


### Per protocol
ggplot(df_results %>% filter(effect == 'Per-Protocol'), aes(x = a1c_lookback, y = exp_est)) + 
  geom_blank(data = df_results %>% filter(effect == 'Intent-To-Treat')) + 
  facet_grid(outcome ~ bmi_lookback_exp, scales = 'free_y') +
  geom_line(aes(color = weights_exp)) +
  geom_point(aes(color = weights_exp), size = 4) + 
  scale_x_continuous(breaks = seq(0,24,6)) +
  scale_color_manual(values = gg_color_hue(5)[c(1,2,3,5)], labels = function(x) parse(text = x)) + 
  labs(x = 'Blood Glucose Labs Lookback Time (Months)',
       y = expression(paste('Estimated Discrete Hazard Ratio (', e^hat(psi)[PP], ')')),
       color = 'Inverse Probability Weights',
       title = 'Effect of Bariatric Surgery on Microvascular Outcomes',
       subtitle = 'Sensitivity to Various Eligibility Lookback Times') + 
  theme(strip.text.y = element_text(size = 11),
        strip.text.x = element_text(size = 12),
        panel.spacing=unit(1,"lines"),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 32),
        plot.subtitle = element_text(size = 24))

ggsave('figures/paper1/bariatric_microvascular_sensitivity_pp.png', height = 12/1.2, width = 16/1.2)s


### Combined Figure
p1 <-
  ggplot(df_results %>% filter(effect == 'Intent-To-Treat'), aes(x = a1c_lookback, y = exp_est)) + 
  geom_blank(data = df_results %>% filter(effect == 'Per-Protocol')) + 
  facet_grid(outcome ~ bmi_lookback_exp, scales = 'free_y', labeller = label_wrap_gen(width = 15)) + 
  geom_line(aes(color = weights_exp_f)) +
  geom_point(aes(color = weights_exp_f), size = 4) + 
  scale_x_continuous(breaks = seq(0,24,6)) +
  scale_color_discrete(labels = function(x) parse(text = x), drop = F) + 
  labs(x = 'Blood Glucose Labs Lookback Time (Months)',
       y = expression(paste('Estimated Discrete Hazard Ratio (', e^hat(psi)[ITT], ')')),
       color = 'Inverse Probability Weights',
       title = 'Intention-To-Treat Effect') + 
  theme(strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 16),
        panel.spacing=unit(1,"lines"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 24),
        plot.title = element_text(size = 30))

### Per protocol
p2 <- 
  ggplot(df_results %>% filter(effect == 'Per-Protocol'), aes(x = a1c_lookback, y = exp_est)) + 
  geom_blank(data = df_results %>% filter(effect == 'Intent-To-Treat')) + 
  facet_grid(outcome ~ bmi_lookback_exp, scales = 'free_y', labeller = label_wrap_gen(width = 15)) + 
  geom_line(aes(color = weights_exp_f)) +
  geom_point(aes(color = weights_exp_f), size = 4) + 
  scale_x_continuous(breaks = seq(0,24,6)) +
  scale_color_discrete(labels = function(x) parse(text = x), drop = F) + 
  labs(x = 'Blood Glucose Labs Lookback Time (Months)',
       y = expression(paste('Estimated Discrete Hazard Ratio (', e^hat(psi)[PP], ')')),
       color = 'Inverse Probability Weights',
       title = 'Per-Protocol Effect'
  ) + 
  theme(strip.text.y = element_text(size = 16),
        strip.text.x = element_text(size = 16),
        panel.spacing=unit(1,"lines"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 24),
        plot.title = element_text(size = 30))



p1/p2 + 
  plot_layout(guides = 'collect') &
  plot_annotation(title = 'Effect of Bariatric Surgery on Microvascular Outcomes',
                  subtitle = 'Sensitivity to Various Eligibility Lookback Times',
                  theme = theme(plot.title = element_text(size = 36, hjust = 0.5),
                                plot.subtitle = element_text(size = 32, hjust = 0.5))) 


ggsave('figures/paper1/bariatric_microvascular_sensitivity_combined.png', height = 24/1.1, width = 16/1.1)
ggsave('figures/paper1/bariatric_microvascular_sensitivity_combined.pdf', height = 24/1.1, width = 16/1.1)

p1 + p2 + 
  plot_layout(guides = 'collect') &
  plot_annotation(title = 'Effect of Bariatric Surgery on Microvascular Outcomes',
                  subtitle = 'Sensitivity to Various Eligibility Lookback Times',
                  theme = theme(plot.title = element_text(size = 36, hjust = 0.5),
                                plot.subtitle = element_text(size = 32, hjust = 0.5))) 


ggsave('figures/paper1/bariatric_microvascular_sensitivity_combined_wide.png', height = 12/1.1, width = 32/1.1)



### Boostrap results for final (3,12) lookback analysis
df_itt <-  future_map_dfr(dir('data/final_itt_bootstraps/', full.names = T), ~read_csv(.x, col_types = cols()))
df_pp <- future_map_dfr(dir('data/final_pp_bootstraps/', full.names = T), ~read_csv(.x, col_types = cols()))
df_results %>% 
  filter(bmi_lookback == 3, a1c_lookback == 12, grepl('R', weights))  %>% 
  group_by(outcome, effect) %>% 
  summarise('hr' = exp(estimate)) %>% 
  ungroup() %>% 
  left_join(
    bind_rows(df_itt, df_pp) %>% 
      group_by(outcome, effect) %>% 
      summarise('lower' = exp(mean(estimate, na.rm = T) + qnorm(0.025) * sd(estimate, na.rm = T)),
                'upper' = exp(mean(estimate, na.rm = T) + qnorm(0.975) * sd(estimate, na.rm = T)))
  ) %>% 
  mutate('text' = paste0(sprintf('%0.3f', hr), ' (', sprintf('%0.3f', lower), ', ', sprintf('%0.3f', upper), ')')) %>% 
  select(outcome, effect, text) %>% 
  pivot_wider(names_from = 'effect', 
              values_from = 'text') %>% 
  kable(format = 'latex')



### ITT CIs
df_bootstraps <- future_map_dfr(dir('data/microvascular_bootstrap_results/', full.names = T), ~read_csv(.x, col_types = cols()))
df_boot <- 
  df_bootstraps %>% 
  group_by(outcome, weights, bmi_lookback, a1c_lookback, effect) %>% 
  summarise('n' = n(),
            'mean' = mean(exp_est),
            'var' = var(exp_est),
            'lower' = quantile(exp_est, 0.025),
            'upper' = quantile(exp_est, 0.975)) %>% 
  mutate('width' = upper - lower) %>% 
  ungroup() %>% 
  inner_join(df_results, by = c('outcome', 'weights', 'bmi_lookback', 'a1c_lookback', 'effect')) %>% 
  mutate('weights_exp' = gsub('x', '%*%', weights),
         'bmi_lookback_exp' = fct_reorder(paste('BMI Lookback:', bmi_lookback, ifelse(bmi_lookback == 1, 'Month', 'Months')), bmi_lookback), 
  )

ggplot(df_boot %>% filter(weights == 'W^A x W^R')  , aes(x = a1c_lookback, y = exp_est)) + 
  facet_grid(outcome ~ bmi_lookback_exp) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), color = 'deepskyblue', lwd = 1) +
  geom_line(lty = 2, color = 'navy') +
  geom_point(size = 3, color = 'navy') + 
  labs(x = 'Blood Glucose Labs Lookback Time (Months)',
       y = expression(paste('Estimated Discrete Hazard Ratio (', e^hat(psi)[ITT], ')')),
       title = expression(paste('Variance of IPW Estimator With Weights ',  W^A %*% W^R)),
       subtitle = '100 Bootstrap Replicates') + 
  scale_x_continuous(breaks = seq(0,24,6)) +
  theme(strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        panel.spacing=unit(1,"lines"),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 32),
        plot.subtitle = element_text(size = 24))
ggsave('figures/paper1/ipwr_ci.png', height = 12/1.2, width = 16/1.2)  


ggplot(df_boot , aes(x = a1c_lookback, y = width)) + 
  facet_grid(outcome ~ bmi_lookback_exp) + 
  geom_line(aes(color = weights_exp)) +
  geom_point(aes(color = weights_exp), size = 3) + 
  scale_color_manual(values = gg_color_hue(5)[c(1,2,4)], labels = function(x) parse(text = x)) + 
  scale_x_continuous(breaks = seq(0,24,6)) +
  labs(x = 'Blood Glucose Labs Lookback Time (Months)',
       y = 'Width of 95% Confidence Interval',
       color = 'Inverse Probability Weights',
       title = '95% Confidence Interval Width of IPW Estimators',
       subtitle = '100 Bootstrap Replicates') + 
  theme(strip.text.y = element_text(size = 11),
        strip.text.x = element_text(size = 12),
        panel.spacing=unit(1,"lines"),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 32),
        plot.subtitle = element_text(size = 24))

ggsave('figures/paper1/ci_widths.png', height = 12/1.2, width = 16/1.2)  

ggplot(df_boot , aes(x = a1c_lookback, y = var)) + 
  facet_grid(outcome ~ bmi_lookback_exp, scales = 'free_y') + 
  geom_line(aes(color = weights_exp)) +
  geom_point(aes(color = weights_exp), size = 3) + 
  scale_color_manual(values = gg_color_hue(5)[c(1,2,4)], labels = function(x) parse(text = x)) + 
  scale_x_continuous(breaks = seq(0,24,6)) +
  labs(x = 'Blood Glucose Labs Lookback Time (Months)',
       y = expression(paste('Var(', e^hat(psi)[ITT], ')')),
       color = 'Inverse Probability Weights',
       title = 'Variance of IPW Estimators',
       subtitle = '100 Bootstrap Replicates') + 
  theme(strip.text.y = element_text(size = 11),
        strip.text.x = element_text(size = 12),
        panel.spacing=unit(1,"lines"),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 32),
        plot.subtitle = element_text(size = 24))

ggsave('figures/paper1/boot_var.png', height = 12/1.2, width = 16/1.2)  


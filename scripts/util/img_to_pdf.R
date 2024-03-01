library(tidyverse)
library(magick)
### Convert images into one PDF of plots
img_to_pdf <- function(fig_dir, outfile) {
  img_files <-  dir(fig_dir, full.names = T)
  img_files <- img_files[grepl('\\.png', img_files)]
  img_pdf <- reduce(map(img_files, image_read), c)
  image_write(img_pdf, format = 'PDF', paste0(fig_dir, '/', outfile))
}

# img_to_pdf(fig_dir = 'figures/eda', outfile = 'tte.pdf')
# img_to_pdf(fig_dir = 'figures/missing_data', outfile = 'plots.pdf')
# img_to_pdf(fig_dir = 'figures/simulations/complete_case_sims/eda/', outfile = 'simulation_settings.pdf')
# img_to_pdf(fig_dir = 'figures/simulations/robust_variance_sims/v1/eda/', outfile = 'simulation_settings.pdf')
# img_to_pdf(fig_dir = 'figures/simulations/robust_variance_sims/v1/bias/', outfile = 'simulation_results_bias.pdf')
# img_to_pdf(fig_dir = 'figures/simulations/robust_variance_sims/v1/coverage/', outfile = 'simulation_results_coverage.pdf')
# img_to_pdf(fig_dir = 'figures/simulations/robust_variance_sims/v2/eda/', outfile = 'simulation_settings.pdf')
# img_to_pdf(fig_dir = 'figures/simulations/robust_variance_sims/v2/bias/', outfile = 'simulation_results_bias.pdf')
# img_to_pdf(fig_dir = 'figures/simulations/robust_variance_sims/v2/coverage/', outfile = 'simulation_results_coverage.pdf')
# img_to_pdf(fig_dir = 'figures/simulations/robust_variance_sims/v2/incidence/', outfile = 'event_rates.pdf')
# img_to_pdf(fig_dir = 'figures/simulations/robust_variance_sims/v2/correlation/', outfile = 'correlation_eda.pdf')
# img_to_pdf(fig_dir = 'figures/simulations/missing_data/eda/', outfile = 'weight_trajectory_plots.pdf')
# img_to_pdf(fig_dir = 'figures/simulations/missing_data/results/', outfile = 'bias_plots.pdf')
# img_to_pdf(fig_dir = 'figures/simulations/missing_data/results/tables', outfile = 'result_tables.pdf')
img_to_pdf(fig_dir = 'figures/paper1/', outfile = 'paper1_figs.pdf')


# Missing Eligibity Data in Target Trial Emulation
Code for analysis/simulations related to missing eligibility data in electronic health record based observational studies using target trial emulation. For more information, see our manuscript:

Benz, L., Mukherjee, R., Dahabreh, I., Wang, R., Arterburn, D., Fischer, H., Lee, C., Shortreed, S., and Hauense, S. "Adjusting for Selection Bias Due to Missing Eligibility Criteria in Emulated Target Trials." _Under Review_ , 2024. [(Pre-Print)](https://arxiv.org/abs/2406.16830)

## R Scripts (`scripts/`)

### Data (`scripts/data`)

* __build_complete_case_T2DM_population.R__: Build dataset which will define the notion of a complete case population for simulations.
* __clean_weights.R__: Cleans weights for all measures. 
* __microvascular_dataset.R__: Creates dataset for replicating (O'Brien, 2018)
* __microvascular_dataset_tte.R__: Creates datasets necessary for replicating (O'Brien, 2018) in TTE format
* __save_parquet.R__: Script to save out some of the larger .sas7bdat files in .parquet format for fast reading

### Utility Functions (`scripts/util/`)

* __helpers.R__: Useful helper functions
* __img_to_pdf.R__: Script to combine all .png image files into 1 PDF
* __fit_models.R__: Fits and saves out model objects

### Exploratory Data Analyses (`scripts/eda`)
* __f31_figures.R__: Figures/table code

### Simulations (`scripts/simulations`)
* __helpers.R__: Helper functions for simulations

#### Missing Data Simulations (`scripts/simulations/missing_data`)
Simulations based on missing data in the eligibility criteria

* __weight_trajectories_function.R__: Functions to generate simulated weight trajectory functions
* __eda_plot_trajectories.R__: Function to make some plots to explore the weight trajectories
* __generate_data.R__: Functions to generate simulated data
* __specify_inputs.R__: Where inputs for simulations are specified
* __run_simulation.R__: Wrapper to run the simulations
* __compute_truth.R__: Compute true ATE
* __fit_analysis_models.R__: Fit analysis model to estimate the ATE
* __analysis.R__: Analyze results
* __explore_sim_params.R__: Scipt to better understand underlying data generation process for selction simulation settings
* __scratch.R__: Scratch pad of snippets useful for better understand underlying data generation process for selction simulation settings
* __boostrap_variance_sim.R__: Simulation for testing the validity of the boostrap for variance.
* __save_boostrap_simulated_data.R__: Save out data files for bootstrap variance simulation
* __analyze_boostrap_variance.R__: Analyze results of boostrap validity

### Analaysis (`scripts/analysis`)
* __microvascular_survival.R__: Analyze survival for data application mirroring (O'Brien 2018)
* __microvascular_tte_pp.R__: Target Trial Emulation for data application mirroring (O'Brien 2018) (PP)
* __microvascular_tte_itt.R__: Target Trial Emulation for data application mirroring (O'Brien 2018) (ITT)
* __boostrap_microvascular_iit.R__: Show boostrapped variance over full grid for 100 replicates of ITT
* __plot_HR.R__: Function to plot discrete hazard ratios.
* __plot_missing_data.R__: Function to summarize missing data as a function of lookbacks.
* __bootstrap_ITT__: Final boostrap (1000 replicates) for ITT based on lookbacks of 3 (BMI) and 12 (A1c)
* __bootstrap_PP__: Final boostrap (1000 replicates) for PPe based on lookbacks of 3 (BMI) and 12 (A1c)

## Data (`data/`)
Due to privacy concerns, the raw/derived datasets can not be stored locally or in GitHub. 

* __simulations/__: Folder of simulation inputs and outputs
* __microvascular_results/__: Results for effect of bariatric surgery on microvascular outcomes from sensitivity analysis.


## Figures (`figures/`)
Figures saved out from various analyses

## Shell (`shell/`)
.sh files for transferring files to and from the cluster.

## Jobs (`jobs/`)
.sh files for batch jobs on the cluster 

* __missing_data_sims.sh__: Job file for missing data simulations
* __microvascular_pp.sh__: Run analysis for target trial emulation for microvascular outcomes (PP)
* __microvascular_itt.sh__: Run analysis for target trial emulation for microvascular outcomes (ITT)
* __boot_microvasular_itt.sh__: Compute boostrapped variance over full grid for 100 replicates of ITT
* __boot_final_itt.sh__: Compute boostrapped variance for final choice of lookbacks for 1000 replicates of ITT
* __boot_final_pp.sh__: Compute boostrapped variance for final choice of lookbacks for 1000 replicates of PP
* __data_prep_bootsrap_variance_sims.sh__: Data prep for bootstrap variance simulations
* __bootstrap_variance_sim.sh__: Bootstrap variance job file
* __boot_variance_all.sh__: Wrapper for all boot variance jobs

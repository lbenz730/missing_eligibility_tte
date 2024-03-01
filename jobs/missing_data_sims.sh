#!/usr/bin/bash

#SBATCH -c 12 ## number of cores
#SBATCH -t 0-12:00 ## amount of time in D-HH:MM
#SBATCH -p fasse_bigmem  ## Partition to submit to
#SBATCH --mem=450000 ## memory pool for all cores
#SBATCH -o logs/missing_data/log.stdout_%a ## STDOUT
#SBATCH -e logs/missing_data/log.stderr_%a ## STDERR
#SBATCH --account=haneuse_lab

module load R/4.2.2-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER

cd $HOME/target_trial_emulation/

Rscript scripts/simulations/missing_data/run_simulation.R $SLURM_ARRAY_TASK_ID
#!/usr/bin/bash

#SBATCH -c 48 ## number of cores
#SBATCH -t 0-01:00 ## amount of time in D-HH:MM
#SBATCH -p fasse_bigmem  ## Partition to submit to
#SBATCH --mem=200000 ## memory pool for all cores
#SBATCH -o logs/bootstrap_valid/log.stdout_%a ## STDOUT
#SBATCH -e logs/bootstrap_valid/log.stderr_%a ## STDERR
#SBATCH --account=haneuse_lab
#SBATCH --array=1-1000

module load R/4.2.2-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER

cd $HOME/target_trial_emulation/

Rscript scripts/simulations/missing_data/bootstrap_variance_sim.R $SLURM_ARRAY_TASK_ID $*
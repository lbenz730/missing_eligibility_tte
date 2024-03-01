#!/usr/bin/bash

#SBATCH -c 1 ## number of cores
#SBATCH -t 0-1:30 ## amount of time in D-HH:MM
#SBATCH -p fasse ## Partition to submit to
#SBATCH --mem=184000 ## memory pool for all cores
#SBATCH -o logs/boot_PP/log.stdout_%a ## STDOUT
#SBATCH -e logs/boot_PP/log.stderr_%a ## STDERR
#SBATCH --account=haneuse_lab
#SBATCH --array=1-4000

module load R/4.2.2-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER

cd $HOME/target_trial_emulation/

Rscript scripts/analysis/bootstrap_PP.R $SLURM_ARRAY_TASK_ID
#!/bin/bash
#SBATCH --job-name=procrm_calibrate_priors
#SBATCH --output=procrm_calibrate_priors.txt
#SBATCH --error=procrm_calibrate_priors.txt

#SBATCH --cpus-per-task=4
#SBATCH --time=01:30:00
#SBATCH --mail-user=emily.alger@icr.ac.uk
#SBATCH --mail-type=ALL

source ~/.bashrc
mamba activate mamba-emily

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun Rscript /home/ealger/revision_calibrate_priors/procrm_sim_HPC.R

#!/bin/bash
#SBATCH --job-name=calibrate_priors
#SBATCH --output=calibrate_priors.txt
#SBATCH --error=calibrate_priors.txt

#SBATCH --cpus-per-task=4
#SBATCH --time=01:30:00
#SBATCH --mail-user=emily.alger@icr.ac.uk
#SBATCH --mail-type=ALL

source ~/.bashrc
mamba activate mamba-emily

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun Rscript /home/ealger/revision_calibrate_priors/efftox_sim_HPC.R

#!/bin/bash
#SBATCH -A matgen_g
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 04:00:00
#SBATCH -n 4
#SBATCH --gpus-per-task=1

source activate openmm
srun python submit.py ../mixed_sims/other $SCRATCH/mixed_sims/other

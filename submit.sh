#!/bin/bash
#SBATCH -A matgen_g
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 00:05:00
#SBATCH -n 4
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-task=1

source activate openmm
srun python production_sim.py
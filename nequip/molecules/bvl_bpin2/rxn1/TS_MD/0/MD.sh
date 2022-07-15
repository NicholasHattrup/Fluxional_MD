#!/bin/bash
#SBATCH -A orca
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -t 12:00:00
#SBATCH -p RM-shared
#SBATCH --job-name MD_0

module load orca


/opt/packages/orca/orca_5_0_1_linux_x86-64_shared_openmpi411/orca MD_0.inp

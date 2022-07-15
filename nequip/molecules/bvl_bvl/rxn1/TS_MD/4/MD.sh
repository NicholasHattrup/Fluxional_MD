#!/bin/bash
#SBATCH -A CHE220025P
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -t 4:00:00
#SBATCH -p RM-shared
#SBATCH --job-name MD_4

module load orca


/opt/packages/orca/orca_5_0_1_linux_x86-64_shared_openmpi411/orca MD_4.inp

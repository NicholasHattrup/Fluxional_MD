#!/bin/bash
#SBATCH -A orca
#SBATCH -p RM-shared
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=16

module load orca

$(which orca) get_bvl_bvl_rxn1_TS.inp

#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH --mem=2g
# by Valeria Añorve-Garibay
# code to simulate any evolutionary dynamic based on SLiM scripts
# TYPE defines the type of evolutionary scenario simulated.
# The three options are "neutral_evolution", "stabilizing_selection" and "directional_selection" which defines 
# if the trait is evolving under neutral evolution, stabilizing selection or directional selection, respectively.
# hsq stands for heritability value, assign accordingly
# path corresponds to the desired output dir
# nrun assigns a run ID for the specific replicate, e.g. sim_``1``.slim based on SLURM_ARRAY_TASK_ID
# w represents the strength of selection, assign accordingly. we used w={1,2,3,4,5} to represent cases going from stronger to weaker stabilizing selection
# QTL_sd represents the standard deviation of the QTLs distribution (normal with mean 0 and sd QTL_sd)
# e.g. slurm command
# neutral evolution with h2 = 1
# sbatch --array=1 sim_any_evo_dynamic.sh neutral_evolution 1.0 neutral_evolution/
# stabilizing selection with h2 = 1 and w = 1
# sbatch --array=1 sim_any_evo_dynamic.sh stabilizing_selection 1.0 1.0 stabilizing_selection/
# directional selection with h2 = 1 and QTLs sd = 0.25
# sbatch --array=1 sim_any_evo_dynamic.sh directional_selection 1.0 0.25

module load slim/4.0.1-kymgtmu

TYPE=$1

if [ "$TYPE" == "neutral_evolution" ]; then
  echo "Running code for $TYPE"
  nrun=${SLURM_ARRAY_TASK_ID}
  path=$3
  slim -d hsq=$2 -d nrun=${nrun} -d path="'${path}'" neutral_evolution.slim > ${path}/OUT.${nrun}

elif [ "$TYPE" == "stabilizing_selection" ]; then
  echo "Running code for $TYPE"
  nrun=${SLURM_ARRAY_TASK_ID}
  path=$4
  slim -d hsq=$2 -d w=$3 -d nrun=${nrun} -d path="'${path}'" stabilizing_selection.slim > ${path}/OUT.${nrun}

elif [ "$TYPE" == "directional_selection" ]; then
  echo "Running code for $TYPE"
  nrun=${SLURM_ARRAY_TASK_ID}
  path=$4
  slim -d hsq=$2 -d QTL_sd=$3 -d nrun=${nrun} -d path="'${path}'" directional_selection.slim > ${path}/OUT.${nrun}
  
fi

#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH --mem=2g
module load slim/4.0.1-kymgtmu
# code to simulate any evolutionary dynamic based on SLiM scripts
# TYPE defines the type of evolutionary scenario simulated.
# The three options are "neutral_evolution", "stabilizing_selection" and "directional_selection" which defines 
# if the trait is evolving under neutral evolution, stabilizing selection or directional selection, respectively.

TYPE=$1

if [ "$TYPE" == "neutral_evolution" ]; then
  echo "Running code for $TYPE"
  nrun=${SGE_TASK_ID}
  path=$3
  slim -d hsq=$2 -d nrun=${nrun} -d path="'${path}'" neutral_evolution.slim > ${path}/OUT.${nrun}

elif [ "$TYPE" == "stabilizing_selection" ]; then
  echo "Running code for $TYPE"
  nrun=${SGE_TASK_ID}
  path=$4
  slim -d hsq=$2 -d w=$3 -d nrun=${nrun} -d path="'${path}'" stabilizing_selection.slim > ${path}/OUT.${nrun}

elif [ "$TYPE" == "directional_selection" ]; then
  echo "Running code for $TYPE"
  nrun=${SGE_TASK_ID}
  path=$4
  slim -d hsq=$2 -d QTL_sd=$3 -d nrun=${nrun} -d path="'${path}'" directional_selection.slim > ${path}/OUT.${nrun}
  
fi

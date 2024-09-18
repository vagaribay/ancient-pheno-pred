#!/bin/bash
# code to simulate any evolutionary dynamic based on SLiM scripts
TYPE=$1

if [ "$TYPE" == "neutral_evolution" ]; then
  echo "Running code for $TYPE"
  #nrun=${SGE_TASK_ID}
  nrun=$4
  path=$3
  slim -d hsq=$2 -d nrun=${nrun} -d path="'${path}'" neutral_evolution.slim > ${path}/OUT.${nrun}

elif [ "$TYPE" == "stabilizing_selection" ]; then
  echo "Running code for $TYPE"
  #nrun=${SGE_TASK_ID}
  nrun=$5
  path=$4
  slim -d hsq=$2 -d w=$3 -d nrun=${nrun} -d path="'${path}'" stabilizing_selection.slim > ${path}/OUT.${nrun}

elif [ "$TYPE" == "directional_selection" ]; then
  echo "Running code for $TYPE"
  #nrun=${SGE_TASK_ID}
  nrun=$5
  path=$4
  slim -d hsq=$2 -d QTL_sd=$3 -d nrun=${nrun} -d path="'${path}'" directional_selection.slim > ${path}/OUT.${nrun}
  
fi

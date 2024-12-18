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

sbatch --array=1 sim_any_evo_dynamic.sh neutral_evolution 0.5 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/neutral_evolution/

sbatch --array=1 sim_any_evo_dynamic.sh stabilizing_selection 0.5 1.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/w_1/
sbatch --array=1 sim_any_evo_dynamic.sh stabilizing_selection 0.5 5.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/w_5/

sbatch --array=1 sim_any_evo_dynamic.sh directional_selection 0.5 0.25 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/directional_selection/sd_025/
sbatch --array=1 sim_any_evo_dynamic.sh directional_selection 0.5 0.0025 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/directional_selection/sd_00025/





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


sbatch --array=1-100 sim_any_evo_dynamic.sh neutral_evolution 0.5 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/neutral_evolution/mid/


sbatch --array=1 sim_any_evo_dynamic.sh stabilizing_selection 0.5 1.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/w_1/
sbatch --array=1 sim_any_evo_dynamic.sh stabilizing_selection 0.5 5.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/w_5/

sbatch --array=1 sim_any_evo_dynamic.sh directional_selection 0.5 0.25 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/directional_selection/sd_025/
sbatch --array=1 sim_any_evo_dynamic.sh directional_selection 0.5 0.0025 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/directional_selection/sd_00025/



sbatch --array=1-100 sim_any_evo_dynamic.sh stabilizing_selection 1.0 2.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/total/2/
sbatch --array=1-100 sim_any_evo_dynamic.sh stabilizing_selection 1.0 3.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/total/3/
sbatch --array=1-100 sim_any_evo_dynamic.sh stabilizing_selection 1.0 4.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/total/4/
sbatch --array=1-100 sim_any_evo_dynamic.sh stabilizing_selection 1.0 5.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/total/5/

sbatch --array=1-100 sim_any_evo_dynamic.sh stabilizing_selection 0.5 1.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/mid/1/
sbatch --array=1-100 sim_any_evo_dynamic.sh stabilizing_selection 0.5 2.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/mid/2/
sbatch --array=1-100 sim_any_evo_dynamic.sh stabilizing_selection 0.5 3.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/mid/3/
sbatch --array=1-100 sim_any_evo_dynamic.sh stabilizing_selection 0.5 4.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/mid/4/
sbatch --array=1-100 sim_any_evo_dynamic.sh stabilizing_selection 0.5 5.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/mid/5/




sbatch --array=33,55,8 sim_any_evo_dynamic.sh directional_selection 1.0 0.25 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/directional_selection/total/025/

sbatch --array=24,34,6,68,8,90,98,101-105 sim_any_evo_dynamic.sh directional_selection 1.0 0.0025 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/directional_selection/total/00025/

sbatch --array=93 sim_any_evo_dynamic.sh directional_selection 0.5 0.25 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/directional_selection/mid/025/
sbatch --array=1-100 sim_any_evo_dynamic.sh directional_selection 0.5 0.025 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/directional_selection/mid/0025/
sbatch --array=1,16,26,101-110 sim_any_evo_dynamic.sh directional_selection 0.5 0.0025 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/directional_selection/mid/00025/


sbatch --array=69,101 sim_any_evo_dynamic.sh stabilizing_selection 1.0 1.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/total/1/
sbatch --array=48 sim_any_evo_dynamic.sh stabilizing_selection 1.0 5.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/total/5/



sbatch --array=76,52 sim_any_evo_dynamic.sh stabilizing_selection 0.5 1.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/mid/1/
sbatch --array=91 sim_any_evo_dynamic.sh stabilizing_selection 0.5 4.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/mid/4/
sbatch --array=42 sim_any_evo_dynamic.sh stabilizing_selection 0.5 5.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/mid/5/


sbatch --array=15,30,66,101-105 sim_any_evo_dynamic.sh neutral_evolution 1.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/neutral_evolution/total/
sbatch --array=19,2,30,34,44,6 sim_any_evo_dynamic.sh neutral_evolution 0.5 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/neutral_evolution/mid/



sbatch --array=56 sim_any_evo_dynamic.sh stabilizing_selection 1.0 1.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stabilizing_selection/total/1/


sbatch --array=22,34,43,49,51,7,79,85,92,93,95,96,101-110 sim_any_evo_dynamic.sh stabilizing_selection 0.5 1.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stab_big_gsize/mid/1/
sbatch --array=86,101 sim_any_evo_dynamic.sh stabilizing_selection 0.5 5.0 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/stab_big_gsize/mid/5/


sbatch --array=11-120 sim_any_evo_dynamic.sh stabilizing_selection 0.8 7.28 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/realistic_params/HT/
sbatch --array=11-120 sim_any_evo_dynamic.sh stabilizing_selection 0.7 6.61 /users/vanorveg/data/vanorveg/ancient-pheno-pred/output/realistic_params/BMI/













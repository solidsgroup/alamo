#!/usr/bin/env bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem-per-cpu=10000
#SBATCH --job-name="propellant_study"
#SBATCH --array=0-5
#SBATCH --output="%x-%A_%a-log.txt"
#SBATCH --mail-user=mungerct@iastate.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --constraint="nova25|nova24"
echo "======================================================"
echo " Job '$SLURM_JOB_NAME'"
echo " Job ID: $SLURM_JOB_ID"
echo " Array Task ID: $SLURM_ARRAY_TASK_ID"
echo " Submitted by: $SLURM_JOB_USER"
echo " Running on node(s): $SLURM_NODELIST"
echo " Simulation with plot file=${plot_file_name} started"
echo " Start time: $(date)"
echo "======================================================"
module load openmpi_hpc
# -------------------------------------------------
# Define array of zeta values (non-integers allowed)
# -------------------------------------------------
ind=00
num=39
AP_arr=(
    "I/I_${ind}_AP.xyzr"
    "J/J_${ind}_AP.xyzr"
    "K/K_${ind}_AP.xyzr"
    "R/R_${ind}_AP.xyzr"
    "S/S_${ind}_AP.xyzr"
    "T/T_${ind}_AP.xyzr"
)
plot_file_arr=(
    "I_${ind}_novoid_${num}"
    "J_${ind}_novoid_${num}"
    "K_${ind}_novoid_${num}"
    "R_${ind}_novoid_${num}"
    "S_${ind}_novoid_${num}"
    "T_${ind}_novoid_${num}"
)

# Get value corresponding to this array task
# zeta=${zeta_values[$SLURM_ARRAY_TASK_ID]}

AP=${AP_arr[$SLURM_ARRAY_TASK_ID]}
plot_file_name=${plot_file_arr[$SLURM_ARRAY_TASK_ID]}

# echo "Running simulation with zeta = ${zeta}"
# Optional: still use IC index if needed

i=${SLURM_ARRAY_TASK_ID}
# srun --mpi=pmix ./bin/alamo-2d-g++ input.scpspheres_kodga_validation \
#     plot_file="output.A_0_novoid_zeta_${zeta}" \
#     propellant.fullfeedback.phi.zeta=${zeta}
srun --mpi=pmix ./bin/alamo-2d-g++ input.scpspheres_kodga_validation \
    plot_file="output.${plot_file_name}" \
    phi.ic.psread.file.name="../papers/RegressionWithVoidsFullFeedback/results/new_datasets/${AP}" \
    pf.eta.ic.type=constant \
    pf.eta.ic.constant.value=1.0

echo "======================================================"
echo " Simulation with plot file=${plot_file_name} completed"
echo " End time: $(date)"
echo "======================================================"

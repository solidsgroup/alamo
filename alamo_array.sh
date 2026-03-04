#!/usr/bin/env bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem-per-cpu=10000
#SBATCH --job-name="propellant_study"
#SBATCH --array=0-7
#SBATCH --output="%x-%A_%a-log.txt"
#SBATCH --mail-user=mungerct@iastate.edu
#SBATCH --mail-type=BEGIN,END,FAIL
echo "======================================================"
echo " Job '$SLURM_JOB_NAME'"
echo " Job ID: $SLURM_JOB_ID"
echo " Array Task ID: $SLURM_ARRAY_TASK_ID"
echo " Submitted by: $SLURM_JOB_USER"
echo " Running on node(s): $SLURM_NODELIST"
echo " Start time: $(date)"
echo "======================================================"
module load openmpi_hpc
# -------------------------------------------------
# Define array of zeta values (non-integers allowed)
# -------------------------------------------------
ind=00
num=6
AP_arr=(
    "A_vf_11/A_${ind}_AP.xyzr"
    "B_vf_6/B_${ind}_AP.xyzr"
    "C_vf_6/C_${ind}_AP.xyzr"
    "D_vf_7/D_${ind}_AP.xyzr"
    "E_vf_7/E_${ind}_AP.xyzr"
    "F_vf_6/F_${ind}_AP.xyzr"
    "G_vf_14/G_${ind}_AP.xyzr"
    "H_vf_13/H_${ind}_AP.xyzr"
)
eta_arr=(
    "A_vf_11/A_${ind}_void.xyzr"
    "B_vf_6/B_${ind}_void.xyzr"
    "C_vf_6/C_${ind}_void.xyzr"
    "D_vf_7/D_${ind}_void.xyzr"
    "E_vf_7/E_${ind}_void.xyzr"
    "F_vf_6/F_${ind}_void.xyzr"
    "G_vf_14/G_${ind}_void.xyzr"
    "H_vf_13/H_${ind}_void.xyzr"
)
plot_file_arr=(
    "A_vf_11_${ind}_void_${num}"
    "B_vf_6_${ind}_void_${num}"
    "C_vf_6_${ind}_void_${num}"
    "D_vf_7_${ind}_void_${num}"
    "E_vf_7_${ind}_void_${num}"
    "F_vf_6_${ind}_void_${num}"
    "G_vf_14_${ind}_void_${num}"
    "H_vf_13_${ind}_void_${num}"
)

# Get value corresponding to this array task
# zeta=${zeta_values[$SLURM_ARRAY_TASK_ID]}

AP=${AP_arr[$SLURM_ARRAY_TASK_ID]}
eta=${eta_arr[$SLURM_ARRAY_TASK_ID]}
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
    pf.eta.ic.psread.file.name="../papers/RegressionWithVoidsFullFeedback/results/new_datasets/${eta}"

echo "======================================================"
echo " Simulation with plot file=${plot_file_name} completed"
echo " End time: $(date)"
echo "======================================================"

#!/usr/bin/env bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem-per-cpu=15000
#SBATCH --job-name="propellant_study"
#SBATCH --array=3-4
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
# zeta_values=(5.0_um 10.0_um 15.0_um 20.0_um 25.0_um 30.0_um)
AP_arr=(D_vf_7/D_00_AP.xyzr E_vf_7/E_00_AP.xyzr F_vf_6/F_00_AP.xyzr G_vf_14/G_00_AP.xyzr H_vf_13/H_00_AP.xyzr)
eta_arr=(D_vf_7/D_00_void.xyzr E_vf_7/E_00_void.xyzr F_vf_6/F_00_void.xyzr G_vf_14/G_00_void.xyzr H_vf_13/H_00_void.xyzr)
plot_file_arr=(D_vf_7_0_void_1 E_vf_7_0_void_1 F_vf_6_0_void_1 G_vf_14_0_void_1 H_vf_13_0_void_1)

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

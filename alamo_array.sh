#!/usr/bin/env bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem-per-cpu=15000
#SBATCH --job-name="zeta_study"
#SBATCH --array=0-5
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
zeta_values=(5.0, 10.0, 15.0, 20.0, 25.0, 30.0)

# Get value corresponding to this array task
zeta=${zeta_values[$SLURM_ARRAY_TASK_ID]}

echo "Running simulation with zeta = ${zeta}"

# Optional: still use IC index if needed
i=${SLURM_ARRAY_TASK_ID}

srun --mpi=pmix ./bin/alamo-2d-g++ input.scpspheres_kodga_validation \
    plot_file="output.A_0_novoid_zeta_${zeta}" \
    propellant.fullfeedback.phi.zeta=${zeta}

echo "======================================================"
echo " Simulation with zeta=${zeta} completed"
echo " End time: $(date)"
echo "======================================================"

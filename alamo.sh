#!/usr/bin/env bash

#SBATCH --time=500:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name="SCP_half_sandwich_TPG_mass_momentum"
#SBATCH --output="%x-%j-log.txt"
#SBATCH --mail-user=mungerct@iastate.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --constraint=nova25

JOB_NAME="SCP_half_sandwich_TPG_mass_momentum"
PLOT_FILE="output.${JOB_NAME}"

echo "======================================================"
echo " Job '$SLURM_JOB_NAME' (ID: $SLURM_JOB_ID) is starting"
echo " Submitted by: $SLURM_JOB_USER"
echo " Running on node(s): $SLURM_NODELIST"
echo " Start time: $(date)"
echo "======================================================"

module load openmpi_hpc

srun --mpi=pmix ./bin/alamo-2d-6species-hdf5-g++ input.sandwich_rocfire_flame_diffuse  plot_file="${PLOT_FILE}"

echo "======================================================"
echo " Job finished at: $(date)"
echo "======================================================"

#!/usr/bin/env bash

#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=64G
#SBATCH --job-name="airfoil_supersonic_2"
#SBATCH --output="%x-%j-log.txt"
#SBATCH --mail-user=mungerct@iastate.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --constraint=nova25

JOB_NAME="airfoil_supersonic_2"
PLOT_FILE="output.${JOB_NAME}"

echo "======================================================"
echo " Job '$SLURM_JOB_NAME' (ID: $SLURM_JOB_ID) is starting"
echo " Submitted by: $SLURM_JOB_USER"
echo " Running on node(s): $SLURM_NODELIST"
echo " Start time: $(date)"
echo "======================================================"

module load openmpi_hpc

srun --mpi=pmix ./bin/sfi-2d-hdf5-g++ input.flow_airfoil_subsonic plot_file="${PLOT_FILE}"
# srun --mpi=pmix ./bin/sfi-2d-g++ input.flow_airfoil_supersonic restart=./output.flow_airfoil_supersonic/27580cell restart_cell=./output.flow_airfoil_supersonic/27580cell plot_file="${PLOT_FILE}"

echo "======================================================"
echo " Job finished at: $(date)"
echo "======================================================"

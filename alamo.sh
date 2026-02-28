#!/usr/bin/env bash

#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem-per-cpu=15000
#SBATCH --job-name="A_1_void_m3"
#SBATCH --output="%x-%j-log.txt"
#SBATCH --mail-user=mungerct@iastate.edu
#SBATCH --mail-type=BEGIN,END,FAIL
# #SBATCH --constraint=nova25

JOB_NAME="A_1_void_M3"
PLOT_FILE="output.${JOB_NAME}"

echo "======================================================"
echo " Job '$SLURM_JOB_NAME' (ID: $SLURM_JOB_ID) is starting"
echo " Submitted by: $SLURM_JOB_USER"
echo " Running on node(s): $SLURM_NODELIST"
echo " Start time: $(date)"
echo "======================================================"

module load openmpi_hpc

srun --mpi=pmix ./bin/alamo-2d-g++ \
    input.scpspheres_kodga_validation \
    plot_file="${PLOT_FILE}"

echo "======================================================"
echo " Job finished at: $(date)"
echo "======================================================"

#!/usr/bin/env bash

#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem-per-cpu=15000
#SBATCH --job-name="arb_geo_Kohga2006_4"
#SBATCH --output="%x-%j-log.txt"
#SBATCH --mail-user=mungerct@iastate.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --constraint=nova25

JOB_NAME="arb_geo_Kohga2006_4"
PLOT_FILE="output.${JOB_NAME}"

echo "======================================================"
echo " Job '$SLURM_JOB_NAME' (ID: $SLURM_JOB_ID) is starting"
echo " Submitted by: $SLURM_JOB_USER"
echo " Running on node(s): $SLURM_NODELIST"
echo " Start time: $(date)"
echo "======================================================"

module load openmpi_hpc

srun --mpi=pmix ./bin/alamo-2d-g++ input.scpspheres_arbitary_geometry plot_file="${PLOT_FILE}"

echo "======================================================"
echo " Job finished at: $(date)"
echo "======================================================"

#!/usr/bin/env bash
#SBATCH --time=14:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=10000
#SBATCH --job-name="alamo"
#SBATCH --array=2-10
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

# IC file index from array task ID
i=${SLURM_ARRAY_TASK_ID}

echo "Running simulation with IC file ${i}"

srun --mpi=pmix ./bin/alamo-2d-g++ input.scpspheres_kodga_validation \
    phi.ic.psread.file.name="all_datasets/C/setC_xyzrs/setC_${i}_AP.xyzr" \
    plot_file="output.scpspheres_novoids_set${i}"

echo "======================================================"
echo " Simulation ${i} completed"
echo " End time: $(date)"
echo "======================================================"

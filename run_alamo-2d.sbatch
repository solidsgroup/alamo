#!/usr/bin/env bash

#SBATCH --account=brunnels
#SBATCH --chdir=/work/brunnels/jstocker/alamo/reg_with_voids
#SBATCH --constraint=nova25
#SBATCH --cpus-per-task=1
#SBATCH --job-name=alamo-2d-reg_with_voids
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jstocker@iastate.edu
#SBATCH --mem-per-cpu=1G
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --output=%x-%j.log
#SBATCH --partition=nova
#SBATCH --time=24:00:00

HOST="$(hostname)"
START_TIME=$(date +%s)

INPUT_FILE="/work/brunnels/jstocker/alamo/reg_with_voids/input.scpspheres_arbitrary_with_voids"
OUTPUT_DIR="/work/brunnels/jstocker/alamo/reg_with_voids/output.reg_with_voids_2d-arb_geo$SLURM_JOB_ID"
LOG_FILE="$(dirname $OUTPUT_DIR)/log-$(basename "$OUTPUT_DIR").txt"

echo "=================================================================================="
echo "Starting Slurm job $SLURM_JOB_ID on $HOST"
echo "Job name: $SLURM_JOB_NAME"
echo "Start time: $START_TIME"
echo "Output directory: $OUTPUT_DIR"
echo "=================================================================================="

# Load modules
module purge
module load openmpi_hpc
# Run the simulation with output redirected to log file
echo ""
echo "Running Alamo simulation..."
echo "Simulation output will be written to: $LOG_FILE"

# Write Slurm metadata to log file
{
    echo "=================================================================================="
    echo "SLURM JOB METADATA"
    echo "=================================================================================="
    echo "Job ID: $SLURM_JOB_ID"
    echo "Job Name: $SLURM_JOB_NAME"
    echo "Node List: $SLURM_JOB_NODELIST"
    echo "Partition: $SLURM_JOB_PARTITION"
    echo "Account: $SLURM_JOB_ACCOUNT"
    echo "Number of Nodes: $SLURM_JOB_NUM_NODES"
    echo "Number of Tasks: $SLURM_NTASKS"
    echo "CPUs per Task: $SLURM_CPUS_PER_TASK"
    echo "Memory per CPU: $SLURM_MEM_PER_CPU"
    echo "Start Time: $START_TIME"
    echo "Host: $HOST"
    echo "=================================================================================="
    echo ""
    mkdir -p "$OUTPUT_DIR"
    # Run the actual simulation
    srun --mpi=pmix ../bin/alamo-2d-g++ $INPUT_FILE plot_file="$OUTPUT_DIR"
} > "$LOG_FILE" 2>&1
SIMULATION_EXIT_CODE=$?

SIMULATION_END_TIME=$(date +%s)
SIMULATION_ELAPSED=$(( SIMULATION_END_TIME - START_TIME ))

echo ""
echo "=================================================================================="
echo "Simulation completed with exit code: $SIMULATION_EXIT_CODE"
echo "Simulation runtime: $((SIMULATION_ELAPSED / 60)) minutes ($SIMULATION_ELAPSED seconds)"
echo "=================================================================================="
echo ""

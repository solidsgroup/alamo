#!/usr/bin/env bash
#SBATCH --partition=nova
#SBATCH --constraint=zen5
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=64G
#SBATCH --job-name=sandwich-no-thermal
#SBATCH --output=sandwich-no-thermal-%j.log
set -euo pipefail
cd /work/brunnels/mungerct/alamo
module purge
module load openmpi gcc eigen
export OMP_NUM_THREADS=1
exec mpirun -np 64 ./bin/lowmach-2d-g++ input.lm.ap_htpb_sandwich_neo_hookean plot_file=output.lm.ap_htpb_sandwich_neo_hookean_${SLURM_JOB_ID}

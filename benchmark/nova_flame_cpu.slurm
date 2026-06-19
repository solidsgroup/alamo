#!/bin/bash
#SBATCH -A brunnels
#SBATCH -J flame_cpu
#SBATCH --partition=nova
#SBATCH --nodes=1
#SBATCH --ntasks=16                 # MPI ranks = domain decomposition
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=flame_cpu.%j.out
#SBATCH --error=flame_cpu.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jackplum@iastate.edu
# ===========================================================================
# CPU baseline run of the Flame center-bore case on NOVA (for the CPU-vs-GPU
# comparison). Builds the CPU binary in-job if missing.
#
#   sbatch benchmark/nova_flame_cpu.slurm
# ===========================================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

INPUT="${INPUT:-input_nova_centre_bore}"

module purge 2>/dev/null || true
module load gcc 2>/dev/null || module load gcc/12 2>/dev/null || true
module load openmpi 2>/dev/null || module load openmpi4 2>/dev/null || true

# --- Build CPU binary if missing (profiled, optimized) ----------------------
CPU_BIN="$(ls -t bin/alamo-2d*-g++ 2>/dev/null | grep -v cuda | head -1 || true)"
if [ -z "${CPU_BIN}" ]; then
  echo "No CPU binary -- building (./configure --comp=g++ --profile && make)"
  ./configure --comp=g++ --profile
  make -j"${SLURM_CPUS_ON_NODE:-16}"
  CPU_BIN="$(ls -t bin/alamo-2d*-g++ | grep -v cuda | head -1)"
fi
echo "CPU binary: ${CPU_BIN}  (${SLURM_NTASKS} ranks)"

# async_out gives compute/IO overlap; NOVA's MPI provides MPI_THREAD_MULTIPLE.
OVERRIDES="amrex.async_out=1 amrex.async_out_nfiles=0 \
           tiny_profiler.device_synchronize_around_region=1"

echo "=== launching ${SLURM_NTASKS} MPI rank(s) ==="
time srun --mpi=pmix -n "${SLURM_NTASKS}" "${CPU_BIN}" "${INPUT}" ${OVERRIDES}
echo "=== done ==="

#!/bin/bash
#SBATCH -A brunnels
#SBATCH -J flame_cpu_3d
#SBATCH --partition=nova
#SBATCH --nodes=1
#SBATCH --ntasks=16                 # MPI ranks = domain decomposition
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=flame_cpu_3d.%j.out
#SBATCH --error=flame_cpu_3d.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jackplum@iastate.edu
# ===========================================================================
# CPU baseline 3D run of the Flame case on NOVA (for the GPU-vs-CPU-node
# crossover comparison). Builds the CPU binary in-job if missing unless
# BUILD_CPU_IF_MISSING=0 is exported by the one-shot sweep driver.
# This is the 3D sibling of nova_flame_cpu.slurm (--dim 3 + input_3d_flame).
#
#   sbatch benchmark/nova_flame_cpu_3d.slurm
# ===========================================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

INPUT="${INPUT:-input_3d_flame}"
BUILD_CPU_IF_MISSING="${BUILD_CPU_IF_MISSING:-1}"

module purge 2>/dev/null || true
module load gcc 2>/dev/null || module load gcc/12 2>/dev/null || true
module load openmpi 2>/dev/null || module load openmpi4 2>/dev/null || true

# --- Build CPU binary if missing (profiled, optimized) ----------------------
CPU_BIN="$(ls -t bin/alamo-3d*-g++ 2>/dev/null | grep -v cuda | head -1 || true)"
if [ -z "${CPU_BIN}" ]; then
  if [ "${BUILD_CPU_IF_MISSING}" != "1" ]; then
    echo "No bin/alamo-3d*-g++ CPU binary found and BUILD_CPU_IF_MISSING=${BUILD_CPU_IF_MISSING}."
    echo "The Phase-3 one-shot should have submitted alamo_cpu_build_3d before this job."
    exit 1
  fi
  echo "No CPU binary -- building (./configure --comp=g++ --dim 3 --profile --get-eigen && make)"
  ./configure --comp=g++ --dim 3 --profile --get-eigen
  make -j"${SLURM_CPUS_ON_NODE:-16}"
  CPU_BIN="$(ls -t bin/alamo-3d*-g++ | grep -v cuda | head -1)"
fi
echo "CPU binary: ${CPU_BIN}  (${SLURM_NTASKS} ranks)"

# async_out gives compute/IO overlap; NOVA's MPI provides MPI_THREAD_MULTIPLE.
OVERRIDES="amrex.async_out=1 amrex.async_out_nfiles=0 \
           tiny_profiler.device_synchronize_around_region=1"

echo "=== launching ${SLURM_NTASKS} MPI rank(s) ==="
time srun --mpi=pmix -n "${SLURM_NTASKS}" "${CPU_BIN}" "${INPUT}" ${OVERRIDES}
echo "=== done ==="

#!/bin/bash
#SBATCH -A brunnels
#SBATCH -J flame_cpu_3d_a2
#SBATCH --partition=nova
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=flame_cpu_3d_a2.%j.out
#SBATCH --error=flame_cpu_3d_a2.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jackplum@iastate.edu
# The partition defaults to nova; pass --partition=scavenger to sbatch to use
# the scavenger queue (or use phase_a2_sweep.sh --scavenger --submit):
#   sbatch benchmark/nova_flame_cpu_3d_a2.slurm
#   sbatch --partition=scavenger benchmark/nova_flame_cpu_3d_a2.slurm
# ===========================================================================
# A2 CPU baseline 3D Flame+Elastic run on NOVA (GPU_ROADMAP_V2 Phase A2).
#
# Provides the CPU-node reference for the combined flame+elastic crossover.
# 64 MPI ranks matches the Phase-3 R3 CPU baseline (06-22 run), ensuring a
# fair comparison.  8-hour limit accounts for the added elastic MLMG cost
# that was absent in Phase-3.
#
# Inputs: input_3d_centre_bore_128_a2, input_3d_centre_bore_256_a2
#
# Env knobs:
#   INPUT                 input file (default input_3d_centre_bore_128_a2)
#   BUILD_CPU_IF_MISSING  1 = build in-job if no CPU binary found (default 1)
#   PARTITION             nova | scavenger  (informational; actual partition set
#                         at sbatch time via --partition or --scavenger sweep flag)
# ===========================================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

INPUT="${INPUT:-input_3d_centre_bore_128_a2}"
BUILD_CPU_IF_MISSING="${BUILD_CPU_IF_MISSING:-1}"
ASYNC_OUT="${ASYNC_OUT:-0}"
TIME_CMD="${TIME_CMD:-/usr/bin/time -v}"
PARTITION="${SLURM_JOB_PARTITION:-${PARTITION:-nova}}"

module purge 2>/dev/null || true
module load gcc 2>/dev/null || module load gcc/12 2>/dev/null || true
module load openmpi 2>/dev/null || module load openmpi4 2>/dev/null || true

CPU_BIN="$(ls -t bin/alamo-3d*-g++ 2>/dev/null | grep -v cuda | head -1 || true)"
if [ -z "${CPU_BIN}" ]; then
  if [ "${BUILD_CPU_IF_MISSING}" != "1" ]; then
    echo "No bin/alamo-3d*-g++ CPU binary found and BUILD_CPU_IF_MISSING=${BUILD_CPU_IF_MISSING}."
    exit 1
  fi
  echo "No CPU binary -- building (./configure --comp=g++ --dim 3 --profile --get-eigen && make)"
  ./configure --comp=g++ --dim 3 --profile --get-eigen
  make -j"${SLURM_CPUS_ON_NODE:-64}"
  CPU_BIN="$(ls -t bin/alamo-3d*-g++ | grep -v cuda | head -1)"
fi
echo "CPU binary: ${CPU_BIN}  (${SLURM_NTASKS} ranks)"
echo "Partition: ${PARTITION}  Input: ${INPUT}"

OVERRIDES="amrex.async_out=${ASYNC_OUT} amrex.async_out_nfiles=0 \
           tiny_profiler.device_synchronize_around_region=1"

echo "=== launching A2 CPU (${SLURM_NTASKS} ranks) ==="
${TIME_CMD} srun --mpi=pmix -n "${SLURM_NTASKS}" "${CPU_BIN}" "${INPUT}" ${OVERRIDES}
echo "=== done ==="

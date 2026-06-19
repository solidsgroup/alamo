#!/bin/bash
#SBATCH -A brunnels
#SBATCH -J flame_gpu
#SBATCH --partition=nova
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --output=flame_gpu.%j.out
#SBATCH --error=flame_gpu.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jackplum@iastate.edu
# NOTE: GPU count / type / gres are set dynamically below via --export or by
# editing GPU_TYPE/NGPU, then re-emitted into the job by the wrapper logic.
# Because #SBATCH lines can't read variables, set the gres on the command line:
#   sbatch --gres=gpu:a100:1 --ntasks=1 --ntasks-per-node=1 benchmark/nova_flame_gpu.slurm
# or use the matching values below (defaults: 1x A100).
# ===========================================================================
# Multi-GPU (1 MPI rank per GPU) run of the Flame center-bore case on NOVA.
# Build the binaries first with benchmark/build_alamo_nova.sh.
#
# Examples:
#   sbatch --gres=gpu:a100:1 --ntasks=1 benchmark/nova_flame_gpu.slurm
#   GPU_TYPE=h200 sbatch --gres=gpu:h200:2 --ntasks=2 benchmark/nova_flame_gpu.slurm
#   MODE=fast sbatch --gres=gpu:a100:1 --ntasks=1 benchmark/nova_flame_gpu.slurm
# ===========================================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

# --- Selectors (override via --export=ALL,VAR=...) --------------------------
GPU_TYPE="${GPU_TYPE:-a100}"        # a100 (sm_80) or h200 (sm_90)
INPUT="${INPUT:-input_nova_centre_bore}"
# MODE: bench = profiled, accurate per-region holdups (device syncs, slower).
#       fast  = max-speed production: pure device arena, no profiler sync.
MODE="${MODE:-bench}"

case "${GPU_TYPE}" in
  a100) ARCH=80 ;;
  h200) ARCH=90 ;;
  *) echo "Unknown GPU_TYPE=${GPU_TYPE} (use a100 or h200)"; exit 1 ;;
esac
NRANKS="${SLURM_NTASKS:-1}"

# --- Modules (match build_alamo_nova.sh) ------------------------------------
module purge 2>/dev/null || true
module load cuda 2>/dev/null || module load cuda/12 2>/dev/null || true
module load gcc 2>/dev/null || module load gcc/12 2>/dev/null || true
module load openmpi 2>/dev/null || module load openmpi4 2>/dev/null || true

# --- Pick the matching-arch binary ------------------------------------------
GPU_BIN="$(ls -t bin/alamo-2d*cuda${ARCH}*-g++ 2>/dev/null | head -1 || true)"
if [ -z "${GPU_BIN}" ]; then
  echo "No bin/alamo-2d*cuda${ARCH}*-g++ found -- run benchmark/build_alamo_nova.sh first."
  exit 1
fi
echo "GPU binary: ${GPU_BIN}  (GPU_TYPE=${GPU_TYPE} sm_${ARCH}, ${NRANKS} rank(s)/GPU, mode=${MODE})"

# --- Runtime overrides ------------------------------------------------------
# async_out : background-thread IO so the GPU runs step N+1 during step N write
if [ "${MODE}" = "fast" ]; then
  # pure device arena (no host page-migration), no profiler syncs = fastest.
  OVERRIDES="amrex.async_out=1 amrex.async_out_nfiles=0 amrex.the_arena_is_managed=0"
else
  # managed arena + per-region device sync for accurate TinyProfiler holdups.
  OVERRIDES="amrex.async_out=1 amrex.async_out_nfiles=0 \
             amrex.the_arena_is_managed=1 \
             tiny_profiler.device_synchronize_around_region=1"
fi

# Optional Nsight Systems trace: NSYS_PREFIX="nsys profile --force-overwrite true -o flame_gpu_${SLURM_JOB_ID} --"
NSYS_PREFIX="${NSYS_PREFIX:-}"

echo "=== launching ==="
time ${NSYS_PREFIX} srun --mpi=pmix -n "${NRANKS}" --gpus-per-task=1 \
     "${GPU_BIN}" "${INPUT}" ${OVERRIDES}
echo "=== done ==="

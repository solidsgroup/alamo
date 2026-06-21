#!/bin/bash
#SBATCH -A brunnels
#SBATCH -J flame_gpu_3d_multi
#SBATCH --partition=nova
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --output=flame_gpu_3d_multi.%j.out
#SBATCH --error=flame_gpu_3d_multi.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jackplum@iastate.edu
# NOTE: GPU count / type / gres are set dynamically below via --export or by
# editing GPU_TYPE/NGPUS, then re-emitted into the job by the wrapper logic.
# Because #SBATCH lines can't read variables, set the gres + tasks on the command
# line to match NGPUS (1 MPI rank per GPU):
#   NGPUS=4 sbatch --gres=gpu:a100:4 --ntasks=4 --ntasks-per-node=4 benchmark/nova_flame_gpu_3d_multi.slurm
# or use the matching values below (defaults: 4x A100, 1 rank/GPU).
# ===========================================================================
# Multi-GPU (1 MPI rank per GPU) 3D run of the Flame case on NOVA.
# This is the multi-GPU 3D sibling of nova_flame_gpu.slurm, parameterized for N
# GPUs for the roadmap-3.4 strong/weak scaling study (problem-size x GPU-count).
# Build the 3D binaries first with benchmark/build_alamo_nova_3d.sh.
#
# Examples:
#   NGPUS=2 sbatch --gres=gpu:a100:2 --ntasks=2 benchmark/nova_flame_gpu_3d_multi.slurm
#   GPU_TYPE=h200 NGPUS=4 sbatch --gres=gpu:h200:4 --ntasks=4 benchmark/nova_flame_gpu_3d_multi.slurm
#   MODE=fast NGPUS=4 sbatch --gres=gpu:a100:4 --ntasks=4 benchmark/nova_flame_gpu_3d_multi.slurm
# ===========================================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

# --- Selectors (override via --export=ALL,VAR=...) --------------------------
GPU_TYPE="${GPU_TYPE:-a100}"        # a100 (sm_80) or h200 (sm_90)
INPUT="${INPUT:-input_3d_flame}"
# NGPUS: number of GPUs = number of MPI ranks (1 rank per GPU, roadmap 3.4).
NGPUS="${NGPUS:-4}"
# MODE: bench = profiled, accurate per-region holdups (device syncs, slower).
#       fast  = max-speed production: pure device arena, no profiler sync.
MODE="${MODE:-bench}"

case "${GPU_TYPE}" in
  a100) ARCH=80 ;;
  h200) ARCH=90 ;;
  *) echo "Unknown GPU_TYPE=${GPU_TYPE} (use a100 or h200)"; exit 1 ;;
esac
# 1 rank per GPU: prefer the SLURM-allocated task count, else fall back to NGPUS.
NRANKS="${SLURM_NTASKS:-${NGPUS}}"

# --- Modules (match build_alamo_nova_3d.sh) ---------------------------------
module purge 2>/dev/null || true
module load cuda 2>/dev/null || module load cuda/12 2>/dev/null || true
module load gcc 2>/dev/null || module load gcc/12 2>/dev/null || true
module load openmpi 2>/dev/null || module load openmpi4 2>/dev/null || true

# --- Pick the matching-arch binary ------------------------------------------
GPU_BIN="$(ls -t bin/alamo_gpu-3d*cuda${ARCH}*-g++ 2>/dev/null | head -1 || true)"
if [ -z "${GPU_BIN}" ]; then
  echo "No bin/alamo_gpu-3d*cuda${ARCH}*-g++ found -- run benchmark/build_alamo_nova_3d.sh first."
  exit 1
fi
echo "GPU binary: ${GPU_BIN}  (GPU_TYPE=${GPU_TYPE} sm_${ARCH}, ${NRANKS} rank(s) = ${NRANKS} GPU(s), 1 rank/GPU, mode=${MODE})"

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

# Optional Nsight Systems trace: NSYS_PREFIX="nsys profile --force-overwrite true -o flame_gpu_3d_multi_${SLURM_JOB_ID} --"
NSYS_PREFIX="${NSYS_PREFIX:-}"

echo "=== launching ${NRANKS} rank(s), 1 GPU/rank ==="
time ${NSYS_PREFIX} srun --mpi=pmix --ntasks="${NRANKS}" --gpus-per-task=1 \
     "${GPU_BIN}" "${INPUT}" ${OVERRIDES}
echo "=== done ==="

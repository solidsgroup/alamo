#!/bin/bash
#SBATCH -A brunnels
#SBATCH -J flame_gpu_3d
#SBATCH --partition=nova
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=72
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=flame_gpu_3d.%j.out
#SBATCH --error=flame_gpu_3d.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jackplum@iastate.edu
# NOTE: GPU count / type / gres are set dynamically below via --export or by
# editing GPU_TYPE/NGPU, then re-emitted into the job by the wrapper logic.
# Because #SBATCH lines can't read variables, set the gres on the command line:
#   sbatch --gres=gpu:v100:1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=72 benchmark/nova_flame_gpu_3d.slurm
# ==========================================================================
# Single-GPU (1 MPI rank per GPU) 3D run of the Flame center-bore case on NOVA.
# Build the 3D binaries first with benchmark/build_alamo_nova_3d.sh.
# This is the 3D sibling of nova_flame_gpu.slurm (3D binary glob + center-bore inputs).
#
# Examples:
#   sbatch --gres=gpu:v100:1 --ntasks=1 --cpus-per-task=32 benchmark/nova_flame_gpu_3d.slurm
#   GPU_TYPE=v100 MODE=fast sbatch --gres=gpu:v100:1 --ntasks=1 --cpus-per-task=32 benchmark/nova_flame_gpu_3d.slurm
#   sbatch --partition=scavenger --gres=gpu:a100:1 --ntasks=1 --cpus-per-task=32 benchmark/nova_flame_gpu_3d.slurm
# ==========================================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

GPU_TYPE="${GPU_TYPE:-v100}"        # v100 (sm_70), a100 (sm_80) or h200 (sm_90)
INPUT="${INPUT:-input_3d_centre_bore_128}"
MODE="${MODE:-bench}"
ASYNC_OUT="${ASYNC_OUT:-0}"
PROFILE_KIND="${PROFILE_KIND:-summary}"
TIME_CMD="${TIME_CMD:-/usr/bin/time -v}"

case "${GPU_TYPE}" in
  v100) ARCH=70 ;;
  a100) ARCH=80 ;;
  h200) ARCH=90 ;;
  *) echo "Unknown GPU_TYPE=${GPU_TYPE} (use v100, a100 or h200)"; exit 1 ;;
esac
NRANKS="${SLURM_NTASKS:-1}"

module purge 2>/dev/null || true
module load cuda 2>/dev/null || module load cuda/12 2>/dev/null || true
module load gcc 2>/dev/null || module load gcc/12 2>/dev/null || true
module load openmpi 2>/dev/null || module load openmpi4 2>/dev/null || true

GPU_BIN="$(ls -t bin/alamo_gpu-3d*cuda${ARCH}*-g++ 2>/dev/null | head -1 || true)"
if [ -z "${GPU_BIN}" ]; then
  echo "No bin/alamo_gpu-3d*cuda${ARCH}*-g++ found -- run benchmark/build_alamo_nova_3d.sh first."
  exit 1
fi
echo "GPU binary: ${GPU_BIN}  (GPU_TYPE=${GPU_TYPE} sm_${ARCH}, ${NRANKS} rank(s)/GPU, mode=${MODE})"

if [ "${MODE}" = "fast" ]; then
  OVERRIDES="amrex.async_out=${ASYNC_OUT} amrex.async_out_nfiles=0 amrex.the_arena_is_managed=0"
else
  OVERRIDES="amrex.async_out=${ASYNC_OUT} amrex.async_out_nfiles=0              amrex.the_arena_is_managed=1              tiny_profiler.device_synchronize_around_region=1"
fi

if [ -z "${NSYS_PREFIX:-}" ]; then
  if [ "${PROFILE_KIND}" = "nsys" ]; then
    NSYS_PREFIX="nsys profile --force-overwrite true --trace=cuda,nvtx,osrt,mpi --sample=none --cpuctxsw=none -o flame_gpu_3d_${SLURM_JOB_ID}"
  else
    NSYS_PREFIX=""
  fi
fi

echo "=== launching ==="
${TIME_CMD} ${NSYS_PREFIX} srun --mpi=pmix -n "${NRANKS}" --gpus-per-task=1      "${GPU_BIN}" "${INPUT}" ${OVERRIDES}
echo "=== done ==="

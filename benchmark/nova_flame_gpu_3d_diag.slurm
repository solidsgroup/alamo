#!/bin/bash
#SBATCH -A brunnels
#SBATCH -J flame_gpu_3d_diag
#SBATCH --partition=nova
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
#SBATCH --output=flame_gpu_3d_diag.%j.out
#SBATCH --error=flame_gpu_3d_diag.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jackplum@iastate.edu
#
# Single-GPU diagnostic run for Phase-3 CUDA error 700 failures.
#
# This intentionally runs the smallest 3D Flame case for one step with
# CUDA_LAUNCH_BLOCKING enabled and writes all runtime artifacts into a per-job
# directory. Use benchmark/phase3_nova_diag.sh to submit the common variants.
#
# Example:
#   DIAG_TOOL=blocking INPUT=input_3d_flame_128 GPU_TYPE=a100 \
#     sbatch --gres=gpu:a100:1 --ntasks=1 benchmark/nova_flame_gpu_3d_diag.slurm
#   DIAG_TOOL=sanitizer INPUT=input_3d_flame_128 GPU_TYPE=a100 \
#     sbatch --gres=gpu:a100:1 --ntasks=1 benchmark/nova_flame_gpu_3d_diag.slurm

set -euo pipefail

ROOT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${ROOT_DIR}"

GPU_TYPE="${GPU_TYPE:-a100}"          # a100 (sm_80) or h200 (sm_90)
INPUT="${INPUT:-input_3d_flame_128}"
DIAG_TOOL="${DIAG_TOOL:-blocking}"    # blocking | sanitizer
MANAGED_ARENA="${MANAGED_ARENA:-1}"
ASYNC_OUT="${ASYNC_OUT:-0}"

case "${GPU_TYPE}" in
  a100) ARCH=80 ;;
  h200) ARCH=90 ;;
  *) echo "Unknown GPU_TYPE=${GPU_TYPE} (use a100 or h200)"; exit 1 ;;
esac

module purge 2>/dev/null || true
module load cuda 2>/dev/null || module load cuda/12 2>/dev/null || true
module load gcc 2>/dev/null || module load gcc/12 2>/dev/null || true
module load openmpi 2>/dev/null || module load openmpi4 2>/dev/null || true

GPU_BIN="$(ls -t bin/alamo_gpu-3d*cuda${ARCH}*-g++ 2>/dev/null | head -1 || true)"
if [ -z "${GPU_BIN}" ]; then
  echo "No bin/alamo_gpu-3d*cuda${ARCH}*-g++ found -- run benchmark/build_alamo_nova_3d.sh first."
  exit 1
fi
GPU_BIN="$(cd "$(dirname "${GPU_BIN}")" && pwd)/$(basename "${GPU_BIN}")"

if [[ "${INPUT}" = /* ]]; then
  INPUT_PATH="${INPUT}"
else
  INPUT_PATH="${ROOT_DIR}/${INPUT}"
fi
if [ ! -f "${INPUT_PATH}" ]; then
  echo "Input file not found: ${INPUT_PATH}"
  exit 1
fi

safe_input="$(basename "${INPUT_PATH}" | tr -c 'A-Za-z0-9_.-' '_')"
RUN_DIR="${DIAG_DIR:-${ROOT_DIR}/phase3_diag_${SLURM_JOB_ID:-manual}_${DIAG_TOOL}_${safe_input}}"
mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"

finish() {
  local status=$?
  echo
  echo "=== diagnostic artifact directory ==="
  echo "${RUN_DIR}"
  echo
  echo "=== top-level artifacts ==="
  ls -lah || true
  echo
  for bt in Backtrace.*; do
    [ -e "${bt}" ] || continue
    echo "=== ${bt} ==="
    sed -n '1,220p' "${bt}" || true
  done
  exit "${status}"
}
trap finish EXIT

export CUDA_LAUNCH_BLOCKING="${CUDA_LAUNCH_BLOCKING:-1}"
export CUDA_DEVICE_WAITS_ON_EXCEPTION="${CUDA_DEVICE_WAITS_ON_EXCEPTION:-0}"

OVERRIDES="max_step=1 stop_time=1e-12_s plot_file=output_diag \
           amr.plot_int=-1 amr.thermo.plot_int=-1 \
           amrex.async_out=${ASYNC_OUT} amrex.async_out_nfiles=0 \
           amrex.the_arena_is_managed=${MANAGED_ARENA} \
           tiny_profiler.device_synchronize_around_region=1"

echo "=== diagnostic configuration ==="
echo "repo          = ${ROOT_DIR}"
echo "run_dir       = ${RUN_DIR}"
echo "gpu_bin       = ${GPU_BIN}"
echo "input         = ${INPUT_PATH}"
echo "gpu_type      = ${GPU_TYPE} sm_${ARCH}"
echo "diag_tool     = ${DIAG_TOOL}"
echo "managed_arena = ${MANAGED_ARENA}"
echo "CUDA_LAUNCH_BLOCKING=${CUDA_LAUNCH_BLOCKING}"
echo "overrides     = ${OVERRIDES}"
echo

case "${DIAG_TOOL}" in
  blocking)
    echo "=== launching with CUDA_LAUNCH_BLOCKING ==="
    time srun --mpi=pmix -n 1 --gpus-per-task=1 \
      "${GPU_BIN}" "${INPUT_PATH}" ${OVERRIDES}
    ;;
  sanitizer)
    SANITIZER="$(command -v compute-sanitizer || true)"
    if [ -z "${SANITIZER}" ]; then
      SANITIZER="$(command -v cuda-memcheck || true)"
    fi
    if [ -z "${SANITIZER}" ]; then
      echo "No compute-sanitizer or cuda-memcheck found in PATH after loading CUDA modules."
      exit 1
    fi
    echo "=== launching under ${SANITIZER} ==="
    time srun --mpi=pmix -n 1 --gpus-per-task=1 \
      "${SANITIZER}" --tool memcheck --error-exitcode 99 \
      "${GPU_BIN}" "${INPUT_PATH}" ${OVERRIDES}
    ;;
  *)
    echo "Unknown DIAG_TOOL=${DIAG_TOOL} (use blocking or sanitizer)"
    exit 2
    ;;
esac

echo "=== done ==="

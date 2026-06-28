#!/bin/bash
#SBATCH -A brunnels
#SBATCH -J flame_gpu_3d_diag
#SBATCH --partition=nova
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --output=flame_gpu_3d_diag.%j.out
#SBATCH --error=flame_gpu_3d_diag.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jackplum@iastate.edu
#
# Single-GPU diagnostic / profiling run.  DIAG_TOOL selects the mode:
#
#   blocking  -- CUDA_LAUNCH_BLOCKING=1, 1 step  (correctness / CUDA error location)
#   sanitizer -- compute-sanitizer memcheck,  1 step
#   ncu       -- Nsight Compute kernel profiles, 2 steps, elastic enabled
#                outputs to benchmark/profiles_alpha1/ncu_<jobid>/
#   nsys      -- Nsight Systems timeline, 10 steps, elastic enabled
#                outputs to benchmark/profiles_alpha1/nsys_<jobid>/
#
# Example (correctness, nova partition):
#   DIAG_TOOL=blocking INPUT=input_3d_centre_bore_128 GPU_TYPE=v100 sbatch --gres=gpu:v100:1 --ntasks=1 --cpus-per-task=32 benchmark/nova_flame_gpu_3d_diag.slurm
#
# Example (A1 profiling, A100, nova partition):
#   DIAG_TOOL=ncu  INPUT=input_3d_flame_128 GPU_TYPE=a100 sbatch --gres=gpu:a100:1 --ntasks=1 --cpus-per-task=32 benchmark/nova_flame_gpu_3d_diag.slurm
#   DIAG_TOOL=nsys INPUT=input_3d_flame_128 GPU_TYPE=a100 sbatch --gres=gpu:a100:1 --ntasks=1 --cpus-per-task=32 benchmark/nova_flame_gpu_3d_diag.slurm
#
# Example (scavenger partition -- better queue times, preemptible):
#   DIAG_TOOL=nsys INPUT=input_3d_flame_128 GPU_TYPE=a100 sbatch --partition=scavenger --gres=gpu:a100:1 --ntasks=1 --cpus-per-task=32 benchmark/nova_flame_gpu_3d_diag.slurm

set -euo pipefail

ROOT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${ROOT_DIR}"

GPU_TYPE="${GPU_TYPE:-v100}"          # v100 (sm_70), a100 (sm_80) or h200 (sm_90)
INPUT="${INPUT:-input_3d_centre_bore_128}"
DIAG_TOOL="${DIAG_TOOL:-blocking}"    # blocking | sanitizer | ncu | nsys
PROFILE_DIR="${PROFILE_DIR:-${ROOT_DIR}/benchmark/profiles_alpha1}"
MANAGED_ARENA="${MANAGED_ARENA:-1}"
ASYNC_OUT="${ASYNC_OUT:-0}"

case "${GPU_TYPE}" in
  v100) ARCH=70 ;;
  a100) ARCH=80 ;;
  h200) ARCH=90 ;;
  *) echo "Unknown GPU_TYPE=${GPU_TYPE} (use v100, a100 or h200)"; exit 1 ;;
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

export CUDA_LAUNCH_BLOCKING="${CUDA_LAUNCH_BLOCKING:-0}"
export CUDA_DEVICE_WAITS_ON_EXCEPTION="${CUDA_DEVICE_WAITS_ON_EXCEPTION:-0}"

OVERRIDES_SHORT="max_step=1 stop_time=1e-12_s plot_file=output_diag \
                 amr.plot_int=-1 amr.thermo.plot_int=-1 \
                 amrex.async_out=${ASYNC_OUT} amrex.async_out_nfiles=0 \
                 amrex.the_arena_is_managed=${MANAGED_ARENA} \
                 tiny_profiler.device_synchronize_around_region=1"

# ncu/nsys profiling: use the elastic-enabled input variant so no CLI overrides
# are needed to initialize the elastic solver.  The _elastic input has elastic.interval=1
# and elastic.tstart=0.0 so the solve fires on step 1.
PROFILE_INPUT="${ROOT_DIR}/${INPUT/_128/_128_elastic}"
if [ ! -f "${PROFILE_INPUT}" ]; then
    PROFILE_INPUT="${INPUT_PATH}"
    echo "WARN: elastic input ${PROFILE_INPUT} not found; falling back to ${INPUT_PATH} (elastic may not run)"
fi

OVERRIDES_PROFILE="max_step=2 stop_time=1e99 amr.plot_int=-1 amr.thermo.plot_int=-1 \
                   amrex.async_out=0 amrex.the_arena_is_managed=${MANAGED_ARENA} \
                   tiny_profiler.device_synchronize_around_region=0 \
                   elastic.solver.verbose=0"

OVERRIDES_NSYS="max_step=10 stop_time=1e99 amr.plot_int=-1 amr.thermo.plot_int=-1 \
                amrex.async_out=0 amrex.the_arena_is_managed=${MANAGED_ARENA} \
                tiny_profiler.device_synchronize_around_region=0 \
                elastic.solver.verbose=0"

echo "=== diagnostic configuration ==="
echo "repo          = ${ROOT_DIR}"
echo "run_dir       = ${RUN_DIR}"
echo "gpu_bin       = ${GPU_BIN}"
echo "input         = ${INPUT_PATH}"
echo "gpu_type      = ${GPU_TYPE} sm_${ARCH}"
echo "diag_tool     = ${DIAG_TOOL}"
echo "managed_arena = ${MANAGED_ARENA}"
echo "profile_dir   = ${PROFILE_DIR}"
echo

case "${DIAG_TOOL}" in
  blocking)
    export CUDA_LAUNCH_BLOCKING=1
    echo "=== launching with CUDA_LAUNCH_BLOCKING ==="
    time srun --mpi=pmix -n 1 --gpus-per-task=1 \
      "${GPU_BIN}" "${INPUT_PATH}" ${OVERRIDES_SHORT}
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
      "${GPU_BIN}" "${INPUT_PATH}" ${OVERRIDES_SHORT}
    ;;

  ncu)
    NCU="$(command -v ncu 2>/dev/null || true)"
    if [ -z "${NCU}" ]; then
      echo "ncu not found in PATH -- load the nsight-compute module or check PATH" >&2
      exit 1
    fi
    echo "ncu = ${NCU}  ($(${NCU} --version 2>&1 | head -1))"
    mkdir -p "${PROFILE_DIR}"
    OUT_BASE="${PROFILE_DIR}/ncu_${SLURM_JOB_ID:-manual}"
    mkdir -p "${OUT_BASE}"
    echo "=== ncu profiles will be written to ${OUT_BASE}/ ==="
    echo "=== profiling input: ${PROFILE_INPUT} ==="

    # AMReX launches every GPU kernel through one wrapper (amrex::launch_global),
    # so ncu only ever sees the function name "launch_global" -- matching on
    # --kernel-name-base function never selects the physics kernels (that bug left
    # every ncu_*/ empty while the job still exited 0). Fix: don't match by name at
    # all -- profile a bounded WINDOW of launches (--launch-skip/--launch-count) and
    # read the *demangled* names in post; they carry the enclosing function (Fapply,
    # Fsmooth, interpolation, ...). (NVTX --nvtx-include was tried and never matched
    # AMReX's deeply-nested push/pop ranges -- "No kernels were profiled".)
    #
    # Metric set + section files. Two gotchas on this ncu (Nsight Compute 2025.1.1):
    #   - there is NO set named "default" anymore (sets are basic/detailed/full/
    #     roofline/...). "--set default" enables nothing -> "No metrics to collect
    #     found in sections" -> "No kernels were profiled". Use basic; override with
    #     NCU_SET=full (roofline + full memory workload) or NCU_SET=detailed.
    #   - pass --section-folder (the install's sections dir) in case ncu can't find
    #     it on its own; fall back to explicit --metrics if even that is missing.
    NCU_DIR="$(dirname "$(readlink -f "${NCU}")")"
    SECTION_DIR="$(find "${NCU_DIR}/.." -maxdepth 4 -name 'SpeedOfLight.section' -printf '%h\n' 2>/dev/null | head -1 || true)"
    if [ -n "${SECTION_DIR}" ]; then
      METRIC_ARGS=(--set "${NCU_SET:-basic}" --section-folder "${SECTION_DIR}")
      echo "ncu sections  = ${SECTION_DIR}"
    else
      METRIC_ARGS=(--metrics "${NCU_METRICS:-gpu__time_duration.sum,sm__throughput.avg.pct_of_peak_sustained_elapsed,gpu__compute_memory_throughput.avg.pct_of_peak_sustained_elapsed,sm__warps_active.avg.pct_of_peak_sustained_active,launch__registers_per_thread,launch__waves_per_multiprocessor}")
      echo "ncu sections  = NOT FOUND near ${NCU_DIR} -- using explicit --metrics"
    fi
    echo "=== ncu --list-sets ==="
    "${NCU}" --list-sets ${SECTION_DIR:+--section-folder "${SECTION_DIR}"} 2>&1 | head -20 || true

    # One pass: profile launches [skip, skip+count) with no name/NVTX filter. The
    # elastic solve dominates the stream, so a small window catches Fapply (+ the
    # other hot kernels). Identify them by demangled name in the summary below.
    NCU_SKIP="${NCU_SKIP:-0}"
    NCU_COUNT="${NCU_COUNT:-60}"
    OUT_REP="${OUT_BASE}/elastic_kernels"
    echo "--- ncu: profiling launches [skip ${NCU_SKIP}, count ${NCU_COUNT}], set=${NCU_SET:-basic} ---"
    srun --mpi=pmix -n 1 --gpus-per-task=1 \
      "${NCU}" \
        --target-processes all \
        "${METRIC_ARGS[@]}" \
        --launch-skip "${NCU_SKIP}" \
        --launch-count "${NCU_COUNT}" \
        --print-kernel-base demangled \
        --export "${OUT_REP}" \
        --force-overwrite \
        "${GPU_BIN}" "${PROFILE_INPUT}" ${OVERRIDES_PROFILE} || true

    echo "=== ncu complete; profiles in ${OUT_BASE}/ ==="
    ls -lh "${OUT_BASE}/" || true
    if [ ! -s "${OUT_REP}.ncu-rep" ]; then
      echo "=== ncu: NO report written (${OUT_REP}.ncu-rep) -- see warnings above ===" >&2
      exit 1
    fi
    echo "=== per-kernel SoL summary (grep Fapply / Fsmooth / interpolation) ==="
    "${NCU}" -i "${OUT_REP}.ncu-rep" --csv --page raw --print-kernel-base demangled \
      --metrics gpu__time_duration.sum,sm__throughput.avg.pct_of_peak_sustained_elapsed,gpu__compute_memory_throughput.avg.pct_of_peak_sustained_elapsed,sm__warps_active.avg.pct_of_peak_sustained_active \
      2>/dev/null | head -120 || true
    ;;

  nsys)
    NSYS="$(command -v nsys 2>/dev/null || true)"
    if [ -z "${NSYS}" ]; then
      echo "nsys not found in PATH -- load the nsight-systems module or check PATH" >&2
      exit 1
    fi
    echo "nsys = ${NSYS}  ($(${NSYS} --version 2>&1 | head -1))"
    mkdir -p "${PROFILE_DIR}"
    OUT_BASE="${PROFILE_DIR}/nsys_${SLURM_JOB_ID:-manual}"
    mkdir -p "${OUT_BASE}"
    echo "=== nsys report will be written to ${OUT_BASE}/ ==="

    echo "=== profiling input: ${PROFILE_INPUT} ==="
    # --trace=cuda only: osrt requires nsys process injection which conflicts with
    # OpenMPI's UCX plugin dlopen (boost raises "Expected shared object name, found
    # a path delimiter" and aborts). CUDA trace is sufficient for A1 metrics.
    # OMPI_MCA_opal_common_ucx_opal_mem_hooks=0: belt-and-suspenders UCX hook disable.
    OMPI_MCA_opal_common_ucx_opal_mem_hooks=0 \
    srun --mpi=pmix -n 1 --gpus-per-task=1 \
      "${NSYS}" profile \
        --trace=cuda,nvtx \
        --cuda-memory-usage=false \
        --output "${OUT_BASE}/nsys_report" \
        --force-overwrite true \
        "${GPU_BIN}" "${PROFILE_INPUT}" ${OVERRIDES_NSYS}

    echo "=== nsys stats summary ==="
    "${NSYS}" stats \
      --report cuda_api_sum,cuda_gpu_kern_sum,osrt_sum \
      --format csv \
      --output "${OUT_BASE}/nsys_stats" \
      "${OUT_BASE}/nsys_report.nsys-rep" 2>/dev/null || \
    "${NSYS}" stats \
      --report cuda_api_sum,cuda_gpu_kern_sum \
      --format csv \
      "${OUT_BASE}/nsys_report.nsys-rep" | tee "${OUT_BASE}/nsys_stats_stdout.csv" || true

    echo "=== nsys complete; report in ${OUT_BASE}/ ==="
    ls -lh "${OUT_BASE}/" || true
    ;;

  *)
    echo "Unknown DIAG_TOOL=${DIAG_TOOL} (use blocking | sanitizer | ncu | nsys)"
    exit 2
    ;;
esac

echo "=== done ==="

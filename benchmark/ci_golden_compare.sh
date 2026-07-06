#!/usr/bin/env bash
# CI correctness gate for chamber-gpu.
#
# Two tripwires, re-armed on every push to chamber-gpu:
#   1. Golden compare: runs benchmark/baseline_suite.py in "check" mode against
#      benchmark/baseline_references/{canonical_step1,canonical_step2,
#      eta_expression_step1} and fails the build on any mismatch.
#   2. NaN-flag assertion smoke: a short Flame run. Flame.cpp's advance kernel
#      sets a Util::DeviceErrorFlag (Util::SetDeviceError) on any NaN/Inf in K,
#      rho, cp, L, eta, alpha, mdot, heatflux, or temperature, and
#      Util::AbortIfDeviceError aborts the run (nonzero exit) if that flag is
#      ever set. A clean exit therefore proves the tripwire did not fire.
#
# Mode is selected with GOLDEN_MODE=cpu|gpu (default cpu). The CPU mode runs on
# any standard runner with no GPU. The GPU mode additionally requires a CUDA
# toolchain/runner and is expected to be invoked only from a runner-gated job
# (see .github/workflows/chamber-gpu-correctness.yml); it is not exercised on
# standard CI runners.
#
# Usage:
#   GOLDEN_MODE=cpu  benchmark/ci_golden_compare.sh   # default, no GPU needed
#   GOLDEN_MODE=gpu  benchmark/ci_golden_compare.sh   # requires nvcc + GPU
set -euo pipefail
cd "$(dirname "$0")/.." || exit 1

GOLDEN_MODE="${GOLDEN_MODE:-cpu}"
DIM="${DIM:-2}"
COMP="${COMP:-g++}"
BUILD_JOBS="${BUILD_JOBS:-$(nproc 2>/dev/null || echo 4)}"
NAN_SMOKE_INPUT="${NAN_SMOKE_INPUT:-input}"
NAN_SMOKE_MAX_STEP="${NAN_SMOKE_MAX_STEP:-2}"

echo "=== ci_golden_compare.sh: mode=${GOLDEN_MODE} dim=${DIM} comp=${COMP} ==="

run_nan_flag_smoke() {
  local binary="$1"
  local out_dir="benchmark/ci_nan_smoke_${GOLDEN_MODE}"
  rm -rf "${out_dir}"
  mkdir -p "${out_dir}"

  echo "--- NaN-flag assertion smoke (${GOLDEN_MODE}): ${binary}"
  # A short run that completes with exit 0 demonstrates that
  # Util::AbortIfDeviceError never observed a set device-error flag, i.e. the
  # NaN/Inf tripwires in Flame.cpp's advance kernel did not fire. If the flag
  # were ever set, the run aborts with a nonzero exit code and this script
  # fails.
  if ! mpiexec -np 1 "${binary}" "${NAN_SMOKE_INPUT}" \
      "max_step=${NAN_SMOKE_MAX_STEP}" \
      "stop_time=1e99_s" \
      "amr.plot_int=-1" \
      "amr.thermo.plot_int=1" \
      "amr.thermo.int=1" \
      "elastic.solver.verbose=0" \
      "elastic.print_model=0" \
      "plot_file=${out_dir}/plot" \
      >"${out_dir}/run.log" 2>&1; then
    echo "NaN-flag assertion smoke FAILED: ${binary} aborted; see ${out_dir}/run.log" >&2
    tail -n 80 "${out_dir}/run.log" >&2 || true
    return 1
  fi
  echo "NaN-flag assertion smoke OK: ${binary} completed without tripping the device-error flag"
}

if [ "${GOLDEN_MODE}" = "cpu" ]; then
  # CPU-resident correctness leg. No fast-math tricks are in play on the CPU
  # build, so a plain g++ build is the correctness baseline for this leg.
  echo "--- configure (CPU, dim=${DIM}, comp=${COMP})"
  ./configure --dim="${DIM}" --comp="${COMP}"
  echo "--- make -j${BUILD_JOBS}"
  make -j"${BUILD_JOBS}"

  # Select the binary this configure/make just produced by its deterministic
  # name -- a wildcard sort here once picked a stale alamo-2d-perf-clang++
  # over the freshly built alamo-2d-g++ and gated against week-old code.
  CPU_BIN="${CPU_BIN:-bin/alamo-${DIM}d-${COMP}}"
  if [ ! -x "${CPU_BIN}" ]; then
    echo "Expected CPU binary ${CPU_BIN} missing after build" >&2
    exit 2
  fi
  echo "cpu binary: ${CPU_BIN}"

  export CPU_BIN
  echo "--- golden compare (baseline_suite.py check, profiles=cpu)"
  python3 benchmark/baseline_suite.py check --profiles=cpu

  run_nan_flag_smoke "${CPU_BIN}"

elif [ "${GOLDEN_MODE}" = "gpu" ]; then
  # Device-resident correctness leg: no-fast-math (--cuda-fp strict) CUDA
  # build, exercising the bit/tolerance-strict golden compare against the
  # gpu_strict baseline references. Requires nvcc and a CUDA-capable GPU; only
  # intended to run from a runner that actually has both (see the
  # golden-gpu job's `if:` guard in chamber-gpu-correctness.yml).
  if ! command -v nvcc >/dev/null 2>&1; then
    echo "GOLDEN_MODE=gpu requested but nvcc is not on PATH; this leg must run on a CUDA runner" >&2
    exit 2
  fi

  ARCH="${ARCH:-$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '. ')}"
  if [ -z "${ARCH}" ]; then
    echo "Could not detect CUDA compute capability; set ARCH explicitly" >&2
    exit 2
  fi

  echo "--- configure (GPU, dim=${DIM}, comp=${COMP}, cuda=${ARCH}, cuda-fp=strict)"
  ./configure --dim="${DIM}" --comp="${COMP}" --cuda "${ARCH}" --cuda-fp strict
  echo "--- make -j${BUILD_JOBS}"
  make -j"${BUILD_JOBS}"

  # Deterministic just-built name (see CPU-leg comment on stale-binary risk).
  GPU_STRICT_BIN="${GPU_STRICT_BIN:-bin/alamo_gpu-${DIM}d-nofast-cuda${ARCH}-${COMP}}"
  if [ ! -x "${GPU_STRICT_BIN}" ]; then
    echo "Expected no-fast-math GPU binary ${GPU_STRICT_BIN} missing after build" >&2
    exit 2
  fi
  echo "gpu_strict binary: ${GPU_STRICT_BIN}"

  export GPU_STRICT_BIN
  echo "--- golden compare (baseline_suite.py check, profiles=gpu_strict)"
  python3 benchmark/baseline_suite.py check --profiles=gpu_strict

  run_nan_flag_smoke "${GPU_STRICT_BIN}"

else
  echo "Unknown GOLDEN_MODE='${GOLDEN_MODE}' (expected 'cpu' or 'gpu')" >&2
  exit 2
fi

echo "=== ci_golden_compare.sh (${GOLDEN_MODE}): PASS ==="

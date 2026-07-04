#!/usr/bin/env bash
# Run a short Flame case on CPU and GPU and compare thermo.dat scalars.
set -euo pipefail
cd "$(dirname "$0")/.." || exit 1

INPUT="${1:-input}"
MAXSTEP="${2:-1}"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT="${OUT:-benchmark/golden_${STAMP}}"
NP="${NP:-1}"
ABS_TOL="${ABS_TOL:-1e-8}"
REL_TOL="${REL_TOL:-1e-6}"

mkdir -p "${OUT}"

if [ -z "${CPU_BIN:-}" ]; then
  CPU_BIN="$(find bin -maxdepth 1 -type f -executable -name 'alamo-2d-*' ! -name '*cuda*' | sort | tail -1)"
fi

if [ -z "${GPU_BIN:-}" ]; then
  LOCAL_ARCH="${CUDA_ARCH:-}"
  if [ -z "${LOCAL_ARCH}" ] && command -v nvidia-smi >/dev/null 2>&1; then
    LOCAL_ARCH="$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '. ' || true)"
  fi
  if [ -n "${LOCAL_ARCH}" ]; then
    GPU_BIN="$(find bin -maxdepth 1 -type f -executable -name "alamo_gpu-2d*cuda${LOCAL_ARCH}*" | sort | tail -1)"
  fi
  if [ -z "${GPU_BIN:-}" ]; then
    GPU_BIN="$(find bin -maxdepth 1 -type f -executable -name 'alamo_gpu-2d*cuda*' | sort | tail -1)"
  fi
fi

if [ -z "${CPU_BIN:-}" ] || [ -z "${GPU_BIN:-}" ]; then
  echo "Missing binary. CPU_BIN='${CPU_BIN:-}' GPU_BIN='${GPU_BIN:-}'" >&2
  exit 2
fi

COMMON_OVERRIDES=(
  "max_step=${MAXSTEP}"
  "stop_time=1e99_s"
  "amr.plot_int=-1"
  "amr.thermo.plot_int=1"
  "amr.thermo.int=1"
  "elastic.solver.verbose=0"
  "elastic.print_model=0"
)

run_one() {
  local label="$1"
  local binary="$2"
  local plot="${OUT}/out_${label}"
  local log="${OUT}/${label}.log"

  echo "--- ${label}: ${binary}"
  mpiexec -np "${NP}" "${binary}" "${INPUT}" "${COMMON_OVERRIDES[@]}" "plot_file=${plot}" >"${log}" 2>&1
  if [ ! -s "${plot}/thermo.dat" ]; then
    echo "${label} did not write ${plot}/thermo.dat; see ${log}" >&2
    exit 3
  fi
}

echo "=== Flame CPU/GPU golden compare ==="
echo "input=${INPUT} max_step=${MAXSTEP} out=${OUT}"
echo "cpu=${CPU_BIN}"
echo "gpu=${GPU_BIN}"

run_one cpu "${CPU_BIN}"
run_one gpu "${GPU_BIN}"

python3 benchmark/compare_thermo.py \
  "${OUT}/out_cpu/thermo.dat" \
  "${OUT}/out_gpu/thermo.dat" \
  --abs-tol "${ABS_TOL}" \
  --rel-tol "${REL_TOL}" | tee "${OUT}/thermo_compare.txt"

echo "Artifacts in ${OUT}"

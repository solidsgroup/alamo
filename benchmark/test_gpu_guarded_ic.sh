#!/usr/bin/env bash
# Verify guarded host-loop IC paths abort clearly in CUDA builds.
set -euo pipefail
cd "$(dirname "$0")/.." || exit 1

INPUT="${1:-input}"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT="${OUT:-benchmark/guarded_ic_${STAMP}}"
mkdir -p "${OUT}"

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

if [ -z "${GPU_BIN:-}" ]; then
  echo "Missing GPU binary. Set GPU_BIN=bin/alamo_gpu-2d-..." >&2
  exit 2
fi

LOG="${OUT}/guarded_ic.log"
OVERRIDES=(
  "max_step=1"
  "stop_time=1e-12_s"
  "plot_file=${OUT}/out_guarded_ic"
  "allow_unused=1"
  "amr.plot_int=-1"
  "amr.thermo.plot_int=-1"
  "pf.eta.ic.type=laminate"
  "pf.eta.ic.laminate.thickness=0.01_m"
  "pf.eta.ic.laminate.orientation=1 0"
  "pf.eta.ic.laminate.eps=1e-5_m"
  "pf.eta.ic.laminate.mollifier=gaussian"
  "pf.eta.ic.laminate.singlefab=1"
  "pf.eta.ic.laminate.invert=0"
  "elastic.solver.verbose=0"
  "elastic.print_model=0"
)

echo "=== CUDA guarded IC negative test ==="
echo "input=${INPUT} out=${OUT}"
echo "gpu=${GPU_BIN}"

set +e
mpiexec -np 1 "${GPU_BIN}" "${INPUT}" "${OVERRIDES[@]}" >"${LOG}" 2>&1
status=$?
set -e

if [ "${status}" -eq 0 ]; then
  echo "Expected guarded IC abort, but run succeeded. See ${LOG}" >&2
  exit 1
fi

if ! grep -q "IC::Laminate is not supported in CUDA builds" "${LOG}"; then
  echo "Run failed, but not with the expected guarded-IC message. See ${LOG}" >&2
  exit 1
fi

echo "Guarded IC abort observed as expected. Log: ${LOG}"

#!/usr/bin/env bash
# Build a local-workstation CUDA binary for quick chamber-gpu smoke tests.
# The NOVA helper still builds sm_80/sm_90 for A100/H200; this script targets
# whatever GPU is installed on the current machine, e.g. sm_86 for an RTX A1000.
set -euo pipefail
cd "$(dirname "$0")/.."

LOCAL_CUDA="${LOCAL_CUDA:-$(pwd)/.local/cuda}"
if ! command -v nvcc >/dev/null 2>&1 && [ -x "${LOCAL_CUDA}/bin/nvcc" ]; then
  export CUDA_HOME="${LOCAL_CUDA}"
  export CUDA_PATH="${LOCAL_CUDA}"
  export PATH="${CUDA_HOME}/bin:${PATH}"
  export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${CUDA_HOME}/lib:${CUDA_HOME}/lib/x86_64-linux-gnu:${CUDA_HOME}/lib/cuda/lib64:${LD_LIBRARY_PATH:-}"
export LIBRARY_PATH="${CUDA_HOME}/lib/stubs:${CUDA_HOME}/lib64/stubs:${LIBRARY_PATH:-}"
fi

COMP="${COMP:-g++}"
DIM="${DIM:-2}"
PROFILE="${PROFILE:-0}"
CUDA_FP="${CUDA_FP:-fast}"
BUILD_JOBS="${BUILD_JOBS:-$(nproc 2>/dev/null || echo 8)}"
INPUT="${INPUT:-input_nova_centre_bore}"
SMOKE="${SMOKE:-1}"
MAX_STEP="${MAX_STEP:-1}"
OUTPUT_BINARY="${OUTPUT_BINARY:-}"
SOURCE_BINARY="${SOURCE_BINARY:-}"
if [ -z "${SOURCE_BINARY}" ]; then echo "SOURCE_BINARY is required" >&2; exit 1; fi
BUILD_MARKER=""
for p in "${SOURCE_BINARY}" "${OUTPUT_BINARY}"; do
  if [ -n "${p}" ]; then case "${p}" in /*) ;; *) echo "binary paths must be absolute" >&2; exit 1;; esac; fi
done
if [ -n "${OUTPUT_BINARY}" ] && [ -e "${OUTPUT_BINARY}" ]; then echo "OUTPUT_BINARY already exists" >&2; exit 1; fi
if [ -n "${SOURCE_BINARY}" ]; then BUILD_MARKER="$(mktemp)"; fi
PRE_BUILD_STAT="$(stat -c '%s:%Y:%i' "${SOURCE_BINARY}" 2>/dev/null || true)"
trap 'rm -f "${BUILD_MARKER}"' EXIT

if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "nvidia-smi not found; pass ARCH explicitly, e.g. ARCH=86 $0" >&2
  exit 1
fi

ARCH="${ARCH:-$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')}"
if [ -z "${ARCH}" ]; then
  echo "Could not detect CUDA compute capability; pass ARCH explicitly, e.g. ARCH=86 $0" >&2
  exit 1
fi

CONFIG_FLAGS=(--comp="${COMP}" --dim "${DIM}" --cuda "${ARCH}" --cuda-fp "${CUDA_FP}")
if [ "${PROFILE}" = "1" ]; then
  CONFIG_FLAGS+=(--profile)
fi

echo "=== local CUDA build ==="
echo "gpu:   $(nvidia-smi --query-gpu=name,compute_cap --format=csv,noheader | head -1)"
echo "arch:  sm_${ARCH}"
echo "comp:  ${COMP}"
echo "fp:    ${CUDA_FP}"
echo "jobs:  ${BUILD_JOBS}"
echo "nvcc:  $(command -v nvcc || echo '<not found>')"

./configure "${CONFIG_FLAGS[@]}"
make -j"${BUILD_JOBS}"

if [ -z "${SOURCE_BINARY}" ]; then
  echo "Set SOURCE_BINARY to the exact executable emitted by this build; refusing wildcard/latest selection." >&2
  exit 1
fi
GPU_BIN="${SOURCE_BINARY}"
case "${GPU_BIN}" in /*) ;; *) echo "SOURCE_BINARY must be absolute" >&2; exit 1;; esac
if [ ! -x "${GPU_BIN}" ]; then
  echo "SOURCE_BINARY is not an executable: ${GPU_BIN}" >&2
  exit 1
fi
POST_BUILD_STAT="$(stat -c '%s:%Y:%i' "${GPU_BIN}" 2>/dev/null || true)"
if [ -n "${BUILD_MARKER}" ] && { [ ! "${GPU_BIN}" -nt "${BUILD_MARKER}" ] || [ "${PRE_BUILD_STAT}" = "${POST_BUILD_STAT}" ]; }; then
  echo "SOURCE_BINARY was not emitted/relinked by this invocation" >&2; exit 1
fi

if [ -n "${OUTPUT_BINARY}" ]; then
  case "${OUTPUT_BINARY}" in /*) ;; *) echo "OUTPUT_BINARY must be absolute" >&2; exit 1;; esac
  if [ -e "${OUTPUT_BINARY}" ]; then
    echo "OUTPUT_BINARY already exists; refusing overwrite" >&2
    exit 1
  fi
  mkdir -p "$(dirname "${OUTPUT_BINARY}")"
  cp --reflink=auto "${GPU_BIN}" "${OUTPUT_BINARY}"
  chmod +x "${OUTPUT_BINARY}"
  GPU_BIN="${OUTPUT_BINARY}"
fi

echo "binary: ${GPU_BIN}"
echo "binary_sha256: $(sha256sum "${GPU_BIN}" | awk '{print $1}')"

if [ "${SMOKE}" = "1" ]; then
  OUT="output_local_cuda_smoke_sm${ARCH}"
  echo "=== one-step smoke test ==="
  "${GPU_BIN}" "${INPUT}" max_step="${MAX_STEP}" stop_time=1e-12_s \
    plot_file="${OUT}" amr.plot_int=-1 amr.thermo.plot_int=-1
fi

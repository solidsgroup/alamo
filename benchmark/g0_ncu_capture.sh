#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

# On NOVA: ncu is in PATH after loading the cuda module.
# Locally: falls back to the nsight-compute install.
NCU="${NCU:-$(command -v ncu 2>/dev/null || echo "${ROOT}/.local/nsight-compute/usr/bin/ncu")}"
# Default to 3D A100 binary (sm_80) for NOVA; override for local 2D runs.
GPU_BIN="${GPU_BIN:-$(ls -t "${ROOT}"/bin/alamo_gpu-3d*cuda80*-g++ 2>/dev/null | head -1 || echo "${ROOT}/bin/alamo_gpu-2d-cuda86-g++")}"
INPUT="${INPUT:-input_3d_flame_128}"
OUT="${OUT:-benchmark/profiles_alpha1/g0_ncu}"

if [ ! -x "${NCU}" ]; then
    echo "ncu not found at ${NCU}; install Nsight Compute or set NCU=/path/to/ncu" >&2
    exit 1
fi

COMMON_ARGS=(
    "${INPUT}"
    "max_step=2"
    "stop_time=1e99"
    "amr.plot_int=-1"
    "amr.thermo.plot_int=-1"
    "elastic.type=static"
    "elastic.on=1"
    "elastic.interval=1"
    "elastic.tstart=0.0"
    "elastic.solver.verbose=0"
    "elastic.print_model=0"
    "elastic.solver.max_iter=100"
    "elastic.solver.nriters=20"
)

run_ncu() {
    local label="$1"
    local kernel_regex="$2"
    local plot_dir="${OUT}_${label}_plot"
    "${NCU}" \
        --target-processes all \
        --set default \
        --kernel-name-base function \
        --kernel-name "regex:${kernel_regex}" \
        --launch-count 4 \
        --csv \
        --page raw \
        --print-kernel-base demangled \
        --export "${OUT}_${label}" \
        --force-overwrite \
        mpiexec -np 1 "${GPU_BIN}" "${COMMON_ARGS[@]}" "plot_file=${plot_dir}"
}

mkdir -p "$(dirname "${OUT}")"
run_ncu "flame_advance"    ".*Flame.*Advance.*|.*flame.*advance.*"
run_ncu "elastic_fapply"   ".*Elastic.*Fapply.*|.*elastic.*fapply.*"
run_ncu "elastic_diagonal" ".*Elastic.*Diagonal.*|.*elastic.*diagonal.*"
run_ncu "newton_prepare"   ".*Newton.*prepareForSolve.*|.*newton.*prepare.*"
run_ncu "operator_interp"  ".*interpolat.*"

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

# AMReX launches all GPU work via amrex::launch_global, so ncu only sees the
# function name "launch_global" -- name-regex matching never hits the physics
# kernels (leaves the export empty while ncu still exits 0). Target the
# TinyProfiler NVTX region instead. Note: Nsight Compute 2025.x has no "default"
# set (it's basic/detailed/full/...), so "--set default" collects nothing; use
# basic (override NCU_SET=full|detailed). --section-folder guards section lookup.
NCU_DIR="$(dirname "$(readlink -f "${NCU}")")"
SECTION_DIR="$(find "${NCU_DIR}/.." -maxdepth 4 -name 'SpeedOfLight.section' -printf '%h\n' 2>/dev/null | head -1 || true)"
if [ -n "${SECTION_DIR}" ]; then
    METRIC_ARGS=(--set "${NCU_SET:-basic}" --section-folder "${SECTION_DIR}")
    echo "ncu sections = ${SECTION_DIR}"
else
    METRIC_ARGS=(--metrics "${NCU_METRICS:-gpu__time_duration.sum,sm__throughput.avg.pct_of_peak_sustained_elapsed,gpu__compute_memory_throughput.avg.pct_of_peak_sustained_elapsed,sm__warps_active.avg.pct_of_peak_sustained_active,launch__registers_per_thread,launch__waves_per_multiprocessor}")
    echo "ncu sections = NOT FOUND near ${NCU_DIR} -- using explicit --metrics"
fi

NCU_FAILS=0
run_ncu() {
    local label="$1"
    local nvtx="$2"
    local plot_dir="${OUT}_${label}_plot"
    "${NCU}" \
        --target-processes all \
        "${METRIC_ARGS[@]}" \
        --nvtx --nvtx-include "${nvtx}" \
        --launch-count 4 \
        --csv \
        --page raw \
        --print-kernel-base demangled \
        --export "${OUT}_${label}" \
        --force-overwrite \
        mpiexec -np 1 "${GPU_BIN}" "${COMMON_ARGS[@]}" "plot_file=${plot_dir}" || true
    if [ -s "${OUT}_${label}.ncu-rep" ]; then
        echo "    ncu: ${label} -> $(du -h "${OUT}_${label}.ncu-rep" | cut -f1)"
    else
        echo "!!! ncu: ${label} produced NO report -- check 'No kernels were profiled' / 'No metrics' warnings above" >&2
        NCU_FAILS=$((NCU_FAILS + 1))
    fi
}

mkdir -p "$(dirname "${OUT}")"
run_ncu "elastic_fapply"   "*Fapply*/"
run_ncu "operator_interp"  "*interpolation*/"
run_ncu "elastic_diagonal" "*Diagonal*/"
run_ncu "bottom_bicgstab"  "*bicgstab*/"
[ "${NCU_FAILS}" -eq 0 ] || { echo "ncu: ${NCU_FAILS} pass(es) produced no data" >&2; exit 1; }

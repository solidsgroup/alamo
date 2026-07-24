#!/usr/bin/env bash
# Run the fullfeedback AP monopropellant case at several pressures for one
# candidate (pre_exponential, activation_temperature) pair, and measure each
# run's steady-state regression rate. Used both standalone and as the
# per-evaluation step inside scripts/optimize_fullfeedback.py.
#
# Usage:
#   scripts/run_pressure_sweep.sh <pre_exponential> <activation_temperature_K> \
#       <workdir> <results_csv> <pressure_MPa> [<pressure_MPa> ...]
#
# Environment:
#   LOWMACH_BIN   path to the lowmach binary
#                 (default: /home/mungerct/research/alamo/bin/lowmach-2d-hdf5-clang++)
#   TEMPLATE      path to the templated input
#                 (default: input.lm.ap_monopropellant_fullfeedback.template,
#                 resolved relative to this script's repo root)
#   KEEP_OUTPUT   if set to 1, do not delete each run's plotfile directory
#                 after measuring (useful for debugging)
#
# Writes <results_csv> with header "pressure_mpa,reg_rate_cm_s" and one row
# per requested pressure (reg_rate_cm_s is empty if that run never reached a
# steady-burning window). All pressures run concurrently.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

LOWMACH_BIN="${LOWMACH_BIN:-/home/mungerct/research/alamo/bin/lowmach-2d-hdf5-clang++}"
TEMPLATE="${TEMPLATE:-${REPO_ROOT}/input.lm.ap_monopropellant_fullfeedback.template}"
KEEP_OUTPUT="${KEEP_OUTPUT:-0}"

if [[ $# -lt 5 ]]; then
    echo "usage: $0 <pre_exponential> <activation_temperature_K> <workdir> <results_csv> <pressure_MPa> [<pressure_MPa> ...]" >&2
    exit 2
fi

PRE_EXPONENTIAL="$1"; shift
ACTIVATION_TEMPERATURE="$1"; shift
WORKDIR="$1"; shift
RESULTS_CSV="$1"; shift
PRESSURES=("$@")

if [[ ! -x "${LOWMACH_BIN}" ]]; then
    echo "error: lowmach binary not found or not executable: ${LOWMACH_BIN}" >&2
    echo "       set LOWMACH_BIN to override" >&2
    exit 1
fi
if [[ ! -f "${TEMPLATE}" ]]; then
    echo "error: template not found: ${TEMPLATE}" >&2
    exit 1
fi

mkdir -p "${WORKDIR}"
echo "pre_exponential=${PRE_EXPONENTIAL} activation_temperature=${ACTIVATION_TEMPERATURE}" \
    > "${WORKDIR}/params.txt"

pids=()
for p in "${PRESSURES[@]}"; do
    run_dir="${WORKDIR}/P${p}"
    mkdir -p "${run_dir}"
    input_file="${run_dir}/input"
    python3 "${SCRIPT_DIR}/render_input.py" \
        --template "${TEMPLATE}" \
        --pressure-mpa "${p}" \
        --pre-exponential "${PRE_EXPONENTIAL}" \
        --activation-temperature "${ACTIVATION_TEMPERATURE}" \
        --out "${input_file}" > "${run_dir}/render.log"

    (
        cd "${run_dir}"
        "${LOWMACH_BIN}" input "plot_file=${run_dir}/output" \
            > "${run_dir}/run.log" 2>&1
    ) &
    pids+=($!)
    echo "started P=${p} MPa (pid $!) -> ${run_dir}"
done

failures=0
for pid in "${pids[@]}"; do
    wait "${pid}" || failures=$((failures + 1))
done
if [[ "${failures}" -gt 0 ]]; then
    echo "warning: ${failures} sim(s) exited nonzero; check run.log in each P<pressure> dir" >&2
fi

{
    echo "pressure_mpa,reg_rate_cm_s"
    for p in "${PRESSURES[@]}"; do
        run_dir="${WORKDIR}/P${p}"
        rate="$(python3 "${SCRIPT_DIR}/regression_rate.py" --rate-only --unit cm/s \
            "${run_dir}/output" 2>>"${run_dir}/rate.log" || true)"
        if [[ "${rate}" == "nan" || -z "${rate}" ]]; then
            echo "${p},"
        else
            echo "${p},${rate}"
        fi
        if [[ "${KEEP_OUTPUT}" != "1" ]]; then
            rm -rf "${run_dir}/output"
        fi
    done
} > "${RESULTS_CSV}"

echo "Wrote ${RESULTS_CSV}"

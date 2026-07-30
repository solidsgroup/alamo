#!/usr/bin/env bash
# Run the AP monopropellant case at several pressures for one candidate
# (rate_multiplier, activation_temperature) pair, and measure each run's
# steady-state regression rate. Used both standalone and as the
# per-evaluation step inside scripts/optimize_ap_regression.py.
#
# Usage:
#   scripts/run_pressure_sweep.sh <rate_multiplier> <activation_temperature_K> \
#       <workdir> <results_csv> <pressure_MPa> [<pressure_MPa> ...]
#
# Environment:
#   LOWMACH_BIN   path to the lowmach binary
#                 (default: /home/mungerct/research/alamo/bin/lowmach-2d-clang++)
#   TEMPLATE      path to the templated input
#                 (default: input.lm.ap_monopropellant.template,
#                 resolved relative to this script's repo root)
#   RATE_UNIT     unit forwarded to regression_rate.py --unit (default: mm/s,
#                 matching AP_reg_rate.csv)
#   KEEP_OUTPUT   if set to 1, do not delete each run's plotfile directory
#                 after measuring (useful for debugging)
#   MIN_TIME_LOW  seconds of simulated time to exclude from the start of the
#                 run at the lowest requested pressure, to skip its startup
#                 transient before measuring the regression rate (default:
#                 0.0, no exclusion)
#   MIN_TIME_HIGH same as MIN_TIME_LOW but for the highest requested
#                 pressure (default: 0.0). Pressures between the lowest and
#                 highest requested pressure get a --min-time linearly
#                 interpolated between MIN_TIME_LOW and MIN_TIME_HIGH.
#
# Writes <results_csv> with header "pressure_mpa,reg_rate_<unit>" and one row
# per requested pressure (empty if that run never reached a steady-burning
# window). All pressures run concurrently.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

LOWMACH_BIN="${LOWMACH_BIN:-/home/mungerct/research/alamo/bin/lowmach-2d-clang++}"
TEMPLATE="${TEMPLATE:-${REPO_ROOT}/input.lm.ap_monopropellant.template}"
RATE_UNIT="${RATE_UNIT:-mm/s}"
RATE_UNIT_SUFFIX="${RATE_UNIT//\//_}"
KEEP_OUTPUT="${KEEP_OUTPUT:-0}"
MIN_TIME_LOW="${MIN_TIME_LOW:-0.0}"
MIN_TIME_HIGH="${MIN_TIME_HIGH:-0.0}"

if [[ $# -lt 5 ]]; then
    echo "usage: $0 <rate_multiplier> <activation_temperature_K> <workdir> <results_csv> <pressure_MPa> [<pressure_MPa> ...]" >&2
    exit 2
fi

RATE_MULTIPLIER="$1"; shift
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
echo "rate_multiplier=${RATE_MULTIPLIER} activation_temperature=${ACTIVATION_TEMPERATURE}" \
    > "${WORKDIR}/params.txt"

P_MIN="$(printf '%s\n' "${PRESSURES[@]}" | sort -g | head -1)"
P_MAX="$(printf '%s\n' "${PRESSURES[@]}" | sort -g | tail -1)"

# Linearly interpolate --min-time between MIN_TIME_LOW (at P_MIN) and
# MIN_TIME_HIGH (at P_MAX) for this pressure.
min_time_for() {
    awk -v p="$1" -v pmin="${P_MIN}" -v pmax="${P_MAX}" \
        -v tlo="${MIN_TIME_LOW}" -v thi="${MIN_TIME_HIGH}" \
        'BEGIN {
            if (pmax == pmin) { t = tlo }
            else { t = tlo + (thi - tlo) * (p - pmin) / (pmax - pmin) }
            printf "%.10g", t
        }'
}

pids=()
for p in "${PRESSURES[@]}"; do
    run_dir="${WORKDIR}/P${p}"
    mkdir -p "${run_dir}"
    input_file="${run_dir}/input"
    python3 "${SCRIPT_DIR}/render_input.py" \
        --template "${TEMPLATE}" \
        --pressure-mpa "${p}" \
        --rate-multiplier "${RATE_MULTIPLIER}" \
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
    echo "pressure_mpa,reg_rate_${RATE_UNIT_SUFFIX}"
    for p in "${PRESSURES[@]}"; do
        run_dir="${WORKDIR}/P${p}"
        min_time="$(min_time_for "${p}")"
        rate="$(python3 "${SCRIPT_DIR}/regression_rate.py" --rate-only --unit "${RATE_UNIT}" \
            --min-time "${min_time}" \
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

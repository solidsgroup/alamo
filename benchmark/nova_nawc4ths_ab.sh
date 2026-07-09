#!/usr/bin/env bash
# Submit NAWC Motor No. 6 quarter-domain GPU comparison runs on NOVA.
#
# References:
#   - benchmark/NOVA_SLURM_RUNBOOK.md for NOVA account/module/GRES assumptions
#   - benchmark/build_alamo_nova_3d.sh for the 3D CUDA build job
#   - benchmark/nova_flame_gpu_3d_a2.slurm for the 3D GPU launch path
#
# Dry run, default:
#   bash benchmark/nova_nawc4ths_ab.sh
#
# Submit both cases using an existing 3D GPU build:
#   bash benchmark/nova_nawc4ths_ab.sh --submit
#
# Submit a build first, then run both cases after the build succeeds:
#   bash benchmark/nova_nawc4ths_ab.sh --build --submit
#
# Useful knobs:
#   GPU_TYPE=a100|v100|h200      default: a100
#   PARTITION=nova|scavenger     default: nova
#   ARCHES="80"                  default: arch matching GPU_TYPE
#   MODE=bench|fast              default: bench
#   ELASTIC_EXTRA_ARGS="..."     appended only to the elastic case
#   NOELASTIC_EXTRA_ARGS="..."   appended only to the no-elastic case
#   COMMON_EXTRA_ARGS="..."      appended to both cases
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

GPU_TYPE="${GPU_TYPE:-a100}"
PARTITION="${PARTITION:-nova}"
ACCOUNT="${ACCOUNT:-brunnels}"
CPUS_PER_TASK="${CPUS_PER_TASK:-8}"
MEM="${MEM:-64G}"
TIME_LIMIT="${TIME_LIMIT:-04:00:00}"
MODE="${MODE:-bench}"
ASYNC_OUT="${ASYNC_OUT:-0}"
PROFILE_KIND="${PROFILE_KIND:-summary}"
COMMON_EXTRA_ARGS="${COMMON_EXTRA_ARGS:-}"
ELASTIC_EXTRA_ARGS="${ELASTIC_EXTRA_ARGS:-elastic.solver.verbose=2}"
NOELASTIC_EXTRA_ARGS="${NOELASTIC_EXTRA_ARGS:-}"
RUN_SCRIPT="${RUN_SCRIPT:-benchmark/nova_flame_gpu_3d_a2.slurm}"
BUILD_SCRIPT="${BUILD_SCRIPT:-benchmark/build_alamo_nova_3d.sh}"

SUBMIT=0
DO_BUILD=0
CASE_FILTER="both"

usage() {
  sed -n '1,32p' "$0"
  cat <<USAGE

Options:
  --submit          Submit jobs with sbatch. Without this, print commands only.
  --build           Submit the 3D NOVA build first and depend run jobs on it.
  --elastic-only    Submit/print only the elastic-enabled case.
  --noelastic-only  Submit/print only the disabled-elastic case.
  --help            Show this help.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --submit) SUBMIT=1 ;;
    --build) DO_BUILD=1 ;;
    --elastic-only) CASE_FILTER="elastic" ;;
    --noelastic-only) CASE_FILTER="noelastic" ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
  shift
done

case "${GPU_TYPE}" in
  v100) DEFAULT_ARCH=70 ;;
  a100) DEFAULT_ARCH=80 ;;
  h200) DEFAULT_ARCH=90 ;;
  *) echo "Unknown GPU_TYPE=${GPU_TYPE} (use v100, a100, or h200)" >&2; exit 2 ;;
esac
ARCHES="${ARCHES:-${DEFAULT_ARCH}}"

ELASTIC_INPUT="nawc6_input_4ths.in"
NOELASTIC_INPUT="nawc6_input_4ths_noelastic.in"
for f in "${ELASTIC_INPUT}" "${NOELASTIC_INPUT}" "${RUN_SCRIPT}" "${BUILD_SCRIPT}"; do
  [[ -f "${f}" ]] || { echo "Missing required file: ${f}" >&2; exit 2; }
done

print_build_command() {
  printf 'ACCOUNT=%q BUILD_PARTITION=%q ARCHES=%q sh %q\n' \
    "${ACCOUNT}" "${PARTITION}" "${ARCHES}" "${BUILD_SCRIPT}"
}

submit_build() {
  echo "=== submitting 3D NOVA build (${BUILD_SCRIPT}, ARCHES=${ARCHES}) ===" >&2
  local output jobid
  output="$(ACCOUNT="${ACCOUNT}" BUILD_PARTITION="${PARTITION}" ARCHES="${ARCHES}" sh "${BUILD_SCRIPT}")"
  printf '%s\n' "${output}" >&2
  jobid="$(printf '%s\n' "${output}" | sed -n 's/.*submitted job \([0-9][0-9]*\).*/\1/p' | tail -1)"
  if [[ -z "${jobid}" ]]; then
    echo "ERROR: could not parse build job id from ${BUILD_SCRIPT} output" >&2
    exit 1
  fi
  printf '%s\n' "${jobid}"
}

append_common_extra() {
  local case_extra="$1"
  if [[ -n "${COMMON_EXTRA_ARGS}" && -n "${case_extra}" ]]; then
    printf '%s %s' "${COMMON_EXTRA_ARGS}" "${case_extra}"
  elif [[ -n "${COMMON_EXTRA_ARGS}" ]]; then
    printf '%s' "${COMMON_EXTRA_ARGS}"
  else
    printf '%s' "${case_extra}"
  fi
}

emit_case() {
  local case_id="$1"
  local input="$2"
  local extra_args="$3"
  local dependency="${4:-}"
  local job="nawc4ths_${case_id}"
  local sbatch_opts=(
    --partition="${PARTITION}"
    --nodes=1
    --gres="gpu:${GPU_TYPE}:1"
    --ntasks=1
    --ntasks-per-node=1
    --cpus-per-task="${CPUS_PER_TASK}"
    --mem="${MEM}"
    --time="${TIME_LIMIT}"
    --job-name="${job}"
    --output="${job}.%j.out"
    --error="${job}.%j.err"
  )
  [[ -n "${dependency}" ]] && sbatch_opts+=(--dependency="afterok:${dependency}")

  if [[ "${SUBMIT}" -eq 1 ]]; then
    command -v sbatch >/dev/null 2>&1 || { echo "ERROR: sbatch not found" >&2; exit 1; }
    echo "=== submitting ${case_id}: input=${input} gpu=${GPU_TYPE} mode=${MODE} ==="
    INPUT="${input}" GPU_TYPE="${GPU_TYPE}" MODE="${MODE}" ASYNC_OUT="${ASYNC_OUT}" \
      PROFILE_KIND="${PROFILE_KIND}" PARTITION="${PARTITION}" EXTRA_ARGS="${extra_args}" \
      sbatch "${sbatch_opts[@]}" "${RUN_SCRIPT}"
  else
    printf 'INPUT=%q GPU_TYPE=%q MODE=%q ASYNC_OUT=%q PROFILE_KIND=%q PARTITION=%q EXTRA_ARGS=%q sbatch' \
      "${input}" "${GPU_TYPE}" "${MODE}" "${ASYNC_OUT}" "${PROFILE_KIND}" "${PARTITION}" "${extra_args}"
    printf ' %q' "${sbatch_opts[@]}" "${RUN_SCRIPT}"
    printf '\n'
  fi
}

BUILD_JOBID=""
if [[ "${DO_BUILD}" -eq 1 ]]; then
  if [[ "${SUBMIT}" -eq 1 ]]; then
    BUILD_JOBID="$(submit_build)"
    echo "Build job id: ${BUILD_JOBID}"
  else
    echo "# Build command:"
    print_build_command
    echo "# Run commands below should use --dependency=afterok:<build_jobid> after submission."
  fi
fi

ELASTIC_ARGS="$(append_common_extra "${ELASTIC_EXTRA_ARGS}")"
NOELASTIC_ARGS="$(append_common_extra "${NOELASTIC_EXTRA_ARGS}")"

case "${CASE_FILTER}" in
  both)
    emit_case "elastic" "${ELASTIC_INPUT}" "${ELASTIC_ARGS}" "${BUILD_JOBID}"
    emit_case "noelastic" "${NOELASTIC_INPUT}" "${NOELASTIC_ARGS}" "${BUILD_JOBID}"
    ;;
  elastic)
    emit_case "elastic" "${ELASTIC_INPUT}" "${ELASTIC_ARGS}" "${BUILD_JOBID}"
    ;;
  noelastic)
    emit_case "noelastic" "${NOELASTIC_INPUT}" "${NOELASTIC_ARGS}" "${BUILD_JOBID}"
    ;;
  *) echo "Internal error: bad CASE_FILTER=${CASE_FILTER}" >&2; exit 2 ;;
esac

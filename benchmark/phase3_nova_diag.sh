#!/usr/bin/env bash
#
# phase3_nova_diag.sh -- submit focused NOVA diagnostics for Phase-3 CUDA 700.
#
# This does not run the full scaling sweep. It submits one tiny single-GPU V100
# diagnostic job per requested tool so the first illegal access can be localized.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ALAMO_DIR="${ALAMO_DIR:-${ROOT_DIR}}"

GPU_TYPE="${GPU_TYPE:-v100}"
INPUT="${INPUT:-input_3d_centre_bore_128}"
TOOLS="${TOOLS:-blocking sanitizer}"  # also: ncu nsys
MANAGED_ARENA="${MANAGED_ARENA:-1}"
DEPENDENCY="${DEPENDENCY:-}"
GPU_NODE_CPUS="${GPU_NODE_CPUS:-32}"
DRY_RUN=0
SCAVENGER=0

usage() {
    cat <<'EOF'
Usage: bash benchmark/phase3_nova_diag.sh [--dry-run] [--scavenger]

Submits focused single-GPU diagnostics for the Phase-3 CUDA error 700 failure.
Run from the NOVA checkout after the 3D GPU binary has been built.

Env knobs:
  ALAMO_DIR      repo root / working tree (default: repo root inferred)
  GPU_TYPE       v100|a100|h200 (default: v100)
  INPUT          input file (default: input_3d_centre_bore_128)
  TOOLS          "blocking", "sanitizer", "ncu", "nsys", or any combo (default: "blocking sanitizer")
  --scavenger    submit to the scavenger partition instead of nova
  MANAGED_ARENA  1|0 passed to amrex.the_arena_is_managed (default: 1)
  DEPENDENCY     optional Slurm dependency, bare job id or full afterok:...
  GPU_NODE_CPUS  CPUs on one NOVA GPU node (default: 72)
EOF
}

for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
        --scavenger) SCAVENGER=1 ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "error: unknown argument '$arg' (try --help)" >&2
            exit 2
            ;;
    esac
done

case "${GPU_TYPE}" in
    v100|a100|h200) ;;
    *) echo "error: GPU_TYPE must be v100, a100 or h200 (got ${GPU_TYPE})" >&2; exit 2 ;;
esac
if ! [[ "${GPU_NODE_CPUS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "error: GPU_NODE_CPUS must be a positive integer (got ${GPU_NODE_CPUS})" >&2
    exit 2
fi

dependency_opt() {
    local dependency="${1:-}"
    local deptext=""
    if [[ -z "$dependency" ]]; then
        return 0
    fi
    case "$dependency" in
        afterok:*|afterany:*|afternotok:*|singleton|expand:*)
            deptext="$dependency"
            ;;
        *)
            deptext="afterok:$dependency"
            ;;
    esac
    printf '%s' "--dependency=${deptext}"
}

cd "${ALAMO_DIR}"

echo "============================================================"
echo " Phase-3 NOVA CUDA-700 diagnostics"
echo " repo          = ${ALAMO_DIR}"
echo " GPU_TYPE      = ${GPU_TYPE}"
echo " INPUT         = ${INPUT}"
echo " TOOLS         = ${TOOLS}"
echo " MANAGED_ARENA = ${MANAGED_ARENA}"
echo " GPU_NODE_CPUS = ${GPU_NODE_CPUS}"
echo " SCAVENGER     = ${SCAVENGER}"
echo " DRY_RUN       = ${DRY_RUN}"
if [[ -n "${DEPENDENCY}" ]]; then
    echo " DEPENDENCY    = ${DEPENDENCY}"
fi
echo "============================================================"

for tool in ${TOOLS}; do
    case "${tool}" in
        blocking|sanitizer|ncu|nsys) ;;
        *) echo "error: unknown diagnostic tool '${tool}'" >&2; exit 2 ;;
    esac

    sbatch_opts=(--nodes=1 --gres="gpu:${GPU_TYPE}:1" --ntasks=1 --ntasks-per-node=1 --cpus-per-task="${GPU_NODE_CPUS}")
    if [[ "${SCAVENGER}" -eq 1 ]]; then
        sbatch_opts+=(--partition=scavenger)
    fi
    depopt="$(dependency_opt "${DEPENDENCY}")"
    if [[ -n "${depopt}" ]]; then
        sbatch_opts+=("${depopt}")
    fi

    echo
    echo "Diagnostic: ${tool}"
    printf 'Command: INPUT=%q GPU_TYPE=%q DIAG_TOOL=%q MANAGED_ARENA=%q sbatch' \
      "${INPUT}" "${GPU_TYPE}" "${tool}" "${MANAGED_ARENA}"
    printf ' %q' "${sbatch_opts[@]}"
    printf ' benchmark/nova_flame_gpu_3d_diag.slurm\n'

    if [[ "${DRY_RUN}" -eq 0 ]]; then
        INPUT="${INPUT}" GPU_TYPE="${GPU_TYPE}" DIAG_TOOL="${tool}" MANAGED_ARENA="${MANAGED_ARENA}" \
          sbatch "${sbatch_opts[@]}" benchmark/nova_flame_gpu_3d_diag.slurm
    fi
done

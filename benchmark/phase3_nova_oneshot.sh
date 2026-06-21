#!/usr/bin/env bash
#
# phase3_nova_oneshot.sh -- one-shot NOVA driver for Phase 3.
#
# Run this on the NOVA login node from the repo root (or anywhere with the repo
# checked out). It will:
#   1) print the local Phase-3 memory budget table,
#   2) submit the 3D NOVA build job,
#   3) submit the Phase-3 sweep matrix with an afterok dependency on the build.
#
# The actual work is done by:
#   - benchmark/build_alamo_nova_3d.sh
#   - benchmark/phase3_scaling_sweep.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ALAMO_DIR="${ALAMO_DIR:-${ROOT_DIR}}"
SWEEP_MODE="${SWEEP_MODE:-strong}"
GPU_TYPE="${GPU_TYPE:-a100}"
SIZES="${SIZES:-128 256 512}"
GPUS="${GPUS:-1 2 4 8}"
OUT="${OUT:-commands.txt}"
BUILD_ONLY="${BUILD_ONLY:-0}"
DRY_RUN=0

usage() {
    cat <<'EOF'
Usage: bash benchmark/phase3_nova_oneshot.sh [--dry-run] [--build-only]

This is the Phase-3 NOVA launcher. It prints the memory budget table, submits
the 3D build job, then submits the strong/weak scaling sweep after the build
completes successfully.

Env knobs:
  ALAMO_DIR   repo root / working tree (default: repo root inferred from script)
  SWEEP_MODE  strong|weak (default: strong)
  GPU_TYPE    a100|h200 (default: a100)
  SIZES       size list for strong scaling (default: "128 256 512")
  GPUS        GPU list (default: "1 2 4 8")
  OUT         sweep command output file (default: commands.txt)
  BUILD_ONLY  set to 1 to submit the build and stop
EOF
}

for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
        --build-only) BUILD_ONLY=1 ;;
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

if [[ "$SWEEP_MODE" != "strong" && "$SWEEP_MODE" != "weak" ]]; then
    echo "error: SWEEP_MODE must be 'strong' or 'weak' (got '$SWEEP_MODE')" >&2
    exit 2
fi

cd "$ALAMO_DIR"

echo "============================================================"
echo " Phase-3 NOVA one-shot launcher"
echo " repo        = $ALAMO_DIR"
echo " sweep mode  = $SWEEP_MODE"
echo " GPU_TYPE    = $GPU_TYPE"
echo " SIZES       = $SIZES"
echo " GPUS        = $GPUS"
echo " OUT         = $OUT"
echo " BUILD_ONLY  = $BUILD_ONLY"
echo " DRY_RUN     = $DRY_RUN"
echo "============================================================"

echo
echo "=== Phase-3 memory budget (local preflight) ==="
if command -v python3 >/dev/null 2>&1; then
    python3 "${ROOT_DIR}/benchmark/phase3_memory_budget.py"
else
    echo "WARN: python3 not found; skipping the local memory budget preflight."
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo
    echo "=== Dry run only ==="
    echo "Would submit:"
    echo "  sh ${ROOT_DIR}/benchmark/build_alamo_nova_3d.sh"
    echo "  DEPENDENCY=<build-jobid> MODE=${SWEEP_MODE} GPU_TYPE=${GPU_TYPE} \\"
    echo "    SIZES='${SIZES}' GPUS='${GPUS}' OUT='${OUT}' \\"
    echo "    bash ${ROOT_DIR}/benchmark/phase3_scaling_sweep.sh --submit"
    exit 0
fi

BUILD_LOG="$(mktemp -t phase3_nova_build.XXXXXX.log)"
trap 'rm -f "$BUILD_LOG"' EXIT

echo
echo "=== Submitting 3D build job ==="
sh "${ROOT_DIR}/benchmark/build_alamo_nova_3d.sh" 2>&1 | tee "$BUILD_LOG"

BUILD_JOBID="$(grep -oE 'submitted job [0-9]+' "$BUILD_LOG" | tail -1 | awk '{print $3}')"
if [[ -z "$BUILD_JOBID" ]]; then
    echo "error: could not parse build job id from ${BUILD_LOG}" >&2
    exit 1
fi

echo
echo "Build job id: ${BUILD_JOBID}"

if [[ "$BUILD_ONLY" -eq 1 ]]; then
    echo "BUILD_ONLY=1 set; stopping after build submission."
    exit 0
fi

echo
echo "=== Submitting Phase-3 sweep matrix (afterok build) ==="
DEPENDENCY="$BUILD_JOBID" \
MODE="$SWEEP_MODE" \
GPU_TYPE="$GPU_TYPE" \
SIZES="$SIZES" \
GPUS="$GPUS" \
OUT="$OUT" \
bash "${ROOT_DIR}/benchmark/phase3_scaling_sweep.sh" --submit

echo
echo "All jobs staged."
echo "Build job:  ${BUILD_JOBID}"
echo "Sweep file: ${OUT}"

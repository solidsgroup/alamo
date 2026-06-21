#!/usr/bin/env bash
#
# phase3_nova_oneshot.sh -- one-shot NOVA driver for Phase 3.
#
# Run this on the NOVA login node from the repo root (or anywhere with the repo
# checked out). It will:
#   1) print the local Phase-3 memory budget table,
#   2) submit the 3D NOVA build job,
#   3) submit one serialized CPU build job,
#   4) submit the Phase-3 sweep matrix with afterok dependencies on the builds.
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
    echo "  sbatch --dependency=afterok:<gpu-build-jobid> <cpu-build-job>"
    echo "  GPU_DEPENDENCY=<gpu-build-jobid> CPU_DEPENDENCY=<cpu-build-jobid> \\"
    echo "    MODE=${SWEEP_MODE} GPU_TYPE=${GPU_TYPE} \\"
    echo "    SIZES='${SIZES}' GPUS='${GPUS}' OUT='${OUT}' \\"
    echo "    bash ${ROOT_DIR}/benchmark/phase3_scaling_sweep.sh --submit"
    exit 0
fi

BUILD_LOG="$(mktemp -t phase3_nova_build.XXXXXX.log)"
CPU_BUILD_SCRIPT="$(mktemp -t phase3_nova_cpu_build.XXXXXX.sbatch)"
trap 'rm -f "$BUILD_LOG" "$CPU_BUILD_SCRIPT"' EXIT

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

cat > "$CPU_BUILD_SCRIPT" <<'EOF'
#!/bin/bash
#SBATCH -A brunnels
#SBATCH -J alamo_cpu_build_3d
#SBATCH --partition=nova
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=alamo_cpu_build_3d.%j.out
#SBATCH --error=alamo_cpu_build_3d.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jackplum@iastate.edu
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

module purge 2>/dev/null || true
module load gcc 2>/dev/null || module load gcc/12 2>/dev/null || true
module load openmpi 2>/dev/null || module load openmpi4 2>/dev/null || true

echo "=== configuring 3D CPU build ==="
./configure --comp=g++ --dim 3 --profile --get-eigen
echo "=== building 3D CPU binary ==="
make -j"${SLURM_CPUS_PER_TASK:-16}"
echo "=== CPU binaries ==="
ls -lh bin/alamo-3d*-g++ 2>/dev/null | grep -v cuda
EOF

echo
echo "=== Submitting 3D CPU build job (afterok GPU build) ==="
CPU_BUILD_JOBID="$(sbatch --parsable --dependency="afterok:${BUILD_JOBID}" "$CPU_BUILD_SCRIPT")"
CPU_BUILD_JOBID="${CPU_BUILD_JOBID%%;*}"
if [[ -z "$CPU_BUILD_JOBID" ]]; then
    echo "error: could not parse CPU build job id" >&2
    exit 1
fi
echo "CPU build job id: ${CPU_BUILD_JOBID}"

echo
echo "=== Submitting Phase-3 sweep matrix (afterok GPU/CPU builds) ==="
GPU_DEPENDENCY="$BUILD_JOBID" \
CPU_DEPENDENCY="$CPU_BUILD_JOBID" \
MODE="$SWEEP_MODE" \
GPU_TYPE="$GPU_TYPE" \
SIZES="$SIZES" \
GPUS="$GPUS" \
OUT="$OUT" \
bash "${ROOT_DIR}/benchmark/phase3_scaling_sweep.sh" --submit

echo
echo "All jobs staged."
echo "GPU build job: ${BUILD_JOBID}"
echo "CPU build job: ${CPU_BUILD_JOBID}"
echo "Sweep file:    ${OUT}"

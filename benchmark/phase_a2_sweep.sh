#!/usr/bin/env bash
#
# phase_a2_sweep.sh -- Phase-A2 combined flame+elastic crossover sweep.
#
# Roadmap: GPU_ROADMAP_V2.md task A2 -- "Combined flame+elastic A100
# baseline-of-record."  Re-runs the Phase-3 size matrix (128^3 + 256^3) with
# elastic.type=static enabled.  Single A100 only; no multi-GPU (multi-GPU is
# a confirmed Phase-3 negative, parked until D1).
#
# This script is staged for NOVA: by default it is a dry run that only emits
# the matrix of sbatch commands into $OUT.  Pass --submit on a SLURM host to
# actually submit.
#
# Inputs:  input_3d_centre_bore_128_a2, input_3d_centre_bore_256_a2
# Scripts: benchmark/nova_flame_gpu_3d_a2.slurm, benchmark/nova_flame_cpu_3d_a2.slurm
#
# Flags:
#   --submit           sbatch each command (requires SLURM; default: dry run)
#   --dry-run          write $OUT only, submit nothing (default)
#   --scavenger        use the scavenger partition instead of nova
#
# Env knobs:
#   GPU_TYPE=a100        GPU type for GRES (default a100)
#   GPU_NODE_CPUS=8      CPU cores on one NOVA GPU node
#   CPU_RANKS=64         MPI ranks for the CPU baseline (match Phase-3 R3)
#   SIZES="128 256"      Grid sizes to sweep (default both)
#   OUT=a2_commands.txt  Output file for the sbatch command matrix
#   DEPENDENCY=...       Optional SLURM dependency applied to every job
#   GPU_DEPENDENCY=...   Optional dependency for GPU jobs only
#   CPU_DEPENDENCY=...   Optional dependency for CPU jobs only
#
set -euo pipefail

GPU_TYPE="${GPU_TYPE:-a100}"
GPU_NODE_CPUS="${GPU_NODE_CPUS:-8}"
CPU_RANKS="${CPU_RANKS:-64}"
SIZES="${SIZES:-128 256}"
OUT="${OUT:-a2_commands.txt}"
DEPENDENCY="${DEPENDENCY:-}"
GPU_DEPENDENCY="${GPU_DEPENDENCY:-${DEPENDENCY}}"
CPU_DEPENDENCY="${CPU_DEPENDENCY:-${DEPENDENCY}}"
PARTITION="${PARTITION:-nova}"

SCRIPT_GPU="benchmark/nova_flame_gpu_3d_a2.slurm"
SCRIPT_CPU="benchmark/nova_flame_cpu_3d_a2.slurm"

SUBMIT=0
for arg in "$@"; do
    case "$arg" in
        --submit)     SUBMIT=1 ;;
        --dry-run)    SUBMIT=0 ;;
        --scavenger)  PARTITION=scavenger ;;
        -h|--help)
            sed -n '2,45p' "$0" | sed 's/^# \{0,1\}//'
            echo
            echo "Usage: bash benchmark/phase_a2_sweep.sh [--dry-run|--submit] [--scavenger]"
            echo "  --dry-run (default): write \$OUT, submit nothing."
            echo "  --submit           : sbatch each command (requires SLURM)."
            echo "  --scavenger        : submit to scavenger partition instead of nova."
            exit 0
            ;;
        *)
            echo "error: unknown argument '$arg' (try --help)" >&2
            exit 2
            ;;
    esac
done

case "$GPU_TYPE" in
    v100|a100|h200) ;;
    *) echo "error: GPU_TYPE must be v100, a100, or h200 (got '$GPU_TYPE')" >&2; exit 2 ;;
esac

echo "============================================================"
echo " Phase-A2 combined flame+elastic crossover sweep"
echo " STAGED FOR NOVA -- this is a no-op without SLURM."
echo "   GPU_TYPE      = $GPU_TYPE"
echo "   GPU_NODE_CPUS = $GPU_NODE_CPUS"
echo "   CPU_RANKS     = $CPU_RANKS"
echo "   SIZES         = $SIZES"
echo "   PARTITION     = $PARTITION"
echo "   OUT           = $OUT"
if [[ -n "$DEPENDENCY" ]]; then
    echo "   DEPENDENCY    = $DEPENDENCY"
fi
if [[ "$SUBMIT" -eq 1 ]]; then
    echo "   ACTION        = SUBMIT (sbatch each command)"
else
    echo "   ACTION        = DRY RUN (write commands only; submit nothing)"
fi
echo "============================================================"

dependency_opt() {
    local dependency="${1:-}"
    if [[ -z "$dependency" ]]; then
        return 0
    fi
    case "$dependency" in
        afterok:*|afterany:*|afternotok:*|singleton|expand:*)
            printf '%s' "--dependency=${dependency}"
            ;;
        *)
            printf '%s' "--dependency=afterok:${dependency}"
            ;;
    esac
}

input_for_size() {
    case "$1" in
        128) echo "input_3d_centre_bore_128_a2" ;;
        256) echo "input_3d_centre_bore_256_a2" ;;
        *) echo "error: unknown size $1 for A2 (use 128 or 256)" >&2; exit 2 ;;
    esac
}

emit_gpu() {
    local input="$1" size="$2"
    local cpus_per_task="${GPU_NODE_CPUS}"
    local mem_gb
    case "$size" in
        128) mem_gb=32G ;;
        256) mem_gb=64G ;;
        *) mem_gb=64G ;;
    esac
    local depopt
    depopt="$(dependency_opt "${GPU_DEPENDENCY}")"
    local sbatch_opts=(--partition="${PARTITION}" --nodes=1 --gres="gpu:${GPU_TYPE}:1" --ntasks=1 --ntasks-per-node=1 --cpus-per-task="${cpus_per_task}" --mem="${mem_gb}")
    [[ -n "$depopt" ]] && sbatch_opts+=("$depopt")

    local sbatch_opts_text=""
    printf -v sbatch_opts_text ' %q' "${sbatch_opts[@]}"

    local cmd="INPUT=${input} GPU_TYPE=${GPU_TYPE} NGPUS=1 PARTITION=${PARTITION} sbatch${sbatch_opts_text} ${SCRIPT_GPU}"
    printf '# A2 size=%s gpu x1 (%s) partition=%s\n%s\n\n' "$size" "$GPU_TYPE" "$PARTITION" "$cmd" >> "$OUT"

    if [[ "$SUBMIT" -eq 1 ]]; then
        if command -v sbatch >/dev/null 2>&1; then
            echo "submitting: $cmd"
            INPUT="${input}" GPU_TYPE="${GPU_TYPE}" NGPUS=1 PARTITION="${PARTITION}" \
                sbatch "${sbatch_opts[@]}" "${SCRIPT_GPU}"
        else
            echo "WARN: --submit requested but 'sbatch' not found; skipping: $cmd" >&2
        fi
    fi
}

emit_cpu() {
    local input="$1" size="$2"
    local mem_gb
    case "$size" in
        128) mem_gb=32G ;;
        256) mem_gb=64G ;;
        *) mem_gb=64G ;;
    esac
    local depopt
    depopt="$(dependency_opt "${CPU_DEPENDENCY}")"
    local sbatch_opts=(--partition="${PARTITION}" --nodes=1 --ntasks="${CPU_RANKS}" --cpus-per-task=1 --mem="${mem_gb}")
    [[ -n "$depopt" ]] && sbatch_opts+=("$depopt")

    local sbatch_opts_text=""
    printf -v sbatch_opts_text ' %q' "${sbatch_opts[@]}"

    local cmd="INPUT=${input} BUILD_CPU_IF_MISSING=0 PARTITION=${PARTITION} sbatch${sbatch_opts_text} ${SCRIPT_CPU}"
    printf '# A2 size=%s cpu-node baseline (%s ranks) partition=%s\n%s\n\n' "$size" "$CPU_RANKS" "$PARTITION" "$cmd" >> "$OUT"

    if [[ "$SUBMIT" -eq 1 ]]; then
        if command -v sbatch >/dev/null 2>&1; then
            echo "submitting: $cmd"
            INPUT="${input}" BUILD_CPU_IF_MISSING=0 PARTITION="${PARTITION}" \
                sbatch "${sbatch_opts[@]}" "${SCRIPT_CPU}"
        else
            echo "WARN: --submit requested but 'sbatch' not found; skipping: $cmd" >&2
        fi
    fi
}

: > "$OUT"
{
    echo "# Phase-A2 combined flame+elastic crossover command matrix"
    echo "# generated by benchmark/phase_a2_sweep.sh"
    echo "# GPU_TYPE=$GPU_TYPE  CPU_RANKS=$CPU_RANKS  SIZES='$SIZES'  PARTITION=$PARTITION"
    echo "# scripts: nova_flame_gpu_3d_a2.slurm, nova_flame_cpu_3d_a2.slurm"
    echo "#"
} >> "$OUT"

for size in $SIZES; do
    input="$(input_for_size "$size")"
    emit_gpu "$input" "$size"
    emit_cpu "$input" "$size"
done

echo
echo "Wrote command matrix to: $OUT"
if [[ "$SUBMIT" -eq 0 ]]; then
    echo "(dry run -- nothing submitted; re-run with --submit on a SLURM host)"
fi

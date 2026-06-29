#!/usr/bin/env bash
#
# phase3_scaling_sweep.sh -- Phase-3 crossover sweep driver (roadmap 3.5).
#
# Sweeps problem size x GPU count to find where the 3D GPU build beats the CPU
# node. This script is staged for NOVA: by default it is a dry run that only
# emits the matrix of sbatch commands into commands.txt -- it submits nothing
# and is a no-op without SLURM. Pass --submit (on a SLURM host) to actually sbatch.
#
# Inputs:  input_3d_centre_bore_128, input_3d_centre_bore, input_3d_centre_bore_512
# Scripts: nova_flame_gpu_3d.slurm, nova_flame_gpu_3d_multi.slurm, nova_flame_cpu_3d.slurm
#
# Current NOVA inventory notes (2026-06-21):
#   - wide CPU nodes: 96 logical CPUs, about 378 GB RAM, no GPU GRES
#   - GPU nodes:      72 logical CPUs, 2x V100 GRES, about 189 GB RAM
#
# Env knobs:
#   MODE=strong|weak   strong (default) = fixed size, vary GPU count;
#                      weak             = size grows with GPU count.
#   GPU_TYPE=...       GPU type label passed through to the SLURM scripts
#                      (default v100 on the currently observed NOVA GPU nodes).
#   GPUS_PER_NODE=2    GPU slots per NOVA GPU node.
#   GPU_NODE_CPUS=8    CPUs requested per GPU node, split across GPU ranks.
#   CPU_RANKS=64       CPU ranks for the full wide-node baseline.
#   SIZES="128 256 512"   grid sizes for the strong sweep.
#   SMALL_GPUS="1"     GPU counts for the smaller resolutions.
#   LARGE_GPUS="1 2"   GPU counts for the largest resolution.
#   GPUS="1 2"         GPU counts to sweep in weak mode.
#   OUT=commands.txt      output file for the emitted command matrix.
#   DEPENDENCY=...        optional SLURM dependency applied to every submitted job
#                         (bare jobid or full sbatch dependency syntax).
#   GPU_DEPENDENCY=...    optional dependency for GPU jobs only; defaults to
#                         DEPENDENCY.
#   CPU_DEPENDENCY=...    optional dependency for CPU jobs only; defaults to
#                         DEPENDENCY.
#
set -euo pipefail

MODE="${MODE:-strong}"
GPU_TYPE="${GPU_TYPE:-v100}"
GPUS_PER_NODE="${GPUS_PER_NODE:-2}"
GPU_NODE_CPUS="${GPU_NODE_CPUS:-8}"
CPU_RANKS="${CPU_RANKS:-64}"
SMALL_GPUS="${SMALL_GPUS:-1}"
LARGE_GPUS="${LARGE_GPUS:-1 2}"
SIZES="${SIZES:-128 256 512}"
GPUS="${GPUS:-1 2}"
OUT="${OUT:-commands.txt}"
DEPENDENCY="${DEPENDENCY:-}"
GPU_DEPENDENCY="${GPU_DEPENDENCY:-${DEPENDENCY}}"
CPU_DEPENDENCY="${CPU_DEPENDENCY:-${DEPENDENCY}}"

SCRIPT_GPU_SINGLE="benchmark/nova_flame_gpu_3d.slurm"
SCRIPT_GPU_MULTI="benchmark/nova_flame_gpu_3d_multi.slurm"
SCRIPT_CPU="benchmark/nova_flame_cpu_3d.slurm"

SUBMIT=0
SCAVENGER=0
for arg in "$@"; do
    case "$arg" in
        --submit) SUBMIT=1 ;;
        --dry-run) SUBMIT=0 ;;
        --scavenger) SCAVENGER=1 ;;
        -h|--help)
            sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'
            echo
            echo "Usage: bash benchmark/phase3_scaling_sweep.sh [--dry-run|--submit] [--scavenger]"
            echo "  --dry-run (default): write \$OUT, submit nothing."
            echo "  --submit           : sbatch each command (requires SLURM)."
            echo "  --scavenger        : use the scavenger partition instead of nova."
            exit 0
            ;;
        *)
            echo "error: unknown argument '$arg' (try --help)" >&2
            exit 2
            ;;
    esac
done

if [[ "$MODE" != "strong" && "$MODE" != "weak" ]]; then
    echo "error: MODE must be 'strong' or 'weak' (got '$MODE')" >&2
    exit 2
fi
if ! [[ "$GPUS_PER_NODE" =~ ^[1-9][0-9]*$ ]]; then
    echo "error: GPUS_PER_NODE must be a positive integer (got '$GPUS_PER_NODE')" >&2
    exit 2
fi
if ! [[ "$GPU_NODE_CPUS" =~ ^[1-9][0-9]*$ ]]; then
    echo "error: GPU_NODE_CPUS must be a positive integer (got '$GPU_NODE_CPUS')" >&2
    exit 2
fi
if ! [[ "$CPU_RANKS" =~ ^[1-9][0-9]*$ ]]; then
    echo "error: CPU_RANKS must be a positive integer (got '$CPU_RANKS')" >&2
    exit 2
fi
case "$GPU_TYPE" in
    v100|a100|h200) ;;
    *) echo "error: GPU_TYPE must be v100, a100, or h200 (got '$GPU_TYPE')" >&2; exit 2 ;;
esac

echo "============================================================"
echo " Phase-3 scaling/crossover sweep (roadmap 3.5)"
echo " STAGED FOR NOVA -- this is a no-op without SLURM."
echo "   MODE          = $MODE"
echo "   GPU_TYPE      = $GPU_TYPE"
echo "   GPUS_PER_NODE = $GPUS_PER_NODE"
echo "   GPU_NODE_CPUS = $GPU_NODE_CPUS"
echo "   CPU_RANKS     = $CPU_RANKS"
echo "   SMALL_GPUS    = $SMALL_GPUS"
echo "   LARGE_GPUS    = $LARGE_GPUS"
echo "   SIZES         = $SIZES"
echo "   GPUS          = $GPUS"
echo "   OUT           = $OUT"
if [[ -n "$DEPENDENCY" ]]; then
    echo "   DEPENDENCY    = $DEPENDENCY"
fi
if [[ -n "$GPU_DEPENDENCY" && "$GPU_DEPENDENCY" != "$DEPENDENCY" ]]; then
    echo "   GPU_DEPENDENCY = $GPU_DEPENDENCY"
fi
if [[ -n "$CPU_DEPENDENCY" && "$CPU_DEPENDENCY" != "$DEPENDENCY" ]]; then
    echo "   CPU_DEPENDENCY = $CPU_DEPENDENCY"
fi
echo "   SCAVENGER     = $SCAVENGER"
if [[ "$SUBMIT" -eq 1 ]]; then
    echo "   ACTION        = SUBMIT (sbatch each command)"
else
    echo "   ACTION        = DRY RUN (write commands only; submit nothing)"
fi
echo "============================================================"

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

emit() {
    local input="$1" ngpus="$2" script="$3" tag="$4" kind="$5"
    local dependency="" depopt=""
    local mem_gb
    case "$input" in
        *128*) mem_gb=32G ;;
        *256*|*512*) mem_gb=64G ;;
        *) mem_gb=64G ;;
    esac
    local sbatch_opts=()
    if [[ "$kind" == "gpu" ]]; then
        local tasks_per_node="$ngpus"
        local cpus_per_task=$((GPU_NODE_CPUS / tasks_per_node))
        if [[ "$cpus_per_task" -lt 1 ]]; then
            cpus_per_task=1
        fi
        dependency="$GPU_DEPENDENCY"
        sbatch_opts+=(--nodes=1 --gres="gpu:${GPU_TYPE}:${ngpus}" --ntasks="${ngpus}" --ntasks-per-node="${tasks_per_node}" --cpus-per-task="${cpus_per_task}" --mem="${mem_gb}")
    else
        dependency="$CPU_DEPENDENCY"
        sbatch_opts+=(--nodes=1 --ntasks="${CPU_RANKS}" --cpus-per-task=1 --mem="${mem_gb}")
    fi
    if [[ "${SCAVENGER}" -eq 1 ]]; then
        sbatch_opts+=(--partition=scavenger)
    fi

    depopt="$(dependency_opt "$dependency")"
    if [[ -n "$depopt" ]]; then
        sbatch_opts+=("$depopt")
    fi

    local sbatch_opts_text=""
    if [[ "${#sbatch_opts[@]}" -gt 0 ]]; then
        printf -v sbatch_opts_text ' %q' "${sbatch_opts[@]}"
    fi

    local cmd
    if [[ "$kind" == "gpu" ]]; then
        cmd="INPUT=${input} NGPUS=${ngpus} GPU_TYPE=${GPU_TYPE} sbatch${sbatch_opts_text} ${script}"
    else
        cmd="INPUT=${input} BUILD_CPU_IF_MISSING=0 sbatch${sbatch_opts_text} ${script}"
    fi
    printf '# %s
%s
' "$tag" "$cmd" >> "$OUT"
    if [[ "$SUBMIT" -eq 1 ]]; then
        if command -v sbatch >/dev/null 2>&1; then
            echo "submitting: $cmd"
            if [[ "$kind" == "gpu" ]]; then
                INPUT="${input}" NGPUS="${ngpus}" GPU_TYPE="${GPU_TYPE}"                     sbatch "${sbatch_opts[@]}" "${script}"
            else
                INPUT="${input}" BUILD_CPU_IF_MISSING=0                     sbatch "${sbatch_opts[@]}" "${script}"
            fi
        else
            echo "WARN: --submit requested but 'sbatch' not found; skipping: $cmd" >&2
        fi
    fi
}

weak_size_for_gpus() {
    case "$1" in
        1) echo 128 ;;
        2) echo 512 ;;
        4) echo 512 ;;
        8) echo 512 ;;
        *) echo 512 ;;
    esac
}

input_for_size() {
    case "$1" in
        128) echo "input_3d_centre_bore_128" ;;
        256) echo "input_3d_centre_bore" ;;
        512) echo "input_3d_centre_bore_512" ;;
        *) echo "input_3d_centre_bore_512" ;;
    esac
}

gpus_for_size() {
    case "$1" in
        512) echo "$LARGE_GPUS" ;;
        *) echo "$SMALL_GPUS" ;;
    esac
}

: > "$OUT"
{
    echo "# Phase-3 scaling sweep command matrix"
    echo "# generated by benchmark/phase3_scaling_sweep.sh"
    echo "# MODE=$MODE GPU_TYPE=$GPU_TYPE GPUS_PER_NODE=$GPUS_PER_NODE GPU_NODE_CPUS=$GPU_NODE_CPUS CPU_RANKS=$CPU_RANKS SIZES='$SIZES' GPUS='$GPUS'"
    echo "# inputs: input_3d_centre_bore family; scripts: nova_flame_*_3d*.slurm"
    echo "#"
} >> "$OUT"

if [[ "$MODE" == "strong" ]]; then
    for size in $SIZES; do
        input="$(input_for_size "$size")"
        for n in $(gpus_for_size "$size"); do
            if [[ "$n" -le 1 ]]; then
                emit "$size" "$input" "$n" "$SCRIPT_GPU_SINGLE" "strong size=${size} gpu x${n}" "gpu"
            else
                emit "$size" "$input" "$n" "$SCRIPT_GPU_MULTI" "strong size=${size} gpu x${n}" "gpu"
            fi
        done
        emit "$size" "$input" "$CPU_RANKS" "$SCRIPT_CPU" "strong size=${size} cpu-node baseline" "cpu"
    done
else
    for n in $GPUS; do
        size="$(weak_size_for_gpus "$n")"
        input="$(input_for_size "$size")"
        if [[ "$n" -le 1 ]]; then
            emit "$size" "$input" "$n" "$SCRIPT_GPU_SINGLE" "weak gpu x${n} size=${size}" "gpu"
        else
            emit "$size" "$input" "$n" "$SCRIPT_GPU_MULTI" "weak gpu x${n} size=${size}" "gpu"
        fi
    done
    big_size="$(weak_size_for_gpus 8)"
    emit "$big_size" "$(input_for_size "$big_size")" "$CPU_RANKS" "$SCRIPT_CPU" "weak cpu-node baseline size=${big_size}" "cpu"
fi

echo
echo "Wrote command matrix to: $OUT"
if [[ "$SUBMIT" -eq 0 ]]; then
    echo "(dry run -- nothing submitted; re-run with --submit on a SLURM host)"
fi

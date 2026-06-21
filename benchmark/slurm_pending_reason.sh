#!/usr/bin/env bash
#
# slurm_pending_reason.sh -- explain why pending Slurm jobs are waiting.
#
# Usage:
#   bash benchmark/slurm_pending_reason.sh
#   bash benchmark/slurm_pending_reason.sh <jobid> [<jobid> ...]
#   bash benchmark/slurm_pending_reason.sh --user <name>
#   bash benchmark/slurm_pending_reason.sh --all
#
# By default this prints all pending jobs for the current user. It does not
# submit, cancel, or modify anything.

set -euo pipefail

USER_NAME="${USER:-}"
SHOW_ALL=0
JOBIDS=()

usage() {
    cat <<'EOF'
Usage:
  bash benchmark/slurm_pending_reason.sh
  bash benchmark/slurm_pending_reason.sh <jobid> [<jobid> ...]
  bash benchmark/slurm_pending_reason.sh --user <name>
  bash benchmark/slurm_pending_reason.sh --all

Defaults:
  - No arguments: diagnose the current user's pending jobs.
  - Explicit job IDs: diagnose only those jobs.

This is read-only. It prints the scheduler reason, priority, dependency, and
resource request summary for each pending job.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        --user)
            shift
            if [[ $# -eq 0 ]]; then
                echo "error: --user requires a username" >&2
                exit 2
            fi
            USER_NAME="$1"
            ;;
        --all)
            SHOW_ALL=1
            ;;
        *)
            JOBIDS+=("$1")
            ;;
    esac
    shift
done

if ! command -v squeue >/dev/null 2>&1; then
    echo "error: squeue not found. Run this on a Slurm login node." >&2
    exit 1
fi

if [[ "${#JOBIDS[@]}" -eq 0 ]]; then
    if [[ "$SHOW_ALL" -eq 1 ]]; then
        mapfile -t JOBIDS < <(squeue -h -t PD -o '%i' 2>/dev/null | sort -u)
    else
        mapfile -t JOBIDS < <(squeue -h -u "$USER_NAME" -t PD -o '%i' 2>/dev/null | sort -u)
    fi
fi

if [[ "${#JOBIDS[@]}" -eq 0 ]]; then
    if [[ "$SHOW_ALL" -eq 1 ]]; then
        echo "No pending jobs found."
    else
        echo "No pending jobs found for user '$USER_NAME'."
    fi
    exit 0
fi

get_field() {
    local key="$1"
    local line="$2"
    local token
    for token in $line; do
        case "$token" in
            "${key}="*)
                printf '%s' "${token#*=}"
                return 0
                ;;
        esac
    done
    return 1
}

for jid in "${JOBIDS[@]}"; do
    echo "================================================================"
    echo "Job: $jid"

    summary="$(squeue -j "$jid" -h -o '  %i %t %P %j %u %M %D %R' 2>/dev/null || true)"
    if [[ -n "$summary" ]]; then
        echo "Summary:"
        echo "$summary"
    fi

    job_line="$(scontrol show job -o "$jid" 2>/dev/null || true)"
    if [[ -z "$job_line" ]]; then
        echo "  (scontrol show job returned no data)"
        continue
    fi

    reason="$(get_field Reason "$job_line" || true)"
    priority="$(get_field Priority "$job_line" || true)"
    partition="$(get_field Partition "$job_line" || true)"
    qos="$(get_field QOS "$job_line" || true)"
    dependency="$(get_field Dependency "$job_line" || true)"
    gres="$(get_field Gres "$job_line" || true)"
    tres_per_node="$(get_field TresPerNode "$job_line" || true)"
    num_nodes="$(get_field NumNodes "$job_line" || true)"
    num_cpus="$(get_field NumCPUs "$job_line" || true)"
    features="$(get_field Features "$job_line" || true)"
    nodelist="$(get_field NodeList "$job_line" || true)"
    exc_nodelist="$(get_field ExcNodeList "$job_line" || true)"
    submit_time="$(get_field SubmitTime "$job_line" || true)"
    eligible_time="$(get_field EligibleTime "$job_line" || true)"
    start_time="$(get_field StartTime "$job_line" || true)"

    echo "Details:"
    [[ -n "$reason" ]] && echo "  Reason        : $reason"
    [[ -n "$priority" ]] && echo "  Priority      : $priority"
    [[ -n "$partition" ]] && echo "  Partition     : $partition"
    [[ -n "$qos" ]] && echo "  QOS           : $qos"
    [[ -n "$dependency" ]] && echo "  Dependency    : $dependency"
    [[ -n "$gres" ]] && echo "  GRES          : $gres"
    [[ -n "$tres_per_node" ]] && echo "  TresPerNode   : $tres_per_node"
    [[ -n "$num_nodes" ]] && echo "  RequestedNodes: $num_nodes"
    [[ -n "$num_cpus" ]] && echo "  RequestedCPUs : $num_cpus"
    [[ -n "$features" ]] && echo "  Features      : $features"
    [[ -n "$nodelist" && "$nodelist" != "(null)" ]] && echo "  NodeList      : $nodelist"
    [[ -n "$exc_nodelist" ]] && echo "  ExcludedNodes : $exc_nodelist"
    [[ -n "$submit_time" ]] && echo "  SubmitTime    : $submit_time"
    [[ -n "$eligible_time" ]] && echo "  EligibleTime  : $eligible_time"
    [[ -n "$start_time" ]] && echo "  StartTime     : $start_time"

    case "$reason" in
        Dependency)
            echo "  Interpretation: waiting for an upstream job or explicit dependency to finish."
            ;;
        Priority)
            echo "  Interpretation: lower priority than jobs currently ahead in the queue."
            if command -v sprio >/dev/null 2>&1; then
                echo "  Priority breakdown:"
                sprio -j "$jid" 2>/dev/null || true
            fi
            ;;
        Resources)
            echo "  Interpretation: scheduler has not found a node/GPU allocation that fits yet."
            ;;
        ReqNodeNotAvail|ReqNodeNotAvail,_Reserved_For_High_Priority|BadConstraints)
            echo "  Interpretation: the requested node type/features are not currently available."
            ;;
        QOSMaxGRESPerUser|QOSMaxCpuPerUserLimit|QOSMaxGRESPerJob|AssocMaxGRESPerJob|AssocGrpGRES)
            echo "  Interpretation: blocked by a QOS/account resource limit."
            ;;
        JobHeldUser|JobHeldAdmin)
            echo "  Interpretation: the job is held."
            ;;
        *)
            echo "  Interpretation: see the reason string above; this is the scheduler's current blocker."
            ;;
    esac
done

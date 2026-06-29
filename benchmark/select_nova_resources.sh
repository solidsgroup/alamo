#!/usr/bin/env bash
# ============================================================================
# select_nova_resources.sh -- read-only NOVA slurm state probe + scheduling
#                              decision (GPU type x CPU x RAM).
#
# Source this from another script and call the function:
#   source benchmark/select_nova_resources.sh
#   select_nova_resources
#   echo "$SELECTED_GPU_TYPE $SELECTED_CPUS $SELECTED_MEM_GB $SELECTED_NODE"
#
# Or run standalone on the NOVA login node to see the decision without
# submitting anything:
#   bash benchmark/select_nova_resources.sh [--explain]
#
# Decision policy
# ----------------
# GPU candidates, most-preferred first: h200, a100. (a100-pcie is NEVER
# considered -- its CUDA compute capability on NOVA has not been confirmed,
# so we can't be sure which binary it would need.)
#
# GPU job resource profiles, most-preferred first: (PREF_CPUS, PREF_MEM_GB),
# (PREF_CPUS, MIN_MEM_GB), (MIN_CPUS, PREF_MEM_GB), (MIN_CPUS, MIN_MEM_GB).
# MIN_CPUS/MIN_MEM_GB are hard floors -- this script never selects below them.
#
# For each (gpu, cpu, mem) pair in that priority order, the script asks: is
# there a node in $PARTITION, in state IDLE or MIXED (not DRAIN/DOWN/MAINT),
# that RIGHT NOW has >=1 free GPU of that type, >=cpu free cores, and >=mem
# free RAM? The first match wins -- "available right now" beats "more
# preferred but would have to wait".
#
# If nothing matches at any profile on either GPU, the script falls back to
# the single most-preferred profile (h200 / PREF_CPUS / PREF_MEM_GB) and lets
# it queue normally in slurm -- not guaranteed to start immediately, but
# guaranteed to eventually run once resources free up.
#
# This module is READ-ONLY: it only calls sinfo/scontrol, never sbatch or
# scancel.
# ============================================================================

PARTITION="${PARTITION:-nova}"

MIN_CPUS="${MIN_CPUS:-8}"
PREF_CPUS="${PREF_CPUS:-8}"
MIN_MEM_GB="${MIN_MEM_GB:-32}"
PREF_MEM_GB="${PREF_MEM_GB:-64}"

# GPU candidates, most-preferred first. Allow callers to override the search
# order with NOVA_GPU_CANDIDATES_OVERRIDE="a100" or "h200 a100".
if [[ -n "${NOVA_GPU_CANDIDATES_OVERRIDE:-}" ]]; then
    read -r -a NOVA_GPU_CANDIDATES <<< "${NOVA_GPU_CANDIDATES_OVERRIDE}"
else
    NOVA_GPU_CANDIDATES=(h200 a100)
fi

# CPU/RAM profiles, most-preferred first.
NOVA_PROFILE_CPUS=("${PREF_CPUS}" "${PREF_CPUS}" "${MIN_CPUS}" "${MIN_CPUS}")
NOVA_PROFILE_MEM_GB=("${PREF_MEM_GB}" "${MIN_MEM_GB}" "${PREF_MEM_GB}" "${MIN_MEM_GB}")

_nova_log() { echo "$@" >&2; }

# Pull "KEY=value" out of a scontrol/squeue space-separated token line.
_nova_field() {
    local key="$1" line="$2" token
    for token in $line; do
        case "${token}" in
            "${key}="*) printf '%s' "${token#*=}"; return 0 ;;
        esac
    done
    return 1
}

# Extract the integer count for gres type $2 out of a gres/gresused string
# like "gpu:h200:4,gpu:a100:0" or "gpu:h200:2(IDX:0,1)". Prints 0 if absent.
_nova_gres_count() {
    local str="$1" gtype="$2" match
    if [[ -z "${str}" || "${str}" == "(null)" ]]; then
        echo 0
        return
    fi
    match="$(grep -oE "gpu:${gtype}:[0-9]+" <<< "${str}" | head -1)"
    if [[ -n "${match}" ]]; then
        echo "${match##*:}"
    else
        echo 0
    fi
}

# Extract an integer TRES count from a comma-separated CfgTRES/AllocTRES string.
# Supports keys like "gres/gpu" and "gres/gpu:h200". Prints 0 if absent.
_nova_tres_count() {
    local str="$1" key="$2" token
    if [[ -z "${str}" || "${str}" == "(null)" ]]; then
        echo 0
        return
    fi
    for token in ${str//,/ }; do
        case "${token}" in
            "${key}="*)
                printf '%s' "${token#*=}" | sed 's/[^0-9].*$//'
                return 0
                ;;
        esac
    done
    echo 0
}

# Is a scontrol State=... string usable right now (idle/mixed, not
# draining/down/maintenance/unresponsive)?
_nova_state_ok() {
    local state="${1:-}"
    case "${state}" in
        IDLE|IDLE+*)
            return 0
            ;;
        MIXED|MIXED+*)
            case "${state}" in
                *DRAIN*|*DOWN*|*MAINT*|*NOT_RESPONDING*|*FAIL*) return 1 ;;
                *) return 0 ;;
            esac
            ;;
        *)
            return 1
            ;;
    esac
}

# Populates NOVA_NODE_NAMES + the NOVA_NODE_* associative arrays below by
# probing every node in $PARTITION. Returns 1 (without raising, since this
# file is meant to be sourced) if sinfo is unavailable or returns nothing.
_nova_probe_nodes() {
    NOVA_NODE_NAMES=()
    declare -gA NOVA_NODE_STATE=()
    declare -gA NOVA_NODE_CPU_TOTAL=()
    declare -gA NOVA_NODE_CPU_ALLOC=()
    declare -gA NOVA_NODE_MEM_FREE_MB=()
    declare -gA NOVA_NODE_GPU_TOTAL=()   # key "node|gputype"
    declare -gA NOVA_NODE_GPU_USED=()    # key "node|gputype"

    if ! command -v sinfo >/dev/null 2>&1; then
        _nova_log "warn: sinfo not found -- not running on a Slurm login node?"
        return 1
    fi

    mapfile -t NOVA_NODE_NAMES < <(sinfo -p "${PARTITION}" -N -h -o "%N" 2>/dev/null | sort -u)
    if [[ "${#NOVA_NODE_NAMES[@]}" -eq 0 ]]; then
        _nova_log "warn: sinfo returned no nodes for partition '${PARTITION}'."
        return 1
    fi

    local node line state cputot cpualloc realmem allocmem freemem gres gresused cfgtres alloctres gtype total used
    for node in "${NOVA_NODE_NAMES[@]}"; do
        line="$(scontrol show node -o "${node}" 2>/dev/null || true)"
        [[ -z "${line}" ]] && continue

        state="$(_nova_field State "${line}")"
        cputot="$(_nova_field CPUTot "${line}")"
        cpualloc="$(_nova_field CPUAlloc "${line}")"
        realmem="$(_nova_field RealMemory "${line}")"
        allocmem="$(_nova_field AllocMem "${line}")"
        freemem="$(_nova_field FreeMem "${line}")"
        gres="$(_nova_field Gres "${line}")"
        gresused="$(_nova_field GresUsed "${line}")"
        cfgtres="$(_nova_field CfgTRES "${line}")"
        alloctres="$(_nova_field AllocTRES "${line}")"

        NOVA_NODE_STATE["${node}"]="${state}"
        NOVA_NODE_CPU_TOTAL["${node}"]="${cputot:-0}"
        NOVA_NODE_CPU_ALLOC["${node}"]="${cpualloc:-0}"

        # Prefer the scheduler's own FreeMem (most accurate); fall back to
        # RealMemory - AllocMem when FreeMem is not reported ("N/A" or unset).
        if [[ -n "${freemem}" && "${freemem}" != "N/A" && "${freemem}" =~ ^[0-9]+$ ]]; then
            NOVA_NODE_MEM_FREE_MB["${node}"]="${freemem}"
        elif [[ -n "${realmem}" && "${realmem}" =~ ^[0-9]+$ && -n "${allocmem}" && "${allocmem}" =~ ^[0-9]+$ ]]; then
            NOVA_NODE_MEM_FREE_MB["${node}"]=$(( realmem - allocmem ))
        else
            NOVA_NODE_MEM_FREE_MB["${node}"]=""
        fi

        for gtype in "${NOVA_GPU_CANDIDATES[@]}"; do
            total="$(_nova_tres_count "${cfgtres}" "gres/gpu:${gtype}")"
            if [[ "${total}" -eq 0 ]]; then
                total="$(_nova_gres_count "${gres}" "${gtype}")"
            fi

            used="$(_nova_tres_count "${alloctres}" "gres/gpu:${gtype}")"
            if [[ "${used}" -eq 0 ]]; then
                used="$(_nova_gres_count "${gresused}" "${gtype}")"
            fi
            if [[ "${used}" -eq 0 && "${alloctres}" != "(null)" ]]; then
                used="$(_nova_tres_count "${alloctres}" "gres/gpu")"
            fi
            NOVA_NODE_GPU_TOTAL["${node}|${gtype}"]="${total}"
            NOVA_NODE_GPU_USED["${node}|${gtype}"]="${used}"
        done
    done
    return 0
}

# Find the first node (in NOVA_NODE_NAMES order) that right now has >=1 free
# GPU of type $1, >=$2 free CPUs, and >=$3 GB free RAM, and is in a usable
# state. Prints the node name and returns 0 on a match, else returns 1.
_nova_find_node() {
    local gtype="$1" need_cpus="$2" need_mem_gb="$3"
    local need_mem_mb=$(( need_mem_gb * 1024 ))
    local node state cpu_total cpu_alloc cpu_free mem_free gpu_total gpu_used gpu_free

    for node in "${NOVA_NODE_NAMES[@]}"; do
        state="${NOVA_NODE_STATE[${node}]:-}"
        _nova_state_ok "${state}" || continue

        gpu_total="${NOVA_NODE_GPU_TOTAL[${node}|${gtype}]:-0}"
        [[ "${gpu_total}" -lt 1 ]] && continue
        gpu_used="${NOVA_NODE_GPU_USED[${node}|${gtype}]:-0}"
        gpu_free=$(( gpu_total - gpu_used ))
        [[ "${gpu_free}" -lt 1 ]] && continue

        cpu_total="${NOVA_NODE_CPU_TOTAL[${node}]:-0}"
        cpu_alloc="${NOVA_NODE_CPU_ALLOC[${node}]:-0}"
        cpu_free=$(( cpu_total - cpu_alloc ))
        [[ "${cpu_free}" -lt "${need_cpus}" ]] && continue

        mem_free="${NOVA_NODE_MEM_FREE_MB[${node}]:-}"
        [[ -z "${mem_free}" ]] && continue
        [[ "${mem_free}" -lt "${need_mem_mb}" ]] && continue

        echo "${node}"
        return 0
    done
    return 1
}

# Main entry point. Populates SELECTED_GPU_TYPE, SELECTED_CPUS,
# SELECTED_MEM_GB, SELECTED_NODE, SELECTED_IMMEDIATE (1 = a free node was
# found right now, 0 = falling back to queueing, -1 = manual FORCE_* override).
select_nova_resources() {
    SELECTED_GPU_TYPE=""
    SELECTED_CPUS=""
    SELECTED_MEM_GB=""
    SELECTED_NODE=""
    SELECTED_IMMEDIATE=0

    if [[ -n "${FORCE_GPU_TYPE:-}" || -n "${FORCE_CPUS:-}" || -n "${FORCE_MEM_GB:-}" ]]; then
        SELECTED_GPU_TYPE="${FORCE_GPU_TYPE:-h200}"
        SELECTED_CPUS="${FORCE_CPUS:-${PREF_CPUS}}"
        SELECTED_MEM_GB="${FORCE_MEM_GB:-${PREF_MEM_GB}}"
        SELECTED_NODE="(forced)"
        SELECTED_IMMEDIATE=-1
        _nova_log "select_nova_resources: FORCE_* override -> gpu=${SELECTED_GPU_TYPE} cpus=${SELECTED_CPUS} mem=${SELECTED_MEM_GB}G"
        return 0
    fi

    if ! _nova_probe_nodes; then
        SELECTED_GPU_TYPE="${NOVA_GPU_CANDIDATES[0]}"
        SELECTED_CPUS="${PREF_CPUS}"
        SELECTED_MEM_GB="${PREF_MEM_GB}"
        SELECTED_NODE="(unknown -- cluster probe failed)"
        SELECTED_IMMEDIATE=0
        _nova_log "select_nova_resources: cluster probe failed -- defaulting to preferred profile (${SELECTED_GPU_TYPE}/${PREF_CPUS}/${PREF_MEM_GB}G); will queue normally."
        return 0
    fi

    local gtype idx cpus mem_gb node
    for gtype in "${NOVA_GPU_CANDIDATES[@]}"; do
        for idx in "${!NOVA_PROFILE_CPUS[@]}"; do
            cpus="${NOVA_PROFILE_CPUS[$idx]}"
            mem_gb="${NOVA_PROFILE_MEM_GB[$idx]}"
            if node="$(_nova_find_node "${gtype}" "${cpus}" "${mem_gb}")"; then
                SELECTED_GPU_TYPE="${gtype}"
                SELECTED_CPUS="${cpus}"
                SELECTED_MEM_GB="${mem_gb}"
                SELECTED_NODE="${node}"
                SELECTED_IMMEDIATE=1
                _nova_log "select_nova_resources: immediate match -> ${gtype} on ${node} (${cpus} cpus, ${mem_gb}G free now)"
                return 0
            fi
        done
    done

    # Nothing immediately available anywhere -- queue at the single most
    # preferred profile and let slurm schedule it normally.
    SELECTED_GPU_TYPE="${NOVA_GPU_CANDIDATES[0]}"
    SELECTED_CPUS="${NOVA_PROFILE_CPUS[0]}"
    SELECTED_MEM_GB="${NOVA_PROFILE_MEM_GB[0]}"
    SELECTED_NODE="(none free now -- will queue)"
    SELECTED_IMMEDIATE=0
    _nova_log "select_nova_resources: no immediately-available node at any profile -- submitting preferred profile (${SELECTED_GPU_TYPE}/${SELECTED_CPUS}/${SELECTED_MEM_GB}G) to queue normally."
    return 0
}

# Human-readable dump of every node's current free CPU/RAM/GPU -- for
# debugging the decision above.
_nova_explain() {
    _nova_probe_nodes || { echo "Could not probe cluster state." >&2; return 1; }
    printf "%-16s %-12s %10s %10s %s\n" "NODE" "STATE" "CPU_FREE" "MEM_FREE_G" "GPU_FREE"
    local node state cpu_total cpu_alloc cpu_free mem_free_mb mem_free_g gpu_str gtype t u
    for node in "${NOVA_NODE_NAMES[@]}"; do
        state="${NOVA_NODE_STATE[${node}]:-?}"
        cpu_total="${NOVA_NODE_CPU_TOTAL[${node}]:-0}"
        cpu_alloc="${NOVA_NODE_CPU_ALLOC[${node}]:-0}"
        cpu_free=$(( cpu_total - cpu_alloc ))
        mem_free_mb="${NOVA_NODE_MEM_FREE_MB[${node}]:-}"
        if [[ -n "${mem_free_mb}" ]]; then
            mem_free_g=$(( mem_free_mb / 1024 ))
        else
            mem_free_g="?"
        fi
        gpu_str=""
        for gtype in "${NOVA_GPU_CANDIDATES[@]}"; do
            t="${NOVA_NODE_GPU_TOTAL[${node}|${gtype}]:-0}"
            [[ "${t}" -lt 1 ]] && continue
            u="${NOVA_NODE_GPU_USED[${node}|${gtype}]:-0}"
            gpu_str="${gpu_str}${gpu_str:+, }${gtype}:$((t - u))/${t}"
        done
        printf "%-16s %-12s %10s %10s %s\n" "${node}" "${state}" "${cpu_free}/${cpu_total}" "${mem_free_g}" "${gpu_str:-none}"
    done
}

# Standalone mode: only runs when this file is executed directly, not when
# sourced -- so sourcing it never mutates the caller's shell options.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    set -eo pipefail
    select_nova_resources
    if [[ "${1:-}" == "--explain" ]]; then
        echo
        echo "=== Full node availability (partition=${PARTITION}) ==="
        _nova_explain
        echo
    fi
    echo "SELECTED_GPU_TYPE=${SELECTED_GPU_TYPE}"
    echo "SELECTED_CPUS=${SELECTED_CPUS}"
    echo "SELECTED_MEM_GB=${SELECTED_MEM_GB}"
    echo "SELECTED_NODE=${SELECTED_NODE}"
    echo "SELECTED_IMMEDIATE=${SELECTED_IMMEDIATE}"
fi

#!/usr/bin/env bash
# ===========================================================================
# phase0_capture.sh -- local driver for benchmark/phase0_capture.slurm.
#
# Campaign PLAN §12 asks for "one command composes and submits a job; a second
# collects and renders on return". This is that pair. Hand-running the profiler
# at each phase boundary is the failure mode where the figure set gets captured
# twice and then abandoned.
#
# Everything here goes over the ControlMaster tunnel (~/.ssh/config Host `nova`).
# Nothing in this script runs a simulation locally -- kermit is correctness-only
# (campaign PLAN §2).
#
# ACTIONS
#   dryrun            Print every command the job would run. No ssh, no sbatch.
#   push              rsync the working tree's benchmark/ + src/ + decks to NOVA.
#   build             Submit the NOVA build job for BOTH binary variants.
#   submit [DIM ...]  Submit one capture job per dimension (default: 2 3).
#   status            squeue for this user, plus pending-reason for any held job.
#   collect <REMOTE>  rsync one remote artifact dir back under results/figures/.
#
# USAGE
#   bash benchmark/phase0_capture.sh dryrun
#   bash benchmark/phase0_capture.sh push
#   bash benchmark/phase0_capture.sh build
#   bash benchmark/phase0_capture.sh submit 2
#   bash benchmark/phase0_capture.sh status
#   bash benchmark/phase0_capture.sh collect /work/brunnels/jackplum/alamo/benchmark/_phase0_2d_a100_123456
# ===========================================================================
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

HOST="${HOST:-nova}"
REMOTE_DIR="${REMOTE_DIR:-/work/brunnels/jackplum/alamo}"
GPU_TYPE="${GPU_TYPE:-a100}"
LOCAL_RESULTS="${LOCAL_RESULTS:-${ROOT}/docs/agent_plans/20260727-phase0-baseline/results/figures}"

usage () { sed -n '2,32p' "${BASH_SOURCE[0]}"; exit "${1:-0}"; }
ACTION="${1:-}"; shift || true
[ -z "${ACTION}" ] && usage 0

ssh_nova () { ssh -o BatchMode=yes -o ConnectTimeout=15 "${HOST}" "$@"; }

case "${ACTION}" in

  dryrun)
    DRYRUN=1 DIM="${DIM:-2}" GPU_TYPE="${GPU_TYPE}" bash benchmark/phase0_capture.slurm
    ;;

  push)
    echo "=== rsync -> ${HOST}:${REMOTE_DIR}"
    rsync -az --info=stats1 \
      --exclude '.git' --exclude 'bin' --exclude 'ext' --exclude '__pycache__' \
      --exclude 'benchmark/_*' --exclude 'benchmark/baseline_runs' \
      ./benchmark ./src ./input* ./configure ./Makefile \
      "${HOST}:${REMOTE_DIR}/"
    echo "NOTE: bin/ and ext/ are deliberately not pushed -- build on NOVA."
    ;;

  build)
    echo "=== submitting NOVA build (both variants, both dims)"
    ssh_nova "cd ${REMOTE_DIR} && VARIANTS='profile plain' DIMS='2 3' ARCHES='80 90' sh benchmark/build_alamo_nova.sh"
    ;;

  submit)
    DIMS=("$@"); [ "${#DIMS[@]}" -eq 0 ] && DIMS=(2 3)
    for d in "${DIMS[@]}"; do
      echo "=== submit DIM=${d} GPU_TYPE=${GPU_TYPE}"
      ssh_nova "cd ${REMOTE_DIR} && sbatch --parsable --gres=gpu:${GPU_TYPE}:1 \
                  --export=ALL,DIM=${d},GPU_TYPE=${GPU_TYPE} \
                  benchmark/phase0_capture.slurm"
    done
    echo "Record the job ids. Do not block on the queue (campaign §2)."
    ;;

  status)
    ssh_nova 'squeue -u $USER -o "%.10i %.12j %.8T %.10M %.20R" 2>/dev/null'
    echo "--- pending reasons (if any)"
    ssh_nova "cd ${REMOTE_DIR} && bash benchmark/slurm_pending_reason.sh 2>/dev/null | head -20" || true
    ;;

  collect)
    REMOTE="${1:-}"; [ -z "${REMOTE}" ] && { echo "collect needs a remote artifact dir"; usage 2; }
    NAME="$(basename "${REMOTE}")"
    mkdir -p "${LOCAL_RESULTS}/${NAME}"
    echo "=== rsync ${HOST}:${REMOTE} -> ${LOCAL_RESULTS}/${NAME}"
    # Traces are large; CSVs and logs are what the figure set is built from.
    rsync -az --info=stats1 "${HOST}:${REMOTE}/" "${LOCAL_RESULTS}/${NAME}/"
    echo "=== collected:"
    find "${LOCAL_RESULTS}/${NAME}" -maxdepth 2 -type f | head -40
    ;;

  *) echo "unknown action '${ACTION}'"; usage 2 ;;
esac

#!/usr/bin/env bash
# ============================================================================
# collect_phaseA_to_local.sh -- tar the Phase-A NOVA outputs and scp them home.
#
# Run this ON NOVA from the alamo repo root after the A1 (ncu/nsys diag) and A2
# (combined crossover) jobs have finished. It bundles:
#   - benchmark/profiles_alpha1/   (ncu_*/ and nsys_*/ reports + CSVs -- A1)
#   - flame_gpu_3d_diag.*.{out,err} (A1 diag SLURM logs; .err carries ncu/nsys
#                                    tool errors -- e.g. why a profile came out empty)
#   - flame_gpu_3d_a2.*.{out,err}   (A2 GPU crossover logs: TinyProfiler + time -v)
#   - flame_cpu_3d_a2.*.{out,err}   (A2 CPU baseline logs)
# into a timestamped tarball, then scp's it to the local box's Downloads dir.
#
# Usage:
#   bash benchmark/collect_phaseA_to_local.sh            # tar + scp
#   bash benchmark/collect_phaseA_to_local.sh --no-send  # tar only, print path
#
# Env knobs:
#   DEST_USER=jackplum            ssh user on the local box
#   DEST_HOST=10.24.220.162       local box IP
#   DEST_DIR=Downloads            target dir (relative to remote $HOME)
#   ROOT_DIR=$PWD                 alamo repo root (default: cwd)
#   OUT=phaseA_outputs_<ts>.tgz   tarball name
# ============================================================================
set -euo pipefail

DEST_USER="${DEST_USER:-jackplum}"
DEST_HOST="${DEST_HOST:-10.24.220.162}"
DEST_DIR="${DEST_DIR:-Downloads}"
ROOT_DIR="${ROOT_DIR:-$PWD}"
TS="$(date +%Y%m%d_%H%M%S)"
OUT="${OUT:-phaseA_outputs_${TS}.tgz}"

SEND=1
case "${1:-}" in
    --no-send) SEND=0 ;;
    -h|--help) sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    "") ;;
    *) echo "error: unknown argument '$1' (try --help)" >&2; exit 2 ;;
esac

cd "${ROOT_DIR}"
if [ ! -d benchmark ]; then
    echo "error: run from the alamo repo root (no ./benchmark here: ${ROOT_DIR})" >&2
    exit 1
fi

# Collect existing targets only; nullglob so missing patterns don't get tarred
# as literal names.
shopt -s nullglob
TARGETS=()
[ -d benchmark/profiles_alpha1 ] && TARGETS+=("benchmark/profiles_alpha1")
TARGETS+=( flame_gpu_3d_diag.*.out flame_gpu_3d_diag.*.err \
           flame_gpu_3d_a2.*.out   flame_gpu_3d_a2.*.err \
           flame_cpu_3d_a2.*.out   flame_cpu_3d_a2.*.err )
shopt -u nullglob

if [ "${#TARGETS[@]}" -eq 0 ]; then
    echo "error: nothing to collect -- no profiles_alpha1/ dir and no *_a2/diag .out logs found in ${ROOT_DIR}" >&2
    exit 1
fi

echo "=== Phase-A collection ==="
echo "  repo      = ${ROOT_DIR}"
echo "  tarball   = ${OUT}"
echo "  including :"
printf '    %s\n' "${TARGETS[@]}"

tar czf "${OUT}" "${TARGETS[@]}"
echo "  wrote $(du -h "${OUT}" | cut -f1) -> ${ROOT_DIR}/${OUT}"

if [ "${SEND}" -eq 0 ]; then
    echo "(--no-send: tarball written, nothing sent)"
    exit 0
fi

echo "=== scp -> ${DEST_USER}@${DEST_HOST}:${DEST_DIR}/ ==="
scp "${OUT}" "${DEST_USER}@${DEST_HOST}:${DEST_DIR}/"
echo "=== done: ${DEST_DIR}/${OUT} on ${DEST_HOST} ==="
echo "    unpack locally with:  tar xzf ~/${DEST_DIR}/${OUT}"

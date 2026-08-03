#!/usr/bin/env bash
# benchmark/ledger_append.sh -- append one row to benchmark/metrics_ledger.csv.
#
# PLAN v2 §14.2: one CSV, appended at every phase boundary, keyed
#   (phase, commit_sha, tree_hash, case, n_ranks, device, build_config)
# with columns T0 through T7 plus wall time. Purpose is regression detection
# ACROSS phases, which is why the key includes tree_hash: a row keyed on
# commit_sha alone cannot identify a dirty-tree measurement (§5.1), and the
# first baseline was captured from 574 dirty files.
#
# Unmeasured targets are written as the literal `NA`, never as 0 and never
# blank. A blank reads as zero in every spreadsheet; a zero for T2 page faults
# or T4 allocations is a PASS. Recording an unrun target as a pass is the
# false-pass failure the whole target set was rewritten to prevent.
#
# USAGE
#   PHASE=0 CASE=input_copy DEVICE=a100 RANKS=1 BUILD=plain-cuda80 \
#     T0=0.3525 T0_SD=0.0041 WALL=31.726 \
#     bash benchmark/ledger_append.sh
#
#   Provenance (commit_sha, tree_hash) is read from benchmark/_pushed_rev.txt
#   when present, so a row always matches the tree that was actually captured.
#   Override with SHA= / TREE_HASH= for a local (unpushed) measurement.
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"

LEDGER="${LEDGER:-benchmark/metrics_ledger.csv}"
REV="benchmark/_pushed_rev.txt"

[ -f "${LEDGER}" ] || { echo "missing ${LEDGER}" >&2; exit 2; }

read_rev () { awk -F= -v k="$1" '$1==k{print $2; exit}' "${REV}" 2>/dev/null; }

SHA="${SHA:-$(read_rev local_head)}"
SHA="${SHA:-$(git rev-parse HEAD)}"
TREE_HASH="${TREE_HASH:-$(read_rev tree_hash)}"
TREE_HASH="${TREE_HASH:-UNPUSHED}"

: "${PHASE:?set PHASE}"
: "${CASE:?set CASE}"
DEVICE="${DEVICE:-unknown}"
RANKS="${RANKS:-1}"
BUILD="${BUILD:-unknown}"

# Every target defaults to NA. See the header comment on why not 0.
for v in T0 T0_SD T1 T2 T3A T3B T3C T4 T5A T5B T5C T6A T6B T7 WALL; do
  eval "${v}=\"\${${v}:-NA}\""
done
NOTES="${NOTES:-}"

if [ "${TREE_HASH}" = "UNPUSHED" ]; then
  echo "WARNING: no tree_hash (benchmark/_pushed_rev.txt absent)." >&2
  echo "  The row will not identify the measured source. Run 'phase0_capture.sh push' first" >&2
  echo "  or pass TREE_HASH= explicitly." >&2
fi

printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
  "${PHASE}" "${SHA}" "${TREE_HASH}" "${CASE}" "${RANKS}" "${DEVICE}" "${BUILD}" \
  "${T0}" "${T0_SD}" "${T1}" "${T2}" "${T3A}" "${T3B}" "${T3C}" "${T4}" \
  "${T5A}" "${T5B}" "${T5C}" "${T6A}" "${T6B}" "${T7}" "${WALL}" \
  "\"${NOTES}\"" >> "${LEDGER}"

echo "appended to ${LEDGER}:"
tail -1 "${LEDGER}"

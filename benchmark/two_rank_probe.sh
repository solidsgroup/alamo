#!/usr/bin/env bash
# two_rank_probe.sh -- Phase 0.5 multi-rank chamber-reduction correctness probe.
#
# WHY THIS EXISTS
# ---------------
# The chamber model reduces regression rate over the burning surface into a single
# scalar that drives a pressure ODE (Flame.cpp:1113-1130 -> Integrator.cpp:1271-1275
# -> Flame.cpp:706). A source read says that reduction is globally allreduced before
# the ODE advances. This script turns that read into runtime evidence: run the same
# deck at -np 1 and -np 2 and require the chamber scalar history to agree.
#
# If it does NOT agree, every multi-rank chamber result is physically wrong and the
# campaign's Phase 2 device-resident-scalar redesign is being built on a broken
# contract. That is why this runs before any measurement.
#
# SCOPE: correctness only. This produces no admissible timing number -- kermit's
# A1000 is shared and power-capped (campaign PLAN 20260727-gpu-memory-strategy §2).
#
# CONCURRENCY: writes only to benchmark/_two_rank_probe_*/ and never touches
# benchmark/baseline_runs/ or bin/. Still: do not run this while status.sh or
# ci_golden_compare.sh is in flight -- they rebuild bin/ underneath a live mpiexec.
#
# USAGE
#   bash benchmark/two_rank_probe.sh cpu 10
#   bash benchmark/two_rank_probe.sh gpu 10
#   DRYRUN=1 bash benchmark/two_rank_probe.sh gpu 10     # print commands, run nothing
#   DECK=input_rod_and_tube_2d REL_TOL=1e-5 bash benchmark/two_rank_probe.sh cpu 5
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

usage () { sed -n '2,28p' "${BASH_SOURCE[0]}"; exit "${1:-0}"; }
case "${1:-}" in -h|--help|help) usage 0 ;; esac

BUILD="${1:-cpu}"
MAX_STEP="${2:-10}"
DECK="${DECK:-input}"
NP_HI="${NP_HI:-2}"
REL_TOL="${REL_TOL:-1e-6}"
ABS_TOL="${ABS_TOL:-1e-8}"
DRYRUN="${DRYRUN:-0}"
# The campaign deck sets no amr.max_grid_size, so AMReX's 2D default (128) leaves
# level 0 as a single 64x64 box that lives entirely on rank 0. Force a multi-box
# layout so both ranks actually own coarse-level work -- otherwise the probe can
# pass vacuously. Same reasoning as MGS in benchmark/local_a100_gate.sh.
MGS="${MGS:-32}"
ARENA_INIT_SIZE="${ARENA_INIT_SIZE:-1073741824}"   # 1 GiB/rank; default 3/4-of-8GB OOMs at init
OUT="${OUT:-$ROOT/benchmark/_two_rank_probe_$(date +%Y%m%d_%H%M%S)_${BUILD}}"

case "$BUILD" in
  cpu) BIN="${BIN:-bin/alamo-2d-g++}" ;;
  gpu) BIN="${BIN:-bin/alamo_gpu-2d-cuda86-g++}"
       # shellcheck source=/dev/null
       [ -f benchmark/local_cuda_env.sh ] && . benchmark/local_cuda_env.sh ;;
  *)   echo "ERROR: build must be 'cpu' or 'gpu' (got '$BUILD')"; usage 2 ;;
esac

[ -x "$BIN" ]  || { echo "ERROR: missing binary $BIN (build it first)"; exit 2; }
[ -f "$DECK" ] || { echo "ERROR: missing deck $DECK"; exit 2; }

# Deck overrides mirror benchmark/baseline_suite.py:184-198 so this probe and the
# golden gate exercise the same configuration. thermo.int=1 forces the chamber
# reduction every step, which is the quantity under test.
ARGS=( "$DECK"
       "max_step=$MAX_STEP"
       "stop_time=1e99_s"
       "amr.plot_int=-1"
       "amr.thermo.plot_int=1"
       "amr.thermo.int=1"
       "amr.max_grid_size=$MGS"
       "elastic.solver.verbose=0"
       "elastic.print_model=0" )
[ "$BUILD" = gpu ] && ARGS+=( "amrex.the_arena_init_size=$ARENA_INIT_SIZE"
                              "amrex.abort_on_out_of_gpu_memory=1" )
# EXTRA: space-separated deck overrides appended last, e.g. EXTRA="elastic.interval=1"
# to force the elastic solve inside a short horizon (the campaign deck's default
# elastic.interval=50 means a 10-step probe never enters MLMG at all).
# shellcheck disable=SC2206
[ -n "${EXTRA:-}" ] && ARGS+=( ${EXTRA} )

run_np () {   # $1 = rank count
  local np="$1"
  local dir="$OUT/np${np}"
  mkdir -p "$dir/plot"
  local cmd=( mpiexec -np "$np" "$ROOT/$BIN" "${ARGS[@]}" "plot_file=$dir/plot" )
  printf '  np=%s: %s\n' "$np" "${cmd[*]}"
  [ "$DRYRUN" = 1 ] && return 0
  local t0 t1 rc
  t0=$(date +%s)
  "${cmd[@]}" >"$dir/run.log" 2>&1
  rc=$?
  t1=$(date +%s)
  if [ $rc -ne 0 ]; then
    echo "  FAIL np=$np exited $rc (wall=$((t1-t0))s) -- see $dir/run.log"
    tail -n 15 "$dir/run.log"
    return 1
  fi
  [ -f "$dir/plot/thermo.dat" ] || { echo "  FAIL np=$np produced no thermo.dat"; return 1; }
  echo "  ok   np=$np wall=$((t1-t0))s  rows=$(( $(wc -l <"$dir/plot/thermo.dat") - 1 ))"
}

mkdir -p "$OUT"
echo "=== two-rank probe: build=$BUILD deck=$DECK max_step=$MAX_STEP np=1 vs np=$NP_HI ==="
echo "    binary=$BIN"
echo "    out=$OUT"
run_np 1     || exit 1
run_np "$NP_HI" || exit 1

if [ "$DRYRUN" = 1 ]; then
  echo "=== DRYRUN: nothing executed, no comparison performed ==="
  exit 0
fi

echo "=== compare np=1 vs np=$NP_HI (abs_tol=$ABS_TOL rel_tol=$REL_TOL) ==="
python3 benchmark/compare_thermo.py \
  "$OUT/np1/plot/thermo.dat" "$OUT/np${NP_HI}/plot/thermo.dat" \
  --abs-tol "$ABS_TOL" --rel-tol "$REL_TOL" | tee "$OUT/compare.txt"
rc="${PIPESTATUS[0]}"

if [ "$rc" -eq 0 ]; then
  echo "PROBE PASS -- chamber history is rank-count invariant at rel_tol=$REL_TOL"
else
  echo "PROBE FAIL (exit $rc) -- do NOT loosen the tolerance to make this pass."
  echo "  Record the first diverging column in the task folder RESULT.md and stop."
fi
echo "artifacts: $OUT"
exit "$rc"

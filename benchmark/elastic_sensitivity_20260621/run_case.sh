#!/usr/bin/env bash
# Elastic-solver sensitivity harness (GPU break-vs-sensitivity study, 2026-06-21).
# Runs ONE elastic solve (max_step=1, interval=1, tstart=0) and prints a one-line
# summary so the orchestrator's context stays clean.
#
# Usage: run_case.sh LABEL BIN NP [OVERRIDE ...]
#   LABEL    : short tag -> log_<LABEL>.log, out_<LABEL>
#   BIN      : path to alamo binary
#   NP       : mpi ranks (1 = no mpiexec)
#   OVERRIDE : each remaining arg is one ParmParse override, passed verbatim as a
#              single argv token. Multi-value ones MUST be one quoted arg, e.g.
#              'amr.n_cell=64 64 64'
set -uo pipefail
ROOT=/home/jackplum/Projects/alamo
cd "$ROOT"
source "$ROOT/benchmark/local_cuda_env.sh" 2>/dev/null || true

LABEL="$1"; BIN="$2"; NP="$3"; shift 3
EXTRA=("$@")
DIR="$ROOT/benchmark/elastic_sensitivity_20260621"
LOG="$DIR/log_${LABEL}.log"

if [ "$NP" -gt 1 ]; then CMD=(mpiexec -np "$NP" "$BIN"); else CMD=("$BIN"); fi

START=$(date +%s)
timeout 900 "${CMD[@]}" "$DIR/input_base" \
  max_step=1 stop_time=1e99_s \
  elastic.interval=1 elastic.tstart=0.0 elastic.solver.verbose=4 elastic.print_model=0 \
  amr.plot_int=-1 amr.thermo.plot_int=-1 amr.thermo.int=1 \
  plot_file="$DIR/out_${LABEL}" \
  "${EXTRA[@]}" > "$LOG" 2>&1
RC=$?
END=$(date +%s)

# ---- classify ----
PARSE_ERR=$(grep -iE 'too many values|ParmParse.*not.*found|unused inputs|table has unused' "$LOG" | head -1)
DIVERGE=$(grep -iE 'MLMG: Failing|failing so|diverg' "$LOG" | head -1)
CUDAERR=$(grep -iE 'unspecified launch|CUDA error|illegal memory|out of memory' "$LOG" | head -1)
SIGABRT=$(grep -iE 'amrex::Abort::SIGABRT|signal 11|segfault|Segmentation' "$LOG" | head -1)
LASTRESID=$(grep -oiE 'resid/bnorm =[ ]*[0-9.eE+-]+' "$LOG" | tail -1)
FINAL=$(grep -E 'Final Iter' "$LOG" | tail -2 | tr '\n' '|')
NRES=$(grep -cE 'Final Iter' "$LOG")   # number of MLMG solves that reached "Final Iter"

VERDICT="UNKNOWN"
if [ -n "$PARSE_ERR" ]; then VERDICT="CONFIG_ERROR";
elif [ -n "$CUDAERR" ]; then VERDICT="CUDA_FAULT";
elif [ -n "$DIVERGE" ]; then VERDICT="DIVERGED";
elif [ "$RC" = "0" ]; then VERDICT="CONVERGED(exit0)";
elif [ "$RC" = "124" ]; then VERDICT="TIMEOUT";
elif [ -n "$SIGABRT" ]; then VERDICT="ABORT(other)";
fi

printf '%-30s rc=%-4s t=%ss  verdict=%s  (Final-Iter blocks: %s)\n' "$LABEL" "$RC" "$((END-START))" "$VERDICT" "$NRES"
printf '    last_resid: %s\n' "${LASTRESID:-<none>}"
printf '    final_iter: %s\n' "${FINAL:-<none>}"
[ -n "$PARSE_ERR" ] && printf '    PARSE: %s\n' "$PARSE_ERR"
[ -n "$DIVERGE" ]  && printf '    DIVERGE: %s\n' "$DIVERGE"
[ -n "$CUDAERR" ]  && printf '    CUDA: %s\n' "$CUDAERR"
exit 0

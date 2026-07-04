#!/usr/bin/env bash
# =============================================================================
# Phase 3 -- microarchitectural efficiency via `perf stat`.
# -----------------------------------------------------------------------------
# IPC, stalled cycles, cache + TLB miss rates and branch misprediction tell you
# WHY the CPU run costs what it does (memory-bound vs compute-bound vs
# front-end bound).  For the GPU binary we stat the HOST side: a healthy GPU
# offload shows low host IPC and high "task-clock" idle (host waiting on device).
#
# perf needs kernel.perf_event_paranoid <= 2 for counters and <= 1 for kernel
# samples.  This box is at 4, so we try; on failure we print the one-liner to
# relax it (run via `! sudo sysctl ...` from the Claude prompt) and continue.
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

if ! have perf; then warn "perf not installed -- skipping phase 3"; exit 0; fi

PARANOID="$(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo 4)"
info "perf_event_paranoid = $PARANOID"
if [ "$PARANOID" -gt 2 ]; then
  warn "paranoid > 2: hardware counters likely blocked."
  warn "To enable:  ! sudo sysctl kernel.perf_event_paranoid=1"
fi

EVENTS="task-clock,context-switches,cpu-migrations,page-faults,cycles,instructions,branches,branch-misses,cache-references,cache-misses,stalled-cycles-frontend,stalled-cycles-backend"

stat_one() {
  local kind="$1" bin="$2" env_prefix="$3"
  if [ ! -x "$bin" ]; then warn "skip $kind: $bin missing"; return; fi
  local out="$RESULTS_DIR/perfstat_${kind}.txt"
  info "perf stat: $kind (steps=$PROFILE_STEPS)"
  rm -rf "$RESULTS_DIR/plt_${kind}_ps"*
  # shellcheck disable=SC2086
  eval "$env_prefix perf stat -o '$out' -e $EVENTS -- \
    '$bin' '$INPUT' max_step=$PROFILE_STEPS \
      amr.plot_int=-1 amr.thermo.plot_int=-1 \
      plot_file='$RESULTS_DIR/plt_${kind}_ps' \
      elastic.solver.verbose=0 ${COMMON_OVERRIDES[*]}" \
    >"$RESULTS_DIR/perfstat_${kind}.runlog" 2>&1
  if grep -q "<not supported>\|<not counted>\|Permission" "$out" 2>/dev/null; then
    warn "$kind: some counters unavailable (paranoid/perms)."
  fi
  [ -s "$out" ] && ok "wrote $out" || warn "no perf output for $kind"
}

stat_one cpu "$BIN_CPU" ""
GPU_ENV_PREFIX=""
[ -f "$GPU_ENV" ] && GPU_ENV_PREFIX="source '$GPU_ENV' >/dev/null 2>&1;"
stat_one gpu "$BIN_GPU" "$GPU_ENV_PREFIX"

py "$ANALYSIS_DIR/lib/parse_perfstat.py" \
   "$RESULTS_DIR/perfstat_cpu.txt" "$RESULTS_DIR/perfstat_gpu.txt" "$RESULTS_DIR"
ok "phase 3 (perf stat) complete"

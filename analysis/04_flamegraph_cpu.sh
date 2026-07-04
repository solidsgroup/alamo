#!/usr/bin/env bash
# =============================================================================
# Phase 4 -- CPU flame graph (perf record -> FlameGraph SVG).
# -----------------------------------------------------------------------------
# Produces an interactive flame graph of where CPU cycles are spent so you can
# see, at a glance, the MLMG / phase-field / I/O call stacks and their width.
#
#   perf record -F <freq> --call-graph dwarf  (dwarf works without frame ptrs;
#       the binaries are -O2 so fp may be omitted -- dwarf unwinding is safest)
#   perf script | stackcollapse-perf.pl | flamegraph.pl  -> SVG
#
# FlameGraph is auto-vendored from GitHub (shallow clone) if absent.
# Sampling needs perf_event_paranoid <= 1; we detect and instruct on failure.
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

if ! have perf; then warn "perf missing -- skipping CPU flame graph"; exit 0; fi

# --- vendor FlameGraph -------------------------------------------------------
FG_DIR="$VENDOR_DIR/FlameGraph"
if [ ! -x "$FG_DIR/flamegraph.pl" ]; then
  info "vendoring FlameGraph -> $FG_DIR"
  if have git && git clone --depth 1 https://github.com/brendangregg/FlameGraph \
        "$FG_DIR" >/dev/null 2>&1; then
    ok "FlameGraph cloned"
  else
    warn "could not fetch FlameGraph (offline?). Will still emit folded stacks."
  fi
fi

PARANOID="$(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo 4)"
if [ "$PARANOID" -gt 1 ]; then
  warn "perf_event_paranoid=$PARANOID blocks stack sampling."
  warn "Enable with:  ! sudo sysctl kernel.perf_event_paranoid=1"
  warn "(also: ! sudo sysctl kernel.kptr_restrict=0 for kernel symbols)"
fi

flame_one() {
  local kind="$1" bin="$2" env_prefix="$3"
  if [ ! -x "$bin" ]; then warn "skip $kind flame: $bin missing"; return; fi
  local data="$RESULTS_DIR/perf_${kind}.data"
  local folded="$RESULTS_DIR/flame_${kind}.folded"
  local svg="$RESULTS_DIR/flamegraph_${kind}.svg"
  info "perf record: $kind (freq=${PERF_FREQ}Hz, steps=$PROFILE_STEPS)"
  rm -rf "$RESULTS_DIR/plt_${kind}_fg"*
  # shellcheck disable=SC2086
  eval "$env_prefix perf record -F $PERF_FREQ --call-graph dwarf -o '$data' -- \
    '$bin' '$INPUT' max_step=$PROFILE_STEPS \
      amr.plot_int=$PLOT_INT amr.thermo.plot_int=$PLOT_INT \
      plot_file='$RESULTS_DIR/plt_${kind}_fg' \
      elastic.solver.verbose=0 ${COMMON_OVERRIDES[*]}" \
    >"$RESULTS_DIR/perf_${kind}.runlog" 2>&1
  if [ ! -s "$data" ]; then
    warn "$kind: perf produced no data (paranoid/perms). See runlog."
    return
  fi
  perf script -i "$data" 2>/dev/null \
    | "$FG_DIR/stackcollapse-perf.pl" 2>/dev/null >"$folded"
  if [ ! -s "$folded" ]; then warn "$kind: empty folded stacks"; return; fi
  ok "folded stacks: $folded ($(wc -l <"$folded") unique stacks)"
  if [ -x "$FG_DIR/flamegraph.pl" ]; then
    "$FG_DIR/flamegraph.pl" --title "alamo $kind flame graph (void=20MPa)" \
        --width 1400 --colors hot "$folded" >"$svg" 2>/dev/null
    [ -s "$svg" ] && ok "flame graph: $svg" || warn "$kind: flamegraph.pl failed"
  fi
  # quick top-of-stack hotspot summary (works even without FlameGraph)
  awk '{n=split($0,a," ");c=a[n];sub(/^.* /,"");
        m=split($0,b,";");print b[m]" "$NF}' "$folded" 2>/dev/null \
    | sort | uniq -c | sort -rn | head -20 \
    > "$RESULTS_DIR/flame_${kind}_top.txt" 2>/dev/null || true
}

flame_one cpu "$BIN_CPU" ""
# GPU host-side flame graph (shows CUDA launch/sync stacks). Real device kernel
# flame charts require nsys (not installed); see README.
GPU_ENV_PREFIX=""
[ -f "$GPU_ENV" ] && GPU_ENV_PREFIX="source '$GPU_ENV' >/dev/null 2>&1;"
flame_one gpu "$BIN_GPU" "$GPU_ENV_PREFIX"

ok "phase 4 (flame graphs) complete"

#!/usr/bin/env bash
# =============================================================================
# run_all.sh -- orchestrate the full CPU/GPU performance analysis suite.
# -----------------------------------------------------------------------------
# Phases:
#   1  wall-clock & efficiency   (post-hoc, from completed run logs)   [always]
#   2  I/O profile               (strace, short instrumented re-runs)
#   3  perf stat                 (microarchitecture counters)
#   4  CPU flame graph           (perf record -> FlameGraph SVG)
#   5  GPU timeline              (nvidia-smi sampling / nsys if present)
#   6  consolidated report       (REPORT.md + index.html)             [always]
#
# Usage:
#   analysis/run_all.sh [--profile-steps N] [--plot-int N]
#                       [--skip-io] [--skip-perf] [--skip-flame] [--skip-gpu]
#                       [--phases 1,2,6]
#
# Env overrides honored: BIN_CPU BIN_GPU INPUT CPU_LOG GPU_LOG PROFILE_STEPS
#                        PLOT_INT PERF_FREQ SAMPLE_MS RESULTS_DIR
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

PHASES="1 2 3 4 5 6"
while [ $# -gt 0 ]; do
  case "$1" in
    --profile-steps) export PROFILE_STEPS="$2"; shift 2;;
    --plot-int)      export PLOT_INT="$2"; shift 2;;
    --phases)        PHASES="${2//,/ }"; shift 2;;
    --skip-io)       PHASES="${PHASES/2/}"; shift;;
    --skip-perf)     PHASES="${PHASES/3/}"; shift;;
    --skip-flame)    PHASES="${PHASES/4/}"; shift;;
    --skip-gpu)      PHASES="${PHASES/5/}"; shift;;
    -h|--help)       grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *) err "unknown arg: $1"; exit 2;;
  esac
done

hr
info "ALAMO performance suite"
info "CPU bin : $BIN_CPU"
info "GPU bin : $BIN_GPU"
info "input   : $INPUT"
info "CPU log : $CPU_LOG  $([ -f "$CPU_LOG" ] && echo '(found)' || echo '(MISSING)')"
info "GPU log : $GPU_LOG  $([ -f "$GPU_LOG" ] && echo '(found)' || echo '(MISSING)')"
info "results : $RESULTS_DIR"
info "phases  : $PHASES   profile_steps=$PROFILE_STEPS plot_int=$PLOT_INT"
hr

want() { [[ " $PHASES " == *" $1 "* ]]; }

if want 1; then
  hr; info "PHASE 1: wall-clock & efficiency"
  py "$ANALYSIS_DIR/01_wallclock.py" "$CPU_LOG" "$GPU_LOG" "$RESULTS_DIR" || warn "phase 1 issues"
fi
if want 2; then hr; info "PHASE 2: I/O profile"; bash "$ANALYSIS_DIR/02_io_profile.sh" || warn "phase 2 issues"; fi
if want 3; then hr; info "PHASE 3: perf stat";   bash "$ANALYSIS_DIR/03_perf_stat.sh" || warn "phase 3 issues"; fi
if want 4; then hr; info "PHASE 4: CPU flame";   bash "$ANALYSIS_DIR/04_flamegraph_cpu.sh" || warn "phase 4 issues"; fi
if want 5; then hr; info "PHASE 5: GPU timeline";bash "$ANALYSIS_DIR/05_gpu_timeline.sh" || warn "phase 5 issues"; fi
if want 6; then
  hr; info "PHASE 6: consolidated report"
  py "$ANALYSIS_DIR/06_report.py" "$RESULTS_DIR" || warn "phase 6 issues"
fi

hr
ok "suite complete. Artifacts in: $RESULTS_DIR"
ok "open report: $RESULTS_DIR/index.html"
ls -1 "$RESULTS_DIR" | sed 's/^/    /'

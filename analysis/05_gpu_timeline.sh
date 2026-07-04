#!/usr/bin/env bash
# =============================================================================
# Phase 5 -- GPU device behaviour.
# -----------------------------------------------------------------------------
# Two tiers, best available is used:
#
#   (A) nsys (Nsight Systems), if installed -> full kernel/API timeline +
#       `nsys stats` summaries (kernel time, memops, CUDA API, gaps). This is
#       the real GPU flame-chart source. Not installed on this box.
#
#   (B) Fallback: sample `nvidia-smi` at high rate in the background while the
#       GPU job runs, capturing SM utilization, memory utilization, memory used,
#       power and SM clock as a time series -> CSV + SVG timeline. Combined with
#       the ioctl time from phase 2 this characterizes device occupancy and the
#       host<->device sync overhead that dominates small-kernel workloads.
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

if [ ! -x "$BIN_GPU" ]; then warn "GPU binary missing -- skipping phase 5"; exit 0; fi
[ -f "$GPU_ENV" ] && source "$GPU_ENV" >/dev/null 2>&1

# ---- Tier A: nsys ----------------------------------------------------------
# Extra per-run overrides may be passed via $NSYS_EXTRA_ARGS (space-separated),
# e.g. to change grid/elastic params for a specific comparison.
if have_nsys; then
  info "nsys found ($NSYS_BIN) -- capturing full GPU timeline"
  REP="$RESULTS_DIR/gpu_nsys"
  # --sample/--cpuctxsw none: this box has perf_event_paranoid>=2 (no CPU
  # sampling perm); CUDA API+kernel trace needs no special privilege.
  "$NSYS_BIN" profile -o "$REP" --force-overwrite true --stats=false \
    --trace=cuda,nvtx --sample=none --cpuctxsw=none \
    "$BIN_GPU" "$INPUT" max_step="$PROFILE_STEPS" \
      amr.plot_int="$PLOT_INT" amr.thermo.plot_int="$PLOT_INT" \
      plot_file="$RESULTS_DIR/plt_gpu_nsys" \
      elastic.solver.verbose=0 "${COMMON_OVERRIDES[@]}" ${NSYS_EXTRA_ARGS:-} \
    >"$RESULTS_DIR/gpu_nsys.runlog" 2>&1
  "$NSYS_BIN" stats --force-export=true \
    --report cuda_api_sum --report cuda_gpu_kern_sum --report cuda_gpu_mem_time_sum \
    --format csv --output "$RESULTS_DIR/gpu_nsys_stats" \
    "$REP.nsys-rep" >>"$RESULTS_DIR/gpu_nsys.runlog" 2>&1 || true
  ok "nsys report: $REP.nsys-rep (open in Nsight Systems GUI for flame timeline)"
fi

# ---- Tier B: nvidia-smi sampling fallback ----------------------------------
info "sampling nvidia-smi during GPU run (steps=$PROFILE_STEPS, plot_int=$PLOT_INT)"
CSV="$RESULTS_DIR/gpu_timeline.csv"
SAMPLE_MS="${SAMPLE_MS:-100}"
rm -rf "$RESULTS_DIR/plt_gpu_tl"*

# background sampler (epoch ms, util.gpu, util.mem, mem.used MiB, power W, sm clk)
nvidia-smi --query-gpu=timestamp,utilization.gpu,utilization.memory,memory.used,power.draw,clocks.sm,temperature.gpu \
  --format=csv,nounits -lms "$SAMPLE_MS" > "$CSV" 2>/dev/null &
SMI_PID=$!

"$BIN_GPU" "$INPUT" max_step="$PROFILE_STEPS" \
  amr.plot_int="$PLOT_INT" amr.thermo.plot_int="$PLOT_INT" \
  plot_file="$RESULTS_DIR/plt_gpu_tl" \
  elastic.solver.verbose=0 "${COMMON_OVERRIDES[@]}" \
  > "$RESULTS_DIR/gpu_timeline.runlog" 2>&1
RC=$?

kill "$SMI_PID" 2>/dev/null; wait "$SMI_PID" 2>/dev/null
[ $RC -ne 0 ] && warn "GPU run rc=$RC (see gpu_timeline.runlog)"
[ -s "$CSV" ] && ok "timeline CSV: $CSV" || warn "no nvidia-smi samples captured"

py "$ANALYSIS_DIR/lib/parse_gpu_timeline.py" "$CSV" "$RESULTS_DIR"
ok "phase 5 (GPU timeline) complete"

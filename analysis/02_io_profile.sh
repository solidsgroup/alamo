#!/usr/bin/env bash
# =============================================================================
# Phase 2 -- I/O profile via strace syscall accounting.
# -----------------------------------------------------------------------------
# `strace -f -c` aggregates wall time spent INSIDE each syscall across all
# threads.  That is the most direct, privilege-free measurement of "time locked
# in I/O" available on this box (no eBPF/biosnoop access at paranoid=4).
#
# We run a SHORT job with plotting ENABLED so the VisMF plotfile write path
# (open/write/fsync/close/mkdir/rename) is actually exercised -- a plotting-off
# run has essentially no I/O and would mislead.
#
# The parser classifies syscalls into: file-write, file-read, sync, metadata,
# gpu-driver (ioctl), memory (mmap/munmap/brk) and other; then reports the
# I/O-attributable fraction of in-kernel time for CPU vs GPU.
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

run_one() {
  local kind="$1" bin="$2" extra_env="$3"
  if [ ! -x "$bin" ]; then warn "skip $kind: binary not found ($bin)"; return; fi
  info "strace I/O profile: $kind  (steps=$PROFILE_STEPS, plot_int=$PLOT_INT)"
  local out="$RESULTS_DIR/strace_${kind}.txt"
  local runlog="$RESULTS_DIR/strace_${kind}.runlog"
  rm -rf "$RESULTS_DIR/plt_${kind}_io"*
  # %file,%desc => path- and fd-based syscalls; add ioctl for GPU driver time.
  # shellcheck disable=SC2086
  ( eval "$extra_env" \
    strace -f -c -e "trace=%file,%desc,ioctl,mmap,munmap,brk" -o "$out" \
      "$bin" "$INPUT" \
        max_step="$PROFILE_STEPS" \
        amr.plot_int="$PLOT_INT" amr.thermo.plot_int="$PLOT_INT" \
        plot_file="$RESULTS_DIR/plt_${kind}_io" \
        elastic.solver.verbose=0 "${COMMON_OVERRIDES[@]}" \
    ) >"$runlog" 2>&1
  local rc=$?
  if [ $rc -ne 0 ]; then warn "$kind run exited rc=$rc (see $runlog)"; fi
  if [ -s "$out" ]; then ok "wrote $out"; else warn "no strace summary for $kind"; fi
}

run_one cpu "$BIN_CPU" ""
# GPU run needs the CUDA runtime environment sourced into the child.
GPU_ENV_PREFIX=""
if [ -f "$GPU_ENV" ]; then GPU_ENV_PREFIX="source '$GPU_ENV' >/dev/null 2>&1;"; fi
run_one gpu "$BIN_GPU" "$GPU_ENV_PREFIX"

info "parsing strace summaries"
py "$ANALYSIS_DIR/lib/parse_strace.py" \
   "$RESULTS_DIR/strace_cpu.txt" "$RESULTS_DIR/strace_gpu.txt" "$RESULTS_DIR"
ok "phase 2 (I/O) complete"

#!/usr/bin/env bash
# =============================================================================
# analysis/lib/common.sh  --  shared configuration & helpers for the perf suite
# -----------------------------------------------------------------------------
# Sourced by every script in analysis/.  Holds the single source of truth for
# binary paths, the workload definition, output locations, tool detection and
# logging primitives.  Override anything from the environment, e.g.
#
#     PROFILE_STEPS=120 PLOT_INT=20 ./analysis/run_all.sh
#
# Design notes (architect):
#   * The "production" comparison runs (full 1.5 s, plotting OFF) are measured
#     post-hoc from their /usr/bin/time -v footer in the .log files.
#   * The "instrumented" runs are SHORT, deliberately re-executed under
#     perf / strace / nvidia-smi.  They intentionally ENABLE plotting so that
#     the I/O path (VisMF writes + fsync) is actually exercised and timed --
#     a zero-I/O run tells you nothing about I/O.
#   * Everything degrades gracefully: a missing tool downgrades one artifact,
#     never aborts the suite.
# =============================================================================
set -uo pipefail

# --- locate repo root (analysis/ lives directly under it) --------------------
COMMON_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(dirname "$COMMON_DIR")"
ROOT_DIR="$(dirname "$ANALYSIS_DIR")"
cd "$ROOT_DIR"

# --- workload definition (override via env) ---------------------------------
BIN_CPU="${BIN_CPU:-$ROOT_DIR/bin/alamo-2d-clang++}"
BIN_GPU="${BIN_GPU:-$ROOT_DIR/bin/alamo_gpu-2d-cuda86-g++}"
INPUT="${INPUT:-$ROOT_DIR/input_copy}"
GPU_ENV="${GPU_ENV:-$ROOT_DIR/benchmark/local_cuda_env.sh}"

# Resolve nsys: explicit $NSYS env > project-local extracted toolkit > PATH.
# The box has no system nsys; Phase 5 extracts Nsight Systems CLI under
# .local/nsight (dpkg-deb -x, no root). Newest target-linux-x64 wins.
if [ -z "${NSYS_BIN:-}" ]; then
  if [ -n "${NSYS:-}" ] && [ -x "${NSYS:-}" ]; then
    NSYS_BIN="$NSYS"
  else
    NSYS_BIN="$(ls -1 "$ROOT_DIR"/.local/nsight/opt/nvidia/nsight-systems/*/target-linux-x64/nsys 2>/dev/null | sort -V | tail -1)"
    [ -z "$NSYS_BIN" ] && command -v nsys >/dev/null 2>&1 && NSYS_BIN="$(command -v nsys)"
  fi
fi
export NSYS_BIN
have_nsys() { [ -n "${NSYS_BIN:-}" ] && [ -x "${NSYS_BIN:-}" ]; }

# Physics overrides shared by every run so CPU and GPU solve the IDENTICAL
# problem.  Star geometry (input_copy native) is used because the centre-bore
# (simple_circle) variant fails MLMG at void=20 MPa.
COMMON_OVERRIDES=(
  model_void.kappa=20_MPa
  model_void.mu=20_MPa
  elastic.print_model=0
)

# Production-run artifacts (already executed at full stop_time, plotting off)
CPU_LOG="${CPU_LOG:-$ROOT_DIR/out_cpu_star_20mpa.log}"
GPU_LOG="${GPU_LOG:-$ROOT_DIR/out_gpu_star_20mpa.log}"

# Instrumented short-run knobs
PROFILE_STEPS="${PROFILE_STEPS:-80}"   # >50 so >=1 elastic solve is captured
PLOT_INT="${PLOT_INT:-40}"             # enable plotting for the I/O-aware runs
PERF_FREQ="${PERF_FREQ:-997}"          # Hz, prime to dodge lockstep aliasing

# --- output locations --------------------------------------------------------
RESULTS_DIR="${RESULTS_DIR:-$ANALYSIS_DIR/results}"
VENDOR_DIR="$ANALYSIS_DIR/vendor"
mkdir -p "$RESULTS_DIR" "$VENDOR_DIR"

# --- logging -----------------------------------------------------------------
if [ -t 1 ]; then
  C_RST=$'\e[0m'; C_DIM=$'\e[2m'; C_RED=$'\e[31m'; C_GRN=$'\e[32m'
  C_YEL=$'\e[33m'; C_BLU=$'\e[34m'; C_BLD=$'\e[1m'
else
  C_RST=; C_DIM=; C_RED=; C_GRN=; C_YEL=; C_BLU=; C_BLD=
fi
info()  { printf '%s[*]%s %s\n'  "$C_BLU" "$C_RST" "$*"; }
ok()    { printf '%s[+]%s %s\n'  "$C_GRN" "$C_RST" "$*"; }
warn()  { printf '%s[!]%s %s\n'  "$C_YEL" "$C_RST" "$*" >&2; }
err()   { printf '%s[x]%s %s\n'  "$C_RED" "$C_RST" "$*" >&2; }
hr()    { printf '%s%s%s\n' "$C_DIM" "------------------------------------------------------------" "$C_RST"; }
have()  { command -v "$1" >/dev/null 2>&1; }

# Run a single-rank instrumented job.  We invoke the binary DIRECTLY (no
# mpiexec) so the profiler attaches to the solver, not the MPI launcher.
# $1 = binary, $2 = run tag (output subdir / plot prefix), rest = extra args.
alamo_short_cmd() {
  local bin="$1" tag="$2"; shift 2
  printf '%s %s max_step=%s amr.plot_int=%s amr.thermo.plot_int=%s plot_file=%s %s %s' \
    "$bin" "$INPUT" "$PROFILE_STEPS" "$PLOT_INT" "$PLOT_INT" \
    "$RESULTS_DIR/plt_${tag}" "${COMMON_OVERRIDES[*]}" "$*"
}

py() { python3 "$@"; }

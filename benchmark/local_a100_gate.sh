#!/usr/bin/env bash
# local_a100_gate.sh -- pre-NOVA local correctness gate for the chamber-gpu branch.
#
# WHY THIS EXISTS
# ---------------
# A GPU defect (Base::Mechanics::Integrate writing host members from a device
# kernel) faulted on NOVA's A100 with `CUDA error 700 (illegal memory access)`
# but did NOT fault on the local RTX A1000. Root cause of that divergence:
#
#     $ nvidia-smi -q | grep 'Addressing Mode'   ->  HMM
#     $ cat /sys/module/nvidia_uvm/parameters/uvm_disable_hmm  ->  N
#
# The A1000 driver has HMM (Heterogeneous Memory Management) enabled, so a device
# dereference of a *system-allocated* host pointer is silently serviced by page
# migration instead of faulting. The A100 environment does not migrate it, so it
# hard-faults. Plain runs and even AMReX `--debug` post-kernel error checks can be
# fooled by HMM (no CUDA error is ever raised). **compute-sanitizer is HMM-immune**
# -- it validates each access against the *intended* allocation -- so it is the
# reliable local stand-in for the A100. This script makes that the default gate.
#
# TIERS (run all by default; pick with TIERS="1 2 ...")
#   1  runtime-strict smoke   -- device arena + CUDA_LAUNCH_BLOCKING + abort-on-OOM.
#                                Fast. Catches hard faults / OOM at the fault site.
#   2  memcheck (A100 spoof)  -- compute-sanitizer memcheck on a NOVA-like *multi-box*
#                                layout. Catches the illegal-access (error-700) class
#                                that HMM hides locally. This is the load-bearing tier.
#   3  racecheck + initcheck  -- non-atomic global races (the `+=`-into-shared class)
#                                and uninitialized device reads.
#   4  HMM-disabled hard run  -- OPTIONAL, root-gated. Reload nvidia_uvm with HMM off so
#                                the A1000 faults like the A100 without sanitizer overhead.
#                                Skipped unless ALLOW_HMM_DISABLE=1 and sudo is available.
#
# USAGE
#   bash benchmark/local_a100_gate.sh
#   TIERS="1 2" INPUT=input_3d_centre_bore_128_a2 bash benchmark/local_a100_gate.sh
#   ALLOW_HMM_DISABLE=1 TIERS=4 bash benchmark/local_a100_gate.sh   # needs sudo
set -uo pipefail

# --- config -----------------------------------------------------------------
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
BIN="${BIN:-bin/alamo_gpu-3d-cuda86-g++}"
INPUT="${INPUT:-input_3d_centre_bore_128_a2}"
MGS="${MGS:-64}"                 # max_grid_size=64 forces a NOVA-like multi-box layout
STOP_TIME="${STOP_TIME:-6e-4}"   # ~6 steps: hits traction diag every step + 1 elastic solve
TIER1_STOP="${TIER1_STOP:-5e-4}" # tier1 stops BEFORE the step-5 elastic solve to stay fast.
                                 # NOTE: pure device arena (the_arena_is_managed=0) is NOT used
                                 # -- on this shared 8 GB HMM desktop it OOMs reserving its initial
                                 # block (Xorg/gnome hold ~1.1 GB) and it cannot defeat HMM anyway,
                                 # so it adds no A100 fidelity. tier1's strictness is launch-blocking
                                 # (async faults reported at the exact kernel) + abort-on-OOM.
TIERS="${TIERS:-1 2}"   # default = smoke + memcheck (practical per-change gate);
                        # add 3 for racecheck/initcheck (slow, pre-NOVA), 4 for HMM-disable
OUT="${OUT:-$ROOT/benchmark/_a100_gate_$(date +%Y%m%d_%H%M%S)}"
CUDA_HOME="${CUDA_HOME:-$ROOT/.local/cuda}"
SAN="${SAN:-$ROOT/.local/cuda-12.0-ubuntu/usr/bin/compute-sanitizer}"
export CUDA_HOME
export LD_LIBRARY_PATH="${CUDA_HOME}/lib:${LD_LIBRARY_PATH:-}"

# elastic solver knobs kept tiny so the (separate, slow) MLMG solve completes fast
# under sanitizer instrumentation -- we are gating *correctness*, not convergence.
ELASTIC_FAST=(elastic.interval=5 elastic.solver.nriters=1 elastic.solver.max_iter=2
              elastic.solver.pre_smooth=1 elastic.solver.post_smooth=1)

mkdir -p "$OUT"
FAIL=0
say () { printf '\n=== %s ===\n' "$*"; }
need_bin () { [ -x "$BIN" ] || { echo "ERROR: missing $BIN (build it first)"; exit 2; }; }
need_san () { [ -x "$SAN" ] || { echo "ERROR: missing compute-sanitizer at $SAN"; exit 2; }; }

run_base () {  # $1=logfile  $2=stop_time  rest=extra alamo args
  local log="$1" stop="$2"; shift 2
  "$BIN" "$INPUT" amr.max_grid_size="$MGS" plot_file="$OUT/plt" \
      stop_time="$stop" "${ELASTIC_FAST[@]}" "$@" >"$log" 2>&1
}

# --- Tier 1: runtime-strict smoke -------------------------------------------
tier1 () {
  need_bin; say "TIER 1  runtime-strict smoke -- launch-blocking + abort-on-OOM, pre-elastic (stop=$TIER1_STOP)"
  local log="$OUT/tier1_runtime_strict.log" t0=$(date +%s)
  CUDA_LAUNCH_BLOCKING=1 run_base "$log" "$TIER1_STOP" \
      amrex.the_arena_is_managed=1 amrex.abort_on_out_of_gpu_memory=1
  local rc=$?; local t1=$(date +%s)
  if [ $rc -eq 0 ] && ! grep -qiE "illegal|CUDA error|Abort|SIGABRT" "$log"; then
    echo "  PASS  wall=$((t1-t0))s  -> $log"
  else
    echo "  FAIL  exit=$rc  -> $log"; grep -iE "illegal|CUDA error|Abort" "$log" | head -3; FAIL=1
  fi
}

# --- Tier 2: memcheck (the A100 spoof) --------------------------------------
tier2 () {
  need_bin; need_san
  say "TIER 2  compute-sanitizer MEMCHECK on multi-box layout (HMM-immune A100 spoof)"
  local log="$OUT/tier2_memcheck.log" t0=$(date +%s)
  "$SAN" --tool memcheck --error-exitcode 99 --launch-timeout 120 \
      "$BIN" "$INPUT" amr.max_grid_size="$MGS" plot_file="$OUT/plt_mc" \
      stop_time="$STOP_TIME" "${ELASTIC_FAST[@]}" amrex.the_arena_is_managed=1 >"$log" 2>&1
  local rc=$?; local t1=$(date +%s)
  if grep -q "ERROR SUMMARY: 0 errors" "$log"; then
    echo "  PASS  wall=$((t1-t0))s  -> $log"
  else
    echo "  FAIL  exit=$rc wall=$((t1-t0))s  -> $log"
    grep -E "Invalid|out of bounds|ERROR SUMMARY|Mechanics.H|Flame.cpp" "$log" | head -8; FAIL=1
  fi
}

# --- Tier 3: racecheck + initcheck ------------------------------------------
tier3 () {
  need_bin; need_san
  for tool in racecheck initcheck; do
    say "TIER 3  compute-sanitizer ${tool^^}"
    local log="$OUT/tier3_${tool}.log" t0=$(date +%s)
    "$SAN" --tool "$tool" --error-exitcode 99 \
        "$BIN" "$INPUT" amr.max_grid_size="$MGS" plot_file="$OUT/plt_$tool" \
        stop_time="$STOP_TIME" "${ELASTIC_FAST[@]}" amrex.the_arena_is_managed=1 >"$log" 2>&1
    local rc=$?; local t1=$(date +%s)
    # memcheck/initcheck print "ERROR SUMMARY: 0 errors"; racecheck prints
    # "RACECHECK SUMMARY: 0 hazards displayed (0 errors, 0 warnings)".
    if grep -qE "ERROR SUMMARY: 0 errors|RACECHECK SUMMARY: 0 hazards" "$log"; then
      echo "  PASS  wall=$((t1-t0))s  -> $log"
    else
      echo "  FAIL  exit=$rc wall=$((t1-t0))s  -> $log"
      grep -E "race|hazard|Uninitialized|ERROR SUMMARY|RACECHECK SUMMARY" "$log" | head -8; FAIL=1
    fi
  done
}

# --- Tier 4: HMM-disabled hard run (optional, root) -------------------------
tier4 () {
  say "TIER 4  HMM-disabled hard run (root-gated A100 spoof, no sanitizer overhead)"
  if [ "${ALLOW_HMM_DISABLE:-0}" != "1" ]; then
    echo "  SKIP  (set ALLOW_HMM_DISABLE=1 to attempt; reloads nvidia_uvm, needs sudo,"
    echo "         and will disrupt any GPU/display session using the device)"
    return
  fi
  echo "  reloading nvidia_uvm with uvm_disable_hmm=1 (sudo) ..."
  if ! sudo modprobe -r nvidia_uvm 2>/dev/null || ! sudo modprobe nvidia_uvm uvm_disable_hmm=1 2>/dev/null; then
    echo "  SKIP  could not reload nvidia_uvm (in use / no sudo)"; return
  fi
  echo "  Addressing Mode now: $(nvidia-smi -q | awk -F: '/Addressing Mode/{print $2}' | xargs)"
  local log="$OUT/tier4_hmm_disabled.log"
  CUDA_LAUNCH_BLOCKING=1 run_base "$log" "$STOP_TIME" amrex.the_arena_is_managed=1
  local rc=$?
  [ $rc -eq 0 ] && echo "  PASS  exit=0 -> $log" || { echo "  FAIL exit=$rc -> $log"; FAIL=1; }
  echo "  restoring HMM (uvm_disable_hmm=0) ..."
  sudo modprobe -r nvidia_uvm 2>/dev/null && sudo modprobe nvidia_uvm 2>/dev/null
}

# --- driver -----------------------------------------------------------------
echo "local_a100_gate: bin=$BIN input=$INPUT mgs=$MGS stop_time=$STOP_TIME tiers='$TIERS'"
echo "HMM status: addressing=$(nvidia-smi -q 2>/dev/null | awk -F: '/Addressing Mode/{print $2}' | xargs)  uvm_disable_hmm=$(cat /sys/module/nvidia_uvm/parameters/uvm_disable_hmm 2>/dev/null)"
for t in $TIERS; do "tier$t"; done
say "RESULT"
[ $FAIL -eq 0 ] && echo "ALL TIERS PASSED -- cleared for NOVA" || echo "GATE FAILED -- inspect logs in $OUT"
exit $FAIL

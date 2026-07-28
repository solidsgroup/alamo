#!/usr/bin/env bash
# benchmark/status.sh -- sole source of chamber-gpu branch state.
#
# Replaces prose status docs (docs/archive/CURRENT.md). Run at the start of
# every session per docs/llm/CONVENTIONS.md.
#
# TWO MODES. The distinction is load-bearing; read it before quoting a PASS.
#
#   default (fast)   Per-session orientation. Cheap legs only. Every label
#                    below names the leg it actually runs, because the previous
#                    version printed "a100-sanitizer: PASS" for a run that
#                    executed no sanitizer (campaign PLAN v2 section 5.2).
#   FULL=1           The correctness gate. Tier 2 compute-sanitizer memcheck +
#                    gpu_strict golden compare. Required before any commit that
#                    touches src/, and at every phase exit.
#
# WHAT THE FAST LEGS DO AND DO NOT COVER
# --------------------------------------
#   device-lint            Static pattern lint. Full coverage in both modes.
#   golden-cpu             ci_golden_compare.sh GOLDEN_MODE=cpu. CPU binary
#                          only. Says NOTHING about the GPU path.
#   smoke-pre-elastic      local_a100_gate.sh TIERS=1. Launch-blocking +
#                          abort-on-OOM smoke, MANAGED memory, and it stops
#                          before the step-5 elastic solve (TIER1_STOP=5e-4).
#                          No compute-sanitizer. No elastic. ~2s.
#
# The two slow legs are the ones that catch the defect classes this branch
# actually produces: the error-700 illegal-access class that local HMM hides
# (Tier 2 memcheck), and GPU-vs-CPU numerical divergence (gpu_strict).
#
#   golden-gpu-strict      ci_golden_compare.sh GOLDEN_MODE=gpu. Needs nvcc.
#   a100-memcheck-tier2    local_a100_gate.sh TIERS="1 2". HMM-immune.
#
# COVERAGE STILL MISSING FROM BOTH MODES (campaign PLAN v2 section 5.2 requires
# these; they are not built yet, so do not read a FULL=1 PASS as covering them):
#   - pressure history over enough steps to diverge
#   - at least one regrid
#   - two ranks
#   - at least one checkpoint/restart cycle
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"

FULL="${FULL:-0}"

echo "== chamber-gpu status $(date -Is) =="
echo "branch:  $(git branch --show-current)"
echo "HEAD:    $(git log -1 --format='%h %ci %s')"
echo "dirty:   $(git status --porcelain | wc -l) files"
echo

if [ "${FULL}" = "1" ]; then
  echo "== gates (FULL -- correctness gate) =="
else
  echo "== gates (fast -- orientation only, NOT a correctness gate) =="
fi

run_gate () {
  local name="$1"; shift
  if [ ! -x "$1" ]; then echo "$name: MISSING ($1)"; return; fi
  local log="benchmark/_gate_logs/$(basename "$1" .sh).log"
  mkdir -p benchmark/_gate_logs
  if "$@" >"$log" 2>&1; then echo "$name: PASS"; else echo "$name: FAIL (see $log)"; fi
}

run_gate "device-lint          " benchmark/lint_device_patterns.sh

# Export, rather than prefixing the run_gate call: run_gate is a shell function
# and command-prefix assignments in front of a function are not reliably
# exported to the processes it spawns.
if [ "${FULL}" = "1" ]; then
  export GOLDEN_MODE=gpu
  export TIERS="1 2"
  # ci_golden_compare.sh GOLDEN_MODE=gpu reconfigures and rebuilds (deliberately
  # -- it exists to defeat stale-binary selection) and hard-requires nvcc on
  # PATH. Report a missing toolchain as BLOCKED, not FAIL: FAIL reads as a
  # correctness result, and an unrun leg is not a result at all.
  if ! command -v nvcc >/dev/null 2>&1; then
    for c in .local/cuda/bin .local/cuda-12.6.3-redist/bin; do
      [ -x "$c/nvcc" ] && { PATH="$PWD/$c:$PATH"; export PATH; break; }
    done
  fi
  if command -v nvcc >/dev/null 2>&1; then
    run_gate "golden-gpu-strict    " benchmark/ci_golden_compare.sh
  else
    echo "golden-gpu-strict    : BLOCKED (nvcc not on PATH; leg did not run)"
  fi
  run_gate "a100-memcheck-tier2  " benchmark/local_a100_gate.sh
else
  export GOLDEN_MODE=cpu
  export TIERS=1
  run_gate "golden-cpu           " benchmark/ci_golden_compare.sh
  run_gate "smoke-pre-elastic    " benchmark/local_a100_gate.sh
  echo
  echo "  NOT RUN: golden-gpu-strict, a100-memcheck-tier2."
  echo "  These legs are the correctness gate. A green line above is orientation,"
  echo "  not assurance -- run 'FULL=1 bash benchmark/status.sh' before committing"
  echo "  src/ changes or exiting a phase."
fi

echo
echo "== open task folders =="
ls -d docs/agent_plans/*/ 2>/dev/null | while read -r d; do
  if [ ! -f "${d}results/DONE" ]; then echo "OPEN: $d"; fi
done

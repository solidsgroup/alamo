#!/usr/bin/env bash
# benchmark/status.sh -- sole source of chamber-gpu branch state.
#
# Replaces prose status docs (docs/archive/CURRENT.md). Run at the start of
# every session per docs/llm/CONVENTIONS.md.
#
# Gate invocations below are the real fast-path each script supports, not
# placeholders:
#   device-lint      benchmark/lint_device_patterns.sh   (Phase 3; MISSING until then)
#   golden-compare   benchmark/ci_golden_compare.sh, GOLDEN_MODE=cpu (its default;
#                    no separate fast flag needed -- CPU leg + 3 tiny golden cases
#                    + NaN-smoke measured at ~1m36s wall on kermit, well under budget)
#   a100-sanitizer   benchmark/local_a100_gate.sh, TIERS=1 only (its own "Tier 1"
#                    is documented as the fast runtime-strict smoke, ~2s wall;
#                    Tier 2+ memcheck is the slow load-bearing gate, run separately
#                    before NOVA submission per the script's own header)
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"
echo "== chamber-gpu status $(date -Is) =="
echo "branch:  $(git branch --show-current)"
echo "HEAD:    $(git log -1 --format='%h %ci %s')"
echo "dirty:   $(git status --porcelain | wc -l) files"
echo
echo "== gates =="
run_gate () {
  local name="$1"; shift
  if [ ! -x "$1" ]; then echo "$name: MISSING ($1)"; return; fi
  local log="benchmark/_gate_logs/$(basename "$1" .sh).log"
  mkdir -p benchmark/_gate_logs
  if "$@" >"$log" 2>&1; then echo "$name: PASS"; else echo "$name: FAIL (see $log)"; fi
}
export GOLDEN_MODE="${GOLDEN_MODE:-cpu}"
export TIERS="${TIERS:-1}"
run_gate "device-lint     " benchmark/lint_device_patterns.sh
run_gate "golden-compare  " benchmark/ci_golden_compare.sh
run_gate "a100-sanitizer  " benchmark/local_a100_gate.sh
echo
echo "== open task folders =="
ls -d docs/agent_plans/*/ 2>/dev/null | while read -r d; do
  if [ ! -f "${d}results/DONE" ]; then echo "OPEN: $d"; fi
done

#!/usr/bin/env bash
# record_references.sh -- Phase 1 task 1.E, local leg.
#
# Runs the full local case set (cpu + gpu_strict profiles) and copies the
# resulting bundles into the NAMED reference locations
# benchmark/validate/references/{cpu_strict,gpu_alpha1_local}/, then runs
# compare_validation.py to establish BASELINE_DELTAS.md -- the current
# GPU-vs-CPU physics-error baseline (what passes today, before any v3
# optimization work lands; see physics_budget.md's review checklist, which
# explicitly wants these real deltas to sanity-check the ENGINEERING
# tolerances before 1.F makes them gating).
#
# References get OVERWRITTEN each time this runs -- they represent "the
# current state", not a history (bulky per-run artifacts already accumulate
# under benchmark/validate/runs/, which is where history lives).
#
# NOVA leg (the roadmap's "alpha-1.0 GPU bundle" on an actual A100) is NOT
# done by this script -- it needs benchmark/validate/run_validation_nova.slurm
# submitted by the user (cluster access), then its bundle copied in beside
# gpu_alpha1_local as e.g. references/gpu_alpha1_nova_a100/.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

PROFILES="${PROFILES:-cpu,gpu_strict}"
REFS_DIR="$ROOT/benchmark/validate/references"

echo "=== running local validation (profiles=$PROFILES) ==="
python3 benchmark/validate/run_validation_local.py --profiles "$PROFILES" \
  --runs-dir "$ROOT/benchmark/validate/runs"

LATEST_CPU="$(ls -dt "$ROOT"/benchmark/validate/runs/*_cpu 2>/dev/null | head -1 || true)"
LATEST_GPU_STRICT="$(ls -dt "$ROOT"/benchmark/validate/runs/*_a1000_sm*_strict 2>/dev/null | head -1 || true)"

if [ -z "$LATEST_CPU" ] || [ -z "$LATEST_GPU_STRICT" ]; then
  echo "ERROR: could not find a fresh cpu and gpu_strict bundle under runs/" >&2
  exit 1
fi

echo "cpu bundle:        $LATEST_CPU"
echo "gpu_strict bundle: $LATEST_GPU_STRICT"

rm -rf "$REFS_DIR/cpu_strict" "$REFS_DIR/gpu_alpha1_local"
mkdir -p "$REFS_DIR"
cp -r "$LATEST_CPU" "$REFS_DIR/cpu_strict"
cp -r "$LATEST_GPU_STRICT" "$REFS_DIR/gpu_alpha1_local"

echo "=== comparing references/cpu_strict vs references/gpu_alpha1_local ==="
RC=0
python3 benchmark/validate/compare_validation.py \
  "$REFS_DIR/cpu_strict" "$REFS_DIR/gpu_alpha1_local" \
  --budget benchmark/validate/physics_budget.yaml \
  --out-md "$REFS_DIR/../BASELINE_DELTAS.md" \
  --out-json "$REFS_DIR/../BASELINE_DELTAS.json" || RC=$?

echo "wrote $REFS_DIR/cpu_strict, $REFS_DIR/gpu_alpha1_local, benchmark/validate/BASELINE_DELTAS.md"
exit $RC

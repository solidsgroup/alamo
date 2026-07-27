#!/usr/bin/env bash
# Extract NOVA plotfile metrics in the local environment (which provides yt)
# and run the strict baseline/native physics-budget comparisons.
set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "usage: $0 <copied-nova-runtime-job-directory>" >&2
    exit 2
fi

ROOT=$(git rev-parse --show-toplevel)
JOB_DIR=$(realpath "$1")
EXTRACT="$ROOT/benchmark/validate/extract_metrics.py"
COMPARE="$ROOT/benchmark/validate/compare_validation.py"
BUDGET="$ROOT/benchmark/validate/physics_budget.yaml"
CPU="$ROOT/benchmark/validate/references/cpu_strict"

for ranks in 1 2; do
    base="$JOB_DIR/physics-${ranks}rank/baseline"
    native="$JOB_DIR/physics-${ranks}rank/native"
    for arm in "$base" "$native"; do
        python3 "$EXTRACT" "$arm/canonical_2d_elastic" --budget "$BUDGET"
    done

    python3 "$COMPARE" "$base" "$native" \
        --case canonical_2d_elastic --gate \
        --out-md "$JOB_DIR/physics-${ranks}rank/compare.md" \
        --out-json "$JOB_DIR/physics-${ranks}rank/compare.json"

    set +e
    python3 "$COMPARE" "$CPU" "$base" \
        --case canonical_2d_elastic --gate \
        --out-md "$JOB_DIR/physics-${ranks}rank/baseline-vs-cpu.md" \
        --out-json "$JOB_DIR/physics-${ranks}rank/baseline-vs-cpu.json"
    baseline_cpu_rc=$?
    python3 "$COMPARE" "$CPU" "$native" \
        --case canonical_2d_elastic --gate \
        --out-md "$JOB_DIR/physics-${ranks}rank/native-vs-cpu.md" \
        --out-json "$JOB_DIR/physics-${ranks}rank/native-vs-cpu.json"
    native_cpu_rc=$?
    set -e

    printf 'baseline=%d native=%d\n' "$baseline_cpu_rc" "$native_cpu_rc" \
        | tee "$JOB_DIR/physics-${ranks}rank/cpu-golden-exit.txt"
    test "$baseline_cpu_rc" -eq "$native_cpu_rc"
done

echo "PASS: NOVA 1-rank and 2-rank strict baseline/native physics gates"

if [ "${CHECK_CROSS_RANK:-0}" = 1 ]; then
    for arm in baseline native; do
        python3 "$COMPARE" \
            "$JOB_DIR/physics-1rank/$arm" "$JOB_DIR/physics-2rank/$arm" \
            --case canonical_2d_elastic --gate \
            --out-md "$JOB_DIR/physics-1rank-vs-2rank-$arm.md" \
            --out-json "$JOB_DIR/physics-1rank-vs-2rank-$arm.json"
    done
    echo "PASS: fixed-BoxArray one-rank/two-rank physics gates"
fi

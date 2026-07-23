#!/usr/bin/env bash
set -euo pipefail
set -x

REPO=/home/jackplum/Projects/alamo
TASK="$REPO/docs/agent_plans/20260721-fapply-runtime-optimization"
ART="$TASK/artifacts/a1000-sm86-2e6a8f8f-20260721"
CAND="$ART/candidate-step3"
ORACLE="$CAND/oracle"
VALIDATION="$CAND/validation"
COMBINED="$CAND/combined_2x2"
RESUME="${RESUME:-0}"

GPU2_STRICT="$CAND/isolated/bin/alamo_gpu-2d-strict-sm86-step3-isolated"
GPU2_FAST="$CAND/isolated/bin/alamo_gpu-2d-fast-sm86-step3-isolated"
GPU3_STRICT="$CAND/isolated/bin/alamo_gpu-3d-strict-sm86-step3-isolated"
GPU3_FAST="$CAND/isolated/bin/alamo_gpu-3d-fast-sm86-step3-isolated"

BASE_STRICT_2="$ART/baseline/validation/strict_2d"
BASE_STRICT_3="$ART/baseline/validation/strict_3d"
BASE_FAST_2="$ART/baseline/validation/fast_2d"
BASE_FAST_3="$ART/baseline/validation/fast_3d"
BASE_2X2_3="$ART/step2/correctness_2x2/fast_3d"
MANIFEST_2X2="$TASK/step2_cases_2x2.manifest.yaml"

# `--require-compatible-manifest` currently compares build_flags as a literal
# string, including artifact paths. Use the frozen baseline strings as the
# compilation-identity key for canonical A/B bundles; each candidate manifest
# still records its actual binary path/hash, and the raw build logs plus
# CANDIDATE_SHA256SUMS record the actual candidate output.
FLAGS_STRICT_2="COMP=g++ DIM=2 PROFILE=0 CUDA_FP=strict BUILD_JOBS=4 SMOKE=0 ARCH=86 SOURCE_BINARY=$REPO/bin/alamo_gpu-2d-nofast-cuda86-g++ OUTPUT_BINARY=$ART/baseline/bin/alamo_gpu-2d-strict-sm86-baseline"
FLAGS_STRICT_3="COMP=g++ DIM=3 PROFILE=0 CUDA_FP=strict BUILD_JOBS=4 SMOKE=0 ARCH=86 SOURCE_BINARY=$REPO/bin/alamo_gpu-3d-nofast-cuda86-g++ OUTPUT_BINARY=$ART/baseline/bin/alamo_gpu-3d-strict-sm86-baseline"
FLAGS_FAST_2="COMP=g++ DIM=2 PROFILE=0 CUDA_FP=fast BUILD_JOBS=4 SMOKE=0 ARCH=86 SOURCE_BINARY=$REPO/bin/alamo_gpu-2d-cuda86-g++ OUTPUT_BINARY=$ART/baseline/bin/alamo_gpu-2d-fast-sm86-baseline"
FLAGS_FAST_3="COMP=g++ DIM=3 PROFILE=0 CUDA_FP=fast BUILD_JOBS=4 SMOKE=0 ARCH=86 SOURCE_BINARY=$REPO/bin/alamo_gpu-3d-cuda86-g++ OUTPUT_BINARY=$ART/baseline/bin/alamo_gpu-3d-fast-sm86-baseline"

cd "$REPO"
if [[ "$RESUME" == 0 ]]; then
  test ! -e "$ORACLE"
  test ! -e "$VALIDATION"
  test ! -e "$COMBINED"
else
  test -d "$ORACLE/logs"
  test ! -e "$VALIDATION"
  test -d "$COMBINED/softvoid/logs"
  grep -q 'ci_golden_compare.sh (cpu): PASS' "$ORACLE/logs/ci_golden_cpu.log"
fi
for binary in "$GPU2_STRICT" "$GPU2_FAST" "$GPU3_STRICT" "$GPU3_FAST"; do
  test -x "$binary"
done
for bundle in "$BASE_STRICT_2" "$BASE_STRICT_3" "$BASE_FAST_2" \
  "$BASE_FAST_3" "$BASE_2X2_3"; do
  test -d "$bundle"
done
test "$(nvidia-smi --query-gpu=uuid --format=csv,noheader)" = \
  GPU-ff00e057-b36d-c833-9da1-8f71516e7d73
test -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)"
mkdir -p "$ORACLE/logs" "$COMBINED/softvoid/logs"

if [[ "$RESUME" == 0 ]]; then
  sha256sum "$CAND"/isolated/bin/* "$MANIFEST_2X2" \
    benchmark/validate/cases.manifest.yaml benchmark/validate/physics_budget.yaml \
    tests/ElasticSoftVoid/input tests/ElasticSoftVoid/test \
    >"$CAND/CANDIDATE_SHA256SUMS"
  git diff --check -- src/Operator/Elastic.cpp src/Operator/Elastic.H \
    src/Test/Operator/Elastic.H src/test.cc
  benchmark/lint_device_patterns.sh
  GOLDEN_MODE=cpu benchmark/ci_golden_compare.sh \
    >"$ORACLE/logs/ci_golden_cpu.log" 2>&1
fi
scripts/runtests.py --dim=2 --comp=g++ \
  --sections 2d-serial 2d-parallel -- tests/ElasticSoftVoid \
  >"$ORACLE/logs/runtests_elasticsoftvoid_cpu.log" 2>&1
if [[ ! -e "$ORACLE/a100_gate_isolated" ]]; then
  BIN="$GPU3_FAST" TIERS='1 2' OUT="$ORACLE/a100_gate_isolated" \
    benchmark/local_a100_gate.sh \
    >"$ORACLE/logs/local_a100_gate_isolated_tiers12.log" 2>&1
else
  grep -q 'ERROR SUMMARY: 0 errors' \
    "$ORACLE/a100_gate_isolated/tier2_memcheck.log"
fi

SOFTVOID_COMMON=(
  max_step=2 amr.max_level=2 amr.n_cell=64 64 8
  explicitmesh.lo1=32 32 0 explicitmesh.hi1=95 95 0
  explicitmesh.lo2=72 72 0 explicitmesh.hi2=183 183 0
  pf.eta.ic.expression.constant.w=0.002
  model_void.kappa=0.2_MPa model_void.mu=0.2_MPa
)

run_softvoid () {
  local root="$1" label="$2" binary="$3"
  shift 3
  local plot="$root/$label"
  local logs="$root/logs/$label"
  test ! -e "$plot"
  test ! -e "$logs"
  mkdir -p "$logs"
  "$binary" tests/ElasticSoftVoid/input \
    "${SOFTVOID_COMMON[@]}" "$@" plot_file="$plot" \
    >"$logs/stdout" 2>"$logs/stderr"
  cp "$logs/stdout" "$logs/stderr" "$plot/"
  python3 tests/ElasticSoftVoid/test "$plot"
}

run_softvoid "$ORACLE" strict_single "$GPU2_STRICT"
run_softvoid "$ORACLE" strict_multibox "$GPU2_STRICT" amr.max_grid_size=32
run_softvoid "$ORACLE" fast_single "$GPU2_FAST"
run_softvoid "$ORACLE" fast_multibox "$GPU2_FAST" amr.max_grid_size=32

run_validation () {
  local profile="$1" case_id="$2" binary="$3" bundle="$4" flags="$5"
  python3 benchmark/validate/run_validation_local.py \
    --profiles "$profile" --case "$case_id" --binary "$binary" \
    --bundle-dir "$bundle" \
    --build-command 'bash benchmark/build_alamo_local_gpu.sh' \
    --build-flags "$flags"
}

run_validation gpu_strict canonical_2d_elastic "$GPU2_STRICT" \
  "$VALIDATION/strict_2d" "$FLAGS_STRICT_2"
run_validation gpu_strict centre_bore_3d_128_a2_converged "$GPU3_STRICT" \
  "$VALIDATION/strict_3d" "$FLAGS_STRICT_3"
run_validation gpu_fast canonical_2d_elastic "$GPU2_FAST" \
  "$VALIDATION/fast_2d" "$FLAGS_FAST_2"
run_validation gpu_fast centre_bore_3d_128_a2_converged "$GPU3_FAST" \
  "$VALIDATION/fast_3d" "$FLAGS_FAST_3"

python3 benchmark/validate/compare_validation.py "$BASE_STRICT_2" \
  "$VALIDATION/strict_2d" --case canonical_2d_elastic --gate \
  --require-compatible-manifest
python3 benchmark/validate/compare_validation.py "$BASE_STRICT_3" \
  "$VALIDATION/strict_3d" --case centre_bore_3d_128_a2_converged --gate \
  --require-compatible-manifest
python3 benchmark/validate/compare_validation.py "$BASE_FAST_2" \
  "$VALIDATION/fast_2d" --case canonical_2d_elastic --gate \
  --require-compatible-manifest
python3 benchmark/validate/compare_validation.py "$BASE_FAST_3" \
  "$VALIDATION/fast_3d" --case centre_bore_3d_128_a2_converged --gate \
  --require-compatible-manifest

run_softvoid "$COMBINED/softvoid" fast_single "$GPU2_FAST" \
  elastic.solver.pre_smooth=2 elastic.solver.post_smooth=2
run_softvoid "$COMBINED/softvoid" fast_multibox "$GPU2_FAST" \
  elastic.solver.pre_smooth=2 elastic.solver.post_smooth=2 amr.max_grid_size=32

python3 benchmark/validate/run_validation_local.py \
  --profiles gpu_fast --case centre_bore_3d_128_a2_converged \
  --manifest "$MANIFEST_2X2" --binary "$GPU3_FAST" \
  --bundle-dir "$COMBINED/fast_3d" \
  --build-command 'bash benchmark/build_alamo_local_gpu.sh' \
  --build-flags "COMP=g++ DIM=3 PROFILE=0 CUDA_FP=fast BUILD_JOBS=4 SMOKE=0 ARCH=86 OUTPUT_BINARY=$GPU3_FAST"
python3 benchmark/validate/compare_validation.py "$BASE_2X2_3" \
  "$COMBINED/fast_3d" --case centre_bore_3d_128_a2_converged --gate \
  --out-md "$COMBINED/fast_3d_vs_step2_2x2.md" \
  --out-json "$COMBINED/fast_3d_vs_step2_2x2.json"

nvidia-smi --query-gpu=name,uuid,driver_version,pstate,clocks.current.sm,clocks.current.memory,power.limit,memory.used \
  --format=csv,noheader
nvidia-smi --query-compute-apps=pid,process_name,used_gpu_memory \
  --format=csv,noheader
git diff --check -- src/Operator/Elastic.cpp src/Operator/Elastic.H \
  src/Test/Operator/Elastic.H src/test.cc
set +x
echo STEP3_CORRECTNESS_AND_2X2_STABILITY_PASS

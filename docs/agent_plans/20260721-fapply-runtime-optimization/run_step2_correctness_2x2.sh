#!/usr/bin/env bash
set -euo pipefail
set -x

REPO=/home/jackplum/Projects/alamo
TASK="$REPO/docs/agent_plans/20260721-fapply-runtime-optimization"
ART="$TASK/artifacts/a1000-sm86-2e6a8f8f-20260721"
OUT="$ART/step2/correctness_2x2"
BIN2="$ART/baseline/bin/alamo_gpu-2d-fast-sm86-baseline"
BIN3="$ART/baseline/bin/alamo_gpu-3d-fast-sm86-baseline"
BASE3="$ART/baseline/validation/fast_3d"
MANIFEST="$TASK/step2_cases_2x2.manifest.yaml"

cd "$REPO"
test ! -e "$OUT"
git diff --exit-code -- src/Operator/Elastic.cpp src/Operator/Elastic.H src/Set/Matrix4_Major.H
test "$(nvidia-smi --query-gpu=uuid --format=csv,noheader)" = \
  GPU-ff00e057-b36d-c833-9da1-8f71516e7d73
test -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)"
sha256sum "$BIN2" "$BIN3" tests/ElasticSoftVoid/input \
  input_3d_centre_bore_128_a2 "$MANIFEST" tests/ElasticSoftVoid/test \
  benchmark/validate/physics_budget.yaml benchmark/validate/extract_metrics.py \
  benchmark/validate/compare_validation.py
mkdir -p "$OUT/softvoid/logs"
nvidia-smi --query-gpu=name,uuid,driver_version,pstate,clocks.current.sm,clocks.current.memory,power.limit,memory.used --format=csv,noheader

run_softvoid () {
  local label="$1"
  shift
  local plot="$OUT/softvoid/$label"
  local logs="$OUT/softvoid/logs/$label"
  test ! -e "$plot"
  test ! -e "$logs"
  mkdir -p "$logs"
  "$BIN2" tests/ElasticSoftVoid/input \
    max_step=2 amr.max_level=2 amr.n_cell=64 64 8 \
    explicitmesh.lo1=32 32 0 explicitmesh.hi1=95 95 0 \
    explicitmesh.lo2=72 72 0 explicitmesh.hi2=183 183 0 \
    pf.eta.ic.expression.constant.w=0.002 \
    model_void.kappa=0.2_MPa model_void.mu=0.2_MPa \
    elastic.solver.pre_smooth=2 elastic.solver.post_smooth=2 \
    "$@" plot_file="$plot" >"$logs/stdout" 2>"$logs/stderr"
  cp "$logs/stdout" "$logs/stderr" "$plot/"
  python3 tests/ElasticSoftVoid/test "$plot"
}

run_softvoid fast_single
run_softvoid fast_multibox amr.max_grid_size=32

python3 benchmark/validate/run_validation_local.py \
  --profiles gpu_fast --case centre_bore_3d_128_a2_converged \
  --manifest "$MANIFEST" --binary "$BIN3" \
  --bundle-dir "$OUT/fast_3d" \
  --build-command 'bash benchmark/build_alamo_local_gpu.sh' \
  --build-flags 'COMP=g++ DIM=3 PROFILE=0 CUDA_FP=fast BUILD_JOBS=4 SMOKE=0 ARCH=86 SOURCE_BINARY=/home/jackplum/Projects/alamo/bin/alamo_gpu-3d-cuda86-g++ OUTPUT_BINARY=/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline/bin/alamo_gpu-3d-fast-sm86-baseline'

# Manifest compatibility is intentionally not required here: the smoothing
# override and task-local manifest are the independent variable under test.
python3 benchmark/validate/compare_validation.py "$BASE3" "$OUT/fast_3d" \
  --case centre_bore_3d_128_a2_converged --gate \
  --out-md "$OUT/fast_3d_vs_baseline.md" \
  --out-json "$OUT/fast_3d_vs_baseline.json"

nvidia-smi --query-gpu=name,uuid,driver_version,pstate,clocks.current.sm,clocks.current.memory,power.limit,memory.used --format=csv,noheader
nvidia-smi --query-compute-apps=pid,process_name,used_gpu_memory --format=csv,noheader
git diff --exit-code -- src/Operator/Elastic.cpp src/Operator/Elastic.H src/Set/Matrix4_Major.H
set +x
echo STEP2_CORRECTNESS_2X2_PASS

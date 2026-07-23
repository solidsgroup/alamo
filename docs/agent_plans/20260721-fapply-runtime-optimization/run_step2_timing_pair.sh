#!/usr/bin/env bash
set -euo pipefail
set -x

if [[ $# -ne 1 || ! "$1" =~ ^(2d|3d)$ ]]; then
  echo "usage: $0 2d|3d" >&2
  exit 2
fi

REGIME=$1
REPO=/home/jackplum/Projects/alamo
TASK="$REPO/docs/agent_plans/20260721-fapply-runtime-optimization"
ART="$TASK/artifacts/a1000-sm86-2e6a8f8f-20260721"
OUT="$ART/step2/smoothing_sweep"

cd "$REPO"
git diff --exit-code -- src/Operator/Elastic.cpp src/Operator/Elastic.H src/Set/Matrix4_Major.H
test "$(nvidia-smi --query-gpu=uuid --format=csv,noheader)" = \
  GPU-ff00e057-b36d-c833-9da1-8f71516e7d73
test -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)"

if [[ "$REGIME" == 2d ]]; then
  CASE=2d_conservative
  INPUT=tests/ElasticSoftVoid/input
  FAST="$ART/baseline/bin/alamo_gpu-2d-fast-sm86-baseline"
  PROFILE="$ART/baseline/bin/alamo_gpu-2d-profile-fast-sm86-baseline"
  COMMON=(
    max_step=2 amr.max_level=2 amr.n_cell=64 64 8
    explicitmesh.lo1=32 32 0 explicitmesh.hi1=95 95 0
    explicitmesh.lo2=72 72 0 explicitmesh.hi2=183 183 0
    pf.eta.ic.expression.constant.w=0.002
    model_void.kappa=0.2_MPa model_void.mu=0.2_MPa
    elastic.print_model=0 elastic.print_residual=0
    elastic.solver.verbose=0 elastic.solver.nr_diagnostics=0
    amr.plot_int=-1 amr.thermo.plot_int=-1
  )
else
  CASE=3d_psi
  INPUT=input_3d_centre_bore_128_a2
  FAST="$ART/baseline/bin/alamo_gpu-3d-fast-sm86-baseline"
  PROFILE="$ART/baseline/bin/alamo_gpu-3d-profile-fast-sm86-baseline"
  COMMON=(
    max_step=2 stop_time=1e99_s elastic.interval=1
    amr.max_grid_size=64 amrex.the_arena_is_managed=1
    elastic.solver.nriters=2 elastic.print_model=0
    elastic.solver.verbose=0 elastic.solver.nr_diagnostics=0
    amr.plot_int=-1 amr.thermo.plot_int=-1
  )
fi

test ! -e "$OUT/$CASE"
sha256sum "$FAST" "$PROFILE" "$INPUT"

run_one () {
  local setting="$1" smooth="$2" mode="$3" label="$4" binary="$5"
  local dir="$OUT/$CASE/$setting/$mode/$label"
  test ! -e "$dir"
  mkdir -p "$dir"
  nvidia-smi --query-gpu=name,uuid,driver_version,pstate,clocks.current.sm,clocks.current.memory,power.limit,memory.used --format=csv,noheader >"$dir/gpu_before.csv"
  nvidia-smi --query-compute-apps=pid,process_name,used_gpu_memory --format=csv,noheader >"$dir/compute_before.csv"
  test ! -s "$dir/compute_before.csv"
  /usr/bin/time -f %e -o "$dir/wall_seconds.txt" \
    "$binary" "$INPUT" "${COMMON[@]}" \
    elastic.solver.pre_smooth="$smooth" elastic.solver.post_smooth="$smooth" \
    plot_file="$dir/plot" >"$dir/stdout" 2>"$dir/stderr"
  nvidia-smi --query-gpu=name,uuid,driver_version,pstate,clocks.current.sm,clocks.current.memory,power.limit,memory.used --format=csv,noheader >"$dir/gpu_after.csv"
  nvidia-smi --query-compute-apps=pid,process_name,used_gpu_memory --format=csv,noheader >"$dir/compute_after.csv"
  test ! -s "$dir/compute_after.csv"
  printf '%s %s %s ' "$setting" "$mode" "$label"
  tr -d '\n' <"$dir/wall_seconds.txt"
  printf '\n'
}

# Candidate first and fresh 4/4 drift arm second in every pair. Thus the final
# measured run is the plan-required unmodified 4/4 post-sweep drift check.
for mode in unprofiled profiled; do
  if [[ "$mode" == unprofiled ]]; then binary="$FAST"; else binary="$PROFILE"; fi
  for label in warmup rep1 rep2 rep3 rep4 rep5; do
    run_one 2x2 2 "$mode" "$label" "$binary"
    run_one drift_4x4 4 "$mode" "$label" "$binary"
  done
done

git diff --exit-code -- src/Operator/Elastic.cpp src/Operator/Elastic.H src/Set/Matrix4_Major.H
set +x
echo "STEP2_TIMING_${REGIME^^}_PASS"

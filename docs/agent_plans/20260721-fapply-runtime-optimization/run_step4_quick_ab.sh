#!/usr/bin/env bash
set -euo pipefail

REPO=/home/jackplum/Projects/alamo
TASK="$REPO/docs/agent_plans/20260721-fapply-runtime-optimization"
ART="$TASK/artifacts/a1000-sm86-2e6a8f8f-20260721"
OUT="$ART/candidate-step4/timing/3d_psi_2x2_quick"
BASE="$ART/baseline/bin/alamo_gpu-3d-fast-sm86-baseline"
CAND="$ART/candidate-step4/isolated/bin/alamo_gpu-3d-fast-sm86-step4-isolated"
INPUT=input_3d_centre_bore_128_a2

cd "$REPO"
test ! -e "$OUT"
test "$(nvidia-smi --query-gpu=uuid --format=csv,noheader)" = \
  GPU-ff00e057-b36d-c833-9da1-8f71516e7d73
test -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)"
mkdir -p "$OUT"
sha256sum "$BASE" "$CAND" "$INPUT"

COMMON=(
  max_step=2 stop_time=1e99_s elastic.interval=1
  amr.max_grid_size=64 amrex.the_arena_is_managed=1
  elastic.solver.nriters=2 elastic.print_model=0
  elastic.solver.verbose=0 elastic.solver.nr_diagnostics=0
  elastic.solver.pre_smooth=2 elastic.solver.post_smooth=2
  amr.plot_int=-1 amr.thermo.plot_int=-1
)

run_one () {
  local arm="$1" label="$2" binary="$3"
  local dir="$OUT/$arm/$label"
  mkdir -p "$dir"
  test ! -e "$dir/wall_seconds.txt"
  nvidia-smi --query-compute-apps=pid,process_name,used_gpu_memory \
    --format=csv,noheader >"$dir/compute_before.csv"
  test ! -s "$dir/compute_before.csv"
  /usr/bin/time -f %e -o "$dir/wall_seconds.txt" \
    "$binary" "$INPUT" "${COMMON[@]}" plot_file="$dir/plot" \
    >"$dir/stdout" 2>"$dir/stderr"
  nvidia-smi --query-compute-apps=pid,process_name,used_gpu_memory \
    --format=csv,noheader >"$dir/compute_after.csv"
  test ! -s "$dir/compute_after.csv"
  printf '%s %s %s\n' "$arm" "$label" "$(<"$dir/wall_seconds.txt")"
}

run_one baseline warmup "$BASE"
run_one candidate warmup "$CAND"
for rep in 1 2 3; do
  if (( rep % 2 == 1 )); then
    run_one baseline "rep$rep" "$BASE"
    run_one candidate "rep$rep" "$CAND"
  else
    run_one candidate "rep$rep" "$CAND"
    run_one baseline "rep$rep" "$BASE"
  fi
done

echo STEP4_QUICK_AB_PASS

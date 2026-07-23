#!/usr/bin/env bash
#SBATCH --account=brunnels
#SBATCH --partition=nova
#SBATCH --job-name=fapply_a100_ab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:a100:1
#SBATCH --time=02:00:00
#SBATCH --output=nova-a100-ab.%j.out
#SBATCH --error=nova-a100-ab.%j.err

set -euo pipefail

REPO=/work/brunnels/jackplum/alamo-fapply-20260722
BIN="$REPO/bin/alamo_gpu-3d-profile-cuda80-g++"
INPUT="$REPO/input_3d_centre_bore_128_a2"
MANIFEST_2X2="$REPO/nova-step2-cases-2x2.manifest.yaml"
VALIDATE_EXACT="$REPO/nova-run-validation-exact.py"
PYTHON="$REPO/.nova-validation-venv/bin/python"
OUT="$REPO/nova-results-$SLURM_JOB_ID"
EXPECTED_HEAD=2e6a8f8f5430a58dd6a85b07b08d472e2bd1d5cd

cd "$REPO"
test "$(git rev-parse HEAD)" = "$EXPECTED_HEAD"
test -x "$BIN"
test -f "$INPUT"
test -f "$MANIFEST_2X2"
test -f "$VALIDATE_EXACT"
test -x "$PYTHON"
test ! -e "$OUT"
mkdir -p "$OUT/timing"

module purge 2>/dev/null || true
module load cuda 2>/dev/null || module load cuda/12
module load gcc 2>/dev/null || module load gcc/12
module load openmpi 2>/dev/null || module load openmpi4
module -t list
nvidia-smi -L
nvidia-smi --query-gpu=name,uuid,driver_version,compute_cap,memory.total,power.limit \
  --format=csv,noheader
"$PYTHON" -c 'import numpy, yaml, yt; print("numpy", numpy.__version__, "yaml", yaml.__version__, "yt", yt.__version__)' \
  >"$OUT/validation_python_versions.txt"
sha256sum "$BIN" "$INPUT" "$MANIFEST_2X2" "$VALIDATE_EXACT" \
  benchmark/validate/cases.manifest.yaml benchmark/validate/physics_budget.yaml \
  benchmark/validate/extract_metrics.py benchmark/validate/compare_validation.py \
  >"$OUT/input_sha256s.txt"

COMMON=(
  max_step=2 stop_time=1e99_s elastic.interval=1
  amr.max_grid_size=64 amrex.the_arena_is_managed=1
  elastic.solver.nriters=2 elastic.print_model=0
  elastic.solver.verbose=0 elastic.solver.nr_diagnostics=0
  amr.plot_int=-1 amr.thermo.plot_int=-1
)

run_one () {
  local arm="$1" smooth="$2" label="$3"
  local dir="$OUT/timing/$arm/$label"
  mkdir -p "$dir"
  nvidia-smi --query-gpu=name,uuid,pstate,clocks.current.sm,clocks.current.memory,power.limit,memory.used \
    --format=csv,noheader >"$dir/gpu_before.csv"
  nvidia-smi --query-compute-apps=pid,process_name,used_gpu_memory \
    --format=csv,noheader >"$dir/compute_before.csv"
  /usr/bin/time -f %e -o "$dir/wall_seconds.txt" \
    "$BIN" "$INPUT" "${COMMON[@]}" \
    elastic.solver.pre_smooth="$smooth" elastic.solver.post_smooth="$smooth" \
    plot_file="$dir/plot" >"$dir/stdout" 2>"$dir/stderr"
  nvidia-smi --query-gpu=name,uuid,pstate,clocks.current.sm,clocks.current.memory,power.limit,memory.used \
    --format=csv,noheader >"$dir/gpu_after.csv"
  nvidia-smi --query-compute-apps=pid,process_name,used_gpu_memory \
    --format=csv,noheader >"$dir/compute_after.csv"
  printf '%s %s %s\n' "$arm" "$label" "$(<"$dir/wall_seconds.txt")"
}

run_one baseline_4x4 4 warmup
run_one candidate_2x2 2 warmup
for rep in 1 2 3 4 5; do
  if (( rep % 2 == 1 )); then
    run_one baseline_4x4 4 "rep$rep"
    run_one candidate_2x2 2 "rep$rep"
  else
    run_one candidate_2x2 2 "rep$rep"
    run_one baseline_4x4 4 "rep$rep"
  fi
done

"$PYTHON" "$VALIDATE_EXACT" \
  --profile gpu_fast --case centre_bore_3d_128_a2_converged \
  --binary "$BIN" --bundle-dir "$OUT/validation_4x4" \
  --build-command 'sbatch nova-build.sbatch' \
  --build-flags "HEAD=$EXPECTED_HEAD ARCH=80 PROFILE=1 CUDA_FP=fast smoothing=4/4"
"$PYTHON" "$VALIDATE_EXACT" \
  --profile gpu_fast --case centre_bore_3d_128_a2_converged \
  --manifest "$MANIFEST_2X2" --binary "$BIN" \
  --bundle-dir "$OUT/validation_2x2" \
  --build-command 'sbatch nova-build.sbatch' \
  --build-flags "HEAD=$EXPECTED_HEAD ARCH=80 PROFILE=1 CUDA_FP=fast smoothing=2/2"
"$PYTHON" benchmark/validate/compare_validation.py \
  "$OUT/validation_4x4" "$OUT/validation_2x2" \
  --case centre_bore_3d_128_a2_converged --gate \
  --out-md "$OUT/validation_2x2_vs_4x4.md" \
  --out-json "$OUT/validation_2x2_vs_4x4.json"

bash benchmark/fapply_register_ab.sh "$BIN" >"$OUT/fapply_resources.txt"
find "$OUT" -type f ! -path '*/plot/*' ! -name ARTIFACT_SHA256SUMS -print0 \
  | sort -z | xargs -0 sha256sum \
  >"$OUT/ARTIFACT_SHA256SUMS"
echo NOVA_A100_AB_PASS

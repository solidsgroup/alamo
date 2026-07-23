#!/usr/bin/env bash
#SBATCH --account=brunnels
#SBATCH --partition=interactive
#SBATCH --job-name=fapply_a100_validate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:a100:1
#SBATCH --time=01:00:00
#SBATCH --output=nova-a100-validate.%j.out
#SBATCH --error=nova-a100-validate.%j.err

set -euo pipefail

REPO=/work/brunnels/jackplum/alamo-fapply-20260722
BIN="$REPO/bin/alamo_gpu-3d-profile-cuda80-g++"
MANIFEST_2X2="$REPO/nova-step2-cases-2x2.manifest.yaml"
VALIDATE_EXACT="$REPO/nova-run-validation-exact.py"
PYTHON="$REPO/.nova-validation-venv/bin/python"
OUT="$REPO/nova-results-11751561"
EXPECTED_HEAD=2e6a8f8f5430a58dd6a85b07b08d472e2bd1d5cd

cd "$REPO"
test "$(git rev-parse HEAD)" = "$EXPECTED_HEAD"
test -x "$BIN"
test -f "$MANIFEST_2X2"
test -f "$VALIDATE_EXACT"
test -x "$PYTHON"
test -d "$OUT/timing"
test ! -e "$OUT/validation_4x4"
test ! -e "$OUT/validation_2x2"
for arm in baseline_4x4 candidate_2x2; do
  for rep in 1 2 3 4 5; do
    test -s "$OUT/timing/$arm/rep$rep/wall_seconds.txt"
  done
done

module purge 2>/dev/null || true
module load cuda 2>/dev/null || module load cuda/12
module load gcc 2>/dev/null || module load gcc/12
module load openmpi 2>/dev/null || module load openmpi4
nvidia-smi -L
"$PYTHON" -c 'import numpy, yaml, yt; print("numpy", numpy.__version__, "yaml", yaml.__version__, "yt", yt.__version__)'

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
  | sort -z | xargs -0 sha256sum >"$OUT/ARTIFACT_SHA256SUMS"
echo NOVA_A100_VALIDATION_TAIL_PASS

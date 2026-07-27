#!/usr/bin/env bash
# Cross-rank physics gate with an explicit grid decomposition.  Keeping the
# BoxArray fixed avoids rank-dependent nodal-patch overlap in yt's L2 norms.
#SBATCH --account=brunnels
#SBATCH --partition=nova
#SBATCH --gres=gpu:a100:2
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=00:15:00
#SBATCH --job-name=mlmg_nosync_mpi
#SBATCH --output=mlmg_nosync_mpi.%j.out
#SBATCH --error=mlmg_nosync_mpi.%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"

TASK="$PWD/docs/agent_plans/20260727-mlmg-nosync-native-validation"
OUT="$TASK/artifacts/nova-fixed-grid-mpi/$SLURM_JOB_ID"
mkdir -p "$OUT"

module purge 2>/dev/null || true
module load cuda 2>/dev/null || module load cuda/12 2>/dev/null || true
module load gcc 2>/dev/null || module load gcc/12 2>/dev/null || true
module load openmpi 2>/dev/null || module load openmpi4 2>/dev/null || true
module list 2>&1 | tee "$OUT/modules.txt"
nvidia-smi --query-gpu=name,compute_cap,memory.total,driver_version,uuid \
    --format=csv | tee "$OUT/gpus.csv"
git rev-parse HEAD | tee "$OUT/git-head.txt"

BIN="$PWD/bin/alamo_gpu-2d-profile-nofast-cuda80-g++"
test -x "$BIN"
sha256sum "$BIN" | tee "$OUT/binary-sha256.txt"

COMMON=(
    max_step=2
    amr.max_grid_size=32
    elastic.solver.verbose=5
    elastic.print_model=0
    amr.plot_int=1
    amr.thermo.plot_int=1
    amr.thermo.int=1
    stop_time=1e99_s
    elastic.interval=1
)

run_case () {
    local ranks="$1" arm="$2" enabled=0
    local case_dir="$OUT/physics-${ranks}rank/$arm/canonical_2d_elastic"
    [ "$arm" = "native" ] && enabled=1
    mkdir -p "$case_dir"
    srun --mpi=pmix --exclusive -n "$ranks" --gpus-per-task=1 \
        "$BIN" input "${COMMON[@]}" \
        elastic.solver.no_gpu_sync="$enabled" \
        plot_file="$case_dir/_plot" >"$case_dir/run.log" 2>&1
    test "$(grep -c '^STEP [0-9][0-9]* ends' "$case_dir/run.log")" -eq 2
    grep -q 'AMReX (.*) finalized' "$case_dir/run.log"
    cp "$case_dir/_plot/thermo.dat" "$case_dir/thermo.dat"
    find "$case_dir/_plot" -mindepth 1 -maxdepth 1 -type d \
        \( -name '[0-9]*node' -o -name '[0-9]*cell' \) \
        -exec cp -a '{}' "$case_dir/" ';'
}

run_case 1 baseline
run_case 1 native
run_case 2 baseline
run_case 2 native

for ranks in 1 2; do
    grep -E 'MLMG: (Initial rhs|Initial residual)' \
        "$OUT/physics-${ranks}rank/baseline/canonical_2d_elastic/run.log" \
        >"$OUT/physics-${ranks}rank/baseline.initials"
    grep -E 'MLMG: (Initial rhs|Initial residual)' \
        "$OUT/physics-${ranks}rank/native/canonical_2d_elastic/run.log" \
        >"$OUT/physics-${ranks}rank/native.initials"
    awk '/MLMG: Final Iter\./ { print $4 }' \
        "$OUT/physics-${ranks}rank/baseline/canonical_2d_elastic/run.log" \
        >"$OUT/physics-${ranks}rank/baseline.iterations"
    awk '/MLMG: Final Iter\./ { print $4 }' \
        "$OUT/physics-${ranks}rank/native/canonical_2d_elastic/run.log" \
        >"$OUT/physics-${ranks}rank/native.iterations"
    cmp "$OUT/physics-${ranks}rank/baseline.initials" \
        "$OUT/physics-${ranks}rank/native.initials"
    cmp "$OUT/physics-${ranks}rank/baseline.iterations" \
        "$OUT/physics-${ranks}rank/native.iterations"
done

echo "PASS: fixed-BoxArray one-rank and two-rank A/B runs" \
    | tee "$OUT/RESULT.txt"

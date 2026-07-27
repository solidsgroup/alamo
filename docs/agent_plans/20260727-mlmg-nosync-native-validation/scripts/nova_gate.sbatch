#!/usr/bin/env bash
# Builds the strict A100 binaries, writes raw one-/two-rank plotfiles, runs
# Compute Sanitizer, and records ten-step timings. After copying OUT locally,
# run postprocess_nova_runtime.sh there to extract yt metrics and gate physics.
#SBATCH --account=brunnels
#SBATCH --partition=nova
#SBATCH --gres=gpu:a100:2
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --job-name=mlmg_nosync_gate
#SBATCH --output=mlmg_nosync_gate.%j.out
#SBATCH --error=mlmg_nosync_gate.%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"

TASK="$PWD/docs/agent_plans/20260727-mlmg-nosync-native-validation"
OUT="$TASK/artifacts/nova/$SLURM_JOB_ID"
mkdir -p "$OUT"

module purge 2>/dev/null || true
module load cuda 2>/dev/null || module load cuda/12 2>/dev/null || true
module load gcc 2>/dev/null || module load gcc/12 2>/dev/null || true
module load openmpi 2>/dev/null || module load openmpi4 2>/dev/null || true
module list 2>&1 | tee "$OUT/modules.txt"
nvidia-smi --query-gpu=name,compute_cap,memory.total,driver_version,uuid \
    --format=csv | tee "$OUT/gpus.csv"
git rev-parse HEAD | tee "$OUT/git-head.txt"
git status --short | tee "$OUT/git-status.txt"
git -C ext/AMReX-Codes/amrex describe --always --dirty --tags \
    | tee "$OUT/amrex-version.txt"

./configure --comp=g++ --dim=2 --cuda=80 --cuda-fp=strict --profile --get-eigen \
    >"$OUT/configure-2d.log" 2>&1
make -j16 bin/alamo_gpu-2d-profile-nofast-cuda80-g++ \
    >"$OUT/build-2d.log" 2>&1
./configure --comp=g++ --dim=3 --cuda=80 --cuda-fp=strict --profile --get-eigen \
    >"$OUT/configure-3d.log" 2>&1
make -j16 bin/alamo_gpu-3d-profile-nofast-cuda80-g++ \
    >"$OUT/build-3d.log" 2>&1

BIN2="$PWD/bin/alamo_gpu-2d-profile-nofast-cuda80-g++"
BIN3="$PWD/bin/alamo_gpu-3d-profile-nofast-cuda80-g++"
sha256sum "$BIN2" "$BIN3" | tee "$OUT/binary-sha256.txt"

PHYSICS_ARGS=(
    max_step=2
    elastic.solver.verbose=5
    elastic.print_model=0
    amr.plot_int=1
    amr.thermo.plot_int=1
    amr.thermo.int=1
    stop_time=1e99_s
    elastic.interval=1
)

run_physics () {
    local ranks="$1" arm="$2" enabled=0
    local case_dir="$OUT/physics-${ranks}rank/$arm/canonical_2d_elastic"
    [ "$arm" = "native" ] && enabled=1
    mkdir -p "$case_dir"
    srun --mpi=pmix --exclusive -n "$ranks" --gpus-per-task=1 \
        "$BIN2" input "${PHYSICS_ARGS[@]}" \
        elastic.solver.no_gpu_sync="$enabled" \
        plot_file="$case_dir/_plot" >"$case_dir/run.log" 2>&1
    test "$(grep -c '^STEP [0-9][0-9]* ends' "$case_dir/run.log")" -eq 2
    grep -q 'AMReX (.*) finalized' "$case_dir/run.log"
    cp "$case_dir/_plot/thermo.dat" "$case_dir/thermo.dat"
    find "$case_dir/_plot" -mindepth 1 -maxdepth 1 -type d \
        \( -name '[0-9]*node' -o -name '[0-9]*cell' \) \
        -exec cp -a '{}' "$case_dir/" ';'
}

run_physics 1 baseline
run_physics 1 native
run_physics 2 baseline
run_physics 2 native

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

SAN=$(command -v compute-sanitizer)
srun --mpi=pmix --exclusive -n 1 --gpus-per-task=1 \
    "$SAN" --tool memcheck --error-exitcode 99 --launch-timeout 120 \
    "$BIN2" input max_step=2 elastic.interval=1 \
    elastic.solver.no_gpu_sync=1 elastic.solver.verbose=2 \
    elastic.print_model=0 amr.plot_int=-1 amr.thermo.plot_int=-1 \
    plot_file="$OUT/sanitizer-output" >"$OUT/sanitizer.log" 2>&1
grep 'ERROR SUMMARY' "$OUT/sanitizer.log" | tee "$OUT/sanitizer-summary.txt"
grep -q 'ERROR SUMMARY: 0 errors' "$OUT/sanitizer-summary.txt"

COMMON=(
    explicitmesh.lo1=32 32 0
    pf.eta.ic.expression.constant.w=0.002
    model_void.kappa=0.2_MPa model_void.mu=0.2_MPa
    elastic.print_model=0 elastic.print_residual=0
    elastic.solver.verbose=2 elastic.solver.nr_diagnostics=0
    amr.plot_int=-1 amr.thermo.plot_int=-1
)

run_timed () {
    local case_name="$1" arm="$2" rep="$3" steps="$4"
    local enabled=0 bin="$BIN2"
    local -a case_args
    [ "$arm" = "native" ] && enabled=1
    if [ "$case_name" = "2d" ]; then
        case_args=(
            amr.max_level=2 amr.n_cell=64 64 8
            explicitmesh.hi1=95 95 0
            explicitmesh.lo2=72 72 0 explicitmesh.hi2=183 183 0
        )
    else
        bin="$BIN3"
        case_args=(
            amr.n_cell=64 64 8 explicitmesh.hi1=95 95 15
            elastic.solver.invariant_periodic=0 0 1
        )
    fi
    local dir="$OUT/timing/$case_name"
    mkdir -p "$dir"
    local sample="$arm"
    [ "$rep" = "startup" ] && sample=startup
    local stem="$dir/$sample.rep$rep"
    local t0 t1
    t0=$(date +%s.%N)
    srun --mpi=pmix --exclusive -n 1 --gpus-per-task=1 \
        "$bin" tests/ElasticSoftVoid/input "${COMMON[@]}" "${case_args[@]}" \
        max_step="$steps" elastic.solver.no_gpu_sync="$enabled" \
        plot_file="$stem.output" >"$stem.log" 2>&1
    t1=$(date +%s.%N)
    awk -v a="$t0" -v b="$t1" 'BEGIN { printf "%.9f\n", b-a }' >>"$dir/$sample.wall"
    awk '/MLMG: Timers: Solve =/ { sum += $5; n += 1 }
         END { printf "%.9f\t%d\n", sum+0, n }' "$stem.log" >>"$dir/$sample.mlmg"
    test "$(grep -c '^STEP [0-9][0-9]* ends' "$stem.log")" -eq "$steps"
}

for case_name in 2d 3d; do
    run_timed "$case_name" baseline startup 1
    for rep in 1 2 3; do
        if (( rep % 2 )); then
            run_timed "$case_name" baseline "$rep" 10
            run_timed "$case_name" native "$rep" 10
        else
            run_timed "$case_name" native "$rep" 10
            run_timed "$case_name" baseline "$rep" 10
        fi
    done
    grep -E 'MLMG: (Initial rhs|Initial residual|Final Iter\.)' \
        "$OUT/timing/$case_name/baseline.rep1.log" \
        >"$OUT/timing/$case_name/baseline.trace"
    grep -E 'MLMG: (Initial rhs|Initial residual|Final Iter\.)' \
        "$OUT/timing/$case_name/native.rep1.log" \
        >"$OUT/timing/$case_name/native.trace"
    cmp "$OUT/timing/$case_name/baseline.trace" \
        "$OUT/timing/$case_name/native.trace"
done

echo "PASS: NOVA direct physics (1/2 rank), sanitizer, and 10-step timing" \
    | tee "$OUT/RESULT.txt"

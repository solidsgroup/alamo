#!/usr/bin/env bash
# Startup-aware native MLMG no-sync A/B benchmark.
#
#   ALAMO_CASE=2d|3d|manybox REPS=5 STEPS=10 bash interleaved_native_ab.sh
#
# Each repetition reverses arm order to spread thermal drift. A one-step run,
# which never reaches this deck's first elastic solve, estimates fixed startup.
# The summary reports total wall, wall/completed step, startup-adjusted wall per
# later step, and MLMG's own startup-excluding solver timer.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
TASK="$REPO/docs/agent_plans/20260727-mlmg-nosync-native-validation"
CASE="${ALAMO_CASE:-2d}"
REPS="${REPS:-5}"
STEPS="${STEPS:-10}"
TIMEOUT="${TIMEOUT:-3600}"
COOLDOWN="${COOLDOWN:-1}"
ARENA_INIT_SIZE="${ARENA_INIT_SIZE:-1073741824}"
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
OUT="${OUT:-$TASK/artifacts/interleaved/$CASE/$RUN_ID}"

if (( STEPS < 2 )); then
    echo "STEPS must be at least 2 so startup can be separated" >&2
    exit 2
fi
if (( REPS < 1 )); then
    echo "REPS must be positive" >&2
    exit 2
fi

cd "$REPO"
if [ -n "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ]; then
    echo "refusing to run: another CUDA process is resident" >&2
    exit 1
fi
# shellcheck source=/dev/null
source "$REPO/benchmark/local_cuda_env.sh"

INPUT=tests/ElasticSoftVoid/input
COMMON=(
    amrex.the_arena_init_size="$ARENA_INIT_SIZE"
    explicitmesh.lo1=32 32 0
    pf.eta.ic.expression.constant.w=0.002
    model_void.kappa=0.2_MPa model_void.mu=0.2_MPa
    elastic.print_model=0 elastic.print_residual=0
    elastic.solver.verbose=2 elastic.solver.nr_diagnostics=0
    amr.plot_int=-1 amr.thermo.plot_int=-1
)
case "$CASE" in
    2d)
        BIN=bin/alamo_gpu-2d-cuda86-g++
        CASE_ARGS=(
            amr.max_level=2 amr.n_cell=64 64 8
            explicitmesh.hi1=95 95 0
            explicitmesh.lo2=72 72 0 explicitmesh.hi2=183 183 0
        )
        ;;
    3d)
        BIN=bin/alamo_gpu-3d-cuda86-g++
        CASE_ARGS=(
            amr.n_cell=64 64 8 explicitmesh.hi1=95 95 15
            elastic.solver.invariant_periodic=0 0 1
        )
        ;;
    manybox)
        BIN=bin/alamo_gpu-2d-cuda86-g++
        CASE_ARGS=(
            amr.max_level=2 amr.n_cell=64 64 8 amr.max_grid_size=16
            explicitmesh.hi1=95 95 0
            explicitmesh.lo2=72 72 0 explicitmesh.hi2=183 183 0
        )
        ;;
    *)
        echo "unknown ALAMO_CASE=$CASE" >&2
        exit 2
        ;;
esac
test -x "$BIN"
mkdir -p "$OUT"

run_one () {
    local arm="$1" rep="$2" steps="$3"
    local enabled=0
    [ "$arm" = "native" ] && enabled=1
    local stem="$OUT/$arm.rep$rep"
    local plot="$OUT/output.$arm.rep$rep"
    local attempt=1 t0 t1 wall
    while true; do
        t0=$(date +%s.%N)
        if timeout "$TIMEOUT" "$BIN" "$INPUT" \
            "${COMMON[@]}" "${CASE_ARGS[@]}" \
            max_step="$steps" plot_file="$plot" \
            elastic.solver.no_gpu_sync="$enabled" >"$stem.log" 2>&1; then
            break
        fi
        if (( attempt < 3 )) && grep -q 'Arena out of memory' "$stem.log"; then
            echo "retrying $arm rep $rep after transient arena allocation failure" >&2
            attempt=$((attempt + 1))
            sleep 3
        else
            echo "$arm rep $rep failed" >&2
            tail -30 "$stem.log" >&2
            exit 1
        fi
    done
    t1=$(date +%s.%N)
    wall=$(awk -v a="$t0" -v b="$t1" 'BEGIN { printf "%.9f", b-a }')
    local completed solves solve_time
    completed=$(grep -c '^STEP [0-9][0-9]* ends' "$stem.log" || true)
    solves=$(grep -c 'MLMG: Timers: Solve =' "$stem.log" || true)
    solve_time=$(awk '/MLMG: Timers: Solve =/ { sum += $5 } END { printf "%.9f", sum+0 }' "$stem.log")
    if (( completed != steps )); then
        echo "$arm rep $rep completed $completed/$steps steps" >&2
        tail -30 "$stem.log" >&2
        exit 1
    fi
    printf '%s\t%s\t%s\t%s\n' "$wall" "$completed" "$solves" "$solve_time" >>"$OUT/$arm.tsv"
    sleep "$COOLDOWN"
}

# max_step=1 ends before the first elastic solve in ElasticSoftVoid.
run_one baseline startup 1
startup_wall=$(awk 'NR == 1 { print $1 }' "$OUT/baseline.tsv")
mv "$OUT/baseline.repstartup.log" "$OUT/startup.log"
mv "$OUT/output.baseline.repstartup" "$OUT/output.startup"
rm "$OUT/baseline.tsv"

for rep in $(seq 1 "$REPS"); do
    if (( rep % 2 == 1 )); then
        arms=(baseline native)
    else
        arms=(native baseline)
    fi
    for arm in "${arms[@]}"; do
        run_one "$arm" "$rep" "$STEPS"
    done
done

for arm in baseline native; do
    grep -E 'MLMG: (Initial rhs|Initial residual|Final Iter\.)' \
        "$OUT/$arm.rep1.log" >"$OUT/$arm.rep1.trace"
done
if cmp -s "$OUT/baseline.rep1.trace" "$OUT/native.rep1.trace"; then
    trace_status=identical
else
    trace_status=different
fi

summarize () {
    local arm="$1"
    sort -g -k1,1 "$OUT/$arm.tsv" >"$OUT/$arm.sorted.tsv"
    awk -v arm="$arm" -v startup="$startup_wall" '
        {
            wall[NR]=$1
            wall_step[NR]=$1/$2
            adjusted[NR]=($1-startup)/($2-1)
            solve[NR]=$4
            solve_call[NR]=($3 > 0 ? $4/$3 : 0)
        }
        END {
            mid=int((NR+1)/2)
            # Rows are sorted by total wall; report metrics from the median run.
            printf "%-8s total=%.4fs wall/step=%.4fs startup-adjusted=%.4fs", \
                arm, wall[mid], wall_step[mid], adjusted[mid]
            printf " mlmg_total=%.4fs mlmg/call=%.6fs n=%d\n", \
                solve[mid], solve_call[mid], NR
        }' "$OUT/$arm.sorted.tsv"
}

echo "case=$CASE reps=$REPS steps=$STEPS startup=${startup_wall}s trace=$trace_status"
summarize baseline
summarize native
echo "artifacts=$OUT"

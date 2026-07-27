#!/usr/bin/env bash
# E1: MG coarsening-depth / bottom-solver sweep.  Configuration only -- no
# source change, so every arm runs the same binary.
#
#   mg_sweep.sh <binary> <reps> <arm-name> [extra parmparse args...]
#
# Records, per rep: wall seconds, per-Newton MLMG iteration counts, final
# resid/bnorm, and MLMG Solve/Bottom timers.  Convergence quality is recorded
# alongside wall time on purpose: a wall win bought by taking fewer/looser
# iterations is not a win, and MLMG_HIGH_CONTRAST_FINDINGS.md specifically
# warns that capping max_coarsening_level without a smoother bottom solver
# detonates the BiCGStab bottom.
set -euo pipefail

REPO=/home/jackplum/Projects/alamo
TASK="$REPO/docs/agent_plans/20260727-gpu-optimization-investigation"

BIN="${1:?usage: mg_sweep.sh <binary> <reps> <arm> [args...]}"
REPS="${2:?reps}"
ARM="${3:?arm name}"
shift 3
EXTRA=("$@")

cd "$REPO"
test -x "$BIN"
if [ -n "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ]; then
  echo "refusing to run: another CUDA process is resident" >&2
  exit 1
fi

# shellcheck source=/dev/null
source "$REPO/benchmark/local_cuda_env.sh"

OUT="$TASK/artifacts/e1_mg_sweep/$ARM"
mkdir -p "$OUT"

# Case selection matches the 20260726 task's harness so numbers are comparable.
# 2d  = 2D conservative (exercises the m_conservative_face_flux branch)
# 3ds = golden-suite 3D section; device-resident, unlike the 128^3 psi case
#       which oversubscribes an 8 GB A1000 and measures page migration.
INPUT=tests/ElasticSoftVoid/input
case "${ALAMO_CASE:-2d}" in
  2d)
    COMMON=(
      max_step=2 amr.max_level=2 amr.n_cell=64 64 8
      explicitmesh.lo1=32 32 0 explicitmesh.hi1=95 95 0
      explicitmesh.lo2=72 72 0 explicitmesh.hi2=183 183 0
      pf.eta.ic.expression.constant.w=0.002
      model_void.kappa=0.2_MPa model_void.mu=0.2_MPa
      elastic.print_model=0 elastic.print_residual=0
      elastic.solver.verbose=2 elastic.solver.nr_diagnostics=0
      amr.plot_int=-1 amr.thermo.plot_int=-1
    )
    ;;
  3ds)
    COMMON=(
      max_step=2 amr.n_cell=64 64 8
      explicitmesh.lo1=32 32 0 explicitmesh.hi1=95 95 15
      pf.eta.ic.expression.constant.w=0.002
      model_void.kappa=0.2_MPa model_void.mu=0.2_MPa
      elastic.solver.invariant_periodic=0 0 1
      elastic.print_model=0 elastic.print_residual=0
      elastic.solver.verbose=2 elastic.solver.nr_diagnostics=0
      amr.plot_int=-1 amr.thermo.plot_int=-1
    )
    ;;
  *) echo "unknown ALAMO_CASE=${ALAMO_CASE}" >&2; exit 2;;
esac

printf '%s\n' "${EXTRA[@]}" > "$OUT/arm_args.txt"

run_one() {
  # Declared separately: a single `local a=.. b="$a"` declares both names
  # before assigning, so $tag would still be unset when $d is evaluated.
  local tag="$1"
  local d="$OUT/$tag"
  mkdir -p "$d"
  local t0 t1
  t0=$(date +%s.%N)
  timeout 1200 "$BIN" "$INPUT" "${COMMON[@]}" "${EXTRA[@]}" \
      plot_file="$d/out" > "$d/stdout" 2> "$d/stderr" || echo "NONZERO_EXIT=$?" >> "$d/stderr"
  t1=$(date +%s.%N)
  echo "$t1 - $t0" | bc > "$d/wall_seconds.txt"
  rm -rf "$d/out"
}

run_one warmup
for r in $(seq 1 "$REPS"); do run_one "rep$r"; done

# --- summarise -------------------------------------------------------------
{
  echo "arm: $ARM"
  echo "args: ${EXTRA[*]:-<none>}"
  echo
  printf '%-8s %10s %14s %14s %12s\n' rep wall_s mlmg_iters final_rel_resid bottom_s
  for r in $(seq 1 "$REPS"); do
    d="$OUT/rep$r"
    w=$(cat "$d/wall_seconds.txt")
    it=$(grep -oP 'Final Iter\. \K[0-9]+' "$d/stdout" | paste -sd+ | bc)
    # Line format: "Final Iter. N resid, resid/bnorm = <abs_resid>, <rel_resid>"
    # -- take the SECOND number, which is the relative residual the solver
    # actually gates on.
    rr=$(grep -oP 'resid/bnorm = [0-9.e+-]+, \K[0-9.e+-]+' "$d/stdout" | tail -1)
    bt=$(grep -oP 'Bottom = \K[0-9.e+-]+' "$d/stdout" | paste -sd+ | bc)
    printf '%-8s %10.3f %14s %14s %12s\n' "rep$r" "$w" "${it:-NA}" "${rr:-NA}" "${bt:-NA}"
  done
  echo
  echo -n "median wall: "
  for r in $(seq 1 "$REPS"); do cat "$OUT/rep$r/wall_seconds.txt"; done \
    | sort -n | awk '{a[NR]=$1} END {print (NR%2 ? a[(NR+1)/2] : (a[NR/2]+a[NR/2+1])/2)}'
} | tee "$OUT/summary.txt"

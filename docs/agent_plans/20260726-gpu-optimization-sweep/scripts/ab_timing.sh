#!/usr/bin/env bash
# Alternating A/B wall-time harness for the 2026-07-26 GPU optimization sweep.
#
#   ab_timing.sh <case: 2d|3d> <label> <binA> <binB> [reps]
#
# Runs warmup + N alternating reps per arm (A first, then B, inside each rep so
# thermal drift hits both arms equally), one process at a time, refusing to run
# if another CUDA process is resident.  Wall seconds land in
# artifacts/<label>/<case>/<arm>/<rep>/wall_seconds.txt.
set -euo pipefail

REPO=/home/jackplum/Projects/alamo
TASK="$REPO/docs/agent_plans/20260726-gpu-optimization-sweep"

CASE="${1:?usage: ab_timing.sh 2d|3d label binA binB [reps]}"
LABEL="${2:?label}"
BIN_A="${3:?binA}"
BIN_B="${4:?binB}"
REPS="${5:-5}"

cd "$REPO"
test -x "$BIN_A"
test -x "$BIN_B"
test -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)"

case "$CASE" in
  2d)
    INPUT=tests/ElasticSoftVoid/input
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
    ;;
  3ds)
    # Golden-suite 3D section: device-resident, so it is free of the managed
    # oversubscription noise the 128^3 psi case suffers on an 8 GB A1000.
    INPUT=tests/ElasticSoftVoid/input
    COMMON=(
      max_step=2 amr.n_cell=64 64 8
      explicitmesh.lo1=32 32 0 explicitmesh.hi1=95 95 15
      pf.eta.ic.expression.constant.w=0.002
      model_void.kappa=0.2_MPa model_void.mu=0.2_MPa
      elastic.solver.invariant_periodic=0 0 1
      elastic.print_model=0 elastic.print_residual=0
      elastic.solver.verbose=0 elastic.solver.nr_diagnostics=0
      amr.plot_int=-1 amr.thermo.plot_int=-1
    )
    ;;
  3d)
    INPUT=input_3d_centre_bore_128_a2
    COMMON=(
      max_step=2 stop_time=1e99_s elastic.interval=1
      amr.max_grid_size=64 amrex.the_arena_is_managed=1
      elastic.solver.nriters=2 elastic.print_model=0
      elastic.solver.verbose=0 elastic.solver.nr_diagnostics=0
      amr.plot_int=-1 amr.thermo.plot_int=-1
    )
    ;;
  *) echo "unknown case $CASE" >&2; exit 2;;
esac

OUT="$TASK/artifacts/$LABEL/$CASE"
mkdir -p "$OUT"
sha256sum "$BIN_A" "$BIN_B" "$INPUT" > "$OUT/inputs.sha256"

run_one () {
  local arm="$1" label="$2" binary="$3"
  local dir="$OUT/$arm/$label"
  mkdir -p "$dir"
  /usr/bin/time -f %e -o "$dir/wall_seconds.txt" \
    "$binary" "$INPUT" "${COMMON[@]}" plot_file="$dir/plot" \
    >"$dir/stdout" 2>"$dir/stderr"
  printf '%-4s %-8s %s\n' "$arm" "$label" "$(tr -d '\n' <"$dir/wall_seconds.txt")"
}

for label in warmup $(seq -f 'rep%g' 1 "$REPS"); do
  run_one A "$label" "$BIN_A"
  run_one B "$label" "$BIN_B"
done

python3 - "$OUT" "$REPS" <<'PY'
import sys, pathlib, statistics
out, reps = pathlib.Path(sys.argv[1]), int(sys.argv[2])
res = {}
for arm in ("A", "B"):
    vals = []
    for r in range(1, reps + 1):
        f = out / arm / f"rep{r}" / "wall_seconds.txt"
        vals.append(float(f.read_text().strip()))
    med = statistics.median(vals)
    mad = statistics.median([abs(v - med) for v in vals])
    res[arm] = (med, mad, vals)
    print(f"{arm}: median {med:.4f} s  MAD {mad:.4f}  raw {vals}")
a, b = res["A"][0], res["B"][0]
print(f"B vs A: {100*(b-a)/a:+.3f}%  (negative = B faster)")
PY

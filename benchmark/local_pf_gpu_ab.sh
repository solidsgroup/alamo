#!/usr/bin/env bash
# Local phase-field GPU A/B harness.
#
# Usage:
#   BASE_BIN=/path/to/base/alamo_gpu OPT_BIN=/path/to/opt/alamo_gpu \
#     benchmark/local_pf_gpu_ab.sh [outdir]
#
# This intentionally benchmarks only the flame/phase-field path. It uses the
# existing GPU perf inputs with elastic.type=disable and leaves NOVA for large
# version-scale runs.
set -euo pipefail

cd "$(dirname "$0")/.." || exit 1

BASE_BIN="${BASE_BIN:-}"
OPT_BIN="${OPT_BIN:-}"
if [[ -z "$BASE_BIN" || -z "$OPT_BIN" ]]; then
  echo "usage: BASE_BIN=/path/base OPT_BIN=/path/opt benchmark/local_pf_gpu_ab.sh [outdir]" >&2
  exit 2
fi
if [[ ! -x "$BASE_BIN" ]]; then
  echo "BASE_BIN is not executable: $BASE_BIN" >&2
  exit 2
fi
if [[ ! -x "$OPT_BIN" ]]; then
  echo "OPT_BIN is not executable: $OPT_BIN" >&2
  exit 2
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
OUT="${1:-benchmark/local_ab_${STAMP}}"
REPEATS="${REPEATS:-3}"
P2_STEPS="${P2_STEPS:-30}"
P1_STEPS="${P1_STEPS:-60}"
RUN_NSYS="${RUN_NSYS:-1}"
PERF_GATE="${PERF_GATE:-0}"
GATE_MIN_P2_SPEEDUP="${GATE_MIN_P2_SPEEDUP:-1.10}"

P2_INPUT="tests/GPU/P2_perf_2d_hiRes_AMR3/input"
P1_INPUT="tests/GPU/P1_perf_2d_hiRes_noAMR/input"

mkdir -p "$OUT"

if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi --query-gpu=name,compute_cap,memory.total --format=csv,noheader \
    | head -n 1 | sed 's/^/local_gpu=/' > "$OUT/environment.txt"
else
  echo "local_gpu=<nvidia-smi not found>" > "$OUT/environment.txt"
fi
{
  echo "base=$BASE_BIN"
  echo "opt=$OPT_BIN"
  echo "p2_input=$P2_INPUT"
  echo "p1_input=$P1_INPUT"
  echo "p2_steps=$P2_STEPS"
  echo "p1_steps=$P1_STEPS"
  echo "repeats=$REPEATS"
  echo "perf_gate=$PERF_GATE"
  echo "gate_min_p2_speedup=$GATE_MIN_P2_SPEEDUP"
} >> "$OUT/environment.txt"

WALL="$OUT/wall_repeats.csv"
echo "case,label,repeat,wall_s" > "$WALL"

run_wall_one() {
  local case_name="$1"
  local label="$2"
  local bin="$3"
  local input="$4"
  local max_step="$5"
  local repeat="$6"
  local run_dir="$OUT/wall_${case_name}_${label}_${repeat}"
  mkdir -p "$run_dir"

  /usr/bin/time -f "%e" -o "$run_dir/time.raw" \
    "$bin" "$input" \
      max_step="$max_step" allow_unused=1 \
      plot_file="$run_dir/plot" amr.plot_int=-1 amr.thermo.plot_int=-1 \
      amrex.async_out=0 tiny_profiler.device_synchronize_around_region=0 \
      > "$run_dir/run.log" 2>&1

  printf "%s,%s,%s,%s\n" "$case_name" "$label" "$repeat" "$(cat "$run_dir/time.raw")" >> "$WALL"
}

for repeat in $(seq 1 "$REPEATS"); do
  run_wall_one p2_amr_thermal_on base "$BASE_BIN" "$P2_INPUT" "$P2_STEPS" "$repeat"
  run_wall_one p2_amr_thermal_on opt  "$OPT_BIN"  "$P2_INPUT" "$P2_STEPS" "$repeat"
  run_wall_one p1_noamr_thermal_on base "$BASE_BIN" "$P1_INPUT" "$P1_STEPS" "$repeat"
  run_wall_one p1_noamr_thermal_on opt  "$OPT_BIN"  "$P1_INPUT" "$P1_STEPS" "$repeat"
done

NSYS_BIN="${NSYS_BIN:-}"
if [[ -z "$NSYS_BIN" ]]; then
  if command -v nsys >/dev/null 2>&1; then
    NSYS_BIN="$(command -v nsys)"
  elif [[ -x ".local/nsight/opt/nvidia/nsight-systems/2026.1.3/bin/nsys" ]]; then
    NSYS_BIN=".local/nsight/opt/nvidia/nsight-systems/2026.1.3/bin/nsys"
  elif [[ -x "../alamo/.local/nsight/opt/nvidia/nsight-systems/2026.1.3/bin/nsys" ]]; then
    NSYS_BIN="../alamo/.local/nsight/opt/nvidia/nsight-systems/2026.1.3/bin/nsys"
  fi
fi

if [[ "$RUN_NSYS" == "1" && -n "$NSYS_BIN" ]]; then
  for label in base opt; do
    bin="$BASE_BIN"
    [[ "$label" == "opt" ]] && bin="$OPT_BIN"
    run_dir="$OUT/nsys_p2_${label}"
    mkdir -p "$run_dir"
    /usr/bin/time -f "wall_s %e" -o "$run_dir/time.txt" \
      "$NSYS_BIN" profile --force-overwrite true --trace=cuda,nvtx \
        --sample=none --cpuctxsw=none --output "$run_dir/nsys_report" \
        "$bin" "$P2_INPUT" \
          max_step="$P2_STEPS" allow_unused=1 \
          plot_file="$run_dir/plot" amr.plot_int=-1 amr.thermo.plot_int=-1 \
          amrex.async_out=0 tiny_profiler.device_synchronize_around_region=0 \
          > "$run_dir/run.log" 2>&1
    "$NSYS_BIN" stats --report cuda_api_sum,cuda_gpu_kern_sum --format csv \
      --output "$run_dir/nsys_stats" "$run_dir/nsys_report.nsys-rep" \
      > "$run_dir/nsys_stats_stdout.txt" 2>&1 || true
  done
fi

python3 - "$OUT" "$PERF_GATE" "$GATE_MIN_P2_SPEEDUP" <<'PY'
import csv
import pathlib
import statistics
import sys

out = pathlib.Path(sys.argv[1])
perf_gate = sys.argv[2] == "1"
gate_min_p2_speedup = float(sys.argv[3])
rows = list(csv.DictReader((out / "wall_repeats.csv").open()))
summary_rows = []
print("case,label,n,median_s,mean_s,stdev_s")
for case in sorted({row["case"] for row in rows}):
    for label in ("base", "opt"):
        samples = [float(row["wall_s"]) for row in rows
                   if row["case"] == case and row["label"] == label]
        stdev = statistics.stdev(samples) if len(samples) > 1 else 0.0
        summary_rows.append({
            "case": case,
            "label": label,
            "n": len(samples),
            "median_s": statistics.median(samples),
            "mean_s": statistics.mean(samples),
            "stdev_s": stdev,
        })
        print(f"{case},{label},{len(samples)},{statistics.median(samples):.3f},"
              f"{statistics.mean(samples):.3f},{stdev:.3f}")

with (out / "summary.csv").open("w", newline="") as fh:
    writer = csv.DictWriter(
        fh, fieldnames=["case", "label", "n", "median_s", "mean_s", "stdev_s"])
    writer.writeheader()
    for row in summary_rows:
        writer.writerow(row)

medians = {(row["case"], row["label"]): row["median_s"] for row in summary_rows}
speedups = {}
for case in sorted({row["case"] for row in summary_rows}):
    base = medians.get((case, "base"))
    opt = medians.get((case, "opt"))
    if base and opt:
        speedups[case] = base / opt

with (out / "speedups.csv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=["case", "speedup"])
    writer.writeheader()
    for case, speedup in sorted(speedups.items()):
        writer.writerow({"case": case, "speedup": speedup})

if speedups:
    print("case,speedup")
    for case, speedup in sorted(speedups.items()):
        print(f"{case},{speedup:.3f}")

p2_speedup = speedups.get("p2_amr_thermal_on")
if perf_gate and (p2_speedup is None or p2_speedup < gate_min_p2_speedup):
    found = "missing" if p2_speedup is None else f"{p2_speedup:.3f}"
    print(f"FAIL: p2_amr_thermal_on speedup {found} < {gate_min_p2_speedup:.3f}")
    sys.exit(1)
PY

echo "artifacts=$OUT"

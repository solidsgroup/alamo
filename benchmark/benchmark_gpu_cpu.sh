#!/usr/bin/env bash
# ============================================================================
# benchmark_gpu_cpu.sh -- compare CPU vs GPU wall-clock + holdups for the
# Flame/chamber solver, and emit performance flame graphs for each run.
#
# What it does
#   1. Runs the same input on the CPU binary and (if present) the CUDA binary.
#   2. Captures AMReX TinyProfiler tables (per-region inclusive/exclusive
#      wall-clock -> the holdup breakdown: Advance vs WritePlotFile vs elastic
#      solve vs Regrid).
#   3. Wraps the GPU run in NVIDIA Nsight Systems (nsys) when available for a
#      CPU+GPU timeline, and renders flame graphs from the profiler data.
#   4. Prints a side-by-side wall-clock / holdup comparison.
#
# Usage
#   benchmark/benchmark_gpu_cpu.sh <input> [max_step]
#
# Prereqs for full output (degrade gracefully if missing):
#   - CPU binary built with profiling:   ./configure --profile && make
#   - GPU binary built with profiling:   ./configure --cuda --profile && make
#   - nsys (CUDA toolkit) for GPU flame graphs
#   - Brendan Gregg's FlameGraph (flamegraph.pl on PATH) for CPU flame graphs
# ============================================================================
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

INPUT="${1:?usage: benchmark_gpu_cpu.sh <input> [max_step]}"
MAXSTEP="${2:-50}"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT="benchmark/results_${STAMP}"
mkdir -p "$OUT"

# Locate the most recent CPU and CUDA alamo binaries.
CPU_BIN="$(ls -t bin/alamo-2d-*-* 2>/dev/null | grep -v cuda | head -1)"
GPU_BIN="$(ls -t bin/alamo-2d-cuda* 2>/dev/null | head -1)"

NP="${NP:-1}"
RUN="mpiexec -np ${NP}"
# Common overrides: short run, profiler on, async IO on.
OV="max_step=${MAXSTEP} amrex.async_out=1 tiny_profiler.device_synchronize_around_region=1"

echo "=== Flame CPU/GPU benchmark ==="
echo "input=${INPUT}  max_step=${MAXSTEP}  results=${OUT}"
echo "CPU binary: ${CPU_BIN:-<none>}"
echo "GPU binary: ${GPU_BIN:-<none>}"

run_one () {  # $1=label  $2=binary  $3=extra-launch-prefix
  local label="$1" bin="$2" pre="${3:-}"
  [ -z "$bin" ] && { echo "[$label] no binary, skipping"; return; }
  local log="${OUT}/${label}.log"
  echo "--- running ${label} ---"
  local t0 t1
  t0=$(date +%s.%N)
  # shellcheck disable=SC2086
  ${pre} ${RUN} "${bin}" "${INPUT}" ${OV} plot_file="${OUT}/out_${label}" >"${log}" 2>&1
  t1=$(date +%s.%N)
  echo "${label} wall_clock_s $(echo "$t1 - $t0" | bc)" | tee -a "${OUT}/summary.txt"
  # Extract the TinyProfiler inclusive table (between the header and the next blank-after-table)
  awk '/TinyProfiler total time|Inclusive/{f=1} f{print} /^$/{if(f>1)f=0; if(f)f++}' "${log}" \
      > "${OUT}/${label}.tinyprof.txt" 2>/dev/null
  echo "  log:      ${log}"
  echo "  profiler: ${OUT}/${label}.tinyprof.txt"
}

# ---- CPU run ----
run_one cpu "${CPU_BIN}" ""

# ---- GPU run (under nsys if available) ----
if [ -n "${GPU_BIN}" ]; then
  if command -v nsys >/dev/null 2>&1; then
    run_one gpu "${GPU_BIN}" "nsys profile --force-overwrite true -o ${OUT}/gpu_trace --stats=true"
    echo "  nsys:     ${OUT}/gpu_trace.nsys-rep  (open in Nsight Systems for the GPU/CPU timeline)"
  else
    echo "[gpu] nsys not found -- running without GPU timeline capture"
    run_one gpu "${GPU_BIN}" ""
  fi
fi

# ---- Flame graphs from TinyProfiler tables (flat, cross-platform) ----
for label in cpu gpu; do
  tp="${OUT}/${label}.tinyprof.txt"
  [ -s "$tp" ] || continue
  folded="${OUT}/${label}.folded"
  python3 benchmark/tinyprofiler_to_folded.py "$tp" > "$folded" 2>/dev/null
  if command -v flamegraph.pl >/dev/null 2>&1 && [ -s "$folded" ]; then
    flamegraph.pl --title "Flame ${label} (TinyProfiler)" "$folded" > "${OUT}/${label}.svg"
    echo "  flamegraph: ${OUT}/${label}.svg"
  fi
done

# ---- Comparison ----
echo ""
echo "=== wall-clock summary ==="
cat "${OUT}/summary.txt" 2>/dev/null
echo ""
echo "=== per-region holdup comparison (CPU vs GPU, top regions) ==="
python3 benchmark/compare_tinyprofiler.py "${OUT}/cpu.tinyprof.txt" "${OUT}/gpu.tinyprof.txt" 2>/dev/null \
  || echo "(comparison needs both cpu/gpu profiler tables)"
echo ""
echo "Done. Artifacts in ${OUT}/"

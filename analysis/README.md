# ALAMO CPU vs GPU Performance Suite

A modular, dependency-light profiling suite that produces a full efficiency,
wall-clock and **I/O** breakdown of the CPU and GPU alamo builds, including
**flame graphs**, and rolls everything into a single HTML report.

## Quick start

```bash
# After the two production runs finish (full 1.5 s, plotting off):
#   out_cpu_star_20mpa.log   (CPU)
#   out_gpu_star_20mpa.log   (GPU)
analysis/run_all.sh
xdg-open analysis/results/index.html
```

Run a subset:

```bash
analysis/run_all.sh --phases 1,6          # just wall-clock + report (no re-runs)
analysis/run_all.sh --skip-flame --skip-perf
analysis/run_all.sh --profile-steps 120 --plot-int 30
```

## What it measures

| Phase | Script | Tool | Output |
|---|---|---|---|
| 1 Wall-clock & efficiency | `01_wallclock.py` | `/usr/bin/time -v` footer + log MLMG timers | `wallclock.{md,csv,json}`, `chart_wallclock.svg`, `chart_elastic.svg` |
| 2 **I/O profile** | `02_io_profile.sh` | `strace -f -c` | `io_profile.{md,json}`, `chart_io.svg` |
| 3 Microarchitecture | `03_perf_stat.sh` | `perf stat` | `perfstat.{md,json}` |
| 4 **CPU flame graph** | `04_flamegraph_cpu.sh` | `perf record` + FlameGraph | `flamegraph_cpu.svg`, `flamegraph_gpu.svg` |
| 5 GPU device timeline | `05_gpu_timeline.sh` | `nsys` *or* `nvidia-smi` sampling | `gpu_timeline.{md,json,csv}`, `chart_gpu_timeline.svg` |
| 6 Consolidated report | `06_report.py` | — | `REPORT.md`, `index.html` |

### Design principles (architect notes)

- **Production vs instrumented separation.** Phase 1 reads the *real* full-length
  runs untouched. Phases 2–5 launch *short* re-runs (`PROFILE_STEPS`, default 80
  so ≥1 elastic solve is captured) under profilers — profiling perturbs timing,
  so it never contaminates the headline wall-clock.
- **I/O is measured with plotting ON.** A plotting-off run has ~no I/O. Phases
  2/4 enable `plot_int` so the VisMF write + `fsync` path is real and timed.
  Phase 1 keeps plotting off for a clean compute comparison.
- **Zero heavy deps.** Charts are hand-rolled SVG (no matplotlib). Parsers are
  stdlib Python. The suite runs on a bare login node.
- **Graceful degradation.** A missing tool downgrades one artifact, never aborts
  the run. Everything that *can* be produced, is.
- **Throughput, not just totals.** Reports wall-per-step and per-solve cost so
  CPU/GPU are compared on equal work even if step counts differ.

## Time spent locked in I/O

Two complementary views:

1. **`/usr/bin/time -v`** (phase 1): `non_compute_s` = `wall − (user+sys)`, a
   direct lower bound on time the (single-thread) process was *not* on-CPU —
   i.e. blocked on I/O, `fsync`, MPI or device sync.
2. **`strace -f -c`** (phase 2): exact wall time *inside* each syscall, bucketed
   into `file_write / file_read / sync / metadata / gpu_driver(ioctl) / memory`,
   with the I/O fraction of in-kernel time. This is the privilege-free,
   syscall-accurate "time locked in I/O" number.

## Flame graphs

- **CPU**: real `perf record --call-graph dwarf` → Brendan Gregg FlameGraph SVG
  (auto-vendored to `analysis/vendor/FlameGraph`). Needs sampling permission:

  ```bash
  ! sudo sysctl kernel.perf_event_paranoid=1
  ! sudo sysctl kernel.kptr_restrict=0      # for kernel symbols
  ```

  (This box defaults to `perf_event_paranoid=4`, which blocks sampling; the
  script detects this and prints the fix. `dwarf` unwinding is used because the
  `-O2` binaries omit frame pointers.)

- **GPU device kernel flame chart**: requires **Nsight Systems** (`nsys`), which
  is *not installed here*. Phase 5 auto-uses it if present (`gpu_nsys.nsys-rep`,
  open in the Nsight GUI for the kernel timeline). Without it, phase 5 falls
  back to an `nvidia-smi` utilization/power/clock **timeline** + the CUDA-driver
  `ioctl` time from phase 2, which together characterize device occupancy and
  host-bound stalls. To install:

  ```bash
  # one option: the standalone CLI
  #   https://developer.nvidia.com/nsight-systems  (apt: nsight-systems-cli)
  ```

## Configuration (env vars)

| Var | Default | Meaning |
|---|---|---|
| `BIN_CPU` | `bin/alamo-2d-clang++` | CPU binary |
| `BIN_GPU` | `bin/alamo_gpu-2d-cuda86-g++` | GPU binary |
| `INPUT` | `input_copy` | input deck |
| `CPU_LOG` / `GPU_LOG` | `out_cpu_star_20mpa.log` / `out_gpu_star_20mpa.log` | production logs for phase 1 |
| `PROFILE_STEPS` | `80` | steps for instrumented re-runs |
| `PLOT_INT` | `40` | plot interval for I/O-aware runs |
| `PERF_FREQ` | `997` | perf sampling Hz |
| `SAMPLE_MS` | `100` | nvidia-smi sampling period |

## Layout

```
analysis/
├── run_all.sh              orchestrator
├── 01_wallclock.py         phase 1
├── 02_io_profile.sh        phase 2
├── 03_perf_stat.sh         phase 3
├── 04_flamegraph_cpu.sh    phase 4 (CPU + GPU-host)
├── 05_gpu_timeline.sh      phase 5
├── 06_report.py            phase 6
├── lib/                    parsers + common.sh + svgchart.py
├── vendor/FlameGraph       auto-cloned
└── results/                all artifacts + index.html
```

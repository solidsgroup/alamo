# Task 001 (backfilled): the analysis/ six-phase profiling suite

## Task
Document the `analysis/` directory rewrite as a task record: a dependency-light,
six-phase CPU-vs-GPU profiling suite that produces wall-clock, I/O, microarchitecture,
flame-graph, and GPU-device-timeline breakdowns, rolled into one HTML/MD report.

## Why
Built 2026-06-19 through 2026-06-20 to replace ad hoc one-off profiling commands with
a repeatable harness, so every perf claim in `docs/llm/perf/` and `changelog/` entries
can cite a script + output file instead of a one-off shell transcript. Backing data
for `docs/llm/perf/2026-06-20-gpu-port-report.md` (the coarse/deep/wide campaigns)
and `analysis/results_phase22*` (Phase 2.2 nsys captures) came from this suite or its
direct predecessor scripts.

## What it is
Six phases, each independently skippable, orchestrated by `analysis/run_all.sh`:

| Phase | Script | Tool | Output |
|---|---|---|---|
| 1 Wall-clock & efficiency | `01_wallclock.py` | `/usr/bin/time -v` + log MLMG timers | `wallclock.{md,csv,json}`, charts |
| 2 I/O profile | `02_io_profile.sh` | `strace -f -c` | `io_profile.{md,json}` |
| 3 Microarchitecture | `03_perf_stat.sh` | `perf stat` | `perfstat.{md,json}` |
| 4 CPU flame graph | `04_flamegraph_cpu.sh` | `perf record` + FlameGraph | `flamegraph_cpu.svg` |
| 5 GPU device timeline | `05_gpu_timeline.sh` | `nsys` or `nvidia-smi` sampling | `gpu_timeline.*` |
| 6 Consolidated report | `06_report.py` | — | `REPORT.md`, `index.html` |

Support: `analysis/lib/` (stdlib-only parsers + `common.sh` + hand-rolled `svgchart.py`,
no matplotlib dependency), `analysis/vendor/FlameGraph` (auto-cloned on first flame-graph
run). Full usage and design rationale (production-vs-instrumented run separation,
I/O-with-plotting-on, graceful degradation when a tool is missing): `analysis/README.md`
— this task record does not duplicate it, see that file for the how-to.

## Files involved (created, still untracked in git)
`analysis/{01_wallclock.py,02_io_profile.sh,03_perf_stat.sh,04_flamegraph_cpu.sh,05_gpu_timeline.sh,06_report.py,run_all.sh,README.md}`, `analysis/lib/*`.

## Test command (structural only — no production logs present on this box to run it end-to-end)
```bash
for f in analysis/run_all.sh analysis/0{2,3,4,5}_*.sh analysis/lib/common.sh; do bash -n "$f"; done
python3 -m py_compile analysis/01_wallclock.py analysis/06_report.py analysis/lib/*.py
```

## Acceptance criteria
All scripts parse (`bash -n` / `py_compile`) with no errors. Re-verified during this
backfill (2026-06-22): all pass.

## Write results to: results/001-RESULT.md

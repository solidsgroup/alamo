# Task 003: Memory-budget calculator + scaling/crossover harness + R3 skeleton (roadmap 3.2, 3.5)

## Goal
Produce the analysis tooling that drives the Phase-3 crossover study on NOVA:
(a) a memory-budget calculator that prints the largest 3D grid that fits per GPU,
(b) a scaling/crossover sweep driver, and (c) the R3 crossover-report skeleton.
All must run on a machine WITHOUT a GPU (pure Python / bash; no CUDA, no SLURM submit).

## Context (only what this task needs)
- Roadmap 3.2 (memory budgeting): compute bytes/node and the largest grid that fits on a
  target device. Recon figures: elastic model fab at `Matrix4<3>` ≈ 360 B/node; phase-field
  fields eta+phi+temp ≈ 3 doubles (24 B); displacement vector 3 doubles (24 B); plus
  RHS/residual/stress scratch. Since Phase-3 GPU runs disable elastic (D1), the dominant
  resident cost is the phase-field/thermal fields + AMReX overhead, NOT the elastic model
  fab. Make bytes/node a CLI parameter with a sensible default (e.g. 256 B/node for the
  elastic-disabled phase-field config; document the breakdown) and also accept the
  ~512 B/node elastic-enabled figure for comparison. Account for ghost cells (a multiplier).
- Target devices and memory (use as the built-in device table):
  A1000=8 GB, A100-40=40 GB, A100-80=80 GB, H100-80=80 GB, H200=141 GB.
- Roadmap 3.5 (crossover hunt): sweep problem size × GPU count to find where the GPU build
  beats the CPU node. The sweep drives the task-001 inputs `input_3d_flame_{128,256,512}`
  and the task-002 NOVA scripts (`nova_flame_gpu_3d.slurm`, `nova_flame_gpu_3d_multi.slurm`,
  `nova_flame_cpu_3d.slurm`). Reference these by name; do NOT create them.
- Standing metric set (roadmap, capture identically every run): launches/step,
  kernel-duration avg, `cudaStreamSynchronize` fraction, achieved occupancy + regs (ncu),
  wall/step GPU vs CPU node, golden-compare residual.

## Files to read first
- `docs/agent_plans/20260620-gpu-phase3-regime-scaling/PLAN.md`
- `GPU-OPT-ROADMAP.txt` is on the user's Desktop (NOT in repo); rely on the PLAN summary
  above instead of trying to open it. (If `benchmark/G0_BASELINE_OF_RECORD.md` exists, skim
  it for the existing metric/reporting conventions to match.)

## Files allowed to modify (create only)
- `benchmark/phase3_memory_budget.py` — bytes/node → max grid per device
- `benchmark/phase3_scaling_sweep.sh`  — size×GPU sweep driver (NOVA-staged)
- `benchmark/PHASE3_R3_crossover.md`    — R3 report skeleton (tables to be filled on NOVA)

## Files NOT allowed to modify
- `src/`, `bin/`, repo-root inputs, any existing `benchmark/*` file, files owned by tasks 001/002.

## Implementation steps
1. `phase3_memory_budget.py` (Python 3, stdlib only, argparse):
   - Args: `--bytes-per-node` (default 256), `--ghost-factor` (default 1.3),
     `--max-level` (default 1), `--refine-fill` (fraction of domain refined, default 0.1),
     `--device` (optional filter).
   - For each device in the table, compute the largest cubic-ish base grid N (with the
     wide-shallow aspect, e.g. Nx=Ny=2*Nz) whose total node memory (incl. ghost factor and
     an AMR overhead term for the refined fraction across levels) fits in device RAM with a
     safety headroom (e.g. 85%). Print a table: device | mem | max Nx×Ny×Nz | #nodes | est. bytes.
   - `--help` must work; running with no args prints the default table. Keep it readable and
     documented; this is a planning tool, approximations are fine and should be labeled as such.
2. `phase3_scaling_sweep.sh` (bash, `set -euo pipefail`, NOVA-staged — does NOT submit here):
   - Define a SIZES list (`128 256 512`) mapped to inputs `input_3d_flame_<size>` and a
     GPUS list (`1 2 4 8`). Emit the matrix of `sbatch` commands (single-GPU script for
     GPUS=1, multi-GPU script for GPUS>1, plus the CPU-node baseline) with `INPUT=...`,
     `NGPUS=...`, `GPU_TYPE=...` env, into a `commands.txt` (dry-run by default; `--submit`
     to actually `sbatch`). Strong-scaling = fixed size, vary GPUS; weak-scaling = size grows
     with GPUS. Document both modes via a `MODE=strong|weak` env (default strong).
   - Echo a clear banner that this is staged for NOVA and is a no-op without SLURM.
3. `PHASE3_R3_crossover.md`: the R3 report skeleton per the roadmap — purpose, method,
   a problem-size × hardware table (columns: size, device, #GPUs, wall/step GPU, wall/step
   CPU-node, GPU/CPU ratio, launches/step, sync fraction, occupancy), and a "Crossover /
   verdict" section with the three valid D3 outcomes (WIN @ single | WIN @ scale | NO WIN).
   Leave numeric cells as TODO to be filled on NOVA.

## Invariants (device, concurrency, architecture)
- Pure planning tooling: no GPU, no CUDA, no live SLURM submission by default.
- Reference inputs/scripts owned by tasks 001/002 by name only; never create them.
- bytes/node is parameterized, not hardcoded; label all estimates as approximate.

## Build and test commands
- `python3 benchmark/phase3_memory_budget.py --help` and a no-arg run must succeed and print
  a device table. `bash -n benchmark/phase3_scaling_sweep.sh` and a default (dry-run)
  invocation must succeed and write `commands.txt` without submitting.

## Expected result
A working memory-budget calculator, a dry-run-safe sweep driver, and an R3 skeleton ready
to fill with NOVA numbers.

## Non-goals
- No real measurements (those happen on NOVA). No editing of inputs or SLURM scripts.
- Do not over-engineer the budget model — a documented first-order estimate is the goal.

## Stop conditions
- None expected; if `python3` is unavailable, record it and still write the bash + md files.

## Final report: write results/003-RESULT.md (RESULT template).

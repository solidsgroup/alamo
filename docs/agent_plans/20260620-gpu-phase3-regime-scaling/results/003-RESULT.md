# Task 003 RESULT: Memory-budget calculator + scaling/crossover harness + R3 skeleton

## Summary
Created the Phase-3 planning tooling that drives the NOVA crossover study (roadmap
3.2 + 3.5): a stdlib-only memory-budget calculator, a dry-run-safe size x GPU
scaling sweep driver, and the R3 crossover-report skeleton. All three run on a
machine WITHOUT a GPU: no CUDA, no SLURM submission by default. No existing file
was modified; only the three new files listed below were created.

- `phase3_memory_budget.py`: from a parameterized bytes/node figure (default 256 B/node
  elastic-disabled/D1; pass `--bytes-per-node 512` for the elastic-enabled comparison),
  ghost factor, and a shallow-AMR overhead term, prints the largest wide-shallow 3D
  grid (Nx=Ny=2*Nz, blocking_factor=32 steps) that fits per device for the built-in
  table A1000 (8 GB), A100-40, A100-80, H100-80, H200 (141 GB). All numbers labelled
  APPROXIMATE / planning-only.
- `phase3_scaling_sweep.sh`: emits the `sbatch` command matrix into `commands.txt`
  (dry-run by default; `--submit` to actually sbatch on a SLURM host). Supports
  `MODE=strong|weak`, `GPU_TYPE`, `SIZES`, `GPUS`, `OUT`. References task-001 inputs
  (`input_3d_flame_{128,256,512}`) and task-002 SLURM scripts
  (`nova_flame_gpu_3d.slurm`, `nova_flame_gpu_3d_multi.slurm`, `nova_flame_cpu_3d.slurm`)
  BY NAME ONLY; creates none of them. Prints a clear NOVA-staged / no-op-without-SLURM banner.
- `PHASE3_R3_crossover.md`: R3 skeleton with purpose, method, the standing metric set,
  a problem-size x hardware table (size, device, #GPUs, wall/step GPU, wall/step CPU-node,
  GPU/CPU ratio, launches/step, sync fraction, occupancy), and a Crossover/verdict section
  with the three valid D3 outcomes (WIN @ single | WIN @ scale | NO WIN). Numeric cells = TODO.

## Files changed
Created (new files only):
- `benchmark/phase3_memory_budget.py` (executable)
- `benchmark/phase3_scaling_sweep.sh` (executable)
- `benchmark/archive/PHASE3_R3_crossover.md`

No existing files modified. The generated `commands.txt` is a transient sweep output
and was removed after validation (regenerated on demand).

## Tests run and results
- `python3 benchmark/phase3_memory_budget.py --help` -> PASS (argparse help prints).
- `python3 benchmark/phase3_memory_budget.py` (no args) -> PASS (prints default device
  table; e.g. A1000 -> 256x256x128, A100-80 -> 576x576x288, H200 -> 704x704x352, all "ok").
- `bash -n benchmark/phase3_scaling_sweep.sh` -> PASS (syntax OK).
- `bash benchmark/phase3_scaling_sweep.sh` (default dry-run) -> PASS (writes 15-command
  matrix to commands.txt; submits nothing).
- Extra robustness checks (all PASS): sweep `--help`; `MODE=weak GPU_TYPE=h200` dry-run
  (size grows with GPUS); `--bytes-per-node 512 --device A100-80` (elastic-enabled
  comparison -> 448x448x224); bad-arg handling exits 2 for both tools; unknown-device
  filter exits 2 with a helpful message.
- `git status --short` confirms the three files are untracked (`??`); no existing file shows
  as modified.

## Issues found
None. Python 3.12.3 available locally (stdlib-only tool runs cleanly).

## Deviations from task
None. bytes/node parameterized (not hardcoded); estimates labelled approximate; inputs and
SLURM scripts referenced by name only; dry-run is the default and is a no-op without SLURM.

## Follow-up needed
- On NOVA: run `phase3_scaling_sweep.sh --submit` once task-001 inputs and task-002 SLURM
  scripts exist, then fill the R3 table and pick the D3 verdict.
- The bytes/node default (256) and ghost/AMR overhead model are first-order estimates; the
  lead's actual measured bytes/node from the 3D smoke run can refine the `--bytes-per-node`
  default if it differs materially.
- `weak_size_for_gpus` caps at 512 (largest provided input) for GPUS>=4; if larger inputs
  are added later, extend the map.

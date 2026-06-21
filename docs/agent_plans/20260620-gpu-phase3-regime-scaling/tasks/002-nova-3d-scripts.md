# Task 002: NOVA 3D build + single/multi-GPU launch scripts (roadmap 3.4)

## Goal
Create 3D siblings of the existing NOVA build/run scripts so we can build a 3D GPU binary
for A100 (sm_80) and H200 (sm_90) on NOVA and launch single- and multi-GPU 3D Flame runs.
These are staged for NOVA; they are NOT run here (no NOVA, single local GPU).

## Context (only what this task needs)
- Existing 2D NOVA scripts (READ as templates, do not modify):
  - `benchmark/build_alamo_nova.sh` — clones `chamber-gpu`, primes AMReX on login node,
    submits a SLURM CPU build job that runs `./configure --comp=g++ --cuda <arch> --profile`
    for each arch in `ARCHES="80 90"`, producing `bin/alamo_gpu-2d-profile-cuda{80,90}-g++`.
  - `benchmark/nova_flame_gpu.slurm` — loads cuda/gcc/openmpi modules; selects binary via
    `ls -t bin/alamo_gpu-2d*cuda${ARCH}*-g++`; runs `srun --gpus-per-task=1`; honors
    `GPU_TYPE=a100|h200`.
  - `benchmark/nova_flame_cpu.slurm` — CPU baseline; binary pattern `bin/alamo-2d*-g++`.
- The 3D build command is `./configure --comp=g++ --dim 3 --cuda <arch> --profile`,
  producing `bin/alamo_gpu-3d-profile-cuda{80,90}-g++` (VERIFIED locally for sm_86).
- The 3D production input is `input_3d_flame` (created by task 001; reference by name).
- Multi-GPU model: 1 rank per GPU (roadmap 3.4). Use `--ntasks` = #GPUs, `--gpus-per-task=1`.

## Files to read first
- `benchmark/build_alamo_nova.sh`, `benchmark/nova_flame_gpu.slurm`, `benchmark/nova_flame_cpu.slurm`
- `docs/agent_plans/20260620-gpu-phase3-regime-scaling/PLAN.md`

## Files allowed to modify (create only)
- `benchmark/build_alamo_nova_3d.sh`        — 3D build driver (mirror build_alamo_nova.sh, --dim 3)
- `benchmark/nova_flame_gpu_3d.slurm`       — single-GPU 3D run
- `benchmark/nova_flame_gpu_3d_multi.slurm` — multi-GPU 3D run (1 rank/GPU; NGPUS env, default 4)
- `benchmark/nova_flame_cpu_3d.slurm`       — CPU baseline 3D run (for the GPU-vs-CPU-node crossover)

## Files NOT allowed to modify
- The existing 2D `benchmark/*nova*` / `*.slurm` files (create siblings, don't edit).
- `src/`, `bin/`, repo-root inputs, files owned by tasks 001/003.

## Implementation steps
1. Copy `build_alamo_nova.sh` → `build_alamo_nova_3d.sh`; add `--dim 3` to the configure
   line; keep `ARCHES="80 90"`; update expected binary names to `bin/alamo_gpu-3d-profile-cuda{80,90}-g++`
   and any echo/log text to say 3D. Do not change the login-node/AMReX-priming logic.
2. Copy `nova_flame_gpu.slurm` → `nova_flame_gpu_3d.slurm`; change the binary glob to
   `bin/alamo_gpu-3d*cuda${ARCH}*-g++`; set the input to `${INPUT:-input_3d_flame}`; keep
   `srun --gpus-per-task=1`, the module loads, and the `GPU_TYPE=a100|h200` arch selection.
3. Create `nova_flame_gpu_3d_multi.slurm` from the single-GPU one but parameterized for N GPUs:
   `NGPUS=${NGPUS:-4}`, `#SBATCH --gpus-per-node` / `--ntasks=${NGPUS}`, `srun --ntasks=${NGPUS}
   --gpus-per-task=1`. Add a comment that this is for the roadmap-3.4 strong/weak scaling study.
4. Copy `nova_flame_cpu.slurm` → `nova_flame_cpu_3d.slurm`; input `input_3d_flame`; if it
   builds a CPU binary on the fly, add `--dim 3`; binary pattern `bin/alamo-3d*-g++`.
5. Make the `.sh` executable (chmod +x). Keep all SBATCH partition/account placeholders
   identical to the 2D versions (do not invent cluster-specific values).

## Invariants (device, concurrency, architecture)
- 1 rank per GPU. Do NOT change AMReX/MPI assumptions from the 2D scripts.
- Reference inputs by name (`input_3d_flame`); do not create input files (task 001 owns those).
- Binary selection by glob, never hardcode a single arch's filename.
- Preserve the existing scripts' module loads, srun flags, and env-var conventions verbatim
  except for the documented 3D changes.

## Build and test commands
- No build, no SLURM submission. Validate with `bash -n <script>` (syntax) on every `.sh`/`.slurm`
  and confirm `grep -n "dim 3" benchmark/build_alamo_nova_3d.sh` and
  `grep -n "input_3d_flame" benchmark/nova_flame_gpu_3d*.slurm` both hit.

## Expected result
Four NOVA-ready 3D scripts that mirror the proven 2D workflow, parameterized for A100/H200
and for single- and multi-GPU runs.

## Non-goals
- No actual NOVA submission. No input authoring (task 001). No sweep driver (task 003).
- Do not retune SBATCH resources beyond GPU count for the multi-GPU script.

## Stop conditions
- If any of the three 2D template scripts is missing, STOP and record which.

## Final report: write results/002-RESULT.md (RESULT template).

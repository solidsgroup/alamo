# 002-RESULT: NOVA 3D build + single/multi-GPU launch scripts (roadmap 3.4)

## Summary
Created four NOVA-ready 3D scripts as siblings of the proven 2D workflow, parameterized
for A100 (sm_80) / H200 (sm_90) and for single- and multi-GPU 3D Flame runs. Each mirrors
its 2D template verbatim except for the documented 3D changes (`--dim 3`, 3D binary glob
`*-3d-*`, `input_3d_flame` as the default input, and a 1-rank-per-GPU multi-GPU variant).
No existing 2D script was modified; no build, no GPU, no SLURM submission was performed.
Validated with `bash -n` on all four files plus the two task-named grep checks (both hit).

## Files changed
Created only (none edited):
- `benchmark/build_alamo_nova_3d.sh`        (chmod +x) — 3D build driver; mirrors
  `build_alamo_nova.sh`, configure line + build-job loop use `--dim 3`, keeps
  `ARCHES="80 90"`, expected binaries `bin/alamo_gpu-3d-profile-cuda{80,90}-g++`,
  sbatch file `build_alamo_gpu_3d.sbatch`, 3D-suffixed job name/log/echo text. Login-node
  AMReX-priming logic, modules, account/partition/email placeholders unchanged.
- `benchmark/nova_flame_gpu_3d.slurm`       — single-GPU 3D run; binary glob
  `bin/alamo_gpu-3d*cuda${ARCH}*-g++`, `INPUT=${INPUT:-input_3d_flame}`,
  `srun --gpus-per-task=1`, GPU_TYPE=a100|h200 arch selection, module loads + MODE
  (bench/fast) overrides preserved verbatim.
- `benchmark/nova_flame_gpu_3d_multi.slurm` — multi-GPU 3D run (1 rank/GPU);
  `NGPUS=${NGPUS:-4}`, `NRANKS=${SLURM_NTASKS:-${NGPUS}}`,
  `srun --ntasks="${NRANKS}" --gpus-per-task=1`; comment noting roadmap-3.4
  strong/weak scaling study. Built from the single-GPU script, only GPU-count
  parameterization added (no SBATCH resource retuning beyond GPU count).
- `benchmark/nova_flame_cpu_3d.slurm`       — CPU baseline 3D run; in-job build uses
  `./configure --comp=g++ --dim 3 --profile`, binary pattern `bin/alamo-3d*-g++`,
  `INPUT=${INPUT:-input_3d_flame}`, `--ntasks=16` and modules unchanged.

## Tests run and results
- `bash -n benchmark/build_alamo_nova_3d.sh`        → OK
- `bash -n benchmark/nova_flame_gpu_3d.slurm`       → OK
- `bash -n benchmark/nova_flame_gpu_3d_multi.slurm` → OK
- `bash -n benchmark/nova_flame_cpu_3d.slurm`       → OK
- `grep -n "dim 3" benchmark/build_alamo_nova_3d.sh`            → hits (lines 16, 78, 106)
- `grep -n "input_3d_flame" benchmark/nova_flame_gpu_3d*.slurm` → hits (single + multi)
- `ls -l` confirms `build_alamo_nova_3d.sh` is executable (`-rwxrwxr-x`); the three
  `.slurm` files are non-executable (matches 2D `.slurm` convention; SLURM does not
  require the exec bit and the 2D `.slurm` templates are likewise non-exec).
- `git status --porcelain benchmark/` shows only the four new files as untracked (`??`);
  no existing 2D `*nova*`/`*.slurm` file modified.

## Issues found
- None. All three 2D template scripts were present (no stop condition triggered).

## Deviations from task
- None functionally. Minor mechanical choices within the task's intent:
  - Multi-GPU `srun` uses `--ntasks="${NRANKS}"` where `NRANKS=${SLURM_NTASKS:-${NGPUS}}`,
    so it honors the SLURM-allocated task count when present and falls back to `NGPUS`
    otherwise (the 2D single-GPU script reads `NRANKS=${SLURM_NTASKS:-1}` similarly). This
    satisfies "1 rank per GPU, `--ntasks=${NGPUS}` --gpus-per-task=1" while staying robust
    if the submitter sets `--ntasks` on the sbatch command line (as the header comments show).
  - Generated sbatch artifact named `build_alamo_gpu_3d.sbatch` (2D used
    `build_alamo_gpu.sbatch`) to avoid clobbering the 2D driver's generated file.

## Follow-up needed
- `input_3d_flame` is referenced by name but owned by task 001; these scripts will only
  run end-to-end once that input exists on `chamber-gpu`.
- On NOVA, the module names (`cuda`/`gcc`/`openmpi`) and SBATCH account/partition
  placeholders are inherited verbatim from the 2D scripts and may need site-specific
  editing — unchanged here per the task's "do not invent cluster-specific values" rule.
- Actual NOVA submission, scaling sweeps (task 003), and the 3D build/smoke are out of
  scope (lead owns the build; task 003 owns the sweep driver).

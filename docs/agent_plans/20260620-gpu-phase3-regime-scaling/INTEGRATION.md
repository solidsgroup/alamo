# Integration: GPU Phase 3 - Regime Scaling

Lead: main thread. Workers: 3 file-disjoint Phase 3 tasks in the shared `chamber-gpu`
tree. Base branch: `chamber-gpu`.

## Completed artifacts

- **3.1 readiness report**: `benchmark/PHASE3_3D_READINESS.md`
  - Confirms the 3D CUDA build runs end-to-end on the local A1000.
  - Records the 3D smoke result and the register-count note from the build log.
  - Reaffirms D1: elastic remains CPU-resident; GPU runs use `elastic.type = disable`.

- **Task 001 - 3D inputs**
  - `input_3d_flame`
  - `input_3d_flame_128`
  - `input_3d_flame_256`
  - `input_3d_flame_512`
  - All are wide-shallow 3D inputs with analytic ICs, six-face BCs, and elastic disabled.

- **Task 002 - NOVA 3D scripts**
  - `benchmark/build_alamo_nova_3d.sh`
  - `benchmark/nova_flame_gpu_3d.slurm`
  - `benchmark/nova_flame_gpu_3d_multi.slurm`
  - `benchmark/nova_flame_cpu_3d.slurm`
  - The scripts keep the 2D workflow shape, switch to `input_3d_flame`, and target
    the 3D binary names and `--dim 3` build path.

- **Task 003 - planning tools**
  - `benchmark/phase3_memory_budget.py`
  - `benchmark/phase3_scaling_sweep.sh`
  - `benchmark/PHASE3_R3_crossover.md`
  - These are GPU-free planning artifacts only; the sweep driver is dry-run by default.

## Validation

- `bash -n benchmark/build_alamo_nova_3d.sh`
- `bash -n benchmark/nova_flame_gpu_3d.slurm`
- `bash -n benchmark/nova_flame_gpu_3d_multi.slurm`
- `bash -n benchmark/nova_flame_cpu_3d.slurm`
- `python3 benchmark/phase3_memory_budget.py`
- `bash benchmark/phase3_scaling_sweep.sh`

The memory-budget tool prints a device table. The sweep driver writes `commands.txt`
and does not submit anything unless `--submit` is passed on a SLURM host.

## Merge order

1. Task 001 inputs.
2. Task 002 NOVA scripts.
3. Task 003 planning tools.
4. `benchmark/PHASE3_3D_READINESS.md` and the plan docs.

All three worker tasks are additive and write disjoint file sets. No shared source
files were modified.

## Remaining execution step

Run the NOVA sweep on the target hardware:

- `bash benchmark/phase3_scaling_sweep.sh --submit`

after the task 001 inputs and task 002 scripts are available on NOVA.

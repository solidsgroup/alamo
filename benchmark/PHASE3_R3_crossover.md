# R3: Phase-3 GPU/CPU Crossover Report (SKELETON)

Status: SKELETON -- numeric cells are TODO, to be filled on NOVA. This file is
staged locally (no GPU here); the A1000's 8 GB cannot reach the saturating regime
(roadmap 3.2), so the real sweep runs on NOVA (A100 / H100 / H200).

Roadmap: Phase 3, sections 3.2 (memory budgeting) and 3.5 (crossover hunt).
Branch: `chamber-gpu`. Elastic is DISABLED on the GPU path (decision D1) via
`elastic.type = disable`.

## Purpose

Find where (if anywhere) the 3D GPU build beats a CPU node, sweeping
**problem size x GPU count**. The strategic question from Phase 2 is whether the
GPU win materializes once the problem reaches a regime that saturates the device
(wide-shallow 3D grids, shallow AMR), and whether multi-GPU scaling extends it.

## Method

- Build: NOVA 3D GPU (`benchmark/build_alamo_nova_3d.sh`, sm_80 A100 / sm_90 H200)
  and the CPU-node baseline. (Scripts owned by task 002; referenced by name.)
- Inputs: `input_3d_flame_{128,256,512}` (wide-shallow, max_level <= 1,
  blocking_factor = 32, large max_grid_size). (Owned by task 001.)
- Memory budget: `python3 benchmark/phase3_memory_budget.py` gives the largest
  grid that fits per device (sanity-check each size against the target GPU before
  launch).
- Sweep driver: `bash benchmark/phase3_scaling_sweep.sh` (MODE=strong|weak)
  emits the `sbatch` matrix; run with `--submit` on NOVA.
  - Strong scaling: fixed size, vary GPUS (1, 2, 4, 8).
  - Weak scaling: size grows with GPUS.
- SLURM scripts: `nova_flame_gpu_3d.slurm` (1 GPU), `nova_flame_gpu_3d_multi.slurm`
  (>1 GPU), `nova_flame_cpu_3d.slurm` (CPU node).

### Standing metric set (capture identically every run, per roadmap)

- launches/step
- kernel-duration avg
- `cudaStreamSynchronize` fraction
- achieved occupancy + registers/thread (ncu; needs perf-counter access on NOVA)
- wall/step GPU vs CPU node
- golden-compare residual (correctness vs CPU reference)

## Results -- problem size x hardware

Fill on NOVA. One row per (size, device, #GPUs). Wall/step in seconds; ratio is
GPU/CPU-node (lower is better for GPU). CPU wall/step is the same-size CPU-node run.

| size | device | #GPUs | wall/step GPU | wall/step CPU-node | GPU/CPU ratio | launches/step | sync frac | occupancy |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 128 | A100-80 | 1 | TODO | TODO | TODO | TODO | TODO | TODO |
| 128 | A100-80 | 2 | TODO | TODO | TODO | TODO | TODO | TODO |
| 128 | A100-80 | 4 | TODO | TODO | TODO | TODO | TODO | TODO |
| 256 | A100-80 | 1 | TODO | TODO | TODO | TODO | TODO | TODO |
| 256 | A100-80 | 2 | TODO | TODO | TODO | TODO | TODO | TODO |
| 256 | A100-80 | 4 | TODO | TODO | TODO | TODO | TODO | TODO |
| 512 | A100-80 | 1 | TODO | TODO | TODO | TODO | TODO | TODO |
| 512 | A100-80 | 2 | TODO | TODO | TODO | TODO | TODO | TODO |
| 512 | A100-80 | 4 | TODO | TODO | TODO | TODO | TODO | TODO |
| 512 | A100-80 | 8 | TODO | TODO | TODO | TODO | TODO | TODO |
| 512 | H200 | 1 | TODO | TODO | TODO | TODO | TODO | TODO |

(Add H100 / H200 rows as runs complete. Note the GPU_TYPE used per row.)

## Memory-budget reference (from phase3_memory_budget.py)

Paste the calculator output here for the record (elastic-disabled default, then
the elastic-enabled comparison). Approximate / planning figures.

```
TODO: paste `python3 benchmark/phase3_memory_budget.py` output
TODO: paste `python3 benchmark/phase3_memory_budget.py --bytes-per-node 512` output
```

## Scaling analysis

- Strong scaling efficiency per size: TODO (speedup vs #GPUs; ideal = linear).
- Weak scaling efficiency: TODO (wall/step should stay flat as size grows with GPUS).
- Where does GPU/CPU ratio cross 1.0 (GPU wins)? size=TODO, #GPUs=TODO, device=TODO.

## Crossover / verdict (decision D3)

Pick exactly one of the three valid outcomes and justify with the table above:

- [ ] **WIN @ single** -- the GPU build beats a CPU node at single-GPU once the
      problem reaches the saturating regime. Smallest winning (size, device): TODO.
- [ ] **WIN @ scale** -- no single-GPU win, but multi-GPU scaling beats the CPU
      node at #GPUs >= TODO for size >= TODO on device TODO.
- [ ] **NO WIN** -- the CPU node beats the GPU build across all tested
      (size x #GPUs); GPU regime scaling does not pay off. Recommend stopping
      Phase-3 perf work and recording the negative result.

Verdict: TODO

Supporting notes / caveats: TODO (e.g. occupancy ceiling, sync-fraction
bottleneck, golden-residual drift, OOM at largest size).

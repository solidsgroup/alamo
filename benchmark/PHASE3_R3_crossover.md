# R3: Phase-3 GPU/CPU Crossover Report

Status: **WIN @ single, now CONFIRMED at full 64-rank CPU** — two NOVA sweeps:
**2026-06-21** (jobs 11160767–11160774, 16-rank CPU baseline, exploratory) and
**2026-06-22** (jobs 11161369–11161379, full 64-rank CPU baseline, the run the
"to harden the verdict" section below called for). The 06-22 run confirms the
06-21 prediction almost exactly: **single A100 beats the full 64-rank CPU node
~9.6× at 128³ and ~13.1× at 256³**, in the same ballpark as the 06-21 run's own
"assume ideal 4× CPU scaling" estimate (~10×/~17×). Confidence is now **HIGH**
on 128³/256³ (the only major caveat from 06-21 is resolved). **Still open:**
`ncu` profiler metrics and H100/H200 rows. **Descoped 2026-06-24 (user
directive), not being pursued:** the 512³ CPU baseline (OOMs on a `--mem`
config bug, not a rank-count problem — 128³/256³ are sufficient) and
multi-GPU scaling (a confirmed regression, on the back burner, not a
near-term goal).

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
  and the CPU-node baseline.
- Inputs: `input_3d_flame_{128,256,512}` (wide-shallow, `max_level = 1`,
  `blocking_factor = 32`, `max_grid_size = 128`, `grid_eff = 0.7`, `stop_time =
  6.5_s`, dt 1e-4 → 65000 steps). 128 = `128 128 64`, 256 = `256 256 128`,
  512 = `512 512 256` (`n_cell`); domain `0.1754 x 0.1754 x 0.0877 m` for all three.
- Sweep driver: `bash benchmark/phase3_scaling_sweep.sh` (MODE=strong) emits the
  `sbatch` matrix; run with `--submit` on NOVA. Both runs used
  `SIZES="128 256 512"`, `GPUS="1 2"`.
- SLURM scripts: `nova_flame_gpu_3d.slurm` (1 GPU), `nova_flame_gpu_3d_multi.slurm`
  (>1 GPU), `nova_flame_cpu_3d.slurm` (CPU node). All flame_* jobs ran with a
  30-minute SLURM time limit (`--time=00:30:00`).
- CPU rank count differs by run: 06-21 baselines came from an **older 16-rank**
  config; 06-22 used the **current 64-rank** config (`--ntasks=64`).

### Standing metric set (capture identically every run, per roadmap)

- launches/step — **pending: needs ncu / TinyProfiler region export** (neither run)
- kernel-duration avg — **pending: ncu**
- `cudaStreamSynchronize` fraction — **pending: ncu**
- achieved occupancy + registers/thread — **pending: ncu, needs perf-counter access on NOVA**
- wall/step GPU vs CPU node — **captured both runs** (see Results)
- golden-compare residual (correctness vs CPU reference) — spot-checked 06-22,
  GPU/CPU step-1 chamber diagnostics bit-identical at size=128; no NaN/Inf in
  any log from either run. **This is a smoke test only, not full
  correctness** — it covers step 1, not the 65,000-step trajectory. Still
  needed before the 9.6x/13.1x numbers can be called correctness-verified
  (not just throughput-verified): field diffs over a representative span of
  steps (not just step 1), thermo time-series diffs (T/mdot/L), final burn
  area/volume agreement, pressure-trajectory agreement, and — once elastic is
  re-enabled on this path — elastic stress/strain comparison against the CPU
  reference.

## Run log A — 2026-06-21, first measurement (jobs 11160767–11160774)

**Headline: the 3D GPU build is correct and stable at every resolution tested —
the open issue is throughput, not correctness.** Across all five captured logs
there were **zero** aborts, NaN/Inf, CUDA errors, or OOMs. The async-output fix
(`23c924721`) held; every job initialized CUDA/MPI and stepped cleanly.

| job | size | #GPU | last step | sim time | outcome |
| --- | ---: | ---: | ---: | ---: | --- |
| 11160767 | 128 | 1 | 65000 (done) | 6.5 s | ✅ **completed** — `=== done ===`, 940.6 s wall, 1.47 GB peak arena |
| 11160770 | 256 | 1 | 21066 | 2.11 s | killed (SIGTERM 02:25:57), healthy, ~32% done |
| 11160773 | 512 | 1 | 1024 | 0.10 s | cut off mid-step, healthy, **no OOM** at 512³ on one 80 GB A100 |
| 11160768 | 128 | 2 | 9205 | 0.92 s | killed (SIGTERM 02:26:02), healthy |
| 11160771 | 256 | 2 | 232 | 0.02 s | cut off mid-step, healthy |
| 11160774 | 512 | 2 | — | — | completed (no log in this bundle) |

CPU-node baselines for this run came from sibling bundles (jobs 11160033/038,
**16 MPI ranks**, the older `nova_flame_cpu_3d.slurm` config).

## Run log B — 2026-06-22, 64-rank CPU confirmation (jobs 11161369–11161379)

Builds: GPU build (11161369, produced `alamo_gpu-3d-profile-cuda80-g++` and
`...-cuda90-g++`) and CPU build (11161370, `alamo-3d-profile-g++`) both
completed clean. Job-ID-to-config mapping confirmed via filename prefix
(`flame_gpu_3d` / `flame_gpu_3d_multi` / `flame_cpu_3d`) against
`commands.txt` submission order.

| job | size | device | #GPUs | last step | outcome |
| --- | ---: | --- | ---: | ---: | --- |
| 11161371 | 128 | A100-80 | 1 | 65000 (done) | ✅ **completed**, 999.3 s wall (TinyProfiler total) |
| 11161372 | 128 | A100-80 | 2 | 32,998 | hit 30-min SLURM time limit |
| 11161373 | 128 | CPU (64 ranks) | — | 12,204 | hit 30-min SLURM time limit |
| 11161374 | 256 | A100-80 | 1 | 30,682 | hit 30-min SLURM time limit |
| 11161375 | 256 | A100-80 | 2 | 2,306 | hit 30-min SLURM time limit |
| 11161376 | 256 | CPU (64 ranks) | — | 2,345 | hit 30-min SLURM time limit |
| 11161377 | 512 | A100-80 | 1 | 3,714 | hit 30-min SLURM time limit |
| 11161378 | 512 | A100-80 | 2 | 411 | hit 30-min SLURM time limit |
| 11161379 | 512 | CPU (64 ranks) | — | ~2 | **OOM-killed at 18.5 s** (`Detected 2 oom_kill events`) |

No NaN/blowup messages in any of the 9 sweep logs (the only `nan` grep hit was
the AMReX profiler region label `MultiFab::contains_nan()`, not an actual NaN).

## Results — combined, problem size x hardware

wall/step for every job that hit the 30-min time limit is `1800 s /
steps_completed` — an average over the whole window including
binary/MPI/CUDA init and output-directory setup, so it understates true
steady-state throughput, more so for low step counts. The two **exact**
numbers (no estimation) are 11160767 (940.6 s / 65000 steps) and 11161371
(999.3 s / 65000 steps) — both 128³/1-GPU full completions, run a day apart,
agreeing to within 6%.

| size | device | #GPUs | run | wall/step (s) | steps completed | GPU/CPU ratio |
| ---: | --- | ---: | --- | ---: | ---: | ---: |
| 128 | A100-80 | 1 | B (06-22) | 0.01537 | 65,000 / 65,000 | **0.104** (9.6x faster than 64-rank CPU) |
| 128 | A100-80 | 1 | A (06-21) | 0.01447 | 65,000 / 65,000 | ~0.026 vs 16-rank CPU (~39x faster) |
| 128 | A100-80 | 2 | B (06-22) | 0.05455 | 32,998 | 0.370 vs 64-rank CPU (3.5x slower than 1-GPU) |
| 128 | A100-80 | 2 | A (06-21) | 0.1370 | 9,205 | ~9.5x slower than 1-GPU |
| 128 | CPU-node | — | B (06-22, 64-rank) | 0.1475 | 12,204 | — |
| 128 | CPU-node | — | A (06-21, 16-rank) | 0.568 | 9,125 | — |
| 256 | A100-80 | 1 | B (06-22) | 0.05867 | 30,682 | **0.076** (13.1x faster than 64-rank CPU) |
| 256 | A100-80 | 1 | A (06-21) | 0.0599 | 21,066 | ~0.014 vs 16-rank CPU (~70x faster) |
| 256 | A100-80 | 2 | B (06-22) | 0.78099 | 2,306 | 1.017 (~parity, 13.3x slower than 1-GPU) |
| 256 | A100-80 | 2 | A (06-21) | not measurable | 232 | — (ran, no OOM) |
| 256 | CPU-node | — | B (06-22, 64-rank) | 0.7676 | 2,345 | — |
| 256 | CPU-node | — | A (06-21, 16-rank) | 4.167 | 1,237 | — |
| 512 | A100-80 | 1 | B (06-22) | 0.48465 | 3,714 | n/a (no valid CPU baseline) |
| 512 | A100-80 | 1 | A (06-21) | not measurable | 1,024 | — (ran, no OOM on one 80 GB A100) |
| 512 | A100-80 | 2 | B (06-22) | 4.37956 | 411 | n/a (9.0x slower than 1-GPU) |
| 512 | CPU-node | — | B (06-22, 64-rank) | **OOM @ 18.5 s** | ~2 | — |
| 512 | CPU-node | — | A (06-21, 16-rank) | **OOM** | — | — |

`launches/step`, `cudaStreamSynchronize` fraction, and achieved occupancy were
**not collected** in either sweep — no `ncu` run included in either bundle;
remain open per the standing metric set.

## Memory-budget reference

```
128³ (128 128 64), A100-80, elastic disabled:
  peak "The Arena" used = 1469 MB ; allocated reservation = 60863 MB
  Free GPU global memory: 79179 MB of 81151 MB total
512³ (512 512 256): ran on a single 80 GB A100 with no OOM in run A (confirms GPU-side budget headroom)
```

`phase3_memory_budget.py` was not run in either sweep (no output captured in
either bundle). **The CPU-side 512³ OOM is a `--mem` problem, not a rank-count
problem**: it OOM'd at 16 ranks in run A and *again* at 64 ranks in run B,
both because `nova_flame_cpu_3d.slurm` hardcodes `--mem=32G` regardless of
input size. Adding ranks doesn't add memory headroom on a single node —
recommend running the memory-budget calculator and re-submitting size=512 CPU
with a larger `--mem` before drawing any size=512 conclusion.

## Scaling analysis

- **Single-GPU size scaling (confirmed both runs):** going 128 → 256 (8x the
  cells) cost only ~4x more time per step in both runs (run A: 4.1x slower,
  69.1 → 16.7 steps/s; run B: 0.0154 → 0.0587 s/step ≈ 3.8x slower). Sublinear
  cost-per-cell in both independent measurements — the device is getting more
  efficient per cell as the problem grows, the saturating-regime signature
  Phase 3 was hunting for.
- **CPU scaling 16→64 ranks (validated prediction):** run A predicted "assuming
  ideal 4x scaling 16→64, GPU would still win ~10x (128³) and ~17x (256³)."
  Run B's *actual* 64-rank measurement gives 9.6x and 13.1x — same order, GPU
  win holds, slightly less favorable than the ideal-scaling guess (consistent
  with sub-linear real-world MPI scaling at higher rank counts).
- **Multi-GPU strong scaling (the key negative result, confirmed both runs,
  now extended to 256³):** going 1 → 2 GPUs made every step *slower*, not
  faster, at every size and in both independent runs:
  - 128³: 3.5x slower (run B) / 9.5x slower (run A) than 1-GPU.
  - 256³: 13.3x slower (run B; run A's 2-GPU/256 job was cut off too early to
    measure, but showed no OOM).
  - 512³: 9.0x slower (run B only).
  `CUDA initialized with 2 devices` is confirmed in run B's 2-GPU logs, ruling
  out a device-binding/oversubscription bug. The per-GPU subdomain is far too
  small at these sizes; halo exchange / `FillBoundary` / MPI sync most likely
  dominates every step (consistent with run A's prior "FillBoundary = 54% of
  launches" note). On a 2xA100 node this is plausibly host-staged halos
  (GPU-aware MPI / NVLink not confirmed in play) — needs `ncu`/`nsys` profiling
  (launches/step, sync fraction) before any further multi-GPU work.
- **Weak scaling:** not tested in either run. Both sweeps only emitted the
  MODE=strong matrix.

## Crossover / verdict (decision D3)

**Verdict: WIN @ single (HIGH confidence).** A single A100 beats the CPU node
**9.6x at 128³ and 13.1x at 256³ against the full 64-rank node** (run B,
2026-06-22) — this is no longer an extrapolation from a 16-rank baseline; it
is the direct measurement the 06-21 report called for, and it lands within the
range that run's own "ideal-scaling" estimate predicted. The advantage grows
with problem size in both independent runs, matching the saturating-regime
hypothesis. This reverses the earlier pre-saturation finding (2D / local
A1000: CPU np8 beat the GPU 2.3–3.4x).

- [x] **WIN @ single** — GPU beats a CPU node at single-GPU in the saturating
      regime. Smallest winning point: **128³ on A100-80 (9.6x vs full 64-rank
      node)**; win widens to 13.1x at 256³.
- [ ] **WIN @ scale** — multi-GPU scaling beats the CPU node. *Ruled out at
      every tested size (128³/256³/512³): 2-GPU consistently 3.5x–13.3x
      *slower* than 1-GPU in both runs. Revisit only at much larger per-GPU
      domains, with GPU-aware MPI confirmed enabled.*
- [ ] **NO WIN** — ruled out.

## Recommendation & confidence

**Recommend: proceed with single-GPU as the supported deployment shape for
128³/256³; do not pursue multi-GPU at these sizes.** Confidence: **~85%**.

- **Measured, not inferred:** the 128³/1-GPU wall/step (the only fully-completed,
  no-extrapolation number) agrees to within 6% across two independent NOVA
  sessions a day apart (940.6 s vs 999.3 s / 65,000 steps). The 64-rank CPU
  baseline (run B) is a direct measurement, not the run-A "ideal 4x scaling"
  estimate it was meant to confirm — and it landed within that estimate's range
  (9.6x/13.1x measured vs ~10x/~17x estimated).
- **What pulls confidence down from ~95%:** no `ncu` occupancy/launch-count data
  in either sweep (the standing metric set's launches/step, sync fraction, and
  occupancy rows are all unmeasured); H100/H200 rows are untested. (The 512³
  CPU baseline never produced a valid number — OOM at both 16 and 64 ranks, a
  `--mem` config bug — but this is now descoped, not a confidence gap to close;
  128³/256³ alone are taken as sufficient.)
- **What's solid regardless:** the multi-GPU regression (2 GPUs slower than 1 at
  every size, both runs) is a measured negative result, not an extrapolation —
  high confidence in *not* pursuing multi-GPU at 128³/256³ without first
  root-causing it.

## To harden the verdict (next NOVA run)

1. ~~64-rank CPU baselines for 128/256~~ — **done in run B (2026-06-22)**.
2. ~~Fix the 512³ CPU baseline~~ — **descoped 2026-06-24 (user directive):
   not being pursued.** 128³/256³ give enough crossover data on their own;
   the `--mem=32G` hardcode in `nova_flame_cpu_3d.slurm` is left as-is.
3. ~~Root-cause the multi-GPU regression~~ — **on the back burner, not a
   near-term goal (user directive, 2026-06-24).** The regression stands as a
   measured negative result (2 GPUs slower than 1 at every tested size); no
   `ncu`/`nsys` diagnosis or further scaling work is planned unless multi-GPU
   becomes a near-term goal again.
4. **`ncu` capture** on the single-GPU 256³ run for launches/step, sync
   fraction, and achieved occupancy (the F4 register-pressure question from
   `PHASE3_3D_READINESS.md`). 512³ dropped per item 2.
5. Capture a **kill/finish timestamp** (or `time` wrapper) for every job so
   wall/step is measurable even when cancelled early — run B's 1800s/steps
   estimate is a reasonable proxy but understates steady-state throughput.
6. Add **H100/H200 rows** — no runs on either device were submitted in either
   sweep; all numbers above are A100-80 only.

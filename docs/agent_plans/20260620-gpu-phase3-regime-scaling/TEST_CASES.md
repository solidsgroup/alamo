# What and why: input_3d_flame_{128,256,512}

These three input decks (plus the canonical `input_3d_flame`, which duplicates
`_256`) are the Phase-3 crossover-hunt sweep — the problem-size axis of
roadmap §3.5 ("sweep problem size and GPU count to locate the point where the
GPU build beats the CPU node, and by how much"). Authored by task 001 of this
folder (`tasks/001-3d-inputs.md` / `results/001-RESULT.md`); this doc is the
standalone what/why summary that sits above that task-level detail.

## What each one is

All three share identical physics (propellant, chamber, thermal blocks copied
from `input`), identical box strategy (wide-shallow: `amr.max_level = 1`,
`amr.blocking_factor = 32`, `amr.max_grid_size = 128`, analytic Expression ICs
only — `IC::BMP` is 2D-only), and identical elastic disposition
(`elastic.type = disable`, no elastic model/bc block — D1 keeps elastic
CPU-resident, so Phase 3 measures the GPU-resident phase-field + thermal path
only). The **only** difference between them is `amr.n_cell`:

| Variant | `amr.n_cell` | Base cells | Role |
|---|---|---:|---|
| `input_3d_flame_128` | 128 128 64  | ~1.0M  | smallest sweep point — sanity/fast-iteration size, fits on the local A1000 |
| `input_3d_flame_256` (= canonical `input_3d_flame`) | 256 256 128 | ~8.4M | mid sweep point |
| `input_3d_flame_512` | 512 512 256 | ~67M | largest sweep point — intended to reach the saturating regime on high-memory NOVA devices (A100 40/80G, H200 141G); exceeds the local A1000's 8G |

Domain is fixed across all three (`geometry.prob_hi = 0.1754m 0.1754m 0.0877m`,
wide in x-y, thin in z) — so the three sizes are a genuine resolution sweep at
fixed physical extent, not just a box-count change.

## Why three sizes (not one)

This is the empirical arm of decision tree **D3** in `docs/llm/ROADMAP.md`
(roadmap §3.5): "does a saturating 3D problem fit in target-GPU memory, and at
the largest fitting size, do phase-field kernels saturate?" A single size
can't answer that — you need at least one point that's clearly
launch-latency-bound (128, cheap to iterate, runs anywhere) and one that's
large enough to plausibly saturate device occupancy (512, NOVA-only). 256 sits
in between as the canonical/production size and a sanity check that scaling
behaves monotonically. The same three-point logic mirrors the earlier 2D
coarse/deep/wide campaigns (`docs/llm/perf/2026-06-20-gpu-port-report.md`
§5-6): bigger, fewer-but-larger boxes is the lever that moves the GPU from
launch-bound to compute-bound, and you need multiple sizes to see the trend,
not just one number.

## What's actually come out of them so far

- **2026-06-20 NOVA attempt** (`analysis/phase3_nova_bundle_20260620_212331/`):
  **failed** — every captured 1-GPU and 2-GPU A100 job for all three sizes
  aborted with CUDA error 700 (illegal memory access) before completing a
  timestep; the 512 CPU-node baseline was OOM-killed. No crossover data. See
  `docs/agent_plans/20260622-analysis-suite-backfill/tasks/002-results-bundles.md`
  for the full breakdown — this bundle is not the source of any of the numbers
  below.
- **2026-06-21 NOVA sweep** (jobs 11160767-11160774, A100-80 sm_80, after
  whatever fixed the above failure): **first real crossover measured**. Single
  A100 beats the CPU node ~39x at 128^3 and ~70x at 256^3 (D3 = WIN @ single).
  Full table and methodology: `benchmark/PHASE3_R3_crossover.md`. No local
  `analysis/` bundle exists for this successful run (the numbers are sourced
  from NOVA job logs directly).
- **Still open** (per `benchmark/PHASE3_R3_crossover.md` and
  `docs/llm/CURRENT.md`): a full 64-rank CPU-node confirmation, `ncu`
  occupancy/register metrics, and GPU counts > 2 (the 512 point and
  multi-GPU weak/strong scaling haven't produced usable numbers yet).

# Alpha-1.0 Baseline of Record

**Tag:** `chamber-gpu-alpha1`
**Commit:** `c00f69086c6d1dc67bdf71e2bd174dbdcc953c85`
**Date:** 2026-06-25
**Branch:** `chamber-gpu` (never merged to master — see `GPU_BRANCH_GUIDE.md`)

---

## What this tag represents

This is the frozen reference point for the GPU Roadmap v2 Phase A instrumentation
campaign (`benchmark/GPU_ROADMAP_V2.md`, task A0). It is the commit at which:

- The multi-box elastic MLMG divergence was fully resolved (one-line `elixir()` fix
  in `Operator<Grid::Node>::interpolation()`, `src/Operator/Operator.cpp:728`).
- Branch DoD items 1–5 are all **DONE** (see `benchmark/PHASE5_BRANCH_DONE.md`).
- Phase-field GPU speedup is **proven**: single A100 vs 64-rank CPU node ≈ 9.6× at
  128³ and 13.1× at 256³ (`benchmark/PHASE3_R3_crossover.md`).
- Combined flame+elastic A100 performance is **unmeasured** — all Phase-3 crossover
  runs used `elastic.type = disable`.

Phase A2 (`benchmark/PHASE_A2_combined_crossover.md`) is the next measurement against
this baseline, with elastic enabled.

---

## Build matrix

Two builds are maintained (see `benchmark/GPU_BRANCH_GUIDE.md`):

| Purpose | Configure flags | Binary pattern |
|---------|----------------|----------------|
| Performance headlines | `--cuda 80 --profile` | `alamo_gpu-3d-profile-cuda80-g++` |
| Correctness / golden-compare | `--cuda 80 --profile` + no-fast-math flag | strict build |
| Multi-arch NOVA build | `--cuda "70 80 90" --profile` | `alamo_gpu-3d-profile-cuda{70,80,90}-g++` |

Build script: `benchmark/build_alamo_nova_3d.sh`

Key configure options (3D GPU build):
```
./configure --comp=g++ --dim 3 --cuda 80 --profile --get-eigen
```
All three NOVA targets (V100/A100/H200):
```
./configure --comp=g++ --dim 3 --cuda "70 80 90" --profile --get-eigen
```
Fast-math is enabled by default in `--profile` builds (`--use_fast_math`).

---

## Canonical input set (alpha-1.0 reference configs)

These are the inputs used for Phase-3 crossover measurements (elastic **disabled**).
A2 will introduce elastic-enabled variants.

| Input file | Grid | Elastic | Use |
|------------|------|---------|-----|
| `input_3d_flame_128` | 128³ | `type=disable` | alpha-1.0 Phase-3 reference |
| `input_3d_flame_256` | 256³ | `type=disable` | alpha-1.0 Phase-3 reference |
| `input_3d_flame_512` | 512³ | `type=disable` | alpha-1.0 Phase-3 reference (CPU OOM) |
| `input_3d_centre_bore_128` | 128³ | `type=disable` | centre-bore geometry reference |
| `input_3d_centre_bore` | 256³ | `type=disable` | centre-bore geometry reference |
| `input_3d_centre_bore_512` | 512³ | `type=disable` | centre-bore geometry reference |

For A2, new elastic-enabled variants (`input_3d_flame_128_elastic`,
`input_3d_flame_256_elastic`) will be added alongside these frozen files.

---

## CPU regression suite status at alpha-1.0

- 118 tests run, 92 verified, **0 failed** on CPU suite.
- CI workflow: `.github/workflows/chamber-gpu-correctness.yml`
- 5 pre-existing CPU-suite failures in Flame lineage (thermal.on=0 null write,
  model_prop query_exactly<2> parse abort) are **not** caused by GPU changes;
  documented in `benchmark/PHASE4_R4_dispatch.md`.

---

## Key prior measurements (do not re-litigate)

| Measurement | Value | Source |
|-------------|-------|--------|
| Phase-field speedup, 128³, 1 A100 vs 64-rank CPU | **9.6×** | `PHASE3_R3_crossover.md` |
| Phase-field speedup, 256³, 1 A100 vs 64-rank CPU | **13.1×** | `PHASE3_R3_crossover.md` |
| Multi-GPU (2 GPU) vs 1 GPU efficiency, 128³ | **5.3%** (loss) | `PHASE3_R3_crossover.md` |
| Combined flame+elastic A100 speedup | **unmeasured** | — |
| Elastic on GPU: correctness | **FIXED** (elixir, commit c00f690) | `GPU_BRANCH_GUIDE.md` |
| Elastic on GPU: performance | **unmeasured** (only 2D local A1000 data) | `PHASE1_ELASTIC_DISPOSITION.md` |

---

## Cross-references

- Forward plan: `benchmark/GPU_ROADMAP_V2.md`
- Branch policy: `benchmark/GPU_BRANCH_GUIDE.md`
- Phase 0–5 DoD (complete): `benchmark/PHASE5_BRANCH_DONE.md`
- Phase-3 crossover results: `benchmark/PHASE3_R3_crossover.md`
- Elastic fix root-cause: `GPU_BRANCH_GUIDE.md` D1 section
- Elixir-race audit (A4): `benchmark/elixir_race_audit.md`
- Static kernel resource dump: `benchmark/G0_BASELINE_OF_RECORD.md`

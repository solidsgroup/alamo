# chamber-gpu Branch Guide

This is the top-level map of the `chamber-gpu` branch: what it is, how to build
it, what is safe to run, and where the supporting evidence lives. It indexes
existing documents rather than duplicating them — follow the links for detail.

Roadmap: `~/Desktop/GPU-OPT-ROADMAP.txt` (Phases 0-5, **complete** — all 5 DoD
items DONE). This guide corresponds to roadmap step **5.3 (Documentation)**; the
companion checklist is `benchmark/PHASE5_BRANCH_DONE.md` (step 5.4).

**Forward plan (alpha-1.0 → beta):** `benchmark/GPU_ROADMAP_V2.md` — the
measurement-driven v2 roadmap (Phases A–D). The headline gap it closes: every
perf number so far is *one solver at a time* (phase-field crossover ran with
elastic disabled; elastic is correctness-fixed but never benchmarked on A100).
Phase A re-baselines combined flame+elastic on A100 and finally unblocks
`ncu`/`nsys` counters on NOVA.

## Branch policy: never merged

`chamber-gpu` is a **permanent, standalone branch**. It is never merged into
`master`. "Done" means the branch satisfies its own definition-of-done
(`PHASE5_BRANCH_DONE.md`), not that it lands on master. The CUDA/GPU work is
isolated behind a separate entry point (`src/alamo_gpu.cc`) and build policy
(`src/GPU/IntegratorPolicy.mk`) so the shared CPU build (`alamo.cc`, all other
integrators) is untouched. See the D4 = ISOLATE decision in
`benchmark/PHASE4_R4_dispatch.md` (Phase 4) for the reasoning.

## Build matrix: which binary for what

| Build | Command | Binary | Use for |
| --- | --- | --- | --- |
| CPU (baseline) | `./configure --comp=g++ --dim=2 && make -j8` | `bin/alamo-2d-g++` | Reference correctness + multi-rank baseline (`np8`) |
| GPU fast (local) | `benchmark/build_alamo_local_gpu.sh` (or `./configure --comp=g++ --dim=2 --cuda local && make`) | `bin/alamo_gpu-2d-cuda86-g++` | Local performance numbers; uses `--use_fast_math` |
| GPU strict / no-fast-math (local) | `./configure --comp=g++ --dim=2 --cuda local --cuda-fp strict && make` | `bin/alamo_gpu-2d-nofast-cuda86-g++` | Correctness gates / golden compare (`--fmad=false`, no fast-math intrinsics) |
| NOVA 3D (A100/H200) | `benchmark/build_alamo_nova_3d.sh` (see also `benchmark/build_alamo_nova.sh` for 2D) | `bin/alamo_gpu-3d-cuda80-*` / `cuda90-*` | Regime-scaling / crossover runs on NOVA |

Full build instructions, NOVA quickstart, async-IO flags, benchmarking scripts
(`benchmark_gpu_cpu.sh`), and the `MODE=fast` vs `MODE=bench` distinction live
in `benchmark/README.md` — read that for the actual commands. Local helper:
`benchmark/build_alamo_local_gpu.sh`; NOVA helpers:
`benchmark/build_alamo_nova.sh`, `benchmark/build_alamo_nova_3d.sh`.

Fast vs strict, in one line: **fast** (`--use_fast_math`) is for performance
headlines; **strict/no-fast-math** (`--cuda-fp strict`) is the only build
correctness claims may be based on. Both binaries coexist (different
`POSTFIX`).

GPU build isolation mechanism (Phase 4.3, already done): `src/GPU/IntegratorPolicy.mk`
declares `ALAMO_GPU_SUPPORTED_INTEGRATORS := flame` and the GPU-clean source
closure; the Makefile's `ifneq (,$(findstring cuda,$(POSTFIX)))` block consumes
it. `src/alamo_gpu.cc` is the Flame-only `main`. Details and the encapsulation
convention for device-captured kernels: `docs/gpu_device_capture_conventions.md`.

## GPU-safe IC/BC matrix

Full matrix: `docs/gpu_safe_ic_bc_matrix.md`.

Summary: only `IC::BMP`, `IC::Constant`, `IC::Expression` (scalar),
`BC::Constant`, and `BC::Operator::Elastic::Constant` are device-safe today.
`IC::PNG`/`PSRead`/`Trig`/`Laminate`, `IC::Expression` (vector), and
`BC::Operator::Elastic::Expression` still use `LoopConcurrentOnCpu` and write
into device-arena fabs from the host; CUDA builds guard these to abort before
reaching the unsafe write (verified by `benchmark/test_gpu_guarded_ic.sh`).
**Canonical chamber GPU inputs must stick to BMP/Constant ICs and Constant
elastic BCs.**

## Elastic disposition (D1): RESOLVED 2026-06-25 — device-elastic works on GPU

Full reports: `benchmark/elastic_sensitivity_20260621/` (root-cause hunt + fix;
see `GPU_ELASTIC_DEBUG_PLAN.md` SOLVED banner and `fix_notes.md`),
`benchmark/PHASE1_ELASTIC_DISPOSITION.md` (original verdict).

**The GPU elastic MLMG divergence is root-caused and FIXED.** Root cause: a GPU
**cross-stream use-after-free race** on the per-box temporary `FArrayBox tmpfab`
in `Operator<Grid::Node>::interpolation()` (`src/Operator/Operator.cpp`). AMReX's
non-tiled GPU `MFIter` cycles iterations across a pool of CUDA streams; `tmpfab`
was freed at iteration end and its device block reused by a later iteration on a
**different stream** while the interpolation/`plus` kernels were still in flight,
corrupting the coarse-grid correction — which the BC-penalty diagonal (~1e17)
then amplified into the 1e20 blow-up. `restriction()` writes its result directly
(no temp), which is exactly why only interpolation broke. **Fix (one line):**
`amrex::Gpu::Elixir tmpfab_eli = tmpfab.elixir();` after `tmpfab.resize(...)` —
defers the temp's free until its own stream completes.

How it was localized (all on the local A1000): (1) box-size sweep — single box
(`max_grid_size=2048`) converges, **any** multi-box decomposition diverged → the
defect is in **cross-box** nodal transfer, not resolution or the bottom solver;
(2) a stock-AMReX `MLNodeLaplacian` reproducer converges multi-box at 2048² → the
defect is in ALAMO's custom transfer overrides, not AMReX/driver/BiCGStab; (3) a
bit-identical CPU↔GPU transfer probe → `interpolation` (not `restriction`); (4) a
pre/post-`nodalSync` bisect → the per-box temp, nondeterministic run-to-run. This
reconciles the whole prior history: the "BiCGStab noise-floor",
"nondeterministic reduction", "deterministic coarse-level operator defect", and
"flaky 2048² resolution threshold" (really a box-count/box-size effect) were all
**downstream of this one race**. The earlier `interpolation()` "live lead" was
the right neighborhood.

Verified end-to-end on the strict/no-fast-math binary: full box sweep
(`max_grid_size` 2048..128 + deep coarsening) **all converge ~1e-9**; the REAL
production operator (`use_psi=1`, 3-material psi-weighted) at multi-box 2048²
**converges in 7 iters to 8.5e-9** (this was the original failing case); an
end-to-end flame+elastic chamber sim (forced multi-box, 300 steps) runs clean —
**10 elastic solves all <1e-8**, no NaN/abort, sane physics (chamber P=5.38 MPa,
mdot=3.65 kg/s at t=0.03s). Single-box behavior unchanged.

Performance footnote (separate, still-open question): where the GPU elastic solve
converged, a fully-subscribed CPU `np8` node was **2.27–3.40× faster** at matched
resolution/iteration count (the original D1 perf data). The *correctness* blocker
is now removed; whether device-elastic is *worth it* on perf grounds is a
distinct question that the original D1 numbers still bear on.

## Crossover (D3): WIN @ single (3D, NOVA)

Full report: `benchmark/PHASE3_R3_crossover.md`.

Summary: on NOVA A100s with the 3D wide-shallow build (`max_level=1`, elastic
disabled), a **single A100 beats a full 64-rank CPU node 9.6× at 128³ and
13.1× at 256³**, with the advantage growing with problem size (sublinear GPU
cost scaling vs near-linear CPU). **D3 verdict: WIN @ single** (HIGH
confidence on 128³/256³ — this is the 64-rank-confirmed figure; an earlier
exploratory 16-rank estimate of ~39×/~70× is superseded and should no longer
be quoted). The 512³ CPU baseline and multi-GPU scaling are **descoped
2026-06-24 (user directive) — not being pursued**: 512³ OOMs on a `--mem`
config bug (128³/256³ are sufficient data); multi-GPU remains a confirmed
negative result (2 GPUs vs 1, 3.5–13.3× slower at every tested size,
non-root-caused) but is on the back burner, not a near-term goal.

Net picture: GPU phase-field wins decisively in the saturating 3D regime. D3 is
independent of D1 — elastic is disabled entirely on the GPU path for the 3D
crossover runs, so this crossover says nothing about whether elastic itself
should run on device (see D1 above, separately reopened by user directive).

## Where the phase reports live

| Phase | Report | File |
| --- | --- | --- |
| 0 | Baseline of record | `benchmark/G0_BASELINE_OF_RECORD.md` |
| 1 | Elastic path disposition (D1) | `benchmark/PHASE1_ELASTIC_DISPOSITION.md` |
| 2 | Phase-field optimization attribution (D2) | `benchmark/phase2_box_sweep/`, `benchmark/phase2_probe_gpu_plot/` (see also `docs/gpu_elastic_device_port_plan.md` for the deprioritized device-elastic plan) |
| 3 | 3D readiness + crossover (D3) | `benchmark/PHASE3_3D_READINESS.md`, `benchmark/PHASE3_NOVA_TESTING_PROGRESS.md`, `benchmark/PHASE3_R3_crossover.md` |
| 4 | Framework dispatch (D4) + CPU regression | `benchmark/PHASE4_R4_dispatch.md`, CPU-regression logs under `benchmark/phase4_cpu_semantics_*/` |
| 5 | Hardening: CI, perf tracking, docs, DoD | `.github/workflows/chamber-gpu-correctness.yml` + `benchmark/ci_golden_compare.sh`; `benchmark/perf_regression_track.py` + `benchmark/PERF_TRACKING.md`; this guide; `benchmark/PHASE5_BRANCH_DONE.md` |

## Correctness tooling quick reference

- Golden CPU-vs-GPU compare on a short run: `benchmark/golden_compare_flame.sh`
- Guarded-IC abort check: `benchmark/test_gpu_guarded_ic.sh`
- Standardized multi-case baseline suite (record/check/report):
  `benchmark/baseline_suite.py`, references in `benchmark/baseline_references/`
- CI wrapper (Phase 5.1): `benchmark/ci_golden_compare.sh` +
  `.github/workflows/chamber-gpu-correctness.yml`

Full command examples for all of the above: `benchmark/README.md`.

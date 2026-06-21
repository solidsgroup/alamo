# Task 001: 3D wide-shallow Flame input files (roadmap 3.3)

## Goal
Author the canonical 3D Flame production input `input_3d_flame` plus grid-size
variants for the Phase-3 crossover sweep. These are the inputs that will run on NOVA
(A100/H200). They must follow the wide-shallow box strategy and use **analytic ICs only**.

## Context (only what this task needs)
- Phase 3 of the GPU optimization roadmap scales the problem 2D→3D to find the
  saturating regime. The box strategy (roadmap 3.3) is: **large base grid + shallow AMR
  + large blocking_factor** (e.g. 32 → 32768 cells/box). Do NOT use deep AMR or tiny boxes.
- `src/IC/BMP.H` is 2D-only — you may NOT use `*.ic.type = bmp`. Use `IC::Expression`.
  Expression syntax (confirmed): `<field>.ic.type = expression`, then
  `<field>.ic.expression.region0 = "<formula>"` with variables `x,y,z,t`; supports
  `sqrt`, `tanh`, `^`, `>`, `and`; constants via `<field>.ic.expression.constant.<name> = <val>`.
- The lead has produced a PROVEN, GPU-VERIFIED 3D template at repo-root `input_3d_smoke`
  (Expression ICs for `pf.eta` and `phi`, full 6-face eta/temp BCs, `elastic.type = disable`,
  NO elastic model/bc block). **Read it first and mirror its IC/BC/elastic structure
  exactly.** Your only job is to scale the grid up to production sizes.
- ELASTIC RULE (verified, critical): use **`elastic.type = disable`** and **OMIT** all
  `model_*` and `elastic.bc.*` lines. Do NOT use `elastic.on = 0` to disable it — that
  still runs the Static solve and crashes on device (CUDA-719). With `disable` present,
  any leftover `model_*`/`elastic.bc.*` entry makes ALAMO's strict parser abort as "unused".
  Per D1 the elastic solve is CPU-resident; Phase 3 measures the GPU-resident
  phase-field + thermal path only.

## Files to read first
- `input_3d_smoke` (repo root) — the proven 3D template; copy its IC/BC blocks verbatim.
- `input` (repo root) — the 2D reference for the physics/propellant/chamber blocks.
- `docs/agent_plans/20260620-gpu-phase3-regime-scaling/PLAN.md` — invariants & contract.

## Files allowed to modify (create only)
- `input_3d_flame`          — canonical production input
- `input_3d_flame_128`      — n_cell 128 128 64 variant
- `input_3d_flame_256`      — n_cell 256 256 128 variant
- `input_3d_flame_512`      — n_cell 512 512 256 variant

## Files NOT allowed to modify
- `input`, `input_3d_smoke`, anything under `src/`, `benchmark/`, or `bin/`.
- Any file owned by tasks 002/003.

## Implementation steps
1. Copy `input_3d_smoke` to `input_3d_flame`. Keep the Expression ICs and 6-face BCs.
2. Set the wide-shallow grid on `input_3d_flame`: `amr.n_cell = 256 256 128`,
   `amr.max_level = 1`, `amr.blocking_factor = 32`, a large `amr.max_grid_size` (e.g. 128),
   `geometry.prob_hi = 0.1754_m 0.1754_m 0.0877_m` (wide x-y, thin z). Keep `amr.node.all = 1`.
3. Keep `elastic.type = disable` (no model/bc block), `thermal.on = 1`, and the full
   propellant/chamber blocks from `input`.
4. Set `plot_file = output_3d_flame`, `amr.plot_int = -1`, `amr.thermo.int = 1`
   (sweeps override step counts on the command line).
5. Create the three size variants `_128 / _256 / _512` identical to `input_3d_flame`
   except for `amr.n_cell` (128 128 64 / 256 256 128 / 512 512 256). Note `_256` will
   duplicate the canonical grid — that is fine; the named variants give the sweep a
   uniform naming scheme.
6. Add a short header comment block to each file (purpose, grid, "Phase 3 / roadmap 3.3").

## Invariants (device, concurrency, architecture)
- Analytic ICs only (no BMP). All six faces (xlo..zhi) have BCs for eta and temp.
- Wide-shallow: max_level ≤ 1, blocking_factor = 32. No deep AMR.
- elastic.type = disable, elastic block omitted (CPU-resident per D1; avoids CUDA-719).
- Do not invent new ParmParse keys — only reuse keys present in `input` / `input_3d_smoke`.

## Build and test commands
- No build, no GPU. Validate by structural diff against `input_3d_smoke`:
  confirm every `*.ic.type` is `expression` (not `bmp`), eta/temp BCs list all of
  xlo,xhi,ylo,yhi,zlo,zhi, and `elastic.type = disable` is present.
  These greps must ALL return nothing: `grep -n "bmp" input_3d_flame*`,
  `grep -nE "elastic\.bc|^model_(prop|void|casing)" input_3d_flame*`.

## Expected result
Four runnable 3D input files at repo root, wide-shallow, analytic ICs, elastic off.

## Non-goals
- No physics tuning / no parameter sweeps of propellant constants.
- No scaling driver scripts (task 003). No SLURM (task 002).

## Stop conditions
- If `input_3d_smoke` is missing, STOP and record it (the lead must produce it first).

## Final report: write results/001-RESULT.md (use the RESULT template in PLAN context).

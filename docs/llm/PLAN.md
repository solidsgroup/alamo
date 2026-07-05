# chamber-gpu — Live Plan

The only live plan. Supersedes `GPU_ROADMAP_V3.md` and `GPU_STRUCTURAL_PLAN_20260703.md`
(both in `docs/archive/`, historical only). Written by REMEDIATION_PLAN Phase 1;
regenerate/edit by hand as tasks close — do not let this exceed 100 lines.

## Current phase & gate

**Phase 3 — elastic `Fapply` structural win**, on `chamber-gpu` /
`chamber-gpu-elastic-opt`. Phase 1 (physics-error-budget validation suite,
`benchmark/validate/`) is the standing gate: **no kernel optimization ships
without a physics-error-budget pass** (strict-build golden compare). Task 3.A
(C1 register edits) already cleared the CPU golden compare; every task below
must clear the same gate, plus an A100 before/after, before it lands.

## Next 3 tasks

1. **3.1 / v3 task 3.A — A100 A/B for the committed C1 edits.**
   Commit `0bb893acc` on `chamber-gpu-elastic-opt` (grad(C) single-`Matrix4`
   temp + boundary `sig` sink in `Operator::Elastic::Fapply`). Static sm_86
   register win confirmed locally; the A100 wall-time + occupancy A/B is the
   only thing still owed. Procedure: `benchmark/PHASE_C1_nova_ab.md`.

2. **3.2 / v3 task 3.F — Launch-bounds sweep on `Fapply`/`Diagonal`.**
   Sweep `__launch_bounds__(256, {1,2,3,4})` on the elastic kernels on A100:
   wall/step + ncu achieved-occupancy + budget gate per point. Statically
   spill-free to at least 80 registers. Use AMReX's existing
   `launch_global<MT, min_blocks>` overload scoped to the elastic kernels only
   — do not patch the global header.

3. **3.2b / v3 task 3.D (first wave) — Cheap kernel surgery in `Fapply`.**
   `src/Operator/Elastic.cpp` (`Fapply`, mirrored in `Diagonal`): (a) load
   `DDW(i,j,k)` once into a local instead of 3 separate loads; (b) column-
   restrict the `Cgrad_d x gradu` contractions (27 FMA instead of 81, x3);
   (c) hand-unroll `Matrix4<3,Major> x Matrix3` with direct `data[]` indexing
   instead of the 45-way if/else in `Matrix4_Major.H:552-565`. Judge by A100
   wall time + ncu executed-instructions, not register count. Bit-exact-able —
   CPU golden compare first, then the budget gate.

## Pointers

- **Status:** `benchmark/status.sh` (Phase 2 of this remediation) — run it, do
  not read a prose status file.
- **Device bug classes to avoid regressing:** `docs/llm/BUG_PATTERNS.md`.
- **NOVA/SLURM procedure:** `benchmark/NOVA_SLURM_RUNBOOK.md`.
- **Historical detail (do not re-read in normal sessions):** `docs/archive/`.

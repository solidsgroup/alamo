# chamber-gpu — Live Plan

The only live plan. Supersedes `GPU_ROADMAP_V3.md` and `GPU_STRUCTURAL_PLAN_20260703.md`
(both in `docs/archive/`, historical only). Written by REMEDIATION_PLAN Phase 1;
regenerate/edit by hand as tasks close — do not let this exceed 100 lines.

## Current phase & gate

**Phase 3 — elastic `Fapply` structural win**, on `chamber-gpu` /
`chamber-gpu-elastic-opt`. Phase 1 (physics-error-budget validation suite,
`benchmark/validate/`) is the standing gate: **no kernel optimization ships
without a physics-error-budget pass** (strict-build golden compare). Task 3.A
(C1 register edits) has now cleared both the CPU golden compare and the A100
before/after (255→244 reg/thread, ~14.5% Fapply wall win, occupancy flat —
see task 3.2); every task below must clear the same two gates before it
lands.

## Next 3 tasks

1. **3.2 / v3 task 3.F — Launch-bounds sweep on `Fapply`/`Diagonal`.**
   Sweep `__launch_bounds__(256, {1,2,3,4})` on the elastic kernels on A100:
   wall/step + ncu achieved-occupancy + budget gate per point. Statically
   spill-free to at least 80 registers. Use AMReX's existing
   `launch_global<MT, min_blocks>` overload scoped to the elastic kernels only
   — do not patch the global header. Now the designed next step, not just
   next-in-queue: task 3.1's A100 A/B (DONE 2026-07-07, see
   `docs/agent_plans/20260707-a100-ab-c1/results/RESULT.md`) found registers
   dropped 255→244 but occupancy stayed flat (~12%), i.e. still short of the
   2-block/SM threshold this lever targets.

2. **3.2b / v3 task 3.D (first wave) — Cheap kernel surgery in `Fapply`.**
   `src/Operator/Elastic.cpp` (`Fapply`, mirrored in `Diagonal`): (a) load
   `DDW(i,j,k)` once into a local instead of 3 separate loads; (b) column-
   restrict the `Cgrad_d x gradu` contractions (27 FMA instead of 81, x3);
   (c) hand-unroll `Matrix4<3,Major> x Matrix3` with direct `data[]` indexing
   instead of the 45-way if/else in `Matrix4_Major.H:552-565`; (d) in
   `Diagonal`, hoist the `DDW(i,j,k)` load above the per-component p-loop
   (currently re-fetched at Elastic.cpp:860/868, 2-3x per node). Judge by A100
   wall time + ncu executed-instructions, not register count. Bit-exact-able —
   CPU golden compare first, then the budget gate.

## Backlog (post next-3)

- **3.I — Fuse Newton convergence norms.** `Solver/Nonlocal/Newton.H` calls
  `MultiFab::norm0` per level x per component (lines ~419, ~572, and
  `FieldNorm0` ~816), each a separate device reduction + stream sync + MPI
  allreduce, multiplied by line-search backtracks. Replace with one fused
  `ReduceOps` pass per field (pattern: `Flame.cpp:616-654`). Convergence-
  semantics-critical: tier 3, CPU golden compare + budget gate.
- **Arena policy A/B (Phase 5).** No explicit arena selection exists in src/;
  managed-memory default is implicated in the multi-GPU loss. Measure device
  arena on A100 vs managed (A1000 needs managed for 8 GB).

## Pointers

- **Status:** `benchmark/status.sh` (Phase 2 of this remediation) — run it, do
  not read a prose status file.
- **Device bug classes to avoid regressing:** `docs/llm/BUG_PATTERNS.md`.
- **NOVA/SLURM procedure:** `benchmark/NOVA_SLURM_RUNBOOK.md`.
- **Historical detail (do not re-read in normal sessions):** `docs/archive/`.

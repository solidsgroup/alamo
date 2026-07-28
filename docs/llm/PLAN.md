# chamber-gpu — Live Plan

The only live plan. Supersedes `GPU_ROADMAP_V3.md` and `GPU_STRUCTURAL_PLAN_20260703.md`
(both in `docs/archive/`, historical only). Written by REMEDIATION_PLAN Phase 1;
regenerate/edit by hand as tasks close — do not let this exceed 100 lines.

## Current phase & gate

**Active work is the chamber-gpu-mem memory-strategy campaign, Phase 0** —
`docs/agent_plans/20260727-gpu-memory-strategy/PLAN.md` (v2.0) and its
deliverable record `results/PHASE0_v2.md`. Phases 1-3 of that campaign are
blocked on preconditions P1-P4, and the branch currently carries a live 3D
GPU correctness defect (`NOTES.md` N11, folder
`docs/agent_plans/20260728-elastic-psi-stencil-oob/`).

**`status.sh` changed 2026-07-28.** Its default run is orientation only and is
explicitly NOT a correctness gate; the sanitizer and gpu_strict legs run under
`FULL=1 bash benchmark/status.sh`. Use `FULL=1` before committing `src/` and at
every phase exit.

The Phase 3 material below is retained as the standing gate definition and as
the elastic-optimization backlog; it is not the active phase.

**Phase 3 — elastic `Fapply` structural win**, on `chamber-gpu` /
`chamber-gpu-elastic-opt`. Phase 1 (physics-error-budget validation suite,
`benchmark/validate/`) is the standing gate: **no kernel optimization ships
without a physics-error-budget pass** (strict-build golden compare). Task 3.A
(C1 register edits) has now cleared both the CPU golden compare and the A100
before/after (255→244 reg/thread, ~14.5% Fapply wall win, occupancy flat —
see task 3.2); every task below must clear the same two gates before it
lands.

GPU performance claims must use multi-step runs (10 steps by default) and
report external wall per step plus a startup calibration. Synchronized solver
or trace-region timers are supporting evidence when available; asynchronous
no-sync region attribution is not authoritative. Two-step decks remain
correctness smokes, not speed evidence.

## Next 3 tasks

1. **3.3 — FApply runtime follow-on. DONE 2026-07-22. RECOMMENDATION REFUTED
   2026-07-28 — DO NOT ADOPT 2/2.**
   Job `11772154` ran both arms on `input_copy` for 800 steps at production
   cadence. 4/4 completed 37 elastic solves (MLMG iters mean 157.41, max 289).
   **2/2 diverged on the first solve** — residual growing ~4.4× per iteration
   to 1.735e+20, then `amrex::Abort::0::MLMG failing so lets stop here`, dead
   at 6.5 s of an expected 431 s. Evidence:
   `docs/agent_plans/20260727-phase0-baseline/results/RESULT.md` §L6.
   The measurements below stand and the task correctly bounded itself to the
   frozen two-step horizon; the horizon simply cannot see this failure mode.
   All decks already set 4/4, so nothing shipped and nothing needs reverting.
   Original text follows.

   Superseded detail: the 2/2 gains were real on the frozen two-step cases
   (A1000 wall -12.8%/-21.0%, MLMG -24.0%/-23.4%; A100 wall -16.56%, MLMG
   -22.54%, FApply -23.30%) and no source change was retained. Step 3
   conservative specialization and Step 4 sequential Cgrad were measured and
   reverted; psi caching failed the memory gate. Full record:
   `docs/agent_plans/20260721-fapply-runtime-optimization/results/RESULT.md`.

2. **3.2 / v3 task 3.F — Launch-bounds sweep on `Fapply`/`Diagonal`.**
   Sweep `__launch_bounds__(256, {1,2,3,4})` on A100: wall/step + ncu
   occupancy + budget gate per point. Occupancy stayed flat (~12%) through 3.1
   and 3.2b — still 1 block/SM. **LOCAL LEG DONE 2026-07-13**: helper
   `src/Operator/ElasticLaunch.H`, knob `ALAMO_ELASTIC_MIN_BLOCKS` (default off
   = bit-identical), 3 sites wired, verifier CONFIRMED, branch
   `launch-bounds-sweep`. Fapply spills hard at min_blocks>=2 (128-reg cap vs
   254 live). REMAINING: A100 wall/ncu sweep (4 arms) + verdict + merge.

3. **3.2b / v3 task 3.D — Fapply/Diagonal kernel surgery. DONE 2026-07-13.**
   Merged cc520b4f8: Fapply exclusive wall -14.5%, MLMG::solve -10.6%,
   occupancy flat (win = spill/replay reduction). Includes Fsmooth 4D launch
   fusion. Record: `docs/agent_plans/20260713-fapply-322b-a100/`.

4. **GPU manual v3 — hostile transferability repair.**
   Active plan: `docs/agent_plans/20260721-gpu-manual-hostile-review/PLAN.md`.
   Separates invariant semantics, port contracts, and corpus examples; adds
   reusable templates, a stateful advisory scanner, and revision-bound per-port
   coverage. Current evidence is only 2/26 file-verified transforms, zero
   transfer-verified. Docs/scripts only; authorizes no source change.

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

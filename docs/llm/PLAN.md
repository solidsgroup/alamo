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

1. **3.3 — FApply runtime follow-on. DONE 2026-07-22.**
   Retain configuration-only 2/2 pre/post smoothing: on A1000 it reduced
   external wall by 12.8%/21.0% and MLMG solve by 24.0%/23.4% in the frozen
   2D-conservative/3D-psi two-step cases. Step 3 conservative specialization
   and Step 4 sequential Cgrad were measured and reverted; psi caching failed
   the memory gate; the target case already has one FApply launch per call, so
   multi-box launch fusion was not justified. Fresh review is clear. Evidence:
   `docs/agent_plans/20260721-fapply-runtime-optimization/results/RESULT.md`.
   NOVA/A100 confirmed the retained configuration: external wall -16.56%, MLMG
   solve -22.54%, and FApply -23.30%, with physics gate PASS. No source change
   is retained; the claim remains limited to the frozen two-step horizon.

2. **3.2 / v3 task 3.F — Launch-bounds sweep on `Fapply`/`Diagonal`.**
   Sweep `__launch_bounds__(256, {1,2,3,4})` on the elastic kernels on A100:
   wall/step + ncu achieved-occupancy + budget gate per point. Rationale:
   occupancy stayed flat (~12%) through both 3.1 and 3.2b — still 1 block/SM.
   **LOCAL LEG DONE 2026-07-13** — helper `src/Operator/ElasticLaunch.H`
   (`launch_global<MT, min_blocks>`, knob `ALAMO_ELASTIC_MIN_BLOCKS`, default
   off = bit-identical), 3 sites wired (Fapply/Diagonal/Fsmooth), CPU golden
   bit-exact, verifier CONFIRMED: commits 4d5289e67+532a757e3 on branch
   `launch-bounds-sweep`; sm_80 ptxas table in
   docs/agent_plans/20260713-launch-bounds-sweep/results/. Fapply spills hard
   at min_blocks>=2 (128-reg cap vs 254 live). REMAINING: A100 wall/ncu sweep
   (4 arms) + verdict + merge.

3. **3.2b / v3 task 3.D — Fapply/Diagonal kernel surgery. DONE 2026-07-13.**
   Merged to chamber-gpu (cc520b4f8) after A100 judgment PASS
   (docs/agent_plans/20260713-fapply-322b-a100/): Fapply exclusive wall
   -14.5% (299.0->255.7 s), MLMG::solve -10.6%, Fapply/launch -22-23%,
   occupancy flat (win = spill/replay reduction, registers 255->254).
   Parity: cell fields bit-identical; node fields within FP-reorder noise
   (strain_zx/zy 2.06e-6 rel = ~5e-9 abs on near-zero shear, adjudicated
   noise). Includes Fsmooth 4D launch fusion (MLMG-inclusive -10.6%).

4. **GPU manual v3 — hostile transferability repair.**
   Active plan: `docs/agent_plans/20260721-gpu-manual-hostile-review/PLAN.md`.
   Separate invariant semantics, port contracts, and corpus examples; add
   reusable scope/closure/inspection/validation/efficiency/harvest templates,
   single-home architecture policies, a stateful advisory scanner, and
   cross-family status. The hostile pass adds revision-bound per-port coverage,
   shape rather than identifier recognizers, a generalized value-dispatch
   contract, host-only numerical-kernel work, and mandatory layout/kernel/
   transfer/resource evidence before baseline efficiency. Current evidence is
   only 2/26 file-verified transforms, with zero transfer-verified; the first
   non-Flame contract instantiation remains an authorized pilot.
   This is docs/scripts-only and authorizes no source or numerical changes.

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

# TASK: rodtube-gate-recalibration

---

## Header

| Field        | Value                                                        |
|--------------|--------------------------------------------------------------|
| Risk tier    | 1 (decks + benchmark references; no src/)                    |
| Model        | opus (root-cause session; steps are tier-1 mechanical)       |
| Verification | full-oracle (benchmark/ci_golden_compare.sh)                 |
| Est. scope   | input_rod_and_tube_2d, input_rod_and_tube_3d, benchmark/baseline_references/rod_and_tube_step2/*, this folder |
| Parallel-safe| no (touches gate references used by status.sh)               |

## Operating rules

1. Read ONLY the files listed in Context budget. docs/archive/ is forbidden.
2. Every step's VERIFY must pass before acting on that step. On failure: STOP,
   report discrepancy, wait.
3. One commit per step unless stated. Message: `<area>: <what> (<task-folder>)`.
4. No scope expansion. New ideas go to NOTES.md in this folder, not into code.
5. If required knowledge is missing, ask the user.

## Context budget

Read first: benchmark/status.sh output, this PLAN.md
Read: input_rod_and_tube_2d, input_rod_and_tube_3d (tolerance blocks),
      benchmark/baseline_suite.py:24-70 (case defs),
      benchmark/baseline_references/rod_and_tube_step2/cpu.json,
      benchmark/ci_golden_compare.sh:27-75
Reference only if step names it: docs/agent_plans/20260718-rodtube-gate-recalibration/results/
Forbidden: docs/archive/*, unrelated task folders

## Objective

Golden gate is red: rod_and_tube_step2/cpu aborts `MLMG failed` at HEAD
(332ecffdd) but completes at origin (d964cfab8). Root cause (established
2026-07-18, this session): NOT a solver regression. The void-recovery operator
(7fbc18cc8 lineage) produces an honest warm-start residual (resid0 == bnorm =
1.05e7) where the old operator inflated resid0 to 3.35e8 (32x bnorm). AMReX
targets tol_rel * max(bnorm, resid0), so origin passed against a 32x looser
bar (180/200 iters, hairline). Both operators plateau at the same near-singular
seam limit (~4e-4 rel); HEAD reaches abs residual 1858 < origin's accepted
3353 yet fails its 32x tighter bar and hits max_iter -> amrex::Abort.
Evidence: with elastic.tol_abs=5000, HEAD step-2 physics is identical to
origin (pressure=970661, mdot=0.12694) and ends 100x MORE converged
nonlinearly (nonlinear_resid 1.96e6 vs origin's accepted 2.17e8; inexact
Newton, outer loop compensates). After this task: decks carry achievable
tolerances with documented rationale, cpu reference regenerated, gate green,
the 14 unpushed commits + recalibration commit pushed.

## Oracle

Command(s): bash benchmark/ci_golden_compare.sh  (exit 0, all cases ok);
            bash benchmark/status.sh  (golden-compare: PASS)
Covers: scalar thermo parity (L_max, area, pressure, disp_*, eta_min) vs
        regenerated references at recalibrated tolerances; run completion.
Does NOT cover: field-level displacement/stress parity (closed separately by
        Step 3 fcompare evidence in results/); GPU profiles (gpu_fast/
        gpu_strict refs go stale with deck edit -- regen deferred until 2D
        CUDA binaries exist; noted in results/RESULT.md).

## Steps

### Step 1 - task folder (this file)
DONE by definition of this commit.

### Step 2 - consumer census
VERIFY: n/a (read-only)
DO: grep -r rod_and_tube across tests/, benchmark/, .github/, slurm.
CHECK: consumer list recorded in results/RESULT.md.
RESULT (done 2026-07-18): live consumers = benchmark/baseline_suite.py +
baseline_references/rod_and_tube_step2/{cpu,gpu_fast,gpu_strict}.json only.
No tests/, no slurm, no CI yaml. 3D deck exists (input_rod_and_tube_3d),
referenced by no live harness.

### Step 3 - field-level parity evidence (kill-switch step)
VERIFY: origin worktree ~/Projects/alamo-origin-verify has built
        bin/alamo-2d-g++ at d964cfab8; fcompare.gnu.ex exists.
DO: run rod_and_tube_step2 with plotfiles on (amr.plot_int=1):
    (a) origin binary, stock deck; (b) HEAD binary, stock deck +
    elastic.tol_abs=5000. fcompare final-step plotfiles, all components.
CHECK: fcompare max rel error small (target <= ~1e-5, judgment on floor
    noise) for displacement/stress components. If fields DIVERGE: STOP --
    tolerance-artifact thesis is dead, fall back to bisecting the 14
    commits. Evidence -> results/fcompare_step3.txt.
RESULT (done 2026-07-18): ADJUDICATED PASS, kill-switch not triggered.
    Cell/flame fields: PLOTFILE AGREE (bit-identical). Node/elastic fields:
    origin<->HEAD differ 3-5% (stress rel 5.3e-2), BUT HEAD self-consistency
    across tol_abs 4000<->5000 is 0.1-0.35% -- 10-40x tighter -- and HEAD
    ends with 100x lower nonlinear residual (1.96e6 vs 2.17e8, rel 0.006 vs
    0.66). The 3-5% is origin's unconverged Newton tail; old reference was
    under-converged. tol_abs floor bracket: fail <=3000 (even max_iter=1000),
    pass >=4000; pick 5000. CONSEQUENCE for Step 5: regenerated cpu.json
    encodes better-converged elastic fields, not merely re-blessed old ones.
    Caveat: no independent truth oracle (no analytic solution); attribution
    rests on self-consistency ratio + nonlinear-residual ordering.

### Step 4 - recalibrate decks
VERIFY: Step 3 passed.
DO: input_rod_and_tube_2d: set elastic.tol_abs=5000 (keep tol_rel=1e-5 as
    secondary); REWRITE the stale comment block at lines ~166-185 ("achievable
    linear tolerance" text describes the pre-recovery operator) to record the
    resid0-honesty + inexact-Newton rationale and this task folder. Mirror in
    input_rod_and_tube_3d (same tol block pattern).
CHECK: HEAD run of gate command exits 0 on the 2D deck.

### Step 5 - regenerate cpu reference + gate green
VERIFY: Step 4 checked.
DO: baseline_suite.py record for rod_and_tube_step2/cpu only (do not touch
    gpu_*.json); run ci_golden_compare.sh.
CHECK: gate exit 0; status.sh golden-compare PASS.

### Step 6 - commit + push
VERIFY: Steps 4-5 green; diff reviewed (decks + cpu.json + task folder only).
DO: one commit `bench: recalibrate rod_and_tube tolerances for honest resid0
    (20260718-rodtube-gate-recalibration)`; push chamber-gpu (14 + this).
CHECK: git push succeeds; origin/chamber-gpu == local HEAD.

## Checkpoints

Tier 1: gates + spot-check diff before the Step 6 commit/push (user
requested push confirmation -- push only with user go-ahead).

## Closeout

- [ ] Oracle passes; status.sh golden-compare PASS
- [ ] results/RESULT.md: what changed, fcompare evidence, deviations
- [ ] changelog/ entry (append-only)
- [ ] touch results/DONE
- [ ] SESSION_LOG.tsv line appended
- [ ] worktree ~/Projects/alamo-origin-verify removed

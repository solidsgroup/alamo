# TASK: phi-quintic-mixing
# Folder: docs/agent_plans/20260723-phi-quintic-mixing/

---

## Header

| Field        | Value |
|--------------|-------|
| Risk tier    | 3 (changes the Flame constitutive material mixture) |
| Model        | Codex root session |
| Verification | partial-oracle |
| Est. scope   | one source file, approximately 8 lines; isolated commit |
| Parallel-safe| no: dirty worktree overlaps Flame and mechanics experiments |

## Operating rules

1. Read only the files listed in the context budget; `docs/archive/` is forbidden.
2. Preserve every pre-existing dirty-worktree change and stage only the quintic formula.
3. Do not alter `eta`, phase-field evolution, pressure loading, tolerances, AMR criteria, or the local `casing_support` experiment.
4. Use a unique candidate output path; do not overwrite the existing baseline.
5. Keep fixed plotting limits and compare MLMG iteration histories before accepting the physics verdict.

## Context budget

Read first: `benchmark/status.sh` output, `docs/llm/PLAN.md`, this `PLAN.md`
Read: `src/Integrator/Flame.cpp:1-25,510-570`,
      `input_rt1s_ideal`,
      `Makefile`,
      `.githooks/*`
Reference only for comparison:
      `output_rt1s_ideal_ncell64_casingAl_void0.5_0.5/`,
      task-owned results and run logs
Forbidden: `docs/archive/*`, unrelated task folders

## Objective

Replace the homogenized chamber's linear `phi` material weights with a clamped
quintic smoothstep while leaving `eta` weights unchanged. Preserve the local
`casing_support` partition and isolate the formula change in its own commit so
the existing and candidate runs form a clean A/B comparison.

## Oracle

Commands: focused 2-D rebuild, `bin/test-2d-g++`, `benchmark/status.sh`, and a
one-second `input_rt1s_ideal` candidate run with a unique output directory.
Covers: compilation, device lint, unit regressions, finite run completion,
MLMG-history comparison, and production-field generation.
Does NOT cover: human VisIt judgment at the prescribed fixed color limits; that
visual verdict remains a required acceptance step.

## Steps

### Step 1 - Implement and isolate the quintic weights

VERIFY:
```bash
git diff -- src/Integrator/Flame.cpp
benchmark/status.sh
```
DO: clamp nodal `phi_avg` to `[0,1]`, compute
`g=p^3(10-15p+6p^2)`, and replace only the three homogenized material weights.
CHECK: `git diff --check`, focused build/test, and device lint pass; staged diff
contains no pre-existing `casing_support` or cell-centering edits.

### Step 2 - Run the frozen A/B case

VERIFY: candidate executable and output path are distinct and the baseline is
readable.
DO: run the same `input_rt1s_ideal` to the same one-second horizon, overriding
only the candidate plot path; capture stdout/stderr and extract MLMG histories.
CHECK: zero exit, finite plot fields, and no meaningful iteration-count jump.

### Step 3 - Visual adjudication

VERIFY: both plotfiles use identical geometry, time, fields, and hierarchy.
DO: compare `P_thetar` at fixed `[-3.324e5, 3.282e5]` limits and compare old
versus quintic `kappa_mix` along the `phi` contour.
CHECK: retain the commit if the line ripple disappears or falls by more than
5x; otherwise retain the formula only as an explicitly approved modeling
practice and promote the lifecycle/staleness experiment.

## Checkpoints

- [x] After plan restatement: user supplied and approved the exact formula, isolation procedure, and pass criteria
- [x] Before commit: formula-only staged diff and build/lint evidence reviewed
- [ ] Before physics verdict: MLMG histories and fixed-scale plots reviewed

## Adversarial review

Review the isolated commit for device-safe clamping, preservation of partition
of unity with `casing_support`, unchanged `eta`, no tolerance/test weakening,
and no accidental resolved-AP/HTPB behavior change.

## Closeout

- [x] Automated oracle passes; status.sh all green
- [x] `results/RESULT.md` records change, evidence, and visual-verdict status
- [ ] `results/DONE` created when the visual verdict is complete
- [ ] One outcome line appended to `docs/llm/SESSION_LOG.tsv`

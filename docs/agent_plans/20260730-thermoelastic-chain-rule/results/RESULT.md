# Thermoelastic chain-rule repair result

## Outcome

The accepted thermoelastic repair was applied to `chamber-gpu` and verified
on the final branch base `21a4c06c368237e2aa1c17d6dd2b2055ada0c53a`.
The constitutive correction, deterministic regression coverage, and the
minimal GPU source-closure update are ready for the task's scoped commit.

The original dirty worktree at `/home/jackplum/Projects/alamo` was not
modified.

## Changes

- Corrected `DW` for `W(F F0^-1)` by right-transforming the elastic stress
  with `F0^-T`.
- Corrected `DDW` by applying the two required `F0^-1` contractions.
- Reused one inverse helper, returned the derived predeformed type from
  `Random()`, and generated a nonsingular perturbation of identity.
- Exported all nine `F0` components in 3-D.
- Added deterministic 2-D/3-D finite-difference derivative,
  stress-free-expansion, and field-layout regression tests.
- Added `InputScraper.cpp` and `OutputLog.cpp` to the Flame CUDA source
  closure. The branch's concurrent development re-merge introduced direct
  references to those translation units, and the official strict CUDA link
  otherwise failed with undefined symbols.

## Verification

| Check | Result |
|---|---|
| Full 2-D test executable | PASS, 0 failures |
| 2-D focused derivatives | `DW = 2.57656e-08`, `DDW = 1.02470e-09` relative error |
| Full 3-D test executable | PASS, 0 failures |
| 3-D focused derivatives | `DW = 1.60469e-08`, `DDW = 1.28394e-09` relative error |
| Free expansion | PASS, normalized stress `0` in 2-D and 3-D |
| `benchmark/status.sh` | PASS: device lint, CPU golden/smoke, A100 sanitizer |
| Local A100 gate tiers 1 and 2 | PASS; memcheck `ERROR SUMMARY: 0 errors` |
| Official strict CUDA build | PASS after completing the source closure |
| Strict GPU valid golden cases | PASS exactly: canonical steps 1/2 and eta-expression step 1 |
| Strict GPU NaN smoke | PASS, two steps and clean AMReX finalization |
| `git diff --check` | PASS |
| Fresh adversarial review | PASS, no material findings |

Evidence is retained in this directory:

- `test-2d.log`
- `test-3d.log`
- `status-final.log`
- `a100-gate.log` and `a100-gate/`
- `gpu-strict-golden-final.log`
- `gpu-strict-nan-smoke/run-final.log`

## Strict-reference adjudication

The strict comparison still reports traction-only differences for
`rod_and_tube_step2`. This is the pre-existing stale GPU reference documented
in `docs/llm/SESSION_LOG.tsv` on 2026-07-26. That case has `F0 = I`, so the
new chain-rule terms reduce to the prior constitutive path and cannot cause
the discrepancy. Three other strict cases compare exactly, the CUDA smoke is
clean, and the A100 runtime and memcheck gates pass. The stale reference is
therefore not a code-relevant failure for this repair.

## Execution notes

The local `chamber-gpu` base advanced from `2a3d20c6` to the corrected
development re-merge `21a4c06c3` while testing was in progress. The update
was preserved, the source closure was reconciled, and all final tests above
were run against the newer base.

The initial baseline A100 check could not run because the isolated worktree
did not yet have the expected CUDA binary. Building that binary resolved the
environmental prerequisite; the subsequent runtime and memcheck gates passed.

# RESULT — rodtube-gate-recalibration (2026-07-18)

## What changed

- `input_rod_and_tube_2d`: `elastic.tol_abs` 1e-8 -> 5000; stale "achievable
  linear tolerance" comment block rewritten with the honest-resid0 +
  inexact-Newton rationale, floor bracket, and long-run dimensional caveat.
- `input_rod_and_tube_3d`: mirrored `tol_abs=5000` with explicit
  "3D floor unmeasured" caveat (no live consumer).
- `benchmark/baseline_references/rod_and_tube_step2/cpu.json`: regenerated.
  Diff vs old reference: ONLY the 4 trac_* boundary tractions
  (trac_xhi_x -60959.1 -> -61097.5; trac_xhi_y 19.05 -> 19.51;
  trac_yhi_x 1.81 -> 2.46; trac_yhi_y -59673.4 -> -59799.0).
  All flame/chamber columns bit-identical to origin (max_rel = 0.0).
- gpu_fast.json / gpu_strict.json NOT regenerated (2D CUDA binaries absent);
  known-stale until rebuilt.

## Why (root cause, evidence)

Golden gate red at HEAD 332ecffdd: rod_and_tube_step2/cpu hit
`amrex::Abort MLMG failed` (max_iter=200, resid plateau 1.8e-4 rel). NOT a
solver regression: void-recovery operator reports honest warm-start resid0
(== bnorm 1.05e7) where the old operator inflated it 32x (3.35e8), so the
old `tol_rel=1e-5` bar (targeted at max(bnorm, resid0)) was accidentally
32x looser -- origin passed at 180/200 iters against it. Seam has an
ABSOLUTE residual floor ~4e3 (fails <=3000 even at max_iter=1000;
passes >=4000).

Evidence (results/fcompare_step3.txt):
- flame/cell fields origin vs HEAD: PLOTFILE AGREE (bit-identical)
- elastic fields: origin<->HEAD 3-5%; HEAD self-consistency across
  tol_abs 4000<->5000 = 0.1-0.35% (10-40x tighter) -> delta is origin's
  unconverged Newton tail (origin accepted at nonlinear_resid_rel=0.66;
  HEAD ends at 0.006, 100x tighter)
- repeatability: gate command twice at HEAD -> thermo.dat byte-identical
- deck-edit run == CLI-override run: thermo.dat byte-identical

User adjudicated 2026-07-18: elastic-field movement is expected accuracy
improvement; recalibrated state adopted as new golden anchor.

## Gate state after

ci_golden_compare.sh: PASS (all 4 cases ok). status.sh: device-lint PASS,
golden-compare PASS, a100-sanitizer PASS.

## Deviations from plan

- Step 3 target "<=1e-5 field parity" not met literally (elastic 3-5%);
  adjudicated PASS via self-consistency + residual-ordering attribution and
  explicit user sign-off. Kill-switch (bisect fallback) not triggered.
- Discovered en route: gate's elastic visibility = trac_* columns only
  (upgrade idea logged in NOTES.md, not implemented).

## Open / follow-ups

- Adversarial review of the recalibration commit requested (tolerance change
  + reference regen in one commit = reviewer target pattern).
- GPU reference regen once 2D CUDA binaries exist; full GPU suite run
  (Track 1b) pending.
- Track 2 ticket: MLMG abort semantics + t=5.58s coarse-cap hypothesis test
  (NOTES.md).
- Worktree ~/Projects/alamo-origin-verify removal at closeout.

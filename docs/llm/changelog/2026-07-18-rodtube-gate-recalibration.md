# rod_and_tube golden-gate recalibration — honest resid0 after void-recovery — 2026-07-18

Branch `chamber-gpu`. Task folder:
`docs/agent_plans/20260718-rodtube-gate-recalibration/` (evidence in
`results/fcompare_step3.txt`, narrative in `results/RESULT.md`).

## TL;DR

| Finding | Class | Status |
| --- | --- | --- |
| rod_and_tube_step2/cpu golden gate red at HEAD (`amrex::Abort MLMG failed`, max_iter=200) | tolerance-calibration artifact, NOT a solver regression | **fixed** — `elastic.tol_abs` 1e-8 → 5000 in both rod_and_tube decks; cpu reference regenerated |
| Old operator inflated warm-start resid0 to 32× bnorm; `tol_rel=1e-5` (targeted at max(bnorm, resid0)) was accidentally a 32× looser bar | latent gate weakness | documented; origin passed at 180/200 iters against the loose bar |
| Seam has an absolute residual floor ~4e3 (fails ≤3000 even at max_iter=1000; passes ≥4000) | conditioning fact | encoded as tol_abs=5000; inexact Newton compensates — end state 100× more converged nonlinearly (1.96e6 vs 2.17e8) |
| Elastic fields moved 3-5% vs old reference | accuracy improvement (old reference under-converged at nonlinear_resid_rel=0.66) | adopted as new golden anchor with user sign-off; flame/cell fields bit-identical |

## Notes

- Gate's elastic visibility = 4 trac_* thermo columns only; field-level
  reference upgrade idea in task NOTES.md.
- gpu_fast/gpu_strict references known-stale until 2D CUDA binaries rebuilt.
- Follow-up ticket material (task NOTES.md): MLMG abort semantics
  (return-best on max_iter), t=5.58s coarse-cap hypothesis test, and the
  observation that the binding convergence metric is not the Fine residual.

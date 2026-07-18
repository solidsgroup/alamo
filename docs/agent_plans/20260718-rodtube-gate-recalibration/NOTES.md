# NOTES (out-of-scope ideas, not implemented)

## Field-level golden reference for rod_and_tube_step2 (awaiting user go-ahead)

Finding from Step 3: the thermo-scalar gate sees elastic state only through
the 4 trac_* boundary integrals. Full-field elastic drift (5% stress) showed
up there as 2e-3..2.6e-1 rel on trac columns -- detectable, but coarse and
partly masked (trac_yhi_x is a near-zero shear component where rel is noise).

Upgrade: archive the final-step node plotfile (00002node) as a reference
artifact and add an fcompare leg to the case in baseline_suite.py
(fcompare.gnu.ex already built in ext/AMReX-Codes/amrex/Tools/Plotfile/).
Would make this deck a true elastic field gate. Cost: ~few MB reference
binary in-repo (or LFS/checksum), one subprocess call in the harness.

## Track 2 (separate ticket, agreed 2026-07-18)

- MLMG hard-abort on max_iter: return-best + let Newton accept decide.
  BEFORE designing: check whether AMReX MLMG has a no-abort/return-best knob.
- Cheap hypothesis test: re-run the t=5.58s coarse-cap chamber deck with
  relaxed elastic.tol_abs -- if it clears 5.58s, the campaign wall is the
  same honest-resid0 artifact class as this gate failure.
- Observed during bracketing: tol_abs=3000 fails even at max_iter=1000
  although the Fine-resid trajectory crosses 3000 by ~iter 110 at 200-iter
  pace. Binding metric is NOT the Fine residual (coarse-AMR-level residual or
  AMReX stall detection). Worth understanding before Track 2 design.

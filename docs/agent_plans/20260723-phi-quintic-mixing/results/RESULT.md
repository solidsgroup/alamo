# Quintic phi mixing result

Status: automated checks passed; fixed-scale VisIt verdict and final commit are pending.

## Change

The homogenized chamber material blend now clamps `phi_avg` to `[0,1]` and
uses

`g(phi) = phi^3 (10 - 15 phi + 6 phi^2)`

for the solid, void, and casing weights. The working-tree implementation keeps
the pre-existing `casing_support` multiplier on `g` and `1-g`, preserves
partition of unity, and leaves the `eta` factors unchanged. The index contains
only the formula/include delta so the eventual commit will not absorb the
overlapping cell-centered-phi or casing-support work.

## Automated evidence

- Focused 2-D build: `make -j4 bin/test-2d-g++ bin/alamo-2d-g++` passed.
- Unit executable: `bin/test-2d-g++` passed with zero failures.
- Device lint: `benchmark/status.sh` passed after the edit and after the run.
- Scoped whitespace check passed for the source and task files.
- Production run:
  `mpiexec -n 4 bin/alamo-2d-g++ input_rt1s_ideal plot_file=output_rt1s_quintic_phi_t1`
  exited zero at step 5000, time 1.0 s, and finalized AMReX normally.
- Final plotfiles exist at `output_rt1s_quintic_phi_t1/05000node` and
  `output_rt1s_quintic_phi_t1/05000cell`.
- Candidate and baseline final headers match in time, domain, two-level AMR
  hierarchy, and cycle metadata (`5000 10000 20000`).
- The run log has no NaN/Inf, abort, divergence, fatal, or failed-convergence
  marker.
- All 49 mechanics solves converged in exactly two Newton iterations. This
  matches `run_rt1s_ideal.log` and the prior cell-centered candidate log.
  MLMG per-cycle verbosity is disabled in these logs, so there are no direct
  MLMG iteration records to compare.

Full log: `results/run_t1.log`.

## Required visual verdict

In VisIt, compare the final candidate against
`output_rt1s_ideal_ncell64_casingAl_void0.5_0.5`:

1. Plot `P_thetar` with fixed limits `[-3.324e5, 3.282e5]`.
2. Compare the old conductivity expression
   `0.162*phi + 70*(1-phi)` against
   `0.162*g + 70*(1-g)`, where
   `g=(clamp(phi,0,1)^3)*(10-15*clamp(phi,0,1)+6*clamp(phi,0,1)^2)`.
3. Accept the mechanism if the `phi ~= 1` line ripple disappears or falls by
   more than 5x. The formula-only change is staged but deliberately uncommitted
   until this verdict.

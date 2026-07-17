# Elastic regression recovery result

## Status

**BLOCKED as of 2026-07-16.**  The paired-stencil regression and several
independent safety/halo defects were repaired, but the final finite-soft-void
chamber does not meet the force-balance oracle in either 2-D or 3-D.  No commit
was made and `DONE` must not be created.

| Repository state | Value |
|---|---|
| Worktree | `/tmp/alamo-elastic-void-continue` |
| Branch | `elastic-void-wip-20260713` |
| HEAD | `a06c6f15d` (`WIP: elastic void robustness -- paired-stencil rewrite`) |
| Trusted comparison parent | `d964cfab8` |
| Current configured target | 2-D production GCC |
| Hook path | `.githooks` |
| Device lint | PASS, 2026-07-16 |

## Completed work

- Same-base attribution proved that the paired `D- G+` rewrite at
  `a06c6f15d`, not the trusted parent, caused false-zero solves, NaNs, and
  mechanics drift.  The full parent operator/stencil stack was restored while
  retaining the parent's existing `elastic.use_psi=0` material path.
- Newton now rejects nonfinite states, restores the baseline after a failed
  line search, never forces final-backtrack acceptance, aborts genuine
  positive-tolerance exhaustion, and preserves `nrtolerance=0` as explicit
  fixed-iteration mode.  Direct unit cases cover its decision boundaries.
- Constant Neumann boundary spacing now remains `amrex::Real`; a direct unit
  test proves nonzero subcell offsets at `dx=0.25`.
- Flame source/model fills now honor periodicity and construct the complete
  allocated model halo.  Elastic computes the diagonal on the two-ring active
  region consumed by its smoother, while generic normalization divides valid
  rows only.  This eliminates the 3-D periodic/mixed-corner false-zero/NaN
  correction.
- `SCPSpheresElastic` now rejects nonfinite mechanics and a nonzero-RHS/zero-
  update false success.  `RubberPlateHole`'s coverage-only fixed iteration
  explicitly opts into `nrtolerance=0` without changing the physical section.
- `ElasticSoftVoid` provides 2-D serial, 2-D MPI, and 3-D sections with
  `elastic.use_psi=0`, `psi_floor=0`, a 0.2 MPa constituent, deterministic AMR
  interfaces, exact mixture checks, realized-modulus checks, finite solver
  diagnostics, periodic-extrusion checks, and plotted force balance.

## Verification outcome

| Gate | Result | Evidence |
|---|---|---|
| Final 2-D GCC build | PASS | `make -j4`; all nine targets linked after closeout cleanup |
| Final 2-D unit executable | PASS | `bin/test-2d-g++`; zero failures |
| Prior 3-D unit executable | PASS | Included new Newton and Constant tests before diagnostic-only cleanup |
| Fresh 2-D soft void, serial | run PASS / check FAIL | `tests/ElasticSoftVoid/output_2026-07-16_15.26.33_kermit_2d-serial`; equilibrium ratio `0.270767` |
| Fresh 2-D soft void, 2 MPI ranks | run PASS / check FAIL | Same test ID, `2d-parallel`; equilibrium ratio `0.270767` |
| Final retained 3-D soft void | run PASS / check FAIL | `output_2026-07-13_19.02.34_kermit_3d-serial`; equilibrium ratio `0.270142` |
| 3-D strict update-tolerance screen | Expected abort | Five iterations end at update `1.40792e-6`, nonlinear residual ratio `0.139951`, tolerance `1e-7` |
| Diff format | PASS | `git diff --check` at closeout |

The fresh 2-D plot realizes `mu_min=0.644443 MPa` and
`kappa_min=0.676235 MPa`; the retained 3-D plot realizes
`mu_min=0.722955 MPa` and `kappa_min=0.760362 MPa`.  The 3-D final field is a
valid periodic extrusion (`rhs_z=0`, `max|disp_z|=8.99e-20`).  Thus the failure
is not a stiff void, psi regularization, MPI decomposition, nonfinite value, or
periodic-z artifact.  Update-based stopping accepts after two iterations while
the nonlinear residual remains about 24% of its initial value and plotted
force imbalance remains about 27% of the applied RHS.

## Unfinished gates

- Localize and repair the residual/Jacobian/boundary/AMR inconsistency exposed
  by the deterministic chamber; the linear material should not require this
  slow nonlinear residual decay.
- Make all three `ElasticSoftVoid` checks pass without increasing the five-
  iteration cap, relaxing tolerances, adding stiffness/psi floors, or forcing
  acceptance.
- Rebuild and run complete final-source 2-D and 3-D GCC suites plus `make test`.
- Obtain an independent Tier-3 diff review.  The attempted closeout reviewers
  were unavailable because the subagent quota was exhausted.
- Restore a usable CUDA toolkit, then run strict/fast 2-D and 3-D correctness,
  comparison, sanitizer/device-lint, and the specified CPU/GPU timing campaign.

See `../NOTES.md` for the chronological almanac, `REVIEW.md` for blocking
findings, and `HANDOFF.md` for exact continuation instructions.

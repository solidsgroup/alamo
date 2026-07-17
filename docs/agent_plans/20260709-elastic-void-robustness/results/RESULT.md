# Result: finite soft void elastic solve

Status: focused Option 1 oracle passed, but full-suite review blocks completion; changes remain deliberately uncommitted.

> **Superseded completion claim (2026-07-13).** The full GCC review found 13 runtime failures and 17 numerical check failures. This focused result is retained as evidence only; see FULL_SUITE_REVIEW_20260713.md and REGRESSION_RECOVERY_HANDOFF.md before treating the work as complete.

## Decision

The physical regression uses a finite, positive soft-void material rather than
a zero-support phase mask:

```text
elastic.use_psi = 0
elastic.psi_floor = 0
model_void.kappa = model_void.mu = 0.2 MPa
```

This is not a psi floor.  The soft material remains part of the elastic
operator, so there are no inactive displacement rows in the tested problem.

## Mechanism and correction

The old nonlinear residual composed collocated central gradient and divergence
operators, while the MLMG `Fapply` used a different direct second-derivative
operator.  Thus the supplied linear operator was not the residual's Jacobian.
A literal central-composition tangent was also not usable: its collocated
gradient has a checkerboard nullspace.

The correction changes both the residual and tangent to the same conservative
paired `D- [ DDW G+ ]` discretization.  `Elastic::Diagonal` is the matching
local impulse diagonal, and the smoother/normalization operate on the same
valid plus first-ghost active rows.  The implementation also supplies the
paired stencil's coefficient/phase-field collars across AMR and MG levels,
keeps initialized C/F anchor data from being overwritten by multi-ghost
exchange, and uses a bounded phase interpolation only for auxiliary support
reads.

For a raw psi-masked solve, a true zero diagonal is now rejected rather than
allowing smoother normalization to create NaN/Inf.  That is a diagnostic for
the unsupported exact-zero-mask formulation, not regularization.

## Validation matrix

| Check | Result |
|---|---|
| Final 2-D CPU release oracle: `scripts/runtests.py --dim=2 --serial --comp=g++ --timeout 180 --no-backspace --no-clean tests/ElasticSoftVoid` | PASS.  AMR level 2; final linear relative residual `8.90303e-06`; final nonlinear relative residual `8.97774e-06`. |
| Final 3-D CPU release oracle: same command with `--dim=3` | PASS.  AMR level 1; final linear relative residual `9.37483e-06`; final nonlinear relative residual `9.44120e-06`. |
| 2-D, two-rank MPI direct smoke, with AMR level 2 and `0.2 MPa` void | PASS; same final nonlinear residual `8.97774e-06`. |
| 2-D void-stiffness sensitivity, `1 MPa` void | PASS; final nonlinear relative residual `9.50446e-06`. |
| 2-D resolution sensitivity, `amr.n_cell='48 48 8'`, AMR level 2 | PASS; final nonlinear relative residual `5.83449e-07`. |
| Property oracle | Confirms no solver failure/non-finite values, residual targets, nontrivial finite displacement/stress/work, positive material contrast, requested AMR level, and exact `use_psi=0`/`psi_floor=0` metadata. |
| Static checks | `git diff --check`, Python compilation of the test oracle, and `benchmark/status.sh` passed (`device-lint: PASS`). |

## Limitations and follow-up scope

- A raw exact-zero psi mask remains a different formulation: it needs an
  inactive-DOF/nullspace policy.  It is not claimed as passing here, and the
  physical regression does not use a floor.
- CPU release builds were exercised in 2-D and 3-D.  CUDA/GPU execution was
  not rerun after the final changes.
- The optional cross-AMR coefficient-average-down path and AMR refinement
  ratio 4 were not accepted as regression scope.  A ratio-4 uniform probe did
  not converge within the configured iteration budget, whereas the supported
  regression uses ratio 2.
- The paired operator requires two displacement/coefficient ghosts.  Newton
  checks that contract; callers of `Fapply` outside that path should preserve
  it explicitly.
- A direct solid-region comparison with a floored/masked reference would
  compare different physical models, so it is not an acceptance oracle for
  Option 1.

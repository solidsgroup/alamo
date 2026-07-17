# Adversarial review

> **Status update, 2026-07-13:** the focused review below is superseded as an
> acceptance decision. The full GCC regression review found blocking
> correctness failures. Read FULL_SUITE_REVIEW_20260713.md and REGRESSION_RECOVERY_HANDOFF.md before
> relying on any earlier “no blocker” or “oracle weakening: no” conclusion.


Scope: final finite-soft-void (Option 1) vector-mechanics path.

## Findings

| Concern | Finding |
|---|---|
| Accidental psi regularization | No.  The regression metadata requires `elastic.use_psi=0` and `elastic.psi_floor=0`; its void is represented by finite positive material coefficients. |
| Residual/Jacobian consistency | Addressed.  Residual assembly, `Fapply`, and the diagonal use the paired conservative `D- G+` stencil; the final Newton residual, not merely update size, is the acceptance condition. |
| Smoother/diagonal row mismatch | Addressed.  `Fapply`, `Elastic::Diagonal`, `Fsmooth`, and normalization use valid nodes plus the first active ghost ring. |
| AMR/MG support reads | Addressed for the supported two-ghost contract.  DDW and psi collars are initialized before restriction; only valid/periodic data are later exchanged so a peer ghost cannot replace an initialized C/F anchor. |
| MPI/periodic behavior | No confirmed blocker.  The two-rank 2-D smoke passed, and the coefficient handoff uses the coarse geometry's periodicity. |
| Boundary conditions | Intentional discretization change.  The shared paired flux requires a one-sided normal boundary derivative with centered tangential closure; retaining the former incompatible boundary derivative would break the exact tangent relationship. |
| Oracle weakening | No.  The test checks reported linear residuals, residual-based Newton convergence, completion status, finite fields, nontrivial response, AMR depth, and material contrast. |

## Remaining non-blocking items

- Exact-zero raw psi masks remain unsupported pending an inactive-DOF/nullspace
  policy.  The zero-diagonal guard is not a substitute for that policy.
- The current regression does not exercise CUDA after the final source update.
- The generic bounded-gradient helper cannot make an undersized displacement
  FAB safe at the anchor itself; the existing two-ghost precondition must be
  retained by any future caller.
- The scalar Newton overload retains a legacy multi-ghost final DDW fill,
  outside the reviewed vector mechanics path.
- Restriction-temporary lifetime under an explicit GPU `NoSyncRegion` is a
  hardening opportunity, not a normal CPU/MPI correctness failure.

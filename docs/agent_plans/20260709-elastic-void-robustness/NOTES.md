# Investigation Notes

## Reproduced Failure

All commands use the existing 2D CPU binary and force an elastic solve at
`t=1e-4` with `max_step=2 elastic.interval=1`.

| Configuration | Result |
|---|---|
| `use_psi=1`, `psi_floor=0`, void `1/1 MPa` | MLMG aborts on Newton step 2 after 200 iterations at relative residual `1.91e-4`. |
| `use_psi=0`, void `1/1 MPa` | Two linear solves converge in 11 and 10 iterations. |
| `use_psi=0`, void `0.2/0.2 MPa`, interface width `1.5e-3 m` | MLMG diverges. |
| `use_psi=0`, void `0.2/0.2 MPa`, interface width `3e-3 m` | Linear solves converge in 8 and 7 iterations. |

The masked formulation scales the already-soft void tangent by eta. Exact-zero
eta therefore leaves the operator with the implementation's `1e-8` offset,
rather than the chosen void modulus. The unmasked formulation uses the finite
void constitutive model directly. The paper requires material discontinuities
to be diffuse and sufficiently resolved; the original width resolves to only
about 2.2 finest cells, while the successful width resolves to about 4.4.

## 3D Evidence

The 3D CUDA binary completed a half-base-resolution extruded rod-and-tube
screen (`32 32 8`, width `6e-3 m`, preserving the 4.4 finest-cells interface
width) with `use_psi=0`, `psi_floor=0`, and a `0.2/0.2 MPa` void. Each of ten
Newton linearizations reached the `1e-5` linear tolerance in 4--11 MLMG
iterations. The full 64x64x16 screen is too close to the 8 GB local GPU limit
for a reliable execution oracle and will need a higher-memory device or CPU
3D build for the final full-resolution run.

## Newton/Operator Mismatch

The nonlinear residual is computed in `Newton::prepareForSolve` as
`b - Divergence(DW(Gradient(u)))`. Its central first derivative followed by a
central divergence produces a two-cell-wide second-derivative stencil.
`Elastic::Fapply` instead uses direct one-cell second derivatives in
`gradgradu`; this is the linear operator given to MLMG.

The mismatch is observable with a uniform tangent and a 1000x smaller load:
MLMG reports a relative linear residual of `2.06e-6`, but the Newton residual
recomputed immediately afterward is `0.399983` of its initial norm. This also
happens with no AMR levels, so it is not a composite-norm artifact. A Newton
line-search can then accept a final forced backtrack and the update-only gate
reports convergence despite a large residual.

The implementation must align the Jacobian and residual discretizations before
the no-floor configuration can be treated as a correct nonlinear solve. The
candidate correction is to make `Elastic::Fapply` (and its smoother diagonal)
the derivative of the existing discrete divergence residual, then add an
oracle that rejects a nonconverged line search rather than accepting it as a
small update.

## Restart Handoff

Read this file and PLAN.md first, then inspect these source anchors:

- `src/Integrator/Flame.cpp:456-465`: psi construction.
- `src/Integrator/Flame.cpp:321-330`: `elastic.use_psi` installation.
- `src/Operator/Elastic.H:145-148`: the independent `m_psi_small=1e-8` mask
  offset.
- `src/Solver/Nonlocal/Newton.H:267-317`: nonlinear residual assembly
  `b - Divergence(DW)`.
- `src/Solver/Nonlocal/Newton.H:389-456`: Newton solve and forced final
  line-search acceptance.
- `src/Operator/Elastic.cpp:513-630`: Jacobian `Fapply`, which uses direct
  second derivatives and a continuum product-rule split.
- `src/Operator/Elastic.cpp:790-900`: hand-derived diagonal that must change
  with any `Fapply` stencil change.

Exact recorded commands and logs:

```bash
# Masked zero-floor failure
bin/alamo-2d-g++ input_rod_and_tube_2d max_step=2 elastic.interval=1 \
  elastic.use_psi=1 elastic.psi_floor=0 \
  model_void.kappa=1_MPa model_void.mu=1_MPa

# No-mask soft-void success with resolved interface
bin/alamo-2d-g++ input_rod_and_tube_2d max_step=2 elastic.interval=1 \
  elastic.use_psi=0 elastic.psi_floor=0 \
  pf.eta.ic.expression.constant.w=0.003 \
  model_void.kappa=0.2_MPa model_void.mu=0.2_MPa

# 3D CUDA confirmation, reduced only for the local 8 GB GPU
bin/alamo_gpu-3d-cuda86-g++ input_rod_and_tube_3d max_step=2 \
  amr.n_cell='32 32 8' elastic.interval=1 elastic.use_psi=0 \
  elastic.psi_floor=0 pf.eta.ic.expression.constant.w=0.006 \
  model_void.kappa=0.2_MPa model_void.mu=0.2_MPa
```

Logs remain under `/tmp/alamo-elastic-void-screen/`; do not treat their
existing prose as authority. Rerun the commands after any source change.

Before editing source, decide and document one of these mutually exclusive
formulations:

1. Preserve `Newton::prepareForSolve`'s discrete `D(DW(Du))` residual and make
   `Fapply`/`Diagonal` its exact tangent. This is the preferred candidate
   because it gives Newton a true Jacobian, but it changes the linear stencil
   from the paper's direct second-derivative product-rule implementation.
2. Preserve the paper's product-rule linear operator and rederive a nonlinear
   residual whose exact derivative is that operator. This is more invasive and
   has no identified generic construction for arbitrary finite-strain models.

Do not merely relax `tol_rel`, increase `nriters`, or retain a psi floor: those
hide the inconsistency and do not meet the objective.

## Step 3 authorization and invariant

On 2026-07-09, the user asked to resume the most recent elastic-solver work.
Treat that as approval to begin the bounded source correction.  The invariant is
that, at every interior node, `Elastic::Fapply(v)` must equal the directional
derivative of the residual assembly in `Newton::prepareForSolve`:
`D_h(DDW(G_h(u)) G_h(v))`, with the *same raw psi weighting* when psi is
installed.  The elastic diagonal must be the diagonal of that same discrete
operator.  No hidden `m_psi_small` regularization may remain in this path.

The initial focused oracle is the existing one-level uniform-coefficient probe:
after one Newton update its recomputed nonlinear residual must fall from the
current `0.399983` relative value to the linear-solve tolerance scale, rather
than merely satisfying the update gate.

## Step 3 trial: central composition is not solver-admissible

The literal exact tangent trial, `D_c(DDW * G_c(v))`, compiled in the 2D
build and its sinusoidal `Fapply` probe agreed with the derived closed form to
`2.93e-15` relative error.  It nevertheless made the one-level uniform case
diverge in MLMG (relative linear residual grew from `2.56e-1` to `1.31e21` in
15 iterations).  This is a formulation result, not a coding discrepancy:
the composition of collocated central first derivatives has a checkerboard
nullspace (`G_c((-1)^i)=0`).  It cannot be used as the solve operator even
though it is the literal derivative of the present residual.

Do not keep that trial implementation.  The next diagnosis must choose a
non-degenerate paired discrete gradient/divergence formulation (and change
both residual and Jacobian together), or derive a nonlinear residual for the
existing one-cell product-rule operator.  The former needs a stencil and
boundary-condition review before another source edit.  An exact-zero
psi-masked domain also has genuine zero rows; it needs an explicit inactive
DOF policy rather than a floor, so the current focused physical oracle remains
`elastic.use_psi=0` with the soft nonzero void model.

## Step 3 paired-stencil specification

The selected correction is the conservative, collocated forward/backward
pair.  Let `Gpair(u)_i` use the forward difference in each coordinate, except
for the backward one-sided difference at a physical high boundary.  Use that
same gradient at *every* node, including boundary-condition rows: a boundary
node's stored stress is also the flux used by its adjacent interior row.  Let
`theta_i` be one without psi or the raw nodal psi average otherwise.  The
residual and tangent are respectively

```
L_i(u) = sum_d [theta_i DW_i.col(d) - theta_(i-e_d) DW_(i-e_d).col(d)] / dx_d
A_i(v) = sum_d [theta_i (DDW_i Gpair(v)_i).col(d)
              - theta_(i-e_d) (DDW_(i-e_d) Gpair(v)_(i-e_d)).col(d)] / dx_d
```

This gives the usual one-cell second-difference in the scalar uniform limit
and removes the collocated-central checkerboard nullspace.  It is the exact
discrete derivative of the chosen residual, including coefficient and raw-psi
variation, without a product-rule surrogate.

This changes a traction/Neumann boundary's tangential derivative from central
to the paired one-sided gradient.  It is required for a one-field, exact
residual/Jacobian formulation; retaining the old boundary gradient would make
the shared boundary flux have two incompatible derivatives.  Stencil selection
must use the global nodal domain, never an MFIter tile edge.

For an interior component `p`, the exact diagonal is the local impulse result:
the current forward gradient has `E0(p,d)=-1/dx_d`, each preceding flux has
`Ed(p,d)=+1/dx_d`, and their weighted `DDW` fluxes are differenced as above.
The boundary diagonal remains the BC evaluation of that same unit impulse.

The legacy `elasticop.small` mask offset changes from an implicit default to an
explicit opt-in regularization with default zero; if supplied, it must weight
residual and Jacobian identically.  A raw exact-zero psi can produce
zero-support rows and disconnected rigid modes; handling those inactive DOFs
or a nullspace is a separate Tier-3 formulation.  This correction must instead
fail fast with a clear message for that unsupported masked case.  The physical
no-floor oracle remains `elastic.use_psi=0` with a finite soft `model_void`.

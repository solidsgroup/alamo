# Phi radial-stress aberration result

## Outcome

Neither proposed solution removes the visible band.

- **Post-regrid consistency:** solving mechanics on the final hierarchy is a
  real correctness requirement for derived stress output, as the historical
  stale-versus-reconstructed mismatch demonstrates. In the valid frozen A/B,
  however, the radial oscillatory metric changes by at most `0.0227%`, far
  below the plan's `80%` selection threshold.
- **Native-face outcome:** the `phi`-following radial transition is already
  present in native coordinate-face normal traction. It is smooth and
  monotone, not a checkerboard ripple.
- **Visualization outcome:** synchronized saved `stress_*`, the exact
  face-average reconstruction, and the VisIt-equivalent radial projection
  agree to at worst `1.88e-5 Pa` in the finest-level target band. Plotting a
  differently centered face quantity changes the curve by up to `0.145 MPa`
  across this steep gradient but does not eliminate it.
- **Visible-aberration outcome:** the color feature is a smooth
  `1.4-1.8 MPa` constitutive radial-stress transition across the diffuse
  material field `phi`. The actual native alternating component is only
  `49-90 Pa`, with zero monotone overshoot.

The target band has `eta = 1` exactly, so `P + p eta I` cannot remove this
feature; it supplies only a constant diagonal offset there.

## Decision

The plan's decision rule selects neither solution as the visual fix.
Accordingly:

- no elastic operator, constitutive model, material interpolation, regrid
  lifecycle, or production output code was changed;
- no temporary restart hook was restored;
- only task-owned analysis scripts and reports were added.

The next diagnostic should target `phi` itself: compare the current
bitmap/diffuse material interpolation with an analytic circular `phi`, then
run a width-versus-resolution convergence study. That is the first route
likely to change this band in the simulation rather than merely rename or
recenter the plotted stress. It was deliberately excluded from this task's
operator-preserving scope.

## Evidence

- Existing pair (invalid for selection):
  [`EXISTING_PAIR.md`](EXISTING_PAIR.md)
- Controlled same-state solve:
  [`POST_REGRID_AB.md`](POST_REGRID_AB.md)
- Native-face classification:
  [`NATIVE_FACE.md`](NATIVE_FACE.md)
- Adversarial review:
  [`REVIEW.md`](REVIEW.md)
- Controlled comparison plot:
  [`controlled-comparison/comparison.png`](controlled-comparison/comparison.png)
- Native lineouts:
  [`native-face/profiles.png`](native-face/profiles.png)

## Verification

- analysis scripts compile with `python3 -m py_compile`;
- manufactured oracle passes and recovers both injected modes within
  floating-point accuracy;
- controlled comparison passes `--require-state-match`;
- all controlled non-mechanics and hierarchy hashes match;
- `bin/test-2d-g++` reports zero failed tests;
- `git diff --check` passes;
- repository device-lint, golden-compare, and A100-sanitizer gates pass.

The task-start and task-end hashes in
[`SOURCE_HASHES.sha256`](SOURCE_HASHES.sha256) are identical. Existing dirty
source changes were preserved without modification.

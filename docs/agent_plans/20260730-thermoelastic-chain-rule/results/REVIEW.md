# Adversarial review

## Verdict

PASS. The fresh read-only reviewer found no material issue in the exact
candidate diff and considered the known stale `rod_and_tube` strict reference
non-blocking for this scoped repair.

## Findings

- For `Fe = F F0^-1`, the corrected stress is exactly
  `Pe F0^-T`.
- The corrected tangent applies one `F0^-1` factor to each deformation
  gradient slot. Major symmetry is preserved by the transformation and by
  the existing `Sym::Major` storage.
- The edited model and tests are valid in both 2-D and 3-D host/device
  compilation paths.
- `Random()` now returns the derived model and constructs `F0 = I + 0.1 R`;
  the perturbation is safely nonsingular for the bounded random matrix used.
- The 3-D field-name and copy layouts now contain every component in the
  expected row-major order.
- The deterministic nonidentity-`F0` derivative tests necessarily detect
  the original missing chain rule. The free-expansion and field-layout tests
  cover the remaining repaired behavior.
- Adding `InputScraper.cpp` and `OutputLog.cpp` is the minimal correct link
  closure for direct references introduced by the branch re-merge; it changes
  no runtime policy.
- No test threshold, golden tolerance, or unrelated behavior was weakened.

## Known external reference

The reviewer confirmed that the `rod_and_tube_step2` traction mismatch does
not implicate this patch: the case uses `F0 = I`, where the repaired
constitutive expressions reduce to the original path. The mismatch was
already documented as a stale GPU reference on 2026-07-26.

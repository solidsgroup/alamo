# Eta-interface ripple root-cause result

## Verdict

The alleged 131-140 kPa eta ripple is not an oscillatory solver error. It is
the steep, resolved elastic-stress transition required to balance the diffuse
pressure load, amplified by using maximum second difference as the primary
roughness metric.

For the production equation

```text
div(P) = -p grad(eta),
```

constant pressure gives the exact discrete balance

```text
div(P + p eta I) = 0.
```

`P` is the constitutive elastic stress and must change by approximately `p`
through the eta band. The pressure-corrected total traction
`Q = P + p eta I`, not raw `P`, is the appropriate quantity for deciding
whether the pressure/interface equilibrium contains an alternating mode.

This satisfies the plan's affirmative H0 criterion. No equilibrium or operator
source change is warranted.

## Evidence

### Frozen production state

- The lagged pressure recovered from the serialized RHS is
  4,316,386.252462386 Pa.
- The production `CellGradientOnNode` is the divergence of the audited eta
  face flux to `5.12e-13 1/m`.
- Two frozen solves reproduce identical Newton histories and bitwise-identical
  cell and node plotfiles.
- The old max-second-difference measure is 146.5-146.9 kPa in raw traction but
  remains 13.2-24.2 kPa after pressure correction, demonstrating that this
  measure labels smooth steep curvature as roughness.
- The locked Nyquist projection is 18.5-19.1 kPa in raw traction and only
  0.548-0.902 kPa after correction, a 95.3-97.1% reduction.
- The accepted composite residual is 2.32284 MPa/m. The final cell eta differs
  slightly from the eta used to assemble the lagged RHS
  (`4.58e-4` relative Linf RHS fit), which bounds how exact an offline
  correction using the later cell plot can be.

### Aligned homogeneous analytical control

The one-level, one-box tanh control uses uniform material, `F0=I`, fixed eta,
and the production pressure source.

| Quantity | x-normal | y-normal |
|---|---:|---:|
| fitted raw slope / expected `-p` relative error | `2.48e-14` | `1.45e-14` |
| corrected traction Linf from a constant | `1.16e-7 Pa` | `9.69e-8 Pa` |
| corrected Nyquist projection | `0 Pa` | `9.24e-9 Pa` |
| elastic jump minus integrated RHS | `1.10e-7 Pa` | `3.73e-8 Pa` |
| post-solve residual Linf | `3.23e-4 Pa/m` | `3.13e-4 Pa/m` |

The x and y raw profiles are transposes to roundoff. The affine and smooth
manufactured oracle detects an injected 1,234.5 Pa alternating mode as
1,222.6 Pa while reporting `2.79e-11 Pa` on the clean field.

### Pressure and grid-phase controls

At pressure factors 0, 0.1, 0.5, 1, and 2, the aligned raw Nyquist amplitude is
respectively 0, 1.678, 8.388, 16.776, and 33.551 kPa: exactly linear in
pressure. The corrected amplitude remains below `6e-9 Pa`.

Normal shifts of 0, 0.25, 0.5, and 0.75 cells leave the aligned corrected
Nyquist projection at roundoff in both x and y. Thus the claimed feature is not
switched by bitmap/grid phase.

At 45 degrees, where a normal traction must be reconstructed by co-locating
adjacent native face columns, the corrected Nyquist amplitude is
2.071-2.107 kPa across the same shifts. Its less-than-2% phase variation and
first-order refinement (orders 1.06 and 1.02) identify a small
orientation/co-location truncation term, not the raw pressure transition.

### Signed-distance radial refinement

The analytic quarter-circle control has exact symmetry, homogeneous material,
and a physical 10-90 eta width of 4, 8, and 16 cells on the 64, 128, and 256
grids.

| Grid | raw Nyquist | corrected Nyquist |
|---:|---:|---:|
| 64 | 30.661 kPa | 365.97 Pa |
| 128 | 16.689 kPa | 85.62 Pa |
| 256 | 6.344 kPa | 17.69 Pa |

The corrected mode converges at orders 2.10 and 2.27. Holding the eta width at
eight cells while refining gives 312.87, 85.62, and 19.92 Pa, also
approximately second order. A small curved-interface truncation term exists,
but it is not the 131-140 kPa feature and decreases normally with resolution.

## Hypothesis disposition

- **H0 selected:** raw elastic stress contains the expected diffuse pressure
  transition; corrected total traction has no material alternating mode in the
  aligned control and converges normally for curved/oblique controls.
- **H1 rejected for the observed feature:** source and pressure face flux are
  algebraically identical on audited rows, and the aligned solve closes both
  local and global balance to roundoff.
- **H2 rejected as the primary cause:** quarter-cell shifts do not switch the
  aligned mode; the diagonal residual is nearly phase-independent.
- **H3 rejected as the primary cause:** a completely uniform material
  reproduces the raw transition and removes it exactly under pressure
  correction.
- **H4 applies only to the small corrected oblique/curved remainder:** it
  converges at first order for the declared co-location reconstruction and
  second order for native radial axis tractions.
- **H5/H6 factorial work was not run:** the plan explicitly skips Steps 6 and 7
  once H0 meets the affirmative criterion. The known post-solve-regrid output
  inconsistency remains a separate defect; the deterministic frozen solve and
  residual checks provide no indication that under-convergence generates the
  raw transition.

## Remedy

1. Keep the conservative mechanics equation unchanged.
2. Plot or compare `P + p eta I` when diagnosing diffuse pressure-interface
   equilibrium. Continue to expose raw `P`, but label it as constitutive
   elastic stress rather than total pressure-balanced traction.
3. Use the locked smooth-detrended Nyquist/even-odd metric as the ripple
   discriminator. Retain maximum second difference only as a historical
   steepness measure.
4. Continue to solve mechanics after any regrid before writing stress output.
5. For quantitative curved-interface work, prefer native staggered columns;
   explicitly label any normal-traction co-location and use the refinement
   behavior above as its error estimate.

The remedy is implemented in the task analysis and plots. No production solver
change was retained. The experimental frozen-restart hook was reverted and the
post-revert source was rebuilt; a fresh control again gave corrected Linf
`3.73e-8 Pa` and corrected Nyquist `5.49e-9 Pa`.

## Verification

Passed:

```bash
python3 docs/agent_plans/20260722-eta-ripple-rootcause/analyze_eta_ripple.py --self-test
python3 docs/agent_plans/20260722-eta-ripple-rootcause/analyze_eta_ripple.py --check-pressure-flux-identity
python3 docs/agent_plans/20260722-eta-ripple-rootcause/analyze_eta_ripple.py \
  --manifest docs/agent_plans/20260722-eta-ripple-rootcause/results/manifest.json \
  --check-complete
python3 docs/agent_plans/20260722-eta-ripple-rootcause/summarize_controls.py
make -j4
bin/test-2d-g++
```

The manifest contains 29 valid cases with commands, hashes, hierarchy, grid
spacing, eta width, pressure, and MPI size. Three stopped/invalid exploratory
attempts are retained and explicitly excluded.

## Artifacts

- `control_summary.json`, `control_summary.csv`, and `control_summary.png`:
  aggregate pressure, phase, and refinement evidence.
- `frozen-final-analysis/metrics.json`: frozen production metrics.
- `planar-analysis/*/metrics.json`: aligned and 45-degree controls.
- `radial-analysis/*/metrics.json`: signed-distance refinement controls.
- `manifest.json`: reproducibility inventory.
- `analyze_eta_ripple.py`, `analyze_planar.py`, `run_planar.sh`,
  `run_radial.sh`, and `summarize_controls.py`: executable analysis and decks.

## Limitations

- The production BMP double-circle was not rerun through a matched analytic
  double-circle geometry sweep. The uniform aligned and analytic radial
  controls are sufficient for the H0 criterion, but they do not quantify the
  sub-kilopascal contribution of the exact production bitmap.
- The frozen analysis corrects with the final plotted eta, while the serialized
  RHS was assembled from the slightly earlier lagged eta.
- No claim is made that a reconstructed 45-degree normal traction is a native
  staggered operator quantity.

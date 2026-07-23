# TASK: eta-ripple-rootcause
# Folder: docs/agent_plans/20260722-eta-ripple-rootcause/

---

## Header

| Field | Value |
|---|---|
| Risk tier | 3 (the investigation may change the mechanics load/operator coupling) |
| Model | Codex root session; keep root-cause reasoning in the main agent |
| Verification | partial-oracle |
| Est. scope | New diagnostic/input scripts and results in this folder; production changes only after a mechanism is demonstrated |
| Parallel-safe | no: the live tree is already dirty and the relevant Flame, Elastic, and Newton files overlap other work |

## Execution status

- Steps 1-2: complete; the frozen rerun is deterministic and the pressure-flux
  identity holds to roundoff.
- Steps 3-5: complete; pressure, x/y/45-degree phase, and fixed-width
  planar/radial refinement controls satisfy the affirmative H0 criterion.
- Steps 6-7: skipped by the plan's H0 decision rule.
- Step 8: complete as an analysis/documentation remedy; no solver change is
  warranted. The temporary restart hook was reverted.
- Final evidence: `results/RESULT.md`.

## Operating rules

1. Preserve the current dirty worktree. Put new decks, scripts, logs, and plots
   under this task folder; do not overwrite the accepted t=1 s artifacts.
2. Use a frozen static-mechanics case for localization. Hold `eta`, `phi`,
   casing support, temperature/F0, pressure, material properties, hierarchy,
   and boundary conditions fixed unless that field is the single variable under
   test.
3. Do not advance burn, thermal, or ballistics during a localization run.
   Replace `traction_from_chamber=1` with the recorded scalar pressure.
4. Do not edit production source until Steps 1-7 identify a mechanism with the
   affirmative-location criterion below. Temporary source A/Bs must be behind
   an explicit experimental knob and reverted after their result is recorded.
5. Change one factor at a time, retain paired input decks, and write a manifest
   containing git hash, executable hash, full command, input hash, MPI size,
   hierarchy, `dx`, measured eta width, and pressure for every run.
6. Treat an eta transition narrower than four finest-grid cells as
   under-resolved. Prefer six to eight cells for stencil comparisons.
7. Do not spend more runs on the legacy product-rule operator until a separate
   task gives it a convergent configuration for this high-contrast case.
8. Every step's CHECK must pass before continuing. At each checkpoint, show the
   experiment table, plots, metric deltas, and source diff (if any) to the user.

## Context budget

Read first: `benchmark/status.sh` output, this `PLAN.md`

Read:

- `docs/agent_plans/20260722-free-circular-casing/results/RIPPLE_DIAGNOSTICS.md`
- `docs/agent_plans/20260722-free-circular-casing/check_stress_ripple.py`
- `docs/agent_plans/20260722-free-circular-casing/results/ripple_no_postsolve_regrid/metrics.json`
- `docs/agent_plans/20260722-free-circular-casing/ripple_no_postsolve_regrid.log`
- `output_ripple_no_postsolve_regrid/05002cell/*`
- `output_ripple_no_postsolve_regrid/05002node/*`
- `output_ripple_no_postsolve_regrid/metadata`
- `input_rt1s_ideal`
- `src/Integrator/Flame.cpp:170-210,250-350,454-810`
- `src/Integrator/Flame.H:100-205`
- `src/Integrator/Base/Mechanics.H:60-150,188-350`
- `src/Integrator/Integrator.cpp:30-45,480-820,1170-1225`
- `src/Integrator/Integrator.H:85-115,450-470`
- `src/Numeric/Stencil.H:748-815,850-925`
- `src/Operator/Elastic.H:90-165`
- `src/Operator/Elastic.cpp:190-330,390-490,620-650`
- `src/Solver/Nonlocal/Newton.H:700-915,1070-1100,1438-1512`
- `src/IC/PSRead.H`

Reference only when a step names it: the new files and results created in this
task folder, focused mechanics unit tests, and the exact source lines changed by
an experimental or production patch.

Forbidden: `docs/archive/*`, unrelated task folders, and the out-of-scope
propellant sweep campaign.

## Objective

Determine whether the remaining approximately 131-140 kPa eta-interface
feature is (a) the expected discrete/continuum response to the diffuse pressure
load, (b) a sampling or diagnostic misclassification of a steep but smooth
profile, or (c) an error from load staggering, eta-dependent material blending,
grid phase/resolution, AMR coupling, or solver under-convergence. The result must
identify a minimal case in which one controlled factor predicts and switches
the unwanted feature on or off, then map that mechanism to the smallest safe
remedy.

The already-established post-solve-regrid inconsistency and nodal-output
amplification remain separate defects. All runs here must solve on their final
hierarchy and assess native face tractions; otherwise the eta investigation is
invalid.

## Competing hypotheses

| ID | Hypothesis | Affirmative signature |
|---|---|---|
| H0 | No numerical ripple remains; the present second-difference metric labels the resolved pressure transition as roughness | Raw traction changes across eta, but the pressure-corrected total traction and smooth-fit residual have no alternating mode and converge normally |
| H1 | Pressure source and elastic face flux become incompatible at an AMR, box, or boundary treatment | The uniform-interior identity passes but fails at actual composite/AMR or boundary rows, or repairing those exact rows removes at least 80% of the phase-aware ripple without changing the physical load |
| H2 | Bitmap/profile sampling or Cartesian grid phase produces aliasing | The amplitude is periodic under quarter-cell shifts or depends strongly on interface angle; matched analytic signed-distance eta removes it |
| H3 | The eta-dependent propellant/void material blend produces the feature | It remains under zero pressure with an imposed mechanical load, follows contrast, or disappears with a uniform model while the same RHS is retained |
| H4 | The feature is an under-resolved diffuse-interface truncation error | It decays at a consistent order when physical width is fixed and `h` is reduced, and is controlled primarily by width in cells |
| H5 | AMR fill/prolongation/reflux or box decomposition produces it | Same-finest-`dx` uniform and AMR solutions disagree, with error localized to box or coarse/fine boundaries |
| H6 | Linear/nonlinear solve under-convergence produces it | The feature changes materially when the accepted composite residual falls by two decades while all fields and discretizations are fixed |

## Oracle and metric contract

The first step creates
`docs/agent_plans/20260722-eta-ripple-rootcause/analyze_eta_ripple.py` and a
machine-readable run manifest. Its self-test is the executable oracle for all
later comparisons.

Commands after Step 1:

```bash
python3 docs/agent_plans/20260722-eta-ripple-rootcause/analyze_eta_ripple.py --self-test
python3 docs/agent_plans/20260722-eta-ripple-rootcause/analyze_eta_ripple.py \
  --manifest docs/agent_plans/20260722-eta-ripple-rootcause/results/manifest.json \
  --check-complete
```

For every interface, use the coordinate-face eta values implied by the
production `CellGradientOnNode` stencil. Lock the sampling coordinates, fit
window in physical units, fit family/order, minimum sample count, and pass
threshold in the self-test before inspecting production results. Set the pass
threshold from affine and smooth manufactured controls, not from the circular
case. Record over `0.1 < eta_face < 0.9`:

1. Each native staggered traction column: `P*e_x` on x faces and `P*e_y` on y
   faces. The aligned planar cases may use the exact native quantities `Pxx`
   and `Pyy`. Any curved-interface normal/tangential projection requires both
   columns at a common location and must be explicitly named, defined, and
   validated as a reconstruction rather than called native traction.
2. Pressure-corrected total traction. In the planar homogeneous control the
   production equation implies a constant quantity of the form
   `P_face + p*eta_face*I`; derive and unit-test the exact sign and indexing
   rather than assuming them from a plot.
3. A fixed checkerboard/Nyquist projection or even/odd pair-difference metric,
   plus the RMS and Linf after subtracting the predeclared smooth fit. Keep max
   second difference only as a secondary historical metric because it is
   phase-sensitive and is large for any legitimate steep transition.
4. The interface composite residual normalized by `p/w_10-90`, along with the
   absolute composite norm reported by the in-code solve. A per-level offline
   residual is diagnostic only and must not be presented as the AMR composite
   solver norm.
5. Integrated RHS, elastic physical-boundary reaction, and their force-balance
   mismatch as secondary global checks. A closed diffuse pressure interface can
   have nearly zero net RHS even when local balance is wrong, so these are not
   H0/H1 decision metrics. Use elastic `P` at a physical boundary unless eta is
   proven to vanish there.
6. Ripple divided by pressure and by the local stress scale; x/y transpose and
   90-degree rotation errors only between cases with matched subcell phase.
7. Identical polar/normal bins and sample counts across paired cases, with axis
   profiles retained only as smoke tests.

The affirmative-location criterion is met only if either:

- the aligned homogeneous planar case proves the exact corrected-flux identity,
  the alleged ripple disappears under the locked alternating-mode metric, and
  the phase/refinement behavior confirms that result; or
- one factor changes the normalized unwanted amplitude by at least 5x (and at
  least 80%), while load, width, grid phase, residual quality, and all other
  fields are controlled, and the predicted scaling repeats at a second
  resolution or pressure.

## Steps

### Step 1 - Freeze a reproducible static baseline and replace the roughness oracle

VERIFY:

```bash
benchmark/status.sh
test -f docs/agent_plans/20260722-free-circular-casing/results/ripple_no_postsolve_regrid/metrics.json
test -f input_rt1s_ideal
```

DO:

- Create either (a) a restart-cell plus restart-node deck that restores the
  synchronized final hierarchy and disables every subsequent phase, thermal,
  ballistic, and regrid advance, or (b) a narrow test harness that loads the
  frozen fields, solves mechanics, and writes immediately. Demonstrate that no
  state mutation occurs between the solve and output; a normal evolve cycle is
  not accepted on assertion alone.
- Freeze the saved `eta`, `phi`, casing support, temperature/F0, and hierarchy.
  Recover the scalar pressure that generated the saved RHS by fitting
  `rhs=-p*CellGradientOnNode(eta)` on nonzero-gradient interior rows and
  cross-checking the original solve log. Do not substitute the final thermo
  pressure because the mechanics coupling is explicitly lagged.
- Start from zero displacement and perform one static mechanics solve on the
  fixed hierarchy. Set `elastic.print_residual=1` to invoke the post-solve
  `compResidual`/reflux path and save `res`; retain `nr_diagnostics` separately
  for iteration history. Define the composite norm, fine-covered-row handling,
  and physical-BC-row handling in the analysis.
- Run the identical deck twice. Require bitwise-identical solution fields on
  the same platform, or record and justify a deterministic numeric tolerance.
- Implement the metric contract above, including plots versus normal distance
  of eta, raw traction, corrected total traction, smooth component, high-pass
  residual, and equilibrium residual.
- Recompute the old 131-140 kPa metric so the new measurements remain traceable
  to the diagnostic report.

CHECK:

```bash
python3 docs/agent_plans/20260722-eta-ripple-rootcause/analyze_eta_ripple.py --self-test
python3 docs/agent_plans/20260722-eta-ripple-rootcause/analyze_eta_ripple.py --case frozen-baseline --check-repeatability
```

STOP if saved stress does not match native-face reconstruction, the repeat is
not deterministic within the declared tolerance, or any regrid occurs after
the solve.

### Step 2 - Audit the discrete pressure balance before changing the solver

VERIFY: Step 1 passes and the baseline raw-face feature is reproduced.

DO:

- Algebraically rewrite every component of the 2-D production
  `CellGradientOnNode(eta)` as the divergence of explicit eta values on the
  same nodal-edge faces used by the conservative elastic flux. For example,
  the x component is the x difference of transverse two-cell averages. Encode
  this identity as a roundoff-level manufactured-field test.
- Form the corresponding pressure flux on every face and add it to the native
  elastic traction with the sign established by `div(P_face)-rhs`.
- Check whether the visible line and alternating residual are in raw elastic
  traction, corrected total traction, or both.
- Restrict this identity to the production case's spatially constant pressure.
  For varying pressure, `div(p*eta*I)` contains an additional `eta*grad(p)` term
  and is not the current equation; do not use that different problem as an H1
  control.

CHECK:

```bash
python3 docs/agent_plans/20260722-eta-ripple-rootcause/analyze_eta_ripple.py --check-pressure-flux-identity
```

Decision:

- If the identity holds and corrected traction is smooth, provisionally select
  H0: the raw elastic stress transition is the field required to balance
  `-p grad(eta)`. Confirm with the aligned planar and phase/refinement checks in
  Steps 3-5 before skipping to Step 8.
- If the identity holds but corrected traction retains a high-frequency mode,
  continue; the source is conservative but another coupling may be responsible.
- If the identity fails on production indexing or at AMR rows, prioritize H1
  and retain the exact failing rows for Steps 3 and 7.

### Step 3 - Build the smallest planar reproduction

VERIFY: Step 2 classifies the pressure-corrected baseline.

DO:

- Create a one-level, one-box, static 2-D case with homogeneous elastic
  material, constant temperature/F0, constant pressure, and analytic tanh eta.
- Begin with a planar interface normal to x, then repeat normal to y and at 45
  degrees. Choose boundary conditions with a known one-dimensional equilibrium
  so corrected total normal traction is constant in the aligned cases.
- With `F0=I` and a constant model, first run a sufficiently small pressure
  series and verify the expected small-strain linear scaling. Then run factors
  0, 0.5, 1, and 2 around the production pressure. Departure from linearity at
  production load is a constitutive/nonlinear result, not by itself evidence
  against H1 or H4.
- Only if Step 2 found a source/flux mismatch, add a temporary matched-face-load
  arm behind a test knob and compare it with the production source on identical
  fields. Do not substitute a superficially different stencil that is
  algebraically identical to the current one.

CHECK: the planar analytical balance, global reaction, and second-order affine
manufactured controls pass; the experiment either reproduces the high-pass
feature or records a clean negative result.

Decision: if the homogeneous planar case reproduces the feature, continue with
source/phase/resolution tests. If it does not, prioritize geometry and material
coupling in Steps 4 and 6.

### Step 4 - Separate eta representation, geometry, orientation, and grid phase

VERIFY: all comparisons use the same physical width and frozen pressure.

DO:

- Translate the planar zero contour by 0, 0.25, 0.5, and 0.75 `h`; repeat for
  x, y, and 45-degree normals.
- Compare the frozen production eta with an analytic signed-distance eta that
  reproduces its actual inner and outer/double-circle zero contours, physical
  width, and subcell phase. Do not use a single circle that changes the grain.
- Where possible, include an equilibrium phase-field tanh profile and the BMP
  profile sampled from the same signed distance. Compare eta itself, its normal
  derivative, and the pressure-corrected traction.

CHECK: report amplitude versus shift and angle, with matched phase for every
x/y rotation comparison.

Decision: periodic phase/angle sensitivity affirms H2. Bitmap-only sensitivity
maps to an analytic/SDF or improved sampling remedy; sensitivity shared by all
profiles maps to the discretization/width study.

### Step 5 - Run two non-confounded width/refinement studies

VERIFY: every case has constant pressure, static fields, at least four cells
across the measured 10-90% width, and the same grid phase.

DO:

- Study A: hold physical eta width fixed and reduce `h` at least twice. Measure
  observed order for corrected-traction high-pass RMS/Linf and composite
  residual.
- Study B: hold width in cells fixed while scaling physical width with `h`.
- Repeat the decisive planar case and one matched signed-distance radial case.
  Keep the outer casing far from the eta band in the planar case.

CHECK: produce convergence tables and log-log plots. A single refinement or a
change in both width and pressure is not accepted as convergence evidence.

Decision: consistent decay in Study A affirms H4 and supplies a mesh/width
requirement for the same continuum problem. Study B characterizes sensitivity
to width in cells, but cannot by itself prove a numerical artifact because its
physical width, and therefore its continuum problem, also changes.

### Step 6 - Factor eta's load, material, and nonlinear roles

VERIFY: use the minimal case from Steps 3-5 and preserve its width, phase, and
solver tolerances.

DO the following factorial controls:

| Arm | Pressure RHS | Eta-dependent material | Operator | Purpose |
|---|---|---|---|---|
| A | on | uniform | current conservative face operator | load only |
| B | off | on | current conservative face operator | coefficient interface under imposed boundary load |
| C | on | on | current conservative face operator | current unmasked coupling |

- For B, impose a nonzero mechanical boundary load; zero RHS plus zero load is
  a null experiment and cannot test material blending.
- Sweep propellant/void modulus contrast geometrically from 1 to the production
  value.
- Do not include `use_psi=1` as a mask-only arm. In Flame it also switches the
  nonlinear residual and stress semantics to the legacy product-rule operator,
  so it is not comparable to the native conservative-face baseline. If later
  evidence requires that operator, give it a separate plan, solvable setup,
  diagnostics, and acceptance criteria.
- Compare linearized/small-load behavior to the finite-strain result via the
  pressure factors from Step 3. Hold F0 fixed so thermal eigenstrain cannot
  masquerade as eta coupling.

CHECK: attribute amplitude to a main effect or interaction with repeatable
pressure/contrast scaling. An arm that fails to converge is inconclusive, not
evidence against its mechanism.

Decision: load-only sensitivity maps to H1/H4; coefficient or contrast
sensitivity maps to H3; nonlinear pressure scaling triggers a focused
constitutive/Newton follow-up before any discretization change.

### Step 7 - Test AMR and convergence only after the uniform mechanism is known

VERIFY: select one uniform-grid case that clearly exhibits the classified
feature and one that does not.

DO:

- Rerun the same uniform grid split into multiple boxes. Move box boundaries
  relative to eta without changing nodes or `dx`. Repeat identical fixed-grid
  fields with at least two MPI rank counts/distribution mappings. Set a numeric
  tolerance for cross-layout comparisons because reduction order may differ;
  require bitwise identity only for identical layouts.
- Compare one-level finest grid with a static AMR hierarchy having identical
  finest `dx` and a forced fine buffer wider than all source/operator stencils.
  Then move the coarse/fine boundary toward and across the eta band.
- Compare composite residual and native face traction after reflux/fill-patch;
  do not use a per-level offline norm as the acceptance oracle.
- Designate the post-solve `elastic.print_residual=1` `compResidual`/reflux norm
  as the H6 metric. State its norm, covered-cell policy, and physical-BC-row
  treatment. Independently tighten linear tolerance and the absolute nonlinear
  update tolerance, record iteration limits/counts, and require this post-solve
  composite residual to fall by two decades if possible.
  Current `nr_convergence` supports `update` and `psi_update`, not residual
  convergence. In the vector Newton path `solnorm` is not accumulated, so the
  logged update metric is effectively an absolute correction, despite the
  `relative` label. Do not claim a residual-mode or relative-update test unless
  a separately reviewed implementation adds it.

CHECK: spatially register the paired solutions and report differences at box,
coarse/fine, physical-boundary, and eta rows separately.

Decision: boundary-local differences affirm H5. A material ripple change after
a two-decade residual reduction affirms H6; otherwise stop spending runs on
tolerances.

### Step 8 - Implement only the remedy selected by the evidence

VERIFY: the affirmative-location criterion is satisfied, the mechanism is
restated in one sentence, and the user approves the proposed source change.

Remedy map:

| Finding | Remedy to prototype | Required regression |
|---|---|---|
| H0: expected diffuse pressure transition | Plot/compare pressure-corrected total traction where appropriate; document raw elastic-stress meaning; do not alter equilibrium | Discrete pressure-flux identity and analytical planar balance |
| H1: composite/AMR or boundary source-flux mismatch | Repair the demonstrated nonuniform-row treatment while preserving the exact uniform-interior face-divergence identity | Constant/affine eta, planar x/y/45°, uniform-versus-AMR balance, CPU/GPU agreement |
| H2: eta sampling/phase | Use matched analytic/SDF eta or an isotropic conservative sampling/prolongation path | Shift/rotation sweep and exact geometry preservation |
| H3: material/psi coupling | Make coefficient averaging/blending face-consistent or revise the mixture/mask formulation demonstrated by the factorial test | Constant-coefficient limit, contrast sweep, energy/equilibrium checks |
| H4: under-resolution | Enforce/document minimum cells across eta and tag a sufficient fine buffer; change stencil only if refinement is impractical | Fixed-width convergence with stated observed order |
| H5: AMR coupling | Correct fill/prolongation/reflux/source synchronization and always solve again after regrid | Uniform-versus-AMR composite comparison and post-regrid solve test |
| H6: under-convergence | Require both update control and a reviewed composite-residual gate | Newton unit tests plus unchanged converged solution under tighter tolerances |

After the minimal regression passes, rerun the frozen circular case, then the
full t=1 s case with a solve on the final hierarchy. Require:

- at least 5x reduction of the unwanted normalized high-pass feature, or H0's
  proof that only the expected corrected balance remains;
- no worse force balance or composite residual;
- x/y transpose symmetry at matched grid phase;
- saved output agreeing with native-face reconstruction;
- no post-solve regrid inconsistency; and
- focused 2-D unit/regression tests passing in debug and production builds.

Do not tune a production threshold to this one geometry. The selected remedy
must pass the planar orientation/phase case and the circular case.

## Checkpoints

- [x] After plan restatement: human confirms the frozen-baseline and
      pressure-corrected-traction approach before new experiment/source work.
- [x] After Step 2: execution continued under the user's explicit standing
      authorization; the feature was classified as expected balance stress.
- [x] After Steps 3-5: execution continued under the user's explicit standing
      authorization; H0 met the affirmative criterion, so factorial/AMR runs
      were skipped.
- [x] Before any experimental source edit: exact knob, expected discriminating
      result, and revert procedure are shown.
- [x] Before production implementation: not applicable; H0 selects no
      production equilibrium change.
- [ ] Before each commit: diff summary, focused oracle output, and proof that
      unrelated dirty files were not included.

## Adversarial review

After implementation and gates pass, use a fresh reviewer with no task context:

> Review the eta-ripple remedy commit on `chamber-gpu`. Assume the diagnosis is
> wrong. Check source/face indexing, sign convention, AMR composite balance,
> physical-boundary rows, cell/nodal ownership, GPU captures, hierarchy
> freshness, convergence semantics, and whether the tests distinguish a steep
> smooth interface from a checkerboard. Report findings only.

Record findings in `results/REVIEW.md` and adjudicate every finding before
commit/merge.

## Closeout

- [x] Root cause is affirmed or H0 is proven with the metric/identity evidence.
- [x] The selected minimal case switches the unwanted feature on/off and its
      scaling repeats at a second pressure or resolution.
- [x] Oracle and focused 2-D regression pass; `benchmark/status.sh` has no new
      failure attributable to this task.
- [x] `results/RESULT.md` records the mechanism, eliminated hypotheses, raw
      metrics, remedy, limitations, and exact reproduction commands.
- [x] Production and experimental diffs are separated; all temporary knobs not
      selected for production are reverted.
- [x] If the task is complete, touch `results/DONE`; append-only changelog and
      session-log rules are followed.

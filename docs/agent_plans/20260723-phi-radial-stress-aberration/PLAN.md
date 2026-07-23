# TASK: phi-radial-stress-aberration
# Folder: docs/agent_plans/20260723-phi-radial-stress-aberration/

---

## Header

| Field         | Value |
|---------------|-------|
| Risk tier     | 3 (mechanics/regrid lifecycle and solver-native stress output) |
| Model         | opus |
| Verification  | partial-oracle |
| Est. scope    | 2 analysis/run scripts first; conditionally 3-5 source/test files, approximately 150-300 lines |
| Parallel-safe | no: the dirty worktree already overlaps `Integrator`, `Flame`, and `Newton` |

## Operating rules

1. Read only the files listed in the Context budget. `docs/archive/` is
   forbidden.
2. Preserve every pre-existing worktree change. Before a source experiment,
   save the exact task-file diff and hashes in `results/baseline/`; never
   revert or commit unrelated changes.
3. Diagnose before retaining production code. Temporary hooks must default
   off, be captured as a patch, and be reverted unless selected by the decision
   table below.
4. Keep `phi`, `eta`, pressure, material parameters, tolerances, boundary
   conditions, hierarchy, and output time identical between paired arms. An
   arm with a state/hash mismatch is invalid.
5. Do not change the elastic operator, constitutive model, diffuse pressure
   source, or interface geometry in this task.
6. Analyze only the `phi` band where `eta` is locally constant. The previously
   resolved `eta` pressure transition is not the target.
7. Use native coordinate-face traction on the two symmetry axes as the primary
   oracle. A polar tensor reconstructed by co-locating staggered columns is
   secondary evidence and must be labeled as such.
8. Maximum second difference is not an oscillation oracle. Use a smooth
   detrend plus even/odd (Nyquist) amplitude, overshoot, exact face residual,
   and paired differences.
9. Every step's CHECK must pass before the next step. On failure, stop and
   record the invalid arm rather than tuning the run until it passes.

## Context budget

Read first: `benchmark/status.sh` output,
`docs/llm/PLAN.md`,
this `PLAN.md`

Read:

- `src/Integrator/Integrator.H:90-205,370-395`
- `src/Integrator/Integrator.cpp:470-530,873-920,1180-1230,1312-1395`
- `src/Integrator/Base/Mechanics.H:60-345`
- `src/Integrator/Flame.H:20-75,120-205`
- `src/Integrator/Flame.cpp:454-675,960-1015`
- `src/Solver/Nonlocal/Newton.H:90-345,365-390,600-690`
- `src/Test/Solver/Nonlocal/Newton.H`
- `src/test.cc`
- `Makefile`

Reference only when the named step calls for it:

- `input_rt1s_ideal`
- `output_ideal_rod_and_tube_free_circular_t1/05000node`
- `output_ripple_no_postsolve_regrid/05002cell`
- `output_ripple_no_postsolve_regrid/05002node`
- `docs/agent_plans/20260722-free-circular-casing/check_stress_ripple.py`
- `docs/agent_plans/20260722-free-circular-casing/results/RIPPLE_DIAGNOSTICS.md`
- `docs/agent_plans/20260722-free-circular-casing/results/ISSUE_INVESTIGATION_REPORT.md`
- `docs/agent_plans/20260722-face-consistent-output-stress/results/RESULT.md`
- `docs/agent_plans/20260722-eta-ripple-rootcause/run_frozen_baseline.sh`
- `docs/agent_plans/20260722-eta-ripple-rootcause/results/CHECKPOINT_BEFORE_TEST_HOOK.md`
- `docs/agent_plans/20260722-eta-ripple-rootcause/results/runs/frozen-final-a/diff.patch`

Forbidden: `docs/archive/*`, unrelated task folders, propellant sweep files

## Objective

Determine whether the radial-stress aberration following the black `phi`
contours is caused by (A) mechanics being solved before the final regrid or
(B) reconstructing and polar-projecting a nodal tensor instead of inspecting
the conservative native face traction. Quantify both contributions on the
same frozen state, retain only fixes selected by the paired evidence, and
produce a VisIt-usable output or documented native-face diagnostic that does
not misrepresent the `phi` interface.

The task must not claim that `P + p eta I` fixes this feature: `eta` is
constant in the target band, so that correction is only a constant offset
there.

## Existing evidence to reproduce, not assume

- The original final level-2 plot differed from an offline face reconstruction
  by approximately `48.6 kPa` in the outer-interface band.
- A solve on a final, non-regridded hierarchy reduced saved-versus-reconstructed
  stress to approximately `2.6e-5 Pa` and reduced the interior face residual
  from approximately `8.16e8` to `2.32e6 Pa/m`.
- The old axis-layer maximum-second-difference values changed very little in
  that comparison. This suggests lifecycle synchronization fixes a real
  consistency defect but may not remove the native `phi` feature.
- Face-consistent output previously reduced the original nodal `phi` mismatch
  by 90.2%, from `0.2108 MPa` to `0.0207 MPa`, but it still assembles
  directionally staggered face columns at nodes.

These observations are provisional for the current visual aberration until the
new oracle measures the same `phi` band and the same radial quantity shown in
VisIt.

## Oracle

Commands:

```bash
python3 docs/agent_plans/20260723-phi-radial-stress-aberration/analyze_phi_radial.py --self-test
python3 docs/agent_plans/20260723-phi-radial-stress-aberration/analyze_phi_radial.py \
  --node output_ideal_rod_and_tube_free_circular_t1/05000node \
  --out docs/agent_plans/20260723-phi-radial-stress-aberration/results/current
python3 docs/agent_plans/20260723-phi-radial-stress-aberration/analyze_phi_radial.py \
  --node output_ripple_no_postsolve_regrid/05002node \
  --out docs/agent_plans/20260723-phi-radial-stress-aberration/results/final-hierarchy
python3 docs/agent_plans/20260723-phi-radial-stress-aberration/compare_arms.py \
  docs/agent_plans/20260723-phi-radial-stress-aberration/results/current/metrics.json \
  docs/agent_plans/20260723-phi-radial-stress-aberration/results/final-hierarchy/metrics.json
```

After a production change is selected:

```bash
make -j4 bin/test-2d-g++ bin/alamo-2d-g++
bin/test-2d-g++
python3 docs/agent_plans/20260723-phi-radial-stress-aberration/analyze_phi_radial.py --self-test
benchmark/status.sh
```

Covers:

- exact reproduction of the production 2-D face stress and face divergence;
- the `phi=0.1-0.9` band with negligible local `eta` gradient;
- saved nodal/polar stress, face-averaged nodal stress, and native axis-face
  traction;
- bottom-axis `Pxx` versus right-axis `Pyy` transpose symmetry;
- detrended even/odd amplitude, monotone-envelope overshoot, and interface
  residual;
- state hashes, hierarchy identity, nonlinear history, and output freshness;
- unchanged non-mechanics fields and repository correctness gates.

Does not cover:

- whether the bitmap-defined `phi` geometry is physically desirable;
- a sharp-interface analytical solution at arbitrary angle;
- the first-order co-location error in an off-axis polar reconstruction;
- 3-D behavior beyond the repository gates.

## Experimental matrix and decision rule

Measure all valid cells of this matrix:

| Mechanics state | Saved nodal/polar `Prr` | Face-averaged nodal `Prr` | Native axis-face normal traction |
|---|---:|---:|---:|
| Existing output after the lifecycle under test | required | required | required |
| Same frozen fields after a final-hierarchy mechanics solve | required | required | required |

Interpretation:

| Result | Diagnosis | Retained action |
|---|---|---|
| Native-face aberration drops after the final-hierarchy solve | Post-regrid stale displacement/equilibrium contributes | Implement solution 1 |
| Native face is smooth but either nodal/polar field is not | Reconstruction/co-location contributes | Implement solution 2 |
| Both effects occur | Mixed lifecycle and output problem | Implement both, in separate commits |
| Native face remains aberrant after the synchronized solve | Neither proposed fix fully explains it | Keep only independently justified consistency/output fixes; stop and report that `phi` discretization/geometry is next |

A proposed fix is selected for the visual aberration only if the paired
detrended even/odd amplitude or overshoot falls by at least 80% without
changing the smooth bulk radial profile outside the `phi` band by more than
the accepted nonlinear-solve tolerance. A post-regrid freshness fix may still
be retained as a correctness fix if it restores face agreement and residual
closure but leaves the visual metric unchanged; the result must explicitly say
that it did not fix the aberration.

## Steps

### Step 1 - Lock the `phi`-band radial-stress oracle

VERIFY:

```bash
benchmark/status.sh
test -d output_ideal_rod_and_tube_free_circular_t1/05000node
test -d output_ripple_no_postsolve_regrid/05002node
git diff --check
```

DO:

- Create `analyze_phi_radial.py` in this task folder by extracting only the
  necessary plotfile reader and exact `FaceStress` formulas from the prior
  diagnostic.
- Resolve components by Header name; honor valid boxes and AMR coverage.
- Define the target band by `0.1 <= phi <= 0.9` and a locally flat-`eta`
  criterion. Report the actual `eta` range and gradient in the selected band.
- On the bottom and right symmetry axes, compare raw native `x`-face `Pxx` and
  `y`-face `Pyy`; do not co-locate tensor columns for the primary metric.
- Add secondary face-averaged nodal and full polar projections using the
  correct center `(0.0877, 0.0877)` and both `Pxy` and `Pyx`.
- Report smooth-detrended even/odd amplitude, overshoot, saved-versus-face
  mismatch, exact face residual, x/y transpose error, and hashes.
- Add a manufactured self-test that is smooth through a radial material jump,
  then injects a known alternating face mode and a known nodal-only mode. The
  oracle must distinguish them.

CHECK:

```bash
python3 docs/agent_plans/20260723-phi-radial-stress-aberration/analyze_phi_radial.py --self-test
python3 docs/agent_plans/20260723-phi-radial-stress-aberration/analyze_phi_radial.py \
  --node output_ideal_rod_and_tube_free_circular_t1/05000node \
  --out docs/agent_plans/20260723-phi-radial-stress-aberration/results/current
```

The self-test must recover each injected amplitude within 2% and report the
clean native-face mode below `1e-8` of the injected amplitude.

### Step 2 - Reanalyze the existing lifecycle pair

VERIFY:

- The analyzer passes Step 1.
- Record plot Header, grid, distribution, component, executable, input, and
  source-diff hashes for both existing artifacts.

DO:

- Analyze the original final output and the existing final-hierarchy output
  with the identical oracle.
- Create `compare_arms.py` and a four-panel plot showing, for both axes:
  saved/polar `Prr`, face-averaged nodal `Prr`, native face normal traction,
  `phi`, and local `eta`.
- Compare non-mechanics fields and hierarchy metadata. If they are not a valid
  pair, label the result historical-only and proceed to Step 3 without selecting
  a fix.
- Record which prior `48.6 kPa` discrepancy is output staleness and whether the
  native `phi`-band mode itself changes.

CHECK:

```bash
python3 docs/agent_plans/20260723-phi-radial-stress-aberration/analyze_phi_radial.py \
  --node output_ripple_no_postsolve_regrid/05002node \
  --out docs/agent_plans/20260723-phi-radial-stress-aberration/results/final-hierarchy
python3 docs/agent_plans/20260723-phi-radial-stress-aberration/compare_arms.py \
  docs/agent_plans/20260723-phi-radial-stress-aberration/results/current/metrics.json \
  docs/agent_plans/20260723-phi-radial-stress-aberration/results/final-hierarchy/metrics.json
```

Write `results/EXISTING_PAIR.md` with a valid/invalid verdict and no production
recommendation unless all pairing checks pass.

### Step 3 - Controlled post-regrid mechanics A/B

VERIFY:

- Human approves the Step 2 metrics and the exact experimental diff.
- The frozen restart restores `phi`, `eta`, temperature, model, RHS, pressure,
  hierarchy, and boundary-condition time. If model deserialization is not
  available in the current tree, restore the prior default-off diagnostic hook
  as an isolated task patch; do not silently rebuild the model from a later
  plotted state.

DO:

- Starting from one frozen checkpoint, produce:
  - Arm A: the restored pre-resolve state analyzed without another regrid;
  - Arm B: exactly one mechanics solve on that same final hierarchy, followed
    by stress derivation and output, with no intervening regrid or physics
    advance.
- Run serial first. Repeat the valid B arm once and require bitwise-identical
  fields and Newton history before considering MPI.
- Hash all frozen non-mechanics fields before and after. Record the solver
  residual and the exact plot-write ordering.
- Apply the Step 1 oracle to both arms.

CHECK:

```bash
python3 docs/agent_plans/20260723-phi-radial-stress-aberration/compare_arms.py \
  docs/agent_plans/20260723-phi-radial-stress-aberration/results/controlled-a/metrics.json \
  docs/agent_plans/20260723-phi-radial-stress-aberration/results/controlled-b/metrics.json \
  --require-state-match
```

The comparison must prove identical `phi`, `eta`, material, pressure, and
hierarchy. Otherwise the A/B is invalid. Write `results/POST_REGRID_AB.md` and
apply the decision table before any production edit.

### Step 4 - Isolate native-face versus nodal/polar output

VERIFY:

- Step 3 has a valid synchronized state.
- Exact offline face stress agrees with `Detail::FaceStress` on the
  manufactured/unit oracle.

DO:

- On both symmetry axes, compare the single native coordinate-face normal
  traction against:
  1. the low/high-face average stored at a node;
  2. the assembled nodal tensor projected to `Prr`;
  3. the VisIt-equivalent polar expression.
- For off-axis context only, co-locate the staggered columns explicitly and
  label the reconstruction and expected order. Do not use it to overrule the
  axis-native result.
- Determine whether the visible line is introduced by low/high face averaging,
  x/y column co-location, polar projection, or is already present on the native
  faces.
- Generate a task-owned CSV and plot that can be independently compared with a
  VisIt lineout.

CHECK:

```bash
python3 docs/agent_plans/20260723-phi-radial-stress-aberration/analyze_phi_radial.py \
  --node docs/agent_plans/20260723-phi-radial-stress-aberration/results/controlled-b/node \
  --out docs/agent_plans/20260723-phi-radial-stress-aberration/results/native-face \
  --write-visit-comparison
```

Write `results/NATIVE_FACE.md` with one of: native-smooth/nodal-rough,
native-rough, mixed, or inconclusive.

### Step 5 - Implement only the selected fix or fixes

VERIFY:

- Human approves `POST_REGRID_AB.md`, `NATIVE_FACE.md`, the selected branch of
  the decision table, and the proposed source-file diff.

DO, solution 1 if selected:

- Refactor the static mechanics solve into a callable operation without
  changing its solver parameters or equations.
- Track whether the hierarchy was regridded after the last valid mechanics
  solve.
- Before stress output, solve once on the final hierarchy only when stale;
  then derive stress from that displacement. Do not solve twice when no regrid
  occurred and do not permit another regrid between solve and write.
- Add an ordering regression test with a forced regrid immediately before an
  output step. Assert final-hierarchy residual closure and saved-versus-face
  agreement.

DO, solution 2 if selected:

- Preserve existing `stress_*` semantics.
- Add an opt-in, clearly named native-face traction diagnostic or the smallest
  face-derived radial output justified by Step 4. Do not call a co-located
  off-axis tensor “native.”
- Reuse `Detail::FaceStress`; do not duplicate constitutive or gradient
  formulas in production.
- Add a unit/manufactured test proving exact axis-normal traction, centering,
  symmetry-plane behavior, and no feedback into displacement or material
  advancement.

CHECK:

```bash
make -j4 bin/test-2d-g++ bin/alamo-2d-g++
bin/test-2d-g++
git diff --check
```

Each solution is a separate commit and must independently pass its relevant
oracle. If neither solution meets the selection rule, retain no production
source change and close with a diagnostic result.

### Step 6 - End-to-end confirmation

VERIFY:

- Selected focused tests and frozen A/B pass.
- The working diff contains only approved task changes.

DO:

- Run a short forced-regrid case that writes immediately after regrid.
- Run the one-second rod-and-tube confirmation only after the frozen and short
  oracles pass.
- Produce side-by-side VisIt-equivalent plots and radial lineouts for raw
  `stress_radial`, native face traction, `phi`, and `eta`. Do not use
  `P + p eta I` to judge the `phi` feature.
- Compare displacement, pressure, temperature, `phi`, `eta`, and bulk radial
  and hoop stress outside the target band.
- Run all repository gates.

CHECK:

```bash
make -j4
bin/test-2d-g++
benchmark/status.sh
```

Acceptance:

- no stale saved-versus-face mismatch after a regrid/output event;
- selected `phi`-band oscillatory metric reduced by at least 80%;
- final face residual is at the solver-accepted scale;
- no material change outside the `phi` band beyond nonlinear tolerance;
- device lint, golden compare, and A100 sanitizer pass.

If lifecycle consistency passes but the native face remains aberrant, report
solution 1 as a correctness fix only and do not claim the stress aberration is
fixed.

## Checkpoints

- [ ] After plan restatement: agent restates the paired matrix and the user
      confirms before any new script or source edit
- [ ] After Step 1: human reviews region selection, metrics, and manufactured
      oracle
- [ ] After Step 2: human reviews whether existing artifacts form a valid pair
- [ ] Before Step 3: exact temporary hook diff and frozen-state hash list
      reviewed
- [ ] After Step 3: human selects or rejects solution 1
- [ ] After Step 4: human selects or rejects solution 2
- [ ] Before each production commit: diff summary, focused oracle, and proof
      that no unrelated dirty changes are included
- [ ] Before the one-second run: short forced-regrid case passes
- [ ] Before closeout: full results and repository gates reviewed

## Adversarial review

After implementation is complete and gates pass, use a fresh reviewer with:

> Review the retained commits for
> `20260723-phi-radial-stress-aberration`. Assume they contain a defect.
> Check mechanics/regrid/output event ordering, duplicate solves, hierarchy
> freshness, pressure and eta time levels, face centering, symmetry parity,
> device-lambda captures, Eigen expression lifetimes, whether the diagnostic
> feeds back into the solution, and whether a steep smooth transition was
> mislabeled as an oscillation. Report findings only.

Store findings in `results/REVIEW.md`; adjudicate every finding before closeout.

## Closeout

- [x] Oracle passes; `benchmark/status.sh` is all green
- [x] `results/RESULT.md` states separately:
      post-regrid consistency outcome, native-face outcome, visible-aberration
      outcome, and any remaining `phi` geometry/discretization limitation
- [x] Experimental hooks are either reverted or covered by production tests
- [x] `results/DONE` created only when no required work remains
- [x] Append-only changelog/version update only if release-worthy
- [x] Session log line appended to `docs/llm/SESSION_LOG.tsv`

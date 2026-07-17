# Elastic regression recovery handoff

## Status and decision

**BLOCKED: do not commit, refresh references, loosen tolerances, raise iteration
caps, or add a coefficient/psi floor.**

The finite-soft-void focused oracle passes, but the full GCC regression run
does not. The task is not complete and the current source changes are
deliberately uncommitted. The next agent should treat the failures as a
solver-correctness investigation, not as a reference-update task.

| Item | Value |
|---|---|
| Worktree | /tmp/alamo-elastic-void-continue |
| Branch | elastic-void-wip-20260713 |
| HEAD/base | a06c6f15d — WIP paired-stencil rewrite |
| Full-run report | report/output_2026-07-13_15.29.01_kermit.json and .html |
| Full-run command | scripts/runtests.py --comp=g++ --no-backspace |
| Current source delta | src/Numeric/Stencil.H; src/Operator/Elastic.H; src/Operator/Elastic.cpp; src/Operator/Operator.cpp; src/Solver/Nonlocal/Newton.H |
| New focused test | tests/ElasticSoftVoid/input and tests/ElasticSoftVoid/test |
| Status gate | benchmark/status.sh: device-lint PASS; worktree is intentionally dirty |

Read FULL_SUITE_REVIEW_20260713.md before inspecting source. It is the
detailed evidence record; this document is the execution handoff.

## What did pass

The selected formulation is a finite, positive soft material in the void:

~~~
elastic.use_psi = 0
elastic.psi_floor = 0
model_void.kappa = model_void.mu = 0.2 MPa
~~~

It does not claim to solve raw exact-zero psi masking. Focused checks passed:

- ElasticSoftVoid 2-D serial: residual 8.97774e-06.
- ElasticSoftVoid 3-D serial: residual 9.44120e-06.
- A direct two-rank 2-D smoke test and two focused sensitivity screens.
- All four Unit sections in the full GCC run.

Those results only establish that this narrow formulation can complete. They
do not validate the revised operator across the existing mechanics suite.

## Full-suite outcome

All 144 sections executed: 131 run passes and 13 run failures. There are 85
passing checks, 17 failing numerical checks, and 42 sections without a check.
There were no skips. The 17 check failures are genuine numerical mismatches;
the similarly named diff.patch files are runner metadata and are irrelevant.

### Runtime failures

| Family / section | Failure signature |
|---|---|
| Eshelby / 3D-serial-4levels | Stops at 20 MLMG iterations, relative residual 2.097364563e-08. |
| Eshelby / 3D-parallel-5levels | Stops at 20 iterations, relative residual 4.256658658e-08. |
| EshelbyFiniteKinematics / 2D-serial | Stops at 20 iterations, relative residual 1.036873858e-08. |
| PlateHole / 2D-serial | Stops at 20 iterations, relative residual 2.023381522e-06. |
| PlateHole / 3D-serial | Stops at 20 iterations, relative residual 2.803395211e-06. |
| PlateHole / 3D-parallel | Stops at 20 iterations, relative residual 2.912295376e-06. |
| Scratch / 2D-serial | Stops at 40 iterations, relative residual 1.952953903e-08. |
| Scratch / 2D-parallel | Stops at 40 iterations, relative residual 1.952953682e-08. |
| FracturePFCZM / 2d-serial-notch | Diverges at iteration 404, relative residual 1.057125712e+20; Linear.H:144 abort. |
| RubberPlateHole / serial-2d | At t=0.5, Newton iteration 2 stalls at 0.241361769 after 150 iterations. |
| RubberPressurizedHole / 2d-serial | At t=2, Newton iteration 3 diverges at iteration 21 to 3.536749885e+20. |
| RubberPressurizedHole / 2d-parallel | Same event in MPI; diverges at iteration 21 to 4.119987487e+20. |
| TopOp / parallel | At step 58, diverges at iteration 88 to 1.273287521e+20; MPI_ABORT/error 6. |

### Numerical check failures

| Group | Sections | Evidence / priority |
|---|---|---|
| Critical hidden failure | SCPSpheresElastic / 2d-parallel-long | At t=0.006 (step 6000), a solve grows from roughly 9.6e-03 to 5.77e3, the next RHS/residual is O(1e23), and plot 10000 onward has NaN displacement, stress, and residual fields. The current checker does not reject it. Start here. |
| Mechanics drift | RubberWithInclusion / serial-2d; Suture / 2D-serial-4levels; ThermoElastic / 2d; TwinGrowth / serial | Completed runs differ from references: e.g. stress error 8.33e-03 versus 1e-04 in RubberWithInclusion; ThermoElastic displacement/stress about 0.038 versus 0.01 tolerance; Suture stress about 0.095553 versus 1e-04. |
| Broad elastic drift | VoronoiElastic / 2d-serial-restart, 2d-serial, 2d-parallel, periodic-2d-serial | Eta is nearly unchanged but displacement differs about 0.867–0.873 and stress about 0.790–0.823, each against 0.01 tolerance. |
| Periodic/AMR drift | VoronoiSimplePeriodic / 2d-amr0-mg0-np1, 2d-amr2-mg0-np1, 2d-amr0-mga-np1, 2d-amr2-mga-np1, 2d-amr0-mg0-np4, 2d-amr2-mg0-np4, 2d-amr0-mga-np4, 2d-amr2-mga-np4 | Model and displacement are relatively close; stress drift is about 0.0668–0.0674 versus 0.05 tolerance. |

Exact stdout/stderr are retained under each test's
output_2026-07-13_15.29.01_kermit_<section> directory. In particular:

~~~
tests/PlateHole/output_2026-07-13_15.29.01_kermit_2D-serial/stdout
tests/SCPSpheresElastic/output_2026-07-13_15.29.01_kermit_2d-parallel-long/stdout
~~~

Keep those artifacts until the cause is known.

## Attribution and baseline boundary

Do not attribute every failure to the latest unstaged edit without a same-base
comparison. The current HEAD was already a WIP paired-stencil rewrite.

- PlateHole 2-D is confirmed pre-existing relative to today's unstaged delta:
  it failed in this worktree before the first modified source timestamp and
  showed the same 20-iteration-cap class in the main checkout.
- The remaining failures, including Eshelby, Rubber, and Scratch, have not
  been compared against a clean same-base build and remain attribution
  unresolved.
- Historical June 29 reports in the main checkout passed many affected
  families, including PlateHole, Eshelby, Rubber, Scratch, Suture,
  ThermoElastic, TwinGrowth, and Voronoi. Thus these are not longstanding
  repository-baseline failures; the WIP stack needs isolation.

Relevant historical reports are:

~~~
/home/jackplum/Projects/alamo/report/output_2026-06-29_12.57.15_kermit.json
/home/jackplum/Projects/alamo/report/output_2026-06-29_13.16.44_kermit.json
/home/jackplum/Projects/alamo/report/output_2026-06-29_13.19.24_kermit.json
/home/jackplum/Projects/alamo/report/output_2026-06-29_13.38.11_kermit.json
~~~

## Source-review findings that require follow-up

1. **Public SetPsi halo contract (P1).** Elastic.cpp allocates a three-cell
   psi collar at lines 60–82. The public legacy overload
   Elastic::SetPsi(int, const MultiFab&) at lines 166–183 copies only the
   caller's grown tile and does not call FillPsiGhosts(). The Field overload
   at lines 187–201 does fill all levels, and the normal vector Newton path
   uses it at Newton.H lines 319–321. Fapply's first active ghost rows use a
   predecessor flux whose cell-to-node psi average can read the third low
   collar. Direct legacy callers can therefore see default 1.0 or stale
   collar data. No in-tree direct caller was found, but this public contract
   is unsafe and is untested.

2. **Masked cross-AMR coefficient restriction (P2).** With
   elastic.solver.average_down_coeffs=1, Elastic.cpp lines 843–994 restricts
   DDW, while psi/theta is filled/coarsened separately. Fapply independently
   multiplies theta and DDW (roughly lines 309–360). A masked C/F interface
   may therefore use an inconsistent effective coefficient. FracturePFCZM
   exercises masking and average-down and is a priority reproducer. Define
   whether theta*DDW, rather than only DDW, must be restricted together.

3. **Displacement FAB precondition (P3).** Numeric::PairedGradient and its
   bounded helpers do not protect an anchor that lies outside the supplied
   displacement Array4. Newton checks a two-ghost precondition at
   Newton.H lines 91–96 and 244–249; direct Apply callers do not have an
   equivalent runtime guard. CUDA was not rerun after the final source edit.

4. **No raw-zero-mask claim.** A zero diagonal is diagnosed rather than
   regularized. That remains a deliberate unsupported formulation requiring
   inactive-DOF/nullspace handling; do not turn the diagnostic into a hidden
   floor.

## Known focused-oracle gaps

ElasticSoftVoid was useful as a smoke/property test but is not an adequate
acceptance oracle:

- Its 2-D run overrides the void modulus to 0.2 MPa, yet its checker only
  requires a material contrast of at least ten; the base 140/10 MPa values
  already meet that threshold. It must assert the expected kappa/mu minima
  and contrast.
- It verifies residuals, finiteness, AMR existence, and nontrivial response,
  but not force balance, symmetry, interface traction, a manufactured
  displacement/stress, or a trusted physical reference.
- It has no automated MPI section; the prior two-rank smoke was manual.
- It does not prove a coarse/fine interface intersects the heterogeneous
  material, and its component checks permit any matching field instead of
  requiring all components.
- Existing unit coverage does not directly exercise PairedGradient,
  StencilWithinArrayBounds, bounded interpolation, impulse diagonal versus
  Fapply, or exact-zero rejection.

## Reproduction commands

Build both dimensions before reproducing, preserving the existing output:

~~~
./configure --dim=2 --comp=g++ && make -j4
./configure --dim=3 --comp=g++ && make -j4
~~~

Fast first reproducers:

~~~
scripts/runtests.py --comp=g++ --no-backspace --no-clean \
  --sections 2d-parallel-long tests/SCPSpheresElastic

scripts/runtests.py --comp=g++ --no-backspace --no-clean \
  --sections 2D-serial tests/PlateHole

scripts/runtests.py --comp=g++ --no-backspace --no-clean \
  --sections 2d-serial-notch tests/FracturePFCZM
~~~

Full validation:

~~~
scripts/runtests.py --comp=g++ --no-backspace --no-clean
~~~

The aggregate Make target was also attempted. Its tab check passed, but
documentation stopped because the local doxygen executable cannot load
libclang-18.so.18. That is an environment dependency; it does not explain the
mechanics test failures. When the environment is repaired, run make test as an
additional gate.

## Required investigation order

1. **Establish clean baselines without disturbing this worktree.** Create a
   separate worktree at the parent/pre-WIP baseline and at a clean copy of
   a06c6f15d. Rebuild and compare the first three reproducers above. Record
   the first commit/configuration that changes each signature.
2. **Stop the SCPSpheres NaN path.** Reduce it to the first failed solve near
   t=0.006. Inspect residual/Jacobian consistency, finite values before and
   after Fapply/Diagonal, and the fixed-iteration/nonlinear acceptance path.
   Add a test that rejects NaN plots and an unconverged solve.
3. **Diagnose the broad linear-solver change.** PlateHole is the quickest
   deterministic signal. Compare the paired operator, diagonal, boundary
   closure, omega, and smoother against a known-good implementation. Do not
   tune nriters or tolerance until an invariant explains the mismatch.
4. **Exercise masked AMR/MPI separately.** Reproduce FracturePFCZM and
   RubberPressurizedHole with/without psi and average_down_coeffs. Test
   whether restricting theta*DDW jointly fixes a demonstrable C/F
   inconsistency.
5. **Repair the public coefficient contract and add low-level coverage.**
   Decide whether the per-level SetPsi overload must fill/finalize all levels,
   be made private/deprecated, or guard Apply until callers do so. Add direct
   tests for halo use, stencil bounds, diagonal impulse agreement, and raw
   zero handling.
6. **Strengthen the finite-soft-void oracle only after correctness is
   restored.** Assert actual material values, add automated MPI and C/F
   heterogeneity coverage, and add a mechanics-sensitive oracle.

## Non-negotiable guardrails

- Do not use elastic.psi_floor, m_psi_small, an artificial stiff void, relaxed
  tolerances, larger iteration limits, forced line-search acceptance, or
  reference refresh to make failures disappear.
- Do not call a run passing because its process exits successfully if its
  plots have NaN or its physical checker fails.
- Keep residual, Fapply, Diagonal, smoother, and boundary closure derived
  from the same discrete operator.
- Keep direct and vector SetPsi behavior coherent across all AMR levels.
- Do not commit until the complete regression gate and a fresh adversarial
  review are clean. The user has not requested a commit.

## Acceptance gates for recovery

1. Every root cause is recorded with a minimal reproducer and an attribution
   boundary.
2. Targeted tests cover each repaired contract, including serial/MPI and
   masked AMR where relevant.
3. scripts/runtests.py --comp=g++ --no-backspace completes with all existing
   run and check sections passing.
4. ElasticSoftVoid remains a genuine finite-soft-void test and passes in
   2-D, 3-D, and its new automated MPI coverage.
5. make test passes once the local doxygen/libclang dependency is usable.
6. A fresh Tier-3 review finds no oracle weakening, hidden regularization,
   stale hierarchy data, or residual/tangent/diagonal mismatch.

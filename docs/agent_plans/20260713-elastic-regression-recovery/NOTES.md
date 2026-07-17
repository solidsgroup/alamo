# Elastic regression recovery notes

## 2026-07-13 plan restatement and guardrails

The user assigned the blocked elastic-void continuation.  Work proceeds from
clean same-base comparisons before any source edit.  The preserved dirty WIP
remains at `a06c6f15d` in `/tmp/alamo-elastic-void-continue`; comparison builds
are isolated in separate worktrees.

The required invariant is stronger than process success: the nonlinear
residual, `Fapply`, `Diagonal`, smoother, and physical boundary closure must be
one discrete operator.  A repair may not use a psi/coefficient floor,
artificial void stiffness, relaxed tolerance, larger iteration cap, forced
line-search acceptance, or refreshed reference to hide a defect.  Exact-zero
masked void remains outside the supported formulation unless an explicit
inactive-DOF/nullspace policy is implemented.

The assignment is treated as authorization for the plan's baseline-first
investigation.  No solver source will be edited until the first deterministic
violation and the expected effect of a candidate repair are recorded here.

## Preserved evidence

| Evidence | Status |
|---|---|
| Full report `report/output_2026-07-13_15.29.01_kermit.json` | Present in WIP worktree |
| Full HTML report with 144 sections | Present in WIP worktree |
| PlateHole retained stdout | Present in WIP worktree |
| SCPSpheresElastic retained stdout | Present in WIP worktree |
| Device lint from `benchmark/status.sh` | PASS at `a06c6f15d`, dirty 11 files |

## Baseline matrix

Results distinguish the parent `d964cfab8`, clean WIP base `a06c6f15d`, and
the preserved unstaged continuation.

| Target | Parent `d964cfab8` | Clean `a06c6f15d` | Dirty continuation | Attribution |
|---|---|---|---|---|
| PlateHole / 2D-serial | PASS/PASS; 14 iterations, `4.756775594e-7` relative | run PASS, check FAIL; solver falsely reports exact zero after one iteration and plots yield `sigxx error nan` | iteration cap, relative residual `2.023381522e-06` | paired WIP is the regression; later ghost fixes change its signature but do not restore correctness |
| SCPSpheresElastic / 2d-parallel-long | run/check PASS through 20,000 steps; 200 finite converged elastic solves | run PASS, check FAIL; displacement is NaN after a nonzero solve collapses to exact zero at step 100 | late divergent solve and NaN plots | clean paired WIP introduces the failure |
| FracturePFCZM / 2d-serial-notch | run/check PASS; 41.18 s simulation | run PASS, check FAIL; crack relative error `1.540865447e-1` after the initial nonzero solve collapses to exact zero | divergence at iteration 404 | clean paired WIP changes mechanics before the later dirty divergence |

### First deterministic attribution

The clean parent `d964cfab8` PlateHole solve and checker both pass.  Its single
linear solve reaches relative residual `4.756775594e-7` in 14 iterations.
The clean paired-stencil commit `a06c6f15d` instead reports an impossible
exact-zero residual after one iteration, exits successfully, and writes a
solution for which the unchanged checker obtains `sigxx error nan`.  The
unstaged continuation no longer has that exact-zero symptom, but stalls at
the 20-iteration cap.  Therefore the paired rewrite is a proven regression;
the later halo/smoother edits alter, but do not repair, its broken operator
hierarchy.

### Discrete defect and expected repair effect

The committed paired operator has two concrete storage inconsistencies.  Its
first low ghost-row flux reaches a cell-to-node psi average three cells beyond
the valid box, while `a06c6f15d` allocates only two psi ghosts.  It also builds
only one ghost ring of diagonal values while the inherited normalization path
divides two rings.  The dirty continuation expands the psi halo and restricts
normalization, explaining why it replaces the clean exact-zero/NaN symptom
with finite nonconvergence.

Those halo repairs cannot make the paired discretization equivalent to the
established operator.  On a uniform, single-level isotropic grid with
`u=(0,x^2 y)`, the mixed forward-gradient/backward-divergence construction
adds the grid-dependent term `h (mu-lambda)` to the x force.  It is therefore
first-order for general Lamé coefficients and changes solid-region mechanics
even when psi, AMR, multigrid, and physical boundaries are absent.  This
violates the task's solid-mechanics invariant and is consistent with the
single-level periodic Voronoi stress drift in the preserved full run.

The minimal source repair is therefore to restore the complete parent
operator stack (`Stencil`, `Elastic`, generic operator smoothing/normalizing,
and nonlinear residual) rather than retain selected pieces of the paired
rewrite.  The expected effect is exact recovery of the parent PlateHole,
SCPSpheresElastic, FracturePFCZM, and uniform-solid behavior while preserving
the already-existing no-psi material-coefficient path.  Independent narrow
changes may then reject nonfinite/failed line searches and report genuine
Newton exhaustion; they must not alter converged mechanics.

## Repair implementation

The paired-stencil changes in `src/Numeric/Stencil.H`,
`src/Operator/Elastic.H`, `src/Operator/Elastic.cpp`,
`src/Operator/Operator.cpp`, and `src/Solver/Nonlocal/Newton.H` were restored
to the exact `d964cfab8` formulation before independent safety edits.  This
preserves the existing Flame `elastic.use_psi=0` path: when disabled, Flame
does not attach psi to the solver and Elastic uses the actual material model
without applying `m_psi_small`.

The independent Newton safety edit removes forced acceptance at the final
line-search backtrack, restores the previous solution before aborting a
failed search, rejects nonfinite residual/correction/trial fields, and aborts
positive-tolerance Newton exhaustion.  The historical `nrtolerance=0`
fixed-iteration behavior remains unchanged.  The scalar line-search
acceptance predicate is exposed under `Solver::Nonlocal::Detail` and covered
by the C++ unit executable for decreases, the permitted 0.01% roundoff band,
zero residual, excessive growth, negative values, NaN, and infinity.

The new `ElasticSoftVoid` regression uses no psi and no psi floor.  Its 2-D
serial and two-rank MPI sections use 0.2 MPa void bulk/shear moduli; its 3-D
serial section uses 0.5 MPa.  The checker asserts the configured void value,
MPI process count, finite/nontrivial mechanics fields, realized material
contrast and AMR, finite linear solves, accepted line-search residuals, and
unscaled-update convergence.  It verifies the exact
`void + psi * (solid - void)` mixture on every AMR level, proves that the 2-D
diffuse interface intersects an internal coarse/fine patch boundary, and
checks parent-formulation displacement and stress witnesses (0.5% tolerance)
so a merely finite but mechanically changed stencil cannot pass.  The legacy
SCPSpheresElastic checker now also
rejects the demonstrated nonzero-RHS/zero-update false success and scans all
final mechanics fields for nonfinite values without changing its reference
data or tolerances.

### Oracle discrimination

The strengthened SCPSpheresElastic checker passes the retained clean-parent
output, including its unchanged eta/displacement/temperature reference
comparison.  Applied to the retained clean `a06c6f15d` output, it fails at
the first bad event with `nonzero nonlinear RHS 395525 accepted a zero Newton
update`, before the legacy reference checker encounters NaN.  This proves the
new oracle detects the original false-success mechanism rather than merely
following refreshed output.

The regression geometry was then matched to the assigned rod-and-tube
chamber (`r1=0.033 m`, `r2=0.046 m`) with explicit clamped 2-D corners.  A
controlled 0.2 MPa screen on the restored operator converged in 8 and 7 MLMG
iterations for that chamber.  The earlier draft's much wider soft annulus
(`r1=0.020 m`, `r2=0.060 m`) diverged on its first linear solve with or
without explicit corner declarations.  The final regression therefore tests
the requested chamber rather than silently weakening the void modulus,
linear tolerance, iteration cap, or acceptance rule.

## Focused verification

| Gate | Result | Key evidence |
|---|---|---|
| 2-D production GCC build | PASS | all 9 binaries linked |
| C++ Unit / 2-D | PASS | line-search predicate covers finite boundary cases, NaN, and Inf |
| PlateHole / 2D-serial | run/check PASS | parent-identical 14 MLMG iterations, `4.756775594e-7` relative |
| FracturePFCZM / 2d-serial-notch | run/check PASS | 400 steps; initial solve 99 MLMG iterations, `1.149237956e-4` relative to bnorm |
| VoronoiSimplePeriodic / all 2-D | 8/8 run/check PASS | serial/MPI, AMR0/2, MG0/automatic all recover parent mechanics |
| ElasticSoftVoid / 2d-serial | run/check PASS | 0.2 MPa; MLMG 8+7, final update `9.71183e-7` |
| ElasticSoftVoid / 2d-parallel | run/check PASS | two active ranks; 4/9/32 grids across AMR levels; final update `9.70531e-7` |
| SCPSpheresElastic / 2d-serial-short | run/check PASS | 1,000 steps, 10 finite elastic solves, strengthened oracle |

An initial nonfinite guard used AMReX's zero-argument `contains_nan()` and
`contains_inf()`, which intentionally inspect every allocated ghost cell.
MLMG leaves unused outer correction ghosts undefined, so that draft produced
false aborts in five Voronoi sections.  Restricting residual checks to valid
nodal equations and correction checks to the physical/periodic domain plus
the active periodic ghost shell removed all five false aborts without changing
an iteration or residual signature.  The guard synchronizes all GPU streams
before inspecting those values.

An adversarial follow-up then tested a stronger-looking alternative that
reconstructed the outer correction ghosts with a nodal FillPatch before the
finite check.  That experiment passed the narrow soft-void and Voronoi cases,
but failed the existing `EshelbyFiniteKinematics/2D-serial` numerical oracle:
the nonlinear residual plateaued near `3.6` and required 32 Newton iterations,
whereas the unchanged parent path reduced the residual to about `0.007` and
converged in 30.  MLMG already supplies meaningful active correction ghosts;
the reconstruction overwrote them and changed the mechanics.  The helper and
its calls were removed.  The final guard is therefore read-only, and the
latest exact-source rerun passes Unit, EshelbyFiniteKinematics, PlateHole,
both 0.2 MPa ElasticSoftVoid sections, and all eight 2-D Voronoi sections.

## Full-suite compatibility repair: explicit fixed nonlinear iteration

The first final-source full-suite run reached one execution failure in
`RubberPlateHole/serial-2d-coverage`.  Its coverage-only runner overrides
`solver.nriters=1` and `solver.fixed_iter=1`, has no numerical checker, but
inherits `solver.nrtolerance=1e-5` from the physical case.  The hardened
Newton contract now correctly aborts after that one update because the
measured metric is `0.05`, not below the positive tolerance.  This is an
internally contradictory coverage configuration, not a mechanics regression.

The invariant for the test repair is: a positive nonlinear tolerance requests
convergence and must abort on exhaustion; an intentionally fixed number of
nonlinear updates must opt into the existing `nrtolerance=0` mode.  The
expected change is limited to the unchecked coverage section: adding its
explicit `solver.nrtolerance=0` override lets it execute one update and exit,
while the ordinary `serial-2d` physics section retains `1e-5` unchanged.

## 3-D periodic chamber defect and source-fix invariant

The strengthened-oracle review found that the retained 3-D chamber is not the
z-invariant periodic extrusion specified by its input.  In
`ElasticSoftVoid/output_2026-07-13_17.57.41_kermit_3d-serial/00002node`, the
first z plane has `psi` in `[0.105225, 0.875]` while an interior plane has
`psi` in `[0.210451, 1]`.  The corresponding shear modulus minimum is
`15.1789 MPa` on that first plane versus `29.8579 MPa` in the interior.  This
creates `|rhs_z|=3.64465e8` and a dominant `|disp_z|=4.81584e-4` even though
the chamber expression has no z dependence.  The prior 3-D mechanics witness
therefore encoded a periodic-ghost artifact and is invalid.

The first source violation is explicit in `Flame::UpdateModel`: `phi_mf`,
`eta_mf`, and `temp_mf` call zero-argument `FillBoundary()`, which does not
copy across periodic domain faces.  `CellGradientOnNode(eta)` and
`CellToNodeAverage(eta/temp)` then read those unfilled z ghosts at the domain
edge.  The derived nodal `model_mf` is likewise finished through
`Util::RealFillBoundary`, whose current implementation ignores its Geometry
argument and also calls zero-argument `FillBoundary()`.

Before changing this source, the required invariant is: a field that is
constant in a periodic direction must remain constant through every ghost
read used to construct `psi`, body force, and material coefficients; the
derived model ghosts consumed by Newton must have the same periodicity.  The
minimal expected repair is confined to `Flame::UpdateModel`: pass
`geom[lev].periodicity()` when filling the three source fields and the derived
model, while preserving the existing multi-ghost behavior.  The repaired
3-D regression must assert z-invariance of `psi`, `model_mu/kappa`, and the
in-plane mechanics, plus negligible `rhs_z` and `disp_z`, before any CPU/GPU
performance result is accepted.  No 3-D mechanics reference will be updated
until that independent symmetry condition passes.

## Soft-void oracle audit and deterministic AMR redesign

An independent audit rejected the first `ElasticSoftVoid` acceptance oracle
despite its passing result.  The configured endpoint moduli were below
1 MPa, but the plotted diffuse mixture was not: the retained 2-D case reached
only `mu=3.890 MPa`, `kappa=4.154 MPa`, and the contaminated 3-D case reached
only `mu=15.179 MPa`, `kappa=16.231 MPa`.  The checker also counted exposed
faces of every chopped fine-grid box, including same-level neighbors, so its
reported heterogeneous coarse/fine boundary was a false positive.  In 3-D,
level 1 covered the full domain and had no internal coarse/fine boundary at
all.  Its 0.5% scalar-max mechanics tolerance would also have accepted the
discarded nodal-ghost reconstruction that failed Eshelby.

The replacement invariant is observable rather than configurational: both
plotted effective moduli must actually fall below 1 MPa; only the outer
boundary of the union of same-level fine boxes counts as a coarse/fine
interface; that union must be internal in x/y and intersect intermediate
`psi`; and the final plotted equilibrium residual must be bounded relative to
the plotted RHS.  A 2-D width screen found that `w=0.0015 m` did not converge,
whereas `w=0.002 m` converged at the unchanged tolerances and realized
`mu=0.644443 MPa`, `kappa=0.676235 MPa`; its second Newton solve backtracked to
`alpha=1/256`.  The regression now uses centered explicit nested meshes whose
outer boundaries cut the annulus, plus `w=0.002 m` and a 0.2 MPa constituent
in both dimensions.  New mechanics references will be derived only after the
periodic source repair, symmetry checks, realized-modulus check,
topology-aware coarse/fine check, and equilibrium check all pass.

## Constant Neumann cell-size truncation invariant

Source review found an independent boundary-condition regression in
`BC::Constant::FillBoundary`: `m_geom.CellSize()` returns `amrex::Real`, but
the three cached spacings are declared `int`.  For the chamber's sub-unit
physical spacing this truncates every spacing to zero, so any nonzero Constant
Neumann condition silently produces a zero-offset ghost value.  A repository
input audit found no existing nonzero Constant Neumann value, which explains
why the full suite did not expose it.

Before changing this source, the required invariant is: for a unit physical
domain with four cells (`dx=0.25`), a Constant Neumann value of 2 on a
one-cell x ghost must change the copied interior value by exactly `2*0.25`.
The minimal expected repair is to preserve the spacing as `amrex::Real` (or
equivalent inferred floating type); a direct CPU/GPU-safe C++ unit probe will
exercise both lower and upper x faces at nonzero data and would fail with the
integer declarations.

The one-line type repair and direct unit probe then passed in the 3-D GCC unit
executable.  The probe uses four cells on `[0,1]` (`dx=0.25`), interior value
10, and checks lower/upper x ghosts `9.5` and `9.25` for Neumann values 2 and
3.  It uses managed storage plus an explicit stream synchronization, so the
same probe is valid in a CUDA unit binary.  This scalar BC path is distinct
from the elastic displacement/traction BC; current Flame inputs use only zero
Neumann data, explaining why chamber mechanics did not already fail here.

## 3-D explicit coarse/fine correction-ghost failure

The first exact 3-D run after the periodic and deterministic-mesh changes
failed before its first Newton update.  On the explicit level-1 box
`[32:95]x[32:95]x[0:15]` over a `64x64x8` base grid, the first MLMG solve
reported initial RHS/residual `2.616043465e8` and final iteration-1 residual
zero, after which the hardened guard reported a non-finite Newton correction.
The pre-solve plot is finite and already proves the periodic data repair:
`rhs_z` is identically zero, `psi_min=0.00374074`, `mu_min=0.722955 MPa`, and
`kappa_min=0.760362 MPa`.  Thus this is not the former z-ghost/model artifact
and the redesigned case does realize sub-MPa material.

This mesh introduces true internal fine-level x/y boundaries that the earlier
3-D case lacked.  MLMG owns valid nodal correction rows, while some allocated
fine-level coarse/fine ghost storage may be intentionally undefined.  The
current non-finite guard checks every allocated correction point lying inside
the physical domain, which conflates those coarse/fine ghosts with valid rows.
Before changing acceptance logic, the next diagnostic will separately test
valid correction rows and the active ghost shell.  A valid-row failure means
the linear solve is broken; a ghost-only failure must be narrowed to exactly
the ghost values actually consumed by the update/residual path.  It may not be
hidden by disabling all non-finite checks or by accepting a contaminated trial.

The split valid-row check proved that the correction itself is invalid, not
merely an unused outer ghost.  Existing operator probes localized the source:
all diagonal valid entries are finite and nonzero, while the explicit fine
level has 1,584 NaNs in the grown diagonal region both before and after nodal
synchronization.  `Elastic::Diagonal` intentionally computes valid nodes plus
one ghost ring, but `Operator<Grid::Node>::normalize` divides two ghost rings.
That division reads the undefined second-ring diagonal, creates NaNs in the
normalized correction, and the subsequent MLMG operations propagate them to
valid rows while misleadingly reducing the residual norm to zero.  The same
failure occurs with a 10 MPa void and with adaptive refinement, excluding low
stiffness and explicit box placement as causes.

The repair invariant is the AMReX nodal-normalize contract: normalize only
valid algebraic rows, then synchronize/fill the correction through the
operator's boundary machinery.  AMReX's `MLNodeLaplacian::normalize` likewise
iterates the valid MultiFab region rather than dividing allocated ghosts.  The
minimal expected change is therefore the final `nghost` argument of
`a_x.divide` from 2 to 0; the finite guards remain enabled to prove that this
does not merely move the non-finite value out of sight.

That normalize-only candidate did not change the failure: valid correction
rows were still non-finite after one iteration.  Inspection of `Fsmooth`
explains why.  The custom Jacobi smoother applies and combines `Ax`, `Dx`, and
`Rx` over two correction ghost rings, so it requires a diagonal on the same
active two-ring region even when BiCGStab normalization itself is restricted
to valid rows.  The coefficient/model field already carries two ghosts, and
`Diagonal` uses an analytic impulse diagonal at each point; no third-ring
coefficient read is needed on the no-psi path.

The refined repair invariant is therefore: every diagonal entry read by the
smoother inside the physical or periodic active domain must be computed from
the same operator, while storage outside nonperiodic physical boundaries is
not an algebraic row.  The next minimal change grows the diagonal work box by
`a_diag.nGrowVect()` and grows only periodic domain directions by the same
amount.  The valid-only normalize change remains because normalization has no
reason to divide physical or coarse/fine ghosts.  The unchanged finite-row
guard will discriminate whether this restores the actual solve.

Growing the diagonal exposed, but did not itself repair, the upstream hole:
the fine grown diagonal then contained 3,216 NaNs and valid correction rows
still failed.  Their topology is the mixed corner between periodic z ghosts
and internal x/y coarse/fine ghosts.  `Flame::UpdateModel` computes the model
only on `grownnodaltilebox() & physical_domain`.  A subsequent periodic
`FillBoundaryAndSync` copies from valid regions; it cannot populate a z ghost
whose x/y coordinates are themselves ghost coordinates.  The diagonal now
correctly asks for those two-ring coefficients and reveals that they never
existed.

Flame's allocation contract already anticipates the needed direct build:
`model` and nodal `phi` carry two ghosts, while cell-centered `eta` and `temp`
carry three because `CellToNodeAverage` over the model's grown box reaches one
cell farther.  The next source invariant is therefore to evaluate the
constitutive blend over the complete allocated model grown box, including
mixed periodic/coarse-fine corners, and then synchronize shared nodal values.
Removing the physical-domain intersection from that model-build box uses the
existing documented halo contract; it does not extrapolate or regularize a
coefficient.  The finite diagonal/correction guards remain the acceptance
test.

## Closeout cleanup invariant

The split valid-row correction check was diagnostic-only: it localized the
first 3-D failure before the diagonal/model-halo repairs, but it now duplicates
the retained active-domain correction check on every Newton iteration.  The
closeout source removes only that first full-field scan.  The second check,
which covers valid rows plus the periodic/physical active ghost shell, remains
enabled; no acceptance, convergence, or finite-value behavior changes.

## Deterministic 2-D line-search oracle correction

Replaying the strengthened checker against both retained deterministic 2-D
artifacts failed only because it required at least one `alpha < 1` step.  That
condition came from the earlier adaptive width screen, whose second solve used
`alpha=1/256`; the final explicit nested mesh accepts every residual-decreasing
step at `alpha=1`.  Backtracking is conditional algorithm behavior, not a
physical invariant of this chamber.  The regression therefore continues to
require line search to be enabled, finite/consistent accepted updates, and a
nonincreasing nonlinear residual, but no longer requires this particular mesh
to trigger a rejection.  The line-search acceptance boundary and Newton
termination states remain covered directly in `src/Test/Solver/Nonlocal/Newton.H`.
An end-to-end state-rollback regression remains a future coverage item.

## Mixed-corner repair result

The three-part correction-halo repair succeeds at the failure it was designed
to fix:

1. `Operator<Grid::Node>::normalize` divides valid algebraic rows only.
2. `Elastic::Diagonal` computes every ghost entry read by the two-ring custom
   smoother, limited to the physical or periodic active domain.
3. `Flame::UpdateModel` evaluates its constitutive blend over the complete
   allocated grown model box, after periodic source fills, and synchronizes the
   resulting nodal model with the level periodicity.

On the final retained 3-D case at
`tests/ElasticSoftVoid/output_2026-07-13_19.02.34_kermit_3d-serial`, the first
two linear solves now perform real work instead of reporting a false zero:
46 iterations at `9.824901275e-6` and 29 iterations at `9.992775719e-6`.
Both Newton corrections are finite and accepted with `alpha=1`.  The final
plot is a genuine periodic extrusion: `rhs_z=0`,
`max|disp_z|=8.99364e-20`, and the per-grid z spread is exactly zero for
`psi`, `model_mu`, `model_kappa`, `rhs_x`, and `rhs_y` (in-plane displacement
spread is below `4.3e-19`).  Realized minima are `mu=0.722955 MPa` and
`kappa=0.760362 MPa`.  This closes the earlier periodic-field, mixed-corner
model, undefined diagonal, and false-zero correction defects.

The diagnostic valid-row correction scan used to localize that sequence was
removed at closeout.  The broader active-domain nonfinite correction guard
remains.  Existing `ALAMO_ML_*` coefficient/Fapply/diagonal probes in
`Elastic.cpp` predate this recovery (they are present in parent `d964cfab8`)
and were not introduced or removed here.

## Final nonlinear-equilibrium blocker

The repaired finite-soft-void case still stops on an update criterion far
before equilibrium.  This is now reproduced in both dimensions and across
MPI decomposition:

| Case | Process result | Final accepted update | Nonlinear residual ratio | Plotted `max|res|/max|rhs|` | Oracle |
|---|---|---:|---:|---:|---|
| Fresh 2-D serial, `output_2026-07-16_15.26.33_kermit_2d-serial` | Complete | `1.10789e-5` after 2 iterations | `0.236769` | `0.270767` | FAIL |
| Fresh 2-D, 2 MPI ranks, `output_2026-07-16_15.26.33_kermit_2d-parallel` | Complete | `1.10786e-5` after 2 iterations | `0.236768` | `0.270767` | FAIL |
| Retained 3-D serial, `output_2026-07-13_19.02.34_kermit_3d-serial` | Complete | `1.21842e-5` after 2 iterations | `0.236027` | `0.270142` | FAIL |

For the fresh 2-D serial plot, `max|rhs|=2.616043465e8` and
`max|res|=7.083374394e7`; realized minima are `mu=0.644443 MPa` and
`kappa=0.676235 MPa`.  The two fresh 2-D runs pass every checker condition
before force balance: completion metadata, finite/nonzero fields, linear
residuals, line-search consistency, explicit rectangular AMR hierarchy,
topology-aware coarse/fine interface intersection, exact material-mixture
identity on every level, and a positive sub-MPa realized void.  Their
identical serial/MPI signatures exclude decomposition as the cause.

A 3-D screen with the original five-iteration cap and a stricter positive
`elastic.solver.nrtolerance=1e-7` is also decisive.  It remained finite and
reduced the nonlinear residual monotonically, but exhausted all five
iterations and correctly aborted:

| Iteration | Update | Nonlinear residual ratio | MLMG iterations |
|---:|---:|---:|---:|
| 1 | `9.01532e-5` | `0.264205` | 46 |
| 2 | `1.21842e-5` | `0.236027` | 29 |
| 3 | `1.68282e-6` | `0.194880` | 122 |
| 4 | `1.16766e-6` | `0.164039` | 143 |
| 5 | `1.40792e-6` | `0.139951` | 161 |

The screen took 472.21 seconds and ended with the intended positive-tolerance
exhaustion error (`final metric=1.40792e-6`, tolerance `1e-7`).  A cheap 2-D
diagnostic with `nrtolerance=0` forced exactly the existing five iterations;
it likewise ended at update `1.03835e-6` and nonlinear residual ratio
`0.137469` (1.52 seconds).  The zero-tolerance run is diagnostic evidence only,
not an acceptable production setting.  These screens prove that lowering the
update tolerance cannot satisfy the `0.05` equilibrium gate within the
existing cap, while raising the cap would violate the task guardrail.

The chamber material models are linear, so the slow residual decay is a
strong signal to check residual/Jacobian/boundary/AMR consistency before
adding any residual-based stopping mode.  A residual mode would currently
abort too; it must not be used to conceal the mismatch.  The next investigation
should localize the maximum residual by interior, physical boundary, and
coarse/fine boundary, then compare `prepareForSolve`'s residual with `Fapply`
of the first correction on the cheap 2-D case.

## Closeout verification and environment

After the diagnostic scan removal, the configured 2-D production GCC build
completed all nine links and `bin/test-2d-g++` passed with zero failures,
including the Newton line-search/termination and Constant Neumann tests.
`benchmark/status.sh` on 2026-07-16 reports device-lint PASS.  The 3-D unit
binary passed before the diagnostic-only cleanup; the cleanup was then rebuilt
and tested in 2-D.  `git diff --check` is the final documentation/source
format gate.

One harness invocation placed the positional test directory after
`--sections`; because that option consumes a variable-length list, it started
unrelated Dendrite and ThermoElastic sections.  They were interrupted and are
not evidence.  Ignore report
`output_2026-07-16_15.24.46_kermit.{json,html}` and the corresponding killed
test directories.  Correct syntax places the test directory first, as in:

~~~
scripts/runtests.py tests/ElasticSoftVoid --comp=g++ --no-backspace \
  --no-clean --sections 2d-serial 2d-parallel
~~~

Local hardware remains an Intel Xeon w5-2545 (12 cores/24 threads) and an
NVIDIA RTX A1000 (8,188 MiB, compute capability 8.6, driver 595.71.05).
`nvcc` is not on `PATH`, and no toolkit was found under `/usr/local` or `/opt`,
so CUDA validation and all CPU/GPU timing work remain unstarted.  No commit was
made and no `results/DONE` marker exists.
## Residual localization and controlled modulus/thermal/AMR matrix (2026-07-16)

The retained final plots were inspected with the non-mutating
`benchmark/elastic_void_matrix.py` diagnostic harness.  The maximum residual
is on the level-1 outer fine-union row in both retained cases, not on a
physical boundary or in the interior:

| Plot | max field/value | level/index | coordinate (m) | classification | max rhs | ratio |
|---|---:|---|---|---|---:|---:|
| 2-D `output_2026-07-16_15.26.33.../00002node` | `res_x=+7.083374394e7` | 1 / `(64,24,0)` | `(0.1315500,0.0767375)` | outer coarse/fine-union | `2.616043465e8` | `0.270766694` |
| 3-D `output_2026-07-13_19.02.34.../00002node` | `res_y=-7.067038657e7` | 1 / `(40,0,6)` | `(0.0986625,0.0438500,0.0164625)` | outer coarse/fine-union | `2.616043465e8` | `0.270142249` |

The harness was expanded as a new diagnostic only; the existing
`tests/ElasticSoftVoid/test` oracle and retained outputs were not changed.
It stages all products under `/tmp`, runs the production `alamo` binaries,
records every Newton/MLMG line, and reports the residual maximum's level,
coordinate, and topology.  The matrix uses the retained 64x64x8 base mesh,
0.002 m interface width, explicit fine boxes, 0.2 MPa soft or 10 MPa hard
void modulus, thermal on/off, AMR on/off, and both 2-D and 3-D.

Corrected matrix report: `/tmp/alamo-elastic-void-matrix-20260716c/matrix.json`.
Completed thermal-on/AMR-on soft cases reproduce the retained behavior:
2-D has MLMG `11,11`, Newton nonlinear ratios `0.264885 -> 0.236769`, and
plotted ratio `0.270767`; 3-D has MLMG `46,29`, ratios `0.264205 -> 0.236027`,
and plotted ratio `0.270142`.  Hard 10 MPa cases remain finite and linear
MLMG-converged but still stop at nonlinear ratios about `0.236` and plotted
ratios `0.281393` (2-D) and `0.280606` (3-D) with AMR/thermal on.

With AMR off and thermal on, hard cases complete two Newton updates with
MLMG 6/6 (2-D) or 7/7 (3-D), nonlinear ratios about `0.21484`/`0.21486`,
and plotted ratios about `0.23425`; the soft 2-D case fails line search after
8 backtracks, while the soft 3-D case fails MLMG after 200 iterations at
relative residual `2.60934e19`.  Thermal off aborts in the phase-field kernel
with a non-finite value before an elastic solve for all AMR-on/off modulus
cases.  These are diagnostic outcomes, not acceptance results.

The plotted ratio is `max_{valid plotted nodes, components}|res| /
max_{valid plotted nodes, components}|rhs|`: it measures remaining global
discrete force imbalance against the applied elastic RHS.  It is independent
of MLMG's per-linear-solve relative residual, which is only the error in the
current Jacobian correction solve.  Thus MLMG values near `1e-5` can coexist
with a nonlinear force ratio near `0.24` when the assembled correction does
not remove the nonlinear residual consistently; update-based Newton stopping
at a small step does not establish equilibrium.

## Cross-worktree audit and discrete invariant correction (2026-07-16)

The earlier statement that Flame's chamber material is linear was incorrect.
`Flame` instantiates `Model::Solid::Finite::NeoHookeanPredeformed`; both `DW`
and `DDW` depend on the deformation gradient.  A future one-correction linear
oracle must therefore use an actual `Model::Solid::Linear::*` model rather than
infer linearity from the small displacement in this regression.

The cross-worktree/history audit found two independent repairs that must be
preserved: coefficient hierarchy resynchronization after every Newton
relinearization, and the GPU interpolation temporary's `Elixir` lifetime.
Coarsening caps, bottom-solver substitutions, extra smoothing, diagonal
inflation, modulus/psi floors, and the rejected paired `D- G+` rewrite did not
repair equilibrium.  A literal centered-gradient/centered-divergence tangent
was also already tested: it matched its manufactured operator to
`2.93e-15`, but its checkerboard nullspace made MLMG grow from `0.256` to
`1.31e21` in 15 iterations.  It is not a viable consistency repair.

The source-level invariant to enforce before any production repair is:

~~~
J_h(u) v = d R_h(u + epsilon v) / d epsilon at epsilon = 0
~~~

on every active algebraic row, with the same physical-boundary closure and
AMR composite-row policy.  Current code violates this even on a uniform,
single-level interior.  Newton forms `R_h` by a centered divergence of stress
computed from a centered displacement gradient.  `Elastic::Fapply` instead
uses one-cell Hessians plus an explicit coefficient-gradient product-rule
term.  In one dimension with constant coefficient `C`, the first construction
linearizes to

~~~
C (v[i+2] - 2 v[i] + v[i-2]) / (4 h^2),
~~~

whereas `Fapply` uses

~~~
C (v[i+1] - 2 v[i] + v[i-1]) / h^2.
~~~

The report's non-divergence-form concern is therefore confirmed as a concrete
discrete residual/Jacobian mismatch.  The retained C/F maxima are an
amplification/localization of that mismatch, not proof that C/F interpolation
alone is causal: the hard, no-AMR controls retain force ratios near `0.234`.
The current C/F `reflux` replaces coarse C/F/covered residual rows with a
full-weighted fine residual; it is not a conservative flux-register
correction.  A separate source audit also found that AMR RHS/model ghost
freshness is not fully established before `TimeStepBegin` mechanics, so ghost
lifecycle remains a secondary test target.

Current smoother facts also correct the publication-era report: `Fsmooth` is
weighted Jacobi with default `omega=2/3`, two sweeps per call, and the
soft-void input requests four pre- and four post-smoothing calls.  Smoother
tuning can affect linear robustness but cannot make a mismatched Jacobian the
derivative of the residual.  Material interpolation is arithmetic by design
in Flame and is asserted by the existing regression; harmonic averaging has
not been tested and must not silently replace the phase-field free-energy
mixture.  It is a possible face-flux homogenization choice only after a
compatible conservative formulation exists.

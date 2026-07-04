# GPU Elastic Fix Notes

Date: 2026-06-21

## ✅ FIXED 2026-06-25 — the actual root cause and one-line fix

**Root cause:** a GPU **cross-stream use-after-free race** on the per-box
temporary `FArrayBox tmpfab` in `Operator<Grid::Node>::interpolation()`
(`src/Operator/Operator.cpp`). AMReX's non-tiled GPU `MFIter` runs iterations
across a pool of CUDA streams; `tmpfab` (created/destroyed each iteration) was
freed at iteration end and its device block reused by a later iteration on a
**different stream** while the interpolation `ParallelFor`/`plus` kernels were
still in flight. The corrupted coarse-grid correction was then amplified by the
BC-penalty diagonal (~1e17) into the 1e20 blow-up. `restriction()` writes its
output directly (no temp), which is why it was always clean and only
interpolation broke.

**Fix:** `amrex::Gpu::Elixir tmpfab_eli = tmpfab.elixir();` immediately after
`tmpfab.resize(...)` — defers the temp's free until its own stream completes.

**How localized (this session):** box-size sweep (single box converges, any
multi-box decomposition diverged → cross-box); stock-AMReX `MLNodeLaplacian`
reproducer converges multi-box (→ ALAMO custom code, not AMReX/BiCGStab);
bit-identical-input CPU↔GPU transfer probe (→ `interpolation`, not
`restriction`); pre/post-`nodalSync` bisect (→ per-box temp, nondeterministic
run-to-run).

**Reconciles the entries below — all were downstream symptoms of this race:**
the BiCGStab noise-floor, the "nondeterministic reduction", the "deterministic
coarse-level operator defect", and the flaky 2048² threshold (a box-count/box-
size effect, not a clock/thermal artifact). The `Fapply`/`m_ddw_mf`/`m_diag`
exonerations below were all CORRECT — the bug was never in operator application
or coefficients; it was in the transfer temp's lifetime. Everything below is
retained as investigation history.

**Verified:** full box sweep all converge ~1e-9; real production operator
(use_psi=1, 3-material) multi-box 2048² converges 7 iters → 8.5e-9; end-to-end
flame+elastic chamber sim (forced multi-box, 300 steps) clean — 10 elastic
solves all <1e-8, sane physics (P=5.38 MPa, mdot=3.65 kg/s).

## STATUS CORRECTION 2026-06-22 (late — read before trusting conclusions below)

Two conclusions recorded below are now corrected. Current ordered priorities live
in `GPU_ELASTIC_DEBUG_PLAN.md` ("CURRENT PRIORITIES").

1. **The "Fapply itself is correct at both mglev=0 and mglev=1" conclusion
   (Fapply Single-Kernel Probe section) is OVER-CLAIMED.** That probe used only
   constant and linear displacement fields, which lie in the elastic operator's
   near-null-space (exact `A·u ≈ 0`). It cannot detect a defect that manifests
   only on rough/high-frequency inputs — which is exactly the failing correction.
   `Fapply` is NOT exonerated. The discriminating test (high-frequency analytic
   field, GPU-vs-CPU at mglev=1) is Priority 1 in `GPU_ELASTIC_DEBUG_PLAN.md`.
2. **The ret==9 acceptance-policy hypothesis is now fully REFUTED.** The redo on
   the correct AMReX tree (`ext/AMReX-Codes/amrex`,
   `log_ret9patch_coarsen1_2048_uniform_i60.log`, 13:15) still diverges
   148 → 1.75e20. Combined with the smoother-only GPU stall (zero reductions,
   deterministic), the standing cause is a **deterministic GPU-specific
   coarse-level (mglev≥1) operator-application defect**, not the bottom-solver
   return-code policy and not a nondeterministic reduction. Do not resume the
   bicgstab-acceptance or coefficient threads.

## Current Narrowing

The smallest failing hierarchy is the uniform 2048x2048 GPU case with
`elastic.max_coarsening_level=1`, i.e. one fine level plus one coarse/bottom MG
level. `elastic.max_coarsening_level=0` did not diverge before the harness
timeout.

Default bottom solver behavior is the sharpest trigger observed so far:

| Case | Bottom solver | Result | Final `resid/bnorm` | Notes |
|---|---|---|---:|---|
| `coarsen1_2048_uniform_i60` | default `bicgstab` | explosive divergence after 12 iterations | `2.624074357e+20` | Pre-fix baseline. |
| `coarsen1_2048_uniform_i60_transferfix` | default `bicgstab` | explosive divergence after 13 iterations | `2.648627461e+21` | After explicit 2D nodal transfer branches. Trajectory changed, failure remains. |
| `coarsen1_2048_uniform_i60_bv5` | default `bicgstab` | explosive divergence after 12 iterations | `4.118904137e+20` | `elastic.solver.verbose=5`; first V-cycle bottom `UP=2.906744463e+05`, then fine `UP before smooth=1.025065514e+11` and fine `resid/bnorm=128.1298774`. |
| `coarsen1_2048_uniform_i12_bottom_cg_bv5` | `cg` | failed at `max_iter=12` | `0.2374600709` | Bottom CG hits its internal 201-iteration cap but stays bounded; first V-cycle bottom `UP=2.858193165e+07`, fine `resid/bnorm=1.925044027`. |
| `coarsen1_2048_uniform_i12_bottom_smoother_bv5_serial` | `smoother` | failed at `max_iter=12` | `3.031934833` | Serial run; first concurrent attempt exhausted AMReX arena. First V-cycle bottom `UP=2.564384674e+07`, fine `resid/bnorm=2.647172004`. |
| `coarsen1_2048_uniform_i12_bicg_btol5e2_bv5` | `bicgstab`, `bottom_tol_rel=0.05` | failed at `max_iter=12` | `4.361997343e+09` | Relaxing bottom tolerance does not remove the first-cycle amplification. |
| `coarsen1_2048_uniform_i12_bicg_final0_bv5` | `bicgstab`, `final_smooth=0` | explosive divergence after 12 iterations | `1.076040848e+20` | Bottom final smoothing is not the primary trigger; coarse-to-fine correction still amplifies strongly. |
| `coarsen1_2048_uniform_i12_bicg_cordiag1` | `bicgstab`, `ALAMO_ML_COR_DIAG=1` | explosive divergence after 12 iterations | `1.668333904e+20` | First interpolation call: coarse correction after bottom `norminf=8.80301e-07`, `l2=5.32357e-04`; fine before interpolation `norminf=2.38418e-10`, `l2=6.80096e-08`; fine after interpolation/sync `norminf=8.78737e-07`, `l2=6.20160e-04`; AMReX fine `UP before smooth=1.023226618e+11`. Correction magnitude stays bounded, but subsequent fine residual explodes. |
| `coarsen1_2048_uniform_i12_cg_cordiag1` | `cg`, `ALAMO_ML_COR_DIAG=1` | failed at `max_iter=12` | `1.502793872` fine `resid/bnorm` | First interpolation call: coarse `norminf=1.08308e-08`, `l2=1.89467e-06`; fine before `norminf=1.43066e-09`, `l2=4.79534e-07`; fine after/sync `norminf=1.1906e-08`, `l2=2.34079e-06`; fine `UP before smooth=1.990960594e+09`. Fine residual remains O(1). |
| `coarsen1_2048_uniform_i12_smoother_cordiag1` | `smoother`, `ALAMO_ML_COR_DIAG=1` | failed at `max_iter=12` | `1.254453504` fine `resid/bnorm` | First interpolation call: coarse `norminf=1.08308e-08`, `l2=1.89467e-06`; fine before `norminf=1.43066e-09`, `l2=4.79534e-07`; fine after/sync `norminf=1.20608e-08`, `l2=2.85432e-06`; fine `UP before smooth=2.003429621e+09`. Fine residual remains O(1). |
| `coarsen1_2048_uniform_i60_bottom_cg` | `cg` | failed at `max_iter=60` | `0.09265647885` | No explosive divergence. |
| `coarsen1_2048_uniform_i60_bottom_smoother` | `smoother` | failed at `max_iter=60` | `0.7369654712` | No explosive divergence. |
| `coarsen1_2048_uniform_i12_bicg_resdiag12` | default `bicgstab`, `ALAMO_ML_COR_DIAG=1`, `ALAMO_ML_RESID_DIAG=12` | explosive divergence after 12 iterations | `5.062253594e+19` fine `resid/bnorm` | First bottom correction on `mglev=1`: valid-region `norminf=8.80301e-07`, `l2=5.32357e-04`, but `contains_nan=1` and `contains_inf=1`. After interpolation/sync on fine `mglev=0`: `norminf=8.79118e-07`, `l2=6.71301e-04`, `contains_nan=0`, `contains_inf=0`. The first fine operator application on this correction gives `A*x norminf=1.02367e+11`, matching `UP: Norm before smooth=1.023665476e+11`. |
| `coarsen1_2048_uniform_i2_cg_resdiag4` | `cg`, `ALAMO_ML_COR_DIAG=1`, `ALAMO_ML_RESID_DIAG=4` | failed at `max_iter=2` as expected | `1.606864833` fine `resid/bnorm` | First bottom correction on `mglev=1`: `norminf=1.08308e-08`, `l2=1.89467e-06`, `contains_nan=0`, `contains_inf=0`. First fine operator application after interpolation gives `A*x norminf=1.91323e+09`, with fine `resid/bnorm` remaining O(1). |
| `coarsen1_2048_uniform_i2_smoother_resdiag4` | `smoother`, `ALAMO_ML_COR_DIAG=1`, `ALAMO_ML_RESID_DIAG=4` | failed at `max_iter=2` as expected | `1.590027202` fine `resid/bnorm` | First bottom correction on `mglev=1`: `norminf=1.08308e-08`, `l2=1.89467e-06`, `contains_nan=0`, `contains_inf=0`. First fine operator application after interpolation gives `A*x norminf=1.25921e+09`, with fine `resid/bnorm` remaining O(1). |

## Code Change Under Test

`src/Operator/Operator.cpp` now has explicit `AMREX_SPACEDIM == 2` branches in
custom nodal `restriction` and `interpolation`. The previous code used the 3D
stencil logic unconditionally in 2D and could reference `k - 1` / `k + 1` in
transfer stencils. This is still a correctness fix candidate, but it is not
sufficient to resolve the 2048x2048 two-level divergence.

## Latest Diagnostic Result: 2026-06-22

`ALAMO_ML_RESID_DIAG=N` now instruments `correctionResidual` around the operator
application that produces AMReX's `Norm before smooth` values. It records the
correction `x`, RHS `b`, raw `A*x`, and final residual before/after sync.

The first bad amplification is now localized to applying the fine-level operator
to the interpolated `bicgstab` bottom correction:

- before interpolation, the default `bicgstab` bottom correction on `mglev=1` has
  finite valid-region reductions but reports `contains_nan=1` and
  `contains_inf=1` somewhere in the `MultiFab` storage/ghost region;
- interpolation plus fine sync removes those invalid values from the fine
  correction as detected by `contains_nan=0`, `contains_inf=0`;
- the fine correction magnitude is still only O(`1e-3`) in L2, but the immediate
  fine operator application produces `A*x norminf=1.02367e+11`, which is the
  observed first `UP: Norm before smooth` blow-up;
- `cg` and `smoother` produce clean bottom corrections (`contains_nan=0`,
  `contains_inf=0`) around `1.9e-06` L2 and fine `A*x` around `1e9`, with fine
  residual ratios staying O(1).

## Diagnostic Update: 2026-06-22 Ghost/Valid Classification

Added explicit valid-region versus grown-region NaN/Inf flags to the existing
`ALAMO_ML_COR_DIAG` / `ALAMO_ML_RESID_DIAG` output in `src/Operator/Operator.cpp`.
Also added an env-gated A/B experiment, `ALAMO_ML_ZERO_CRSE_GHOSTS=1`, which
zeros the coarse correction ghosts before interpolation while preserving valid
nodal data.

New short runs used the uniform 2048x2048, two-level case with
`elastic.solver.max_iter=2`:

| Case | Bottom solver / env | First bottom correction flags | First fine correction after interpolation | First `UP: Norm before smooth` | Iteration 1 `resid/bnorm` | Interpretation |
|---|---|---|---:|---:|---:|---|
| `coarsen1_2048_uniform_i2_bicg_ghostdiag` | default `bicgstab` | valid clean, grown dirty: `valid_contains_nan=0`, `valid_contains_inf=0`, `grown_contains_nan=1`, `grown_contains_inf=1` | `l2=7.13073e-04`, clean | `1.019629288e+11` | `120.5511062` | Invalids are outside the valid region, but the valid correction still triggers the huge fine operator response. |
| `coarsen1_2048_uniform_i2_bicg_ghostzero` | `bicgstab`, `ALAMO_ML_ZERO_CRSE_GHOSTS=1` | ghost-zeroed copy is fully clean | `l2=6.07807e-04`, clean | `1.022714492e+11` | `133.1733782` | Zeroing invalid coarse ghosts does **not** remove the amplification. Ghost invalids are real, but not the primary cause of the first blow-up. |
| `coarsen1_2048_uniform_i2_cg_ghostdiag` | `cg` | fully clean | `l2=2.42014e-06`, clean | `1.917575955e+09` | `2.109561211` | Bounded path; bottom correction is about 295x smaller in L2 than `bicgstab`. |
| `coarsen1_2048_uniform_i2_smoother_ghostdiag` | `smoother` | fully clean | `l2=3.80679e-06`, clean | `2.309435902e+09` | `2.979791965` | Bounded path; same order as `cg`. |

Conclusion: the previous `contains_nan=1` / `contains_inf=1` on the
`bicgstab` bottom correction is confined to grown/ghost storage, not valid nodal
data. Cleaning those ghosts changes the interpolated correction slightly but does
not change the first bad amplification. The primary defect is therefore not
"invalid ghost values are read by interpolation." The sharper target is the valid
`bicgstab` bottom correction itself: it is hundreds of times larger than the
`cg`/`smoother` correction and maps through the fine operator to an O(`1e11`)
response.

## Diagnostic Update: 2026-06-22 Spatial Classification

Added `ALAMO_ML_SPATIAL_DIAG=N` in `src/Operator/Operator.cpp`. It is env-gated
like the existing norm diagnostics and records per-component min/max coordinates,
max absolute value location, domain-boundary/box-boundary/interior max absolute
values, total absolute-value fractions on boundaries, and a nearest-neighbor
jump metric. It is wired into:

- coarse correction before interpolation (`crse_after_bottom`);
- fine correction before/after interpolation;
- `correctionResidual` before `apply` and immediately after `A*x`.

Short two-iteration runs:

| Case | Bottom solver | First bottom correction spatial summary | First fine correction after interpolation | First fine `A*x` spatial summary | Fine residual path |
|---|---|---|---|---|---|
| `coarsen1_2048_uniform_i2_bicg_spatialdiag` | default `bicgstab` | interior extrema, not domain boundary; comp0 max abs `8.80301e-07` at coarse `(793,465)`, comp1 max abs `8.73755e-07` at `(507,750)`; box-boundary sum fraction ~`0.060` | comp0 max abs `8.80298e-07` at fine `(1330,930)`, comp1 `8.73801e-07` at `(886,1500)`; neighbor jump is O(max) (`0.998`/`1.000` of max) | first two fine `A*x` maxima are interior, O(`1e7`) before interpolation-bottom correction and O(`3e7`-`6e7`) in early correction residual calls; after interpolation, AMReX reports fine `resid/bnorm=129.808` then `12520.198` | abort at `max_iter=2`, final `resid/bnorm=12520.19781` |
| `coarsen1_2048_uniform_i2_cg_spatialdiag` | `cg` | interior extrema, clean; comp0/1 max abs ~`1.083e-08`, about 81x smaller in norminf than `bicgstab`; box-boundary sum fraction ~`0.05`-`0.06` | fine max abs ~`1.19e-08`, domain-boundary max effectively zero; neighbor jumps O(max), but absolute magnitude remains small | fine `A*x` maxima are interior and O(`4e6`-`9e6`) in early correction residual calls; fine `resid/bnorm` stays bounded (`1.640` then `0.475`) | abort at `max_iter=2`, final `resid/bnorm=0.4753410389` |

Interpretation: the bad `bicgstab` correction energy is not concentrated on the
physical domain boundary, and not solely explained by box/tile edges. The first
coarse correction extrema are interior for both `bicgstab` and `cg`, but the
default `bicgstab` bottom correction is ~80x larger in norminf and hundreds of
times larger in L2 than the bounded path. Bilinear interpolation of both paths
creates nearest-neighbor jumps O(max), but only the much larger `bicgstab`
correction maps to the explosive fine residual. The next target is therefore the
bottom `bicgstab` correction quality/failure policy and whether its returned
valid correction is numerically compatible with the fine operator after a failed
bottom solve, not a simple boundary-node or ghost-read defect.

## Diagnostic Update: 2026-06-22 Fapply Single-Kernel Probe (mglev=0 vs mglev=1)

Added `ALAMO_ML_FAPPLY_PROBE=<mglev>` to `src/Operator/Elastic.cpp`. On the
first `Fapply` call at the requested `mglev`, it overrides the input
displacement field with a known analytic field (constant, then a linear ramp
scaled to a physically-representative ~1e-6 m magnitude instead of raw domain
coordinates), captures `gradgradu.norm()` and the resulting `f` interior/
boundary max-abs, and aborts — letting a single MG-level Fapply evaluation be
checked against the analytic answer (`f=0` in the interior) without running a
full solve.

Raw-coordinate-scale linear field (pre-fix, round-off swamped, not meaningful
on its own): `mglev=0` gave `interior_max_abs=1.18685`, `mglev=1` gave
`0.296887` — `mglev=0` four times worse, which looked at first like the GPU
second-derivative stencil might be specifically broken at the fine level.

Corrected physically-scaled linear field, same uniform 2048x2048
`max_coarsening_level=1` case, both levels now rerun on the current build:

| mglev | `gradgradu.norm()` interior_max | `f` interior_max_abs |
|---|---:|---:|
| 0 | `3.06209e-14` | `9.04966e-06` |
| 1 | `7.65521e-15` | `2.26242e-06` |

`mglev=0` is ~4x worse than `mglev=1` in both metrics, matching the `1/dx^2`
round-off-amplification prediction exactly (mglev=0's `dx` is half of
mglev=1's, so `1/dx^2` is 4x larger for the same fixed-magnitude floating-point
error in `u`). Both values are negligible (1e-5/1e-6) against the problem's
RHS scale (~3.4e7).

Conclusion: **Fapply itself is correct at both mglev=0 and mglev=1** for this
analytic test. The earlier 4x disparity was round-off amplification, not a
GPU-specific second-derivative defect. This closes the "broken Fapply kernel"
hypothesis; it does not explain the bicgstab bottom-solve divergence, which
remains the active target below.

## CORRECTION: 2026-06-22 Wrong AMReX Source Tree

The repo has two vendored AMReX trees: `ext/amrex` and `ext/AMReX-Codes/amrex`.
`configure` defaults to `--build-amrex-fork=AMReX-Codes`, and
`.make/Makefile.pre.conf` confirms `AMREX_TARGET =
ext/AMReX-Codes/amrex/2d-nofast-cuda86-g++-26.06` — i.e. the headers actually
installed and compiled into every local GPU binary come from
`ext/AMReX-Codes/amrex/Src/...`, not `ext/amrex/Src/...`.

The "ret==9 Acceptance-Policy Patch Test" below was applied to
`ext/amrex/Src/LinearSolvers/MLMG/AMReX_MLMG.H`, which is **not** the tree
used by the build. Confirmed via `strings` on the rebuilt binary: none of the
edited code/markers were present. **The "negative result" (divergence got
worse) recorded below is invalid** — the patch was never actually compiled
in; the rebuild just reproduced the unpatched baseline with normal run-to-run
GPU nondeterminism (consistent with the existing nondeterminism finding, not
with the patch). Treat the section below as void pending a redo on
`ext/AMReX-Codes/amrex/Src/...`.

## Diagnostic Update: 2026-06-22 ret==9 Acceptance-Policy Patch Test (NEGATIVE) — VOID, WRONG TREE, SEE CORRECTION ABOVE

Tested the hypothesis that AMReX's bottom-bicgstab "accept partial progress"
policy (`ext/amrex/Src/LinearSolvers/MLMG/AMReX_MLMG.H`,
`actualBottomSolve()`: `if (ret != 0 && ret != 9) setVal(cor[...], 0.0);`) was
letting a GPU-nondeterminism-corrupted, non-converged bicgstab correction
(`ret==8` relabeled `ret==9` when `rnorm < rnorm0`) through to the fine-level
interpolation, producing the observed amplification.

Patch: changed the condition to `if (ret != 0)`, removing the `ret==9`
exemption so *any* non-strictly-converged bicgstab bottom correction is zeroed,
identical treatment to the hard-breakdown codes 1-4. Rebuilt
`bin/alamo_gpu-2d-nofast-cuda86-g++` (`./configure --comp=g++ --dim=2 --cuda
local --cuda-fp strict && make`).

Rerun of the same uniform 2048x2048, `max_coarsening_level=1`, default
`bicgstab`, `max_iter=60` case:

- before patch: explosive divergence after 12 iterations, `resid/bnorm =
  2.624074357e+20`.
- after patch: explosive divergence after **13** iterations, `resid/bnorm =
  5.910612351e+27` — *worse*, not better.

**Conclusion: the `ret==9` partial-acceptance policy is NOT the cause.**
Zeroing every non-strictly-converged bicgstab correction does not prevent or
even delay the explosive blowup; it changes the trajectory slightly (one extra
iteration) but the magnitude got two orders of magnitude worse. This rules out
"bad correction sneaks through under the ret==9 exemption" as the mechanism.

Implication: either (a) bicgstab is reporting `ret==0` (strict convergence) on
GPU at the point of divergence — i.e. the convergence test itself is being
satisfied by a numerically wrong state, not skipped by the ret==9 path — or
(b) the blowup is not rooted in the bottom bicgstab return-code handling at
all, and lies elsewhere (smoother, restriction/interpolation, or the residual
computation/V-cycle bookkeeping around the bottom solve). The `ret != 0` patch
in `AMReX_MLMG.H:1637` should be reverted (or left as a harmless stricter
policy) before any further upstream consideration; it is not a fix.

Next angle: instrument bicgstab's *own* reported `ret` and `rnorm`/`rnorm0` at
the point it returns for the first divergent V-cycle (does it actually claim
`ret==0`?), and separately check whether the explosive growth still appears
under `bottom_solver=smoother` (which bypasses bicgstab/this code path
entirely) at the same `max_iter=60` setting — `fix_notes.md`'s earlier
`coarsen1_2048_uniform_i60_bottom_smoother` run did NOT explode (capped at
`resid/bnorm=0.737`), which argues for (a) — bicgstab specifically, likely
claiming false convergence — over (b).

## Diagnostic Update: 2026-06-22 bicgstab ret/rnorm Instrumentation (correct tree)

Lesson learned: this repo vendors *two* AMReX trees, `ext/amrex` and
`ext/AMReX-Codes/amrex`. `configure` defaults to
`--build-amrex-fork=AMReX-Codes`, and `.make/Makefile.pre.conf` /
`.make/Makefile.post.conf` confirm the actual installed/compiled headers come
from `ext/AMReX-Codes/amrex/Src/...` (with a real Make dependency on
`AMREX_SRC` that triggers an amrex rebuild+reinstall when those files change).
`ext/amrex` is unused by the current build. Any future AMReX-side patch/probe
must go in `ext/AMReX-Codes/amrex/Src/...`; verify with `strings
bin/<binary> | grep <marker>` before trusting a "negative" result.

Added `ALAMO_ML_BICGSTAB_RET_DIAG=N` to
`ext/AMReX-Codes/amrex/Src/LinearSolvers/MLMG/AMReX_MLCGSolver.H`
(`solve_bicgstab`), logging the raw (pre-relabel) `ret`, `iter`, `rnorm`,
`rnorm0`, `rnorm/rnorm0`, and whether the correction will be accepted, for the
first N bottom-solve calls. Verified present in the rebuilt binary via
`strings`.

Rerun: uniform 2048x2048, `max_coarsening_level=1`, default `bicgstab`,
`max_iter=20`, `N=20`:

- Every single bottom-solve call returns `raw_ret=8` (hit `bottom_maxiter=201`
  without meeting `eps_rel=1e-4`) and is relabeled/accepted as `9`. **None
  report `ret==0`** — bicgstab is never falsely claiming strict convergence.
- Each call's own local quality is sane: `rnorm/rnorm0` ranges `0.0019`-`0.28`
  across all 20 calls — a real, bounded local reduction every time, not a
  blown-up or garbage correction.
- `rnorm0` (the bottom problem's RHS/initial-residual norm, i.e. essentially
  the restricted fine residual) grows by a strikingly consistent **~40x-150x
  per outer MLMG iteration**: `1.16e-9 -> 1.70e-7 -> 1.49e-5 -> 6.48e-4 ->
  2.42e-2 -> 0.603 -> 26.2 -> 637 -> 32035 -> 673128` (values shown are the
  first bicgstab call's `rnorm0` per outer iteration). `Fine resid/bnorm`
  grows in lockstep: `159.98 -> 15094 -> 636519 -> 23223291 -> ...`.
- Critically, **iteration 1 is already diverging** (`Fine resid/bnorm =
  159.98`, far above 1) — the instability is present from the very first
  V-cycle, not something that builds up gradually over many iterations.

**Conclusion: the bottom-bicgstab return-code/acceptance policy is not the
locus of the bug.** Every accepted correction is a locally reasonable,
bounded-ratio result. The ~100x-per-V-cycle amplification must be happening in
how that locally-sane correction is handled on the fine level — interpolation,
post-bottom smoothing, or the residual/restriction bookkeeping around the
coarse-fine boundary — not in bicgstab's internals. This is consistent with
the standing `ALAMO_ML_COR_DIAG` finding (a numerically modest `l2~5e-4`
correction produces `Norm before smooth=1.02e11` immediately after
interpolation onto the fine level): the explosion happens *after* the bottom
solve returns a sane result, on the way back up.

## Diagnostic Update: 2026-06-22 CPU-side `SetUniform(true)` diagnostic flag (separate from the GPU interpolation bug above)

Investigated a separate report (low void-stiffness CPU elastic failures) by
re-reading `src/Integrator/Base/Mechanics.H` against `src/Operator/Elastic.cpp`.
Found `src/Integrator/Base/Mechanics.H:194` had been left as
`elastic_op.SetUniform(true); // TEMP DIAGNOSTIC - revert before merging`
(an **uncommitted** local working-tree edit, confirmed via `git diff master`;
the master/last-committed value is `SetUniform(false)`). `SetUniform(true)`
gates off the `grad(C)·grad(u)` coefficient-gradient term in
`Elastic.cpp:419` (`if (!m_uniform) { ... }`). For a masked elastic operator
with a low-stiffness void region, the stiffness field is sharply nonuniform
between solid and void, so this term is exactly the one that matters — solving
with `SetUniform(true)` silently drops it and solves the wrong PDE. This
applies identically to both the CPU and GPU code paths since both share
`Elastic.cpp`.

Action taken: reverted to `elastic_op.SetUniform(false)`. Rebuilt
`bin/alamo-2d-g++` (CPU target; build config backed up/restored around the
rebuild since `.make/Makefile.pre.conf` was pointed at the GPU target for
another concurrent investigation).

Validation run (CPU, `bin/alamo-2d-g++`, `np=4`), ad-hoc input based on
`run_stability_sweep.sh`'s `VD_4_30` config (`psi_floor=0.05`,
`model_void.kappa=4_MPa`, `model_void.mu=3_MPa`, propellant/casing
`kappa=162_MPa`, `mu=113.6_MPa`):

- `stop_time=0.05_s` (500 steps): **PASS** — 0 MLMG failures, 0 NaN, exit 0,
  ran to completion.

Re-run at `stop_time=1.5_s` with void stiffness relaxed further to
`model_void.kappa=model_void.mu=0.2_MPa` (still with the `SetUniform(false)`
fix in place): **fails — Newton stall, not a crash.** Clean (1 Newton
iteration per elastic solve) through `t≈0.49s`. At the elastic solve
triggered at `t=0.52s` (`step=5200`... actually logged at the solve starting
around `step≈5226`, `interval=25`), Newton iteration count climbs instead of
terminating:

- `NR iteration` relative `norm(ddisp)`: `1.3e-5 -> 2.0e-5 -> 3.1e-5 -> 4.7e-5
  -> 7.2e-5 -> 1.1e-4 -> 1.7e-4 -> 2.6e-4 -> 3.9e-4 -> 6.0e-4 -> 8.9e-4 ->
  1.30e-3 -> 1.79e-3 -> 2.22e-3 -> 2.42e-3 -> 2.33e-3 -> 2.16e-3 -> 2.09e-3 ->
  2.13e-3 -> 2.20e-3 -> 2.23e-3 -> 2.22e-3` (iterations 1-22).
- Roughly geometric growth (~1.3-1.5x/iter) through iteration ~15, then it
  does **not** cross into NaN/overflow territory — it flattens into a tight
  limit cycle around `2.1e-3 - 2.4e-3` and stays there with no further trend
  for 6+ iterations. `nrtolerance=1e-5` is never approached.
- Each Newton iteration costs a full MLMG solve (~150 V-cycles, ~5s wall each
  at this point vs. ~18 V-cycles/~0.5s earlier in the run when converging
  cleanly in 1 iteration) — so this single elastic solve is on track to burn
  through all `nriters=200` (per `Newton.H:361-366`, exceeding `m_nriters`
  does **not** abort; the loop just breaks and the simulation continues with
  the unconverged displacement field) at a cost of ~15-20 minutes for this one
  solve alone.
- No NaN, no MLMG `Abort`/exception — this is a genuine Newton stall/limit
  cycle at low void stiffness, not the GPU-side interpolation blow-up
  documented above. Matches the user's prior report of manually killing an
  earlier run of this same config because "the Newton solver wasn't
  converging."
- Implication: `SetUniform(false)` (the correct setting) is necessary but not
  sufficient — once void stiffness drops far enough (here `0.2_MPa`, vs. the
  `elastic_stability_sweep.md` campaign's previously-found safe floor of
  `~4-6_MPa`), the masked elastic operator's conditioning is bad enough that
  Newton itself stalls regardless of the coefficient-gradient term being
  present. This is consistent with `elastic_stability_sweep.md`'s standing
  finding that `rod_and_tube`-style geometries stall (rather than blow up)
  as void/psi_floor are relaxed too far.

**Caveat found along the way (unrelated to either bug above):** the stock
`run_stability_sweep.sh` input files (`input_rod_and_tube_new`,
`input_anchor_neutral_new`, etc.) set `propellant.type = homogenize` but still
carry the old resolved-model blocks (`model_ap`/`model_htpb`) instead of the
required `model_prop`/`model_casing`. This hits `query_exactly<2>`'s parse
abort immediately (0 of the 5 `{lambda,mu,E,nu,kappa}` keys present) — this is
the same "model_prop query_exactly<2> parse abort" bug already catalogued as
one of the two non-dispatch chamber-lineage Flame bugs blocking DoD item 4. It
currently blocks `run_stability_sweep.sh` outright; validation above used a
scratch copy (`/tmp/test_setuniform_input`) with `model_prop`/`model_casing`
added by hand, not an edit to the tracked input files.

## Diagnostic Update: 2026-06-22 Bottom-Smoothing Gap Closed (correct tree)

Added `ALAMO_ML_BOTTOM_SMOOTH_DIAG=N` to `actualBottomSolve()` in
`ext/AMReX-Codes/amrex/Src/LinearSolvers/MLMG/AMReX_MLMG.H`, bracketing the one
checkpoint not covered by existing probes: the bottom correction `x`
immediately as `bottomSolveWithCG` returns (post ret==9 accept/zero handling)
versus after the `nuf`/`nub` post-bottom smoothing sweeps that run before `x`
ever reaches `interpolation()`.

Combined run with `ALAMO_ML_BICGSTAB_RET_DIAG`, `ALAMO_ML_BOTTOM_SMOOTH_DIAG`,
and `ALAMO_ML_COR_DIAG` all enabled, same uniform 2048x2048
`max_coarsening_level=1` default-`bicgstab` case, `max_iter=10`:

- **Bottom smoothing (`n_smooth=8`) does essentially nothing**:
  `ratio_post_over_pre` is `0.99991`-`1.00004` on every call. The smoother is
  not amplifying or meaningfully damping; `crse_after_bottom` (seen by
  `interpolation()`) is, for practical purposes, pure bicgstab/CG output.
- Growth trace across outer iterations is now fully bracketed and
  self-consistent: `crse_after_bottom` norminf `8.89e-7 -> 5.64e-5 -> 2.01e-3
  -> 8.65e-2` (~40x-60x per iteration); `fine_after_interp_sync` tracks it
  almost 1:1 (interpolation itself is not where the blowup happens); `Fine
  resid/bnorm` grows in lockstep `138.9 -> 4922 -> 263376 -> ...`.
- **Iteration 1 alone is already broken**: starting from a correction with
  `norminf~8.87e-7` (genuinely tiny, not itself anomalous), the very first
  V-cycle already ends at `resid/bnorm=138.9`. Every later iteration is the
  same unstable map re-applied to an already-larger input, not new behavior
  building up gradually.

**Conclusion: the explosion is fully localized to one already-identified
operation** — applying the fine-level operator to the freshly-interpolated
correction during iteration 1 (the `A*x norminf=1.02e11` on `x` with
`norminf~8.79e-7` from the earlier `ALAMO_ML_RESID_DIAG` probe). `Fapply` was
independently proven correct on an analytic linear field at both `mglev=0` and
`mglev=1`. Bottom smoothing is now proven inert. Interpolation does not
amplify (`crse_after_bottom` -> `fine_after_interp_sync` is ~1:1). So the
defect is not in any single kernel being "wrong" — it is that bicgstab's
bottom correction, despite a locally sane convergence ratio in its own residual
norm, is spatially rough enough (relative to this operator's
stiffness/`dx^2` scale) that the genuinely-correct, genuinely-stiff fine
operator turns a tiny correction into a residual 138x larger than bnorm. This
points at the *coarse-grid representation/consistency* of the correction
(Galerkin/operator-coarsening quality, or the coarse operator's coefficients)
rather than at any individual transfer or smoothing kernel.

Ruled out so far: Fapply (both levels), bottom-smoothing amplification,
interpolation/restriction amplification, bicgstab false-convergence
(`ret==0`), and the `ret==9` partial-acceptance exemption is unconfirmed
pending a redo on the correct tree (the original test was void, see
CORRECTION section above).

## Next Target

Do not broaden the sweep yet. The next diagnostic target is bottom `bicgstab`
failure handling and correction acceptance policy:

- inspect the AMReX bottom `bicgstab` path to determine what state is returned
  after `MLMG: Bottom solve failed` and whether ALAMO/AMReX should accept that
  correction for coarse-to-fine interpolation;
- identify whether the failed bottom solve can be configured to return a bounded
  correction like `cg`/`smoother`, reject the failed correction, or fall back to a
  safer bottom solver for this nodal elastic operator;
- keep `ALAMO_ML_ZERO_CRSE_GHOSTS=1` and `ALAMO_ML_SPATIAL_DIAG=N` as
  negative-control/diagnostic knobs, not fixes;
- only revisit boundary or ghost-read hypotheses if a later run shows the
  correction acceptance policy is clean.

## Diagnostic Update: 2026-06-23 Bottom-Smoothing Reproduction — closes the "n==0 ghost-staleness" sub-thread from GPU_ELASTIC_DEBUG_PLAN.md's 2026-06-23 code-review pass

`GPU_ELASTIC_DEBUG_PLAN.md`'s same-day code-review pass on the ghost-fill path
re-flagged `interpolation()`'s unguarded raw `Array4` ghost read as a live
suspect, conditioned on whether `n == (ret==0) ? nub : nuf` could be `0` for
this case (which would mean the bottom-solve producer leaves `x`'s ghosts
unrefreshed before `interpolation()` reads them). That question is
**answered: no.** Reran `ALAMO_ML_BOTTOM_SMOOTH_DIAG=20` on the standard
uniform 2048² `max_coarsening_level=1` default-bicgstab case
(`log_coarsen1_2048_uniform_i12_bicg_bottomsmoothdiag_cc.log`, `max_iter=12`,
independent rerun from the 2026-06-22 `i10` run above): every one of 20 calls
reports `ret=9` (bottom solve "failed" but accepted) and therefore always
takes the `nuf` branch with `n_smooth=8` — never `0`. `ratio_post_over_pre`
stays in `0.999919`-`1.000054` across all 20 calls (essentially reproducing
the 2026-06-22 `0.99991`-`1.00004` range to within rerun noise), confirming
(a) the post-bottom smoother is inert (not an amplifier, as already
established) and (b) it always runs with `n=8`, so per the already-established
"`Fsmooth` always re-fills/syncs ghosts every call" fact, `interpolation()` is
**not** reading stale/unrefreshed ghosts via this mechanism for the default
configuration. The `n==0` ghost-staleness scenario the code-review pass raised
is a configuration that doesn't occur on this failure path; do not re-open it
without a concrete case where `ret==0` and `nub` is small.

This does **not** resolve the separately-observed `grown_contains_nan=1` /
`grown_contains_inf=1` finding in the coarse correction's ghost region
(`GPU_ELASTIC_DEBUG_PLAN.md`'s "Code-review pass" section) — that is a claim
about ghost *validity* (NaN/Inf contamination), not ghost *staleness* or
*amplitude*, and is a distinct, still-open question. Given the post-smooth
ratio≈1 result here, if that NaN/Inf observation is confirmed real (not an
artifact of the diagnostic itself), it would have to be present already in the
pre-smooth state too (since smoothing barely changes the field) — worth
checking directly with `compute-sanitizer --tool initcheck` on this exact case
per `GPU_ELASTIC_DEBUG_PLAN.md` Priority 2, targeting the bottom-solve-to-
interpolation handoff specifically, rather than re-deriving it from amplitude
diagnostics.

Standing next step is unchanged from the 2026-06-22 "Next Target" above:
bottom `bicgstab` failure handling / correction acceptance policy is still the
most direct unexplained mechanism (a tiny, non-anomalous correction at outer
iteration 1 already produces `resid/bnorm=138.9` on the genuinely-correct fine
operator). This run's full growth trace (`pre_smooth_norminf` across the 10
outer iterations, 2 calls/iteration): `2.1e-6, 8.8e-7 | 7.0e-5, 5.9e-5 | 5.1e-3,
3.9e-3 | 0.13, 0.11 | 5.7, 4.5 | 240, 222 | 7418, 5523 | 2.38e5, 2.00e5 |
4.35e6, 3.68e6 | 3.95e8, 3.26e8` — geometric blowup (~30-100x per outer
iteration) starting immediately, consistent with and reproducing the
2026-06-22 trace shape.

## Diagnostic Update: 2026-06-23 CPU Newton-Stall Full Exploration

Followed up the "CPU-side `SetUniform(true)`" entry above by fully exploring the
low-void-stiffness Newton stall (the "Side Thread" in `GPU_ELASTIC_DEBUG_PLAN.md`,
where the why/what-to-do summary now lives). Harness:
`scratchpad/run_stall.sh` runs the exact stall config `/tmp/test_setuniform_input`
(double_circle annular grain, `use_psi=1`/`psi_floor=0.05`, fixed elastic load
`traction=1.0_MPa` with `traction_from_chamber=0`, `interval=25`, `nriters=200`,
`nrtolerance=1e-5`, MLMG `max_iter=200`) on CPU `bin/alamo-2d-g++`, varying one
knob per run. Clean regime: 1 NR iter/solve, ~15-19 MLMG V-cycles to
resid/bnorm~5e-9.

| Run | Override (vs baseline) | Stall onset (step / t) | Max NR | MLMG in stall | Outcome |
|---|---|---|---|---|---|
| A | baseline (void 0.2, default coarsening) | 5200 / 0.520 | 200 | Final Iter 122, resid/resid0=9.2e-9 (reaches 1e-8) | STALL, exhausts nriters, limit cycle ~2.19e-3 |
| B | `max_coarsening_level=0` | 5200 / 0.520 | 200 | Final Iter 49, resid/resid0=9.5e-9 | STALL, IDENTICAL limit cycle ~2.19e-3 (just cheaper iters) |
| C | `bottom_solver=cg` | 5200 / 0.520 | 200 | Final Iter 127-130, resid/resid0=9.3e-9 | STALL, identical to default bicgstab |
| D | `model_void.kappa=mu=1_MPa` | 5350 / 0.535 | 200 | — | STALL, later onset |
| E | `model_void.kappa=mu=2_MPa` | 5475 / 0.5475 | 200 | — | STALL, later onset |
| F | `model_void.kappa=4 mu=3_MPa`, stop=0.55 | none <=5500 | 2 | n/a | CLEAN to t=0.55 (FALSE NEGATIVE — see Fext) |
| Fext | `model_void.kappa=4 mu=3_MPa`, stop=0.65 | 5600 / 0.560 | 200 | — | STALL (same limit-cycle signature) |

### Key conclusions
1. **Open question resolved.** A stalled solve EXHAUSTS `m_nriters=200` (never
   escapes, never NaNs), then `Newton::solve` returns 0.0 and the sim silently
   continues (B: hit NR 200 at step 5200, advanced through 5201-5225, started a
   new stalling solve at 5225). NR norm grows geometrically to ~2.42e-3 (iter 15)
   then dead-flat limit cycles at ~2.19e-3.
2. **Mechanism corrected: undamped-Newton overshoot, NOT inexact linear solve.**
   During the stall the MLMG linear solve still reaches `tol_rel=1e-8` (A/B/C all
   show resid/resid0 ~9e-9), and A (full MG, 122 V-cyc) and B (coarsen0, 49
   V-cyc) reach the IDENTICAL nonlinear limit cycle. The stall is independent of
   linear-solver path/quality => it is the full Newton step overshooting on the
   soft-void large-strain tangent. Refutes the side-thread's "compounding effect
   #1 (inexact linear solve)".
3. **`max_coarsening_level=0` (the previously-recommended first try) and any
   bottom-solver change are REFUTED.** Neither touches the nonlinear limit cycle.
4. **Flooring void stiffness only DELAYS onset, does NOT fix it** (+125-150 steps
   per stiffness step: 0.2->5200, 1.0->5350, 2.0->5475, 4/3->5600). The earlier
   "4/3 MPa clean to t=0.55" was a false negative (run stopped before its later
   onset); extended to t=0.65 it stalls at t=0.56. The previously-recorded
   "4-6 MPa safe floor" is therefore not a true floor in this regime.
5. **Silent-failure code cluster** (Newton.H `Set::Vector` overload + caller):
   `Newton.H:397` `return 0.0` on exhaustion == ideal-convergence return value;
   `Mechanics.H:207` discards the return; `Newton.H:378-379` commented out so the
   logged "relative norm(ddisp)" is the ABSOLUTE max increment (norm0, m) and
   `nrtolerance=1e-5` is really a 10-micron absolute step floor (the `Set::Scalar`
   overload is genuinely relative -> the two are inconsistent); `Newton.H:366`
   `if(nriter==m_nriters)break;` is dead code.

## Diagnostic Update: 2026-06-23/24 `compute-sanitizer initcheck` — concrete root-cause candidate, strongest lead to date

Ran `compute-sanitizer --tool initcheck` on `bin/alamo_gpu-2d-nofast-cuda86-g++`
against the standard minimal-failing case (uniform 2048², `max_coarsening_level=1`,
default bicgstab), `elastic.solver.max_iter=1`, `elastic.solver.bottom_max_iter=5`
(tight caps needed purely to keep sanitizer instrumentation overhead tractable —
a full-iteration run hung past 580s without producing the first MLMG print line
even at `max_iter=2`; this run still exercises the exact bottom-solve ->
normalize -> interpolation handoff at least once).
`sanitizer_initcheck_2048_coarsen1_i1_bmi5_cc.log`.

**Result: `ERROR SUMMARY: 296460 errors`, all "Uninitialized __global__ memory
read of size 8 bytes", all 10000 captured instances sharing the IDENTICAL host
call chain** (not noise/multiple unrelated sites):

```
Operator.cpp:427  Operator<Grid::Node>::normalize()
  -> MultiFab::divide (MultiFab.cpp:1387)
  -> [uninitialized read in amrex::Divide kernel, AMReX_FabArrayUtility.H:1336]
called from MLCGSolver.H:182 (solve_bicgstab)
  <- bottomSolveWithCG (MLMG.H:1699)
  <- actualBottomSolve (MLMG.H:1615)
  <- mgFcycle (MLMG.H:1485)
  <- oneIter (MLMG.H:1306)
```

**Mechanism (code-read, not yet independently confirmed by a targeted
single-call probe):**
- `Operator<Grid::Node>::normalize(amrlev, mglev, a_x)` (`Operator.cpp:419-428`)
  does `a_x.divide(*m_diag[amrlev][mglev], 0, getNComp(), 2)` — explicitly
  including **2 ghost layers** — then immediately
  `a_x.FillBoundaryAndSync(Geom(amrlev,mglev).periodicity())`. The
  FillBoundaryAndSync call strongly suggests the design intent is "the ghost
  region computed by the divide is scratch; the sync call immediately after
  fixes it up" — which is true for interior/periodic-neighbor ghosts but not
  for physical-boundary ghosts (see below).
- `m_diag[amrlev][mglev]`'s ghost cells are explicitly zeroed once
  (`diagfab.setVal(0.0)`, `Operator.cpp:315`, inside `Diagonal()`) and never
  written again — so every `normalize()` call divides by a known-zero in the
  ghost region, by construction.
- **The domain is non-periodic** (`geometry.is_periodic = 0 0 0` in
  `input_base`). `FillBoundaryAndSync(periodicity)` only exchanges data
  between periodic images / overlapping valid regions of neighboring boxes —
  it does **not** fill ghost cells that sit at the *physical* domain boundary,
  since there is no periodic donor there. Whatever the divide computes at
  those specific ghost nodes (uninitialized numerator / zero denominator)
  survives `normalize()` unchanged.
- That corrupted vector is one of BiCGStab's own internal Krylov vectors
  (`MLCGSolver.H:182`), which is fed straight back into `Fapply` on the next
  operator application — `Fapply`'s nodal stencil legitimately reads ghost
  cells near the domain boundary (that is how BC-extension layers work for a
  nodal operator), so it directly ingests the corrupted values.

**Why this is consistent with everything found so far:**
- Explains the standing `grown_contains_nan=1`/`grown_contains_inf=1`
  ghost-region finding (`GPU_ELASTIC_DEBUG_PLAN.md`) — this is *where that NaN
  comes from*, not just where it's observed.
- Explains why `cg` only stalls while `bicgstab` explodes: both solvers call
  the same `normalize()` on their own internal vectors (same defect exposure),
  but BiCGStab's two-term mixed recurrence is far more sensitive to a
  corrupted search direction than CG's.
- Explains why the smoother-only bottom-solver path is unaffected by *this*
  defect (it has its own separate, already-documented stall): `Fsmooth` never
  calls `normalize()` at all.
- Consistent with `Fapply` being independently proven correct (sine probe,
  2026-06-23): the defect is upstream of `Fapply` — in the data it's fed, not
  in the stencil evaluation itself.

**Not yet done (do this before treating this as confirmed root cause, not
just strongly-supported hypothesis):**
1. A targeted single-call probe: dump `a_x`'s ghost-cell values immediately
   before and after one `normalize()` call near a physical (non-periodic)
   domain boundary, to directly confirm NaN/Inf/garbage survives at those
   specific nodes and nowhere else (interior/periodic ghosts should be clean).
2. Check whether `normalize()`'s ghost-inclusive divide is itself necessary —
   the natural candidate fix is restricting the divide to the valid region
   only (`nghost=0`), since the ghost region is going to be
   overwritten/attempted-overwritten by `FillBoundaryAndSync` immediately
   after anyway for the cells that call actually fixes; for the physical-
   boundary cells it doesn't fix, a proper BC-aware ghost fill (not just
   `FillBoundary(periodicity)`) would be needed if those ghosts must carry a
   meaningful value at all for this operator. Do not patch and declare victory
   without rerunning the full 2048² case to confirm the explosion is gone, not
   just the sanitizer error.
3. Confirm whether AMReX's other (non-Alamo) nodal linops define `normalize()`
   the same way (ghost-inclusive divide) — if so, this would be expected to
   also affect any non-periodic nodal MLMG use case generically, which would
   make this an upstream-relevant finding, not Alamo-specific.

## Candidate fix #2 tested and REFUTED (2026-06-24): restricting normalize()'s divide to nghost=0

Implemented "not yet done" item 2 above: changed `Operator.cpp:424`'s
`a_x.divide(*m_diag[amrlev][mglev],0,getNComp(),2)` to
`a_x.divide(*m_diag[amrlev][mglev],0,getNComp(),0)` (valid-region-only divide,
no ghost layers touched). Rebuilt `bin/alamo_gpu-2d-nofast-cuda86-g++`
(confirmed via `Operator.cpp.o`/binary mtimes that the rebuild actually picked
up the change) and reran the exact pre-fix baseline case
(`coarsen1_2048_uniform_i60`'s overrides: `amr.n_cell=2048 2048
elastic.use_psi=0 elastic.max_coarsening_level=1 elastic.solver.max_iter=60`,
uniform `kappa=162_MPa mu=113.6_MPa` for prop/void/casing) via `run_case.sh`.

**Result: unchanged explosive divergence.** Diverges at the same iteration
count (12) with `resid/bnorm=1.276639791e+20` — same order of magnitude as the
pre-fix baseline's `2.624074357e+20`
(`postfix_coarsen1_2048_uniform_i60` rc=6, `rc=6 t=164s verdict=DIVERGED`).
The compute-sanitizer `initcheck` finding (296,460 uninitialized-read errors
in `normalize()`'s ghost-inclusive divide on a non-periodic domain) is real
and the analysis of *why* it occurs is presumably still correct, but it is
**not the cause of the 2048² explosion** — fixing it left the failure mode
completely unchanged. The ghost-inclusive divide may still be worth fixing
upstream/defensively (it's still touching uninitialized memory), but it is no
longer a candidate explanation for this bug.

**Patch was reverted** (back to `,2)` / ghost-inclusive divide) and the binary
rebuilt again to match, so the working tree is back at the pre-2026-06-24-fix
baseline. The investigation is back to item 1 above (the `Fapply`/`Fsmooth`
coarse-level defect hypothesis from the 2026-06-22 entries, and the still-open
"not yet done" instrumentation step there) as the live lead — the
`normalize()` ghost divide is now a closed, refuted side-branch.

**Correction:** "item 1" above (`Fapply` coarse-level correctness) was
*already cleared* on 2026-06-23 via the `ALAMO_ML_FAPPLY_PROBE_MODE=sine`
high-frequency differential (see `GPU_ELASTIC_DEBUG_PLAN.md` "CURRENT
PRIORITIES", "Priority 1 is DONE") — bit-correct to roundoff at both
`mglev=0` and `mglev=1`. Re-stating it as the live lead above was stale; see
the next entry for what was actually checked next.

## m_diag (the normalize()/Fsmooth divisor) tested and CLEARED (2026-06-24)

Identified a real, previously-unchecked gap: `ALAMO_ML_COEFF_DIAG` had
verified the *stiffness tensor* `m_ddw_mf` (used by `Fapply`) is correctly
restricted across MG levels, but `m_diag` (computed by
`Elastic<SYM>::Diagonal()`, `Elastic.cpp:779-886` — note this overrides, and
is the one actually used in place of, the unrelated host-loop checkerboard
`Diagonal()` in the generic `Operator.cpp` base class — the latter is dead
code for `Elastic` and was a red herring not worth pursuing) had never been
independently checked GPU-side. `m_diag` is what both `normalize()`'s divide
and `Fsmooth`'s Jacobi update (`x = ... / diagfab(i,j,k,n)`) divide by, so a
wrong value here would corrupt every bottom-solver variant identically
(matching the standing "all of bicgstab/cg/smoother fail the same way"
observation).

Added `ALAMO_ML_DIAG_PROBE` (same pattern as `ALAMO_ML_COEFF_DIAG`;
`Elastic.cpp`, top of file) dumping per-`(amrlev,mglev)` sum/min/max/frob_norm
plus separate nan/inf/zero counts, over both the **valid** region and the
**grown** (ghost-inclusive) region, called from `Diagonal()` both
**presync** (right after the compute loop) and **postsync** (after
`FillBoundaryAndSync`+`nodalSync`) — specifically to catch the same
non-periodic-physical-boundary-ghost mechanism already confirmed in
`normalize()`, applied here instead.

Ran on the standard `2048², max_coarsening_level=1` case
(`log_diagprobe2_2048_coarsen1.log`):

```
mglev=0 valid : min=-1.16446e+17 max=1           nan=0 inf=0 zero=0     frob_norm=3.396e+20
mglev=0 grown : min=-1.16446e+17 max=1.5e+08      nan=0 inf=0 zero=21246 frob_norm=3.495e+20 (postsync)
mglev=1 valid : min=-2.91116e+16 max=1           nan=0 inf=0 zero=0     frob_norm=4.273e+19
mglev=1 grown : min=-2.91116e+16 max=1.5e+08      nan=0 inf=0 zero=10641 frob_norm=4.520e+19 (postsync)
```

**Verdict: `m_diag` is correct, both interior and ghost. CLEARED.**
- Valid-region magnitudes are large (BC-penalty-style terms at physical
  boundary nodes, `ALAMO_ELASTIC_OP_BC_EVAL` branch — large by design for an
  essential-BC penalty, not a bug) but scale **exactly** as expected between
  levels: `mglev=1`'s min/max is `mglev=0`'s `/4`, matching the 2× DX
  coarsening of a second-derivative-based diagonal term to roundoff
  precision. Same self-consistency argument that cleared `m_ddw_mf`.
- The only "bad" entries anywhere, at any region, are **zero — never NaN or
  Inf**. These are benign: the compute loop only writes out to
  `validbox().grow(1) & domain`, so the outer (2nd) ghost layer at a physical
  boundary is genuinely never written and stays at its default-zero value,
  both before and after `FillBoundaryAndSync` (which, as established for
  `normalize()`, cannot reach non-periodic physical-boundary ghosts either).
  This is unwritten padding, not corruption.
- This also explains, retroactively, why the candidate fix #2 above (removing
  `normalize()`'s ghost-inclusive divide) had **zero effect**: there was
  never any genuine garbage (NaN/Inf) in those ghost cells to divide by in
  the first place — just legitimate zeros — so restricting the divide to
  `nghost=0` couldn't have changed anything. The two refutations are
  mutually consistent, not contradictory.

**Implication:** every Alamo-side candidate checked so far is now clean:
`Fapply` (bit-correct), `m_ddw_mf` coefficients (bit-correct restriction),
`m_diag` (correct, this entry), `Fsmooth`'s ghost-refill (structurally always
re-syncs, 2026-06-23), the `n==0` ghost-staleness path (closed, `ret=9`/`nuf=8`
always). The live leads still open, per `GPU_ELASTIC_DEBUG_PLAN.md`'s
"Next:" paragraph (just before "Original priority list"), are: (a)
`Operator<Grid::Node>::interpolation()` (`Operator.cpp:709-768`)'s raw,
unchecked `Array4` read one ghost layer into `crse` with no `FillBoundary` of
its own immediately before — still the most concrete untraced Alamo-side
transfer code; or (b) AMReX's own `bottomSolveWithCG`/`MLCGSolver.H`
internals (Saxpy/dotxy/norm on GPU) producing an already-anomalous
correction before any Alamo code touches it — increasingly looking like an
upstream, not Alamo-specific, GPU code path given how much Alamo-side code
has now been cleared.

## Sanity check (2026-06-24): is dirty `ext/AMReX-Codes/amrex` confounding any of this?

`ext/AMReX-Codes/amrex` is its own nested git repo (not a submodule of the
outer repo, gitignored by the outer repo) and was found dirty:
```
 M Src/LinearSolvers/MLMG/AMReX_MLCGSolver.H
 M Src/LinearSolvers/MLMG/AMReX_MLMG.H
```
Read both diffs in full. Both are purely additive: a `getenv`-gated block
(`ALAMO_ML_BICGSTAB_RET_DIAG`, `ALAMO_ML_BOTTOM_SMOOTH_DIAG`) that prints
diagnostics and is a complete no-op when the env var is unset — the
surrounding solver math/control-flow (`ret==0||ret==8` accept check,
`linop.smooth(...)` call) is untouched, just wrapped. Confirmed these env
vars were NOT set in any run referenced in this entry (`env | grep
ALAMO_ML` was empty). **Not a confound.**

## Baseline non-reproduction (2026-06-24): REPORT.md's 512²/1024² convergence does not currently hold

To rule out the uncommitted `Elastic.cpp`/`Operator.cpp`/`Newton.H` diagnostic
edits as the cause of an unexpectedly bad result (`coarsesweep_512_mcl1`
diverging, which the original 2026-06-21 study never tested but seemed
surprising given 512² convergence elsewhere), ran:
```
git stash push --include-untracked -- src/Operator/Elastic.cpp src/Operator/Operator.cpp src/Solver/Nonlocal/Newton.H
./configure ... (unchanged, same POSTFIX) && make
```
then reran the exact REPORT.md-style "uniform operator" config
(`use_psi=0`, prop=void=casing=162MPa/113.6MPa) at 512² with DEFAULT
(uncapped) `max_coarsening_level`, i.e. exactly the config REPORT.md says
converges in 9 iterations. Result on clean HEAD (`7d9b12336`, last commit
touching these files is `1da57132e` from 2026-06-20, predating REPORT.md):
**DIVERGED** (`cleanhead_sanity_512_default`, 29 iterations →
3.04e20). Also reran the documented canonical `coarsen1_2048_uniform_i60`
case on clean HEAD: still diverges (`cleanhead_coarsen1_2048_i60`,
8.84e28 — same failure mode, sanity-confirms clean HEAD behaves like the
working tree for the already-known-bad case). `git stash pop` to restore the
diagnostic instrumentation afterward; rebuilt again before continuing.

**Conclusion: source code is not the variable.** Both the working tree and
clean HEAD diverge at 512² where REPORT.md (3 days earlier, same source,
binaries explicitly "NOT rebuilt" per its own Method section) documented
clean 9-iteration convergence. Likely candidates not yet checked: CUDA
toolkit/driver version drift, GPU clock/thermal/ECC state, or some other
machine-level change unrelated to this branch's commits. This is a loose
end that should be investigated separately — it means the "2048² is the
resolution cliff" framing from REPORT.md may not currently hold at all (see
the coarse-size sweep immediately below, which shows divergence at every
tested resolution down to 512² once `max_coarsening_level=1`, and even at
512² with default depth).

## Part 1 experiment matrix (2026-06-24): coarse-grid-size and depth sweeps

Harness: `run_case.sh`, binary `bin/alamo_gpu-2d-nofast-cuda86-g++`, uniform
operator (`use_psi=0`, prop=void=casing=162MPa/113.6MPa), `elastic.solver.max_iter=60`
unless noted. All commands of the form:
```
run_case.sh LABEL bin/alamo_gpu-2d-nofast-cuda86-g++ 1 \
  'amr.n_cell=N N' elastic.max_coarsening_level=MCL \
  elastic.use_psi=0 model_{prop,void,casing}.kappa=162_MPa model_{prop,void,casing}.mu=113.6_MPa \
  elastic.solver.max_iter=60
```

**A. Coarse-grid-size sweep, `max_coarsening_level=1` fixed (2 levels only):**

| n_cell | verdict | iters to abort | final resid/bnorm |
|---|---|---|---|
| 512  | DIVERGED | 20 | 1.52e21 |
| 768  | DIVERGED | 15 | 1.58e21 |
| 1024 | DIVERGED | 15 | 1.08e20 |
| 1536 | DIVERGED | 14 | 2.13e21 |
| 2048 | DIVERGED | 12 | 1.83e21 |

**No cliff — every tested resolution diverges once `max_coarsening_level` is
capped at 1.** This contradicts the "resolution is the trigger" framing;
capping depth at 1 level is itself sufficient to trigger the failure
regardless of absolute size, at least for this uniform-operator
configuration. (Caveat: see baseline non-reproduction above — these numbers
are self-consistent with each other but don't match REPORT.md's older
default-depth numbers at the same resolutions.)

**B. Depth sweep, `amr.n_cell=2048 2048` fixed, default-depth + mcl=0 sanity
re-check:**

| max_coarsening_level | verdict | notes |
|---|---|---|
| 0 (smoother-only, no bottom solve) | converging, slow | `resid/bnorm` 0.0047→7.0e-5 monotonically over 24+ iters in <600s before this session's wrapper timeout cut it off; matches existing notes that mcl=0 avoids the bug entirely by avoiding the bottom-solve code path |
| 1 | DIVERGED | 12 iters → 1.83e21 (same as table A) |
| default (uncapped) | DIVERGED | 16 iters → 1.86e20 — same failure as mcl=1, not obviously different in severity |
| 2, 3 | not run | deprioritized once the bottom_tol_abs finding (below) made the bottom-solve-specific mechanism clear; mcl=0's clean (if slow) behavior already brackets "no bottom solve, no problem" vs "any bottom solve, problem" |

512² with DEFAULT depth also diverges (`sanity_512_default`/`cleanhead_sanity_512_default`,
20-29 iters → ~1-3e20) — see baseline non-reproduction note; this means even
the smallest resolution tested now fails once a real bottom solve runs,
default-depth or not.

## Part 2: root-cause via existing `bottom_verbose` knob (no new instrumentation needed) — see GPU_ELASTIC_DEBUG_PLAN.md "CURRENT PRIORITIES" for the full writeup

Discovered `src/Solver/Nonlocal/Linear.H:270-272` already wires
`elastic.solver.verbose>4` to `MLMG::setBottomVerbose`, which (with zero
source patch) turns on AMReX's own `MLCGSolverT::solve_bicgstab` per-iteration
trace (`AMReX_MLCGSolver.H:201-207,242-248`, `verbose>2` branch) — never
previously used in this investigation; every prior probe was a black-box
before/after check, never an inside-the-loop look at the bottom Krylov solve.
Ran with `elastic.solver.verbose=5 elastic.solver.max_iter=3` on the
canonical `2048²`/`mcl=1` case (`log_bicgstab_inner_trace_2048_mcl1.log`):
every one of 6 bicgstab calls across the 3 outer iterations starts from an
`error0` already at the FP noise floor (1.16e-9 down to 4.17e-11) relative to
the operator's actual coefficient scale (O(1e11)-O(1e17), per the
already-cleared `m_diag`/`m_ddw_mf` probes), runs the full 200-iteration cap
without ever satisfying a relative-tolerance break, and reports a
deceptively small final relative error (0.004-0.24) while the **fine-level**
residual explodes 131.8 -> 23133.8 -> 953962.9 over those same 3 outer
iterations. The bottom solve "looks fine" by every metric it reports on
itself; the corruption is in the absolute magnitude/direction of the
correction it hands back, not anything its own convergence check would
catch.

**Confirmation test (not just inference): `elastic.solver.bottom_tol_abs=1e-3`
eliminates the explosion.** This makes BiCGStab's existing (otherwise dead,
at these residual magnitudes) early-return check
(`AMReX_MLCGSolver.H:145`, `if (rnorm0==0 || rnorm0<eps_abs) return ret;`)
actually fire, skipping the noise-grinding entirely.
- `confirm_bottomtolabs_2048_mcl1_full` (max_iter=200): `Fine resid/bnorm`
  decreases smoothly and monotonically every iteration, 0.97 -> 0.32 over the
  full 200-iteration budget, zero instability. Run ends in
  `amrex::Abort("MLMG failed")` only because 200 iterations isn't enough
  budget to reach the requested relative tolerance from this starting point
  ("Failed to converge after 200 iterations. resid, resid/bnorm =
  10802778.5, 0.3207274092") — an ordinary insufficient-budget exit, not a
  numerical blow-up. (Note: `run_case.sh`'s `LASTRESID` regex grabs the
  absolute residual `10802778.5` from this final summary line, not the
  `0.32` relative value — a harness display quirk to be aware of when
  reading `verdict=UNKNOWN` cases; check the log directly.)
- `confirm_bottomtolabs_512_mcl1` (max_iter=60): same clean monotonic decline,
  0.72 -> 0.10, same insufficient-budget (not blow-up) exit.
- **Caveat — does not generalize to default-depth:**
  `confirm_bottomtolabs_512_default` (default uncapped `max_coarsening_level`,
  same `bottom_tol_abs=1e-3`) still DIVERGES (34 iters -> 4.57e20). Expected:
  `bottom_tol_abs` only gates the bottom-level CG/BiCGStab call; default deep
  coarsening has many additional intermediate-level smoothing/transfer steps
  that never go through that early-return at all. If the same
  near-zero-residual GPU-sensitivity mechanism generalizes to intermediate
  `Fsmooth`/restriction/interpolation (untested), it would need a different,
  more invasive mitigation than this one ParmParse knob.

**Not done (AMReX-internal per-iteration Krylov-vector instrumentation, e.g.
inside `solve_bicgstab`'s loop dumping `p`/`v`/`r`/`t` themselves):
unnecessary** — the existing `bottom_verbose` knob already exposed the
relevant per-iteration relative-error trace without any source patch, and
the `bottom_tol_abs` confirmation test already validates the mechanism
directly. Could still be done as a follow-up to see exactly which GPU
reduction (`dotxy` in `rho`/`rhTv`/`tvals`) first decorrelates from the CPU
equivalent, but is no longer necessary to establish the root cause.

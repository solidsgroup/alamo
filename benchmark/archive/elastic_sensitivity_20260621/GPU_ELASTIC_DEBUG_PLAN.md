# GPU Elastic MLMG Divergence Debug Plan

> ## ✅ SOLVED 2026-06-25 — READ THIS FIRST (supersedes everything below)
>
> **Root cause:** a GPU **cross-stream use-after-free race** on the per-box
> temporary `FArrayBox tmpfab` in `Operator<Grid::Node>::interpolation()`
> (`src/Operator/Operator.cpp`). AMReX's non-tiled GPU `MFIter` cycles iterations
> across a pool of CUDA streams; `tmpfab` was freed at iteration end and its
> device block reused by a later iteration on a **different stream** while the
> interpolation/`plus` kernels were still in flight → corrupted coarse-grid
> correction → amplified by the BC-penalty diagonal (~1e17) into the 1e20
> blow-up. `restriction()` writes directly (no temp), so only interpolation broke.
>
> **Fix (one line):** `amrex::Gpu::Elixir tmpfab_eli = tmpfab.elixir();` right
> after `tmpfab.resize(...)` in `interpolation()`.
>
> **How localized:** (1) box sweep — single box converges, any multi-box
> decomposition diverged → cross-box; (2) stock-AMReX `MLNodeLaplacian`
> reproducer converges multi-box → ALAMO custom code, not AMReX/BiCGStab;
> (3) bit-identical CPU↔GPU transfer probe → `interpolation`, not `restriction`;
> (4) pre/post-`nodalSync` bisect → the per-box temp (nondeterministic run-to-run).
> Everything below (noise-floor BiCGStab, nondeterministic reduction,
> deterministic coarse-operator defect, the "2048² cliff") was a **downstream
> symptom** of this race. The `interpolation()` "live lead" was the right spot.
>
> **Verified end-to-end** (strict/no-fast-math binary): full box sweep
> (mgs 2048..128 + deep coarsening) all converge ~1e-9; real production operator
> (use_psi=1, 3-material) multi-box 2048² converges 7 iters → 8.5e-9; end-to-end
> flame+elastic chamber sim (forced multi-box, 300 steps) clean, 10 elastic
> solves all <1e-8, sane physics. The text below is kept for investigation
> history only.

Date: 2026-06-21

## CURRENT PRIORITIES (2026-06-24, latest update — READ FIRST)

**Strongest lead to date, with a validated workaround: BiCGStab's bottom
solve is being handed a residual already at the double-precision noise floor
(~1e-9 to 1e-11 absolute) and grinds it for the full 200-iteration cap with
no absolute-tolerance escape hatch — this is what explodes the V-cycle, not a
ghost/race/coefficient defect.**

Root-cause method: `elastic.solver.verbose=5` (>4) propagates to
`MLMG::setBottomVerbose`, which (with zero source changes, just an existing
ParmParse knob already wired in `src/Solver/Nonlocal/Linear.H:270-272`) turns
on AMReX's own per-iteration BiCGStab trace
(`MLCGSolverT::solve_bicgstab`'s `verbose>2` branch, `ext/AMReX-Codes/amrex/
Src/LinearSolvers/MLMG/AMReX_MLCGSolver.H:201-207,242-248`) — this had never
been used in this investigation; every prior probe only checked Alamo-side
state before/after the bottom solve as a black box. On the canonical
minimal-failing case (2048², `max_coarsening_level=1`, uniform operator,
`use_psi=0`), this trace shows:
```
MLMG: Initial residual (resid0) = 33682118.18
MLCGSolver_BiCGStab: Initial error (error0) = 1.156999015e-09   <- already ~roundoff
MLCGSolver_BiCGStab: Final: Iteration 201 rel. err. 0.0334       <- looks "fine" by its OWN metric
MLMG: Iteration 1 Fine resid/bnorm = 131.8     <- but the FINE level already exploded
MLMG: Iteration 2 Fine resid/bnorm = 23133.8   <- and keeps exploding, ~180x per outer iter
MLMG: Iteration 3 Fine resid/bnorm = 953962.9
```
Every one of the 6 bicgstab calls in this 3-outer-iteration run starts from
an absurdly tiny `error0` (1.16e-9, 4.17e-11, 1.37e-7, 3.10e-8, 2.68e-5,
3.17e-6) and *all six* hit the 200-iteration cap (`Final: Iteration 201`)
without ever satisfying a relative-tolerance break — yet each reports a
small-looking final relative error (0.004–0.24). The bottom solve's own
diagnostics never flag a problem; the damage is being done in absolute terms
while looking converged in relative terms, then injected into the fine level
via interpolation.

**Mechanism (consistent with every earlier finding):** when the restricted
coarse residual is already at the FP noise floor relative to the operator's
actual dynamic range (the BC-penalty diagonal entries are O(1e11)–O(1e17) per
the already-cleared `m_diag`/`m_ddw_mf` probes), BiCGStab's `rho`/`alpha`/
`omega` updates (`AMReX_MLCGSolver.H:158,182-185,218-226`) divide by
near-zero dot products. GPU vs. CPU reduction-order non-associativity (a
completely ordinary, well-known FP effect — not a bug) gets catastrophically
amplified by these near-zero denominators, producing a numerically
uncorrelated-with-the-true-solution correction that nonetheless "passes" the
relative-error check because everything is being measured relative to an
already-microscopic baseline. This explains: why `Fapply`/coefficients/
`m_diag` are all independently correct (none of them are wrong — the inputs
to the bottom solve are fine, it's iterating on noise that was never worth
solving); why CPU doesn't show this (different, evidently more benign,
reduction order — or simply less sensitive to it at this conditioning); why
no input lever (psi_floor, void stiffness, dt, fast-math) fixed it (none of
those change whether the bottom-level residual is noise-floor-tiny); and
why it is deterministic per-run on GPU yet differs from CPU (fixed GPU
reduction order, different fixed CPU reduction order — not a race).

**Validated workaround (confirmation test, not just a hypothesis check):**
setting `elastic.solver.bottom_tol_abs=1e-3` (so BiCGStab's existing
early-return check at `AMReX_MLCGSolver.H:145`, `if (rnorm0==0 || rnorm0 <
eps_abs) return ret;`, actually fires instead of being dead code at these
tiny `rnorm0` values) **completely eliminates the explosion** at 2048²/mcl=1:
`Fine resid/bnorm` now decreases smoothly and monotonically every single
outer iteration (0.97 → 0.32 over the full 200-iteration cap, see
`log_confirm_bottomtolabs_2048_mcl1_full.log`) — it hits `amrex::Abort`
only because 200 iterations isn't enough to reach the requested tolerance
from this starting point (an ordinary "didn't converge in budget" exit, not
a numerical blow-up: `MLMG: Failed to converge after 200 iterations. resid,
resid/bnorm = 10802778.5, 0.3207274092`). Same clean-decline result at
512²/mcl=1 (0.72 → 0.10 over 60 iters, `log_confirm_bottomtolabs_512_mcl1.log`).
**This is the first manipulation in the whole investigation that turns the
explosion into ordinary slow convergence rather than just changing how long
it limps before exploding.**

**Important caveat — does NOT generalize to default/deep coarsening:**
the same `bottom_tol_abs=1e-3` does **not** fix 512² with `max_coarsening_level`
left at its default (uncapped) value — still explodes
(`log_confirm_bottomtolabs_512_default.log`, DIVERGE → 4.57e20). This makes
sense: `bottom_tol_abs` only gates AMReX's *bottom* CG/BiCGStab call.
Default/deep coarsening has many additional intermediate MG levels whose
smoothing/restriction/interpolation never go through that bottom-solver
early-return at all — if the same "near-zero residual + GPU reduction-order
sensitivity" mechanism generalizes to those intermediate-level operations
(not yet confirmed, but the obvious next hypothesis), it would need a
different mitigation (e.g. a residual floor check inside `Fsmooth`/
`Operator<Grid::Node>` generic smoothing, not just the bottom solver).

**Also discovered, unrelated to the above and NOT YET EXPLAINED — flag for a
fresh investigation, do not assume it's stale:** the 2026-06-21 `REPORT.md`
baseline ("uniform operator converges cleanly at 512²/1024², only 2048²
diverges") **does not currently reproduce**, even after (a) `git stash`-ing
all uncommitted `Elastic.cpp`/`Operator.cpp`/`Newton.H` changes back to clean
HEAD (`7d9b12336`) and rebuilding, confirming the Alamo source is
byte-identical to what's been unmodified since commit `1da57132e` (2026-06-20,
predates the REPORT.md run), and (b) confirming the only dirty state in the
vendored `ext/AMReX-Codes/amrex` tree (`AMReX_MLCGSolver.H`,
`AMReX_MLMG.H`) is purely additive, env-var-gated `getenv`/`Print` diagnostic
code with the actual solver math/control-flow completely unchanged (verified
by reading the diffs — see `fix_notes.md`'s "Sanity check" entry below for
the full diff). With identical Alamo source and behaviorally-inert AMReX
diffs, `cleanhead_sanity_512_default` and even `sanity_512_default` (current
working tree) both **diverge** at 512² with default (uncapped) coarsening —
contradicting REPORT.md's documented `conv 9 it`/`conv 22+22` results for
that exact configuration class. Source code is not the variable here, so
suspect environment drift (CUDA toolkit/driver version, GPU clock/thermal/ECC
state) between 2026-06-21 and now, not a regression in this branch's commits.
This needs its own investigation before treating any absolute resolution
threshold (e.g. "2048² is the cliff") as still authoritative — it may no
longer be a clean threshold at all on the current machine state, see the new
`coarsesweep_*`/`sanity_*` logs in this directory (512² through 2048² *all*
diverge under the current environment with `max_coarsening_level=1`, and
512² also now diverges even at default depth).

**Not yet done (deprioritized in favor of the finding above, still worth
doing):** (1) confirm whether the noise-floor mechanism generalizes to
intermediate-level `Fsmooth`/restriction/interpolation under default/deep
coarsening (would explain the default-depth-512² counterexample above); (2)
box-tiling/`amr.max_grid_size` sweep at the 2048²/mcl=1 case (planned, not
run — lower priority now that a mechanism + workaround exists that doesn't
depend on box decomposition at all); (3) root-cause the REPORT.md
non-reproduction (separate from this bug). Raw logs for everything above:
`coarsesweep_{512,768,1024,1536,2048}_mcl1`, `sanity_{512,2048}_default`,
`cleanhead_*`, `depthsweep_2048_mcl{0,2,3}` (mcl0/2/3 launched but not all
completed within this session's time budget — mcl0 is slow-but-clean per
existing notes, mcl2/3 not yet run), `bicgstab_inner_trace_2048_mcl1`,
`confirm_bottomtolabs_*`.

---

## Goal

Find the first incorrect GPU operation responsible for the 2048² elastic MLMG
divergence.

The target is not to prove that the solve diverges; that is already established.
The target is to localize the failure to a specific multigrid phase, level,
kernel family, or runtime condition so it can be fixed or reported upstream with
a minimal reproducer.

## CURRENT PRIORITIES (2026-06-23/24 update — READ FIRST; supersedes the "narrowed target" notes scattered below)

**2026-06-24 LATEST — m_diag (the normalize()/Fsmooth divisor) tested and
CLEARED.** Added `ALAMO_ML_DIAG_PROBE` (mirrors `ALAMO_ML_COEFF_DIAG`) to dump
`m_diag` min/max/frob_norm/nan/inf/zero, valid vs grown region, pre- and
post-sync, in `Elastic<SYM>::Diagonal()`. Result on the standard 2048² case:
no NaN/Inf anywhere; valid-region magnitudes are large (BC-penalty terms at
physical boundaries) but scale exactly by the expected `/4` between mglev=0
and mglev=1; the only "bad" ghost entries are benign unwritten zero padding
(2nd ghost layer at physical boundaries, never written by design, same
FillBoundaryAndSync-can't-reach-non-periodic-ghosts limitation as
`normalize()` — but with nothing actually garbage there). `m_diag` is
correct. Combined with the already-cleared `Fapply` (2026-06-23) and
`m_ddw_mf` coefficients, essentially all Alamo-side coarse-level state is now
verified clean. Full writeup: `fix_notes.md`, "m_diag tested and CLEARED".
**Remaining live leads:** (a) `Operator<Grid::Node>::interpolation()`
(`Operator.cpp:709-768`)'s unchecked raw `Array4` ghost read — still the most
concrete untraced Alamo-side code; (b) AMReX's own `bottomSolveWithCG`
internals (Saxpy/dotxy/norm) — increasingly the more likely culprit given
how much Alamo-side code has been ruled out.

**2026-06-24 LATE — normalize()-divide candidate fix tested and REFUTED.**
Patched `Operator.cpp:424` to `,0)` (valid-region-only divide, no ghosts) per
item 2 of the sanitizer finding's "not yet done" list below, rebuilt, and
reran the standard `coarsen1_2048_uniform_i60` repro. Result: identical
explosive divergence (12 iterations, `resid/bnorm=1.28e20`, same order as the
`2.62e20` pre-fix baseline). The sanitizer's uninitialized-read finding is
real but **is not the cause of the 2048² explosion**. Patch reverted, binary
rebuilt back to baseline. Full writeup: `fix_notes.md`, "Candidate fix #2
tested and REFUTED" entry. The live lead is back to the pre-existing item 1
(`Fapply`/`Fsmooth` coarse-level instrumentation, not yet done) — do not
re-attempt the ghost-divide fix without new evidence.

**2026-06-24 — concrete root-cause candidate found via Priority 2
(`compute-sanitizer initcheck`), strongest lead to date. Read
`fix_notes.md`'s "compute-sanitizer initcheck — concrete root-cause
candidate" entry in full before doing anything else.** Summary: every one of
296,460 flagged errors (10000 captured instances, all identical call chain)
is an uninitialized `__global__` memory read inside
`Operator<Grid::Node>::normalize()` (`Operator.cpp:419-428`)'s
`a_x.divide(*m_diag[...], 0, getNComp(), 2)` — a divide that includes 2 ghost
layers, called from AMReX's internal bicgstab (`MLCGSolver.H:182`) on its own
Krylov vectors. `m_diag`'s ghosts are deliberately zeroed and never re-synced
(divide-by-zero in the ghost region, every call, by construction); the
domain is non-periodic (`geometry.is_periodic = 0 0 0`), so the
`FillBoundaryAndSync(periodicity)` immediately following the divide does NOT
clean up the physical-boundary ghost cells (no periodic donor exists there) —
only interior/periodic ghosts get fixed. The corrupted vector is fed straight
into the next `Fapply` call, which legitimately reads boundary ghost cells as
part of the nodal stencil. This explains the standing
`grown_contains_nan=1`/`grown_contains_inf=1` ghost-region finding (shows
*where* the NaN originates), why `cg` only stalls while `bicgstab` explodes
(same defect, different sensitivity to a corrupted search direction), and why
the smoother-only path is unaffected by this specific defect (`Fsmooth` never
calls `normalize()`). NOT yet independently confirmed by a targeted
before/after ghost-value probe at a physical boundary, and NOT yet patched —
see `fix_notes.md` "Not yet done" list (single-call probe, candidate fix of
restricting the divide to valid-region-only, and whether this is an
Alamo-specific `normalize()` override or shared with other AMReX nodal
linops). Do this before resuming Priority 3/5 below.

**Priority 1 is DONE as of 2026-06-23 — `Fapply` is now genuinely cleared at
both levels.** `ALAMO_ML_FAPPLY_PROBE_MODE=sine` was added to
`src/Operator/Elastic.cpp`'s existing `ALAMO_ML_FAPPLY_PROBE` diagnostic: it
overrides displacement component 0 with `amp*sin(ax*i)*sin(ay*j)` (wavelength
~4 cells, `amp=1e-6`, all other components zero) and computes the matching
closed-form expected `f` analytically in the same kernel (uniformly-sampled
sinusoids are exact eigenfunctions of the standard central-difference
second-derivative and mixed-derivative stencils used here, so this is an exact
ground truth, not an approximation, when the operator is uniform/psi-disabled
— this harness's reproduction baseline). This is the actual discriminating
test the WARNING below used to call for; the old constant/linear probe could
not see a near-null-space-invisible defect, this field is not in the
near-null-space.

Result, GPU, `bin/alamo_gpu-2d-nofast-cuda86-g++`, exact documented minimal-
failing case (2048², `elastic.max_coarsening_level=1`, uniform
prop=void=casing, `use_psi=0`):
- `mglev=0` (fine): `interior_max_abs=1.16e11`, `|F-expected|=0.0193`,
  **relative=1.66e-13** (roundoff).
- `mglev=1` (coarse, the suspect level): `interior_max_abs=2.91e10`,
  `|F-expected|=1.80e-3`, **relative=6.20e-14** (roundoff).
- Sanity control at a known-good small case (512², same `max_coarsening_level=1`):
  `mglev=1` relative=2.33e-14 (roundoff, as expected).

**Conclusion: `Fapply` is bit-correct on GPU at the coarse level, even under a
genuinely high-frequency rough input, at the exact resolution/hierarchy where
the full multi-iteration solve explodes.** The defect is NOT in per-node
stencil evaluation at either level. This re-narrows the standing hypothesis to
the *other* half of the "coarse-level operator application" framing: the
generic `Operator<Grid::Node>` nodal transfer (restriction/interpolation) path,
ghost-fill/sync around the bottom solve, or AMReX's own bottom-`bicgstab`
internal vector ops (AXPY/dot/norm) on GPU — not `Elastic::Fapply` itself.
Reproduce via:
```
ALAMO_ML_FAPPLY_PROBE=1 ALAMO_ML_FAPPLY_PROBE_MODE=sine \
benchmark/elastic_sensitivity_20260621/run_case.sh LABEL \
  bin/alamo_gpu-2d-nofast-cuda86-g++ 1 \
  'amr.n_cell=2048 2048' elastic.use_psi=0 elastic.max_coarsening_level=1 \
  model_prop.kappa=162_MPa model_prop.mu=113.6_MPa \
  model_void.kappa=162_MPa model_void.mu=113.6_MPa \
  model_casing.kappa=162_MPa model_casing.mu=113.6_MPa
```
Note: building this required reconfiguring the live build from the in-progress
Phase 3 3D-fast-math setup to 2D-no-fast-math
(`./configure --comp=g++ --dim=2 --cuda local --cuda-fp strict && make`,
matching the recipe already in `fix_notes.md`'s "Rebuilt
`bin/alamo_gpu-2d-nofast-cuda86-g++`" entries) — reconfigure back to 3D before
resuming Phase 3 work.

### Code-review pass on the ghost-fill path (2026-06-23, before running the sanitizer)

Traced every ghost-fill/transfer call on the V-cycle's down-sweep -> bottom ->
up-sweep path that touches the bottom-level correction, since that is exactly
what `interpolation()` reads:

- `Operator<Grid::Node>::restriction()` (`src/Operator/Operator.cpp:525-648`)
  reads `fine`'s ghost cells (stencil reaches `i±1` in coarse index = up to
  `±1` fine cell beyond its own box) and ends with
  `crse.FillBoundary(...); nodalSync(...)`. The `fine` it reads
  (`rescor[amrlev][mglev]`) is produced by `correctionResidual()`
  (`Operator.cpp:1001-1058`), which itself ends with
  `resid.FillBoundaryAndSync(...)` before `restriction` is called from
  `MLMG.H:1373/1383`. This leg looks consistent — ghosts are filled
  immediately before they're read.
- `Operator<Grid::Node>::Fsmooth()` (`Operator.cpp:350-417`) runs 2 internal
  Jacobi sub-iterations per call and **calls `x.FillBoundary(...);
  nodalSync(...)` at the end of both**, including the last one. So contrary to
  an initial hypothesis, AMReX's generic `smooth(..., skip_fillboundary, niter)`
  wrapper (`AMReX_MLNodeLinOp.cpp:169-179`) skipping `applyBC` before the
  *first* of its `niter` calls to `Fsmooth` does not leave stale ghosts after
  the *last* one — Alamo's `Fsmooth` re-fills and re-syncs every time
  regardless. This refutes a "smoother leaves stale ghosts" theory at the
  `Fsmooth` level specifically.
- **`Operator<Grid::Node>::interpolation()` (`Operator.cpp:709-768`) is the
  one place in this leg with no FillBoundary of its own immediately before
  the read.** It pulls `crsefab = (*cmf)[mfi]` and indexes
  `cdata(I, J, K, n)` / `cdata(I+1, J, K, n)` / etc. directly off the raw
  `Array4` with no bounds check against the box's *valid* region — at every
  fine node on a box edge, `I+1` (or `J+1`/`K+1`) lands in `crse`'s *ghost*
  cell, one layer deep. It trusts that whoever last wrote `crse` (the bottom
  solve) left that ghost layer correctly filled and synced. For
  `bottom_solver=smoother`, that producer is `actualBottomSolve()`
  (`AMReX_MLMG.H:1564-1568`) calling the same `smooth()` analyzed above — so
  by the `Fsmooth`-always-refills finding, this particular leg's ghosts
  *should* be fresh. For `bicgstab`/`cg`, the producer is
  `bottomSolveWithCG()` (not yet traced) followed by `n` post-smooth
  iterations (`MLMG.H:1675`, `n` = `nub` or `nuf`) — if `n==0` for some
  configuration, `x`'s ghosts at the moment `interpolation()` reads them are
  whatever `bottomSolveWithCG()` itself left, which has **not** been audited
  for whether it leaves a valid/synced ghost layer. This is the most concrete
  remaining ghost-fill suspect: an `n==0` (or `bottomSolveWithCG`-internal)
  path that returns a correction without a final ghost sync, consumed
  unguarded by `interpolation()`'s raw `Array4` read.
- This matches the standing instrumentation finding already on record
  (`grown_contains_nan=1` / `grown_contains_inf=1` in the coarse correction's
  ghost region) and explains *why* `ALAMO_ML_ZERO_CRSE_GHOSTS=1` didn't fix the
  blowup: zeroing a bad ghost still feeds the stencil a wrong (just
  differently-wrong) value — the fix has to be filling it correctly, not
  blanking it.
- Next code-level step if the sanitizer doesn't localize this directly: trace
  `bottomSolveWithCG()` (`AMReX_MLMG.H`, not yet read) for whether it leaves
  `x`'s ghost cells filled/synced on return, and check the actual `nub`/`nuf`
  values this harness's `elastic.solver.*` inputs resolve to (zero would make
  this the live path, not just a theoretical one).

**UPDATE 2026-06-23 (later) — the `n==0` ghost-staleness scenario above is
CLOSED, narrowing this further.** Reran `ALAMO_ML_BOTTOM_SMOOTH_DIAG=20` on the
standard default-bicgstab case (`fix_notes.md`, "Bottom-Smoothing
Reproduction" entry): every call has `ret=9`, so `n` always takes the `nuf=8`
branch — never `0`. Combined with the already-established "`Fsmooth` always
re-fills/syncs ghosts" fact, this means `interpolation()` is NOT reading
stale/unrefreshed ghosts via an `n==0` skip on this failure path; that specific
ghost-staleness mechanism does not occur here and should not be re-investigated
without a concrete `ret==0`/small-`nub` case. `ratio_post_over_pre≈1.0` on
every call also reconfirms (this was first shown 2026-06-22) that the
post-bottom smoother is not the amplifier — whatever inflates the correction
is already present when `bottomSolveWithCG` returns, before any Alamo
transfer/smoothing code runs on it.

This does **NOT** clear the separate `grown_contains_nan=1`/`grown_contains_inf=1`
ghost-region finding below, which is about ghost *validity* (NaN/Inf), not
staleness or amplitude — that remains open and is the right target for the
sanitizer run.

Next: Priority 3 (minimal AMReX-only nodal-MLMG reproducer) and Priority 5
(box-size sweep) both target the transfer/ghost path now implicated; Priority 2
(`compute-sanitizer`, confirmed present at
`.local/cuda-12.0-ubuntu/usr/bin/compute-sanitizer` — the "not found" note
below is stale) is still the fastest path to ground truth and should run on
the restriction/interpolation kernels in `src/Operator/Operator.cpp`'s
`Operator<Grid::Node>` transfer path specifically, not `Elastic.cpp` — and
should specifically target the NaN/Inf ghost-region claim now that the
staleness/amplitude sub-questions are closed, since AMReX's own
`bottomSolveWithCG`/BiCGStab (not Alamo's transfer code) is now the leading
suspect for where an already-anomalous correction originates.

---

### Original priority list (2026-06-22, Priority 1 now superseded above)

The investigation has localized the failure well but then **circled**: the most
recent `fix_notes.md` entry exonerates `Fapply` using a test that cannot detect
the suspected defect, and redirects at a cause already ruled out. A fresh agent
should NOT resume the bicgstab-acceptance-policy or coefficient threads. Work the
priorities below in order.

### Corrected root-cause status (what is actually ruled OUT)
- **BiCGStab acceptance policy: REFUTED.** The ret==9 patch, redone on the
  CORRECT AMReX tree (`ext/AMReX-Codes/amrex`, log `log_ret9patch_..._i60`,
  2026-06-22 13:15), still diverges 148 → 1.75e20. bicgstab never falsely
  converges (`raw_ret=8`, local `rnorm/rnorm0` 0.002–0.28 every call).
- **Nondeterministic GPU reduction: REFUTED as the bottom-level cause.**
  `bottom_solver=smoother` does zero reductions/dot-products/atomics yet still
  STALLS on GPU (resid/bnorm 0.74) while converging on CPU (0.18) for the
  bit-identical restricted system. The defect is deterministic.
- **Coarse coefficients: REFUTED.** `ALAMO_ML_COEFF_DIAG` shows the uniform-field
  restriction is bit-identical across MG levels.

### Standing hypothesis (corrected, then re-narrowed 2026-06-23)
A **deterministic, GPU-specific defect in the coarse-level (mglev≥1) operator
application** — `Elastic<SYM>::Fapply` and/or the generic `Operator<Grid::Node>`
nodal transfer/ghost handling — that produces a wrong stencil result on
**rough/high-frequency** inputs but is invisible on smooth ones. The two observed
symptoms (explosive bicgstab, stalling cg/smoother) are almost certainly the SAME
coarse-operator defect under different solvers; treat them as one problem.

> **UPDATE 2026-06-23 — `Fapply` IS now cleared, with a test that can actually
> see the defect.** The WARNING immediately below was correct that the old
> constant/linear probe couldn't detect a rough-input defect. A new
> `ALAMO_ML_FAPPLY_PROBE_MODE=sine` high-frequency closed-form probe (see
> "CURRENT PRIORITIES" above) re-ran the exact same test the WARNING demanded,
> at the exact documented minimal-failing hierarchy, and found `Fapply`
> correct to roundoff (relative diff ~1e-13/1e-14) at both `mglev=0` and
> `mglev=1`. The defect is in the transfer/ghost path or AMReX's bottom solver,
> not `Fapply`. The paragraph below is kept for history; do not re-open the
> "Fapply might still be the cause" thread without new evidence.

> **(historical) WARNING — do not trust the `Fapply` "exoneration" in `fix_notes.md`.** That
> probe used only constant and linear displacement fields, which lie in the
> operator's near-null-space (exact answer ≈ 0). It cannot see a defect that only
> manifests on rough inputs — which is exactly the failing case. `Fapply` is NOT
> cleared.

### Priority 1 — High-frequency `Fapply` differential (THE discriminating test)
Feed `Fapply` at `mglev=1` a known high-frequency analytic field (e.g.
`u = sin(kx)·sin(ky)` scaled to ~1e-6 m, where `A·u` has a closed form) and diff
**GPU vs CPU** at that level, per-node. Extend `ALAMO_ML_FAPPLY_PROBE` in
`src/Operator/Elastic.cpp` to support a sinusoidal field, not just constant/linear.
This confirms or kills the standing hypothesis; everything else is downstream.

### Priority 2 — Unblock the sanitizer (likely the fastest ground truth)
Phase 2 was filed "blocked" without checking the bundled toolkit. First run
`find .local -name 'compute-sanitizer*'` (it ships with the CUDA toolkit). If
present, run `--tool initcheck` then `--tool memcheck` on the two-level
`max_coarsening_level=1` case. The coarse correction already shows
`grown_contains_nan=1 / grown_contains_inf=1` in its ghost region — `initcheck`
is built to localize exactly that (uninitialized/stale device-ghost read). Note:
`ALAMO_ML_ZERO_CRSE_GHOSTS=1` does NOT exonerate ghost-handling — zeroing ghosts
still feeds the stencil wrong values, so the blowup persisting under it is
expected, not a clean negative.

### Priority 3 — Minimal AMReX reproducer (Phase 6.2 — bisects the biggest fork)
Build a standalone uniform-coefficient nodal-MLMG case (no ALAMO physics). If it
diverges on GPU → upstream AMReX bug report. If it converges → the defect is in
ALAMO's custom `Operator<Grid::Node>` transfer / `Elastic::Fapply`, not AMReX.
Nothing else resolves "ALAMO operator vs upstream AMReX" as cleanly, and it has
not been attempted.

### Priority 4 — Quantify the avoidance path (`max_coarsening_level=0`)
`max_coarsening_level=0` (single MG level, no bottom solve) already CONVERGES at
2048² on GPU. Complete the wall-clock-vs-CPU comparison (CPU ref 1401.7 s) that
was listed but never finished. If the goal is literally "elastic runs on GPU"
(the 2026-06-22 D1 reversal), this avoidance path may already satisfy it — see
the strategy fork below — and de-risks the whole track independent of the fix.

### Priority 5 — Box-size sweep (Phase 1.2, now diagnostic — was skipped early)
Box decomposition controls how the coarse operator's per-box owned/ghost regions
are assembled, so it is directly relevant to the current ghost/coarse-operator
hypothesis. It was skipped "because a trigger was found," but it is now MORE
diagnostic, not less. Run it if Priorities 1–3 do not localize the defect.

### Strategy fork (flag for the user / lead; do not silently pick)
- **(A) Avoidance** (`max_coarsening_level=0`): makes GPU elastic run today, no
  fix, at the cost of MG efficiency. Needs only the Priority-4 perf number.
- **(B) Root-cause-and-fix** (Priorities 1–3, 5): required for a correct
  general-purpose GPU MLMG.

The `docs/gpu_elastic_device_port_plan.md` "three hard gates" (no-fast-math
convergence at high res, fair N-rank CPU comparison, ncu on A100) gate (B). The
divergence hunt here is what blocks that plan's Gate 1.

---

## Side Thread: CPU-Side Newton Stall at Low Void Stiffness (2026-06-22)

Separate from the GPU MLMG divergence above, but related (same elastic
operator, same high-contrast-coefficient pathology). Raw run data is in
`fix_notes.md` ("CPU-side `SetUniform(true)` diagnostic flag" and the
"Newton-stall full exploration" sections). This is the why/what-to-do writeup.

### FULLY EXPLORED 2026-06-23 — headline conclusions (supersedes the original
### "potential solutions / recommended first try" notes kept below for history)

1. **Reproduced and characterized.** Config `/tmp/test_setuniform_input`
   (double_circle annular grain; `model_void.kappa=mu=0.2_MPa`;
   `model_prop`/`model_casing`=162/113.6 MPa => contrast 810:1; `use_psi=1`,
   `psi_floor=0.05`; fixed elastic load `traction=1.0_MPa`,
   `traction_from_chamber=0`; `elastic.interval=25`; `nriters=200`,
   `nrtolerance=1e-5`; MLMG `max_iter=200`, `pre/post_smooth=4`). CPU
   `bin/alamo-2d-g++`. Solves are clean (1 NR iter, ~15-19 MLMG V-cycles to
   resid/bnorm~5e-9) through step 5175 / t=0.5175; the solve at step 5200
   (t=0.52) stalls. The trigger is GEOMETRY EVOLUTION at fixed load (the
   annular void grows/changes shape as the grain burns back), not load growth.

2. **The open question is resolved: the stall EXHAUSTS `m_nriters` (200).** It
   neither escapes nor NaNs. NR relative-norm grows geometrically to a peak
   (~2.42e-3 at iter 15) then settles into a dead-flat limit cycle at
   ~2.19e-3 and stays there to iteration 200. On exhaustion `Newton::solve`
   falls through to `return 0.0` and the simulation SILENTLY CONTINUES to the
   next step and the next (also-stalling) solve. Verified empirically: run B
   hit NR iter 200 on the step-5200 solve, then advanced through steps
   5201-5225 and began a new solve at step 5225, with no NaN/abort/error.

3. **Mechanism CORRECTED: this is a pure undamped-Newton overshoot
   (nonlinear), NOT an inexact-linear-solve problem.** During the stall the
   MLMG linear solve STILL reaches `tol_rel=1e-8` ("Final Iter. 122,
   resid/resid0=9.2e-9") — it is exact, merely ~6x costlier than the clean ~19
   V-cycles. Decisive evidence: the full-multigrid path (run A, 122 V-cycles)
   and the single-level `max_coarsening_level=0` path (run B, 49 V-cycles)
   reach the IDENTICAL nonlinear limit cycle (~2.19e-3). The limit cycle is
   therefore independent of linear-solver path/quality. This REFUTES the
   original "compounding effect #1 (inexact linear solve)" below; only
   "compounding effect #2 (soft void => large strain => the per-step tangent
   linearization is a poor local model => the full Newton step overshoots)"
   stands. There is no line search / damping / under-relaxation in
   `Newton.H` — the full step is always applied (`Newton.H:388-389`).

4. **The recommended first fix `max_coarsening_level=0` is REFUTED on CPU.**
   It does not break the limit cycle (run B reaches the same ~2.19e-3); it only
   makes each Newton iteration cheaper so the run burns through `nriters`
   faster. (Note this is the OPPOSITE of its effect on the GPU divergence,
   where the coarse level is *defective* so removing it helps. On CPU the
   coarse level is correct and removing it gives up useful MG acceleration for
   no convergence benefit.)

5. **Tightening the linear solve cannot help, confirmed.** `bottom_solver=cg`
   (run C) stalls at the IDENTICAL step 5200/t=0.52 as default bicgstab, with the
   MLMG still reaching `tol_rel=1e-8` (Final Iter 127-130, resid/resid0=9.3e-9).
   No linear-solver/bottom-solver knob touches the nonlinear limit cycle.

6. **Flooring void stiffness only DELAYS the stall — it does not prevent it.**
   The stall ONSET marches later in the burn as the void stiffens (the evolving
   annular geometry has to burn further before conditioning crosses the
   overshoot threshold), at ~+125-150 steps per stiffness step:
   - 0.2 MPa: stall onset step 5200 (t=0.520)
   - 1.0 MPa: stall onset step 5350 (t=0.535)   [run D]
   - 2.0 MPa: stall onset step 5475 (t=0.5475)  [run E]
   - 4.0/3.0 MPa: stall onset step 5600 (t=0.560) [run Fext]
   All reach NR iter 200 (exhaust) at their onset, same limit-cycle signature.
   The earlier "4/3 MPa clean to t=0.55" was a FALSE NEGATIVE: it just hadn't
   reached its (later) onset; extending the run to t=0.65 stalls it at t=0.56.
   So the previously-recorded "4-6 MPa safe floor" is NOT a true floor in this
   regime — for any finite void softness the evolving annular geometry
   eventually crosses the overshoot threshold. Flooring stiffness only buys
   time; the real fix is on the nonlinear (Newton) side.

### Silent-failure code cluster (fix regardless of the numerical remedy chosen)
All in `src/Solver/Nonlocal/Newton.H` `solve(Set::Field<Set::Vector>...)`
(the overload Flame uses) and its caller `Mechanics.H:207`:
- `Newton.H:397` `return 0.0;` on `m_nriters` exhaustion returns the SAME value
  an ideal zero-residual convergence returns — a stall is actively misreported
  as perfect convergence.
- `Mechanics.H:207` discards the return value entirely, so even that signal is
  unused; a stalled timestep is indistinguishable from success except for the
  climbing "NR iteration N" log lines.
- `Newton.H:378-379` are commented out, so `solnorm` stays 0 and `relnorm =
  cornorm`: the value logged as "relative norm(ddisp)" is actually the ABSOLUTE
  max nodal increment (`norm0`, meters). `nrtolerance=1e-5` is thus a 10-micron
  absolute step floor, not a relative tolerance. The `Set::Scalar` overload
  (`Newton.H:458-459`) DOES compute `solnorm` and is genuinely relative — the
  two overloads are inconsistent.
- `Newton.H:366` `if (nriter == m_nriters) break;` is dead code (unreachable
  inside `for (nriter=0; nriter<m_nriters; ...)`).

### Recommended remedies, ranked (corrected by this exploration)
1. **Damped / line-search Newton** (code change; THE fix). The A-vs-B-vs-C
   identical-limit-cycle-via-different-linear-paths result is the textbook
   signature of full-step overshoot, which a residual-decrease backtracking line
   search (accept a step only if it reduces `compResidual`'s norm; else halve)
   or simple under-relaxation directly targets. Required for cases needing a
   genuinely soft/near-fluid void, since flooring (below) cannot make any finite
   softness safe over a full burn. Not implemented here (needs a CPU rebuild;
   the shared `.make` build config is currently pointed at the GPU target, and a
   regression test run was using `bin/alamo-2d-g++` during this work).
2. **Make the `m_nriters` cap a hard, visible failure** (small code change; do
   regardless of the numerical fix). Replace the silent `return 0.0` with a
   non-convergence signal and have `Mechanics.H:207` act on it (warn / reduce
   dt / abort). Today a stall is indistinguishable from success in the log.
3. **Floor `model_void` stiffness** (input-only) — MITIGATION ONLY, not a fix.
   Stiffer void delays onset (~+125-150 steps per stiffness step) but does not
   prevent it (4/3 MPa stalls at t=0.56). Use it to push the stall past the
   simulation's end time if the physics tolerates a stiffer void, but do not
   treat any value as permanently safe.
4. **Continuation / load-stepping** (code change): ramp void stiffness or load
   over sub-steps so each Newton solve starts near its solution.
5. **True void exclusion (cut-cell / embedded boundary)** (large change):
   removes the coefficient jump at its source.

Levers proven NOT to help (do not retry): `max_coarsening_level=0` (run B), any
bottom-solver change (`bottom_solver=cg`, run C), and stiffness flooring as a
permanent fix (runs D/E/Fext). The stall is nonlinear, so no linear-side knob
and no finite stiffness floor resolves it.

### Original analysis (RETAINED FOR HISTORY — partially superseded above)

### What was found
- `src/Integrator/Base/Mechanics.H:194` had an uncommitted `SetUniform(true)`
  diagnostic flag left in, which drops the `grad(C)·grad(u)` term in
  `Elastic.cpp:419` — wrong whenever the stiffness field is nonuniform (i.e.
  always, for a masked solid/void/casing operator). Reverted to
  `SetUniform(false)`; validated clean (0 MLMG failures) at
  `model_void.kappa=model_void.mu=4_MPa/3_MPa`, `stop_time=0.05s`.
- Pushing void stiffness down to `0.2_MPa` (vs. the previously-known-safe
  floor of ~4-6 MPa from `elastic_stability_sweep.md`) at `stop_time=1.5s`
  reproduces a **genuine Newton stall**, not a crash: NR relative-norm grows
  geometrically (~1.3-1.5x/iter) from `1.3e-5` up through iteration ~15, then
  flattens into a tight limit cycle around `2.1e-3`-`2.4e-3` and stays there
  indefinitely. `nrtolerance=1e-5` is never approached. Each iteration costs a
  full MLMG solve that itself gets much more expensive while stalling (~150
  V-cycles/~5s vs. ~18/~0.5s when converging cleanly). Manually killed at
  iteration ~33-34 of 200 once the stall pattern was confirmed; the run was
  never allowed to reach the iteration cap, so whether it eventually escapes,
  NaNs, or just exhausts `m_nriters` was not observed.

### Why this is likely happening
Two compounding effects, not one:

1. **Stiffness contrast makes each Newton step's own linear solve inexact.**
   At `0.2_MPa` void vs. `162_MPa` solid/casing, the contrast is ~810:1 (vs.
   ~40:1 at the known-safe 4-6 MPa floor). MLMG cost rising from ~18 to ~150
   V-cycles per step is the standard multigrid signature of a coefficient
   jump that geometric coarsening/smoothing doesn't resolve well — feeding an
   under-converged linear solve into Newton breaks the usual quadratic
   convergence.
2. **The void region is soft enough to be genuinely large-strain.** At
   `0.2_MPa` the void deforms far more than at `4_MPa` for the same load,
   pushing the local response into the part of the NeoHookean curve where the
   per-step tangent linearization is a poor local model. Each Newton step then
   overcorrects rather than converging — consistent with the observed
   growth-then-plateau (bounce, not blowup) signature.

Both are consistent with `rod_and_tube` already being flagged as the
"binding"/stall geometry in `elastic_stability_sweep.md` — this is the same
known failure mode, reproduced now that `SetUniform(false)` is correctly
restored (i.e. `SetUniform(false)` is necessary but not sufficient at this
stiffness contrast).

### Potential solutions, roughly by effort/risk
- **Floor `model_void` stiffness** (cheapest): don't let it drop below
  ~1-2 MPa if the physics allows it. Punts on cases that need a near-fluid
  void.
- **Tighten the linear solve for this regime**: more bottom-solve iterations,
  a direct/robust bottom solver instead of coarsening through the
  high-contrast field, or `max_coarsening_level=0` (same avoidance knob
  already known to work around the *GPU* divergence above — worth checking
  whether it also stabilizes this CPU stall, since both are coefficient-jump
  pathologies on the same operator).
- **Damped/line-search Newton**: scale the step when it would increase the
  residual, directly targeting the observed overshoot-and-bounce.
- **Continuation/load-stepping**: ramp void stiffness (or load) onto the soft
  value over several pseudo-steps within a timestep instead of jumping
  straight from a converged stiffer state to the fully-soft one.
- **Make the `m_nriters` cap a hard failure**: `Newton.H` currently breaks
  silently past `m_nriters` and continues with an unconverged displacement
  field rather than aborting or flagging the timestep — worth fixing
  regardless of which numerical fix is chosen, since a stall is otherwise
  indistinguishable from quiet success in the log.
- **True void exclusion (cut-cell/embedded boundary)** instead of soft-material
  masking: removes the coefficient jump at its source rather than tuning
  around it. Biggest change, addresses the root cause directly.

Recommended first try: `max_coarsening_level=0` on the CPU build, since the
mechanism is already understood from the GPU side and it's a one-line input
change with no code risk. Damped Newton is the next thing to reach for if
that's insufficient.

---

## Current Facts

- The failing case is a single static elastic solve, not a time-integration
  failure.
- Timestep changes do not affect the first-solve residual or divergence.
- Fast math is not the cause; fast and no-fast builds both fail at 2048².
- Conditioning is not the cause; a genuinely uniform operator converges on CPU
  but diverges on GPU at 2048².
- GPU convergence is normal at 512² and 1024².
- GPU failure appears at 2048², where multigrid depth and GPU concurrency are
  larger.
- Identical GPU commands can produce different V-cycle trajectories and abort
  iterations, indicating nondeterministic GPU behavior.

## Progress Update: 2026-06-21

Phase 1 has localized the smallest failing hierarchy.

- `elastic.max_coarsening_level=0` creates one MG level. The uniform 2048² GPU
  case did not diverge before the 900 s harness timeout; residual decreased
  monotonically to `resid/bnorm = 3.371003866e-05`.
- `elastic.max_coarsening_level=1` creates two MG levels and is sufficient to
  trigger the explosive failure. With the default bottom solver, the uniform
  2048² GPU case diverged after 12 outer MLMG iterations with
  `resid/bnorm = 2.624074357e+20`.
- The failure is therefore not a deep-hierarchy-only issue. It appears at the
  first fine/coarse transition and bottom-solve path.
- `compute-sanitizer` and `cuda-memcheck` are not installed on the current
  workstation PATH or under the checked CUDA locations, so Phase 2 is currently
  blocked locally unless the tools are installed.

Additional bottom-solver comparison after rebuilding
`bin/alamo_gpu-2d-nofast-cuda86-g++`:

- default `bicgstab`: explosive divergence, `resid/bnorm` reaches `1e20+`;
- default `bicgstab` with `elastic.solver.verbose=5`: the first V-cycle has a
  finite bottom correction residual norm (`UP: Norm after bottom =
  2.906744463e+05`), but coarse-to-fine correction immediately amplifies it to
  `UP: Norm before smooth = 1.025065514e+11` on the fine level and then fine
  `resid/bnorm = 128.1298774`. By outer iteration 12 the run fails with
  `resid/bnorm = 4.118904137e+20`.
- `elastic.solver.bottom_solver=cg`: no explosive divergence, but reaches the
  60-iteration cap at `resid/bnorm = 0.09265647885`;
- `elastic.solver.bottom_solver=cg` with `elastic.solver.verbose=5` and
  `elastic.solver.max_iter=12`: bottom CG also reaches its internal
  201-iteration cap, but remains non-explosive. First V-cycle bottom
  `UP = 2.858193165e+07`, fine `resid/bnorm = 1.925044027`; at iteration 12,
  fine `resid/bnorm = 0.2374600709`.
- `elastic.solver.bottom_solver=smoother`: no explosive divergence, but reaches
  the 60-iteration cap at `resid/bnorm = 0.7369654712`.
- `elastic.solver.bottom_solver=smoother` with `elastic.solver.verbose=5` and
  `elastic.solver.max_iter=12`: first V-cycle bottom `UP = 2.564384674e+07`,
  fine `resid/bnorm = 2.647172004`; at iteration 12, fine `resid/bnorm =
  3.031934833`. This remains bounded compared with `bicgstab`.
- `elastic.solver.bottom_tol_rel=0.05` for `bicgstab` does not fix the
  instability; it still diverges by 12 outer iterations with final
  `resid/bnorm = 4.361997343e+09`.
- `elastic.solver.final_smooth=0` for `bicgstab` does not fix the instability;
  it still diverges after 12 outer iterations with final `resid/bnorm =
  1.076040848e+20`.
- Do not run multiple 2048² GPU diagnostics concurrently on the local RTX A1000;
  concurrent `cg`/`smoother` tests exhausted the AMReX arena. Use serial runs.

Code insight found during inspection:

- `src/Operator/Operator.cpp` custom nodal `restriction` and `interpolation`
  used 3D stencil logic unconditionally in a 2D build, including `k - 1` /
  `k + 1` accesses in transfer stencils. A patch under test adds explicit
  `AMREX_SPACEDIM == 2` branches matching 2D full-weighting restriction and
  bilinear nodal interpolation.
- That patch changes the default `bicgstab` trajectory but does not remove the
  two-level divergence. It remains a correctness fix candidate, not a complete
  fix for the current failure.

Current narrowed target:

- instrument the two-level case at `elastic.max_coarsening_level=1`;
- compare default `bicgstab` against `cg` and `smoother`;
- the first bad amplification now appears after a finite bottom `bicgstab`
  correction is returned and interpolated/corrected back to the fine level.
- Temporary instrumentation in `Operator<Grid::Node>::interpolation` guarded by
  `ALAMO_ML_COR_DIAG=N` logs coarse correction and fine correction field norms for
  the first `N` interpolation calls. The first `bicgstab` diagnostic run
  (`coarsen1_2048_uniform_i12_bicg_cordiag1`) recorded:
  - AMReX bottom residual norm: `UP: Norm after bottom = 2.906744463e+05`;
  - coarse correction field after bottom: `norminf = 8.80301e-07`,
    `l2 = 5.32357e-04`;
  - fine correction before interpolation: `norminf = 2.38418e-10`,
    `l2 = 6.80096e-08`;
  - fine correction after interpolation and nodal sync: `norminf = 8.78737e-07`,
    `l2 = 6.20160e-04`;
  - AMReX fine residual norm immediately before smoothing:
    `UP: Norm before smooth = 1.023226618e+11`;
  - first fine residual ratio: `resid/bnorm = 129.5263407`;
  - failure after 12 iterations: `resid/bnorm = 1.668333904e+20`.
- The correction field itself does not show explosive magnitude growth during
  interpolation/sync. The new narrowed target is the fine-level operator/residual
  application after the interpolated correction is applied, or a mismatch between
  the correction representation and the operator stencil/ghost state.
- `cg` and `smoother` comparison runs with the same instrumentation
  (`coarsen1_2048_uniform_i12_cg_cordiag1` and
  `coarsen1_2048_uniform_i12_smoother_cordiag1`) produced nearly identical
  first-call correction norms: coarse `l2 = 1.89467e-06`, fine-before
  `l2 = 4.79534e-07`, and fine-after/sync `l2 = 2.34079e-06` for `cg` or
  `2.85432e-06` for `smoother`. Their fine residual ratios stayed O(1), ending
  at `1.502793872` for `cg` and `1.254453504` for `smoother` at iteration 12.
- Compared with those bounded cases, default `bicgstab` returns a first coarse
  correction about 281 times larger in L2 (`5.32357e-04` versus `1.89467e-06`)
  and immediately produces fine `resid/bnorm = 129.5263407`. This points to the
  bottom `bicgstab` correction quality/representation and its interaction with
  the fine operator, not to interpolation magnitude blow-up alone.

Spatial diagnostic update (2026-06-22):

- `ALAMO_ML_SPATIAL_DIAG=N` now records per-component extrema locations,
  domain-boundary/box-boundary/interior maxima, boundary absolute-value
  fractions, and nearest-neighbor jump metrics at the existing correction and
  residual checkpoints.
- `coarsen1_2048_uniform_i2_bicg_spatialdiag` versus
  `coarsen1_2048_uniform_i2_cg_spatialdiag` shows the first bottom-correction
  extrema are interior in both paths, not physical-boundary-only and not solely
  box-edge-localized.
- The discriminant remains correction quality/magnitude: `bicgstab` returns an
  interior valid correction about 80x larger in norminf and hundreds of times
  larger in L2 than `cg`. Interpolation creates O(max) nearest-neighbor jumps in
  both paths, but only the larger failed-`bicgstab` correction maps to explosive
  fine residuals (`resid/bnorm = 129.808` at iteration 1 versus `1.640` for
  `cg`).
- Next target: bottom `bicgstab` failure policy and the validity of reusing its
  returned correction after `MLMG: Bottom solve failed`; boundary/ghost reads are
  now negative-control checks, not the primary hypothesis.

## Primary Hypothesis

> **REVISED 2026-06-22 — see CURRENT PRIORITIES at top.** The original "GPU
> concurrency / nondeterministic reduction" framing below is SUPERSEDED: the
> smoother-only path (zero reductions) is deterministic and still fails on GPU,
> so the defect is a deterministic coarse-level operator-application bug, not a
> race/reduction. The text below is retained as original hypothesis history.

The failure is a GPU concurrency defect in the nodal MLMG path: race, stale
ghost read, missing synchronization, invalid memory access, uninitialized read,
or unsafe atomic/reduction behavior.

Resolution is likely acting as a proxy for one or more of:

- deeper multigrid hierarchy,
- larger number of boxes/tiles,
- different coarse-level geometry,
- higher kernel concurrency,
- tile-boundary or ghost-cell interaction.

## Reproduction Baseline

Use the existing harness:

```bash
benchmark/elastic_sensitivity_20260621/run_case.sh LABEL BIN NP OVERRIDE...
```

Recommended baseline case:

```bash
benchmark/elastic_sensitivity_20260621/run_case.sh \
  nofast_2048_uniform \
  bin/alamo_gpu-2d-nofast-cuda86-g++ \
  1 \
  'amr.n_cell=2048 2048' \
  elastic.use_psi=0 \
  model_prop.kappa=162_MPa \
  model_prop.mu=113.6_MPa \
  model_void.kappa=162_MPa \
  model_void.mu=113.6_MPa \
  model_casing.kappa=162_MPa \
  model_casing.mu=113.6_MPa
```

The essential requirement is a uniform elastic operator: `use_psi=0` and equal
propellant, void, and casing stiffness.

## Phase 1: Runtime Knob Localization

Goal: determine whether the failure is tied to multigrid depth, box geometry, or
general fine-level GPU execution.

### 1.1 Max-Coarsening Sweep

Run the 2048² uniform GPU case with:

```text
elastic.max_coarsening_level=0
elastic.max_coarsening_level=1
elastic.max_coarsening_level=2
elastic.max_coarsening_level=3
elastic.max_coarsening_level=4
elastic.max_coarsening_level=5
elastic.max_coarsening_level=6
elastic.max_coarsening_level=7
elastic.max_coarsening_level=-1
```

Record for each run:

- converged/diverged,
- final residual,
- first 10 `resid/bnorm` values,
- abort iteration if divergent,
- number of AMR/MG levels printed by MLMG.

Decision logic:

- If shallow coarsening converges but deeper coarsening diverges, focus on
  restriction, prolongation, coarse solve, or bottom transition.
- If `max_coarsening_level=0` still diverges, focus on fine-level operator
  application, smoother, residual, nodal BC handling, or ghost updates.
- If only a specific level count fails, build a minimal reproducer with that
  hierarchy.

Status:

- Stop broad coarsening sweeps for now. The stop condition has been met:
  `max_coarsening_level=1` is the smallest failing hierarchy observed so far.
- Use `elastic.solver.max_iter=60` for follow-up classification runs.
- Prioritize the first coarse transition and bottom-solve path over deeper
  hierarchy characterization.

### 1.2 Box Size / Decomposition Sweep

Run the 2048² uniform GPU case with several box sizes:

```text
amr.max_grid_size=32
amr.max_grid_size=64
amr.max_grid_size=128
amr.max_grid_size=256
amr.max_grid_size=512
amr.max_grid_size=1024
```

Record:

- converged/diverged,
- number of boxes per level,
- first 10 residual ratios,
- abort iteration.

Decision logic:

- Strong dependence on box size suggests tile-boundary, ghost-fill, launch
  geometry, or cross-block accumulation issues.
- No dependence on box size suggests a level/depth issue or an internal kernel
  race independent of decomposition.

### 1.3 Smoother Count Sweep

Run with:

```text
elastic.solver.pre_smooth=1 elastic.solver.post_smooth=1
elastic.solver.pre_smooth=2 elastic.solver.post_smooth=2
elastic.solver.pre_smooth=4 elastic.solver.post_smooth=4
elastic.solver.pre_smooth=8 elastic.solver.post_smooth=8
```

Decision logic:

- If more smoothing delays but does not prevent divergence, the smoother may be
  amplifying an earlier bad correction.
- If smoother count changes the first bad iteration dramatically, instrument the
  smoother path first.

## Phase 2: Sanitizer Evidence

Goal: get source-level evidence of memory/race/synchronization errors.

Use the smallest failing uniform case from Phase 1. Keep iteration limits low
enough for sanitizer runtime to be practical.

Current local status:

- `compute-sanitizer` was not found after sourcing
  `benchmark/local_cuda_env.sh`.
- `cuda-memcheck` was not found after sourcing `benchmark/local_cuda_env.sh`.
- No matching sanitizer binaries were found under `/usr/local` or `/opt`.
- If these tools are installed later, run them on the two-level
  `max_coarsening_level=1` default-`bicgstab` case first.

### 2.1 Memcheck

```bash
compute-sanitizer --tool memcheck \
  bin/alamo_gpu-2d-nofast-cuda86-g++ INPUT_AND_OVERRIDES
```

Looks for:

- out-of-bounds global memory access,
- invalid shared/local memory access,
- misaligned access,
- illegal address.

Success criterion:

- A source line, kernel name, or stack trace identifying the first invalid
  access.

### 2.2 Racecheck

```bash
compute-sanitizer --tool racecheck \
  bin/alamo_gpu-2d-nofast-cuda86-g++ INPUT_AND_OVERRIDES
```

Looks for:

- shared-memory races,
- unsafe intra-block communication.

Important limitation:

- Racecheck does not catch all global-memory races or all atomic-ordering bugs.
  A clean racecheck run does not exonerate the code.

### 2.3 Initcheck

```bash
compute-sanitizer --tool initcheck \
  bin/alamo_gpu-2d-nofast-cuda86-g++ INPUT_AND_OVERRIDES
```

Looks for:

- uninitialized device-memory reads.

### 2.4 Synccheck

```bash
compute-sanitizer --tool synccheck \
  bin/alamo_gpu-2d-nofast-cuda86-g++ INPUT_AND_OVERRIDES
```

Looks for:

- invalid barrier usage,
- divergent synchronization,
- warp-level synchronization misuse.

## Phase 3: Debug Build

Goal: make sanitizer and profiler output actionable.

Build a no-fast GPU binary with:

```text
-lineinfo
--fmad=false
```

Avoid full `-G` initially unless needed. Full CUDA debug mode can change
scheduling enough to hide races.

If sanitizer output is not source-resolved, add a second debug binary with:

```text
-G
-O0 or reduced optimization
```

Record for every binary:

- git commit,
- AMReX commit/version,
- CUDA version,
- NVIDIA driver version,
- GPU model,
- host compiler version,
- complete CUDA flags,
- complete AMReX flags.

## Phase 4: Multigrid Phase Instrumentation

Goal: identify the first phase and level where GPU state becomes wrong.

This is now the primary next phase on the current workstation because sanitizer
tools are unavailable locally.

Instrument the MLMG solve around these phases if accessible:

- before pre-smooth,
- after pre-smooth,
- after residual computation,
- after restriction,
- before and after coarse/bottom solve,
- after prolongation/correction,
- after post-smooth,
- at end of V-cycle.

For each level and field, compute:

- `min`,
- `max`,
- `norminf`,
- `L2 norm`,
- NaN count,
- Inf count,
- optionally a deterministic checksum.

Fields of interest:

- solution/displacement,
- RHS,
- residual,
- correction/error,
- operator coefficients/model,
- psi if enabled, though the uniform case should keep psi disabled.

Success criterion:

- Identify the first phase where one of these diagnostics jumps, becomes
  nondeterministic, or differs strongly from CPU.

Immediate instrumentation target:

- case: uniform 2048² GPU, `elastic.max_coarsening_level=1`,
  `elastic.solver.max_iter=60`;
- variants: default `bicgstab`, `bottom_solver=cg`, and
  `bottom_solver=smoother`;
- levels: fine `mglev=0` and coarse/bottom `mglev=1`;
- priority checkpoints: after bottom solve on coarse `cor[0][1]`, immediately
  after coarse-to-fine interpolation/correction on fine `cor[0][0]`, and after
  post-smoothing.

Expected diagnostic value:

- If `bicgstab`, `cg`, and `smoother` have matching pre-bottom state but only
  `bicgstab` returns explosive corrections, inspect AMReX bottom `bicgstab` and
  ALAMO's operator behavior on the bottom level.
- If states differ before the bottom solve, inspect restriction, coarse
  operator coefficient transfer, boundary masks, or bottom-level BC setup.
- If bottom state is finite but the fine correction explodes after
  interpolation, inspect ALAMO's custom nodal interpolation/correction path.

## Phase 5: CPU/GPU Differential Trace

Goal: separate "wrong algorithmic state" from "different but valid floating point
path."

Run the same uniform 2048² case on CPU and GPU. Dump phase diagnostics from Phase
4 for both.

Compare by:

- multigrid level,
- phase,
- outer iteration,
- field.

Decision logic:

- If CPU/GPU match until a specific GPU phase, inspect that phase's kernel.
- If GPU traces differ run-to-run before they differ from CPU in norm, inspect
  nondeterministic reductions or atomics.
- If CPU and GPU differ immediately in coefficients/RHS, the bug is in setup,
  ghost fill, BC initialization, or model transfer, not the solver iteration.

## Phase 6: AMReX-Focused Isolation

Goal: determine whether the bug is in ALAMO setup or upstream AMReX nodal MLMG.

### 6.1 Confirm ALAMO Inputs to MLMG

Before `solver.solve(...)`, verify:

- RHS norms are finite and identical CPU/GPU where expected,
- coefficients/model fields are finite,
- boundary condition arrays are initialized,
- nodal ghost cells have expected values,
- no stale data exists outside valid boxes if kernels read grown boxes.

### 6.2 Minimal AMReX Reproducer

If the uniform 2048² case still fails with simple coefficients, build a minimal
AMReX-only reproducer:

- same geometry,
- same nodal operator type,
- same BC type,
- same multigrid depth,
- uniform coefficients,
- synthetic RHS matching the ALAMO norm scale if needed.

Success criterion:

- Reproducer fails without ALAMO physics. This makes the issue suitable for an
  AMReX bug report or upstream patch.

If the AMReX-only reproducer does not fail, the issue is likely in ALAMO's
operator setup, boundary-condition wrapper, field initialization, or GPU-captured
data used by the operator.

## Phase 7: Profiler Use

Goal: map failing phases to concrete kernels and launch geometry.

Use Nsight Systems first if available:

```bash
nsys profile --trace=cuda,nvtx ...
```

Then use Nsight Compute on narrowed kernels:

```bash
ncu --set basic --kernel-name REGEX ...
```

Useful information:

- kernel names around the first bad phase,
- grid/block dimensions,
- shared memory usage,
- atomics,
- achieved occupancy,
- unusual replay or memory warnings.

Profiler output alone is not proof of correctness bugs, but it helps map phase
instrumentation to real kernels.

## Expected Fix Categories

The final fix is likely one of:

- add or move a GPU synchronization between dependent kernels,
- fix ghost-fill ordering before grown-box reads,
- replace unsafe accumulation with deterministic or correctly atomic
  accumulation,
- fix nodal tile overlap writes,
- fix coarse-level boundary or mask initialization,
- fix dimension-specific nodal transfer stencils in ALAMO's custom 2D MLMG
  restriction/interpolation path,
- avoid or fix the default GPU bottom `bicgstab` path for this custom nodal
  elastic operator,
- avoid an AMReX GPU nodal MLMG path with a known upstream defect,
- update AMReX if the issue is already fixed upstream.

## Deliverables

Produce these artifacts as the investigation proceeds:

1. `coarsening_sweep.md`
   Table of convergence versus `elastic.max_coarsening_level`.

2. `box_sweep.md`
   Table of convergence versus `amr.max_grid_size`.

3. `sanitizer_results.md`
   Sanitizer commands, logs, and first actionable findings.

4. `phase_trace.md`
   First-bad-phase diagnostics with CPU/GPU comparison.

5. `minimal_reproducer/`
   AMReX-only reproducer if the issue localizes upstream.

6. `fix_notes.md`
   Patch hypothesis, validation matrix, and before/after results.

## Stop Conditions

Stop broad sweeping once any of these occurs:

- sanitizer reports a concrete invalid access, race, uninitialized read, or sync
  error;
- phase instrumentation identifies a single first-bad phase and level;
- box-size sweep identifies a tile/box geometry trigger;
- coarsening sweep identifies a specific MG-depth trigger;
- a minimal AMReX reproducer fails.

At that point, switch from characterization to code inspection and patching in
the narrowed component.


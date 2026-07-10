# FlowSource1D cutoff + source-term investigation — status as of 2026-07-10

This file documents the full investigation into why `cutoff` (the parameter that
suppresses/relaxes the fluid solver toward the solid state in near-solid regions)
misbehaves when combined with the mass/momentum/energy source terms and the
`lagrange` no-penetration penalty, after the EOS/Gas-model merge (commit
`44a995177`, PR #274). Work is paused here to be picked up later — see "Where to
pick this back up" at the bottom.

Branch: `test/mass-source-test`. Nothing in this document has been committed yet
except the two fixes noted as COMMITTED below.

## Background / already-fixed, already-committed issues

These were fixed and committed earlier in this effort and are NOT part of the
open problem below:

1. **`Gas::ComputeLocalFractions` 0/0 at vacuum cells** (commit `f6bbd289e`).
   At `eta -> 0` (exactly), local density/moles were exactly 0, causing
   `mass_fraction = density/0 = nan`. Fixed by defaulting to a uniform
   `1/nspecies` composition when density/moles are non-positive.
2. **Unfloored velocity divide in `Hydro::RHS`** (commit `e1844bab3`).
   `v = Mx_fluid/density_fluid` was unfloored (old code had `+small`); restored
   the floor to avoid inf/nan velocities leaking out of the `eta->0` region.

With just these two fixes, the three FlowSource1D tests (mass/momentum/energy
source term) pass with `cutoff` inactive (default, or `cutoff=0.0`, or
`lagrange=0.0`). The remaining work below is entirely about making `cutoff`
(state suppression in solid regions) coexist with active source terms.

## The open problem (as posed by the user)

The user wants to use a `cutoff` value up to `0.5` (as they do in other,
non-source-term simulations) *together with* active mass/momentum/energy source
terms and the `lagrange` no-penetration penalty, without crashing or producing
wrong physics. This currently does not work for the momentum-source test case at
any `cutoff > ~0.01`.

## Chronology of root causes found (all confirmed with reproducible evidence)

### 1. Lagrange stiffness (resolved by the user, not by code change)
`lagrange=150` (as introduced in this branch) combined with the sharp
`eps=0.05` interface gives an effective spring stiffness that violates the fixed
explicit timestep's stability limit. **User found `lagrange=10.0` fixes this.**
No code change; purely an input-file tuning issue. `tests/FlowSource1D/input`
currently has (or had) `lagrange = 10.0` for testing, though the working-copy
input was reset to `lagrange = 0.0` by a linter/external edit most recently —
restore to `10.0` when resuming source+lagrange testing.

### 2. Old cutoff+source combination was simply never exercised before
Checked out the pre-merge working commit `0e21e0ed722826cbb828829e28442602a3f1f46d`
in a scratch git worktree (`ext/` symlinked in, `./configure --dim=2 &&
make realclean` needed because the bundled AMReX version differs). Confirmed:
- The old `Advance()` had the exact same unconditional hard reset
  (`if (eta<cutoff) { rho_new=rho_solid; ... }`) as pre-refactor current code.
  Byte-for-byte in `Regrid()`/`TagCellsForRefinement()` — unchanged by the merge.
- The old `tests/FlowSource1D/input` never set `cutoff` or `lagrange` at all, and
  never set `solver.type` (defaulted to **Roe**, not HLLC).
- Once `lagrange=10` and `solver.type=hllc` are forced to match the current
  branch's actual config, the **old code runs cutoff=0.1/0.2/0.5 cleanly** (to
  `stop_time=0.5`, full AMR, `amr.regrid_int=100`). So the old code's simple
  hard-reset mechanism *is* compatible with cutoff+source+lagrange in principle
  — something the EOS merge changed broke it, not the cutoff mechanism itself.

### 3. Root cause A: flux "* eta" does not suppress `inf`/`nan` (FIXED, uncommitted)
Traced via targeted `Util::ParallelMessage` debug prints (added, tested, then
removed each time) directly bracketing the `Advance()` state right before/after
`Regrid()`, and inside `HLLC.H`'s exception path.

**Mechanism:**
- At `eta ~ 0` cells, the fluid-state extraction
  `(state - (1-eta)*state_solid)/(eta+small)` used to build Riemann states is
  already large for tiny `eta` — true in *both* old and new code.
- Old code's pressure formula was a simple closed-form polynomial
  (`(gamma-1)*(E-KE) - pref`) — large but always finite.
- New code instead calls `gas.ComputeT`/`ComputeP`/`gamma`, which can genuinely
  overflow to `inf`/`nan` for the same extreme inputs.
- The resulting flux is scaled by `* eta` (~0) to suppress it at solid cells —
  but `0 * inf = nan`, not `0`. The suppression silently fails.
- That `nan`/huge flux corrupts the **fine-level** conserved state over the ~100
  steps between regrids (invisible until then). When AMR next averages fine
  data down onto the coarse level (a normal regrid operation), the corruption
  jumps onto the coarse level — this is why the crash always coincided with a
  `"Regridding on levelN"` message, and why disabling regrid
  (`amr.regrid_int=100000`) hid the crash entirely.
- Confirmed precisely: at `time=0.05` (STEP 100), `rho_old(lev=0,i=56)=1.225`
  (clean). Regrid fires. At `time=0.0505` (next `Advance()`, before any RHS
  call), `rho_old(i=56)=815.254`, `rho_old(i=57)=-2.87e11` — corruption already
  present, written by the regrid/average-down step itself.

**Fix implemented** (in `src/Integrator/Hydro.cpp`, `RHS()`, flux computation):
skip the Riemann solve entirely (return an explicit zero `Flux`) when **both**
the current cell and the relevant neighbor have `eta < small`, instead of
relying on `* eta` to suppress a possibly-`inf` result:
```cpp
Set::Scalar eta_xlo = invert ? 1.0-eta_patch(i-1,j,k)*eta_patch(i-1,j,k) : eta_patch(i-1,j,k);
// ... eta_xhi, eta_ylo, eta_yhi similarly ...
flux_xlo = (eta < small && eta_xlo < small) ? Solver::Local::Riemann::Flux(0.0,0.0,0.0,0.0) :
    riemannsolver->Solve(state_xlo_fluid, state_x_fluid, gas, molef, i, j, k, 0, small) * eta;
// ... same pattern for flux_ylo, flux_xhi, flux_yhi ...
```
This fix is real and directly verified against the specific mechanism above.
It alone was **not** sufficient to make the momentum test pass at all nonzero
`cutoff` — see root causes B and C below, found afterward.

### 4. Root cause B: `Ldot0` viscous term not tapered by `source_taper` (FIXED, uncommitted)
After root cause A's fix, the **momentum** source test (`u0=-1`) still failed,
now at `cutoff=0.0001` (the smallest tested nonzero value), with a *different*
signature: `eta` at the failing cell was `3.16e-8` — just *above* the
`eta<small` (`1e-8`) threshold from fix A, so that fix didn't apply. Traced via
debug prints in the `RHS()` flux `catch(...)` block:
- `rho(i,j)=1.225` (clean), but `Mx(i,j)=-484727` (already large), and the
  fluid-extracted velocity `u=7.3e11` (exploded via the same `eta`-division
  amplification as before, now just barely above the `small` cutoff).
- Crucially, `Source(Mx)=-1.13e7` was still huge **even though `source_taper`
  should have zeroed `mdot0`/`qdot0`/the `lagrange` penalty term** at this
  `eta` (`eta < cutoff`).
- Root cause: `Ldot0` (the viscous forcing term companion to `div_tau`,
  computed from `mu`, `hess_eta`, and `(u - u0)`) was being added to `Source`
  **completely unscaled** — no `eta`, no `source_taper` gating at all (unlike
  `div_tau`, which *is* scaled by `* eta` where it's added to the RHS). Since
  `Ldot0` depends directly on the already-exploded `u`, this created a positive
  feedback loop: bigger `u` -> bigger `Ldot0` -> bigger `Source` -> bigger `M`
  next step -> even bigger `u` (since `eta` stays tiny) -> ... runaway.

**Fix implemented** (`src/Integrator/Hydro.cpp`, `RHS()`, right after the
`Ldot0`/`div_tau` accumulation loop):
```cpp
// Ldot0 was being added to Source unscaled, unlike div_tau (scaled by "* eta"
// where it's added to the RHS below). It depends directly on the fluid-
// extracted velocity u, which is amplified without bound as eta -> 0. Without
// tapering, this term can keep injecting momentum into deep-cutoff cells even
// when source_taper has correctly zeroed mdot0/qdot0/the lagrange penalty.
Ldot0 *= source_taper;
```
**Verified effect:** re-ran the full 3-case suite across `cutoff` in
`{0.0, 0.0001, 0.001, 0.01, 0.1, 0.5}`. **Mass and energy now pass at every
tested cutoff.** Momentum still failed at all nonzero cutoff with
`relax_rate=0`.

### 5. Root cause C: residual flux leakage into cutoff-suppressed cells needs active correction (PARTIALLY FIXED, uncommitted)
Even with `source_taper` correctly zeroing all *local* source/penalty/viscous
forcing in a cutoff-suppressed cell, that cell still exchanges **ordinary
advective flux** with its active neighbor (the flux magnitude is scaled by the
neighbor's own `eta`, which is small-but-nonzero just above the interface).
Traced via debug prints again at a momentum-test crash (`cutoff=0.0001`,
`relax_rate=0`, STEP 102, right at the first regrid):
- `eta(i)=0.00107` (active, `> cutoff`): `rho(i)=1.253` clean.
- `eta(i+1)=4.7e-5`, `eta(i+2)=2.1e-6` (both `< cutoff=0.0001`, correctly
  source-suppressed): `rho(i+1)=3.1e20`, `rho(i+2)=-2.4e30` — corrupted.

So flux leakage alone, with **nothing actively correcting the state back
toward solid** (`relax_rate=0` means the RHS relaxation term introduced
earlier in this effort is a total no-op), is enough for cutoff-suppressed
cells to accumulate real momentum/energy and eventually blow up via the same
`eta`-division amplification, given enough steps.

**Fix (partial): use a nonzero `relax_rate`.** Re-tested the momentum case with
`relax_rate=100`:

| cutoff | momentum, relax_rate=0 | momentum, relax_rate=100 |
|---|---|---|
| 0.0001 | FAIL | **PASS** |
| 0.001  | FAIL | **PASS** |
| 0.01   | FAIL | **PASS** |
| 0.1    | FAIL | FAIL |
| 0.5    | FAIL | FAIL |

Tried `relax_rate` in `{100, 500, 1000, 2000}` at `cutoff=0.1` — **all fail**,
always at the exact same cell (`i=202, lev=2`). So this is a genuinely
different, third failure mode, not fixable by tuning `relax_rate`.

## Open problem: cutoff reaching into the actively-forced interface region (UNRESOLVED)

`eta.ic.expression.region0 = "0.5*tanh((4.0-x)/0.05) + 0.5"`. Solving for where
`eta=0.1`: `x = 4 + 0.05*atanh(2*0.1-1) = 4 - 0.055 = 3.945`. That's only 0.055
away from the interface center (`x=4`, `eta=0.5`, where `|grad_eta|` peaks and
`lagrange` pushes hardest toward `u0=-1`). So `cutoff=0.1` reaches **into**
the actively-forced region, not just the solid tail.

Relaxing (at any rate) a cell that the Lagrange penalty is simultaneously
forcing toward `u0=-1` creates a direct tug-of-war: the relaxation term pulls
velocity toward the solid state (`u->0`-ish) while `lagrange*source_taper`
(which is *not* fully zero yet at `eta=0.1` unless `cutoff_taper` widens the
taper band) or residual momentum diffusion keeps pushing toward `u0`. This
looks like a physical/numerical conflict inherent to overlapping the cutoff
region with the actively-forced region — not a bug fixable by better-tuning a
single scalar rate.

**Not yet tried:** making the cutoff/relaxation mechanism `grad_eta`-aware —
i.e., suppress/skip the relaxation (not just scale it by `eta`) wherever
`|grad_eta|` exceeds some threshold, regardless of the `eta` value itself, so a
nominally larger `cutoff` can be requested without conflicting with the
actively-forced interface band. This was proposed to the user but not
implemented or tested.

**Also not yet tried:** whether `cutoff_taper > 0` (smoothing the transition
in `eta`-space, independent from `relax_rate`) changes this specific
`i=202`-class failure. Earlier testing of `cutoff_taper` combined with the old
*instantaneous* hard-reset mechanism did not help (crashed even earlier, at
STEP 427 instead of 478) — but that was tested *before* the `Ldot0` taper fix
(root cause B) and *before* switching to the `relax_rate`-based relaxation
mechanism. Worth retrying `cutoff_taper` now that B is fixed.

## Current uncommitted diff (as of pause point)

- `src/Integrator/Hydro.H`: added `Set::Scalar cutoff_taper=NAN;` and
  `Set::Scalar relax_rate=NAN;` fields.
- `src/Integrator/Hydro.cpp`:
  - `Parse()`: added `cutoff_taper` (default `0.0`) and `relax_rate` (default
    `100.0`) parameters.
  - `Advance()`: removed the old instantaneous hard reset entirely (no state
    overwrite happens there anymore — only `eta`/vorticity/dt-max calc
    remains).
  - `RHS()`: added `source_taper` (ramps mass/momentum/energy source terms and
    the `lagrange` penalty from 0 at `eta=cutoff` to full strength at
    `eta=cutoff+cutoff_taper`; hard step to 0 below `cutoff` when
    `cutoff_taper=0`). Added `Ldot0 *= source_taper` (root cause B fix). Added
    a `clamp_weight = 1-source_taper` relaxation source term that subtracts
    `relax_rate*clamp_weight*(state - solid_state)` from `Source` for
    mass/momentum/energy (this is the mechanism `relax_rate` controls). Added
    the `eta<small` flux-skip logic (root cause A fix) in the flux computation
    block.
- `src/Solver/Local/Riemann/HLLC.H`: trivial, no net functional change (a
  removed/re-added comment block from debug iteration — verify with `git diff`
  before committing, should be a no-op or near-no-op vs the merge-base).
- `tests/FlowSource1D/input`: currently has `lagrange=0.0`, `cutoff=0.0`,
  `cutoff_taper=0.0`, `relax_rate=0.0` (reset to inert defaults by a
  linter/external edit at the end of the session — **restore
  `lagrange=10.0`, and pick a real `cutoff`/`relax_rate` when resuming**, per
  whatever configuration is being tested).

All temporary debug instrumentation (added and removed multiple times during
this investigation) has been cleaned up — confirmed via
`grep -rn "TEMP DEBUG\|TEMP TEST\|BLOWUP\|PREADV" src/` returning nothing.

## Verification commands used throughout

Single-case run (no `run_1d_tests.sh` wrapper, direct binary invocation, useful
for fast iteration with parameter overrides):
```
./bin/hydro-2d-g++ tests/FlowSource1D/input \
  m0.ic.expression.region0="10.0" u0.ic.expression.region0="-1.0" q.ic.expression.region0="0.0" \
  cutoff=<X> cutoff_taper=<Y> relax_rate=<Z> stop_time=0.5 \
  plot_file=<scratch path>
```
Full 3-case suite (mass/momentum/energy), matching `run_1d_tests.sh`'s
overrides:
```
# mass:      u0=0.0   q=-0.001
# momentum:  u0=-1.0  q=0.0      <-- most sensitive case, use this for iteration
# energy:    u0=0.0   q=-10.0
```
`make -j8 ./bin/hydro` to rebuild; a full clean rebuild
(`rm obj/obj-2d-g++/Integrator/Hydro.cpp.o ...`) is worth doing before trusting
a "success" result if anything seems inconsistent with prior runs — one round
of this investigation produced a false "all pass" result that turned out to be
from a stale/incorrectly-verified binary state; always re-verify with a fresh
build before reporting results.

## Where to pick this back up

1. Restore `tests/FlowSource1D/input`'s `lagrange=10.0` (currently reset to
   `0.0`).
2. Decide on the `grad_eta`-aware cutoff idea (make relaxation strength depend
   on `|grad_eta|`, not just `eta`) vs. simply documenting that `cutoff` must
   stay below ~`0.01`–`0.05` for this interface width (`eps=0.05`) when
   `lagrange`/source terms are active.
3. Re-verify the full 3-case suite across a `cutoff` sweep with a **fresh,
   clean rebuild** before trusting any "pass" result (see note above about the
   earlier false-positive).
4. Once a working configuration is settled, commit incrementally following the
   plan's git workflow (see `cryptic-humming-boole.md`): each logical fix
   (Ldot0 taper, flux-skip, relax_rate mechanism) as its own commit, trailer
   `Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>`, no push/PR unless
   asked.

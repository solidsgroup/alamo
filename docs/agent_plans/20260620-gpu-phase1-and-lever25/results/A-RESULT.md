# Result: Task A — Lever 2.5 field packing

## Summary
**blocked** (architectural — no source change was made; analysis disqualified the
pack before any code was written, per "do not force a fragile change to
manufacture a win"). This is the documented negative-result outcome the task
explicitly allows.

## Files changed
None. `git diff --stat src/` is empty. `src/Integrator/Flame.cpp`,
`src/Integrator/Flame.H` are untouched.

## Measurements (baseline only — no "after" state exists)
Captured before any implementation attempt, `wide_512_bf32_mgs128` box-sweep
config (`amr.n_cell="64 64 64" amr.max_level=3 amr.blocking_factor=32
amr.max_grid_size=128`, gpu_fast binary `bin/alamo_gpu-2d-cuda86-g++`, np=1,
max_step=5):

- Total GPU kernel launches (5 steps): **63,505** → 12,701/step
- `BC::Constant::FillBoundary` device-kernel instances: **34,108 (53.7%** of all
  GPU kernel launches) — matches the task's documented 56% within capture-config
  noise.
- FillBoundary GPU **execution time**: 35.1 ms of 233 ms total kernel time
  (**15.1%**) — these are cheap kernels; the 54% launch-count share does not
  translate 1:1 into wall-time share. The lever's available win is bounded by
  host-side launch/sync overhead, not device occupancy.
- `cudaLaunchKernel`: 63,505 calls, avg 6.87 µs, total 436 ms (44.3% of API time)
- `cudaStreamSynchronize`: 24,237 calls, avg 9.9 µs, total 240 ms (24.4%)
- Box-sweep gate baseline (`wide_512_bf32_mgs128`, recorded before any change):
  CPU np8 = 1.418 s, GPU fast = 1.835 s, GPU/CPU = 1.294 — `benchmark/phase2_box_sweep/summary.csv`
  shows **5/5 ok** (all cases bit/tolerance-correct).
- Gate state confirmed **9/9 ok** via `CPU_NP=8 GPU_FAST_NP=1 GPU_STRICT_NP=1
  python3 benchmark/baseline_suite.py check` (strict bit-identical) before any
  attempt — this was the only golden-compare run performed, since no code was
  changed; nothing to re-gate.
- No "after" launches/step, FillBoundary share, or wall/step number exists —
  implementation was not attempted (see below).

## What was packed and why (and what was left unpacked)
**Nothing was packed.** Before writing any code, I traced every consumer of
`eta_mf` and `psi_mf` (the task's designated "low-risk, same-BC" pack) and found
three independent, mutually-reinforcing blockers that make this pack neither
low-risk nor mechanically simple:

1. **`Newton.H`/`Elastic.cpp` read psi at an implicit component 0.**
   `Flame::Parse` calls `value.solver.setPsi(value.psi_mf)`
   (`Flame.cpp:295`), storing `Set::Field<Set::Scalar>* m_psi` in `Newton.H:60`.
   Every downstream read — `Newton.H:210` (`m_elastic->SetPsi(i, *(*m_psi)[i])`),
   `Newton.H:256` (`(*m_psi)[lev]->array(mfi)`), and `Elastic.cpp:145-146`
   (`Elastic<SYM>::SetPsi`'s `a_psi.array(mfi)`) — index `(i,j,k)` with **no
   explicit component**, i.e. component 0 only. Packing eta(comp0)+psi(comp1)
   would make `setPsi` silently read **eta's** values as psi — a severe,
   silent correctness bug, not a build error. The only viable component
   ordering is **psi=comp0, eta=comp1** — but that breaks point 2 below. Fixing
   this generally means editing `Newton.H`/`Elastic.cpp`/`Elastic.H`, which are
   solver-core framework files, explicitly outside this task's allowed-files
   list (only Flame.cpp/Flame.H and Flame-local BC construction are permitted).

2. **IC initializers write to an implicit component 0.** `Flame::Initialize`
   calls `ic_eta->Initialize(lev, eta_mf)` (`Flame.cpp:315`) and
   `Flame::Regrid` calls `ic_eta->Initialize(lev, eta_0_mf, time)`
   (`Flame.cpp:915`); `IC::BMP::Add` (and Laminate/Constant/Expression/PNG/
   PSRead — every `IC::*::Add`) writes via `a_field[lev]->array(mfi)` indexed
   `field(i,j,k)`, again implicit component 0. If eta is comp0 (required to
   keep IC init correct without touching ≥5 IC files outside Flame), psi must
   be comp1 — directly contradicting point 1's requirement that psi be comp0
   for `setPsi` to read it correctly. **There is no single component ordering
   that satisfies both constraints** without editing files outside this task's
   scope (`IC/*.H`, `Solver/Nonlocal/Newton.H`, `Operator/Elastic.cpp/.H`).

3. **`bc_eta` (ncomp=1) is shared by `eta_old_mf` and `eta_0_mf`, not just
   eta/psi.** `Flame::Parse` constructs one `BC::Constant bc_eta` with
   `ncomp=1` (`Flame.cpp:118`) and reuses the same pointer for `eta_mf`,
   `eta_old_mf`, `psi_mf`, and `eta_0_mf` (`Flame.cpp:122-128`). `eta_0_mf` is
   registered with the default `evolving=true`, so it goes through the
   framework's generic per-step `FillPatch` (`Integrator.cpp:1218`) using
   `bc_eta` just like eta/psi. `BC::Constant::FillBoundary` has a hard runtime
   assertion (`Constant.cpp:68`: `Util::Assert(a_in.nComp() == m_ncomp)`).
   Packing eta+psi into one 2-component fab requires `bc_eta.m_ncomp == 2` to
   pass that fab's FillBoundary — but then `eta_old_mf`/`eta_0_mf` (genuinely
   1-component) would fail the same assertion on every regrid/timestep. This
   forces constructing a **second**, dedicated `BC::Constant` instance (ncomp=2,
   cloned from the same `pf.eta.bc` input keys) solely for the packed fab —
   workable in principle (within Flame's own setup, so technically in-scope),
   but it is the one piece of the three that is *not* itself disqualifying.

Point 3 alone would be a manageable, contained change. Points 1 and 2 are not:
each requires edits outside the explicitly allowed file list, and there is no
component ordering that resolves both simultaneously without those edits. A
genuinely "low-risk, same-BC" eta+psi pack does not exist in this codebase as
it stands — the framework's implicit-component-0 convention (used by `setPsi`,
every `IC::*::Add`, and several other generic consumers found while tracing
`cell.fab_array`/`physbc_array`) is incompatible with silently widening any
existing 1-component evolving fab to 2 components.

The cross-BC temp pack (step 4 of the task) was not attempted: it strictly
required step 2 (eta+psi) to be green first, and step 2 never reached an
implementation that passed even a basic design review.

## Tests run and results
- `CPU_NP=8 GPU_FAST_NP=1 GPU_STRICT_NP=1 python3 benchmark/baseline_suite.py check`
  → **9/9 ok**, strict bit-identical (run once, pre-implementation; confirms the
  tree was green before this task started and remains green since no source was
  touched).
- `python3 benchmark/phase2_box_sweep.py` → **5/5 ok** (`benchmark/phase2_box_sweep/summary.csv`),
  `wide_512_bf32_mgs128`: cpu=1.418s, gpu=1.835s, gpu/cpu=1.294.
- nsys capture (`wide_512_bf32_mgs128` config, gpu_fast binary, max_step=5):
  succeeded on retry after one transient OOM (free GPU mem dropped to 500 MiB
  momentarily, almost certainly Track B's concurrent elastic job; retried once
  GPU memory recovered to 6.6 GB free, per the plan's OOM-retry guidance).
  Reports saved at `/tmp/lever25_nsys/baseline.nsys-rep` (+ `.sqlite`,
  `_stats_cuda_api_sum.csv`, `_stats_cuda_gpu_kern_sum.csv`); copies of the two
  CSVs at `/tmp/A_baseline_cuda_api_sum.csv` and `/tmp/A_baseline_cuda_gpu_kern_sum.csv`.

## Issues found
- See "What was packed" above — the three architectural blockers are the
  finding. No build, link, or runtime errors were hit because no code was
  written; the blockers were identified by static tracing of every consumer of
  `eta_mf`/`psi_mf` (`grep -rn "psi_mf|setPsi" src/`) before committing to an
  implementation, specifically to avoid discovering this mid-refactor across
  three rebuilds.
- One transient GPU OOM during nsys capture, consistent with the plan's
  documented single-shared-GPU contention risk; resolved by retry once free
  memory recovered. No action needed beyond what the plan already prescribes.

## Deviations from task
- Did not implement the same-BC eta+psi pack described in the task as
  "low-risk." Investigation (not present in the task's own framing) found it
  is not low-risk: it requires either editing `Newton.H`/`Elastic.cpp` (psi's
  implicit-comp0 consumer) or every `IC::*::Add` (eta's implicit-comp0
  initializer), both outside the allowed-files list, with no component
  ordering that avoids both. This is reported as new information for whoever
  scopes a future attempt, rather than silently working around it.
- Did not attempt the cross-BC temp pack (step 4), since it depended on step 2
  succeeding first.
- Did not produce an "after" nsys capture or wall/step delta, since no
  implementation exists to measure.

## Recommendation
**Do not pursue Lever 2.5 for eta+psi (or eta+psi+temp) under the current
framework conventions.** Reasons, combining the architectural finding with the
task's own quantitative bar:
1. Architectural: no compliant (in-scope-files-only) component ordering exists
   for eta+psi — see blockers 1 and 2 above.
2. Even if blockers 1/2 were waived (i.e., `Newton.H`/`Elastic.cpp`/IC files
   were back in scope for a future, larger task), the back-of-envelope ceiling
   is unfavorable: FillBoundary is 53.7% of launch *count* but only 15.1% of
   GPU kernel *time* (these are cheap kernels); removing one of ~5 evolving
   scalar BC fabs (eta+psi merge) saves roughly 1/5 of the 34,108 FillBoundary
   launches (~6,800 launches, ~10.7% of all GPU launches), i.e. roughly
   6,800 × 6.9 µs ≈ 47 ms of host-side `cudaLaunchKernel` submission cost
   across 5 steps (~9.4 ms/step) against a 367 ms/step GPU box-sweep baseline
   (`wide_512_bf32_mgs128`) — about **2.5%**, comfortably under the task's 5%
   stop threshold, and that ignores that a fair fraction of launch-submission
   time overlaps async GPU execution rather than sitting on the wall-clock
   critical path.
3. Both points independently support stopping; together they make this a
   clean, confident negative result rather than a marginal call.

No revert is needed since no commit-worthy change exists — the working tree's
`src/` is byte-identical to the task's starting point. If a future task wants
to revisit this, the prerequisite is widening scope to include `Newton.H`,
`Elastic.cpp`/`.H`, and the `IC::*::Add` family so a consistent component
convention can be chosen — at which point it is a different, larger task than
"Lever 2.5" as scoped here.

# GPU Acceleration Roadmap v3 — Validate, Then Optimize (beta → release candidate)

Branch: `chamber-gpu` (never merged — see `benchmark/GPU_BRANCH_GUIDE.md`).
Supersedes: GPU Roadmap v2 (now at `benchmark/archive/GPU_ROADMAP_V2.md`).
v2's Phase A is **measured** and its B1 gate **passed**; v3 absorbs v2's open items
(see §2) and adds the priorities this arc is actually about.

Created: 2026-06-30.
Status: **Phase 1 ACTIVE (physics validation suite).** Nothing downstream ships
until Phase 1 exists.

> ## The v3 inversion
>
> v2's rule was: *"No kernel-level optimization without a counter that justifies
> it."* That rule got us the counters. It told us **where** the time goes
> (`Operator::Elastic::Fapply` = 74.6% of GPU kernel time, elastic solve ≈ 95% of
> wall). It did **not** tell us whether the optimizations we are about to make
> keep the answer correct.
>
> The remaining work is no longer "make CUDA compile" or even "make it fast." The
> bugs left are **physics-validity, async-lifetime, and performance-architecture**
> bugs — the dangerous class, because the code still runs and still produces
> numbers. A loosened tolerance, a skipped solve, a register-shaved kernel can all
> look faster and stable while silently shifting the pressure history, the burn
> area, or the stress localization that the chamber model exists to predict.
>
> So v3's governing rule is:
>
> > **No optimization ships without a physics-error budget that proves the coupled
> > burn / pressure / stress behavior is preserved within a stated tolerance.**
>
> Counters prove a change is *faster*. The error budget proves it is still *the
> same simulation*. Phase 1 builds that budget as a one-command suite. Everything
> else is gated behind it.

---

## 1. Why a v3 roadmap

Three forces converged at the end of v2:

1. **The strategic target collapsed onto one kernel.** Phase A measured the
   combined flame+elastic A100 workload: elastic MLMG is ≈95% of evolve time,
   Fapply is 74.6% of GPU kernel time at ~255 registers/thread (≈12.5%
   occupancy), and flame is 0.20%. *Elastic is the project now.* Flame is frozen
   unless a correctness bug appears. (`benchmark/archive/PHASE_A_FINDINGS.md`.)

2. **The cheap wins became dangerous.** The Phase C0 input levers (relax
   `tol_rel/abs` 1e-8 → 1e-6, `interval` 50 → 100, tune `bottom_*`) showed a
   ~2.5× wall-progress gain on NOVA — but with **no stress-field parity check**
   in the bundle. A 2.5× "win" that quietly changes the stress field is not a
   win. We cannot keep stacking tolerance, bottom-solver, and interval changes
   and calling the result valid. (`benchmark/archive/PHASE_C0_input_ab.md`.)

3. **The headline number is not yet defensible.** The 22.9× combined speedup is
   measured against a 64-rank CPU baseline that leaves **32 ranks idle** during
   elastic (`Operator.cpp:459` load imbalance). The GPU win is real, but the
   ratio is inflated; it should not be published until the CPU baseline is fair.

v2 was *measurement-driven*. v3 is *validation-driven*: build the instrument that
defines "correct enough," fix the measurement so the numbers are honest, **then**
spend effort on the kernel rewrite that the data has earned.

---

## 2. What v2 settled, and what v3 inherits

**Settled by v2 — cite, do not re-litigate** (full record in
`benchmark/SUCCESS_BOOK.md` + `benchmark/archive/`):

- Phase-field wins ~10× on GPU in the saturating 3D regime (9.6× @128³, 13.1×
  @256³, single A100 vs 64-rank CPU). `archive/PHASE3_R3_crossover.md`.
- Elastic runs correctly on GPU multi-box (2D **and** 3D). Two device-only fault
  classes were root-caused: the cross-stream `tmpfab` use-after-free
  (`Operator.cpp:728`, `elixir()`, commit `c00f69086`) and the chained-Eigen
  `F.inverse().transpose()` fault (`NeoHookean.H`).
- Combined flame+elastic A100 baseline measured: 22.9× per step at 256³, elastic
  fraction `f ≈ 0.95`. `archive/PHASE_A_FINDINGS.md`. **(Caveats: CPU baseline
  unfair — see 2.A; runs were SLURM-preempted — see 2.B.)**
- Branch DoD 5/5; CPU regression suite green. `archive/PHASE5_BRANCH_DONE.md`.

**Open items v2 left on the table — v3 absorbs each into a numbered task:**

| v2 item | v2 status | Lands in v3 as |
| --- | --- | --- |
| A1 — Fapply Speed-of-Light re-export (ncu) | partial (report in newer format, not yet read) | **2.C** |
| A3 — lock metric set + GPU-CI perf regression | not done | **2.D** |
| A4 — sibling GPU-only-fault audit | not done | **4.B / 4.C** |
| B3 — fair CPU baselines | not done | **2.A** |
| C0 — input levers (tol/interval/bottom) | promising, unvalidated | **2.B** (now gated on physics budget) |
| C1 — Fapply register/occupancy edits | staged on `chamber-gpu-elastic-opt`, no A100 A/B | **3.A** |
| C2 — Fapply AoS→SoA + uniform fast-path | not started | **3.C / 3.D** |
| C4 — operator reuse across solves | not started | **3.E** |
| D1/D2 — multi-GPU / H100 envelope | parked | **5.A / 5.B** |

Nothing is lost; it is re-sequenced behind the validation gate.

---

## 3. Guiding principles

Carried forward from v2, plus the v3 additions:

1. **Validate before you trust a speedup.** Every optimization — input lever or
   kernel rewrite — passes the Phase 1 physics-error budget before it counts as a
   win. *(new, governing)*
2. **Measure before you optimize.** A counter still justifies every kernel change
   (v2's rule, retained).
3. **Two builds, two purposes.** Perf headlines use the fast-math binary;
   correctness/validation claims use **only** the strict / no-fast-math binary.
   Every table labels which: `cpu`, `gpu_strict`, `gpu_fast`. Never mix.
4. **Combined physics is the unit of measurement.** Flame-only and elastic-only
   numbers are diagnostic sub-measurements, not the headline.
5. **Engineering quantities over time, not just final fields.** A skipped or
   loosened solve can match the final field and still phase-lag the pressure
   trajectory. The budget compares *trajectories*.
6. **Single-GPU is the supported shape** at 128³/256³ until 3.x makes the Fapply
   path healthy. Multi-GPU is a parked, measured negative (halo-bound).
7. **Reuse the harness, don't fork it.** Phase 1 extends
   `golden_compare_flame.sh` / `compare_thermo.py` / `baseline_suite.py`, not a
   parallel rewrite. Avoid building "two ALAMOs."

---

## 4. Phases & tasks

Five phases, gated in sequence. Each task: **Goal · Steps · Done-when · Depends ·
Artifact.** IDs are stable references (`1.A`, `1.B`, …, `5.B`).

The gate chain: **Phase 1 (the instrument) → Phase 2 (honest numbers) → Phase 3
(the kernel win) → Phase 4 (hardening) → Phase 5 (envelope).** Phase 4 audits may
run opportunistically in parallel; everything else is strictly ordered.

---

### Phase 1 — Physics validation suite & error budget  ⟵ THE PRIORITY

Build the instrument that defines "correct enough" and runs from **one command**
on either the local A1000 (sm_86) or a NOVA A100 (sm_80), emitting organized,
machine-diffable output bundles for storage and comparison. Until this exists, no
Phase 2/3 result is trustworthy.

Design home: `benchmark/validate/`.

---

**1.A — Physics-quantity registry & error budget**
- **Goal:** a single written contract listing every observable the suite tracks,
  how it is extracted, and the tolerance class it must satisfy. This *is* the
  definition of "correct enough."
- **Steps:** enumerate the validation observables and bin each into a tolerance
  class:
  - **CORRECTNESS** (must match CPU-strict to ~1e-6 relative, no-fast-math):
    final per-component field L2/L∞ norms (η/phase field, temperature, displacement,
    stress components), residual at convergence.
  - **ENGINEERING-TRAJECTORY** (must match within a stated physical % over the
    *whole* run, with bounded phase-lag): chamber pressure history, chamber
    volume, burn area, total mdot, burn-front / interface displacement, max & mean
    von Mises stress, max principal stress, elastic energy.
  - **SOLVER-HEALTH** (regression-tracked, not pass/fail on absolute value):
    Newton iteration count per solve, MLMG V-cycle count, bottom-solver
    iterations, residual-history monotonicity.
  For each observable record: name, definition, extraction source
  (`thermo.dat` column / TinyProfiler region / plotfile field-norm), tolerance
  class, and the numeric threshold.
- **Done-when:** `benchmark/validate/physics_budget.md` (+ a machine-readable
  `physics_budget.yaml` the comparator reads) lists every observable with source
  and threshold; reviewed against the chamber model's engineering intent (does a
  passing budget actually mean the burn/pressure/stress prediction is preserved?).
- **Depends:** none. **Artifact:** `benchmark/validate/physics_budget.{md,yaml}`.

**1.B — Single-entry validation runner**
- **Goal:** one script that runs the whole suite on whatever hardware it finds,
  with zero hand-editing between the A1000 and the A100.
- **Steps:** build `benchmark/validate/run_validation.sh` that:
  1. **auto-detects hardware** (parse `nvidia-smi` / `SLURM_*`; map sm_86→local
     A1000, sm_80→A100, sm_90→H200) and selects the matching binary
     (`bin/alamo_gpu-*`), defaulting to the **strict / no-fast-math** build for
     validation; `--fast` opt-in for a perf-labelled smoke row only.
  2. selects a **device-appropriate case set**: 2D + small-3D cases that fit the
     A1000's 8 GB locally; the 256³ chamber case on A100. Cases live in a manifest
     so the same script scales up/down by hardware.
  3. for each case, runs the CPU-strict reference **and** the GPU candidate (or,
     with `--compare-to <ref-bundle>`, the GPU candidate vs a stored reference),
     under managed arena where the A1000 needs it.
  4. is idempotent and re-runnable; never clobbers a prior bundle.
  - Build on the existing harness: this wraps/extends `golden_compare_flame.sh`,
    `compare_thermo.py`, and the `baseline_suite.py` record/check/report flow —
    do not reinvent the runner. The new parts are hardware auto-detect, the
    elastic/stress observables (1.A) the current thermo-only compare lacks, and
    the bundle schema (1.C).
- **Done-when:** `bash benchmark/validate/run_validation.sh` completes on the
  local A1000 and (via a committed `.slurm` wrapper the user submits) on a NOVA
  A100, producing one output bundle per case with no per-host edits.
- **Depends:** 1.A, 1.C (schema). **Artifact:** `benchmark/validate/run_validation.sh`
  + `benchmark/validate/cases.manifest`.

**1.C — Organized output schema (storage + comparison)**
- **Goal:** the user's core requirement — outputs so well-defined and organized
  that two runs, weeks apart on different hardware, drop into a comparator with no
  fixups.
- **Steps:** define a deterministic bundle layout:
  ```
  benchmark/validate/runs/<UTCstamp>_<gitSHA>_<host>_<device>/
    manifest.json     # host, GPU, driver, git SHA, build flags, binary, input hashes, np
    <case>/
      metrics.json    # every observable from physics_budget, as {value, class}
      thermo.dat      # raw scalar history (copied)
      region_times.txt# TinyProfiler regions (elastic vs phase-field split)
      field_norms.json# per-component L2/L∞ of final plotfile
      run.log         # full stdout
  ```
  `manifest.json` makes every bundle self-describing and reproducible; the
  per-case `metrics.json` is the canonical comparison surface. Use UTC timestamps
  + short git SHA + hostname + device tag so bundles sort chronologically and
  never collide across machines.
- **Done-when:** the schema is documented in `benchmark/validate/README.md`; a
  produced bundle validates against it; `metrics.json` carries every 1.A
  observable with its tolerance class tagged.
- **Depends:** 1.A. **Artifact:** `benchmark/validate/README.md` (schema spec).

**1.D — Comparator & verdict report**
- **Goal:** turn two bundles into a pass/fail physics-budget verdict a human and a
  CI job can both read.
- **Steps:** build `benchmark/validate/compare_validation.py`:
  - `compare_validation.py <bundleA> <bundleB>` → a per-observable table {name,
    class, A, B, Δ, tolerance, PASS/FAIL} + a roll-up verdict (any CORRECTNESS
    fail ⇒ overall FAIL; ENGINEERING fails surfaced with magnitude; SOLVER-HEALTH
    deltas reported, non-gating).
  - **Trajectory comparison, not endpoint-only**: for ENGINEERING observables,
    compare the full time series — max relative deviation over the burn **and** a
    phase-lag / cross-correlation metric — so a stale-but-converged field can't
    pass on its final value alone.
  - Emit both human Markdown (`compare_report.md`) and machine JSON
    (`compare_report.json`).
- **Done-when:** running the comparator on a known-good vs a deliberately-detuned
  bundle correctly returns PASS and FAIL respectively, with the offending
  observable named.
- **Depends:** 1.A, 1.C. **Artifact:** `benchmark/validate/compare_validation.py`.

**1.E — Golden references (CPU-strict + alpha-1.0 GPU)**
- **Goal:** populate the budget with real reference values so future changes have
  something to be checked against.
- **Steps:** run the suite to record (a) the **CPU-strict golden bundle** — the
  trusted reference for every CORRECTNESS observable — and (b) the current
  **alpha-1.0 GPU bundle** (tip of `chamber-gpu`, strict build). Store both under
  `benchmark/validate/references/` with their manifests. Commit the small
  `metrics.json`/`manifest.json`; keep bulky plot/thermo artifacts out of git
  (gitignore + a documented external location).
- **Done-when:** `references/cpu_strict/` and `references/gpu_alpha1/` exist and
  the comparator reports the GPU-vs-CPU deltas, establishing the *current*
  physics-error baseline (what passes today, before any v3 change).
- **Depends:** 1.B, 1.D. **Artifact:** `benchmark/validate/references/` +
  `benchmark/validate/BASELINE_DELTAS.md`.

**1.F — Make validation the gate**
- **Goal:** wire the suite in so no Phase 2/3 change can land without passing it.
- **Steps:** add `run_validation.sh --gate <candidate-bundle>` returning non-zero
  on any CORRECTNESS regression (or any ENGINEERING regression beyond budget);
  document it in `benchmark/READ_FIRST_NEXT_STEP.md` as a required pre-commit step
  for hot-path edits; add it as a *gated, non-blocking* job in
  `.github/workflows/chamber-gpu-correctness.yml` that runs on a `[self-hosted,
  cuda]` runner and prints `skipped: no GPU` otherwise.
- **Done-when:** the gate command exists, is documented as mandatory for Phase
  2/3, and the CI job is present and skip-safe on CPU-only runners.
- **Depends:** 1.D, 1.E. **Artifact:** gate mode + workflow job + READ_FIRST update.

---

### Phase 2 — Honest measurement

Make the numbers defensible: fix the unfair CPU baseline, run the Phase C input
levers as *isolated, budget-gated* A/B rows, finish the device counters, and lock
the metric set into CI. Every row here passes the Phase 1 gate.

---

**2.A — Fair CPU baseline**
- **Goal:** retire the inflated 22.9× by occupying all 64 CPU ranks during the
  elastic solve. **Steps:** the operator exposes only 32 boxes, idling half the
  node (`Operator.cpp:459`); re-block so every rank carries elastic work
  (`max_grid_size` / box count), same physics, same output cadence, same
  convergence criteria. Re-run the 256³ combined baseline. **Done-when:** a table
  reporting both `project baseline (22.9×)` and `fair baseline (corrected)`,
  node-hours / GPU-hours alongside step time. **Depends:** Phase 1.
  **Artifact:** `benchmark/PHASE_2A_fair_cpu_baseline.md`.

**2.B — Clean isolated A/B for the Phase C input levers**
- **Goal:** attribute the C0 "2.5×" to specific levers, each gated on physics
  parity — replacing the stacked, unvalidated C0 result. **Steps:** run as
  *separate* rows, never pre-stacked: `{baseline · tol-only · bottom-only ·
  interval-only · tol+interval · all-tuned}`, via `benchmark/phase_c_elastic_ab.sh`
  at 256³/A100. For each row capture: per-solve wall, V-cycle count,
  bottom-solver iters, Fapply time, full step time — **and** the Phase 1 verdict
  (run the bundle through `compare_validation.py` vs the CPU-strict golden).
  **Done-when:** the reviewer's A/B table {case → wall, V-cycles, stress/physics
  delta, PASS/FAIL}; a documented set of levers that are *both* faster and
  budget-clean (those ship; the rest are dropped). **Depends:** Phase 1, 2.A.
  **Artifact:** `benchmark/PHASE_2B_input_levers_ab.md`.

**2.C — Complete the device counter set (Fapply Speed-of-Light)**
- **Goal:** the per-kernel memory-vs-compute split that gates the Phase 3 kernel
  rewrite. **Steps:** finish the NOVA-side `ncu` Speed-of-Light re-export for
  `Fapply`/`Diagonal`/`Newton::prepareForSolve` (the v2 A1 report is in a newer
  ncu format the local 2022.4.1 can't open). Capture: SoL memory %, SoL compute %,
  achieved occupancy, registers/thread, local-memory loads/stores (spill traffic),
  DRAM throughput, dominant warp-stall reason. Use the `--nvtx-include
  "Operator::Elastic::Fapply()/"` + explicit `--metrics`/`--section-folder` form
  (the empty-ncu root cause from the v2 era — see
  `benchmark/archive/`). **Done-when:** a counter table that states plainly
  whether Fapply is memory-bound or compute/occupancy-bound — the fact that
  decides 3.C (SoA layout) vs 3.D (register/math reduction) priority.
  **Depends:** Phase 1. **Artifact:** `benchmark/PHASE_2C_fapply_sol.md`.

**2.D — Lock the standing metric set into CI**
- **Goal:** make every run emit the metric set automatically and catch
  regressions. **Steps:** extend `benchmark/perf_regression_track.py` to also
  record the elastic/phase-field split **and** the Phase 1 budget verdict; add the
  gated GPU perf-regression job to `chamber-gpu-correctness.yml` (`[self-hosted,
  cuda]`, skip-safe). **Done-when:** a perf-regression row is produced by CI on a
  GPU runner (or verified skip-safe on CPU-only), CSV schema includes the elastic
  split and the budget pass/fail. **Depends:** 1.F, 2.A. **Artifact:** updated
  `perf_regression_track.py` + workflow job + first CI row.

---

### Phase 3 — The structural kernel win (elastic `Fapply`)

The real performance prize: raise Fapply's ~12.5% occupancy. Order of attack is
**set by 2.C** — if memory-bound, lead with 3.C; if occupancy/register-bound,
lead with 3.D. **No kernel change lands without an A100 before/after AND a Phase 1
gate pass** (strict-build golden compare). Keep the generic operator as the
CPU/reference path; specialize only the GPU hot path.

---

**3.A — Land the staged register edits with an A100 A/B**
- **Goal:** measure and ship (or reject) the two bit-identical edits already
  staged on worktree branch `chamber-gpu-elastic-opt`: grad(C) accumulated one
  direction at a time (one live `Matrix4` derivative temp instead of three) + the
  boundary-only `sig`-sink that drops a `Matrix4·Matrix` product from every
  interior node. **Steps:** rebuild 3D strict + fast on A100; capture occupancy +
  registers/thread + wall/step before/after (2.C harness); run the Phase 1 gate.
  **Done-when:** occupancy up + wall down + budget PASS ⇒ ship; else a documented
  "no win" with the counter. **Depends:** 2.C. **Artifact:**
  `benchmark/PHASE_3A_fapply_register_edits.md`.
  - **Supporting evidence (2026-07-02 code audit,
    `benchmark/GPU_AUDIT_20260702.md` §5):** direct kernel reading confirms
    grad(C) (`Elastic.cpp:615-624`) is the dominant register-pressure driver —
    three live `Set::Matrix4` temps (63 doubles) per cell, and the
    `if (!m_uniform)` guard around it never actually skips the branch for
    chamber sims (`Mechanics.H:198` forces `SetUniform(false)` unconditionally)
    — i.e. the C1 edit targets a hotspot that is *always* live, not a rare
    path. The audit also found two adjacent, not-yet-staged opportunities:
    `DDW(i,j,k)` is loaded from global memory twice per cell
    (`Elastic.cpp:613` and `:629`) when psi is active — cache it once; and the
    in-kernel boundary branch (`Elastic.cpp:537-670`) is the fuller motivation
    for task 3.B below, not just the C1 `sig`-sink partial fix.

**3.B — Interior/boundary kernel separation**
- **Goal:** stop every interior node from carrying boundary logic and full DDW
  state. **Steps:** split Fapply so boundary-node handling is its own launch;
  interior nodes drop the boundary branch and its live variables. **Done-when:**
  measured register/occupancy/wall improvement on the interior kernel, budget
  PASS. **Depends:** 2.C, 3.A. **Artifact:** `benchmark/PHASE_3B_interior_boundary_split.md`.

**3.C — Uniform/nonuniform material split + interface mask**
- **Goal:** the >90% of nodes that are interior-uniform (3 piecewise-constant
  materials) should skip the 6 neighbor `Matrix4` loads and the grad(C) branch
  entirely. **Steps:** precompute a thin material-interface mask once per solve;
  uniform-interior nodes take a fast path; if register needs diverge sharply,
  compile/launch separate uniform vs nonuniform kernels. If 2.C says memory-bound,
  this is the lead task. **Done-when:** measured wall reduction on the elastic
  region, budget PASS. **Depends:** 2.C, 3.A. **Artifact:**
  `benchmark/PHASE_3C_uniform_fastpath.md`.

**3.D — Replace full `Matrix4` math in the hot path (the big one)**
- **Goal:** get Fapply below ~96–128 registers/thread by not materializing general
  tensor objects in a stencil kernel. **Steps:** exploit the symmetry of the
  material tangent (`DDW`); write the needed contractions as explicit scalar
  arithmetic instead of constructing full `Set::Matrix4` temporaries; precompute
  compact per-node coefficients where the tangent is reused across MLMG
  iterations. Accept an inelegant kernel — keep the readable version as the
  CPU/reference path. **Done-when:** registers/thread under target, occupancy and
  wall improved on A100, budget PASS (strict golden compare clean). **Depends:**
  2.C, 3.A (and informs 3.C). **Artifact:** `benchmark/PHASE_3D_matrix4_scalarize.md`.

**3.E — Operator reuse across solves**
- **Goal:** stop rebuilding the operator every solve. **Steps:**
  `Mechanics::TimeStepBegin` calls `Elastic::define()` (~2.4 GB alloc) +
  `prepareForSolve` per solve (`Mechanics.H:200-209`); only coefficients change.
  Hoist construction; re-`SetModel` only. **Done-when:** per-solve setup time
  drops, budget PASS. **Depends:** 3.A. **Artifact:**
  `benchmark/PHASE_3E_operator_reuse.md`.

**3.F — Re-tune launch/block config at the new register footprint**
- **Goal:** the best block size at 255 regs/thread is not the best at <128.
  **Steps:** after 3.D lands, sweep block size / `__launch_bounds__` against
  achieved occupancy and wall. **Done-when:** a chosen launch config with a
  counter-backed justification. **Depends:** 3.D. **Artifact:** note appended to
  `PHASE_3D_matrix4_scalarize.md`.

---

### Phase 4 — Correctness hardening (sibling-bug audits + sanitizers)

Two device-only fault classes already bit elastic (elixir UAF, chained-Eigen
fault). Both usually have siblings. This phase finds them and locks in
sanitizer gates so the Phase 3 rewrites can't reintroduce them. 4.A–4.D may run in
parallel with Phase 2/3 (they are independent correctness insurance).

---

**4.A — Sanitizer gate matrix**
- **Goal:** `compute-sanitizer memcheck` is necessary but not sufficient. **Steps:**
  run `memcheck` + `racecheck` + `initcheck` across the matrix {single-box ·
  multi-box · AMR · restart} × {thermal-on/off · elastic-on/off} × {strict ·
  fast-smoke}. Reuse `benchmark/local_a100_gate.sh` (HMM-immune local gate) as the
  driver. **Done-when:** a results grid; every cell green or a filed defect.
  **Artifact:** `benchmark/PHASE_4A_sanitizer_matrix.md`.

**4.B — Async-lifetime audit beyond `FArrayBox`**
- **Goal:** the elixir audit found the obvious temp-fab case; widen the pattern.
  **Steps:** audit every device-lambda capture of a temporary or stack-owned
  object — transient `Array4`s from short-lived fabs, coefficient arrays,
  `std::vector/array`, Eigen temporaries, and **every `[this]` capture** in a
  device path (the traction-diagnostic bug proved `this` is dangerous). For each:
  elixir-needed / lifetime-guaranteed / fix. **Done-when:** an audit table, each
  hit with a verdict, fixes verified by re-running the multi-box elastic + a 3D
  smoke under 4.A. **Artifact:** extend `benchmark/archive/elixir_race_audit.md`
  → `benchmark/PHASE_4B_async_lifetime_audit.md`.

**4.C — Host-access-from-device audit**
- **Goal:** the `inLaunchRegion()` IC/BC guards protect against kernel launches,
  not against host writes into device-arena memory. **Steps:** for every guarded
  path ask: does it write through an `Array4` on the host? which arena owns the
  data? does the guard still trigger for small test cases? does it differ under
  managed / device / HMM memory? (A small local grid can pass and still be invalid
  on NOVA — see the HMM-masking note in `benchmark/LOCAL_A100_SPOOFING.md`.) **Done-when:** an audit
  table with an arena-ownership verdict per guarded path. **Artifact:**
  `benchmark/PHASE_4C_host_access_audit.md`.

**4.D — Nodal/cell-centered indexing audit**
- **Goal:** the restart node-fab OOB bug (node-centered fab sized with cell dims)
  likely has siblings on the high-side boundary / restart path. **Steps:** test
  restart-after-regrid; restart with elastic-on/thermal-off and the reverse;
  multi-level node fields; high-side boundary nodes; checkpoint↔plotfile roundtrip
  for every registered component. **Done-when:** each path passes the Phase 1 gate
  or has a filed defect. **Artifact:** `benchmark/PHASE_4D_nodal_indexing_audit.md`.

**4.E — Deep-AMR: supported-or-abort**
- **Goal:** no silent-wrong-results state. The fast GPU regime is uniform /
  `max_level=1`; deep subcycling AMR is launch-poison but must not quietly give
  wrong answers if a user enables it. **Steps:** either deep AMR passes the Phase
  1 validator, or it aborts with a clear message. Document the "recommended fast
  GPU regime" vs the "supported correctness regime." **Done-when:** deep AMR is
  green-or-aborts, documented. **Artifact:** `benchmark/PHASE_4E_amr_support.md`.

**4.F — Solver-tuning safety guard**
- **Goal:** detect when a loosened/skipped elastic solve has gone stale or
  under-resolved. **Steps:** watch for nonlinear residual not decreasing, stress
  discontinuities after a skipped solve, pressure/burn-rate phase-lag, local
  stress extrema drifting more than global norms, `psi_floor` stiffness jumps at
  the interface. Tie these to the Phase 1 budget's trajectory + SOLVER-HEALTH
  checks so 2.B can't pass a physically-stale config. **Done-when:** the guard
  conditions are encoded as budget checks. **Artifact:** folded into
  `physics_budget.md` + `PHASE_4F_solver_safety.md`.

**4.G — Fix `DeviceErrorFlag` stream-pool race**
- **Goal:** close a live, unfixed defect found by the 2026-07-02 code audit
  (`benchmark/GPU_AUDIT_20260702.md` §2): `DeviceErrorFlag::value()`
  (`Util.H:68-78`) syncs only the *current* CUDA stream
  (`Device::streamSynchronize()`), but `MFIter` cycles boxes across a pool of
  streams — the same mechanism behind the historical `interpolation()` `tmpfab`
  UAF. `SetDeviceError` writes from earlier-streamed boxes may still be in
  flight when `AbortIfDeviceError` reads the flag right after the MFIter loop
  (`Newton.H:121-124,228-231` — unconditionally compiled, always live in
  production; `Elastic.cpp:675,889` — `AMREX_DEBUG`-gated). Writes are atomic
  so there's no corruption; the failure mode is a **silent false negative** — a
  NaN/invalid-`kinvar` error on a non-last-streamed box goes undetected on
  multi-box runs. **Steps:** add `amrex::Gpu::streamSynchronizeAll()` before
  the flag read, cleanest inside `Util::AbortIfDeviceError` itself
  (`Util.H:148-151`) so all four call sites are covered by one edit — the same
  idiom already used for the `dsol_mf` UAF fix at `Newton.H:459`. Cost is one
  full-device sync per `prepareForSolve`/diagonal check (once per Newton
  iteration) — negligible. **Done-when:** the fix lands, a multi-box GPU test
  with a deliberately-injected device error (e.g. NaN model coefficient on a
  non-first-streamed box) reliably aborts instead of silently passing; existing
  `tests/GPU/` 9/9 stay green. **Depends:** none (independent, parallel-ok).
  **Artifact:** fix in `src/Util/Util.H` + note in
  `benchmark/GPU_AUDIT_20260702.md`.

**4.H — Wire up GPU CI so it actually runs**
- **Goal:** close the CI gap found by the same audit (§4): both `golden-gpu`
  and `phase1-budget-gate` jobs in
  `.github/workflows/chamber-gpu-correctness.yml` are gated on a
  `has-cuda-runner` repo topic that is never set, so only the `golden-cpu` leg
  executes — the workflow reports green with zero GPU signal, silently, not
  even in a visible "skipped: no GPU" way. **Steps:** either set the
  `has-cuda-runner` topic on a real self-hosted CUDA runner, or switch the
  jobs' `runs-on:` selector to the `[self-hosted, cuda]` label directly per
  task 1.F's documented target state, with an explicit skip-safe message on
  CPU-only runners. **Done-when:** the GPU legs either execute on a real
  runner or print a visible skip, never a silent pass-by-omission. **Depends:**
  overlaps 1.F/2.D scope — coordinate rather than duplicate if those are
  active. **Artifact:** updated `.github/workflows/chamber-gpu-correctness.yml`.

---

### Phase 5 — Scaling envelope (deferred)

Only after single-GPU Fapply is healthy (Phase 3 done). Nice-to-have envelope, not
gating for the release candidate.

---

**5.A — Multi-GPU revisit at large per-GPU domains** *(conditional)*
- **Goal:** find the per-GPU domain size where halo cost amortizes and 2 GPUs beat
  1. **Steps:** currently a measured loss (2 GPUs 3.5–13.3× slower, halo-bound);
  revisit only with much larger per-GPU subdomains + GPU-aware MPI/NVLink
  confirmed, using 2.C's sync-fraction methodology. **Artifact:**
  `benchmark/PHASE_5A_multigpu.md`.

**5.B — H100/H200 rows**
- **Goal:** extend the crossover beyond A100-80. **Steps:** re-run the 2.A/2.B
  matrix on H100/H200 (`build_alamo_nova_3d.sh` already targets `sm_90`).
  **Artifact:** rows added to `PHASE_2A_fair_cpu_baseline.md`.

---

## 5. Task summary

| ID | Phase | Task | Depends | Gate/conditional | Artifact |
| --- | --- | --- | --- | --- | --- |
| 1.A | 1 | Physics-quantity registry & error budget | — | — | `validate/physics_budget.{md,yaml}` |
| 1.B | 1 | Single-entry validation runner | 1.A,1.C | — | `validate/run_validation.sh` |
| 1.C | 1 | Organized output schema | 1.A | — | `validate/README.md` |
| 1.D | 1 | Comparator & verdict report | 1.A,1.C | — | `validate/compare_validation.py` |
| 1.E | 1 | Golden references (CPU-strict + alpha-1) | 1.B,1.D | — | `validate/references/` |
| 1.F | 1 | Make validation the gate | 1.D,1.E | **gates Phases 2–3** | gate mode + CI job |
| 2.A | 2 | Fair CPU baseline | Ph1 | — | `PHASE_2A_fair_cpu_baseline.md` |
| 2.B | 2 | Clean isolated input-lever A/B | Ph1,2.A | budget-gated | `PHASE_2B_input_levers_ab.md` |
| 2.C | 2 | Fapply Speed-of-Light counters | Ph1 | **orders Phase 3** | `PHASE_2C_fapply_sol.md` |
| 2.D | 2 | Lock metric set into CI | 1.F,2.A | — | `perf_regression_track.py` + job |
| 3.A | 3 | Staged register edits + A100 A/B | 2.C | A/B + budget | `PHASE_3A_fapply_register_edits.md` |
| 3.B | 3 | Interior/boundary kernel split | 2.C,3.A | A/B + budget | `PHASE_3B_interior_boundary_split.md` |
| 3.C | 3 | Uniform fast-path + interface mask | 2.C,3.A | if memory-bound | `PHASE_3C_uniform_fastpath.md` |
| 3.D | 3 | Scalarize Matrix4 hot path | 2.C,3.A | if occupancy-bound | `PHASE_3D_matrix4_scalarize.md` |
| 3.E | 3 | Operator reuse across solves | 3.A | — | `PHASE_3E_operator_reuse.md` |
| 3.F | 3 | Re-tune launch/block config | 3.D | — | note in 3.D |
| 4.A | 4 | Sanitizer gate matrix | — | parallel-ok | `PHASE_4A_sanitizer_matrix.md` |
| 4.B | 4 | Async-lifetime audit | — | parallel-ok | `PHASE_4B_async_lifetime_audit.md` |
| 4.C | 4 | Host-access-from-device audit | — | parallel-ok | `PHASE_4C_host_access_audit.md` |
| 4.D | 4 | Nodal/cell indexing audit | — | parallel-ok | `PHASE_4D_nodal_indexing_audit.md` |
| 4.E | 4 | Deep-AMR supported-or-abort | Ph1 | — | `PHASE_4E_amr_support.md` |
| 4.F | 4 | Solver-tuning safety guard | Ph1 | — | `PHASE_4F_solver_safety.md` |
| 4.G | 4 | Fix `DeviceErrorFlag` stream-pool race | — | parallel-ok | fix in `src/Util/Util.H` |
| 4.H | 4 | Wire up GPU CI (`has-cuda-runner` gap) | — | parallel-ok | updated workflow file |
| 5.A | 5 | Multi-GPU at large domains | 2.C | if domain large enough | `PHASE_5A_multigpu.md` |
| 5.B | 5 | H100/H200 rows | 2.A,2.B | — | rows in 2.A |

Critical path: **1.A → 1.{B,C,D} → 1.E → 1.F → 2.A → (2.B ∥ 2.C) → 3.A → 3.{C/D}**.
Phase 4 audits run anytime.

---

## 6. Exit criteria — beta → release candidate

The branch is a release candidate when:

1. **Phase 1 done:** a one-command physics validation suite exists, runs on both
   the A1000 and the A100, and emits diffable bundles; the CPU-strict golden and
   alpha-1.0 GPU references are recorded; the gate is wired and mandatory.
2. **2.A done:** a *fair* CPU baseline replaces the inflated 22.9× headline; both
   numbers are reported with node/GPU-hours.
3. **2.B done:** the Phase C input levers are attributed individually and each
   shipped lever has a budget-PASS — no stacked, unvalidated tuning.
4. **2.C done:** the Fapply Speed-of-Light counters exist and have chosen the
   Phase 3 order (memory- vs occupancy-bound).
5. **Phase 3 executed *or* explicitly closed:** at least the staged register edits
   (3.A) measured on A100 with a budget pass; the bigger kernel rewrites (3.C/3.D)
   either landed with measured wins or documented "no win" against the counter.
6. **Phase 4 baseline:** the sanitizer matrix (4.A) is green and the four audits
   (4.B–4.E) have a verdict table — no known silent-wrong-results path.

Multi-GPU (5.A) and H100/H200 (5.B) are **beyond the release candidate** —
envelope extension, not gating.

---

## 7. Cross-references

- Entry point / process: `benchmark/READ_FIRST_NEXT_STEP.md`
- Build matrix / branch policy / IC-BC safety: `benchmark/GPU_BRANCH_GUIDE.md`
- Wins log (settled fact, do not re-litigate): `benchmark/SUCCESS_BOOK.md`
- Standing metric set + perf-regression harness: `benchmark/PERF_TRACKING.md`
- Frozen alpha-1.0 reference: `benchmark/ALPHA1_BASELINE.md`
- **History archive (superseded v2 plan, Phase A–C findings, Phase 1–5 reports,
  test-suite fixes, old debug harnesses):** `benchmark/archive/README.md`
- Validation suite home (Phase 1 deliverable): `benchmark/validate/`
- Full-`src/` code audit (2026-07-02): `benchmark/GPU_AUDIT_20260702.md` — live
  `DeviceErrorFlag` defect (task 4.G), CI gap (task 4.H), Fapply perf anatomy
  supporting task 3.A, dormant landmines for anyone extending the GPU closure
- Structural plan (2026-07-03): `benchmark/GPU_STRUCTURAL_PLAN_20260703.md` —
  candidate new tasks 3.G (stored-stencil fine-level operator), 3.H (Fsmooth
  fusion), PF.1–PF.5 (phase-field launch/sync structural fixes incl. landing the
  measured `codex/gpu-pf-structural-speedups` branch), and the multi-GPU
  prerequisite list that amends 5.A's framing

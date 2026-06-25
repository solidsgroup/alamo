# Phase 5.4: chamber-gpu Branch Definition-of-Done

Roadmap: `~/Desktop/GPU-OPT-ROADMAP.txt`, Phase 5, "Branch definition-of-done"
checklist (`chamber-gpu` is never merged to master; this is the branch's own
done-bar, not a merge gate). See `benchmark/GPU_BRANCH_GUIDE.md` for the
navigable index this checklist points into.

Status legend: **DONE** (evidence exists and was verified), **PARTIAL** (some
evidence exists, gap is named), **PENDING** (work tracked but not yet landed —
points at the expected artifact).

## Checklist

### 1. Golden compare passes (no-fast-math) on coarse config

**Status: DONE** (descoped to coarse-only 2026-06-24, explicit user
directive — see note below)

- Coarse: **DONE**. `benchmark/G0_BASELINE_OF_RECORD.md` records
  `baseline_suite.py check` passing for `canonical_step1`, `canonical_step2`,
  `eta_expression_step1` against `gpu_strict` (no-fast-math) at local/coarse
  resolution. Harness: `benchmark/golden_compare_flame.sh`,
  `benchmark/baseline_suite.py`.
- **Saturating-3D golden compare: descoped, no longer part of this item.**
  Originally PENDING — the 3D NOVA crossover runs
  (`benchmark/PHASE3_R3_crossover.md`, 128³/256³, D3 = WIN @ single) report
  wall-clock and stability but never ran a no-fast-math golden compare at
  scale. Per explicit user directive (2026-06-24): CPU-vs-GPU correctness
  parity at scale is not the current problem; performance and stability (the
  elastic cross-box transfer defect under D1) were. **The elastic defect is
  RESOLVED 2026-06-25** (GPU cross-stream race on a per-box temp in
  `interpolation()`, fixed with `tmpfab.elixir()`; see `GPU_BRANCH_GUIDE.md`
  D1 section). Not planned unless a correctness regression is suspected at scale.

### 2. No host-loop writes into device-arena fabs reachable from any documented-supported input

**Status: DONE**

- `docs/gpu_safe_ic_bc_matrix.md` enumerates every IC/BC path and its status.
  Supported (device kernel) paths: `IC::BMP`, `IC::Constant`, `IC::Expression`
  (scalar), `BC::Constant`, `BC::Operator::Elastic::Constant`.
  Guarded-unsupported paths (`IC::PNG`/`PSRead`/`Trig`/`Laminate`,
  `IC::Expression` vector, `BC::Operator::Elastic::Expression`) still use
  `LoopConcurrentOnCpu`, but CUDA builds **abort before** entering the
  host-loop write — verified by `benchmark/test_gpu_guarded_ic.sh`.
  "Documented-supported input" = BMP/Constant IC + Constant elastic BC per the
  matrix; that surface has zero host-loop device-arena writes.
- Caveat: this is "guarded," not "ported" — the unsupported paths are fenced
  off, not device-implemented. That is sufficient for this DoD item (no
  *reachable* unsafe write under supported inputs) but is a known gap if a
  future input needs PNG/Trig/Laminate ICs on GPU.

### 3. Device aborts/NaN detection active

**Status: DONE**

- Roadmap step 0.1 called for routing device `Util::Abort` through
  `amrex::Abort()`/`__trap()` instead of a silent no-op, and restoring NaN
  detection via a per-cell scratch-fab flag reduced once per kernel launch.
  `benchmark/G0_BASELINE_OF_RECORD.md` records the corrected-build status as
  satisfied (G0 correctness/build requirements passed); the guarded-IC abort
  behavior in item 2 is direct evidence the device abort path fires (CUDA
  builds abort, not silently corrupt, on unsupported IC/BC entry).
- Evidence pointer for the CI-enforced form of this gate: the NaN-flag
  assertion smoke described in task 004
  (`.github/workflows/chamber-gpu-correctness.yml` +
  `benchmark/ci_golden_compare.sh`) — see item 5 below for that artifact's
  landing status.

### 4. CPU regression suite green for all integrators

**Status: DONE — fixed and committed 2026-06-23 (commit `2cacb50dd`).**

- **UPDATE 2026-06-23:** Both Flame defects are fixed and committed
  (`src/Integrator/Flame.cpp`, `src/Integrator/Flame.H`, commit `2cacb50dd`).
  Bug #1 (`L_mf` null write at `thermal.on=0`): `L_mf` is now registered
  unconditionally (`Flame.cpp:130-135`); `L` (always consumed by the eta
  evolution) is validated regardless of `thermal_on`, while `K`/`rho`/`cp`
  (purely thermal, legitimately NAN for burn-rate-only propellant models) are
  only validated when `thermal_on`. Bug #2 (`model_prop` arity abort): the
  ctor and `UpdateModel` now branch on `homogenized` — homogenized configs
  parse `model_prop`/`model_void`/`model_casing` as before; resolved configs
  (e.g. `SCPSpheresElastic`) now parse and blend `model_ap`/`model_htpb` by
  the species field `phi`, restoring the original mesoscale rule-of-mixtures
  behavior instead of aborting at parse time.
- Re-ran `scripts/runtests.py --dim=2` after rebuilding `bin/alamo-2d-g++`:
  **118 run, 92 verified, 0 failed** (was 113 run, 89 verified, 5 failed). The
  3 `SCPSandwich` + 2 `SCPSpheresElastic` cases that previously
  segfaulted/aborted now run to completion and verify clean. No new failures
  introduced (de-virtualization conclusion from 2026-06-21 stands).
- **Reproducibility note:** hit an unrelated local environment break before
  this re-run — a stale user-level `numpy==2.4.6` (in
  `~/.local/lib/python3.12/site-packages`, ahead of the apt-installed
  `python3-numpy==1.26.4` on `sys.path`) was incompatible with the
  apt-installed `yt`/`matplotlib`/`pandas` (built against numpy 1.x), so
  every yt-based test's "Checking result" step failed with an `ImportError`
  regardless of source changes. Fixed via
  `pip3 install --user --break-system-packages "numpy<2"`. Anyone re-running
  this suite locally should check for the same shadowing if "Checking
  result" fails universally.

- **UPDATE 2026-06-21 (task 001 result landed,
  `docs/agent_plans/20260621-gpu-phase4-5/results/001-RESULT.md`):** the 2D
  suite is **113 run, 89 verified, 5 failed**. The de-virtualization concern is
  resolved **negative** — every genuine `Solid`-dispatch / elastic integrator
  passed (Eshelby, PlateHole, Rubber*, UniaxialTension, Solid, ThermoElastic,
  VoronoiElastic, DynamicBar, …). The 5 failures are 2 flame test cases, both
  confirmed chamber-lineage Flame bugs (master = merge-base `4e80a3e68` runs
  both cleanly), **neither de-virtualization**:
  - `SCPSandwich` (×3): SIGSEGV at `Flame.cpp:677` (`L_out(i,j,k)=L`) — `L_mf`
    registered only under `if(thermal.on)` (Flame.cpp:180) but written
    unconditionally; crashes when `thermal.on=0`.
  - `SCPSpheresElastic` (×2): parse-time abort (`ParmParse.H:1297 query_exactly`
    on `elastic.model_prop` requiring exactly 2 of {lambda,mu,E,nu,kappa}; input
    supplies none) — aborts in the Flame ctor before any solve.
  This item flips to DONE once those two defects are fixed (register `L_mf`
  unconditionally / guard the write; resolve the `model_prop` arity contract).
  See `benchmark/PHASE4_R4_dispatch.md`.
- Tracked by sibling task 001 (`docs/agent_plans/20260621-gpu-phase4-5/tasks/001-phase4-1-cpu-regression.md`),
  which runs `scripts/runtests.py` over all 53 integrators (2D required, 3D
  best-effort) on `chamber-gpu` to confirm the de-virtualization of
  `Solid`/`BC`/`Operator` (removing `virtual` for nvcc) did not change any
  CPU integrator relying on runtime `Solid*` dispatch.
- Evidence pointer once landed: `docs/agent_plans/20260621-gpu-phase4-5/results/001-RESULT.md`
  and `benchmark/PHASE4_R4_dispatch.md` (task 002, which folds the 4.1 verdict
  into the R4 framework-dispatch report). Until those land, this item has no
  pass/fail data and must not be marked DONE.
- Partial circumstantial evidence: `benchmark/README.md` "Status" section
  states the chamber/Flame port is "bit-for-bit identical to before the port"
  for the center-bore CPU run, and `src/GPU/IntegratorPolicy.mk` scopes the
  GPU build to `flame` only so other integrators' *build* is unaffected — but
  build-unaffected is not the same as a verified-green *test* suite, hence
  PENDING rather than DONE.

### 5. R3 win (or documented no-win + recommendation) recorded on the branch

**Status: DONE**

- `benchmark/PHASE3_R3_crossover.md` records **D3 = WIN @ single**
  (medium-high confidence): single A100 beats a 16-rank CPU node ~39× at 128³
  and ~70× at 256³, advantage growing with size (saturating-regime signature).
  "WIN @ scale" (multi-GPU) is explicitly **not** claimed — current 2-GPU data
  at small per-GPU domains is strongly negative (~9.5× slower, halo-bound) and
  the report flags this as open, not resolved.
- The only caveat on this DONE: the CPU baseline used is 16-rank, not the
  64-rank full node the current SLURM script requests; the report itself
  argues the margin (39-70×) is too large for a 4× CPU rank increase to
  overturn the verdict, and lists the 64-rank confirmation as future
  hardening, not a blocker to recording the win.

## Summary table

| # | Item | Status | Evidence |
| --- | --- | --- | --- |
| 1 | Golden compare (no-fast-math): coarse | DONE | `G0_BASELINE_OF_RECORD.md`. Saturating-3D golden compare descoped 2026-06-24 (user directive) — not required |
| 2 | No reachable host-loop device-arena writes | DONE | `docs/gpu_safe_ic_bc_matrix.md`, `benchmark/test_gpu_guarded_ic.sh` |
| 3 | Device aborts/NaN detection active | DONE | `benchmark/G0_BASELINE_OF_RECORD.md`; CI form pending in `chamber-gpu-correctness.yml` |
| 4 | CPU regression suite green, all integrators | DONE | commit `2cacb50dd`; `scripts/runtests.py --dim=2` = 118 run, 92 verified, 0 failed |
| 5 | R3 win / no-win recorded | DONE | `benchmark/PHASE3_R3_crossover.md` (D3 = WIN @ single) |

## Overall branch status

**5 of 5 DONE.** Item 4 flipped to DONE 2026-06-23 (commit `2cacb50dd`): both
chamber-lineage Flame defects (a `thermal.on=0` null write, a `model_prop`
parse-arity abort) are fixed, and `scripts/runtests.py --dim=2` is fully
green (0 failed). Item 1 flipped to DONE 2026-06-24 by descoping the
saturating-3D golden compare (user directive — not the current priority).
Per branch policy this checklist governs "branch done," never "merged" —
`chamber-gpu` stays permanently unmerged regardless of its completion state.

## Next actions (updated 2026-06-24)

All 5 checklist items are DONE; the checklist itself is closed. Remaining
work on the branch is tracked outside this DoD list: the GPU elastic
cross-box transfer defect (D1) is **RESOLVED 2026-06-25** — root-caused to a
GPU cross-stream race on a per-box temp in `interpolation()`, fixed with one
line (`tmpfab.elixir()`), verified end-to-end
(`benchmark/elastic_sensitivity_20260621/GPU_ELASTIC_DEBUG_PLAN.md` SOLVED
banner). The Newton-damping fix port (currently only in worktree
`/tmp/alamo-newton-damping`) is still open. 512³ CPU baseline and multi-GPU scaling are
explicitly descoped/back-burnered (see `docs/llm/CURRENT.md`) and are not
DoD blockers.

# Physics Validation Suite — `benchmark/validate/`

**Status: SCAFFOLD (2026-06-30).** This directory is the home of **GPU Roadmap v3
Phase 1** — the physics-error budget and the one-command validation suite. The
plan and task breakdown (1.A–1.F) are in `benchmark/GPU_ROADMAP_V3.md` §4. This
README is the schema spec (roadmap task **1.C**); the rest is built per the tasks
below.

## Why this exists (the v3 governing rule)

> No optimization ships — input lever or kernel rewrite — without a physics-error
> budget that proves the coupled burn / pressure / stress behavior is preserved
> within a stated tolerance.

A loosened tolerance, a skipped solve, or a register-shaved kernel can look faster
and stable while silently shifting the pressure history, burn area, or stress
localization the chamber model exists to predict. This suite is the instrument
that catches that. Counters prove a change is *faster*; this suite proves it is
still *the same simulation*.

## Two entry points, not one auto-detecting script

Local and NOVA runs use **separate scripts**, not a single hardware-detecting
`run_validation.sh`. NOVA's Slurm/module environment (account, partition,
GRES names, MPI launch mode — see `benchmark/NOVA_SLURM_RUNBOOK.md`) is
specific enough that folding it into a local runner adds branching without
removing any real duplication — the two environments don't share a launch
path. Both still emit the same bundle schema below.

- `run_validation_local.sh` — local A1000 (sm_86) only. Runs directly via
  `mpiexec` against `bin/alamo*` binaries already built in this checkout.
  Defaults to the **strict / no-fast-math** build for validation; `--fast`
  opt-in for a perf-labelled smoke row.
- `run_validation_nova.slurm` (+ `cases.manifest`'s NOVA-sized case set) —
  submitted by the user via `sbatch` per `NOVA_SLURM_RUNBOOK.md` conventions
  (account `brunnels`, partition `nova`, `GPU_TYPE` env var). Not run by an
  agent — cluster access is the user's.

Both emit one organized output bundle per case, idempotently, with no
per-host edits to the bundle schema itself.

## Output bundle schema (task 1.C — the contract for storage + comparison)

```
benchmark/validate/runs/<UTCstamp>_<gitSHA>_<host>_<device>/
  manifest.json      # host, GPU, driver, git SHA, build flags, binary, input hashes, np
  <case>/
    metrics.json     # every physics_budget observable, as {value, tolerance_class}
    thermo.dat       # raw scalar history (copied)
    region_times.txt # TinyProfiler regions (elastic vs phase-field split)
    field_norms.json # per-component L2/L-inf of the final plotfile
    run.log          # full stdout
```

`<UTCstamp>` (UTC, sortable) + short `<gitSHA>` + `<host>` + `<device>` make every
bundle self-describing, chronologically sortable, and collision-free across
machines. `metrics.json` is the canonical comparison surface.

## Tolerance classes (defined in `physics_budget.md`, task 1.A)

- **CORRECTNESS** — final per-component field L2/L∞ norms + converged residual;
  must match the CPU-strict golden to ~1e-6 relative (no-fast-math). A CORRECTNESS
  fail ⇒ overall FAIL.
- **ENGINEERING-TRAJECTORY** — chamber pressure history, volume, burn area, total
  mdot, burn-front displacement, max/mean von Mises, max principal stress, elastic
  energy; compared over the *whole* run (max relative deviation + phase-lag), not
  just the endpoint, within a stated physical %.
- **SOLVER-HEALTH** — Newton iters/solve, MLMG V-cycles, bottom-solver iters,
  residual-history monotonicity; regression-tracked, non-gating.

## To build (roadmap Phase 1 tasks)

| Task | Deliverable (in this dir) | Status |
|------|---------------------------|--------|
| 1.A | `physics_budget.md` + `physics_budget.yaml` (observable registry + thresholds) | ✅ draft, pending domain review (see checklist in `physics_budget.md`) |
| 1.B | `run_validation_local.sh` + `run_validation_nova.slurm` + `cases.manifest` (separate local/NOVA entry points, see above) | TODO |
| 1.C | **this `README.md`** (bundle schema spec) | ✅ scaffold |
| 1.D | `compare_validation.py` (per-observable verdict, trajectory compare) | ✅ done (`--selftest` passes: known-good→PASS, detuned→FAIL) |
| 1.E | `references/{cpu_strict,gpu_alpha1_local}/` + `BASELINE_DELTAS.md` (local leg only; NOVA leg needs user to submit `run_validation_nova.slurm`) | ✅ local leg done (2026-07-02) — `canonical_2d_elastic` bit-identical PASS; `centre_bore_3d_128_a2` CORRECTNESS FAIL (~1% stress/strain field norm divergence, GPU reduction non-associativity), but all ENGINEERING-TRAJECTORY observables PASS. NOVA leg (`gpu_alpha1_nova_a100`) still TODO. |
| 1.F | `--gate` mode + non-blocking GPU CI job | wired in `.github/workflows/chamber-gpu-correctness.yml` (`phase1-budget-gate` job); `--gate` mode exists in `compare_validation.py`. Job still inert on real CI — gated behind an unset `has-cuda-runner` repo topic (roadmap task 4.H) — but will now find `references/cpu_strict` if that gate is ever opened. |

Build on the existing harness — do not fork it: `benchmark/golden_compare_flame.sh`,
`benchmark/compare_thermo.py`, `benchmark/baseline_suite.py`. The new parts are
hardware auto-detect, the elastic/stress observables the current thermo-only
compare lacks, the trajectory comparison, and this bundle schema.

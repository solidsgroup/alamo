# NOTES — gpu-memory-strategy

Scope-creep sink and glaring-defect log. Nothing here authorizes a code change.
Per PLAN.md §1, a glaring defect is logged here with evidence *first*, then gets
its own tier-3 task folder from `docs/llm/TASK_TEMPLATE.md`.

## Candidate defects (evidence logged, not yet actioned)

| # | Site | Observation | Evidence | Status |
|---|---|---|---|---|
| N1 | `src/Integrator/Flame.cpp:1127` | `reduce_data.value(reduce_op)` forces a device sync **per box, per level**, inside the `MFIter` loop, to accumulate `chamber.{volume,area,mdot}` on the host | Source read 2026-07-27 | Phase 2 target, not a defect per se — the comment at `:1112` shows it was a deliberate race fix. Device-resident accumulator is the replacement. |
| N2 | `src/Integrator/Flame.cpp:697-701` | Five separate `ParallelDescriptor::ReduceReal{Max,Min}` calls back to back (`thermo_max_temp`, `thermo_mdot_max`, `thermo_heatflux_max`, `thermo_L_max`, `thermo_eta_min`) | Source read 2026-07-27 | Batchable into one allreduce. Phase 2. |
| N3 | `src/Integrator/Integrator.cpp:1246`, `:1263` | `MFIter mfi(grids[...], dmap[...], true)` — tiling hardcoded `true` in the integrate path. If not `TilingIfNotGPU()`-guarded, CPU tiling runs on GPU as pure overhead | Source read 2026-07-27 | **Unverified.** Check the guard before calling it a defect. Phase 3 lead. |
| N4 | `benchmark/baseline_suite.py` / gpu_strict | `rod_and_tube_step2` GPU golden reference is stale (pre-existing, carried from earlier work) | Prior sessions | Blocks trusting the gpu_strict leg. Fix or quarantine before Phase 1 leans on it. |
| N5 | `benchmark/status.sh`, `benchmark/ci_golden_compare.sh` | Not concurrency-safe: shared log path, shared `baseline_runs/` tree, shared `bin/alamo-2d-g++`. Two overlapping runs → spurious RED (`canonical_step1/cpu failed with 131`, SIGQUIT) | Observed 2026-07-27; solo re-run `EXIT=0` | Harness hygiene. Cheap fix: per-run log/dir naming or a lockfile. |
| N6 | `src/Integrator/Integrator.cpp:1229` vs `src/Integrator/Flame.cpp:706` | `IntegrateVariables` (which recomputes and allreduces `chamber.{volume,area,mdot}`) is gated on `thermo.interval`, a **diagnostics-output** knob. `chamber.model.Advance(...)` consumes those scalars **every** step regardless. At `amr.thermo.int > 1` with `variable_pressure = 1`, the pressure ODE integrates stale mass flux. | Source read + `Integrator.cpp:122-123` (default 1), `input:17` (`amr.thermo.int = 1`) — 20260727-phase05-two-rank-probe RESULT §1 | **Latent, not live.** Every current deck sets 1 and the default is 1. Physics-input freshness should not be coupled to an output-interval knob. Own tier-3 folder if actioned; do not fix inline. |

## Deferred ideas

(empty)

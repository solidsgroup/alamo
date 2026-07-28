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
| N7 | `src/Integrator/Hydro.cpp:332` | `Hydro::TimeStepComplete` does `if (dynamictimestep.on) DynamicTimestep_Update(); return;` — the `return` makes lines 334-344 **dead code**: the three `ParallelDescriptor::ReduceRealMax(c_max/vx_max/vy_max)` and the CFL `SetTimestep(new_timestep)` never execute | Source read 2026-07-27 | Not a live bug (nothing wrong runs) but Hydro's own CFL timestep control is disabled in favour of `DynamicTimestep_Update`. Either intentional and the block should go, or a regression. Needs the Hydro assessment (PLAN §18) first — do not delete on sight. |
| N8 | `src/Integrator/Hydro.cpp`, all 10 `MFIter` sites | Tiling hardcoded `true`/`false`; `TilingIfNotGPU()` used 0 times. Elastic.cpp 8/13, Flame.cpp 1/10, Integrator.cpp 0/3 | Source read 2026-07-27 | Same class as N3, systemic in Hydro. CPU tiling on GPU is pure overhead (PLAN §10). Blocked behind the PLAN §18 question of whether Hydro runs on GPU at all. |
| N6 | `src/Integrator/Integrator.cpp:1229` vs `src/Integrator/Flame.cpp:706` | `IntegrateVariables` (which recomputes and allreduces `chamber.{volume,area,mdot}`) is gated on `thermo.interval`, a **diagnostics-output** knob. `chamber.model.Advance(...)` consumes those scalars **every** step regardless. At `amr.thermo.int > 1` with `variable_pressure = 1`, the pressure ODE integrates stale mass flux. | Source read + `Integrator.cpp:122-123` (default 1), `input:17` (`amr.thermo.int = 1`) — 20260727-phase05-two-rank-probe RESULT §1 | **Latent, not live.** Every current deck sets 1 and the default is 1. Physics-input freshness should not be coupled to an output-interval knob. Own tier-3 folder if actioned; do not fix inline. |

## v2 precondition verification — 2026-07-28

PLAN v2 §2 declares P1-P4 blocking. The plan's own accusations were tested
against the tree rather than accepted. Results below; all confirmed.

| Claim | PLAN v2 § | Verdict | Evidence |
|---|---|---|---|
| Main gate self-certified against an inadequate oracle | 5.2 (P3) | **CONFIRMED** | `benchmark/status.sh:33` forced `TIERS=1`. `local_a100_gate.sh:83-86` = launch-blocking smoke with `the_arena_is_managed=1`; `:46` stops before the step-5 elastic solve; `:26` names Tier 2 "the load-bearing tier" and it never ran. Golden was `GOLDEN_MODE=cpu` (`ci_golden_compare.sh:26,84`); the `gpu_strict` path at `:119` was unreachable. |
| Baseline captured from a dirty tree | 5.1 (P2) | **CONFIRMED** | Capture provenance block already self-reports `local_dirty_files=574` against `local_head=cb1a3071e`. |
| Vendored AMReX already defaults `the_arena_is_managed=false` | 3.3 T1 | **CONFIRMED** | `ext/AMReX-Codes/amrex/Src/Base/AMReX_Arena.cpp:59`. |
| `AbortIfDeviceError` = two unconditional full-stream syncs per Flame level per step | 5.3, 5.4 | **CONFIRMED** | `Util/Util.H:156` `Gpu::streamSynchronizeAll()` then `:157` `flag.value()`. Call sites `Flame.cpp:861` and `:915`, both unguarded, per level per `Advance`. |
| Per-step device error flag allocates a `DeviceScalar` with H2D/D2H | 8 (T4) | **CONFIRMED** | `Util/Util.H:77` `amrex::Gpu::DeviceScalar<int> m_flag`; constructed unguarded at `Flame.cpp:747`. |

### Consequences

**P3 actioned.** `status.sh` split into a fast orientation mode and `FULL=1`
(Tier 2 memcheck + `gpu_strict` golden). Fast legs relabelled to name what they
execute — `golden-cpu`, `smoke-pre-elastic` — and the fast mode prints an
explicit "not assurance" notice. The header lists the §5.2 coverage that
*neither* mode has yet (pressure divergence, regrid, two ranks, restart) so a
`FULL=1` PASS is not over-read either. Missing nvcc reports `BLOCKED`, not
`FAIL`: an unrun leg is not a correctness result.

**A prior claim is withdrawn.** `20260727-phase0-baseline/results/RESULT.md`
§N2 concluded "T1 passes" from an empty device-arena crash set. Given
`AMReX_Arena.cpp:59`, the device arm is the AMReX default and the managed arm
was the synthetic one; an empty crash set does not distinguish "no managed
allocations" from "managed allocations that happen to work" (§3.3 T1). The
measurement stands; the T1 conclusion does not.

## Candidate defects, second batch (2026-07-28)

| # | Site | Observation | Evidence | Status |
|---|---|---|---|---|
| N9 | `src/` tree-wide | ALAMO **never selects an arena**. Two occurrences of "arena" in `src/`: one error string (`Util/Util.H:32`) and one call, `src/Test/BC/Constant.H:43` `The_Managed_Arena()`. No `SetArena`, no pinned/async/comms arena, nothing outside `src/Test/`. | `grep -rni arena src/` 2026-07-28 | Reshapes PLAN §7: the data-class-to-arena table is **new plumbing**, not a flip of existing call sites. Phase 1a scope is larger than "mechanical flip" implies. |
| N11 | `src/Operator/Elastic.cpp:461` → `src/Numeric/Stencil.H:1592` | **LIVE 3D GPU correctness defect.** `Operator::Elastic<1>::Diagonal()` reads out of bounds during the first elastic solve. 2,497 × `Invalid __global__ read of size 8 bytes`, then `cudaErrorLaunchFailure (error 719)`. Faulting address is **18,320 bytes *before*** the nearest allocation — reading below the fab, not past its end. | `FULL=1 bash benchmark/status.sh` 2026-07-28, deck `input_3d_centre_bore_128_a2`, log `benchmark/_a100_gate_20260728_125617/tier2_memcheck.log` | **Needs its own tier-3 folder. Do not fix inline.** Detail below. |
| N10 | `Flame.cpp:747,861,915` vs `Elastic.cpp:200-205,391,517` | Elastic wraps both the `DeviceErrorFlag` and the `AbortIfDeviceError` check in `#ifdef AMREX_DEBUG` / `#ifdef ALAMO_GPU`, so they compile out of release. Flame's are **unguarded** and ship in every build. | Source read 2026-07-28 | §5.3 asks whether this "compiles out"; in Elastic it already does, in Flame it does not. **Fix = apply Elastic's existing guard to Flame, NOT delete the sync** — `Util.H:150-155` documents that the sync closes a real cross-stream read race this branch already fixed once. Deleting it reintroduces a known bug. |

### N11 detail — what the first honest gate run found

The P3 fix was validated by running `FULL=1`. It went red on its first
execution. That is the finding: the branch has been shipping a green board while
the load-bearing tier was never run.

```
device-lint          : PASS
golden-gpu-strict    : PASS
a100-memcheck-tier2  : FAIL
```

`golden-gpu-strict` passing is worth recording separately — it independently
re-confirms the §L4 `rod_and_tube_step2` re-recording against a fresh
`--cuda-fp strict` rebuild.

**Call path**, from the memcheck backtrace:

```
Mechanics.H:225  TimeStepBegin -> solver.solve(...)
  MLMG::solve -> MLMG::prepareForSolve
    Operator.cpp:515  Operator<Grid::Node>::prepareForSolve()
      Elastic.cpp:441  Elastic<1>::Diagonal() ParallelFor
        Elastic.cpp:461 (device frame)
          Stencil.H:1592  Interpolate::CellToNodeAverage
```

**Mechanism, high confidence but not yet proven by a fix.** `Stencil.H:1587-1592`
computes `ilo/jlo/klo = (stencil[d] == StencilType::Lo ? 0 : 1)` and reads
`f(i-ilo, j-jlo, k-klo, m)`. With no stencil argument the default is all-Central,
so `ilo=jlo=klo=1` and the read steps to `(i-1,j-1,k-1)` unconditionally.

`Elastic.cpp:461` calls `CellToNodeAverage(psi, i, j, k, 0)` **without** a
stencil argument — while the same lambda computed the boundary-aware stencil at
`:444-445` (`Numeric::GetStencil(i, j, k, stencilbox)`) and correctly passed it
to `Gradient_Diagonal` at `:448`. The psi interpolation ignores the stencil the
surrounding code went to the trouble of building. Reading below the allocation is
consistent with stepping to `i-1,j-1,k-1` at the low corner.

**It is systemic, not a one-line slip.** Of 19 `CellToNodeAverage` call sites in
`src/`, exactly **one** passes a stencil: `Solver/Nonlocal/Newton.H:1411`. All
three in `Elastic.cpp` (`:255` Fapply, `:461` Diagonal, `:669`) omit it, and
`PhaseFieldMicrostructure.cpp:270-273` has `//, sten);` commented out — someone
removed the argument deliberately at some point. The single site that passes it
is evidence the requirement was known.

Do not assume the fix is "pass `sten`". Growing psi's ghost region or clamping
the read are also candidates, and the three Elastic sites may not want the same
answer. That adjudication belongs in the tier-3 folder with a memcheck-clean
gate as its oracle.

**Provenance.** The only dirty `src/` files are `Flame.{cpp,H}`; the fault is in
Elastic/Stencil, which they do not touch, so this is almost certainly
pre-existing on the branch rather than introduced by uncommitted work. Not
proven — no clean-tree rebuild was run to confirm.

**Blocking status.** This is a correctness defect in the elastic path, in 3D,
psi-gated (`if (m_psi_set)`). PLAN v2 §2 makes oracle integrity a precondition
for Phases 1-3; a red oracle is a stronger block than a broken one. Phase 1
should not open over it.

## Deferred ideas

(empty)

# Result: Task 001 — Phase 4.1 CPU-semantics regression on chamber-gpu

## Status
DONE (verdict landed). Suite is **NOT fully green**, but the failures are
**not** caused by the nvcc de-virtualization. Two independent Flame-only bugs,
both confirmed chamber-gpu-introduced. De-virtualization (D4=ISOLATE) is
**validated**.

## Headline

`scripts/runtests.py` (2D, g++), captured in
`benchmark/phase4_cpu_semantics_20260621_030514/dim2.log`:

```
113 tests run
89 tests run and verified
5 tests failed
```

The 5 failures are 5 *modes* of just **2 test cases**, both of which run
`alamo.program = flame` (the `SCP*` test directories were long ago repurposed
to flame configs — they are not Solid-crystal-plasticity tests despite the
name):

| Test | Failing modes | Class |
| --- | --- | --- |
| `tests/SCPSandwich` | serial, serial-coverage, parallel | runtime segfault |
| `tests/SCPSpheresElastic` | 2d-serial-short, 2d-parallel-long | parse-time abort |

## The de-virtualization question — ANSWERED: clean

Task 001 exists to confirm that removing `virtual` from `Solid`/`BC`/`Operator`
(required for nvcc) did not break CPU integrators that rely on runtime `Solid*`
dispatch. It did not. Every genuine Solid-dispatch / elastic integrator in the
suite **passed**:

- `Eshelby`, `EshelbyDynamics`, `EshelbyFiniteKinematics`
- `Inclusion`, `PlateHole`, `RubberPlateHole`, `RubberPressurizedHole`,
  `RubberWithInclusion`
- `UniaxialTension`, `UniaxialTensionPeriodic`, `DynamicBar`
- `Solid`, `ThermoElastic`, `VoronoiElastic`, `CompositeImpact`, `Suture`,
  `TopOp`, `FractureLimestone`, `FracturePFCZM`
- `SCPThermalSandwich`, `SCPThermalVoid` (flame + elastic-disabled, thermal on)

If de-virtualization had corrupted `Solid` vtable dispatch, these would be the
first to fail. They pass. The 2 failing cases contain **no Solid-dispatch path
that the de-virtualization touched** (one has elastic disabled entirely; the
other aborts before any solve). See `benchmark/archive/PHASE4_R4_dispatch.md` for how
this folds into the D4 verdict.

## Root cause #1 — SCPSandwich (3 modes): null write in Flame::Advance

- **Crash:** SIGSEGV at `src/Integrator/Flame.cpp:677`, `L_out(i, j, k) = L;`
- **Why:** `L_out = L_mf.Patch(lev, mfi)` (Flame.cpp:651) is an empty `Array4`
  (`p = 0x0`) because `L_mf` is registered **only inside** the
  `if (value.thermal.on)` block — `RegisterNewFab(value.L_mf, …)` at
  Flame.cpp:180, guard opens at Flame.cpp:150. SCPSandwich sets
  `thermal.on = 0`, so `L_mf` is never allocated, yet the Advance kernel writes
  the mobility to it unconditionally (the sibling thermal read on Flame.cpp:665
  *is* guarded: `thermal_on ? temp(i,j,k) : NAN`; the `L_out` write is not).
- **Introduced, not pre-existing:** master (= merge-base `4e80a3e68`) has **no
  `L_mf` field at all** — it computes `L` locally and never stores it. The
  `L_mf` mobility diagnostic was added on the chamber lineage with a
  thermal-only registration but an unconditional write. The master binary runs
  the identical SCPSandwich flame config to completion (STEP 1/STEP 2, exit 0);
  the chamber-gpu binary segfaults at STEP 1.
- **Trigger is `thermal.on = 0` + the mobility write.** Not units (master's
  no-suffix input crashes the chamber-gpu binary identically), not Solid
  dispatch, not the elastic path.
- **Fix options (owner decision):** register `L_mf` unconditionally (it is the
  Allen–Cahn mobility, needed regardless of thermal — registering it outside
  the thermal block is the natural fix), or guard the `L_out(i,j,k) = L` write
  behind `if (L_out)` / `thermal_on`.

## Root cause #2 — SCPSpheresElastic (2 modes): parse-time arity abort

- **Abort:** `ABORT ./src/IO/ParmParse.H:1297 (query_exactly)` —
  `while reading [model_prop.lambda, model_prop.mu, model_prop.E,
  model_prop.nu, model_prop.kappa]`, "only 2 values are allowed" — raised from
  `IO::ParmParse::queryclass<Integrator::Flame>` inside the Flame constructor
  (`src/Integrator/Flame.cpp:28`). Fails **before any timestep**.
- **Why:** chamber-gpu's Flame elastic model (`elastic.model_prop`, type
  `Model::Solid::Finite::NeoHookeanPredeformed`) parses with `query_exactly`
  requiring **exactly 2** of the 5 elastic constants
  {lambda, mu, E, nu, kappa}. The SCPSpheresElastic input does not define
  `model_prop.*` at all (it specifies per-phase moduli `model_ap`/`model_htpb`),
  so 0 of 2 are found → hard abort. The uncaught `ParmParseException` propagates
  to `std::terminate`, which the test harness reports loosely as "Segfault."
- **Introduced, not pre-existing:** the master binary parses and runs the same
  SCPSpheresElastic input to completion (exit 0). The stricter
  `model_prop` arity contract is a chamber-gpu change (the elastic model header
  was adopted from `origin/gpu`); it is a **parse/schema** change, not a
  runtime-dispatch regression.
- **Fix options (owner decision):** update the stale input to specify
  `elastic.model_prop` with exactly 2 constants, or make `model_prop` optional /
  derivable from the per-phase models. Needs a call on whether `model_prop` is
  intended to be required.

## Method / evidence (reproducible)

- Suite log: `benchmark/phase4_cpu_semantics_20260621_030514/dim2.log` (run
  03:05–03:34, chamber-gpu working-tree CPU binary `bin/alamo-2d-g++`).
- Crash binary: `bin/alamo-2d-g++` (built 03:04 from chamber-gpu working tree).
  Surviving `Backtrace.0` (03:21) symbolizes to `Flame::Advance`.
- Master baseline: built `master` (= `4e80a3e68`, the chamber-gpu merge-base) in
  an isolated `git worktree` reusing the prebuilt `ext/AMReX-Codes/amrex/
  2d-g++-26.06` (no AMReX rebuild, no writes to the shared tree; worktree bumped
  to `-std=c++20` to match the 26.06 headers). Master runs **both** failing
  configs cleanly (exit 0).
- Exact crash lines from a debug (`-O2 -g`, no LTO) chamber-gpu worktree build
  under gdb: SCPSandwich → Flame.cpp:677 (`L_out` `p=0x0`, `thermal_on=false`);
  SCPSpheresElastic → ParmParse.H:1297 `query_exactly` on `model_prop`.

## Scope notes

- **2D only.** The 3D leg (best-effort per the task) was not completed in this
  capture; the 2D suite is the required bar and is the basis for this verdict.
- **No source/build/git state was modified** by this investigation. All builds
  ran in disposable `/tmp` worktrees that reused the prebuilt AMReX read-only.

## Verdict for DoD item 4 ("CPU regression suite green for all integrators")

**RED, but orthogonal to the GPU-port framework-dispatch decision.** The suite
is not green; the two failures are chamber-lineage Flame bugs (one runtime null
write under `thermal.on=0`, one stale/over-strict `model_prop` parse contract),
both surfaced — not caused — by the de-virtualization audit. They must be fixed
(or the inputs updated) before item 4 can flip to DONE, but they are not a
de-virtualization / `Solid`-dispatch regression. See
`benchmark/archive/PHASE4_R4_dispatch.md`.

# TASK: gpu-memory-strategy
# Folder: docs/agent_plans/20260727-gpu-memory-strategy/

# chamber-gpu-mem: Memory Strategy Migration Plan

**Branch:** `chamber-gpu-mem` (branched from `chamber-gpu` @ cc22ff9ee, 2026-07-27)
**Objective:** Move from "runs on GPU" to "device-resident by construction."
**Character:** Campaign strategy. This document is the campaign spine, not a
single task. Each phase spawns its own `docs/agent_plans/YYYYMMDD-<phase>/PLAN.md`
from `docs/llm/TASK_TEMPLATE.md` before any `src/` edit, per CLAUDE.md.

## Header (this document only)

| Field         | Value                                                       |
|---------------|-------------------------------------------------------------|
| Risk tier     | 0 — docs only; no `src/` change is authorized by this file  |
| Model         | opus (campaign design); per-phase folders route their own   |
| Verification  | judgment (strategy doc); per-phase folders carry oracles    |
| Est. scope    | this PLAN.md + NOTES.md; per-phase folders spawned later    |
| Parallel-safe | yes — docs only, disjoint from all `src/` work              |

Per-phase task folders inherit the template's **Operating rules** verbatim.
Restated here because the campaign spans many sessions:

1. Read only what the phase's Context budget lists. `docs/archive/` is forbidden.
2. Every step's VERIFY passes before acting. On failure: STOP, report, wait.
3. One commit per step. Message: `<area>: <what> (<task-folder>)`.
4. No scope expansion — new ideas go to `NOTES.md`, not into code. The §1
   glaring-defect exception is the *only* carve-out, and it routes through a new
   task folder rather than an inline edit.
5. Missing knowledge → ask. Do not wing it.
6. Tier ≥ 2 → stop at each checkpoint and print the checklist.

## Context budget (campaign level)

Read first: `benchmark/status.sh` output, this PLAN.md
Read: `docs/llm/PLAN.md`, `docs/llm/BUG_PATTERNS.md`,
`benchmark/NOVA_SLURM_RUNBOOK.md`, `ext/AMReX-Codes/amrex/Docs/**/GPU.rst`
Reference only if a step names it: `src/Integrator/Flame.{cpp,H}`,
`src/Integrator/Integrator.cpp`, `src/Solver/Nonlocal/Newton.H`,
`src/Operator/Elastic.*`
Forbidden: `docs/archive/*`, unrelated task folders, the ~/Desktop sweep campaign

## Oracle (campaign level)

Command: `benchmark/status.sh` — all three gates green (device-lint,
golden-compare, a100-sanitizer), run **solo** (see §2 concurrency warning).
Covers: CPU golden correctness, device-pattern lint, runtime-strict smoke.
Does NOT cover: T1-T7 (NOVA measurement only), gpu_strict leg while the
`rod_and_tube_step2` reference is stale, multi-rank correctness (Phase 0.5).

## Closeout (campaign level)

- [ ] Every phase folder has `results/DONE`
- [ ] Metrics ledger populated for all phase boundaries (§12)
- [ ] Patterns exported to `gpu_manual` (§13)
- [ ] `results/RESULT.md` here; `changelog/` appended (never edited in place)
- [ ] One `docs/llm/SESSION_LOG.tsv` line per session

---

## 0. Corrections applied to the source draft

Recorded so the deltas stay visible rather than being silently absorbed.

| # | Draft said | Corrected to | Basis |
|---|---|---|---|
| C1 | Branch `chamber-gpu` | Branch `chamber-gpu-mem` | This campaign's branch |
| C2 | "Elastic/MLMG solver work" out of scope | **Algorithmic redesign** of Elastic/MLMG out of scope; memory residency, arena, and sync work inside them is **in scope**. Glaring algorithmic defects surfaced in Flame/Elastic/MLMG get fixed, not ignored | User ruling; elastic is 95% of GPU wall and excluding it entirely would scope the campaign to ~5% of runtime |
| C3 | Kernel temporaries → `The_Async_Arena()` | Evidence-gated, not default. Already A/B'd on this repo and **rejected**: +2.2% 2D / +4.0% 3D | `gpu_amrex_guide_audit_20260727` |
| C4 | T1 (managed off) is a Phase 1 exit gate | T1 is **NOVA/A100-only**. Local A1000 (8 GB) keeps managed; it cannot validate T1 | A1000 needs managed at 8 GB; `the_arena_init_size` default already OOMs intermittently at init |
| C5 | §14 open questions 1-5 all open | Q1-Q3 partly answered from source (§14). Q4 answered. Q5 answered | Source read 2026-07-27, this session |
| C6 | Testing implicitly local | Explicit access + scheduling model added (§2) | User ruling: `ssh -MN` tunnel, agent-authorized runs, SLURM-bound |

---

## 1. Scope

### In scope

- Memory residency and arena discipline for the chamber timestep loop
- Elimination of per-step host/device synchronization, **including inside
  Elastic/MLMG** — sync and residency are memory-strategy questions wherever
  they occur
- Launch configuration and kernel shape, to the extent memory traffic drives it
- Multi-GPU *structural* decisions only (see §10)
- A repeatable profiling harness with a fixed figure set

### Out of scope

- Algorithmic changes to the Flame solver
- **Algorithmic redesign** of the Elastic/MLMG solver — smoother choice, cycle
  structure, coarsening strategy, Newton formulation. Separate thread.
- `multicomponent/FMA` porting (this plan generates patterns for it, does not
  touch it)
- Multi-GPU performance tuning, halo overlap, load balance
- The propellant parameter-sweep campaign (~/Desktop, sims 030-085) — CLAUDE.md
  hard exclusion

### Glaring-defect exception (C2)

If work under this plan surfaces an obvious algorithmic defect in Flame,
Elastic, or MLMG — wrong reduction scope, a stale hierarchy, a convergence test
that cannot converge, an O(n²) walk where O(n) is available — **fix it**, do not
route around it. Conditions:

1. Log it in `NOTES.md` in this folder with evidence before touching code.
2. It gets its own tier-3 task folder from the template. No inline drive-by edits.
3. Correctness gate is mandatory and non-negotiable: device lint + golden
   compare + compute-sanitizer.

"Glaring" means demonstrable from evidence in hand. A suspicion that an
algorithm could be better is scope creep and goes to `NOTES.md` only.

### Governing rule

> Device owns field data permanently. Host copies are transient, explicit,
> created for I/O, and destroyed immediately. No field crosses the bus inside a
> timestep.

Every phase below is a consequence of this rule. When a decision is ambiguous,
resolve it against this sentence.

---

## 2. Access and scheduling model

**Local (kermit, A1000 sm_86, 8 GB, shared and 50 W-capped):**
correctness only — device lint, CPU golden compare, compute-sanitizer,
`res_usage.sh`. Managed arena stays on locally. **Local timing numbers are not
admissible for T5-T7**; the card is shared and power-capped, and HMM masks
host-pointer defects.

**NOVA (A100/H200, SLURM):** all load-bearing measurement. Access is via the
`ssh -MN` ControlMaster tunnel (`~/.ssh/config` Host `nova`, `ControlPath
~/.ssh/cm-%C`, `ControlPersist 12h`). The agent is authorized to submit and
collect jobs over that tunnel.

**Scheduling discipline.** GPU access is queued, not interactive. Therefore:

- **Batch by phase boundary.** One submission campaign per boundary, capturing
  the entire figure set F1-F10 in a single job where the tooling allows, not one
  job per figure.
- **Compose before submitting.** Every deck, flag, and capture command is
  dry-run locally or on the login node first. A job that dies on a typo costs a
  queue slot and hours.
- **Never block on the queue.** Submit, record the job id, do local work, poll.
  `benchmark/slurm_pending_reason.sh` explains a stuck job.
- **Pin the toolkit version** in the harness and verify flag spellings once per
  cluster change (`NOVA_SLURM_RUNBOOK.md` inventory list).
- Existing entry points: `benchmark/build_alamo_nova{,_3d}.sh`,
  `benchmark/nova_flame_gpu{,_3d,_3d_multi,_3d_diag}.slurm`,
  `benchmark/select_nova_resources.sh`.

**Gate concurrency warning.** `benchmark/status.sh` and
`benchmark/ci_golden_compare.sh` are **not** concurrency-safe: they share one
log path, one `benchmark/baseline_runs/` tree, and one `bin/alamo-2d-g++`. Two
overlapping runs produced a spurious RED on 2026-07-27 (`canonical_step1/cpu
failed with 131` = SIGQUIT, binary swapped by a concurrent `make` under a live
`mpiexec`). Run one at a time. If a gate goes red, re-run solo before believing it.

---

## 3. Target definition

The end state is falsifiable. These are the acceptance criteria for the whole
effort.

| # | Target | Measurement | Threshold | Where measurable |
|---|---|---|---|---|
| T1 | No managed allocations | Runs clean with device arena, managed disabled | Binary | **NOVA only** (C4) |
| T2 | No page-fault migration | nsys unified-memory counters | 0 events | NOVA |
| T3 | Field data does not cross the bus per step | nsys memcpy size summary / steps | Tens of bytes, scalars only | NOVA |
| T4 | No steady-state device allocation | `cudaMalloc` count after warmup | 0 | NOVA; local indicative |
| T5 | Arena footprint stable | High-water mark between regrids | Flat | NOVA |
| T6 | No per-step pipeline drain | Timeline idle fraction inside step loop | Under ~5% | NOVA |
| T7 | Memory-bound kernels near roof | Achieved DRAM bandwidth, top kernels | Band set at Phase 0 | NOVA |

T7 has no absolute number until the Phase 0 baseline exists. Set it then, from
the measured roofline, not from a guess.

**Known T6 starting point:** 2D is already characterized as launch-latency
bound — GPU idle 45.6%, median inter-kernel gap 6.3 µs. 3D is Fapply-compute-
bound (86% of GPU time). The two dimensions will not hit T6 by the same route.

### Explicitly not targets

Zero bus traffic is the wrong goal. Three transfers are legitimate and stay:

- Plotfile and checkpoint writes (bulk, amortized, async, pinned)
- Regrid redistribution (infrequent, device-to-device)
- Convergence and diagnostic scalars (few doubles, batched)

The distinction is not *whether* data moves but whether the movement is
scheduled, explicit, visible in source, and outside the inner loop.

---

## 4. Phase overview

| Phase | Name | Character | Gate to exit |
|---|---|---|---|
| 0.5 | Two-rank correctness probe | Correctness | Multi-rank result matches single-rank |
| 0 | Baseline and gap table | Measurement | Gap table complete, T7 threshold set |
| 1 | Arena hygiene | Mechanical | T1, T4, T5 |
| 2 | Sync elimination | Design | T2, T3, T6 |
| 3 | Launch config and kernel shape | Tuning | T7 |
| 4 | Multi-GPU | Deferred | Out of scope here |

Phases 1 through 3 are strictly ordered. Tuning launch configuration before
removing synchronization tunes around a defect scheduled for removal.

---

## 5. Phase 0.5 — Two-rank correctness probe

**Runs first. Not a performance activity.**

The chamber model reduces regression rate over the burning surface to drive a
scalar pressure ODE. That reduction is inherently global. If it reduces
rank-locally, multi-rank chamber results are already wrong, and the error is
physical rather than numerical.

### Source pre-answer (2026-07-27) — downgrades this risk

Read of the actual call chain:

- `Integrator.cpp:1187-1190` — `IntegrateVariables(...)` then
  `TimeStepComplete(...)`, in that order.
- `Integrator.cpp:1224-1275` — `IntegrateVariables` runs per-box `Integrate`,
  then **`ParallelDescriptor::ReduceRealSum` across ranks** for every
  `thermo.extensives[i]`.
- `Flame.cpp:1113-1130` — per-box `ReduceOps` sum of `(dvol, darea, dmdot)`,
  accumulated host-side into `chamber.{volume,area,mdot}`.
- `Flame.cpp:706` — `chamber.model.Advance(timestep, chamber.mdot,
  chamber.volume, chamber.pressure)` inside `TimeStepComplete`, i.e. **after**
  the MPI reduction.

**Conclusion: the reduction is structurally global, not rank-local.** The
draft's headline risk is most likely already handled. Phase 0.5 therefore
shrinks to confirmation, not investigation.

### Remaining questions

1. Confirm `thermo.extensives[]` is actually true for `volume`, `area`, and
   `mass_flux` — the registration at `Flame.cpp:195-197` and `233-235` must set
   the extensive flag, or the allreduce is skipped and the pre-answer above
   collapses.
2. Does the port survive domain decomposition — ghost fill and BC application
   at rank boundaries?
3. Is GPU-aware MPI actually active at runtime? **Prior evidence says no** — the
   multi-GPU loss was root-caused to blocking comm + no GPU-aware MPI +
   `regrid_int=2` + managed arena. Confirm at runtime rather than assume;
   silent host staging of every ghost cell violates the governing rule and
   pollutes the Phase 0 baseline.

### Exit gate

Two-rank run reproduces single-rank chamber pressure history to solver
tolerance. If it does not, that defect is fixed before Phase 1 — not deferred to
Phase 4.

---

## 6. Phase 0 — Baseline and gap table

**Purpose:** Establish what is actually true, so Phases 1-3 are scoped by
evidence rather than by an assumption that the port has typical defects.

### Prerequisites

- **Regression harness exists and is green.** Answered: `benchmark/baseline_suite.py`
  covers `canonical_step1`, `canonical_step2`, `eta_expression_step1`,
  `rod_and_tube_step2`, gated by `benchmark/ci_golden_compare.sh` (cpu +
  gpu_strict legs). Verified green on this branch 2026-07-27, `EXIT=0`.
  **Known gap:** the `rod_and_tube_step2` GPU golden reference is stale
  (pre-existing, unfixed). Fix or quarantine it before Phase 1 relies on the
  gpu_strict leg.
- **NVTX annotation of the step loop.** An unannotated nsys timeline of a
  multiphysics code is unreadable. Named regions per physics stage convert the
  timeline from a screenshot into a diagnostic. Highest-leverage QOL investment
  in the plan. AMReX `BL_PROFILE` regions can emit NVTX given the right build
  configuration — verify the flag name against the **vendored** AMReX
  (`ext/AMReX-Codes/amrex`, 26.06), not the website and not `ext/amrex` (CPU 25.07).

### Activities

1. Capture baseline nsys profile, managed memory in its current state
2. Flip to device arena on NOVA, capture the crash set (do not fix yet — inventory it)
3. Capture Nsight Compute detail on the top kernels by time
4. Compute the footprint budget (§7)
5. Populate the gap table against T1-T7
6. Set the T7 threshold from the measured roofline

### Gap table

| Target | Current | Gap | Phase | Est. effort |
|---|---|---|---|---|
| T1 managed off | | | 1 | |
| T2 zero page faults | | | 1-2 | |
| T3 per-step bytes | | | 2 | |
| T4 steady-state alloc | | | 1 | |
| T5 stable high-water | | | 1 | |
| T6 idle fraction | 2D: 45.6% idle, 6.3 µs median gap | | 2 | |
| T7 bandwidth vs roof | | | 3 | |

### Decision gate — is chamber the right target?

Phase 0 is also the point to ask whether this work should proceed at all.

Optimize chamber if it gates simulation throughput. If chamber-gpu already turns
around the runs the SRM paper needs, the marginal value of a 2× speedup is low
and `multicomponent/FMA` is the better claim on the time. Readiness is not a
reason.

Answer before Phase 1. Cheap here, expensive later.

---

## 7. Phase 1 — Arena hygiene

**Character:** Mechanical, low-risk, high-certainty. Do this regardless of what
the gap table says.

### Arena assignment

| Data class | Arena | Rationale |
|---|---|---|
| MultiFab field data | `The_Arena()`, device, non-managed | Governing rule |
| Kernel temporaries | **`The_Async_Arena()` — evidence-gated, see below** | Stream-ordered free, no sync on destruction |
| I/O staging, reduction landing | `The_Pinned_Arena()` | DMA-direct; pageable D2H forces a synchronous driver staging copy |
| MPI buffers (non-GPU-aware path) | `The_Comms_Arena()` | Isolation from field pool |
| Managed | Nothing (NOVA); required locally on A1000 | Unless individually justified in writing |

**Async arena caveat (C3).** `The_Async_Arena()` was already A/B'd on this repo
and **rejected**: +2.2% 2D, +4.0% 3D. Do not re-adopt it as a default on the
strength of the general argument. Re-test only after Phase 2 removes the
per-step drain — the earlier measurement was taken against a sync-bound
baseline, so the result may not survive, but that is a hypothesis to test, not a
reason to assume.

### Configuration

- `amrex.the_arena_is_managed=0` — permanent on NOVA, not a debug toggle. Local
  A1000 stays managed (C4).
- `amrex.the_arena_init_size` — claim the pool up front to kill fragmentation
  and mid-run allocation stalls. **Express as a fraction of device memory
  divided by ranks-per-GPU, not a hardcoded byte count.** Currently pinned by
  hand locally because the default (3/4 of 8 GB) intermittently OOMs at init —
  that hand-pin is exactly the hardcode this policy replaces.
- `amrex.abort_on_out_of_gpu_memory=1` — guardrail

Parameter names and defaults have moved across AMReX versions. Verify against
the **vendored** 26.06 tree.

### Method — flip per subsystem

Flip the arena after each subsystem, not once at the end. Porting everything and
then flipping produces fifty simultaneous segfaults with no attribution and
unknown coupling.

Managed memory during porting is a scaffold, not a strategy: it keeps unported
code correct-but-slow so the regression harness stays green throughout. Its cost
is that a forgotten host touch is a silent 10-100× penalty with no signal.
Flipping converts an invisible performance bug into a loud correctness bug with
a line number. That conversion is the entire value.

Note the local wrinkle: because the A1000 must stay managed, the flip is only
observable on NOVA. Budget queue turnaround into the per-subsystem cadence, or
the "loud correctness bug with a line number" arrives a day late.

### Triage of resulting crashes

Not every crash is a missing kernel port. Classify before fixing:

| Crash site | Fix |
|---|---|
| Inside the timestep loop | Port the kernel |
| Diagnostic / screen output | Explicit staged copy through pinned |
| Checkpoint / plotfile | Explicit staged copy, async |
| Debug leftovers | Delete |

Common hiding places: `mf[mfi](i,j,k)` in diagnostic code, hand-rolled host
reduction loops, unported BC fill, `MFIter` loops with plain C++ bodies. See
`docs/llm/BUG_PATTERNS.md` for the device bug classes this branch has already
paid for.

### Exit gate

T1, T4, T5 pass on NOVA. Regression harness green. Device lint clean.

---

## 8. Footprint budget

**Do this arithmetic during Phase 0, before committing to Phase 1's `init_size`
policy.**

```
resident = cells * components * (1 + ghost_overhead) * sizeof(Real) * n_multifab_copies
```

Sum across levels and across the solver's working copies. Compare against 40 GB
(A100) or 64 GB (MI210), divided by ranks-per-GPU. **Also run it against 8 GB**
— that number decides which decks remain runnable locally at all, which sets how
much correctness work can happen off-queue.

Two consequences:

- If the budget clears device memory, FP32 storage stops being an optimization
  and becomes a structural requirement — and that changes what Phases 1-3 build.
  Better known now than discovered as an OOM mid-Phase 1.
- On a bandwidth-bound code, storage precision is a *memory strategy* decision,
  not an arithmetic one. FP32 halves resident footprint and halves traffic: the
  same lever pulled twice.

Chamber is likely comfortable. Run the number anyway — the same calculation is a
hard input for `multicomponent/FMA`, where species count multiplies field count.

---

## 9. Phase 2 — Sync elimination

**Character:** Design work. The chamber-specific phase, and probably the largest
win.

### The chamber feedback loop

```
reduce regression rate over burning surface
  -> mass flux
  -> scalar ODE for chamber pressure
  -> pressure feeds Arrhenius mobility and burn rate law
  -> next step
```

The naive implementation is `reduce -> D2H -> host ODE -> H2D -> next kernel`.
That is a full pipeline drain every timestep, and it does not appear as a slow
kernel. It appears as dead timeline — which is why the NVTX-annotated timeline
from Phase 0 is the diagnostic that finds it.

### Confirmed sync sites (source read, 2026-07-27)

Concrete Phase 2 inventory, already located:

| Site | What it costs |
|---|---|
| `Flame.cpp:1127` `reduce_data.value(reduce_op)` | Device sync **per box, per level**, inside the `MFIter` loop — not once per step. Lands the reduction on the host to accumulate into `chamber.*`. Prime T3/T6 target. |
| `Flame.cpp:664-701` | A second `ReduceOps` block plus **five** `ParallelDescriptor::ReduceReal{Max,Min}` calls (`thermo_max_temp`, `thermo_mdot_max`, `thermo_heatflux_max`, `thermo_L_max`, `thermo_eta_min`) — five separate allreduces where one batched call would do. |
| `Integrator.cpp:1271-1275` | Per-variable `ReduceRealSum` loop over `thermo.extensives` — batchable. |
| `Solver/Nonlocal/Newton.H` ~419, ~572, `FieldNorm0` ~816 | `norm0` per level × per component, each a device reduction + stream sync + MPI allreduce, multiplied by line-search backtracks. **Already tracked as PLAN.md backlog item 3.I** — adopt it into this phase rather than duplicating it. |

### Target shape

Chamber pressure lives permanently in a one-element device array. Reduction
lands on device. The scalar update is a single-thread kernel, or folded into the
head of the following kernel. Host observes pressure only when diagnostics print.

A one-thread kernel is a strange-looking object. It costs roughly 5 µs. It
replaces a full drain.

### Build the reduction Allreduce-shaped now

The one multi-GPU decision that cannot be cleanly deferred.

The burning-surface reduction is global (confirmed, §5). Multi-rank, the scalar
update is an `MPI_Allreduce`, not a local reduction. Building Phase 2 as
`device-local reduce -> one-thread kernel` produces a design that is *replaced*,
not extended, at Phase 4.

Build instead: `device-local reduce -> device-buffer Allreduce -> device scalar
update`. On one rank the Allreduce is a no-op. Cost single-GPU is nil;
correctness multi-GPU is by construction.

Caveat: this is only a device-buffer Allreduce if GPU-aware MPI is active. Phase
0.5 Q3 says it currently is not. Either enable it or accept a staged path and
mark it as debt — do not write "device-buffer Allreduce" in the source and let
it silently stage through host.

### Newton solver convergence check

If the Newton convergence test is a host-side comparison on a reduced residual,
that is an additional sync *per Newton iteration*, stacked on the per-step drain.

Options: fixed iteration count with a check every N iterations, or batch
residual norms and test on device.

Confirm whether the CPU Newton damping / line-search work left a host-side
branch in place. A host-side convergence branch is a reasonable CPU design and a
serious GPU defect; exactly the kind of thing that survives a port unexamined.
**Convergence-semantics-critical: tier 3, CPU golden compare + budget gate.**

### Prior art to fold in, not rediscover

`ALAMO_MLMG_NOSYNC` (default-off: `Gpu::NoSyncRegion` + `max_gpu_streams=1`)
already measures **-12.1% 2D / -1.5% 3D** with bit-identical traces, and is
retained but unshipped. It is a Phase 2 deliverable that already exists —
validate and promote it rather than re-deriving it. Its 2D/3D asymmetry is
itself the T6 story.

### Exit gate

T2, T3, T6 pass. Per-step bus traffic is scalars only. Timeline shows no drain
inside the step loop.

---

## 10. Phase 3 — Launch configuration and kernel shape

**Character:** Tuning. Deliberately last — the bottleneck relocates after Phase
2, so anything tuned earlier is tuned against a defect.

### Levers

**Box size.** `max_grid_size` up from CPU-tuned values (32³ is launch-bound on
GPU); 128 or higher is the usual landing zone. Confirm `TilingIfNotGPU()`
throughout — CPU tiling on GPU is pure overhead. **Known lead:**
`Integrator.cpp:1246` and `:1263` construct `MFIter mfi(grids[...], dmap[...],
true)` with tiling hardcoded true, in the integrate path. Check whether that is
`TilingIfNotGPU()`-guarded; if not, it is a defect, not a tuning knob.

*Constraint:* box count is the domain decomposition granularity. Tuning to 256
on a 512³ domain yields 8 boxes — adequate on 8 ranks, poorly balanced on 6,
unusable on 12. Choose box size against a target rank count, or accept a retune
at Phase 4.

*Prior evidence:* wide-shallow AMR beats deep AMR on GPU here.

**Fusion.** Fuse when it removes a full field read/write cycle. Split when
registers spill. Verify with `-Xptxas -v`, or `benchmark/res_usage.sh` as the
root-free `ncu` substitute. A large fused multiphysics kernel is elegant and
frequently occupies a fraction of the SM.

*Prior evidence:* Fapply is 252 regs in 3D vs 94 in 2D. **2D cannot show
register wins** — run register experiments in 3D or they measure nothing.
Non-RDC 252→248 regs produced no occupancy change and was rejected.

**Occupancy limiter.** Identify what actually caps occupancy — registers, shared
memory, or block count — before adjusting anything. Nsight Compute reports this
directly. Note that on this code the historical Fapply wins came from
spill/replay reduction at **flat occupancy** (~12%), so "raise occupancy" is not
automatically the objective.

### Framing

Phase-field and stencil physics are memory-bound. Consequences:

- Optimize bytes moved, not floating-point operations
- FP32 gains come from bandwidth, not ALU throughput
- Fusion pays because it removes traffic
- Hardware selection is an HBM bandwidth question, not a TFLOP question

### Exit gate

T7 pass, against the threshold set at Phase 0.

---

## 11. Phase 4 — Multi-GPU (deferred)

Deferred content:

- Halo exchange overlap and communication hiding
- Load balance and regrid distribution strategy
- Scaling studies, strong and weak

All three require an optimized single-GPU baseline to be measurable.

Carried forward from earlier phases:

- Global Allreduce-shaped chamber reduction (built at Phase 2)
- Rank-aware `init_size` policy (set at Phase 1)
- Box size chosen against target rank count (Phase 3)
- GPU-aware MPI confirmed active (verified at Phase 0.5 — currently believed inactive)

Known starting point: multi-GPU is presently a **loss** (5.3% efficiency), root
-caused to blocking comm + no GPU-aware MPI + `regrid_int=2` + managed arena.
Three of those four are addressed by Phases 1-2, which is the argument for
deferring Phase 4 rather than attacking it now.

---

## 12. Instrumentation and man-in-the-loop tooling

A human looks at a fixed figure set at each phase boundary and can see what
moved. Regenerating figures must be one command, or it will not happen
consistently.

### Prerequisite — NVTX annotation

Annotate the step loop by physics stage. Without it the timeline is unreadable
and the highest-value diagnostic in this plan is unavailable. Verify the AMReX
NVTX build flag against the **vendored** 26.06 tree.

### Standing capture

| Tool | Purpose | Scope |
|---|---|---|
| Nsight Systems | Timeline, transfers, syncs, page faults, idle gaps | Whole run, few steps |
| Nsight Compute | Per-kernel roofline, occupancy, achieved bandwidth, registers | Top kernels by time only |

Nsight Compute serializes and instruments heavily — restrict it to the top few
kernels and a handful of invocations. Locally `ncu` counters may be barred
(`ERR_NVGPUCTRPERM`); `benchmark/res_usage.sh` is the root-free substitute. On
NOVA, confirm counter permission before budgeting a job around `ncu`.

Trace domains: CUDA, NVTX, MPI, OS runtime, unified-memory page-fault counters,
CUDA memory usage. Pin the tool version in the harness and verify flag spellings
once.

### Figure set

Regenerate all of these at every phase boundary. Consistency across phases is
what makes them useful.

| # | Figure | Reads on | Phase relevance |
|---|---|---|---|
| F1 | Annotated timeline, one full timestep | T6 | 0, 2 — primary sync diagnostic |
| F2 | Memcpy summary, bytes and count by direction per step | T3 | 0, 2 |
| F3 | Unified-memory page-fault count | T2 | 0, 1 |
| F4 | Kernel time Pareto, top 10 with cumulative % | Targeting | 0, 3 |
| F5 | Roofline, top kernels vs HBM roof | T7 | 0, 3 |
| F6 | Achieved DRAM bandwidth as % of peak, per kernel | T7 | 3 |
| F7 | Occupancy with limiting resource identified | T7 | 3 |
| F8 | Arena high-water mark over steps | T5 | 1 |
| F9 | GPU idle fraction inside the step loop | T6 | 2 |
| F10 | Wall time per step, phase over phase | Headline | All |

F1 is the one to read first at every phase boundary. Human pattern recognition
on a timeline outperforms any aggregate metric at answering "why is there a gap
here."

F5 and F6 are the pair that matter at Phase 3: for a bandwidth-bound code,
distance from the HBM roof is the target, and compute-side metrics are largely
decorative.

**Capture every figure in both 2D and 3D.** The two are bound by different
things (2D launch-latency, 3D Fapply-compute) and a single-dimension figure set
will mislead.

### Metrics ledger

One CSV, appended at every phase boundary. Keyed by `(phase, commit_sha,
case_name, n_ranks, device)`. Columns: T1-T7 measurements plus wall time per step.

Purpose is regression detection. Phase 3 tuning can quietly undo a Phase 1
property — a fusion change that reintroduces a temporary allocation in the step
loop — and without a ledger this is invisible until much later.

Reuse the existing path-keyed ledger machinery rather than building new.

### Automation

Target: one command produces the full figure set plus the ledger row.

```
make profile CASE=chamber_small PHASE=1
```

Given the SLURM constraint (§2), the realistic shape is *one command composes and
submits a job; a second collects and renders on return*. Hand-running the
profiler at each boundary is the failure mode where the figure set gets captured
twice and then abandoned. Build the automation during Phase 0, when the harness
is being set up anyway.

---

## 13. Pattern capture

Chamber optimization is also template generation for `multicomponent/FMA`.
Capture as work proceeds rather than reconstructing afterward:

| Pattern | Destination |
|---|---|
| Arena assignment by data class | `gpu_manual` |
| Async arena discipline for temporaries (incl. the rejection result) | `gpu_manual` |
| Device-resident scalar with Allreduce-shaped reduction | `gpu_manual` |
| Sync inventory and elimination checklist | `gpu_manual` |
| Flip-per-subsystem migration method | `gpu_manual` |
| Gap table and figure set | `CHAMBER_GPU_MASTER_REPORT.md`, new section |

The multicomponent port should then be a mechanical pass against a pattern set
rather than a fresh design exercise. Genuinely new work there is confined to the
advection kernels and to whatever register pressure the species arrays introduce.

---

## 14. Open questions — status

| # | Question | Status |
|---|---|---|
| 1 | Does chamber pressure update host-side, and where does the reduction land? | **Answered.** Host-side. Per-box device reduce at `Flame.cpp:1113-1125`, `.value()` sync at `:1127`, host accumulate at `:1128-1130`, MPI sum at `Integrator.cpp:1271-1275`, ODE advance at `Flame.cpp:706`. |
| 2 | Is the burning-surface reduction global or rank-local? | **Answered: global** — `ReduceRealSum` runs before `TimeStepComplete`. One residual check: confirm the `extensives` flag is set for volume/area/mass_flux. |
| 3 | Single rank or multi-rank in current use? | Single-rank is the working configuration; multi-GPU is a measured loss (5.3% eff). Multi-rank correctness still gets probed at 0.5. |
| 4 | Available hardware? | **Answered.** Local A1000 8 GB (shared, 50 W-capped, correctness only) + NOVA A100/H200 via SLURM over `ssh -MN` tunnel. Phase 4 is testable on NOVA multi-GPU nodes; scope depends on allocation, not procurement. |
| 5 | Does the regression harness cover chamber? | **Answered: yes** — `baseline_suite.py`, 4 cases, gated by `ci_golden_compare.sh`, green on this branch. Gap: stale `rod_and_tube_step2` gpu_strict reference. |

Q1-Q3 were the ones collapsing the Phase 2 scope estimate. With them answered,
Phase 2's target list is the concrete site inventory in §9, not an open search.

---

## 15. Risk register

| Risk | Phase | Response |
|---|---|---|
| `extensives` flag unset → reduction silently rank-local after all | 0.5 | One-line source check; cheapest item in the plan |
| Footprint exceeds device memory under device arena | 0-1 | Compute budget at Phase 0; FP32 storage becomes structural |
| GPU-aware MPI silently inactive | 0.5 | Believed inactive; verify at runtime, otherwise baseline is polluted and the Phase 2 "device-buffer Allreduce" is a fiction |
| Newton host-side convergence branch survives from CPU design | 2 | Audit explicitly; do not assume the port caught it. Tier 3. |
| Phase 3 tuning silently undoes a Phase 1 property | 3 | Metrics ledger; re-check T4/T5 at Phase 3 exit |
| Box size retune required at Phase 4 | 3 | Choose against target rank count now |
| Chamber is not the throughput bottleneck | 0 | Decision gate before Phase 1 |
| **SLURM turnaround throttles the per-subsystem flip cadence** | 1 | Batch subsystem flips; accept coarser attribution rather than serializing on the queue |
| **Local A1000 cannot validate T1-T7** | all | Local runs are correctness-only by policy; no timing claim ships from kermit |
| **Concurrent gate runs produce spurious RED** | all | One gate at a time; re-run solo before believing a failure |
| Glaring-defect exception becomes scope creep | all | Evidence in `NOTES.md` first, own tier-3 task folder second, no inline edits |

---

## 16. Ship rule

CLAUDE.md hard rule, restated because this campaign is entirely perf work:

> No kernel/perf optimization ships without a correctness pass — device lint,
> golden compare, compute-sanitizer on A1000.

Tier-2 gate ignores mid-run aborts; supplement with a converging-deck memcheck.

---

## Note on specifics

Configuration parameter names, AMReX NVTX integration, and Nsight command-line
flags shift between versions. Treat every such name here as a pointer to verify
against the **vendored** AMReX (`ext/AMReX-Codes/amrex`, 26.06) and the
installed toolkit, not as a literal. Strategy, phase ordering, and target
definitions are version-independent.

**Confidence:** ~90% on strategy, phase ordering, and target set. ~70% on exact
parameter and flag spellings. Source pre-answers in §5, §9, and §14 are
file:line-cited reads from 2026-07-27 and carry higher confidence than the
draft's assumptions, but were not empirically probed — Phase 0.5 still runs.

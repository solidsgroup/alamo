# TASK: phase05-two-rank-probe
# Folder: docs/agent_plans/20260727-phase05-two-rank-probe/

Phase 0.5 of the `chamber-gpu-mem` memory-strategy campaign
(`docs/agent_plans/20260727-gpu-memory-strategy/PLAN.md` §5). Runs before
Phase 0 baseline measurement.

---

## Header

| Field        | Value                                                                 |
|--------------|-----------------------------------------------------------------------|
| Risk tier    | 1 — benchmark/probe scripts + docs; no `src/` change authorized here  |
| Model        | opus (campaign orchestration); probe execution is mechanical           |
| Verification | partial-oracle — `benchmark/compare_thermo.py` + `benchmark/status.sh` |
| Est. scope   | this folder + one probe script under `benchmark/`; ~150 lines          |
| Parallel-safe| no — probe runs use `bin/alamo*` and would race `status.sh` (§2 warning)|

## Operating rules

1. Read ONLY the files listed in Context budget. `docs/archive/` is forbidden.
2. Every step's VERIFY must pass before acting on that step. On failure: STOP,
   report discrepancy, wait.
3. One commit per step. Message: `<area>: <what> (20260727-phase05-two-rank-probe)`.
4. No scope expansion. New ideas go to the campaign `NOTES.md`, not into code.
   A fix discovered here gets its own tier-3 folder (campaign PLAN §1).
5. Missing knowledge → ask.
6. Tier 1: unattended, gates + spot-check diff.

## Context budget

Read first: `benchmark/status.sh` output, this PLAN.md, campaign PLAN §5
Read: `src/Integrator/Flame.cpp:190-240`, `:655-720`, `:1100-1140`,
`src/Integrator/Integrator.cpp:115-130`, `:1180-1290`,
`src/Integrator/Integrator.H:240-250`, `:440-460`,
`benchmark/compare_thermo.py`, `benchmark/baseline_suite.py:20-70`, `input`
Reference only if a step names it: `benchmark/NOVA_SLURM_RUNBOOK.md`,
`ext/AMReX-Codes/amrex/Docs/**/GPU.rst`
Forbidden: `docs/archive/*`, unrelated task folders, the ~/Desktop sweep campaign

## Objective

The chamber model reduces regression rate over the burning surface into a scalar
pressure ODE. If that reduction were rank-local, every multi-rank chamber result
would be physically wrong, and the Phase 2 device-resident-scalar redesign would
be built on a broken contract. A source read (campaign PLAN §5) says the
reduction is structurally global; this task converts that read into runtime
evidence, and answers the two remaining Phase 0.5 questions (domain-decomposition
survival, GPU-aware MPI status) before any measurement or optimization begins.

After this task: a recorded, reproducible verdict on whether 2-rank chamber
pressure history matches 1-rank to solver tolerance, on both CPU and GPU builds,
plus a runtime record of GPU-aware MPI availability on kermit (and the command to
re-check it on NOVA).

## Oracle

Command(s):
- `benchmark/compare_thermo.py <1rank>/thermo.dat <2rank>/thermo.dat --rel-tol 1e-6`
  exit 0 = pass
- `benchmark/status.sh` — three gates green, run **solo**

Covers: chamber scalar history (pressure, volume, area, mass_flux) agreement
between 1 and 2 ranks over the probe horizon; that no rank-boundary ghost/BC
defect perturbs the integrated quantities; existing CPU golden + device lint +
sanitizer gates still green.

Does NOT cover: per-cell field agreement away from the integrated quantities
(thermo.dat is a reduction, so a sign-cancelling field error could hide);
long-horizon divergence past the probe step count; multi-*node* MPI; NOVA
GPU-aware MPI (login-node/compute-node check only, recorded not gated);
performance of any kind (§2 — kermit is correctness-only).

## Steps

### Step 1 — Source confirmations (Q1 residual, staleness check)

VERIFY:
```bash
grep -n "RegisterIntegratedVariable" src/Integrator/Integrator.H
grep -n "RegisterIntegratedVariable(&value.chamber" src/Integrator/Flame.cpp
grep -n "thermo.interval" src/Integrator/Integrator.cpp
```
DO: record in `results/RESULT.md` whether the `extensive` flag is set for
`volume`, `area`, `mass_flux`, and whether `IntegrateVariables`' gate
(`Integrator.cpp:1229`) fires every step for the campaign decks. Both feed
campaign PLAN §14 Q2 and risk-register row 1.
CHECK: `results/RESULT.md` §1 populated with file:line citations.

**Status: DONE 2026-07-27 — see results/RESULT.md §1.**

### Step 2 — GPU-aware MPI runtime probe

VERIFY:
```bash
ompi_info --parsable --all | grep mpi_built_with_cuda_support
```
DO: record the local verdict. Add the same probe as a one-liner to be run on a
NOVA compute node (not the login node — the MPI module differs), and record the
command in `results/RESULT.md` for the Phase 0 NOVA batch.
CHECK: `results/RESULT.md` §2 states local verdict and the NOVA re-check command.

**Status: DONE 2026-07-27 (local leg) — see results/RESULT.md §2.**

### Step 3 — Probe harness

VERIFY:
```bash
test -x bin/alamo-2d-g++ && test -x bin/alamo_gpu-2d-cuda86-g++
git status --porcelain benchmark/ | head
```
DO: add `benchmark/two_rank_probe.sh` — runs the `input` deck for a fixed
`max_step` at `-np 1` and `-np 2`, into per-run output dirs under the scratch
tree (never `benchmark/baseline_runs/`, per the §2 concurrency warning), then
invokes `compare_thermo.py`. Parameterized by build (cpu | gpu) and step count.
CHECK: `bash -n benchmark/two_rank_probe.sh` clean; `--help`/dry-run prints the
exact `mpiexec` lines without executing them.

**Status: DONE 2026-07-27 — see results/RESULT.md §3. Deviation: `MGS` default 32
added after the first run showed the deck's default grid leaves level 0 as one
box on rank 0.**

### Step 4 — CPU two-rank probe

VERIFY: Step 3 CHECK passed; no other gate or benchmark run is in flight.
DO: `benchmark/two_rank_probe.sh cpu <max_step>` on the `input` deck.
CHECK: `compare_thermo.py` exit 0 at `--rel-tol 1e-6`. If it fails, capture the
first diverging column and step in `results/RESULT.md` and STOP — do not tune the
tolerance to make it pass.

**Status: DONE 2026-07-27 — PASS, bit-identical. Legs A (np=2), B (np=4), D
(np=2, `elastic.interval=1`). See results/RESULT.md §4.**

### Step 5 — GPU two-rank probe

VERIFY: Step 4 passed. Both ranks share the single A1000; confirm the run does
not OOM at init (pin `amrex.the_arena_init_size` per campaign PLAN §7).
DO: `benchmark/two_rank_probe.sh gpu <max_step>`.
CHECK: `compare_thermo.py` exit 0 at the same tolerance. Record the tolerance
actually needed; GPU 1-vs-2-rank will not be bit-identical (reduction order), and
the required tolerance is itself a Phase 2 input.

**Status: DONE 2026-07-27 — PASS, bit-identical at `rel_tol=1e-6`; no loosening
needed. Legs C (10 steps) and E (`elastic.interval=1`). See results/RESULT.md §4.**

### Step 6 — Verdict and closeout

DO: write `results/RESULT.md` — verdict per Phase 0.5 exit gate, the tolerance
each leg needed, and any defect found (defects route to their own tier-3 folder,
not fixed here). Append the campaign `NOTES.md` if anything new surfaced.
CHECK: `benchmark/status.sh` green, run solo.

**Status: DONE 2026-07-27 — verdict PASS, see results/RESULT.md §5. Campaign
NOTES.md gained N6 (chamber ODE freshness coupled to `amr.thermo.int`).**

## Checkpoints

Tier 1 — unattended. Report to the user at Step 4/Step 5 failure, or at Step 6.

## Closeout

- [ ] Oracle passes; `status.sh` all green (solo run)
- [ ] `results/RESULT.md`: verdict, evidence, deviations from plan
- [ ] `changelog/` entry (append-only) if a probe script is retained
- [ ] `touch results/DONE`
- [ ] Session log line appended → `docs/llm/SESSION_LOG.tsv`

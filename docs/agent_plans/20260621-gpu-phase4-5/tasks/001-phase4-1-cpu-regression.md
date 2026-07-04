# Task 001: Phase 4.1 — CPU-semantics regression across all integrators

## Goal
Confirm the GPU port's de-virtualization of the `Solid`/`BC`/`Operator`
hierarchies (removing `virtual`) did NOT change any CPU integrator's results.
Run the CPU regression suite on the chamber-gpu branch state and report pass/fail.

## Context
chamber-gpu removed `virtual` from shared model hierarchies for nvcc. Integrators
that dispatch through a runtime `Solid*` are the risk surface: Eshelby*, PlateHole,
Rubber*, SCP*, Inclusion, CrystalPlasticity, DynamicBar, Fracture*, Solid,
CompositeImpact. Phase-field/flow integrators (AllenCahn, CahnHilliard, Dendrite,
Flow*, Heat*, PFC) do not use Solid dispatch and are a secondary check.
`scripts/runtests.py` compares each test against its stored reference, so a PASS
means semantics were preserved.

## Branch guard (DO FIRST — you run in an isolated worktree)
Run `git cat-file -e HEAD:src/alamo_gpu.cc && echo OK_CHAMBER_GPU` and
`git log --oneline -5`. The worktree MUST be based on chamber-gpu (alamo_gpu.cc
exists at HEAD, GPU commits like `fffbec0a2` present). If alamo_gpu.cc is absent
at HEAD, you are on the wrong branch — STOP and report that in RESULT.

## Files allowed to modify
- Only create logs under `benchmark/phase4_cpu_semantics_<timestamp>/`. NO source/Makefile/configure edits.

## Files NOT allowed to modify
- Any `src/**`, `Makefile`, `configure`, `scripts/**`, `tests/**`.

## Implementation steps
1. Branch guard (above).
2. 2D first (minimum bar): `./configure --dim=2 --comp=g++ && make -j8` then
   `scripts/runtests.py --dim=2 --comp=g++ --permissive --no-backspace --timeout=2000`
   (if g++ is unavailable, fall back to `--comp=clang++` and note it). Tee output to a log.
3. Summarize 2D: total tests, passed, failed; list every failure with its name.
4. 3D (best-effort): same with `--dim=3`. If it cannot finish in your time budget,
   record 2D-complete + 3D-partial and list which 3D tests ran.
5. Interpret: for each FAILURE, state whether it is a Solid-dispatch integrator
   (de-virtualization suspect) or a phase-field/flow one (more likely a pre-existing
   chamber physics diff). A clean Solid-integrator pass set = de-virtualization safe.

## Invariants
- Read-only on the branch's source; this is a measurement task.
- Do not "fix" any failing test — report it.

## Build and test commands
See steps 2/4.

## Expected result
Pass/fail counts per dim, the failure list with Solid-vs-phasefield classification,
and a verdict: did de-virtualization preserve CPU semantics? (yes / no / inconclusive-needs-3D).

## Non-goals
Fixing tests; GPU builds; touching CI.

## Stop conditions
Wrong-branch worktree; build failure (report the compiler error — itself a finding).

## Final report: write results/001-RESULT.md
Use the RESULT.md format. Put the pass/fail tables and verdict in Summary.

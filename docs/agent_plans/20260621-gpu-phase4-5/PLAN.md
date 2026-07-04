# Plan: Execute Phase 4 + Phase 5 of GPU-OPT-ROADMAP on chamber-gpu

## User goal
Autonomously execute Phase 4 (framework dispatch decision/generalization) and
Phase 5 (hardening) of `~/Desktop/GPU-OPT-ROADMAP.txt`. Dispatch agents to do it.

## Hard constraint
`chamber-gpu` is **NEVER merged to master** (user directive 2026-06-21). The
branch is the deliverable. This forces the D4 dispatch decision to **ISOLATE**
(contain the device-model in the CUDA build; never touch the shared CPU build for
other integrators). No agent commits to `chamber-gpu` or merges anything; new
files land in the working tree (or a worktree) for the user to review.

## Current architecture (verified this session)
- GPU work lives behind the CUDA build only: `src/alamo_gpu.cc` (Flame-only main)
  + `src/GPU/IntegratorPolicy.mk` (per-integrator GPU-clean source closure) wired
  via Makefile `ifneq (,$(findstring cuda,$(POSTFIX)))`. The CPU/all-integrator
  launcher `alamo.cc` is unchanged.
- **4.3 (de-fork) is ALREADY DONE**: `IntegratorPolicy.mk` replaced the hand-
  curated source list with a principled, build-tracked policy
  (`ALAMO_GPU_SUPPORTED_INTEGRATORS := flame`). Remaining 4.3 work = verify +
  document, not implement.
- Tests: `scripts/runtests.py` over `tests/` (53 integrators), dims 2 & 3.
  `benchmark/phase4_cpu_semantics_regression.sh` wraps it (4.1) — ready, unrun.
- CI: `.github/workflows/` (linux.yml, performance.yml on self-hosted `scooter`).
  None triggers on `chamber-gpu`.
- Golden refs: `benchmark/baseline_references/{canonical_step1,canonical_step2,
  eta_expression_step1}`; harness `benchmark/baseline_suite.py`.
- Prior decisions: D1 = elastic CPU-resident; D3 = WIN @ single (A100 ~39-70× vs
  CPU node). See `benchmark/archive/PHASE1_ELASTIC_DISPOSITION.md`, `PHASE3_R3_crossover.md`.

## Desired architecture
Phase 4/5 closed out for a permanent branch: 4.1 regression data captured; R4
records D4=ISOLATE (pending Runnels ratification) + documents the done de-fork;
encapsulation convention documented; Phase-5 hardening artifacts (correctness CI,
perf-regression tracking, consolidated docs, branch definition-of-done) created.

## Invariants and constraints
- Never merge/commit to `chamber-gpu`; never edit master. Create NEW files; only
  task 003 extends one solely-owned existing doc.
- Do not edit the uncommitted working files: `Makefile`, `configure`,
  `src/Integrator/Flame.H`, `benchmark/archive/PHASE3_R3_crossover.md`.
- No two tasks write the same file (see ownership below).
- D4 = ISOLATE is fixed by the no-merge policy; agents document it, they do not
  re-litigate it, and they flag it needs human (Runnels) ratification.

## Files involved / ownership (disjoint)
- 001: creates `benchmark/phase4_cpu_semantics_<ts>/` logs only (no source edits).
- 002: `benchmark/archive/PHASE4_R4_dispatch.md` (new).
- 003: `docs/gpu_device_capture_conventions.md` (extend; sole owner).
- 004: `.github/workflows/chamber-gpu-correctness.yml` + `benchmark/ci_golden_compare.sh` (new).
- 005: `benchmark/perf_regression_track.py` + `benchmark/PERF_TRACKING.md` (new).
- 006: `benchmark/GPU_BRANCH_GUIDE.md` + `benchmark/archive/PHASE5_BRANCH_DONE.md` (new).

## Build and test commands
- CPU build: `./configure --dim=<2|3> --comp=g++ && make -j8`
- Tests: `scripts/runtests.py --dim=<d> --comp=g++ --permissive --no-backspace --timeout=2000`
- Golden: `CPU_NP=8 python3 benchmark/baseline_suite.py check`

## Task graph
- 001 (serial, heavy, worktree) → produces regression data.
- 002 (serial, after 001) → R4, consumes 001 summary.
- 003, 004, 005, 006 (parallel-safe, main tree) → independent artifacts.

## Parallelization strategy
Dispatch 001 + 003 + 004 + 005 + 006 concurrently (001 is worktree-isolated; the
others create disjoint new files in the main tree and do not build). Dispatch 002
after 001 returns, injecting 001's RESULT summary.

## Integration strategy
Lead reads each RESULT, reviews diffs, writes INTEGRATION.md. No auto-commit; the
artifacts sit in the working tree for the user. 001's worktree only holds test
logs — summary returns via RESULT; worktree is discarded.

## Risks
- Worktree base: 001 must be based on chamber-gpu, not master. Guard: verify
  `git cat-file -e HEAD:src/alamo_gpu.cc` (exists only on chamber-gpu); STOP if not.
- 4.1 3D suite is slow: 2D-complete is the minimum bar; 3D best-effort.
- Test failures may be pre-existing chamber physics diffs, not de-virtualization —
  001 must distinguish (focus on Solid-dispatch integrators).
- GPU CI needs a CUDA runner that may not exist — 004 ships the CPU-side gate and
  documents the GPU leg as runner-gated.

## Rollback plan
All outputs are new untracked files; `rm` them to roll back. No source/build/git
state is modified.

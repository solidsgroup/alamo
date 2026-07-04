# Result 006: Phase 5.3 + 5.4 — doc consolidation + branch definition-of-done

## Status: DONE

## Files created
- `benchmark/GPU_BRANCH_GUIDE.md` (new)
- `benchmark/archive/PHASE5_BRANCH_DONE.md` (new)

No other files were read-modified. Source docs (`docs/gpu_safe_ic_bc_matrix.md`,
`benchmark/README.md`, `benchmark/archive/PHASE1_ELASTIC_DISPOSITION.md`,
`benchmark/archive/PHASE3_R3_crossover.md`, `docs/gpu_device_capture_conventions.md`,
`benchmark/archive/G0_BASELINE_OF_RECORD.md`) and all source/CI files were read-only.

## Summary

### `benchmark/GPU_BRANCH_GUIDE.md`
A navigable index (not a copy) covering: branch-never-merged policy; build
matrix table (CPU / GPU fast / GPU strict-no-fast-math / NOVA 3D, with the
exact build command and binary name for each, pointing to
`benchmark/README.md` for full detail); the GPU-safe IC/BC matrix summary
(linking `docs/gpu_safe_ic_bc_matrix.md`); the elastic disposition D1 =
CPU-resident summary (linking `benchmark/archive/PHASE1_ELASTIC_DISPOSITION.md`); the
crossover D3 = WIN @ single summary (linking `benchmark/archive/PHASE3_R3_crossover.md`);
a phase-report location table (Phases 0-5); and a correctness-tooling quick
reference. Also notes the Phase 4.3 de-fork mechanism
(`src/GPU/IntegratorPolicy.mk`) and points to `benchmark/archive/PHASE4_R4_dispatch.md`
for the D4 = ISOLATE decision, since the branch-isolation policy is load-bearing
context for "why this guide exists."

### `benchmark/archive/PHASE5_BRANCH_DONE.md`
The roadmap's 5-item branch definition-of-done as a checklist, each item with a
DONE/PARTIAL/PENDING status and evidence pointer:

1. **Golden compare (no-fast-math), coarse + saturating — PARTIAL.** Coarse
   passes (`G0_BASELINE_OF_RECORD.md`); no saturating-3D-config golden compare
   has been run yet (`PHASE3_R3_crossover.md` itself lists this metric as
   pending in its standing metric set).
2. **No reachable host-loop device-arena writes — DONE.** Per
   `docs/gpu_safe_ic_bc_matrix.md`, the supported IC/BC surface (BMP/Constant
   IC, Constant elastic BC) has zero host-loop writes; unsupported paths are
   guarded to abort first (verified by `test_gpu_guarded_ic.sh`).
3. **Device aborts/NaN detection active — DONE.** G0 corrected-build status +
   the guarded-IC abort behavior are direct evidence; CI-enforced form is
   pending in task 004's `chamber-gpu-correctness.yml`.
4. **CPU regression suite green, all integrators — PENDING.** This depends on
   sibling task 001's result, which had not landed in
   `docs/agent_plans/20260621-gpu-phase4-5/results/` at the time this task ran
   (results/ only contained this task's own output). Pointed at
   `results/001-RESULT.md` and `benchmark/archive/PHASE4_R4_dispatch.md` (task 002) as
   the evidence to check once available. Explicitly NOT marked DONE on
   circumstantial build-isolation evidence alone.
5. **R3 win recorded — DONE.** `PHASE3_R3_crossover.md` records D3 = WIN @
   single (39-70x at 128-256 cubed on NOVA A100), with "WIN @ scale" honestly
   left unclaimed (multi-GPU currently regresses at small per-GPU domains).

Overall verdict recorded: 3/5 DONE, 1 PARTIAL, 1 PENDING — branch is not yet at
its own done-bar, reaffirmed as "branch done" criteria only, never a merge
gate. The doc explicitly tells the reader to re-check (not rewrite) items 1 and
4 once tasks 001/004/005 land.

## Notes for the lead
- At the time this task ran, `results/` contained no other RESULT files from
  sibling tasks 001-005, so item 4 (CPU regression) and the CI/perf-tracking
  pointers in item 1/3 are necessarily forward-pointers to expected artifacts,
  not verified-complete citations. This matches the task's own instruction
  ("mark status honestly... point to the expected artifact even if still
  landing").
- No commits were made (per instructions). Both files are new and untracked.

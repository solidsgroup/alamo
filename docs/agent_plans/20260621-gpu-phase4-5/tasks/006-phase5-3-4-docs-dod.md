# Task 006: Phase 5.3 + 5.4 — doc consolidation + branch definition-of-done

## Goal
5.3: a single consolidated GPU branch guide that pulls together the GPU-safe IC/BC
matrix, the build matrix (no-fast-math vs fast-math), and the elastic disposition.
5.4: a branch definition-of-done checklist with current status of each item.

## Context
The branch is never merged; these docs are the branch's self-contained record.
Source material already exists (link/summarize it; do NOT edit the originals):
- IC/BC matrix: `docs/gpu_safe_ic_bc_matrix.md`
- Build matrix: `benchmark/README.md` + build scripts (`benchmark/build_alamo_local_gpu.sh`,
  `benchmark/build_alamo_nova_3d.sh`); fast vs strict (no-fast-math) builds.
- Elastic disposition: `benchmark/archive/PHASE1_ELASTIC_DISPOSITION.md` (D1 = CPU-resident).
- Crossover: `benchmark/archive/PHASE3_R3_crossover.md` (D3 = WIN @ single).
- Roadmap Phase 5 "Branch definition-of-done" checklist: `~/Desktop/GPU-OPT-ROADMAP.txt`.

## Files to read first
- `~/Desktop/GPU-OPT-ROADMAP.txt` (Phase 5 + definition-of-done list)
- `docs/gpu_safe_ic_bc_matrix.md`, `benchmark/README.md`
- `benchmark/archive/PHASE1_ELASTIC_DISPOSITION.md`, `benchmark/archive/PHASE3_R3_crossover.md`

## Files allowed to modify
- `benchmark/GPU_BRANCH_GUIDE.md` (new)
- `benchmark/archive/PHASE5_BRANCH_DONE.md` (new)

## Files NOT allowed to modify
- The existing docs being summarized; any source/CI.

## Implementation steps
1. `GPU_BRANCH_GUIDE.md`: a top-level index/guide for the chamber-gpu branch.
   Sections: branch policy (never merged), build matrix (which binary for what:
   CPU, fast CUDA, no-fast-math CUDA, NOVA 3D; how to build each), GPU-safe IC/BC
   matrix (summary + link), elastic disposition D1 (summary + link), crossover D3
   (summary + link), where the phase reports live. Keep it a navigable map, not a
   copy of every doc.
2. `PHASE5_BRANCH_DONE.md`: the roadmap's branch definition-of-done as a checklist,
   each item annotated DONE / PARTIAL / PENDING with one line of evidence and a
   pointer (golden compare; no host-loop device-arena writes; device aborts/NaN
   detection; CPU regression suite [from task 001]; R3 win recorded). Reflect the
   no-merge policy (this is "branch done", not "merged").
3. Mark task-001 / task-004 / task-005 outputs as the evidence for the relevant
   items even if those are still landing (note "see PHASE4_R4_dispatch.md",
   "chamber-gpu-correctness.yml", etc.).

## Invariants
Index/summarize; do not duplicate or edit source docs. Accurate status only.

## Expected result
Two new docs: a navigable branch guide and an honest definition-of-done checklist.

## Non-goals
Editing existing docs; code; CI.

## Stop conditions
None expected.

## Final report: write results/006-RESULT.md

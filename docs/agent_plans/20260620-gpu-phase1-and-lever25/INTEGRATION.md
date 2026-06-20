# Integration: Lever 2.5 + Phase 1 elastic (2026-06-20)

Lead: Opus. Workers: 2× Sonnet, main-tree, ownership isolation. Base: `fe3844f31`.

## Completed tasks
- **Track B — Phase 1 elastic disposition (D1): COMPLETED.** Verdict delivered
  with data. R1 written (`benchmark/PHASE1_ELASTIC_DISPOSITION.md` extended) +
  `results/B-RESULT.md`.
- **Track A — Lever 2.5 field packing: BLOCKED (negative result, by design).**
  No source written; the pack was disqualified at design review before code.
  `results/A-RESULT.md`.

## Blocked / not-done
- **2.5 implementation** — architecturally infeasible within the task's allowed
  files. The framework's **implicit-component-0 convention** makes a "low-risk
  same-BC" eta+psi pack impossible without editing `Newton.H` / `Elastic.cpp/.H`
  (psi read at comp0) AND/OR every `IC::*::Add` (eta written at comp0) — no
  single component ordering satisfies both, and both are out of scope. Even if
  waived, the ceiling is ~2.5% wall/step (below the 5% bar), because FillBoundary
  is **53.7% of launch count but only 15.1% of kernel time** (cheap kernels).
- **ncu occupancy (step 1.3)** — environment-blocked (`ERR_NVGPUCTRPERM`,
  `perf_event_paranoid=4`, no passwordless sudo). No occupancy data obtainable
  here; not a finding about the kernels. Needs root / driver-param change.

## Conflicts / overlap
None. Verified post-run:
- `git diff --stat src/` empty (Track A made zero source edits).
- `git diff --stat input configure Makefile` empty (neither track touched shared
  build/input state).
- Track A wrote only nsys temp (`/tmp/lever25_nsys/`) + box-sweep outputs it owns;
  Track B wrote only `benchmark/phase1_elastic_2048/*`, `input_phase1_elastic`,
  and its private binary `bin/alamo_gpu_phase1-nofast`. No file owned by both.
- One transient GPU OOM during Track A's nsys capture (Track B's concurrent
  elastic job briefly took the 8 GB) — resolved by retry, exactly as the plan
  prescribed. The single-shared-GPU contention was real but non-fatal.

## Architecture-consistency check
- Golden compare green: gate **9/9 ok** (strict bit-identical), box sweep **5/5
  ok** — confirmed by Track A pre-task; unchanged since (no source touched).
- Device-safety invariants intact (no new host-loop device-arena writes).
- Both verdicts are consistent with the roadmap's D1/D2 decision trees and the
  master pivot tree.

## Tests run
- `baseline_suite.py check` → 9/9 ok (strict bit-identical).
- `phase2_box_sweep.py` → 5/5 ok (`wide_512_bf32_mgs128`: cpu 1.418s, gpu 1.835s,
  gpu/cpu 1.294).
- Elastic matched-resolution CPU np8 vs GPU strict (solve regions only):
  CPU 3.40× faster @512², 2.27× faster @1024²; GPU SIGABRTs @2048² (numerical,
  not OOM); CPU completes @2048².

## Recommended merge / commit order
1. **Already committed:** checkpoint `fe3844f31` (validated Phase 0 + 2.3).
2. **Commit (docs only, no code):** the orchestration record
   `docs/agent_plans/20260620-gpu-phase1-and-lever25/**` + report R1
   `benchmark/PHASE1_ELASTIC_DISPOSITION.md`. These capture the D1 verdict and the
   2.5 negative result — the durable deliverables. *Recommend, pending user OK.*
3. **No code to merge from either track** — Track A produced none; Track B
   produced none. Nothing to revert.
4. **Do NOT commit:** run artifacts (`benchmark/phase1_elastic_2048/*.log` +
   plotfiles, `benchmark/phase2_box_sweep/` outputs), the 275 MB private binary
   `bin/alamo_gpu_phase1-nofast` (remove to reclaim space), `input_phase1_elastic`,
   `/tmp/lever25_nsys/`.

## Strategic consequence (roadmap master pivot tree)
- **D1 = CPU-RESIDENT elastic** → "redirect P2/P3 effort to GPU phase-field +
  CPU-elastic split; this is the most likely best-realistic architecture."
- **D2 = residual is a regime problem** (reinforced: 2.5's launches are cheap) →
  "stop grinding launch levers; jump to **Phase 3** (regime scaling: 2D→3D,
  shallow AMR, multi-GPU)." Phase 3 needs NOVA A100/H100 — the A1000's 8 GB
  can't reach the saturating regime.

## Follow-ups
1. **Elastic operator/conditioning root-cause** (the real elastic next step):
   compare CPU vs GPU operator inputs at 2048² before `solver.solve`; minimal
   elastic-only repro; sweep `use_psi`/`psi_floor`/AMR-depth one knob at a time;
   inspect `Fapply`/`Diagonal`/`Newton::prepareForSolve` for nodal-index/ghost/
   mask/precision GPU-vs-CPU differences. **Investigate the 83-vs-126-iter
   nondeterminism** (atomic/reduction order) — it points beyond static conditioning.
2. **ncu occupancy** — needs GPU perf-counter permission (root / driver param).
3. **Lever 2.5** — shelved unless a future, larger task widens scope to
   `Newton.H` + `Elastic.cpp/.H` + `IC::*::Add` to pick a consistent component
   convention. Not worth it on the ~2.5% ceiling alone.
4. **Pivot to Phase 3** for any further GPU win-hunting (needs bigger hardware).

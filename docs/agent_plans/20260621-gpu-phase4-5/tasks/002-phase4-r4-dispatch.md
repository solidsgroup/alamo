# Task 002: Phase 4.2/4.3 — R4 dispatch report (D4 = ISOLATE)

## Goal
Write `benchmark/archive/PHASE4_R4_dispatch.md`: the Phase-4 report recording the D4
framework-dispatch decision and the de-fork status, per the roadmap R0–R4 template.

## Context
- D4 is FIXED by the no-merge policy: the branch never lands on master, so the
  de-virtualization stays contained in the CUDA build and must NOT be generalized
  framework-wide. **D4 = ISOLATE.** Do not re-argue it; justify and record it, and
  flag that formal ratification needs the group (Runnels) per Gate G4.
- 4.3 (de-fork) is ALREADY DONE: `src/GPU/IntegratorPolicy.mk` provides a
  principled, build-tracked, per-integrator GPU-clean policy
  (`ALAMO_GPU_SUPPORTED_INTEGRATORS := flame`, `ALAMO_GPU_MAIN/SOURCES`), consumed
  by the Makefile `findstring cuda` block; `src/alamo_gpu.cc` is the Flame-only
  main. The CPU launcher `alamo.cc` is untouched. Describe this as the realized
  isolation mechanism.
- The 4.1 CPU-semantics regression result will be provided to you in the dispatch
  prompt — fold its verdict in.

## Files to read first
- `~/Desktop/GPU-OPT-ROADMAP.txt` (Phase 4, D4 tree lines ~210-225, report template ~289)
- `src/GPU/IntegratorPolicy.mk`, `src/alamo_gpu.cc`, `Makefile` (lines ~66-80)
- `benchmark/archive/PHASE1_ELASTIC_DISPOSITION.md` (D1), `benchmark/archive/PHASE3_R3_crossover.md` (D3)

## Files allowed to modify
- `benchmark/archive/PHASE4_R4_dispatch.md` (new) ONLY.

## Files NOT allowed to modify
- Everything else.

## Implementation steps
1. Read the sources above.
2. Write R4 with sections: Objective; Method; Data (4.1 regression result +
   de-fork mechanism description); D4 decision-tree verdict = ISOLATE with the
   reasoning chain (no-merge → contained capability → P3 win is single-GPU Flame,
   not a framework need → ISOLATE); 4.3 de-fork status = DONE (describe
   IntegratorPolicy.mk, how to add an integrator); 4.4 pointer (convention doc);
   Recommendation + confidence; Open risks / what needs Runnels ratification (G4).
3. Be explicit that D4 ratification is a HUMAN decision still pending — do not
   claim it is ratified.

## Invariants
Honest status: ISOLATE is recommended/decided-by-policy but NOT yet ratified.

## Expected result
`benchmark/archive/PHASE4_R4_dispatch.md` exists, follows the template, records D4=ISOLATE.

## Non-goals
Code changes; ratifying D4; editing the roadmap.

## Stop conditions
If the 4.1 result is missing from the prompt, write R4 with 4.1 marked "pending".

## Final report: write results/002-RESULT.md

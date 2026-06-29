# R4: Phase-4 Framework-Dispatch Decision (D4)

Status: **D4 = ISOLATE** (recorded; needs human/Runnels ratification). The
de-fork (4.3) is already implemented; the CPU-semantics regression (4.1) is
captured and confirms the isolation holds. Roadmap: Phase 4, framework dispatch
generalization. Branch: `chamber-gpu` (never merged — see
`benchmark/GPU_BRANCH_GUIDE.md`).

## Decision D4 = ISOLATE

The GPU device-model is **contained in the CUDA build only**; the shared,
all-integrator CPU build is never altered to suit the GPU port. This is forced
by the hard branch policy (`chamber-gpu` is never merged to master, user
directive 2026-06-21): there is no merge for which a shared-build generalization
would pay off, so the lowest-risk option — isolate — is the correct one. Agents
document this decision; they do not re-litigate it. **Open:** Runnels
ratification that ISOLATE (vs. a general framework-dispatch refactor upstreamed
to master) is the intended end state.

### How ISOLATE is realized (4.3 de-fork — already DONE)

- `src/GPU/IntegratorPolicy.mk` declares `ALAMO_GPU_SUPPORTED_INTEGRATORS :=
  flame` and computes the GPU-clean source closure. The Makefile selects it only
  under `ifneq (,$(findstring cuda,$(POSTFIX)))` (`SRC_MAIN = $(ALAMO_GPU_MAIN)`,
  `SRC = $(ALAMO_GPU_SOURCES)`), so the CUDA build compiles `src/alamo_gpu.cc`
  (Flame-only main) and the device-clean subset.
- The CPU launcher `src/alamo.cc` (all 53 integrators) and the CPU build are
  unchanged. The de-virtualization (removing `virtual` from `Solid`/`BC`/
  `Operator` for nvcc) lives in shared headers but is behavior-preserving on the
  CPU — see the 4.1 evidence below.

## 4.1 verdict folded in (from `…/results/001-RESULT.md`)

The CPU-semantics regression (`scripts/runtests.py`, 2D) was run on
`chamber-gpu`: **113 run, 89 verified, 5 failed**. The de-virtualization is
**clean** — every genuine `Solid`-dispatch / elastic integrator passed
(`Eshelby`, `EshelbyFiniteKinematics`, `Inclusion`, `PlateHole`, `RubberPlateHole`,
`RubberPressurizedHole`, `RubberWithInclusion`, `UniaxialTension(Periodic)`,
`DynamicBar`, `Solid`, `ThermoElastic`, `VoronoiElastic`, `CompositeImpact`,
`Suture`, `TopOp`, `Fracture*`, `SCPThermalSandwich`, `SCPThermalVoid`).

The 5 failures are 2 flame test cases, and **neither is a dispatch-isolation
violation**:

1. `SCPSandwich` (×3) — runtime SIGSEGV at `Flame.cpp:677` (`L_out(i,j,k)=L`):
   `L_mf` is registered only under `if(thermal.on)` (Flame.cpp:180) but written
   unconditionally; crashes when `thermal.on=0`. master has no `L_mf` field and
   runs the config cleanly → chamber-lineage Flame bug, not de-virtualization.
2. `SCPSpheresElastic` (×2) — parse-time abort (`ParmParse.H:1297 query_exactly`
   on `elastic.model_prop.{lambda,mu,E,nu,kappa}`, exactly-2 required, 0 given);
   aborts in the Flame constructor before any solve. master parses the same
   input fine → chamber-gpu `model_prop` schema/parse change, not dispatch.

**Conclusion for D4:** ISOLATE is validated. The shared CPU build's runtime
dispatch is intact across all integrators; the GPU port did not perturb any
CPU integrator's semantics through the de-virtualization. The two red tests are
independent Flame defects (one runtime, one input contract) that the audit
surfaced and that must be fixed on the branch, but they do not bear on the
framework-dispatch decision.

## What remains (not blockers to recording D4)

- **Ratify ISOLATE** (Runnels): confirm the branch-contained device-model, not a
  master-upstreamed dispatch generalization, is the intended end state.
- **Fix the two Flame defects** so DoD item 4 can flip to DONE (tracked in
  `benchmark/PHASE5_BRANCH_DONE.md` item 4 and `…/results/001-RESULT.md`):
  register `L_mf` unconditionally (or guard the write); resolve the
  `model_prop` arity contract (update input or make optional).
- **3D regression leg** (best-effort) was not completed; 2D is the required bar.

## Evidence

- `docs/agent_plans/20260621-gpu-phase4-5/results/001-RESULT.md` (full 4.1
  analysis, reproduction, master-baseline comparison).
- `benchmark/phase4_cpu_semantics_20260621_030514/dim2.log` (suite log).
- `src/GPU/IntegratorPolicy.mk`, Makefile CUDA scoping (4.3 de-fork).

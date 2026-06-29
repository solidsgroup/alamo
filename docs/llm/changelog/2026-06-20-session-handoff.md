# GPU Optimization Session Handoff — 2026-06-20

Branch `chamber-gpu`. Roadmap: `docs/llm/ROADMAP.md`. Repo: `/home/jackplum/Projects/alamo`.
Project memory (auto-loaded): `~/.claude/projects/-home-jackplum-Projects-alamo/memory/` — esp. `gpu_perf_nsys_findings.md` (updated this session with all perf findings below).

---

## 1. What is DONE this session

### Phase 2.2 — thermo reduction/sync batching (committed)
- **Commit `02772a369`** on `chamber-gpu`: "GPU: batch Flame thermo diagnostics into one fused reduction per step". One file, 45+/26−.
- `Flame::TimeStepComplete` (`src/Integrator/Flame.cpp` ~539–592) previously ran **6 separate `MultiFab::min/max` calls per level every step** (each a device launch + host sync), plus a per-layer `Util::Message` debug print and a per-step `MONITOR thermal.on` print.
- Replaced with **one fused `amrex::ReduceOps<Max,Max,Max,Max,Min>`** accumulated across all levels into a single `ReduceData` (one `.value()` host sync), then 5 scalar `amrex::ParallelDescriptor::ReduceRealMax/Min`.
- **Critical correctness details (preserved):**
  - MPI all-reduce is MANDATORY: `MultiFab::max` all-reduces internally; raw `ReduceData` sees only local boxes. CPU gate runs np8 and the 64³ base grid is a single box → **7/8 ranks empty at level 0**. Without the all-reduce, rank-0's thermo.dat is wrong.
  - 0.0/1.0 floors preserved; min/max are order-independent → output **bit-identical** to old code. `max(0,allreduce(locals)) == allreduce(max(0,local))`.
  - `a_time` param commented out (now unused), avoids `-Werror=unused-parameter`.
- **Validated:** `baseline_suite.py check` **9/9** `fresh vs reference: ok` (cpu np8 + gpu_fast + gpu_strict, 3 cases); `phase2_box_sweep.py` **5/5** correctness ok.
- **Honest perf:** 2.2 removed ~17 syncs/step out of ~4,170 (<0.5%). Correctness-neutral cleanup; negligible as a perf lever. Wall-clock flat (init-dominated at 5 steps).

### Binaries — ALL rebuilt fresh with the refactor
- `bin/alamo-2d-clang++` (CPU, np8) — `./configure --comp=clang++ --dim 2 && make -j`
- `bin/alamo_gpu-2d-cuda86-g++` (fast CUDA, sm_86) — `COMP=g++ CUDA_FP=fast SMOKE=0 ./benchmark/build_alamo_local_gpu.sh`
- `bin/alamo_gpu-2d-nofast-cuda86-g++` (strict/no-fast-math CUDA) — `CUDA_FP=strict SMOKE=0 ./benchmark/build_alamo_local_gpu.sh`
- CUDA build is scoped to Flame → only changed TU recompiles + relink (~2 min). nvcc enforces `--Werror ext-lambda-captures-this`; the new `[=]` reduction lambda captures only local `Set::Patch` values (no `this`), so it passes both fast and strict.

### Measurement + D2 (Phase 2 decision tree, roadmap line 130)
- nsys: `.local/nsight/opt/nvidia/nsight-systems/2026.1.3/target-linux-x64/nsys`; `source benchmark/local_cuda_env.sh` first.
- **Parity config capture** (`analysis/results_phase22/`, new binary, 30 steps, phase-field only, wide_512_bf32_mgs128 overrides): **11,592 cudaLaunchKernel/step, 4,170 cudaStreamSynchronize/step (18.5% of API time).** Launch breakdown: **FillBoundary 56.3%**, AMR interp (CellConservativeLinear) 30.1%, copy/misc ~6%, **Flame's own per-box loops only ~3.8%** (Advance 230 @53µs + UpdateModel + Integrate + TagCells + Mechanics).
- **nsubsteps 2→1 capture** (`analysis/results_phase22_nsub1/`): launches 11,592→**2,026/step (5.72×)**, syncs 4,170→**961/step (4.34×)**. ⇒ **launch count is subcycle-multiplied (nsubsteps^level), not ghost-exchange-code-bound.**
- **D2 verdict:** still launch-bound, but the residual is a **regime problem (Phase 3)**, not fixable by 2D launch-code refactors. Code levers are low-ROI: 2.4 <4% of launches; 2.5 ~25% best case behind a real refactor. Wide-shallow config already at GPU/CPU ≈ 1.0.

---

## 2. Working-tree state (IMPORTANT)

`HEAD = 02772a369`. The working tree still has **uncommitted prior-session GPU work** (NOT mine, validated green by the G0 gate because binaries were built from the working tree):
- `src/Integrator/Flame.cpp` — NaN-tripwire re-arming in `Advance` (`Util::DeviceErrorFlag`/`SetDeviceError`/`AbortIfDeviceError`, Phase 0.1); tag-fusion **2.3** in `TagCellsForRefinement`; an `Integrate` change.
- `src/Util/Util.H` — `DeviceErrorFlag` infrastructure (+51 lines).
- `src/IC/{Expression,Laminate,PNG,PSRead,Trig}.H`, `src/BC/Operator/Elastic/Expression.H`, `src/Model/Propellant/{FullFeedback,Homogenize}.H`, `src/Util/BMP.H` — device-portability guards.
- `configure`, `benchmark/README.md` also modified; lots of untracked benchmark/analysis artifacts.

I committed **only** my 2.2 hunk (via `git apply --cached` of the single extracted TimeStepComplete hunk). Do not assume the rest is committed. The tag-fusion 2.3 is done+validated but uncommitted.

---

## 3. What I was ABOUT to do (user's last instruction)

User said: **"Spin up an agent to implement 2.5 now, spin up a second agent to move to phase 1 elastic."** (This OVERRIDES my recommendation to skip 2.5 — user wants it done anyway for diligence, in parallel with Phase 1.) I had NOT spawned them yet — user interrupted to request this handoff.

### ⚠️ Coordination warning (critical)
Both tracks touch the **same repo/working tree, single `bin/`, single `configure` state**. Running them concurrently in one tree = build/config races. Options:
- Give the **2.5 agent an isolated git worktree** (`isolation: "worktree"`), since it edits field code + rebuilds + re-gates.
- The **Phase 1 agent can start with the EXISTING `bin/alamo_gpu-2d-nofast-cuda86-g++`** (no-fast-math) for step 1.1 — no rebuild needed initially, so it can run in the main tree.
- Or run sequentially.

---

## 4. Plan / scope for each track

### Track A — Lever 2.5 (field packing) [user wants it despite low-ROI]
Goal: cut the 56%-of-launches FillBoundary population by packing co-evolving cell-centered fields into multi-component MultiFabs so one FillBoundary/FillPatch covers them.
- **Field structure** (`Flame.cpp` ~122–191, declared in `Flame.H`): evolving ghost-bearing fabs are `eta_mf`(bc_eta, ghost 3), `psi_mf`(bc_eta, ghost 2), `temp_mf`(bc_temp, ghost 3). `eta_old_mf`(g2)/`temp_old_mf`(g3) are swap partners (`std::swap` at `Flame.cpp:611`), evolving=false. The nghost=0 fabs (mdot, heatflux, L, alpha, temps) generate **no** FillBoundary — packing them saves nothing.
- **Blocker:** `eta` uses `bc_eta`, `temp` uses `bc_temp` (different physical BCs). The high-value eta+temp pack needs a **combined multi-component BC** (eta-cond on comp0, temp-cond on comp1). `BC::Constant` supports component ranges. Differing ghost counts (3 vs 2) → use max.
- **Framework:** `RegisterNewFab(Set::Field&, BC*, ncomp, nghost, name, writeout, evolving)` (`Integrator.cpp:340`). Per-fab `FillPatch` via `physbc_array` (`Integrator.cpp:394`, `1218`). Packing touches: field decls, registration, **every access** (`.array(mfi)`/`.Patch(lev,mfi)`/`[lev]`), BC construction, and the swap logic.
- **Realistic leverage:** ~25% fewer launches best case; even perfect FillBoundary removal only touches 56% and the 30% interp + everything else stays subcycle-multiplied. **High refactor risk.**
- **Must:** re-run `baseline_suite.py check` (9/9) + `phase2_box_sweep.py` (5/5) + re-capture nsys to confirm FillBoundary launch drop and wall/step vs the 5% threshold. Commit separately if green.

### Track B — Phase 1 elastic disposition (D1, roadmap lines 72–112)
- **1.1** Run **no-fast-math** GPU elastic at **2048²** — the case that SIGABRTs under fast-math. Use `bin/alamo_gpu-2d-nofast-cuda86-g++`. Record MLMG convergence + iteration count vs CPU on the identical system.
- **1.2** Branch (D1, line 83): converges ⇒ fast-math×conditioning confirmed, fix = the build-flag split (already have fast vs strict binaries); still diverges ⇒ operator/conditioning investigation (psi formulation, `psi_floor`, masked-operator conditioning at thin features).
- **1.3** Occupancy/register profile of `Fapply`/`Diagonal`/`Newton` with `maxrregcount` removed (ncu: `benchmark/g0_ncu_capture.sh`, `NCU=.local/nsight-compute/usr/bin/ncu`).
- **1.4** Fair elastic bench: GPU elastic vs **fully-subscribed CPU** MLMG at ≥1024², counting only solve regions.
- **Context (memory `gpu_perf_nsys_findings.md`, `chamber_gpu_port.md`):** GPU elastic SIGABRTs (`MLMG failing`) on first elastic solve at high res under `--use_fast_math` (`configure:703/705`); CPU solves the identical config fine; at coarse 64²-base GPU elastic was 1.65× FASTER than CPU (port works). The benchmark runs SKIP elastic via `elastic.tstart=1e9` — to exercise it set `elastic.on=1`, `elastic.type=static`, sane `tstart`, and the memory's stability settings (`psi_floor`, clamped BCs). Prior artifacts exist: `benchmark/phase1_elastic_2048/`, `benchmark/PHASE1_ELASTIC_DISPOSITION.md` (untracked).

---

## 5. Key commands

```bash
# CUDA env (for nsys / running GPU binary)
source benchmark/local_cuda_env.sh

# G0 correctness gate (compares thermo.dat vs benchmark/baseline_references/; ALL columns)
CPU_NP=8 GPU_FAST_NP=1 GPU_STRICT_NP=1 python3 benchmark/baseline_suite.py check

# Box sweep (cpu_np8 vs gpu_fast, 5 cases, 5 steps) -> benchmark/phase2_box_sweep/{README.md,summary.csv}
python3 benchmark/phase2_box_sweep.py

# Rebuild GPU fast / strict (Flame-scoped, ~2 min)
COMP=g++ CUDA_FP=fast   SMOKE=0 ./benchmark/build_alamo_local_gpu.sh
COMP=g++ CUDA_FP=strict SMOKE=0 ./benchmark/build_alamo_local_gpu.sh
# Rebuild CPU
./configure --comp=clang++ --dim 2 && make -j"$(nproc)"

# nsys launch/sync capture (parity config example)
NSYS=.local/nsight/opt/nvidia/nsight-systems/2026.1.3/target-linux-x64/nsys
"$NSYS" profile -o OUT/nsys --force-overwrite true --stats=false --trace=cuda,nvtx --sample=none --cpuctxsw=none \
  bin/alamo_gpu-2d-cuda86-g++ input allow_unused=True max_step=30 stop_time=1e99_s \
  elastic.tstart=1000000000.0 elastic.solver.verbose=0 elastic.print_model=0 \
  amr.plot_int=-1 amr.thermo.plot_int=-1 amr.thermo.int=1 \
  amr.n_cell="64 64 64" amr.max_level=3 amr.blocking_factor=32 amr.max_grid_size=128 \
  amr.grid_eff=0.9 amr.base_regrid_int=1000000 amr.nsubsteps=2 plot_file=OUT/plt
"$NSYS" stats --force-export=true --report cuda_api_sum --report cuda_gpu_kern_sum --format csv --output OUT/stats OUT/nsys.nsys-rep
```
Best box-sweep case = `wide_512_bf32_mgs128` (n_cell 64³, max_level 3, blocking 32, max_grid 128, grid_eff 0.9, nsubsteps 2), GPU/CPU ≈ 1.0.

---

## 6. Task list (TaskCreate IDs this session)
- #5 Measure 2.2 + apply D2 — **completed**
- #6 Implement lever 2.5 (field packing) — **pending** (Track A)
- #7 Re-run gate + sweep + nsys after 2.5 — **pending**
- #8 Commit Phase 2.2 — **completed** (`02772a369`)
- #9 Phase 1 elastic disposition (D1) — **pending** (Track B)

(#1–#4 from earlier: CPU build, GPU build, G0 gate, box sweep — all completed.)

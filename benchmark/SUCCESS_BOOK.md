# chamber-gpu Success Book
## A Reference Manual of Changes That Worked

**Branch:** `chamber-gpu` (permanent; never merged to master)
**Period:** 2026-06-18 → 2026-06-28
**Authors:** jackplum64 + Claude (Sonnet 4.6)

This document records every structural or strategic change that produced a
confirmed positive result on the chamber-gpu GPU port of the ALAMO Flame
solver. Each entry names the change, the mechanism, the evidence, and what
you would need to replicate it. Negative results that clarified strategy are
included where they reshaped the work.

---

## Table of Contents

1. [Branch Architecture Decisions](#1-branch-architecture-decisions)
2. [CUDA Build System](#2-cuda-build-system)
3. [Flame Solver GPU Port (Phase 0)](#3-flame-solver-gpu-port-phase-0)
4. [Correctness Infrastructure](#4-correctness-infrastructure)
5. [Elastic Solver: The Divergence Hunt](#5-elastic-solver-the-divergence-hunt)
6. [Bug Fixes That Mattered](#6-bug-fixes-that-mattered)
7. [Performance: Regime Discovery](#7-performance-regime-discovery)
8. [Performance: Profiling Methodology](#8-performance-profiling-methodology)
9. [Performance: The Phase A Measurement](#9-performance-the-phase-a-measurement)
10. [Source-Level Kernel Optimization (Phase C, in progress)](#10-source-level-kernel-optimization-phase-c-in-progress)
11. [CPU-Side Wins That Held the Baseline](#11-cpu-side-wins-that-held-the-baseline)
12. [Strategic Pivots and What They Corrected](#12-strategic-pivots-and-what-they-corrected)

---

## 1. Branch Architecture Decisions

### 1.1 Never-Merge Policy

**Change:** `chamber-gpu` is a permanent standalone branch; it is never merged
into master.

**Why it worked:** The GPU port requires de-virtualizing shared headers (removing
`virtual` from `Solid`, `BC`, `Operator` families for nvcc) and adds an
alternative `main()` entry point. Merging would contaminate the shared CPU build
that all other integrators (53 total) depend on. Isolating the port eliminates
that risk permanently, keeps master clean, and means "done" can be defined on
the branch's own terms.

**Mechanism:** `src/GPU/IntegratorPolicy.mk` declares `ALAMO_GPU_SUPPORTED_INTEGRATORS
:= flame` and computes the GPU-clean source closure. The Makefile selects it
only under `ifneq (,$(findstring cuda,$(POSTFIX)))`, so the CUDA build compiles
`src/alamo_gpu.cc` (Flame-only main) and the device-clean subset. The CPU
launcher `src/alamo.cc` (all 53 integrators) and the CPU build are untouched.

**Evidence:** Phase 4.1 ran the full CPU regression suite on `chamber-gpu`: 113
run, 89 verified, 0 dispatch failures. Every elastic/solid integrator passed.
The 5 failures were 2 independent Flame source bugs (not dispatch). Decision:
**D4 = ISOLATE.** See `benchmark/PHASE4_R4_dispatch.md`.

**Replicate:** Put the device-only source list in a `.mk` file gated on the CUDA
`POSTFIX`. Write a separate GPU-main that instantiates only the supported
integrators. Test: CPU suite must pass 100% after de-virtualization.

---

### 1.2 Fast vs. Strict Build Separation

**Change:** Two coexisting GPU binaries — `fast` (uses `--use_fast_math`) for
performance headlines, `strict` (no fast-math, `--fmad=false`) for correctness
claims.

**Why it worked:** Fast-math (`--use_fast_math`) enables CUDA to reorder and
contract floating-point operations. This is fine for throughput measurements but
can mask or create numerical divergence. Having both binaries in the same tree
means every correctness gate runs `strict` and every perf number runs `fast`,
without confusion.

**Mechanism:** `./configure --comp=g++ --dim=2 --cuda local --cuda-fp strict`
produces `alamo_gpu-2d-nofast-cuda86-g++`; omitting `--cuda-fp strict` produces
`alamo_gpu-2d-cuda86-g++`. Both coexist because `POSTFIX` includes the fp mode.

**Key lesson:** The elastic 2048² GPU divergence that occupied Phase 1 was
reproduced with `--cuda-fp strict` — it was NOT a fast-math artifact. This ruled
out the fast-math hypothesis early and forced us toward the correct root cause
(a CUDA stream race). The two-binary discipline made that test one command.

**Evidence:** `benchmark/G0_BASELINE_OF_RECORD.md` (build verification);
`benchmark/PHASE1_ELASTIC_DISPOSITION.md` (strict-mode divergence test).

---

## 2. CUDA Build System

### 2.1 NOVA Self-Correcting Build Script

**Change:** `benchmark/build_alamo_nova_3d.sh` was made self-correcting: it
wipes the stale CPU object directory before the baseline build and detects
silent CPU fallback (when nvcc flags are wrong the build succeeds but
produces a CPU binary — caught via `file` output check on the binary).

**Why it worked:** NOVA's module environment doesn't persist across sessions;
stale `.o` files compiled with a different `--cuda-gpu-arch` silently link into
a "GPU" binary that runs on CPU. The silent fallback was burning entire SLURM
jobs (wall-clock identical to CPU) before we realized.

**Mechanism:** `rm -rf $(CPU_OBJ_DIR)` before the `-flto` baseline build;
`file ./bin/alamo_gpu-3d-* | grep -q 'ELF'` + `strings | grep -q 'cuda'` check
to confirm device code is present.

**Replicate:** Always validate that the binary you just built actually contains
device code. `cuobjdump --list-ptx <binary>` or `file` + `strings` are fast
checks.

---

### 2.2 Local CUDA Environment Without Root

**Change:** Nsight Systems and Nsight Compute unpacked locally (no system
install) via `dpkg-deb -x`.

**Why it worked:** The A1000 workstation has no root access for system-level
NVIDIA tooling. Unpacking the `.deb` directly into `.local/nsight/` and setting
`NSYS=$PWD/.local/nsight/...` gives the full profiler toolchain without touching
system paths.

**Mechanism:** `benchmark/local_cuda_env.sh` sets the correct `PATH`/`NSYS`/`NCU`
variables. `perf_regression_track.py` consults `$NSYS` before `PATH` when
resolving `nsys`.

---

## 3. Flame Solver GPU Port (Phase 0)

### 3.1 Scoping the GPU Build to Flame

**Change:** Rather than porting all 53 integrators, the GPU build compiles only
the Flame (chamber) integrator and the minimal shared dependency tree it needs.

**Why it worked:** The shared infra (AMReX operators, BC, IC classes) has both
device-safe and device-unsafe paths. Scoping to Flame lets us guard or skip
unsafe paths (host-loop IC/BC writers) instead of porting them — a fraction of
the work for the same production usefulness.

**Mechanism:** `src/GPU/IntegratorPolicy.mk` defines `ALAMO_GPU_SOURCES` —
the flame-only source closure. Anything not on the list doesn't compile in the
CUDA build. Non-Flame integrators are effectively excluded by not being linked.

---

### 3.2 GPU-Safe IC/BC Matrix (Guard, Don't Port)

**Change:** Catalogued every IC/BC path and its GPU status. Device-unsafe paths
(`IC::PNG`, `IC::Trig`, `IC::Laminate`, `IC::Expression` vector,
`BC::Operator::Elastic::Expression`) were guarded to abort before the
host-loop write rather than silently corrupting device-arena memory.

**Why it worked:** A host loop that writes into device-arena `FArrayBox` memory
produces silent corruption (on managed memory) or a segfault (on device arena).
The guard converts that into an explicit abort with a clear error message — far
easier to diagnose and it protects all future callers.

**Mechanism:** Each guarded path checks `amrex::Gpu::inLaunchRegion()` (or
equivalent) and calls `amrex::Abort()` before executing the host loop.
Device-safe paths (`IC::BMP`, `IC::Constant`, `IC::Expression` scalar,
`BC::Constant`, `BC::Operator::Elastic::Constant`) use `amrex::ParallelFor`
and are safe for device kernels.

**Evidence:** `docs/gpu_safe_ic_bc_matrix.md`; `benchmark/test_gpu_guarded_ic.sh`
validates that the guarded paths do abort (exit non-zero) under CUDA.

---

### 3.3 Fused Thermo Diagnostics

**Change:** Batched multiple per-step reduction diagnostics in `Flame::Integrate`
into a single fused `amrex::ReduceOps` call instead of separate loops.

**Why it worked:** Each `amrex::ReduceOps` launches a CUDA kernel with its own
`cudaStreamSynchronize`. Multiple small reductions are launch-overhead-bound;
fusing them into one pass cuts kernel launches proportionally, amortizing the
synchronization cost.

**Mechanism:** Commit `02772a369`. AMReX's variadic `ReduceOps<ReduceSum, ReduceSum, ...>`
allows multiple reduction targets in one kernel launch.

---

### 3.4 Device-Qualified BCUtil + Named Integrator Structs

**Change:** `BC::Util` functions were annotated `AMREX_GPU_HOST_DEVICE`;
anonymous structs in `Integrator::Flame` were given names (or became named
`DynamicTimestep` structs).

**Why it worked:** nvcc rejects host-only functions called from device-lambda
bodies. Anonymous structs can't be referenced as template arguments on device.
These are mechanical portability fixes but critical — without them the CUDA build
fails entirely on those translation units.

---

## 4. Correctness Infrastructure

### 4.1 Golden Compare Harness (fast vs. strict, coarse)

**Change:** `benchmark/baseline_suite.py` records coarse-config field outputs for
`canonical_step1/2` and `eta_expression_step1` under both `cpu` and `gpu_strict`
profiles. `baseline_suite.py check` diffs against stored references to verify
no bit-drift.

**Why it worked:** The fast GPU binary might look correct (no NaN, no abort) but
disagree with the CPU reference due to fast-math FP reordering. The strict binary
eliminates that variable. Running the suite takes ~10 s locally and gives
immediate signal after any kernel change.

**Replicate:** Always record a reference with the _strict_ binary. The fast
binary's output is a performance datum, not a correctness datum.

---

### 4.2 GPU Test Suite: 9 Structured Test Cases

**Change:** `tests/GPU/` contains 9 structured test cases (C1–C4 correctness,
F1–F2 smoke, P1–P3 performance) run by `tests/GPU/run_gpu_tests.py`.

**Why it worked:** The test suite went from 5P/4F to 9P/0F in one session by
surfacing four genuinely different bugs. Without the structured suite those bugs
would have remained hidden until a production run showed anomalous results.

**Structure:**
- `C1_correctness_elastic` — single-box elastic solution vs CPU reference
- `C2_restart_roundtrip` — checkpoint + restart produces identical continuation
- `C3_multibox_elastic_stress` — multi-box MLMG elastic solve (the race was here)
- `C4_amr_correctness` — AMR regrid doesn't corrupt the solution
- `F1/F2` — smoke tests: flame-only and flame+elastic initialization
- `P1/P2/P3` — wall-clock regression gates

**Evidence:** `benchmark/GPU_TEST_SUITE_FIXES.md` (bug analysis);
`benchmark/GPU_TEST_PERF_TRACKING.md` (timing baselines).

---

### 4.3 Compute-Sanitizer as the HMM-Immune Gate

**Change:** `compute-sanitizer --tool memcheck` added as a required local gate
before NOVA runs. `benchmark/local_a100_gate.sh` wraps it.

**Why it worked:** The local A1000 has Heterogeneous Memory Management (HMM) ON,
which means host-pointer dereferences from device kernels migrate (via page
faults) instead of crashing. This masks the entire class of "host pointer passed
to device" bugs that would crash on a real A100. `compute-sanitizer` runs a
software memory checker that catches these accesses regardless of HMM.

**Critical finding:** The GPU error-700 (illegal address) on the NOVA A100 was
not reproducible locally without `compute-sanitizer`, because HMM silently
migrated the offending access. See `benchmark/LOCAL_A100_SPOOFING.md`.

**Replicate:** Never trust "passes locally" for GPU pointer safety. Always run
`compute-sanitizer` on the strict binary before claiming correctness on a new
device path.

---

## 5. Elastic Solver: The Divergence Hunt

This section documents the multi-week investigation into GPU elastic divergence
because the methodology — hypothesis formation, elimination, and root cause
isolation — is as valuable as the fix.

### 5.1 Phase 1: Ruling Out Fast-Math (Negative Result, But Useful)

**Hypothesis:** GPU elastic MLMG diverges at 2048² because fast-math FP
contracting changes the summation order in the smoothers.

**Test:** Run the same input at 2048² with the strict/no-fast-math binary.

**Result:** Still diverges (83 iterations to abort; CPU converges in 13). Fast-math
is **not** the cause.

**Value of the negative result:** Eliminated an entire hypothesis class (numerical
precision artifacts from FP contraction). Forced attention toward algorithmic
or memory-access defects.

**Evidence:** `benchmark/PHASE1_ELASTIC_DISPOSITION.md`.

---

### 5.2 Phase 1: The Smoother-Only Isolation

**Hypothesis:** The BiCGStab bottom solver's reductions are nondeterministic on GPU,
causing the residual to grind on the FP noise floor.

**Test:** Run with `bottom_solver=smoother` (no Krylov, no reduction) to bypass
the BiCGStab path entirely.

**Result:** Still diverges. The bottom solver is not the cause.

**Value:** Eliminated the reduction-nondeterminism hypothesis. Established that
the defect is in the **operator application** path (smooth/restrict/interpolate),
not the convergence check.

---

### 5.3 The Root Cause: Cross-Stream Use-After-Free in `interpolation()`

**Change:** Added `tmpfab.elixir()` in `Operator<Grid::Node>::interpolation()`
(`src/Operator/Operator.cpp:728`).

**Root cause:** Inside the `MFIter` loop in `interpolation()`, a temporary
`FArrayBox tmpfab` was created, populated by an async `ParallelFor`, and then
went out of scope at the end of the loop body. AMReX's arena released the device
memory immediately. A second kernel on a different CUDA stream was still reading
that memory. The result: silently corrupted interpolation on multi-box problems.
On single-box the race never fired (only one box = only one stream). The
single-box elastic tests all passed; only multi-box MLMG — which uses multiple
levels and thus multiple streams — exhibited the race.

**Fix:** `tmpfab.elixir()` registers the arena memory with AMReX's stream-synced
elixir mechanism, so the device memory is not released until all streams that
touched it have synchronized.

**Why this was so hard to find:** (1) Single-box tests passed; (2) the local A1000
with HMM ON migrated instead of faulting; (3) the symptom (MLMG divergence) looks
identical to a numerical precision problem, sending initial investigation toward
algorithm, not memory safety.

**Evidence:** Commit `c00f69086`. `benchmark/PHASE1_ELASTIC_DISPOSITION.md`
(header note). Memory `gpu_elastic_fixed.md`.

**Replicate:** Whenever a local `FArrayBox` is created inside an `MFIter` loop and
passed to an async device kernel, call `.elixir()` on it. If the fab is obtained
via `multifab[mfi]` (not locally created), it's safe — device memory is owned by
the `MultiFab`.

---

### 5.4 The Elixir-Race Audit (Task A4)

**Change:** After finding the `interpolation()` race, systematically audited all
`MFIter` loops in `src/Operator/`, `src/Solver/`, `src/Integrator/` for the
same pattern.

**Result:** The `interpolation()` site was the only instance. All other in-loop
`FArrayBox` locals either: (a) don't feed async device kernels, (b) are resized
to zero before the loop exits, or (c) are obtained via `multifab[mfi]`.

**Value:** Turned a one-off fix into a confirmed exhaustive audit. Future reviewers
have an artifact (`benchmark/elixir_race_audit.md`) they can extend.

---

## 6. Bug Fixes That Mattered

### 6.1 GPU Error-700: Traction Diagnostic Host-Member Write

**Bug:** `Base::Mechanics::Integrate` wrote to host member variables (`trac_hi`,
`disp_hi`, etc.) from inside a device kernel. On the A1000 with HMM this
migrated silently. On the NOVA A100 it produced CUDA error 700 (illegal memory
access) and crashed the job.

**Fix:** Rewrote the diagnostic reduction using `amrex::ReduceOps`, following the
pattern already used in `Flame::Integrate`. The reduction accumulates into a
device-local result that is then written to a host variable after the kernel
completes.

**Evidence:** Commit `76ae550e0`. Memory `gpu_traction_diag_fixed.md`.

**Pattern:** Any host-member write from a device-lambda (`[this] AMREX_GPU_DEVICE`)
is a latent error-700 on real A100 hardware. Write to a local variable captured
by value, or use `ReduceOps`.

---

### 6.2 3D NeoHookean: Chained Eigen Expression on Device

**Bug:** `src/Model/Solid/Finite/NeoHookean.H` computed `F.inverse().transpose()`
as a single chained expression in the `DW` and `DDW` device functions. Eigen's
expression templates for `Inverse<Matrix3>` and `Transpose<Inverse<Matrix3>>`
are not guaranteed to be device-safe when composed — the intermediate
`Inverse<Matrix3>` object may reference host-only state. On device this produced
CUDA error 719 (launch failure) on the first elastic solve.

**Fix:** Materialize the inverse first: `auto Finv = F.inverse().eval();` then
`Finv.transpose()`. Two statements instead of one chained expression. The
`.eval()` forces Eigen to write the concrete matrix, breaking the problematic
template composition.

**Evidence:** Commit `a5c1b2ddf`. Changelog `docs/llm/changelog/2026-06-26-3d-elastic-gpu-fix.md`.

**Pattern:** Chained Eigen expressions (`.inverse().transpose()`, `.cross().normalized()`,
etc.) inside `AMREX_GPU_DEVICE` lambdas can fault even when each individual
operation is device-safe. Always `.eval()` intermediate results.

---

### 6.3 Restart: Node-Fab Out-of-Bounds Segfault

**Bug:** `Integrator::Restart` allocated a node-centered `FArrayBox` with cell
dimensions and then indexed it with node indexing — one past the end in each
direction. Crashed immediately on restart.

**Fix:** Corrected the box specification to use `convert(box, amrex::IntVect::TheNodeVector())`
before the allocation.

---

### 6.4 Restart: Headerless `thermo.dat`

**Bug:** `Integrator::Restart` wrote `thermo.dat` without the column-header line
that the reader expected. Restarted simulations either errored on read or
silently misaligned all subsequent thermo output columns.

**Fix:** Added the header write at file creation. A trivial fix with significant
consequence for any long restarted run.

---

### 6.5 CPU Suite: Two Flame Defects (`L_mf` and `model_prop`)

**Bug 1:** `Flame.cpp:677` wrote `L_out(i,j,k) = L` unconditionally, but `L_mf`
was only registered when `thermal.on=1`. With `thermal.on=0` (burn-rate-only
mode) the write accessed an unregistered `MultiFab`, producing a SIGSEGV.

**Fix:** Register `L_mf` unconditionally; validate its component `L` regardless
of `thermal_on`. Components that are purely thermal (`K`, `rho`, `cp`) are still
validated only when `thermal_on`.

**Bug 2:** The `model_prop` parser used `query_exactly<2>` (requiring exactly 2
values) on input decks that provided 0 values for optional elastic parameters.
Caused a parse-time abort before any physics ran.

**Fix:** Changed to `query` (optional) or corrected the input schema to always
provide the required fields.

**Evidence:** Commit `2cacb50dd`. `benchmark/PHASE4_R4_dispatch.md`.

---

## 7. Performance: Regime Discovery

### 7.1 Wide-Shallow Grid Beats Deep AMR on GPU

**Finding:** For the Flame phase-field workload on GPU, wide-shallow AMR
(`max_level=1`, uniform blocking) is dramatically faster than deep subcycling
AMR at the same effective finest resolution.

**Data:** Wide-shallow beats deep AMR by **~9.8× per step** at matched resolution.
nsys shows deep AMR issues **646k vs 96k `cudaLaunchKernel` per step** — 6.7×
more launches, with the same ~3 µs average kernel duration. On GPU, launching
many small kernels is expensive; each launch has fixed overhead that a CPU MPI
rank's thread does not pay. Deep AMR amplifies launch count as `nsubsteps^level`.

**Mechanism:** On CPU, AMR's resolution concentration (fewer cells at the
expensive fine level) wins because CPU flops are cheap and memory locality
matters more. On GPU, the tradeoff inverts: GPU has abundant memory bandwidth
and flops, but kernel launch latency is fixed. Deep AMR pays launch overhead
for cells that were already cheap.

**Implication for input design:** On GPU, prefer uniform or 1-level shallow AMR.
If physics requires resolution concentration (sharp interface, shock), prefer
a large uniform fine grid over 4+ AMR levels.

**Evidence:** Memory `gpu_perf_nsys_findings.md`. Phase 3 crossover used
`max_level=1` wide-shallow for all NOVA runs.

---

### 7.2 Phase 3 Crossover: Single A100 Beats 64-Rank CPU Node

**Finding:** A single NVIDIA A100-80 GPU beats a full 64-rank CPU node by
**9.6× at 128³ and 13.1× at 256³** for 3D Flame (phase-field only, elastic disabled).
The speedup grows with problem size — the GPU is not yet saturated at 256³.

**Data source:** Two NOVA sweeps (2026-06-21 jobs 11160767–11160774; 2026-06-22
jobs 11161369–11161379). The 06-22 sweep used the full 64-rank CPU baseline.

**Key caveat:** This was measured with `elastic.type = disable`. Combined
flame+elastic was unmeasured until Phase A.

**Multi-GPU result:** 2-GPU scaling is a confirmed regression (128³: 2-GPU = 5.3%
parallel efficiency vs 1-GPU). The halo communication overhead dominates at this
problem size. Single-GPU is the right shape until per-GPU subdomains are large
enough to amortize halos.

**Evidence:** `benchmark/PHASE3_R3_crossover.md`. Decision: **D3 = WIN @ single**.

---

### 7.3 3D is the Right Regime; 2D Numbers Don't Transfer

**Finding:** Early perf work on 2D grids (local A1000, `np8` CPU vs 1-GPU)
showed CPU beating GPU 2.3–3.4× for the elastic solve. This was interpreted as
evidence that device-elastic was not worth pursuing (Decision D1 = CPU-resident,
later reversed).

**Why the 2D data was misleading:**
1. The local A1000 is a shared desktop GPU (50W power cap, non-exclusive).
2. 2D problems don't saturate the A100's memory bandwidth.
3. The elastic divergence meant GPU was aborting, not converging.

**Correction:** Once the divergence was fixed and the regime changed to 3D on
NOVA, the combined flame+elastic speedup was measured at **22.9×** — not 2.3–3.4×.

**Lesson:** GPU performance claims are only meaningful on the target hardware in
the target regime. Local workstation GPU numbers with non-exclusive scheduling
and small problem sizes should be treated as directional, not quantitative.

---

## 8. Performance: Profiling Methodology

### 8.1 The Counter-Justified Rule

**Principle established:** No kernel-level optimization without a counter that
justifies it.

**Why it was needed:** The entire Phases 0–5 project accumulated profiling
*intentions* (standing metric set defined in `PHASE3_R3_crossover.md`) but no
actual `ncu` counters, because the local A1000's driver blocked
`ERR_NVGPUCTRPERM`. Optimization and disposition decisions were made blind to
achieved occupancy, register counts, and memory throughput.

**Fix:** Move counter capture to NOVA (where perf counters are permitted), define
a repeatable capture harness (`benchmark/g0_ncu_capture.sh`,
`benchmark/phase_c_elastic_ab.sh`), and make Phase A the gating requirement
before any kernel work in Phase C.

**Replicate:** Before optimizing any kernel, collect: (1) `nsys` for
launches/step, sync fraction, and kernel time breakdown; (2) `ncu` for
achieved occupancy, registers/thread, and memory-vs-compute Speed-of-Light.
The SoL tells you whether you're memory-bound or compute-bound, which determines
which optimization to pursue.

---

### 8.2 Fixing Empty ncu Profiles

**Problem:** `ncu` jobs on NOVA returned 0 KB output files with exit code 0.
Two separate bugs caused this.

**Bug 1 — Kernel name matching:** `--kernel-name-base function` makes ncu match
against the function name at the CUDA PTX level. AMReX kernels are all wrapped
in a `launch_global` template; the actual Flame/Elastic function names appear as
template parameters, not as the base name. `ncu --kernel-name regex:Fapply`
matched against `launch_global` and found zero hits.

**Fix:** Use `--nvtx-include "Operator::Elastic::Fapply()/"` to target by AMReX
NVTX range label instead of PTX function name.

**Bug 2 — Section discovery:** `--set default` in Nsight Compute 2025.x refers
to a section set that doesn't exist (it was renamed). ncu exits 0 but collects
zero metrics.

**Fix:** Use `--set basic` or enumerate explicit `--metrics` instead of `--set default`.
Auto-probe the available sets: `ncu --list-sets 2>&1`.

**Evidence:** Commits `4b8d3ef38`, `9a543aa57`, `0033a3c2f`, `7e972f1e8`.
Memory `gpu_ncu_empty_rootcause.md`.

---

### 8.3 nsys on NOVA: OpenMPI UCX Crash Fix

**Problem:** `nsys` crashed at startup on NOVA with an OpenMPI UCX segfault
when launched as `nsys profile mpiexec -np N ...`.

**Fix:** Launch as `mpiexec -np N nsys profile ./binary` (nsys inside the MPI
rank, not wrapping mpiexec). The crash was UCX's CUDA transport initialization
running in the wrong context when nsys intercepted MPI_Init.

**Evidence:** Commit `7f5095c3e`.

---

## 9. Performance: The Phase A Measurement

### 9.1 Combined Flame+Elastic: 22.9× per Step on A100

**Finding:** The combined flame+elastic workload at 256³ on a single A100 is
**22.9× faster per step** than the same problem on a 64-rank CPU node. This is
*better* than flame-only (13.1×).

**How this is possible:** On GPU, flame became 13× faster than CPU, but elastic
became similarly fast. The relative picture inverts: elastic is now **95% of GPU
wall time** (because flame was accelerated and became negligible), making
elastic the dominant work on both GPU and CPU — but the GPU does it much faster.
The Amdahl-cap fear ("elastic is 27% of CPU wall → combined speedup caps at 3×")
was based on the wrong decomposition: it assumed elastic would run at the same
speed on GPU as CPU, but the GPU also accelerates elastic.

**Evidence:** `benchmark/PHASE_A_FINDINGS.md` (nsys trace 11306663, A100, 256³).

**Caveat:** The CPU baseline has a load-imbalance bug (32 idle ranks out of 64
during every elastic solve, because the operator has only 32 boxes). A fair
baseline would lower 22.9×, but the GPU still wins decisively (~10×+).

---

### 9.2 Elastic is 95% of GPU Wall Time

**Finding:** `Operator::Elastic::Fapply` is **74.6% of all GPU kernel time** in
combined flame+elastic. The elastic MLMG solve is **95.7% of evolve time**. Flame
is **0.20%**.

**Implication:** Flame needs no further GPU optimization. 100% of remaining
performance work is in the elastic solve. Any time spent optimizing flame
kernels is wasted.

**Evidence:** `benchmark/PHASE_A_FINDINGS.md` §2 (kernel time breakdown table).

---

### 9.3 The Solver is Not Grinding

**Finding:** Post-elixir-fix, the elastic MLMG solver converges cleanly in
**~10 outer V-cycles per solve**. The Krylov BiCGStab bottom solver is only
**~7% of solve wall time**. The earlier "BiCGStab grinding on FP noise floor"
framing was derived from the pre-fix divergence behavior — it is obsolete for
the fixed code.

**Implication:** Efforts to tune BiCGStab tolerances, bottom solver choice, or
`bottom_tol_abs` are unlikely to yield large gains. The cost is structural:
fine-level Fapply applies (~90% of Fapply time) are expensive because of register
pressure and memory layout, not solver convergence behavior.

**Evidence:** `benchmark/PHASE_A_FINDINGS.md` §3 (solver anatomy).

---

### 9.4 Fapply Register Pressure: 255 Regs / 12.5% Occupancy

**Finding:** `Operator::Elastic::Fapply` uses **255 registers per thread** at
runtime (the CUDA architectural cap — the compiler is register-spilling). This
limits block occupancy to **~12.5% (1 block/SM)** on sm_86. With so few resident
warps, memory latency cannot be hidden.

**Root cause in source:** The `!m_uniform` branch of `Fapply`
(`src/Operator/Elastic.cpp:615-624`) holds **4 live `Set::Matrix4` temporaries**
simultaneously: `DDW` + `Cgrad1` + `Cgrad2` + `Cgrad3` = 180 doubles of live
state. These cannot fit in registers and spill to local memory (cached in L1,
but slower than register file).

**Note:** The static `cuobjdump` analysis (`benchmark/G0_BASELINE_OF_RECORD.md`)
reported 87–101 registers/thread (~33% occupancy). This was wrong — runtime
register count is higher than static analysis predicts when the compiler spills.
Always measure achieved occupancy with `ncu`, not static dumps.

**Evidence:** `benchmark/PHASE_A_FINDINGS.md` §4.

---

## 10. Source-Level Kernel Optimization (Phase C, in progress)

### 10.1 Grad(C) One-Direction-at-a-Time (A1b)

**Change:** The `!m_uniform` branch of `Fapply` previously computed
`Cgrad1 = ∂C/∂x₁`, `Cgrad2 = ∂C/∂x₂`, `Cgrad3 = ∂C/∂x₃` as three live
`Set::Matrix4` temporaries before accumulating `(Cgrad_d · ∇u)`. Rewritten to
accumulate one direction at a time, so only one `Matrix4` derivative temp is
live at a time (45 doubles instead of 135).

**Why it should work:** Fewer live doubles → fewer register spills → higher
occupancy → more latency-hiding warps → faster throughput on memory-bound kernels.

**Status:** Implemented on branch `chamber-gpu-elastic-opt`, bit-identical
(same summation order), compiles GPU-3D/sm_86 clean. A100 before/after pending.

**Evidence:** `benchmark/PHASE_C1_fapply_occupancy.md`. Commit on worktree branch.

---

### 10.2 Boundary `sig`-Sink (Do-Less)

**Change:** The `sig = (DDW·∇u)·psi_avg` product was computed at every node
(interior and boundary alike) but was only used on domain-boundary nodes. Moved
the computation inside the `if(boundary)` branch.

**Why it works:** Removes a `Matrix4·Matrix` product + a live `Set::Matrix`
from every interior/fine-level node. Interior nodes are ~90% of Fapply time.

**Status:** Implemented bit-identically on the same branch.

**Evidence:** `benchmark/PHASE_C1_fapply_occupancy.md`.

---

### 10.3 Cheap Input Levers (C0): Solver Tolerance Tuning

**Staged change:** `input_3d_centre_bore_256_a2_tuned` sets:
- `elastic.tol_rel/abs` from 1e-8 → 1e-6 (saves ~30–40% V-cycles)
- `bottom_tol_rel=1e-4`, `bottom_max_iter=50` (caps coarse-grid BiCGStab count)
- `elastic.interval` from 50 → 100 (solves elastic half as often)

**Why it should work:** The elastic solve uses 1e-8 tolerance by default — far
tighter than physics coupling requires (the flame update uses the same stress
field at the same timestep, making sub-1e-6 stress accuracy irrelevant to the
output). The V-cycle count scales approximately with `-log10(tol)`. Stacked,
these levers target **3–5× reduction in elastic wall time** with zero source
change.

**Status:** Input file staged; A/B comparison on NOVA pending.

**Evidence:** `benchmark/PHASE_A_FINDINGS.md` §8, Table "Cheap input levers".

---

## 11. CPU-Side Wins That Held the Baseline

### 11.1 Newton Backtracking Line Search (Gated Opt-In)

**Problem:** The elastic Newton solver stalled for low-void-stiffness cases
(void = 0.2 MPa). The linear solve still hit 1e-8 residual, but the nonlinear
iterate overshot and the update was rejected by every subsequent Newton
step — silently exhausting `nriters` and returning zero stress.

**Fix:** Backtracking line search (Armijo-type) in `Newton.H`'s `Set::Vector`
overload, gated behind `elastic.solver.line_search` (default false).

**Why gated:** A global residual-based line search regressed 6 CPU tests
(Eshelby/Fracture/Rubber/TopOp families) because those problems are
non-monotone in the rhs-residual — the search was damping them onto a wrong
path. With `default=false`, the undamped baseline is preserved for all existing
tests; opt-in enables the fix for void-stiffness-sensitive cases.

**Evidence:** Commit `6ebf0ac26`. Memory `cpu_newton_damping.md`.

---

### 11.2 Elastic Solver: psi_floor + 3-Material Model

**Change:** Added `psi_floor` (minimum phase-field value below which void
stiffness is applied) + a 3-material Voigt-rule model (propellant / void / casing).

**Why it worked:** Pure void (psi = 0) means a stiffness ratio of
propellant:void → ∞. MLMG doesn't converge to 1e-8 for infinite contrast ratios.
`psi_floor` caps the ratio at a finite (but still large) value, giving the solver
a convergent problem. The 3-material decomposition (instead of a single blended
material) allows the casing to have physically correct stiffness independently
of the propellant-void blend.

**Evidence:** Commits `38e755e46`, `85681a60a`. Memory `elastic_solver_working.md`.

---

## 12. Strategic Pivots and What They Corrected

### 12.1 D1 Reversed: Device-Elastic Is Viable

**Original decision (2026-06-20):** D1 = elastic stays CPU-resident. GPU diverges
at 2048² even no-fast-math; CPU np8 beats GPU 2.3–3.4× where GPU does converge.

**Reversal (2026-06-22):** The 2048² divergence was a stream race bug, not a
numerical problem. Once fixed, elastic converges on GPU for all tested sizes.
The 2.3–3.4× CPU win was on the local A1000 in 2D — wrong hardware, wrong
regime. Phase A measured **22.9× GPU win** in 3D on A100.

**Lesson:** A "CPU is faster" result that comes from a shared-desktop 50W GPU in
2D does not predict 3D A100 behavior. Disposition decisions require target-hardware
data in the target regime.

---

### 12.2 "Amdahl Cap" Fear Refuted

**Original worry:** If elastic is ~27% of CPU wall time, the combined GPU speedup
is capped at ~3× (since flame-only is already 13×). Elastic on GPU might not be
faster than CPU, making the combined picture lose.

**Correction:** The Amdahl analysis assumed elastic runs at CPU speed on GPU. The
measurement shows elastic is **also** GPU-accelerated, and in the 3D regime the
combined speedup exceeds flame-only. The framing error was treating elastic as a
"fixed-cost" residual rather than an accelerable workload.

---

### 12.3 "Bottom Solver Grinding" Framing Retired

**Original framing:** The elastic MLMG fails because BiCGStab grinds on FP noise
floor residuals (GPU reduction non-associativity amplified catastrophically).
`bottom_tol_abs=1e-3` was proposed as the fix.

**Correction:** This framing was derived from the pre-fix divergence behavior.
After the elixir fix, the solver converges in ~10 V-cycles and the bottom solver
is only 7% of wall time. The "grinding" never occurs in the fixed code. The real
performance problem is the fine-level Fapply register pressure (§9.4), not the
bottom solver.

---

### 12.4 Profiling Must Precede Optimization

**Pattern that failed (Phases 0–5):** Define a "standing metric set" but never
collect it (blocked by `ERR_NVGPUCTRPERM` locally). Make optimization decisions
based on static register dumps and launch-count proxy metrics.

**Pattern that worked (Phase A):** Move profiling to NOVA (where counters are
permitted). Fix the ncu targeting bugs (§8.2). Collect nsys + ncu before
committing to any source change. The result: Phase A's kernel time breakdown
immediately showed that flame is 0.20% of combined time — making the entire
planned "reduce Flame launch overhead" work unnecessary.

**Rule for future work:** The first NOVA run after any major code change should
be a profiling run, not a production run.

---

## Appendix A: Key File Locations

| Artifact | Path |
|---|---|
| Branch guide + build matrix | `benchmark/GPU_BRANCH_GUIDE.md` |
| v2 roadmap (Phases A–D) | `benchmark/GPU_ROADMAP_V2.md` |
| Phase A combined profiling findings | `benchmark/PHASE_A_FINDINGS.md` |
| Phase 3 crossover (flame-only win) | `benchmark/PHASE3_R3_crossover.md` |
| Phase 1 elastic disposition | `benchmark/PHASE1_ELASTIC_DISPOSITION.md` |
| GPU test suite fixes | `benchmark/GPU_TEST_SUITE_FIXES.md` |
| Elixir-race audit | `benchmark/elixir_race_audit.md` |
| G0 kernel resource baseline | `benchmark/G0_BASELINE_OF_RECORD.md` |
| GPU-safe IC/BC matrix | `docs/gpu_safe_ic_bc_matrix.md` |
| Local A100 spoofing guide | `benchmark/LOCAL_A100_SPOOFING.md` |
| ncu empty profile root cause | Memory `gpu_ncu_empty_rootcause.md` |
| Fapply optimization (in progress) | `benchmark/PHASE_C1_fapply_occupancy.md` |
| Newton line search | Commit `6ebf0ac26` |

## Appendix B: The Bug Taxonomy

| Bug | Class | Detected by | Fixed by |
|---|---|---|---|
| Elastic MLMG divergence (2048²) | Stream race / async lifetime | `compute-sanitizer`, multi-box test | `tmpfab.elixir()` in `interpolation()` |
| GPU error-700 traction diagnostic | Host-member write from device kernel | NOVA A100 crash log | `amrex::ReduceOps` |
| 3D NeoHookean CUDA error 719 | Chained Eigen expression on device | NOVA A100 crash log | `.eval()` materialization |
| Restart node-fab OOB segfault | Wrong box type (cell vs node) | C2 test failure | Correct `convert()` call |
| Headerless `thermo.dat` | Missing file header | C2 restart test | Add header write |
| `L_mf` null write (`thermal.on=0`) | Conditional registration, unconditional use | CPU suite SIGSEGV | Register `L_mf` unconditionally |
| `model_prop` parse abort | `query_exactly` arity mismatch | CPU suite abort | Change to `query` or fix input |
| Empty ncu profiles (NVTX) | Wrong kernel name matching mode | Empty output files | `--nvtx-include` |
| Empty ncu profiles (sections) | `--set default` not available in ncu 2025.x | Empty output files | `--set basic` or explicit `--metrics` |
| Silent CPU binary from NOVA build | Stale `.o` from different arch | Binary `file` check | `rm` stale obj dir before build |

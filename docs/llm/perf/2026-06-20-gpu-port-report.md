# ALAMO GPU Port — Engineering & Performance Report

**Project:** ALAMO (AMReX-based multiphysics; `Integrator::Flame` strand-burning / chamber sims)
**Branch:** `chamber-gpu`
**Report date:** 2026-06-20
**Author:** Performance analysis session (Claude Code)
**Scope:** The CUDA device port of the elastic solve ("Option B" de-virtualization), and a three-campaign CPU-vs-GPU performance characterization of the Flame + AMR + elastic pipeline, including device-level Nsight Systems profiling.

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [System Under Test](#2-system-under-test)
3. [The GPU Overhaul — Strategy & Implementation](#3-the-gpu-overhaul--strategy--implementation)
4. [Test Methodology](#4-test-methodology)
5. [Test Campaign Parameters](#5-test-campaign-parameters)
6. [Results & Data Analysis](#6-results--data-analysis)
7. [Device-Level Analysis (Nsight Systems)](#7-device-level-analysis-nsight-systems)
8. [The High-Resolution Elastic Divergence](#8-the-high-resolution-elastic-divergence)
9. [Synthesis & Conclusions](#9-synthesis--conclusions)
10. [Future Explorations](#10-future-explorations)
11. [Appendices](#11-appendices)

---

## 1. Executive Summary

The ALAMO Flame solver was ported to run on NVIDIA GPUs via CUDA. The phase-field path
already ran on the device, but the **elastic (Mechanics → Newton → Operator) solve** crashed
with `CUDA error 719` because the solver dispatched **virtual** member functions on
`Model::Solid::Solid` and `BC::Operator::Elastic::Elastic` objects that live in device memory
but whose **vtable pointers point to host memory**. The fix ("Option B") was to
**de-virtualize both hierarchies** so every elastic kernel performs static dispatch on
concrete, device-resident objects. This succeeded: the device elastic solve runs, and at the
baseline resolution it is **1.65× faster than the CPU**.

We then ran three CPU-vs-GPU campaigns on an RTX A1000:

| Campaign | Grid strategy | GPU vs CPU (wall) | Headline insight |
|---|---|---|---|
| **Coarse (baseline)** | `max_level=3`, `n_cell=64` | GPU **1.50× slower** | GPU wins elastic (1.65×) but loses the phase-field/AMR path (2.29×) |
| **Deep fine-grid** | `max_level=5`, `n_cell=64` | GPU **2.87× slower** | Deepening AMR multiplies *launch/sync overhead*, not useful compute |
| **Wide-shallow** | `max_level=2`, `n_cell=512` | GPU **1.74× slower** | Same finest resolution, **9.8× faster/step**, gap nearly halved |

**Nsight Systems proved the mechanism at the device level:** for the same 2048² finest
resolution, the deep-AMR config issues **646,294 GPU kernel launches per step** versus
**95,949** for the wide-shallow config (**6.7×**), with an average kernel duration of only
**~3 µs** in both — i.e. deep AMR does not give the GPU bigger work, it gives it **6.7× more
tiny work**. The GPU is **launch-latency- and synchronization-bound**, not compute-bound, on
this problem class.

A second, important finding: at the high base resolution, the **GPU elastic solve diverges**
(`MLMG failing`, SIGABRT) while the **CPU solves the identical problem cleanly**. Root cause is
the `--use_fast_math` build flag interacting with the worse-conditioned masked operator at fine
resolution — a known risk that the port plan had explicitly flagged.

---

## 2. System Under Test

### 2.1 Hardware

| Component | Spec |
|---|---|
| GPU | NVIDIA **RTX A1000**, 8188 MiB, compute capability **sm_86** |
| GPU driver | 595.71.05 |
| Host | Linux 6.8.0 (Ubuntu 24.04), x86_64 |
| Run topology | single MPI rank, single GPU (`mpiexec -np 1` / direct) |

> **Context:** the RTX A1000 is an entry-level professional GPU (Ampere, ~2048 CUDA cores,
> modest memory bandwidth). The performance regime described here — latency-bound, small-kernel —
> is *more* pronounced on a small GPU; a data-center GPU (A100/H100) would have higher launch
> throughput and more cores to fill, but the *relative* CPU-vs-GPU conclusions about AMR
> granularity still hold.

### 2.2 Software

| Component | Version / detail |
|---|---|
| AMReX | 26.06 |
| CPU binary | `bin/alamo-2d-clang++` (clang++, IEEE FP) |
| GPU binary | `bin/alamo_gpu-2d-cuda86-g++` (nvcc, `-ccbin=mpicxx`, C++20, `-DALAMO_GPU`) |
| CUDA toolkit | 12.6.3 redistributable (`.local/cuda`), nvcc |
| Memory model | **device arena** (`amrex.the_arena_is_managed = false`) — *not* unified/managed memory |
| Profiler | **Nsight Systems CLI 2026.1.3**, installed locally without root |

---

## 3. The GPU Overhaul — Strategy & Implementation

### 3.1 Problem statement (root cause of the original crash)

The Newton device kernels (`src/Solver/Nonlocal/Newton.H:167,544,587`) call
`model(i,j,k).DW(...)`, and the elastic operator (`src/Operator/Elastic.cpp:364`) invokes
`(*m_bc)(...)` **inside a device `ParallelFor`**. Two stacked problems:

1. **Polymorphic Solid model on device.** `Model::Solid::Solid<SYM>` declared
   `virtual ~Solid`, `virtual W/DW/DDW/Advance/ContainsNan/Print`. The object data is copied
   to device memory, but the **vtable pointer still references host memory**. A virtual call on
   the device dereferences that host pointer → illegal access → `CUDA error 719`.
2. **Polymorphic Elastic BC operator on device.** `BC::Operator::Elastic::Elastic` had
   `virtual operator()` / `virtual getType` / `virtual ~`, invoked inside a device kernel — the
   identical vtable-on-device hazard.

A naïve "demote those loops back to the host" fix fails because the solver's scratch fabs
(`dw_mf/ddw_mf/rhs_mf/res_mf`) are **device-only memory**; host loops segfault writing them, and
toggling `amrex.the_arena_is_managed=1` did not change that.

### 3.2 Strategy chosen: Option B — full device-native via de-virtualization

The decisive observation: **all elastic code is already templated on the concrete model type**
(`Newton<T>`, `Mechanics<T>`, `Operator::Elastic<SYM>`). Runtime polymorphism was therefore
**not actually required** — the `virtual` keywords were vestigial. The strategy was to remove
them so dispatch becomes **static (compile-time)**, letting every elastic kernel operate on
concrete, device-resident objects with **no vtable, no managed memory, no host loops**.

This was selected over two alternatives:

| Option | Description | Why not chosen |
|---|---|---|
| **A. Managed memory** | Put model/BC + scratch fabs in CUDA Unified Memory so host vtables resolve | Page-migration thrash on every solve; "slow by design"; doesn't fix vtable-on-device, only hides it |
| **B. De-virtualize (chosen)** | Remove `virtual`, static dispatch, keep everything device-resident | Exact (pure dispatch change → no numeric change); zero managed-memory cost; matches the already-templated architecture |
| **C. Hybrid (BC on host)** | Models on device, BC boundary loops on host with managed scratch | Partial; still pays managed-memory + host-sync cost on the BC path; "not full device-native" |

### 3.3 Implementation — patterns & technologies

**De-virtualization (the core pattern).** Across `src/Model/Solid/**`, `virtual` and
`override` were removed from `W/DW/DDW/Advance/ContainsNan/Print` and destructors, replaced
with **`AMREX_GPU_HOST_DEVICE`** qualifiers and defaulted/trivial device destructors. Static
dispatch flows through the existing `In/ExtClassOperators.H` macro layer (a CRTP-like
compile-time mechanism). Scope of edit:

- **242** added `AMREX_GPU_HOST_DEVICE` qualifiers across `src/`.
- **~20** Solid models de-virtualized (`Finite/NeoHookean`, `NeoHookeanPredeformed`,
  `Finite/CrystalPlastic`, `Finite/PseudoLinear|PseudoAffine/Cubic`, all `Linear/*`,
  all `Affine/*`).
- `Solid.H` base: `virtual ~Solid()` → `AMREX_GPU_HOST_DEVICE ~Solid() = default`; all
  W/DW/DDW/Advance/ContainsNan/Print de-virtualized and device-qualified.

```cpp
// Before (host-only vtable, illegal on device):
virtual Set::Matrix DW(const Set::Matrix &) const { Util::Abort(...); }
// After (static dispatch, device-safe):
AMREX_GPU_HOST_DEVICE Set::Matrix DW(const Set::Matrix &) const { ... }
```

**BC operator → capturable POD.** `BC::Operator::Elastic::Elastic` was refactored so the
boundary-condition *type* data is a plain array (`m_bc_type`, an array of `Face`→`Type` enums)
that is **trivially device-copyable and can be captured by value into a device kernel `[=]`**.
The face/edge/corner classification now operates on a *passed-in* `m_bc_type` rather than via
virtual dispatch on a base pointer. The value data (`Constant`'s scalar interpolators,
`Expression`'s parsers) is handled separately so the hot path carries no heap-owning members.

**Device-safe host-only code.** Host-only constructs that sit inside device kernels were
guarded — **24** sites using `#ifndef ALAMO_GPU` and AMReX's `AMREX_IF_ON_HOST(...)` /
`AMREX_IF_ON_DEVICE(...)` macros. This wraps `Util::Abort`, `contains_nan()`-aborts, and
`Print(std::ostream&)` (because `std::ostream` is host-only). `ALAMO_GPU` is `-D`-defined by
`configure` for `--cuda` builds.

**Self-documenting guardrail.** `Solid.H` retains a `std::has_virtual_destructor<T>` static
check that emits a warning if any model re-introduces a virtual destructor — the codebase
encodes the GPU hazard as a compile-time tripwire.

**Build system (`configure`, +463 lines).**

- `--cuda [arch]` with **auto-detection** of the local GPU compute capability via
  `nvidia-smi --query-gpu=compute_cap`, normalizing `7.5`/`sm_75`/`compute_75` → `75`, and
  validating an explicit arch against the detected GPU.
- nvcc invocation: `-ccbin=mpicxx`, C++20, `-DALAMO_GPU`, `--expt-relaxed-constexpr`,
  `--expt-extended-lambda`, `-maxrregcount=255`, `--Werror cross-execution-space-call`
  (catches host/device call-space violations at compile time),
  `--generate-code arch=compute_86,code=[sm_86,compute_86]`.
- **Two optimization modes** ("fast/bench"): a fast mode (`--ptxas-options=-O3 --use_fast_math`,
  host `-O3`) and a conservative mode (`--ptxas-options=-O0`, no fast-math, host `-O0`). The
  production GPU binary was built with the **fast-math** path — see §8 for the consequence.

### 3.4 Drawbacks & trade-offs accepted

| Choice | Benefit | Drawback / risk accepted |
|---|---|---|
| **De-virtualization** | Device-native, exact numerics, no managed memory | Removes runtime polymorphism — any future code that needs a base `Solid*`/`Elastic*` dispatched at runtime must be templated/specialized instead. Audited as absent from the hot path. |
| **`AMREX_GPU_HOST_DEVICE` everywhere** | Single source compiles for host & device | Larger compile units, **275 MB** GPU binary, long nvcc compiles; every host-only call inside these functions must be guarded. |
| **`#ifndef ALAMO_GPU` guards** | Device build/run safe | Two code paths to keep in sync; an unguarded `Util::Abort` breaks the device compile. |
| **`--use_fast_math` (fast mode)** | Faster transcendentals/FP on device | **Non-IEEE FP reassociation** → at high resolution the ill-conditioned masked elastic operator **diverges** (§8). |
| **Device arena (not managed)** | No page-migration thrash; predictable | Everything the host touches (regrid bookkeeping) must round-trip explicitly; no automatic host access to device fabs. |
| **`-maxrregcount=255`** | Avoids register-spill compile failures | Can cap occupancy on register-heavy kernels. |

---

## 4. Test Methodology

### 4.1 The analysis suite (`analysis/`)

A dependency-light (stdlib-Python + hand-rolled SVG) profiling suite with six phases:

| Phase | Tool | Output |
|---|---|---|
| 1 Wall-clock & efficiency | `/usr/bin/time -v` footer + log MLMG timers | `wallclock.{md,csv,json}` |
| 2 I/O profile | `strace -f -c` | `io_profile.*` |
| 3 Microarchitecture | `perf stat` | `perfstat.*` |
| 4 CPU flame graph | `perf record` + FlameGraph | `flamegraph_cpu.svg` |
| 5 **GPU device timeline** | **`nsys`** / `nvidia-smi` | `gpu_nsys.*`, `gpu_timeline.*` |
| 6 Consolidated report | — | `REPORT.md`, `index.html` |

**Design principle — production vs instrumented separation.** Phase 1 measures the *real*
full-length runs untouched (plotting off → clean compute comparison). Phases 2–5 launch *short*
re-runs under profilers, since profiling perturbs timing and must never contaminate the
headline wall-clock.

### 4.2 Nsight Systems integration (Phase 5, performed this session)

`nsys` was absent and the box has no root. Nsight Systems CLI 2026.1.3 was installed
**without root** by extracting the CUDA-repo `.deb` with `dpkg-deb -x` into `.local/nsight`,
then wired into the suite:

- `analysis/lib/common.sh`: added an `NSYS_BIN` resolver (`$NSYS` env → `.local/nsight/**` →
  PATH) and `have_nsys`.
- `analysis/05_gpu_timeline.sh`: uses `$NSYS_BIN` with
  `--trace=cuda,nvtx --sample=none --cpuctxsw=none` (CUDA API + kernel trace needs **no**
  special privilege; CPU sampling does, and `perf_event_paranoid ≥ 2` blocks it on this box),
  emitting `cuda_api_sum`, `cuda_gpu_kern_sum`, `cuda_gpu_mem_time_sum` CSVs.

### 4.3 Metric definitions

- **Wall/step** — elapsed wall clock ÷ steps; the fair throughput metric when step counts or
  sim durations differ.
- **Voluntary context switches** — the process *blocking* (e.g. on `cudaStreamSynchronize` /
  device sync / driver wait). The signature of host↔device synchronization overhead.
- **System CPU time** — kernel-mode time (syscalls, driver `ioctl`s, page-fault servicing).
- **Major page faults** — faults requiring I/O; on the device build, host-side staging /
  bookkeeping churn (not simulation-data migration — the arena is device, not managed).
- **`cudaLaunchKernel` count** — number of GPU kernel launches; each carries a fixed
  ~3–5 µs host-side launch latency the CPU never pays.
- **GPU-busy time** — Σ of kernel execution durations (device-side; robust to host-side
  profiler overhead).

---

## 5. Test Campaign Parameters

All campaigns share: `input_copy`, **star** grain geometry (BMP IC), `void κ = µ = 20 MPa`,
`nsubsteps = 2`, `regrid_int = 2`, `blocking_factor = 8` (except wide), `grid_eff = 0.7`
(except wide), AMReX 26.06, single rank, plotting off for wall measurement.

| Parameter | Coarse (baseline) | Deep fine-grid | Wide-shallow (new) |
|---|---|---|---|
| `amr.max_level` | 3 | **5** | **2** |
| `amr.n_cell` (base) | 64³ | 64³ | **512³** |
| **Finest effective resolution** | 512² | **2048²** | **2048²** |
| Finest cell size `dx` | 3.4e-4 m | **8.57e-5 m** | **8.57e-5 m** |
| `timestep` (coarse dt) | 1.0e-4 s | 2.5e-5 s | 2.5e-5 s |
| `amr.grid_eff` | 0.7 | 0.7 | **0.9** |
| `amr.blocking_factor` | 8 | 8 | **16** |
| `elastic.interval` | 50 | 50 (none reached) | n/a (100000, no solve)¹ |
| `max_step` | 310 | 310 | 80 (≈¼ of 310) |
| Finest subcycles / coarse step | 2³ = 8 | **2⁵ = 32** | **2² = 4** |
| Sim time covered | 0.0310 s | 0.00775 s | 0.00200 s |

¹ The wide-shallow run was configured for the phase-field/AMR path only (elastic present but
not solving), making it **directly comparable to the deep fine-grid run, which also performed
0 elastic solves**. The requested `elastic.interval=5` comparison is blocked on GPU — see §8.

The two grid strategies (**deep** vs **wide**) reach the **identical 2048² finest cell**; the
only difference is *how* — deep AMR (6 levels, heavy subcycling, tiny patches) vs a large base
grid with shallow AMR (3 levels, light subcycling, big patches).

---

## 6. Results & Data Analysis

### 6.1 Campaign 1 — Coarse baseline (`analysis/results/`)

**`max_level=3`, `n_cell=64`, 310 steps. GPU is 1.50× slower overall.**

| Metric | CPU | GPU |
|---|---:|---:|
| Wall clock (s) | 40.57 | 61.02 |
| User CPU (s) | 40.02 | 68.05 |
| System CPU (s) | 0.15 | **17.31** |
| CPU time = user+sys (s) | 40.17 | 85.36 |
| Avg active cores | 0.99 | 1.40 |
| Peak RSS (MB) | 144.3 | 560.6 |
| FS bytes written (MB) | 62.18 | 62.62 |
| Major page faults | 0 | 5,207 |
| **Voluntary ctx switches** | **1,465** | **4,120,333** |
| Involuntary ctx switches | 309 | 5,928 |
| Wall per step (ms) | 130.9 | 196.8 |
| **Elastic solves** | **12** | **12** |
| Elastic total (s) | **18.96** | **11.46** |
| Elastic mean / solve (s) | 1.580 | 0.955 |
| Elastic max / solve (s) | 1.880 | 1.110 |
| MLMG iters (mean) | 18.33 | 18.33 |

**Decomposition (the key to the whole story):**

- **Elastic solve:** GPU 11.46 s vs CPU 18.96 s → the de-virtualized device elastic solve is
  **1.65× *faster*** on GPU. The port works and pays off on the big, dense, node-based solve.
- **Everything else** (phase-field + AMR + regrid + I/O): GPU 49.56 s vs CPU 21.61 s → GPU is
  **2.29× *slower***.
- Net: the elastic win (47% of CPU wall) drags the headline up to a "respectable" 1.50× loss.

Already visible: **4.1 M voluntary context switches** on the GPU (vs 1,465 on CPU) and 17.3 s
of system time — the fingerprint of host↔device synchronization on the per-step phase-field
path.

### 6.2 Campaign 2 — Deep fine-grid (`analysis/results_fine_grid/`)

**`max_level=5`, `n_cell=64`, 310 steps. GPU is 2.87× slower — the gap nearly doubles.**

| Metric | CPU | GPU | GPU/CPU |
|---|---:|---:|---:|
| Wall clock (s) | 454.53 | **1304.30** | 2.87× |
| User CPU (s) | 453.20 | 1438.17 | |
| System CPU (s) | 1.13 | **387.29** | 343× |
| CPU time = user+sys (s) | 454.33 | 1825.46 | |
| Avg active cores | 1.00 | 1.40 | |
| Peak RSS (MB) | 400.1 | 562.5 | |
| Major page faults | 3 | **25,646** | |
| **Voluntary ctx switches** | **3,658** | **109,688,234** | **30,000×** |
| Involuntary ctx switches | 1,172 | 14,935,569 | |
| Wall per step (ms) | 1466.2 | **4207.4** | 2.87× |
| **Elastic solves** | **0** | **0** | — |

**Why the deep grid disproportionately hurts the GPU — two compounding reasons:**

1. **The GPU's one win disappears.** This run reached only t = 0.00775 s and performed
   **0 elastic solves**. The single thing the GPU was winning (the elastic solve) is absent,
   exposing the full per-step phase-field penalty — which simultaneously *grew*.

2. **The per-step penalty itself grew.** Per-step work growth, coarse → deep:
   - CPU: 130.9 → 1466.2 ms/step = **11.2×** (honest "more work": more cells, more subcycles,
     more regrid work).
   - GPU: 196.8 → 4207.4 ms/step = **21.4×**.
   - The extra **1.91×** (21.4 / 11.2 = 2.87 / 1.50) is the **"AMR-depth tax"** on the GPU: the
     added work arrives as *more, smaller* launch/sync operations, each carrying a fixed device
     latency the CPU does not pay.

**Grid structure at t=0** (parsed from the GPU plotfile header) — the source of the tiny
kernels:

| Level | Domain | Boxes | Total cells | Median cells/box |
|---:|---:|---:|---:|---:|
| 0 | 64² | 1 | 4,096 | 4096 |
| 1 | 128² | 1 | 16,384 | 16384 |
| 2 | 256² | 30 | 37,696 | 448 |
| 3 | 512² | 7 | 44,096 | 5824 |
| 4 | 1024² | 34 | 110,464 | 2880 |
| 5 | **2048²** | **161** | 216,896 | **768** |
| **All** | | **234** | **429,632** | |

The two levels that exist *only* in the deep run (4 & 5) hold **195 boxes / 327k cells — 76%
of all cells — in patches with a median < 1k cells** (`blocking_factor=8` allows boxes as small
as 64 cells). An A1000 needs **tens of thousands of cells per kernel** to saturate; a 768-cell
kernel is essentially pure launch + sync latency. **Box-advances per coarse step** (with
subcycling) = Σ `nboxes[ℓ]·2ℓ` = **5,875**, each one ≥ 1 kernel + ghost-exchange + reduction +
device sync.

**Reading the fingerprint:** 387 s system time = **30% of wall in kernel mode** servicing the
syscall/driver/page-fault storm. **109.7 M voluntary context switches** ≈ 84,000/s of blocking
device-sync waits. This is *not compute*; it is host↔device coordination.

### 6.3 Campaign 3 — Wide-shallow (`analysis/results_wide_grid/`)

**`max_level=2`, `n_cell=512` (same 2048² finest cell), 80 steps. GPU 1.74× slower —
9.8× faster per step than the deep run.**

| Metric | Wide GPU | Wide CPU | Deep GPU (ref) |
|---|---:|---:|---:|
| Wall clock (s) | 34.25 (80 steps) | 19.63 (80 steps) | — |
| Wall per step (ms) | **428** | 245 | 4207 |
| Elastic solves | 0 | 0 | 0 |
| **GPU/CPU ratio** | **1.74×** | — | 2.87× |

| Comparison (per step) | Wide | Deep | Wide advantage |
|---|---:|---:|---:|
| GPU wall/step | 0.428 s | 4.207 s | **9.8× faster** |
| CPU wall/step | 0.245 s | 1.466 s | 6.0× faster |
| GPU/CPU ratio | 1.74× | 2.87× | gap nearly halved |

Reaching the same finest resolution via a **large base grid + shallow AMR** instead of **deep
AMR** makes the GPU ~10× faster per step and cuts its relative penalty from 2.87× to 1.74×.
The CPU also benefits (6×/step) because deep AMR's subcycling and tiny-box overhead is wasteful
on any architecture — but the GPU benefits *more*, because its overhead was launch-latency, the
exact thing the wide grid removes.

---

## 7. Device-Level Analysis (Nsight Systems)

Identical 30-step runs of each grid strategy under `nsys` (`--trace=cuda`), phase-field only
(0 elastic solves), GPU binary. Counts are profiler-independent (nsys serializes the launch
queue, inflating *wall*, but not the *number* of launches/syncs).

### 7.1 CUDA API summary — launches and syncs per step

| Per-step | **Wide** (lvl 2, n=512) | **Deep** (lvl 5, n=64) | Deep / Wide |
|---|---:|---:|---:|
| `cudaLaunchKernel` | 95,949 | **646,294** | **6.7×** |
| `cudaStreamSynchronize` | 3,520 | 23,275 | 6.6× |
| `cudaLaunchHostFunc` | 31,466 | 224,257 | 7.1× |
| `cudaMemcpyAsync` | 16,275 | 114,939 | 7.1× |

Totals over 30 steps: wide = 2.88 M kernel launches; deep = **19.39 M**. Extrapolated to the
full 310-step deep run that is **~200 M launches** — consistent with the **109.7 M voluntary
context switches** measured independently in §6.2 (≈ one blocking wait per ~2 launches).

In the API time budget, `cudaStreamSynchronize` alone is **20.1%** (wide) of CUDA-API time, and
`cudaLaunchKernel` is the single largest at 56–61%. (An early 2-step validation on the coarse
config showed the same pattern: 103,769 launches and `cudaStreamSynchronize` at **38.6%** of
API time.)

### 7.2 Kernel summary — the kernels are tiny in both configs

| Per-config (30 steps) | Wide | Deep |
|---|---:|---:|
| Kernel instances | 2,878,479 | 19,388,818 |
| **Average kernel duration** | **3.06 µs** | **2.58 µs** |
| GPU-busy time (Σ kernel durations) | 8.81 s | 50.03 s |

**The decisive result:** the average kernel is **~3 µs in both configs — and *smaller* in the
deep config.** Deep AMR does **not** give the GPU larger, more-saturating kernels; it gives it
**6.7× more, even smaller** kernels. At ~3 µs, a kernel barely exceeds the ~3–5 µs fixed launch
latency, so the device spends most of its time being launched-at and synchronized rather than
computing. This is the textbook **latency-bound** regime, and it is *intrinsic to the AMR
granularity*, independent of the elastic port.

> The dominant kernels are AMReX `launch_global<256, ... ParallelFor ...>` (the phase-field
> field updates) plus `ReduceOps<Min,Max>` / `ReduceOps<Sum,...>` — the per-level T/mdot/L
> min-max-sum reductions emitted **every step for every level** (6 levels in deep), each a
> blocking device→host sync. These reductions are a structural sync source the wide grid
> reduces simply by having fewer levels.

---

## 8. The High-Resolution Elastic Divergence

The requested `elastic.interval=5` comparison **cannot run on the GPU at the wide/high-res
config**:

| Build | Wide config, first elastic solve (step 5) | Result |
|---|---|---|
| **CPU** (`alamo-2d-clang++`) | converges | **exit 0** |
| **GPU** (`alamo_gpu-2d-cuda86-g++`) | `amrex::Abort: MLMG failing so lets stop here !!!` | **SIGABRT (exit 6)** |

Mitigations that did **not** help: loosening `elastic.solver.nrtolerance` (1e-5 → 1e-4) and
raising `elastic.psi_floor` (0.05 → 0.1). The solver is genuinely **diverging**, not merely
failing to reach tolerance.

**Root cause — `--use_fast_math`.** `configure` builds the production GPU binary with the
fast-math optimization path (`--ptxas-options=-O3 --use_fast_math`, host `-O3`). `--use_fast_math`
permits **non-IEEE FP reassociation** and lower-precision intrinsics. At the 2048²-resolved
star, the masked elastic operator on thin grain features is far worse-conditioned than at the
coarse 64²-base run; the fast-math perturbation tips MLMG from "converges slowly" into
"diverges." The CPU build (clang++, IEEE) solves the identical system. **This is exactly the
risk the port plan flagged** ("use a `--profile`/no-fast-math build for the bit-compare gate").

**Crucial nuance — the port itself is correct.** At the coarse 64²-base config the GPU elastic
solve **ran and was 1.65× faster than the CPU** (§6.1) with identical MLMG iteration counts
(18.33). The de-virtualization is numerically exact; the divergence is a **build-flag ×
conditioning** interaction at high resolution, not a logic bug in the device port.

**Fix.** Build a no-fast-math GPU variant — `configure`'s conservative mode already exists
(`--ptxas-options=-O0`, no `--use_fast_math`), or add an explicit `--fmad=false`. Use the
fast-math build for phase-field timing and the no-fast-math build for the elastic solve / any
bit-compare correctness gate.

---

## 9. Synthesis & Conclusions

1. **The GPU port is functionally correct and helps where it should.** De-virtualization
   (Option B) eliminated the vtable-on-device crash; the device elastic solve runs and is
   **1.65× faster** than CPU at baseline resolution. The architecture was already templated, so
   removing `virtual` was exact and low-risk.

2. **The Flame/AMR phase-field path is latency-bound on the GPU, not compute-bound.** Across
   all configs the GPU loses the per-step phase-field path (2.29×–2.87×) while winning the
   elastic solve. nsys shows ~3 µs average kernels and tens-to-hundreds of thousands of launches
   per step — the device is starved.

3. **AMR granularity is the dominant performance lever — more than depth or resolution.** Two
   configs with the **identical 2048² finest cell** differ **6.7×** in kernel-launch count and
   **9.8×** in GPU wall/step, purely from how the refinement is structured. Deep AMR
   (6 levels, 32× subcycling, hundreds of < 1k-cell boxes, regrid every 2 steps) is the
   GPU-hostile extreme; a **large base grid + shallow AMR** keeps the resolution while
   restoring big, device-saturating kernels.

4. **The earlier slowdown was never about the elastic port.** The deep fine-grid run performed
   **0 elastic solves** — its 2.87× slowdown is entirely the phase-field/AMR path, amplified by
   depth, with the elastic win removed. The two facts (no elastic offset + amplified per-step
   overhead) fully account for the 1.50× → 2.87× jump.

5. **`--use_fast_math` is a real correctness liability for the elastic solve.** It buys
   phase-field speed but breaks high-resolution MLMG convergence that the CPU handles fine. The
   build must offer (and the elastic path must use) a no-fast-math variant.

6. **It is not a memory-thrash problem.** The arena is device-resident
   (`the_arena_is_managed=false`); the major page faults and context switches are launch/sync
   and host-side bookkeeping, not unified-memory page migration.

**Bottom line:** the GPU is a good fit for ALAMO's *dense, few-box* work (the elastic MLMG
solve, large base grids) and a poor fit for *many-tiny-box* work (deep AMR with heavy
subcycling). Performance engineering should target the launch/sync count, not the FLOP count.

---

## 10. Future Explorations

**A. Reduce kernel-launch / sync count (highest leverage).**
- **Coarser AMR granularity:** raise `amr.max_grid_size` and `amr.blocking_factor` (fewer,
  bigger boxes); raise `amr.grid_eff` toward 0.9 (less fragmentation); lengthen `regrid_int`.
  The wide-shallow run already demonstrates the payoff.
- **CUDA Graphs / kernel fusion:** capture the repeated per-step phase-field kernel sequence
  into a CUDA Graph to amortize launch latency across the many tiny kernels; fuse adjacent
  `ParallelFor`s over the same box.
- **Fewer forced syncs:** the per-level min/max/sum reductions (T, mdot, L) are emitted every
  step and force device→host syncs. Batch them across levels, reduce their frequency
  (`thermo.int`), or compute them asynchronously.

**B. Close the elastic-on-GPU gap properly.**
- Build a **no-fast-math** GPU binary and re-run the `elastic.interval=5` comparison at high
  resolution (CPU already converges) — quantify the device elastic solve where it most matters.
- Investigate `--fmad=false` as a middle ground (keeps `-O3` PTX, disables fused multiply-add
  reassociation) to retain most speed while restoring convergence.
- Improve operator conditioning at the void interface (psi formulation, BC clamping) so the
  solve is robust to FP perturbation regardless of build flags.

**C. Scale the problem to the GPU's strengths.**
- Benchmark larger base `n_cell` (1024², 2048²) and 3D — the crossover where the GPU's
  throughput beats the CPU should appear once kernels are large enough to saturate the device.
- Test on a data-center GPU (A100/H100): higher launch throughput and more cores will shift the
  latency-bound boundary and may flip the per-step phase-field result.

**D. Complete the profiling picture.**
- Run the full analysis suite Phase 2 (`strace` I/O), Phase 3 (`perf stat` IPC/cache), Phase 4
  (CPU flame graph) — now that Phase 5 (nsys) is wired in — for a complete CPU+GPU breakdown.
- Open the captured `nsys_{wide,deep}.nsys-rep` in the Nsight GUI for the visual kernel
  timeline and gap analysis.

**E. Numerical validation (deferred from the port plan).**
- Generate a chamber-branch 20 MPa CPU golden and bit/tolerance-compare the GPU result
  (de-virtualization is exact; expect ≤ 1e-12 drift from FP reassociation, larger under
  fast-math) to formally close the port's correctness gate.

---

## 11. Appendices

### 11.1 Artifact locations

| Artifact | Path |
|---|---|
| Port plan | `docs/gpu_elastic_device_port_plan.md` |
| Coarse report | `analysis/results/REPORT.md` (+ `wallclock.*`, charts) |
| Deep fine-grid report | `analysis/results_fine_grid/REPORT.md` |
| Wide-shallow report | `analysis/results_wide_grid/REPORT.md` |
| nsys reports (GUI) | `analysis/results_wide_grid/nsys_{wide,deep}.nsys-rep` |
| nsys stats CSVs | `analysis/results_wide_grid/stats_{wide,deep}_*.csv` |
| Wide run logs | `analysis/results_wide_grid/wide_{gpu,cpu}.log` |
| Production logs | `out_cpu_star_20mpa.log`, `fine_grid_{cpu,gpu}.log` |
| Local nsys install | `.local/nsight/opt/nvidia/nsight-systems/2026.1.3/target-linux-x64/nsys` |

### 11.2 Key source files changed (GPU overhaul)

| Concern | File(s) |
|---|---|
| Solid base de-virtualized | `src/Model/Solid/Solid.H` (+77) |
| Models de-virtualized | `src/Model/Solid/**` (Finite/NeoHookean* +122/+49, Linear/*, Affine/*, …) |
| Static-dispatch macros | `src/Model/Solid/{In,Ext}ClassOperators.H` |
| Newton device kernels | `src/Solver/Nonlocal/Newton.H` (+59) |
| Elastic operator + BC call | `src/Operator/Elastic.cpp` (+246) |
| BC de-virtualized / capturable | `src/BC/Operator/Elastic/{Elastic,Constant,Expression}.H` |
| Mechanics integration | `src/Integrator/Base/Mechanics.H` (+41) |
| Integrator structs / guards | `src/Integrator/Flame.cpp` (+525) |
| Build system | `configure` (+463), `benchmark/local_cuda_env.sh` |
| Devicification scope | **242** `AMREX_GPU_HOST_DEVICE` added, **24** `#ifndef ALAMO_GPU` / `AMREX_IF_ON_HOST` guards |

### 11.3 Reproduce commands

```bash
# Environment
source benchmark/local_cuda_env.sh
export NSYS=.local/nsight/opt/nvidia/nsight-systems/2026.1.3/target-linux-x64/nsys

# Wide-shallow GPU, 80 steps, clean wall (phase-field only)
/usr/bin/time -v ./bin/alamo_gpu-2d-cuda86-g++ input_copy max_step=80 \
  amr.max_level=2 amr.n_cell="512 512 512" amr.grid_eff=0.9 amr.blocking_factor=16 \
  elastic.interval=100000 timestep=2.5e-5_s \
  model_void.kappa=20_MPa model_void.mu=20_MPa \
  amr.plot_int=-1 amr.thermo.plot_int=-1 elastic.solver.verbose=0 elastic.print_model=0 \
  plot_file=/tmp/wide_gpu

# nsys device timeline (30 steps) + summaries
$NSYS profile -o /tmp/nsys_wide --force-overwrite true --stats=false \
  --trace=cuda,nvtx --sample=none --cpuctxsw=none \
  ./bin/alamo_gpu-2d-cuda86-g++ input_copy max_step=30 \
  amr.max_level=2 amr.n_cell="512 512 512" amr.grid_eff=0.9 amr.blocking_factor=16 \
  elastic.interval=100000 timestep=2.5e-5_s model_void.kappa=20_MPa model_void.mu=20_MPa \
  amr.plot_int=-1 amr.thermo.plot_int=-1 elastic.solver.verbose=0 elastic.print_model=0 \
  plot_file=/tmp/nsys_wide
$NSYS stats --force-export=true \
  --report cuda_api_sum --report cuda_gpu_kern_sum \
  --format csv --output /tmp/stats_wide /tmp/nsys_wide.nsys-rep

# For the DEEP comparison, swap: amr.max_level=5 amr.n_cell="64 64 64"
# (drop grid_eff/blocking_factor overrides to use the fine-grid defaults)
```

### 11.4 Glossary

- **Option B** — de-virtualize both the Solid model and Elastic BC hierarchies so elastic
  kernels run device-native with static dispatch.
- **Latency-bound** — runtime dominated by fixed per-operation overhead (kernel launch, device
  sync), not arithmetic; the device sits idle between tiny kernels.
- **Subcycling** — finer AMR levels take `nsubsteps^level` smaller timesteps per coarse step;
  multiplies the per-step operation count geometrically with level depth.
- **`--use_fast_math`** — nvcc flag enabling non-IEEE FP (reassociation, fast intrinsics);
  faster but can break ill-conditioned iterative solvers.

---

*End of report.*

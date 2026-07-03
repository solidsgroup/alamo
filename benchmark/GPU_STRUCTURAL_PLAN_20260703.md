# GPU Structural Performance Plan — elastic, phase field, multi-GPU (2026-07-03)

> **What this is.** A synthesis investigation (2026-07-03) into the *large structural
> code changes* that would speed up the GPU implementation, covering the elastic
> solver, the phase-field (Flame) regression path, and multi-GPU. It consolidates
> the measured evidence in the almanac, adds several new code-level discoveries and
> one new local measurement, and lays out a ranked structural plan.
>
> **What this is not.** It does not supersede `benchmark/GPU_ROADMAP_V3.md` — it
> *feeds* it. Every proposal below maps to a v3 task ID where one exists, and new
> proposals are labeled NEW. The v3 governing rule stands: **no optimization ships
> without a physics-error-budget pass** (`benchmark/validate/`, READ_FIRST §8a).
>
> Written by an autonomous investigation session; four parallel sub-agent
> deep-reads (almanac digest, elastic code path, flame code path, multi-GPU
> evidence) plus direct verification of every load-bearing claim. All file:line
> cites verified against `chamber-gpu` @ `dd4054056`.

---

## 1. Where the time goes (settled evidence, cited not re-argued)

| Fact | Value | Source |
|---|---|---|
| Elastic MLMG share of combined evolve (A100, 256³) | ≈95.7% | `archive/PHASE_A_FINDINGS.md` §1 |
| `Fapply` share of GPU kernel time | 74.6% | PHASE_A §2 |
| `Fapply` registers/thread (runtime ncu, A100) | 255 (arch cap) → ~12.5% occupancy | PHASE_A §4 |
| Flame share of combined evolve | 0.20% | PHASE_A §1 |
| Combined speedup vs 64-rank CPU (inflated baseline) | 22.9× | PHASE_A §6; fair baseline owed (v3 task 2.A) |
| Flame-only GPU speedup, saturating 3D | 9.6× @128³ / 13.1× @256³ | `archive/PHASE3_R3_crossover.md` |
| Multi-GPU 2×A100 | 3.5–13.3× **slower** than 1 GPU | PHASE3_R3_crossover |
| Per solve: outer V-cycles / Fapply calls | ~10 / ~3,504 (fine applies ≈90% of Fapply time) | PHASE_A §3 |
| Fapply arithmetic intensity (source estimate) | ≈0.2 flop/byte — memory-bound *if* occupancy were healthy | `GPU_AUDIT_20260702.md` §5 |

Key structural context: the solve is **matrix-free**. Every `Fapply` call re-derives
the operator action from the stored 4th-order tangent field `m_ddw_mf`
(`Set::Matrix4<3,Sym::Major>` = 45 doubles = **360 B per node**), and `Fapply` is
called O(100)× per Newton iteration (2 Jacobi sub-iterations × pre/post smooth ×
levels × ~10 V-cycles + residuals + bottom solver).

---

## 2. New discoveries from this investigation (2026-07-03)

### 2.1 The 255-register fact is now confirmed *statically and locally* (new measurement)

`cuobjdump --dump-resource-usage` on the current local 3D binary
(`bin/alamo_gpu-3d-cuda86-g++`, sm_86):

```
Operator::Elastic<1>::Fapply   main ParallelFor:  REG:255  STACK:192
Operator::Elastic<1>::Diagonal main ParallelFor:  REG:254  STACK:8 (another variant REG:84–140, STACK:936)
```

`Elastic<1>` = `Set::Sym::Major` (`src/Set/Matrix4.H:12`) — the chamber
instantiation (`NeoHookeanPredeformed::DDW` returns `Matrix4<DIM,Sym::Major>`,
`src/Model/Solid/Finite/NeoHookeanPredeformed.H:32`; Flame derives from
`Base::Mechanics<NeoHookeanPredeformed>`, `src/Integrator/Flame.H:25`).

This **resolves the long-standing 87–101-regs-static vs 255-regs-runtime
discrepancy** (SUCCESS_BOOK §9.4): the old static dump
(`benchmark/g0_cuda_resource_usage.txt`) was the **2D** binary. The 3D `Sym::Major`
kernel hits the 255-register cap *statically*, on sm_86, exactly as ncu measured
on sm_80. Consequence: **register/occupancy work on Fapply can be A/B'd locally
on the A1000 by static resource usage** — no NOVA round-trip needed for that half
of the evidence (wall-clock still needs A100). Tooling:
`benchmark/fapply_register_ab.sh` (added this session).

Provenance caveat: the shipped binary's STACK:192 could **not** be reproduced by
any same-flags rebuild with the current toolchain (nvcc 12.6.3 gives STACK:0 at
REG:255, and even the main tree's own `obj/` object shows STACK:0) — the binary
was built with older/unpinned flags or toolchain. Two actions: (a) treat the
STACK:192 as historical; (b) make the build scripts (local + NOVA) dump
`fapply_register_ab.sh` output into the build log so every binary carries a known
resource profile — the NOVA production builds may still be on the spilling
toolchain, which is free performance if so.

### 2.2 The staged C1 register edits were one command away from being lost — now committed

The two Phase-C1 `Fapply` edits (single live grad(C) `Matrix4` temp; `sig` built
only in the boundary branch) existed **only as uncommitted worktree state** in
`/home/jackplum/Projects/alamo-elastic-opt`; branch `chamber-gpu-elastic-opt` had
zero commits ahead of `chamber-gpu`. All the *other* diffs in that worktree
(NodeBilinear elixir, Newton device-error placement, traction-diagnostic `boxhi`
fix) had already landed in `chamber-gpu` separately — the Elastic.cpp edits were
the only unique content. **Now committed as `0bb893acc` on
`chamber-gpu-elastic-opt`.** Static register A/B (this session, sm_86 3D, fast-math):

The full same-toolchain static A/B matrix is in **§2.2a** (end of file) — headline:
C1 moves 255→244 registers (modest), and the probe series shows the register
count is **allocator-driven, not code-driven**, with the launch-bounds cap
spill-free down to at least 80 registers. That finding re-orders Phase 3 (see
§2.10 and §3).

### 2.3 Current `chamber-gpu` still pays the wasted interior `sig` contraction

`src/Operator/Elastic.cpp:545` computes
`Set::Matrix sig = (DDW(i,j,k) * gradu) * psi_avg;` **unconditionally**, but `sig`
is consumed only inside the boundary branch (`:550-553`). Every interior node —
the overwhelming majority — pays an 81-FMA `Matrix4×Matrix` contraction and
discards it. (Verified directly; this is exactly what the C1 edit removes. Until
3.A lands, production still pays it.)

### 2.4 Two-thirds of the grad(C) contraction FLOPs are discarded, and the biggest contraction routes through a 45-way if/else

New code-level findings sharpening v3 task 3.D:

- `Elastic.cpp:634-636` (the `!m_uniform` branch, **always live** because
  `Mechanics.H:198` calls `SetUniform(false)` unconditionally): each
  `(Cgrad_d * gradu)` computes the full 9-entry matrix product (81 FMA) then keeps
  **one column** (3 entries). 162 of 243 FMA discarded per node. A column-restricted
  contraction is a pure win and can preserve accumulation order (bit-identity
  candidate).
- `Elastic.cpp:626` `f = (DDW(i,j,k)*gradgradu)*psi_avg` uses
  `operator*(Matrix4<3,Major>, Matrix3)` at `src/Set/Matrix4_Major.H:552-565` —
  a naive 4-nested loop, self-flagged `// TODO: improve efficiency of this method`,
  which indexes through the symmetry-collapsing `operator()` — effectively a
  **45-way if/else chain resolved per access, per node, in the hottest kernel in
  the codebase**. The sibling `Matrix4×Matrix` operator (`:537-549`) is already
  hand-unrolled with direct `data[]` access; `Matrix4×Matrix3` needs the same
  treatment.
- `DDW(i,j,k)` (360 B) is loaded from global memory **twice** per node
  (`Elastic.cpp:613` and `:629` region; also the `:545` use) — cache once
  (`GPU_AUDIT_20260702.md` §5 flagged this; still unfixed).
- Per interior node, `Fapply` loads **7 distinct `Matrix4`s (7×360 B = 2,520 B)**
  — center + 6 face-neighbors for the grad(C) central differences — and each
  node's tangent is re-read by up to 6 neighboring threads with zero reuse
  capture. This is the concrete number that motivates the stored-stencil design
  in §3.4.

### 2.5 Phase field: a per-substep kernel + halo exchange that provably does nothing

`Base::Mechanics<MODEL>::Advance` (`src/Integrator/Base/Mechanics.H:399-412`)
launches, **every substep at every level** whenever elastic isn't `Disable`:
a `ParallelFor` calling `model(i,j,k).Advance(dt, eps, sig, time)` over the grown
nodal box, followed by `model_mf[lev]->FillBoundaryAndSync(...)`.

For the chamber model, `Advance` is the **empty inherited no-op**
`src/Model/Solid/Solid.H:64` (`void Advance(...) {}` — neither `NeoHookean` nor
`NeoHookeanPredeformed` overrides it). So this is a full kernel launch + a full
nodal halo exchange + sync, per level per substep, to execute an empty function.
Fix shape: a compile-time trait (e.g. `static constexpr bool advances = false;`
in `Solid`, overridden by models with real kinetics) gating the whole block.

### 2.6 Phase field: up to 3 full-device synchronizations per level per substep for NaN trapping

`Util::AbortIfDeviceError` calls `amrex::Gpu::streamSynchronizeAll()`
(`src/Util/Util.H:156` — the correct 4.G race fix). But `Flame::Advance` invokes
it after the phase-field kernel (`Flame.cpp:835-836`) **and** after the thermal
kernel (`Flame.cpp:892-893`), i.e. two all-stream barriers per level per substep,
plus the `Mechanics::Advance` sync above. These barriers defeat the cross-box
async overlap the MFIter stream pool exists to provide. Error *detection* does not
need this immediacy: the device flag is sticky — it can be read **once per coarse
step** (all levels, after the level loop) with identical abort semantics, just a
slightly later abort point. This is the highest-leverage *flame-only* structural
fix that isn't already on the codex branch (§2.8).

### 2.7 Phase field: `UpdateModel` is not gated by `elastic.interval`, and rebuilds via 3 separate kernels what one would do

`Mechanics.H:185` calls `UpdateModel(...)` **before** the interval gate at
`Mechanics.H:192`. So on every coarse step — including the (interval−1)/interval
steps where the elastic solve is skipped — `Flame::UpdateModel`
(`Flame.cpp:438-578`) still does, per level: 3 `FillBoundary` calls (phi/eta/temp),
a 3-kernel psi rebuild (`MultiFab::Copy` + `mult` + `plus`, `Flame.cpp:457-462`,
for the affine `psi = floor + (1-floor)·eta`), up to 2 `ParallelFor`s rebuilding
the full per-node 3-material model blend, and a `Util::RealFillBoundary` on the
360 B/node model field. With `elastic.interval=50–200` (the C0 levers), ~99% of
these rebuilds produce state nobody consumes. Fixes: (a) gate the model/RHS
rebuild on the same interval predicate as the solve (verify psi's non-elastic
consumers first), (b) fuse the psi rebuild into one `ParallelFor` regardless.

### 2.8 Phase-field structural work already exists on an unmerged branch — with measurements

Branch `codex/gpu-pf-structural-speedups` (worktree
`/home/jackplum/Projects/alamo-pf-gpu-opt`, 5 commits, based on `7e972f1e8` —
**behind** current `chamber-gpu`) already implements: `evolving=false` registration
for static/diagnostic fields (psi, eta_0, L, mdot, alpha, heatflux, laser,
temp_old, eta_old, phi…), skipping their per-substep `FillPatch` **and** their
per-step `average_down` (`Integrator.cpp` now checks `evolving_array[n]`).

Measured locally (A1000, its `benchmark/GPU_TEST_PERF_TRACKING.md` additions):
P2 AMR3 case **1.45× faster** (1.65 s → 1.14 s median), `cudaLaunchKernel` calls
**−70.5%** (108,035 → 31,835), memcpy/memset API calls −72.9%,
`FillPatchTwoLevels` inclusive 0.916 s → 0.208 s. C4 CPU-vs-GPU AMR correctness
PASS. This branch independently validates the "launch/fill churn, not stencil
throughput" theory of flame AMR cost. **Action: rebase onto current `chamber-gpu`,
run the Phase 1 physics-budget gate (which did not exist when it was written), and
land it.** It also directly overlaps two findings above (it sets `eta_old`
non-evolving; it does *not* address §2.5, §2.6, or §2.7).

### 2.9 Multi-GPU: the "halo-bound" diagnosis is too narrow — and partly wrong

Re-examination of what was actually run (`nova_flame_gpu_3d_multi.slurm`,
`input_3d_flame_{128,256,512}`, elastic **disabled** in all multi-GPU rows):

| # | Shortcoming | Status | Evidence |
|---|---|---|---|
| 1 | **All communication is blocking** — zero `FillBoundary_nowait/_finish` anywhere in `src/` | confirmed | grep; hot paths use blocking `FillBoundary`/`FillBoundaryAndSync`/`ParallelCopy` only |
| 2 | **GPU-aware MPI never enabled or verified** — `amrex.use_gpu_aware_mpi` appears nowhere; AMReX defaults it off; NOVA build scripts never check `MPIX_Query_cuda_support` | confirmed | grep of inputs/slurm/build scripts; `AMReX_ParallelDescriptor.cpp` defaults |
| 3 | **`amr.regrid_int` defaults to 2 and is never overridden** (`Integrator.cpp:54`) → ~32,500 collective regrids (tag gather + DistributionMapping rebuild + ParallelCopy of every field) over a 65,000-step run — a size-independent, rank-scaling tax | confirmed | code + inputs |
| 4 | **Benchmark mode ran on managed memory** (`amrex.the_arena_is_managed=1` in the slurm's bench branch) — cross-GPU touches page-migrate instead of explicit transfers; NVLink topology never confirmed | confirmed (impact unprofiled) | slurm file; NOVA_SLURM_RUNBOOK "pending" note |
| 5 | **128³ level 0 is a single 128×128×64 box** (`max_grid_size=128` ≥ domain) — one GPU idles at level 0. The code even warns about exactly this in the *elastic* path (`Operator.cpp:456-461`) but the Flame path is unguarded | confirmed | AMReX `BoxList::maxSize` semantics |
| 6 | Coarse-fine subcycling sync of the small level-1 flame-front patch | plausible, unprofiled | no nsys/ncu exists for any multi-GPU run |
| 7 | Elastic MLMG node-centered comm (~10 `ParallelCopy` sites + `FillBoundaryAndSync` per MG level, no agglomeration/consolidation configured in `Linear.H`) | latent — elastic was off in all multi-GPU rows | code read |

Critically, 256³/512³ decompose fine (4/32 boxes of 128³, halo surface/volume
≈4.7%) yet regressed *as bad or worse* than the degenerate 128³ case — so **box
granularity/halo volume alone cannot be the cause**; the size-independent
blocking/latency taxes (#1–#4) dominate. The v3 framing of 5.A ("revisit at larger
per-GPU domains") is necessary but **not sufficient**: without #1–#4 fixed, larger
domains will still be gated by blocking-comm latency × regrid frequency.

---

## 3. Structural plan — elastic solver (the project's critical path)

Ordering respects v3: 2.C (ncu Speed-of-Light) decides 3.C-vs-3.D lead; everything
gates on the Phase 1 physics budget. Items marked **bit-exact-able** can
additionally be verified by the cheaper CPU golden compare first.

### 3.1 Land the C1 edits (v3 task 3.A) — now durable, ready for A100 A/B
Commit `0bb893acc` on `chamber-gpu-elastic-opt`. Static sm_86 A/B in §2.2/§2.2a;
A100 wall + occupancy A/B still owed. Note: the three C1 evidence docs
(`PHASE_C1_cpu_golden_compare.md`, `PHASE_C1_fapply_occupancy.md`,
`PHASE_C1_nova_ab.md`) were also untracked worktree files — now committed as
`123de00a2` **on the branch**; they are not visible from the main tree until 3.A
lands.

### 3.2 Launch-bounds sweep (v3 task 3.F — **promoted**, see §2.10)
Sweep `__launch_bounds__(256, {1,2,3,4})` on the Fapply (and Diagonal) launches
on A100: wall/step + ncu achieved-occupancy + budget gate per point. Statically
spill-free to at least 80 regs (§2.2a rows 4–5), config-grade cheap, and the
single most promising near-term lever this investigation found. Use AMReX's
existing `launch_global<MT, min_blocks>` overload scoped to the elastic kernels;
do not patch the global header.

### 3.2b Cheap kernel surgery (v3 task 3.D, first wave) — judged by wall time, not registers
In one pass over `Fapply` (and mirrored in `Diagonal`):
1. Load `DDW(i,j,k)` once into a local, use for all three consumers (§2.4).
2. Column-restricted `Cgrad_d × gradu` contractions — 27 FMA instead of 81, ×3 (§2.4).
3. Hand-unroll `Matrix4<3,Major> × Matrix3` with direct `data[]` indexing,
   replacing the 45-way if/else accessor path (`Matrix4_Major.H:552-565`).

§2.2a row 2 shows these will **not** move the static register count — their case
is instruction count and load elimination, so the metric is A100 wall + ncu
executed-instructions, and the expected win is modest. Bit-exact-able variants
first (CPU golden compare), then budget gate.

### 3.3 Interior/boundary kernel split (v3 task 3.B) — payoff revised down
Boundary handling (`Elastic.cpp:550-553` → the 26-way `BC::eval` if-chain in
`src/BC/Operator/Elastic/Elastic.H:186-247`) leaves the interior kernel entirely:
one launch over the shrunk interior box, one small launch over the boundary
shell. **§2.2a row 3 refutes the register rationale** (identical allocation with
the BC path compiled out); the surviving benefits are warp-divergence removal at
boundary tiles and per-kernel launch bounds. Keep on the list, but behind 3.2.

### 3.4 NEW — Stored-stencil operator for the fine levels (the big structural bet; candidate task **3.G**)
**Idea:** stop re-deriving the operator action from the 360 B/node tangent on
every one of the ~100 `Fapply` calls per Newton iteration. For a *fixed* Newton
state, `Fapply` is a linear stencil: `f(x) = Σ_{δ∈N(x)} B_δ(x) · u(x+δ)` with
19-point support in 3D (center + 6 face + 12 edge neighbors — corners are not
touched by either `gradgradu` or the grad(C) differences) and 3×3 blocks. Build
`B_δ(x)` **once per Newton iteration per level** (a single Fapply-strength kernel
that folds DDW, grad(C), psi, and psi-gradient terms into blocks), then every
smoother/residual/bottom application is a block-SpMV.

Arithmetic:
- **Stored bytes/node:** 19 × 9 × 8 B = **1,368 B** (vs 360 B today → ~3.8× more
  coefficient memory).
- **Traffic per apply per node:** today ≈2,520 B of `Matrix4` loads + ~456 B of u;
  stored-stencil ≈1,368 B + 456 B — ~35–40% less, *and* streaming-friendly instead
  of 45-double AoS strided.
- **Registers per thread:** the apply kernel becomes accumulator + one block in
  flight — plausibly <64 regs → 4–8× more resident warps. Since the kernel is
  memory-latency-bound at 12.5% occupancy, this is where the real win lives; the
  combination is the plausible 2–4× on the elastic solve.
- **Amortization:** ~10 V-cycles × (2×(pre+post) smooth + residuals) on the fine
  level ≈ 100 applies per build. Build cost ≈1 apply. Negligible.
- **Bonus:** `Diagonal` becomes "read the diagonal of B₀" — the per-Newton-iter
  full-hierarchy `Diagonal(true)` recompute (`Elastic.cpp:801-947`, with its own
  wasted-contraction pattern at `:883` and STACK:936 kernels) disappears.
- **Cost/risk:** at 256³, fine level ≈ 17M nodes × 1,368 B ≈ **23 GB** (+ MG
  hierarchy ≈ +14%) vs ~7 GB today for `m_ddw_mf`. Fits A100-80 for the current
  case but is the binding constraint. Mitigations, in preference order:
  (a) store stencils **only for the fine one or two MG levels** where ~90% of
  Fapply time lives (PHASE_A §3 bimodality), keep matrix-free coarse levels;
  (b) exploit operator self-adjointness (Major symmetry ⇒ B_{−δ}(x) = B_{+δ}(x−δ)ᵀ)
  to store ~10 of 19 blocks (~720 B/node) at the cost of a transposed gather;
  (c) FP32 coefficient storage with FP64 accumulation — **only** through the
  physics budget, given the documented high-contrast sensitivity
  (`MLMG_HIGH_CONTRAST_FINDINGS.md`).
- **Not bit-exact** (different accumulation order) — this is exactly what the
  Phase 1 budget exists to adjudicate. Keep the current matrix-free path as the
  CPU/reference implementation (v3 principle: specialize only the GPU hot path).
- **Decision input:** if 2.C's SoL says memory-BW-bound even at fixed occupancy,
  option (b)/(c) compression matters more; if latency/occupancy-bound, plain (a)
  suffices. Do 2.C first; it was designed to order exactly this decision.

### 3.5 Uniform-material fast path (v3 task 3.C)
Complementary to 3.4 (and mooted for whichever levels get stored stencils): a
per-box uniformity test (psi min/max plus model-field constancy) selecting a
no-grad(C) kernel variant. Most chamber volume (deep propellant, casing) is
locally uniform. If 3.4 lands fine-level-only, this still helps the matrix-free
coarse levels.

### 3.6 Operator & MLMG persistence across solves (v3 task 3.E)
`Mechanics.H:197-225` constructs a fresh `Operator::Elastic` (full
`define()` — `m_ddw_mf`/`m_psi_mf`/`m_diag` across every amrlev×mglev) **and** a
fresh `amrex::MLMG` every solve, then tears both down. Persist both across
timesteps; invalidate on regrid only; per solve just `SetModel` + re-restrict
coefficients. Removes multi-GB alloc/free churn per solve (and interacts well
with 3.4's larger footprint — allocate once, not per solve). Risk: lifetime
reasoning vs the existing `streamSynchronizeAll` UAF fixes (`Mechanics.H:225`,
`Newton.H:459`); do after 4.B's async-lifetime audit or with it.

### 3.7 NEW — Fsmooth fusion (candidate task **3.H**)
`Operator.cpp:350-417`: each `Fsmooth` = 2 × (`Fapply` + separate axpy-style
MultiFab ops + per-component `ParallelFor` + `FillBoundary` + `nodalSync`) ≈ 16+
launches, re-reading the fields it just wrote. Fuse apply+Jacobi-update into one
kernel (trivial in the stored-stencil world of 3.4: SpMV+update). Halves smoother
global traffic and cuts ~10 launches per smooth call. Bit-exact-able if the
Jacobi arithmetic order is preserved.

### 3.8 Small but real: per-V-cycle `tmpfab` allocation in `interpolation()`
`Operator.cpp:715-728` allocates an `FArrayBox` + elixir per box per level
transition per V-cycle (the historical UAF site). With 3.6's persistent operator,
pre-allocate per-level scratch once. Minor; bundle with 3.6.

---

## 4. Structural plan — phase field (Flame)

Flame is 0.2% of *combined* wall — these matter for flame-only production runs,
for the CPU-comparison story, and for multi-GPU (where launch/sync churn is the
poison). Ordered by measured or mechanically-certain payoff:

1. **Land `codex/gpu-pf-structural-speedups`** (§2.8) — measured 1.45× on the AMR
   case, −70% launches. Rebase onto `chamber-gpu`, pass the Phase 1 gate, land.
   NEW task **PF.1**.
2. **Batch the device-error check to once per coarse step** (§2.6) — removes 2–3
   `streamSynchronizeAll()` per level per substep while keeping the 4.G
   correctness fix (sticky flag, later read). NEW task **PF.2**. Should be
   A/B'd with the same local harness the codex branch used
   (`benchmark/local_pf_gpu_ab.sh`).
3. **Skip the no-op `model.Advance` kernel + `FillBoundaryAndSync`** for models
   without kinetics via a `static constexpr` trait (§2.5). NEW task **PF.3**.
4. **Gate `UpdateModel`/psi/RHS rebuild on the elastic interval; fuse the psi
   rebuild to one kernel** (§2.7). NEW task **PF.4** (verify psi consumers before
   gating; fusion is unconditional).
5. **Fuse phase-field + thermal kernels** (`Flame.cpp:723` and `:866` — same box,
   overlapping state, currently separated by a swap + full-device sync). NEW task
   **PF.5**, lower priority.
6. `eta_old` FillPatch waste: already fixed on the codex branch (part of PF.1).

---

## 5. Multi-GPU — prerequisites before it is worth re-measuring (feeds v3 5.A)

The 5.A framing ("bigger per-GPU domains") is incomplete (§2.9). Before *any*
re-measurement, all of:

1. **Enable + verify GPU-aware MPI**: build against CUDA-aware OpenMPI, set
   `amrex.use_gpu_aware_mpi=1`, and assert `MPIX_Query_cuda_support()==1` in the
   run log (one-line print at startup). Without this, every halo is
   D2H→MPI→H2D.
2. **Stop regridding every 2 steps in benchmarks**: set `amr.regrid_int` ≥ 16 (or
   static grids for the bench) and report regrid cost separately. ~32,500
   collective DM rebuilds per run is a rank-scaling tax unrelated to halos.
3. **Device arena, not managed**, for A100 bench rows (managed was only ever a
   local-A1000 8 GB necessity).
4. **≥4–8 boxes per GPU at level 0** (128³ needs `max_grid_size ≤ 64`); keep the
   single-GPU row at the same box size so the comparison isolates communication.
5. **Comm/compute overlap** in the Flame advance path (`FillBoundary_nowait` /
   compute-interior / `FillBoundary_finish`) — the only *code* change in this
   list; the rest are build/config. Elastic's ~10 blocking `ParallelCopy` sites
   per MG level (`Elastic.cpp:1177-1508`) and unconfigured MLMG
   agglomeration/consolidation (`Linear.H` never calls
   `setAgglomeration/setConsolidation`) are the *next* tier, relevant once
   elastic multi-GPU is attempted at all.
6. **Profile this time**: nsys with NVTX on one rank; capture sync fraction — the
   number that was never collected.

Only after 1–4 (config-level, cheap) does a 2-GPU row measure communication
architecture rather than configuration debt. Expectation management: single-GPU
remains the supported shape (v3 principle 6); multi-GPU matters when the domain
exceeds one GPU's memory (≥512³-class with elastic on, or ~1024³ flame-only).

---

## 6. Sequencing (respecting v3 gates)

```
Phase 1 gate (exists) ─▶ 2.C SoL counters ─▶ 3.A (C1, committed, A/B owed)
                                   │              + 3.F launch-bounds sweep
                                   │                (§2.10 — bundle into the
                                   │                 same A100 session as 3.A)
              ┌────────────────────┼──────────────────────┐
              ▼                    ▼                      ▼
        3.2b cheap kernel    3.3 interior/boundary   PF.1 codex-branch rebase+gate
        surgery (3.D wave 1)      split (3.B,        PF.2 sync cadence
              │                   payoff revised ↓)  PF.3 no-op Advance skip
              ▼                                      PF.4 UpdateModel gating
        3.4 stored-stencil (3.G, the big bet)
        + 3.7 Fsmooth fusion (3.H)
        + 3.6 operator persistence (3.E)
              │
              ▼
        5.A multi-GPU — only with §5 prerequisites 1–4 done
```

Every arrow lands through: strict-build golden compare where bit-exact-able →
`benchmark/validate/` budget gate → A100 A/B for perf claims (fast-math build,
labeled).

---

## 7. Artifacts & state produced by this investigation

- Commits `0bb893acc` (C1 edits) + `123de00a2` (C1 evidence docs, previously
  untracked) on `chamber-gpu-elastic-opt` — now durable.
- `benchmark/fapply_register_ab.sh` — one-command static register/stack A/B for
  the elastic hot kernels (works on binaries or single `.o` files).
- The §2.2a probe matrix (6 same-toolchain cells) — the raw dumps were session
  artifacts; every number is recorded in this file and reproducible via the
  recipe in §2.2a.
- This document, linked from the READ_FIRST almanac table and ROADMAP_V3 §7.

Open items this doc *does not* resolve (unchanged from v3): 2.A fair CPU
baseline; 2.C ncu SoL export; A100 A/B for 3.A; Phase 1 NOVA strict-build gap.

---

## 2.2a — The static register probe matrix (2026-07-03, new measurements)

All cells: same tree (`alamo-elastic-opt` worktree), same toolchain (nvcc
12.6.3), same flags (`--cuda 86 --cuda-fp fast`, 3D, `-maxrregcount=255`,
`-O3 --use_fast_math`), `Elastic<1>` (= `Sym::Major`, the chamber instantiation)
**main `ParallelFor` kernel** of `Fapply` via
`cuobjdump --dump-resource-usage` (`benchmark/fapply_register_ab.sh`):

| # | Variant | REG | STACK (spill frame) |
|---|---|---|---|
| 0 | pre-C1 (`chamber-gpu` tip code) | 255 (cap) | 0 |
| 1 | C1 edits (`0bb893acc`) | 244 | 0 |
| 2 | C1 + derivative-free grad(C) probe (no `Matrix4` temp ever materialized) | 244 | 0 |
| 3 | C1 + boundary path compiled out (simulates 3.B interior kernel) | 244 | 0 |
| 4 | C1 + `__launch_bounds__(256, 2)` → 128-reg ceiling | **128** | **0** |
| 5 | C1 + `__launch_bounds__(256, 3)` → ~85-reg ceiling | **80** | **0** |

(Reference: the *shipped* `bin/alamo_gpu-3d-cuda86-g++` shows 255/**192** — not
reproducible with the current toolchain; see the provenance caveat in §2.1.
`Diagonal` follows the same pattern: 189 uncapped → 125 @minb2 → 80 @minb3, all
STACK:0.)

Probe code: rows 2–3 were temporary edits in the `alamo-elastic-opt` worktree
(marked `NOT LANDABLE`, since reverted — recipe reproducible from this table);
rows 4–5 were a one-line `__launch_bounds__` patch to the *installed* AMReX
header (`ext/AMReX-Codes/amrex/3d-cuda86-g++-26.06/include/AMReX_GpuLaunchGlobal.H`),
also reverted. A landable version scopes the bound per-kernel (AMReX already
provides the `launch_global<MT, min_blocks>` template overload) rather than
patching every kernel in the TU.

### 2.10 What the matrix means (revises §3 priorities)

1. **ptxas fills whatever budget it is given.** Rows 1–3: removing the derivative
   `Matrix4` temp *and* the entire boundary/BC path changed the register count by
   exactly zero — with 255 available, the allocator spends registers on ILP
   scheduling regardless of source-level liveness. Static REG is therefore
   **not** a useful success metric for kernel micro-surgery (wall time and ncu
   counters are); conversely, micro-surgery is not how occupancy improves.
2. **The occupancy lever is the launch bound, and it is spill-free down to at
   least 80 registers.** Rows 4–5: 2× and 3× the resident warps with zero spill
   frame. The cost is reduced per-thread ILP and more re-loads through L1 — for
   a kernel stuck at 12.5% occupancy and latency-bound behavior, this trade is
   usually strongly net-positive, but it is empirical: **sweep
   `__launch_bounds__(256, {1,2,3,4})` on A100 with wall + ncu occupancy +
   budget gate.** This promotes v3 task 3.F from "after 3.D" to **immediately
   after (or bundled with) 3.A** — it is config-grade cheap and may deliver a
   large fraction of what the kernel rewrites were expected to.
3. **3.B's register rationale is refuted; its divergence rationale survives.**
   Row 3 shows the interior kernel allocates identically without the BC path.
   The split's remaining value: warp-divergence removal at boundary tiles, and
   the freedom to give interior/boundary kernels different launch bounds. Expect
   modest wins; do not lead with it.
4. **3.G (stored-stencil) is unaffected** — its case was always memory traffic
   (2,520 B/node of `Matrix4` loads per apply, re-read ~100× per Newton
   iteration) and streaming access patterns, not register count. Rows 4–5
   actually *help* 3.G: a low-register SpMV-style apply kernel plus a high
   launch bound compound.
5. **C1 (3.A) remains worth landing** — bit-identical, CPU-golden-verified, and
   removes real discarded work (§2.3) — but its expected A100 wall delta should
   be sized as "small, possibly noise-level"; the A/B that matters most now is
   the launch-bounds sweep.

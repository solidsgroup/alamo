# AMReX GPU-practice audit — 2026-07-09 (session notes)

Audit of chamber-gpu src/ against AMReX GPU documentation
(amrex-codes.github.io/amrex/docs_html/GPU.html) + modern GPU practice.
Scouted inventory + doc comparison; items NOT in this task's scope are
recorded here for the live plan backlog.

## Already conformant (no action)
- Reductions: ReduceOps/ReduceData used in all hot diagnostics
  (Flame.cpp:620/1060, Elastic.cpp:27/91/687+, Mechanics.H:461+). Matches
  AMReX "fused device reductions" guidance.
- Temp lifetimes: `elixir()` at Operator.cpp:728, NodeBilinear.H:48 (the
  cross-stream UAF fix). AMReX now suggests The_Async_Arena() as the modern
  replacement for Elixir — optional modernization, zero perf delta expected.
- MFIter: `TilingIfNotGPU()` dominant (~60 sites); matches doc ("disable
  tiling on GPU, maximize work per launch").
- Host/device divergence is confined to launch macros + BC virtual-dispatch
  workaround (Elastic.cpp:158/820, Newton.H:14, Operator.cpp:17). Sound.
- Gpu::DeviceVector staging only in I/O paths (StarAftGrain.H:369, BMP.H:103).

## Gaps found, deferred (backlog candidates)
1. **No `__launch_bounds__` anywhere in src/** — AMReX exposes
   `ParallelFor<MY_BLOCK_SIZE>` / `launch_global<MT,min_blocks>`. Already
   live-plan task 3.2; requires A100 ncu occupancy judgment — NOVA off-limits
   this session, so not touched.
2. **Fusible per-component ParallelFor loops**: Linear.H:83-85,
   Operator.cpp:389-392, IC/Constant.H:62-65 — component loop outside the
   kernel means N small launches instead of one 4D launch. AMReX doc: pass
   ncomp to ParallelFor ("component loop is moved to the innermost loop").
   Low-risk, measurable only at small-box/coarse-MG levels (launch-latency
   bound). Backlog candidate.
3. **Per-component `Dot()`/`norminf()` in Operator.cpp:56-57/82-83 and
   Newton norms (live-plan 3.I)** — each is a separate device reduction +
   sync + allreduce. Already tracked as live-plan 3.I; convergence-critical,
   not this session.
4. **Arena policy** — no explicit arena selection in src/; managed default.
   Already live-plan Phase 5 A/B (A100). A1000 (8 GB) needs managed; no
   local action.
5. **`amrex.use_gpu_aware_mpi=1`** — AMReX doc recommends for multi-GPU;
   already implicated in structural plan multi-GPU loss analysis. NOVA-only.

## In scope this session (live-plan 3.2b)
- Fapply: DDW(i,j,k) global load repeated 3x per node (Elastic.cpp:532,613,629)
  — MATRIX4 Major = 45 doubles = 360 B/load in 3D; hoist to one load.
- Fapply: `(Cgrad_d * gradu).col(d)` computes full 3x3 (81 FMA per product)
  then discards 2/3 — column-restrict to 27 FMA x3.
- Matrix4_Major.H:552-565 `Matrix4 x Matrix3`: quadruple loop over branchy
  10/45-way if-else accessor — branch-heavy, register-hungry; hand-unroll
  with direct data[] indexing (modern GPU practice: branch-free constexpr
  indexing in hot kernels).
- Diagonal: DDW re-fetched inside per-component p-loop (Elastic.cpp:860,868).

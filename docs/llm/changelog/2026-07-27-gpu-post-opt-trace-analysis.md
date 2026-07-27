# GPU post-optimization trace analysis — launch-gap ceiling — 2026-07-27

Follow-up to `2026-07-26-gpu-optimization-sweep.md`. No source changes. Analysis
only, re-derived from `final.sqlite` of the post-optimization nsys trace; full
write-up in the addendum of
`docs/agent_plans/20260726-gpu-optimization-sweep/results/RESULT.md`.

## TL;DR

The 2D-conservative case is now launch-latency bound, not compute bound.

| metric | value |
| --- | --- |
| span, first to last kernel | 776.0 ms |
| kernel busy | 421.8 ms |
| GPU idle | 354.2 ms (**45.6%** of span) |
| gaps | 32,771, median 6.34 us, mean 10.8 us |
| share of GPU time in the top 2 kernels | 91.2% |
| explicit memcpy traffic | 268 transfers, ~0 GB, 0.12 ms |
| launches on a single stream | 32,484 of 32,772 |

## Consequences

1. **Further per-kernel tuning has a hard ceiling.** Halving `Fapply` again
   would cut only ~27% of wall. Launch count and host-side per-launch cost are
   the remaining levers.
2. **Coarse MG levels cost ~25% of the window for near-zero work.** 6,970 of
   11,808 `Fapply` launches (59%) run on <= 9 blocks on a 16-SM device: 111 ms
   kernel plus ~86 ms attributed gap. This promotes coarse-level elimination to
   the top of the candidate list, ahead of face-flux fission. It is *not* a
   block-width problem. The `MLMG_HIGH_CONTRAST_FINDINGS.md` warning against
   capping `max_coarsening_level` without a stronger bottom solver still binds.
3. **Effectively single-stream.** AMReX rotates streams per MFIter box
   (`AMReX_MFIter.cpp:378`, `max_gpu_streams` default 4); the observed
   distribution means one box per level, so no intra-rank concurrency exists on
   this deck. Needs checking against the production `max_grid_size`.
4. **The 3D 26 GB pressure is UVM page migration, not explicit copies** —
   consistent with the abandoned 3D local timing, and confirms that trimming the
   unread `m_ddw_mf` face components targets the right mechanism.
5. **New cheap candidate:** level-dependent `Fapply` block width.
   `ALAMO_ELASTIC_FAPPLY_MT` is one compile-time constant; 128 wins because it
   helps the small-grid launches, but the grid=100 tier is 32% of `Fapply` time
   and is already SM-saturated. Block size cannot change results.
6. **To verify:** the five `Gpu::streamSynchronizeAll()` calls in
   `src/Solver/Nonlocal/Newton.H` (lines 777, 860, 1001, 1294, 1300) — each is a
   full device barrier; confirm none sits inside a per-relinearization loop.
   `Fapply` itself is clean in release builds (its `DeviceErrorFlag` is
   `AMREX_DEBUG`-gated at `src/Operator/Elastic.cpp:200-205`).

## Evidence

Derived tables committed as CSV in
`docs/agent_plans/20260726-gpu-optimization-sweep/artifacts/nsys/`:
`final_launch_gap_summary.csv`, `final_kernel_summary.csv`,
`final_fapply_grid_histogram.csv`, `final_stream_distribution.csv`,
`final_memcpy_summary.csv`. The `.nsys-rep` and `.sqlite` binaries they were
derived from are not tracked, matching the 20260721 task's convention.

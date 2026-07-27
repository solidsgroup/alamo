# GPU optimization sweep — result

Branch `chamber-gpu`, local RTX A1000 (sm_86, 16 SM, 8188 MiB).
Source at start: HEAD `14767e4d7` plus the pre-existing uncommitted
`Integrator/Flame.{cpp,H}` chi-rename WIP (untouched by this task).

## Outcome

Three changes measured, two retained, one reverted. Retained set is
**bit-identical to the unmodified source in a strict (`--fmad=false`) CUDA
build**, on both the canonical deck and the sensitivity deck.

| regime | baseline | final | change |
|---|---:|---:|---:|
| 2D conservative wall (median of 5) | 1.82 s | 1.38 s | **-24.18%** |
| 2D GPU kernel time (nsys) | 589.1 ms | 421.8 ms | **-28.40%** |
| 2D kernel launches (nsys) | 74,756 | 32,772 | **-56.16%** |
| 2D `Fapply` kernel time (11,808 calls both arms) | 419.2 ms | 323.7 ms | **-22.78%** |
| 3D (golden-suite section) wall (median of 5) | 8.66 s | 7.81 s | **-9.82%** |

Decks: `tests/ElasticSoftVoid/input` 2D-parallel-shaped args (`max_step=2`,
`amr.max_level=2`, 4/4 smoothing) and its 3D-serial section. Both run
`elastic.use_psi=0`, i.e. the **conservative face-flux** `Fapply` branch, which
is also what the production deck `input_copy` selects
(`Flame.cpp:382` maps `use_psi=0` to conservative face flux).

## Retained changes

### 1. `Operator::Fsmooth` elementwise fusion — `src/Operator/Operator.cpp:356`

Each Jacobi half-sweep staged `Rx = Ax - x*diag` through two whole extra
MultiFabs (`Dx`, `Rx`) and four elementwise kernels (`Copy`, `Multiply`,
`Copy`, `Subtract`) before a fifth kernel consumed `Rx` pointwise. Every one of
those operations is pointwise in `(i,j,k,n)`, so they fold into the update
kernel that already reads the same indices.

Evidence for the target: the 2D nsys trace attributed 75.4 ms (12.8% of GPU
time) and ~42k of 74.8k launches to exactly those four kernels. After fusion
they are gone; the update kernel absorbs 5.7 ms of the removed work.

Isolated measurement: **-18.0%** wall in 2D (2.00 → 1.64 s), **-3.0%** in 3D.
MLMG iteration counts identical (120 iterations, `Final Iter. 37` and `45`
in both arms) — the win is not fewer iterations.

Also removes two MultiFab allocations per `Fsmooth` call, which matters on the
managed arena.

### 2. `Fapply` read-only coefficient access — `src/Operator/Elastic.cpp:199,242`

- `DDW` and `psi` were bound as mutable `Array4` in a kernel that only reads
  them, forcing the compiler to assume they may alias the `F` store. Now
  `const_array`.
- `MATRIX4 const ddw = DDW(i,j,k)` copied a 45-double (3D `Sym::Major`)
  coefficient into registers on a kernel already at 254 registers/thread — and
  the conservative branch never reads it at all. Now bound by reference.

Isolated measurement: **-7.1%** wall in 2D, **-2.0%** in 3D.

### 3. `Fapply` launch width 256 → 128 — `src/Operator/Elastic.cpp:15`

This came out of the trace, not from a guess. Grouping the 2D `Fapply` launches
by grid size:

| grid (blocks) | launches | avg duration | share of Fapply time |
|---:|---:|---:|---:|
| 2 | 3,936 | 30.1 us | 22.5% |
| 5 | 3,034 | 29.4 us | 16.9% |
| 17 | 2,992 | 27.6 us | 15.7% |
| 50 | 1,312 | 80.4 us | 20.0% |
| 52 | 248 | 88.2 us | 4.2% |

Most `Fapply` launches are coarse-MG-level launches occupying **2 to 17 blocks
on a 16-SM device**, each still costing ~30 us. At that size the block width —
not occupancy — decides how many SMs are engaged at all. Halving it doubles the
number of blocks. Occupancy per SM is unchanged (254 regs/thread caps an SM at
~256 resident threads either way), so this is purely about spreading small
launches across the machine.

Isolated measurement: **-9.8%** wall in 2D, **-3.4%** in 3D.
A further halving to 64 was measured and is **not** better (+0.7% vs 128).

Block size cannot change results: every thread writes only its own
`F(i,j,k,.)`, with no reduction or shared state.

## Reverted

| candidate | result | disposition |
|---|---|---|
| `Set::MulCol` in the conservative face-flux branch (compute only the used flux column) | **+0.63%** — within noise | reverted; nvcc already dead-code-eliminates the unused columns of the unrolled `face` loop |
| `ALAMO_ELASTIC_FAPPLY_MT=64` | **+0.71%** vs 128 | not adopted |

## Correctness

| gate | result |
|---|---|
| `benchmark/lint_device_patterns.sh` | PASS (0 non-allowlisted) |
| `GOLDEN_MODE=cpu benchmark/ci_golden_compare.sh` | PASS, 4/4 cases |
| `TIERS=1 benchmark/local_a100_gate.sh` | PASS |
| strict CUDA build vs unmodified source, `input` and `input_rod_and_tube_2d` | `thermo.dat` **bit-identical** |
| `GOLDEN_MODE=gpu` (strict) golden compare | 3/4 ok; `rod_and_tube_step2` FAIL — **pre-existing, not caused by this work** |

The `rod_and_tube_step2/gpu_strict` failure reproduces with **identical deltas**
(`trac_xhi_x abs=1.384000e+02`) using the untouched pre-existing binary
`bin/alamo_gpu-2d-cuda86-g++` (built 2026-07-24, before this task). Cause: the
GPU references were recorded 2026-07-05 and the CPU reference was re-recorded
2026-07-18 during the void-recovery recalibration; the GPU references were never
regenerated. The 2026-07-18 changelog already records this
("gpu_fast/gpu_strict references known-stale until 2D CUDA binaries rebuilt").
The fresh GPU run now agrees with the *CPU* reference. Re-recording the GPU
references is a separate decision and was deliberately not done here.

Note on the fast-math build: with `--use_fast_math` (`-fmad=true`) the fused
`Ax - x*diag` may contract to an FMA, so fast-build results differ from the
old ones at ~1e-9 relative (visible in the MLMG residual trace's last two
digits, same iteration count). Strict builds — the ones the golden gate uses —
are bit-identical.

## Not measured here

The 3D-psi benchmark case `input_3d_centre_bore_128_a2` was started and
abandoned: it has a ~26 GB managed high-water mark on a 7.8 GB card and took
>18 min/rep against the 4.8 min recorded on this machine in the 2026-07-21
task. Local timing of that case measures page migration, not the kernel. It
belongs on NOVA/A100, where it is resident. The `use_psi=1` (non-conservative)
`Fapply` branch is therefore unmeasured by this task; changes 2 and 3 apply to
it, and change 2's register argument is strictly stronger there (the 45-double
`ddw` copy is actually used on that path).

Nsight Compute could not be used: this driver keeps
`NVreg_RestrictProfilingToAdminUsers` at its default, so unprivileged `ncu`
fails with `ERR_NVGPUCTRPERM`. A ready-to-run capture script is at
`scripts/run_root_ncu.sh` (needs `sudo`). All findings above therefore rest on
nsys tracing plus paired wall-clock A/B, not on hardware counters.

## Next candidates, ranked by the evidence collected

1. **Face-flux kernel fission (conservative branch).** Every interior face flux
   is computed twice — once by the thread on each side. Computing fluxes into a
   face-valued scratch field and differencing them in a second kernel halves the
   `Matrix4` coefficient traffic (~1296 → ~648 B/thread in 3D) at the cost of
   9 doubles/node of round-trip. This is the largest remaining structural item
   for the production path.
2. **Shared-memory tiling of `U` in `Fapply`.** The conservative branch reads a
   27-point neighbourhood of `U` per node, issuing ~180 scalar loads for 81
   distinct doubles. An 8x8x4 block plus halo needs 14.4 KB of shared memory.
   Requires a custom launch (AMReX `ParallelFor` gives no tiled block shape).
3. **Coarse-level launch consolidation.** 3,936 of the 2D `Fapply` launches run
   on 2 blocks. Beyond block width, the real fix is not to run the operator on
   levels that small — a bottom-solver/`max_coarsening_level` study. Note the
   standing warning in `MLMG_HIGH_CONTRAST_FINDINGS.md`: never cap `mcl` without
   a smoother bottom solver.
4. **`Fapply` AoS -> SoA for `m_ddw_mf`.** Still the long-deferred item: a
   45-double struct per node means every coefficient load instruction scatters
   across 32 distinct sectors.
5. **Skip the unused face components of `m_ddw_mf` when the conservative path is
   off.** `m_ddw_mf` carries `AMREX_SPACEDIM+1` `Matrix4` components; on the
   `use_psi=1` path components 1..dim are written, coarsened and ghost-exchanged
   but never read. Dropping them cuts the largest array by 75% — directly
   relevant to the 26 GB managed footprint that blocks local 3D work.
6. **Backlog item 3.I (fused Newton norms) is not worth doing on this evidence:**
   `FabArray::norminf` is 0.13% of the 3D profile.

## Evidence paths

- `artifacts/baseline/bin/`, `artifacts/{fsmooth_fusion,expA_constarray,expB_mulcol,mt_sweep}/bin/` — binaries with SHA-256
- `artifacts/*/2d/`, `artifacts/*/3ds/` — per-rep wall times, stdout, stderr
- `artifacts/nsys/baseline.nsys-rep`, `artifacts/nsys/final.nsys-rep`
- `artifacts/strict_ab/` — pristine / candidate strict binaries used for the
  bit-identity comparison
- `scripts/ab_timing.sh`, `scripts/run_root_ncu.sh`

---

## Addendum — post-optimization trace analysis (2026-07-27)

Re-derived from `artifacts/nsys/final.sqlite` after the session above closed.
Five findings; they revise the ranking in the previous section. All numbers are
from the 2D-conservative 4x4 case, single rank, RTX A1000. Derived tables are
committed as CSV next to the trace (the `.nsys-rep`/`.sqlite` binaries are not
tracked, matching the 20260721 task).

### A. The GPU is idle 45% of the solve window — this is now the ceiling

| metric | value |
|---|---|
| span, first to last kernel | 776.0 ms |
| kernel busy | 421.8 ms |
| idle | 354.2 ms (45.6% of span) |
| number of gaps | 32,771 |
| median gap | 6.34 us |
| mean gap | 10.8 us |

A 6-10 us median gap is CUDA launch latency with the host unable to run ahead of
the device. The consequence is a hard limit on further kernel work: even halving
`Fapply` again would cut only ~27% of wall. What is left to attack is launch
*count* and host-side per-launch cost, not per-kernel efficiency.

Source: `artifacts/nsys/final_launch_gap_summary.csv`.

### B. 91.2% of GPU time is two kernels

`Fapply` 323.65 ms / 76.7% over 11,808 launches; the fused `Fsmooth` update
61.0 ms / 14.5% over 10,496 launches. The next largest single entry is 1.1% and
is a one-shot initialization kernel. There is no third kernel worth touching.

Source: `artifacts/nsys/final_kernel_summary.csv`.

### C. Coarse MG levels cost ~25% of the window for near-zero work

`Fapply` launches grouped by grid size:

| grid (blocks) | launches | avg | total |
|---|---|---|---|
| 100 | 1,312 | 78.6 us | 103.2 ms |
| 34 | 2,992 | 26.5 us | 79.2 ms |
| 3 | 3,936 | 16.3 us | 64.1 ms |
| 9 | 3,034 | 15.6 us | 47.3 ms |
| 104 | 248 | 88.8 us | 22.0 ms |
| 36 | 286 | 27.3 us | 7.8 ms |

6,970 of 11,808 `Fapply` launches (59%) run on 9 blocks or fewer — at most 1,152
threads on a 16-SM device. Those launches account for 111 ms of kernel time plus
roughly 86 ms of attributed gap: **~197 ms, about 25% of the 776 ms window.** At
a 16 us duration a launch is dominated by its own overhead, and no block-size
choice changes that.

This promotes candidate 3 from the previous section to first place and reframes
it. The lever is not block width; it is not issuing the operator on levels this
small at all — either capping `max_coarsening_level` **together with** a
stronger bottom solver (the standing warning in `MLMG_HIGH_CONTRAST_FINDINGS.md`
against capping `mcl` bare still applies), or running the bottom levels on the
host. Cheaper to test than face-flux fission and larger in expected effect.

Source: `artifacts/nsys/final_fapply_grid_histogram.csv`.

### D. Effectively single-stream

32,484 of 32,772 launches land on stream 13; streams 14/15/16 take 96 each.
AMReX rotates GPU streams round-robin per MFIter box
(`ext/amrex/Src/Base/AMReX_MFIter.cpp:378`, `max_gpu_streams` default 4), so
this distribution means the benchmark has one box per level and the rotation
never advances. No intra-rank kernel concurrency is available on this deck.
Worth confirming whether the production deck's `max_grid_size` actually yields
more than one box per level; if it does not, multi-stream capacity is unused
everywhere rather than only in the benchmark.

Source: `artifacts/nsys/final_stream_distribution.csv`.

### E. Explicit memcpy traffic is nil

268 transfers, ~0 GB, 0.12 ms total across the whole run. The 3D case's 26 GB
managed high-water therefore manifests as UVM page migration, not as explicit
copies — which is consistent with the abandoned 3D timing attempt above, and
confirms that candidate 5 (dropping unread `m_ddw_mf` face components) targets
the right mechanism.

Source: `artifacts/nsys/final_memcpy_summary.csv`.

### New candidate: level-dependent `Fapply` block width

`ALAMO_ELASTIC_FAPPLY_MT` is a single compile-time constant. 128 won globally
because it helps the 3- and 9-block launches engage more SMs. But the grid=100
tier is 103.2 ms — 32% of all `Fapply` time — and is already SM-saturated, where
the wider block may be better. Dispatching the block width on grid size or
`mglev` instead of one constant is a single A/B with no correctness exposure
(every thread writes only its own `F(i,j,k,.)`; block size cannot change
results).

### Also worth checking

`src/Solver/Nonlocal/Newton.H` has five `amrex::Gpu::streamSynchronizeAll()`
calls (lines 777, 860, 1001, 1294, 1300). `Fapply` itself is clean in release
builds — its `Util::DeviceErrorFlag` is `AMREX_DEBUG`-gated at
`src/Operator/Elastic.cpp:200-205` — but confirm none of the Newton five sits
inside a per-relinearization loop, since each one is a full device barrier.

### Revised ordering

1. Coarse-level elimination (finding C)
2. Face-flux kernel fission
3. `m_ddw_mf` unread-component trim
4. Level-dependent `Fapply` block width
5. Shared-memory tiling of `U`
6. AoS -> SoA for `m_ddw_mf`

Everything from item 2 down competes for a share of the 421.8 ms that is already
only 54% of the wall window; item 1 attacks the idle half directly.

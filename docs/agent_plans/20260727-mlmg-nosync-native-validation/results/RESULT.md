# Native MLMG no-sync validation

Date: 2026-07-27

## Verdict

All six requested checks pass. Retain this only as an opt-in,
geometry-sensitive experiment, never as a general default:

- the ten-step 2D single-box case is substantially faster;
- the ten-step 3D section case is slightly faster on the local RTX A1000 and
  substantially faster on the A100;
- the ten-step many-box case is substantially slower;
- strict baseline/native physics agrees within roundoff on one and two A100s;
- a fixed-BoxArray one-rank/two-rank comparison agrees within roundoff;
- the enabled path is clean under Compute Sanitizer on an A100;
- the current 3D sanitizer deck has a pre-existing stencil error in both arms
  and cannot adjudicate this feature.

The retained option remains default off.

Raw plotfiles, sanitizer logs, and binaries total approximately 1.5 GiB and
remain in the local task artifact directory (with A100 originals in the NOVA
worktree). They are intentionally excluded from version control. This result,
the reproduction scripts, job IDs, binary hashes, compact timing samples, and
physics verdicts are the reviewable evidence package.

## What changed

- `Solver::Nonlocal::Linear` parses `elastic.solver.no_gpu_sync` (default
  false) and passes it to AMReX 26.06 `MLMG::setNoGpuSync`.
- The old prototype's process-global `amrex.max_gpu_streams=1` requirement is
  gone. AMReX itself creates a solver-scoped single-stream/no-sync region.
- `MLMGSyncStateGuard` snapshots both AMReX region flags and restores them if
  `MLMG::solve` throws. AMReX 26.06 restores them only on its normal-return
  path.
- `PrepareMLMG` now calls `MLMG::setThrowException` directly for Alamo's
  documented recoverable-failure mode. Previously `abort_on_fail=0` still
  reached `MPI_Abort`, so the exception-restoration path was unreachable.
- Local harnesses own `amrex.the_arena_init_size=1073741824`; no production
  input deck or NOVA run is given that workstation-specific limit.
- The benchmark guide and live plan now require GPU performance evidence to
  use ten steps by default, report wall/step plus startup calibration, and
  treat synchronized region/solver timers as supporting evidence only.

Source files:

- `src/Solver/Nonlocal/Linear.H`
- `src/Solver/Nonlocal/MLMGSyncStateGuard.H`
- `benchmark/README.md`
- `benchmark/local_a100_gate.sh`
- `docs/llm/PLAN.md`

## Startup-aware local performance

Hardware: RTX A1000 (sm_86). Build: CUDA fast-math, AMReX 26.06. Each run uses
ten completed steps; step one has no elastic solve and steps 2-10 contain 18
MLMG calls total. Arms alternate order. The local-only 1 GiB arena reservation
is identical in both arms. No run passes `amrex.max_gpu_streams`.

| case | reps | one-step startup | baseline wall / step | native wall / step | wall delta | baseline MLMG | native MLMG | MLMG delta |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 2D conservative | 7 | 0.8423 s | 6.2470 s / 0.6247 s | 4.9849 s / 0.4985 s | -20.20% | 5.3782 s | 4.1203 s | -23.39% |
| 3D section | 3 | 0.8554 s | 64.3567 s / 6.4357 s | 63.3035 s / 6.3303 s | -1.64% | 59.2904 s | 58.1824 s | -1.87% |
| 2D many-box (`max_grid_size=16`) | 5 | 0.8286 s | 8.1170 s / 0.8117 s | 10.3462 s / 1.0346 s | +27.46% | 7.1576 s | 9.3873 s | +31.15% |

All paired residual/iteration traces are identical: 7/7 2D, 3/3 3D, and 5/5
many-box. The 3D baseline's first sample was externally disturbed (115.38 s
versus 62.65/64.36 s); the table uses the median, and the exclusive A100 job is
the hardware verdict.

Evidence:

- `artifacts/interleaved/2d/local-full-v2/`
- `artifacts/interleaved/3d/local-full-v2/`
- `artifacts/interleaved/manybox/local-full-v2/`
- `scripts/interleaved_native_ab.sh`

Interpretation: AMReX's native API trades all MLMG box-level stream concurrency
for fewer host-side synchronization gaps. That is valuable for the 2D
single-box launch-bound case, nearly neutral for the compute-bound 3D case,
and harmful when many boxes supply useful stream concurrency.

## Exception/failure behavior

Two checks pass:

1. A CUDA probe deliberately changes both AMReX region flags inside
   `MLMGSyncStateGuard`, throws, and verifies restoration for outer
   `(false,false)` and `(true,true)` states.
2. A real 2D MLMG solve with `max_iter=0`, `abort_on_fail=0`, and no-sync
   enabled now reports nonconvergence, completes step 2, finalizes AMReX, and
   exits zero. Before the direct `MLMG::setThrowException` call, the same
   command exited 6 through `MPI_Abort`.

Evidence:

- `artifacts/exception/probe.log`
- `artifacts/exception/forced_failure_recoverable.log`
- `scripts/mlmg_sync_guard_probe.cpp`

## Local correctness and gates

- CUDA release builds: 2D PASS, 3D PASS.
- CUDA strict/no-fast-math 2D build: PASS.
- CPU unit binaries: 2D PASS (zero failures), 3D PASS (zero failures).
- Device-pattern lint: PASS, zero non-allowlisted findings.
- Strict canonical 2D baseline/native comparison: PASS for all 27
  observables; field differences are approximately `1e-15` relative and
  solver iteration counts match.
- Stored CPU-strict golden: FAIL for both the current baseline and native arm
  in the same way (roughly 1-2% elastic-field drift and 23 versus 26 MLMG
  cycles). This is current-worktree/reference drift, not a no-sync delta.
- `make check`: blocked by 11 existing/current-worktree undocumented inputs,
  including the user's `chi` additions and pre-existing Newton knobs.

Evidence:

- `artifacts/test-2d-cpu.log`
- `artifacts/test-3d-cpu.log`
- `artifacts/device-lint.log`
- `artifacts/physics/local-baseline-vs-native.{md,json}`
- `artifacts/physics/local-{baseline,native}-vs-cpu.{md,json}`
- `artifacts/make-check.log`

The existing `input_3d_centre_bore_128_a2` local sanitizer gate is not a clean
oracle in this worktree: baseline reports 2,498 invalid reads and native reports
3,362, both beginning at `Numeric/Stencil.H:1592` from
`Operator::Elastic::Diagonal`. Both logs are retained; the canonical 2D case
is used for the feature-specific A100 sanitizer gate.

## NOVA/A100

Hardware: two NVIDIA A100-SXM4-80GB GPUs (sm_80), CUDA 12.8.1, driver
580.159.04. The timing binary is the strict/no-fast-math profile build. Each
timing sample completes ten steps, the arms alternate order, and the table
reports the median of three samples. No arena cap or global stream parameter
is present.

| case | one-step startup | baseline wall / step | native wall / step | wall delta | startup-adjusted delta | MLMG delta |
|---|---:|---:|---:|---:|---:|---:|
| 2D conservative | 2.1779 s | 5.7843 s / 0.5784 s | 4.0100 s / 0.4010 s | -30.67% | -49.20% | -52.12% |
| 3D section | 2.1451 s | 12.0226 s / 1.2023 s | 9.5451 s / 0.9545 s | -20.61% | -25.08% | -25.36% |

The startup-adjusted column subtracts the separate one-step wall calibration
from both ten-step medians. Multi-step external wall is authoritative. The
MLMG timer is retained as supporting evidence, but no-sync execution can defer
GPU completion and invalidate inner asynchronous region attribution.

The A100 profile demonstrates that attribution trap directly: TinyProfiler
assigns 5.544 s to 80,694 baseline `Fapply` calls but only 0.215 s to the same
80,694 native calls. The kernels did not become 25.8 times faster; their work
was charged at later synchronization points. End-to-end wall remained stable
across all repetitions: baseline 12.023/11.908/12.063 s and native
9.485/9.545/9.597 s.

The hardware delta is primarily a denominator effect. A quiet paired local 3D
run removed 0.118 s per MLMG call, while the A100 medians removed 0.133 s per
call. Baseline MLMG time was 3.294 s/call on the RTX A1000 versus 0.525 s/call
on the A100, so approximately the same host synchronization cost was 3.6% of
the former and 25.4% of the latter. The reported local median is additionally
weak: its paired wall deltas were -45.1%, +6.9%, and -3.3%, whereas all three
exclusive A100 deltas were between -19.9% and -21.1%.

Job `11771684` passed:

- strict baseline/native physics on one A100: all 27 observables PASS;
- strict baseline/native physics on two ranks/two A100s: all 27 observables
  PASS;
- enabled-path Compute Sanitizer memcheck: `ERROR SUMMARY: 0 errors`;
- 2D and 3D ten-step residual/iteration traces: byte-identical;
- all timing runs completed ten steps.

The first cross-rank diagnostic used AMReX's rank-dependent default BoxArray:
one rank produced one coarse box and two ranks produced two. The validation
extractor integrates nodal plotfile patches independently, so its L2 norms
included a different number of duplicated patch-interface nodes. That
comparison failed even for initial-condition fields while Linf extrema agreed;
it was an invalid decomposition-changing oracle, not a no-sync discrepancy.

Job `11771696` repeated all four strict runs with
`amr.max_grid_size=32`, fixing the BoxArray across rank counts. All four gates
then pass: baseline/native at one rank, baseline/native at two ranks,
one-rank/two-rank baseline, and one-rank/two-rank native. Cross-rank field
differences are roundoff-scale and solver counts match.

Evidence:

- `artifacts/nova-runtime/11771684/`
- `artifacts/nova-runtime/11771684/physics-{1,2}rank/compare.{md,json}`
- `artifacts/nova-runtime/11771684/sanitizer-summary.txt`
- `artifacts/nova-runtime/11771684/timing/`
- `artifacts/nova-fixed-grid-mpi/11771696/`
- `artifacts/nova-fixed-grid-mpi/11771696/physics-1rank-vs-2rank-{baseline,native}.{md,json}`

Build job `11771646` produced the exact strict binaries but stopped before
simulation because NOVA's Python module lacked PyYAML. Runtime jobs therefore
wrote plotfiles directly; metric extraction and comparison ran locally in the
existing `yt` validation environment. The binary SHA-256 values are:

- 2D:
  `b4e310d201f49fc675b5dee915c6748be5c82622e0d1c9d902becc3ef7ded2c7`
- 3D:
  `795f57d08193ec4c83215f00529e251804820dff876b62af70f8af1af9b765b1`

## Six requested checks

| requested check | result |
|---|---|
| Native AMReX API, no global stream knob | PASS |
| 2D / 3D / many-box without global stream config | PASS correctness; mixed performance |
| Exception/failure restoration | PASS |
| A100 Compute Sanitizer + strict physics | PASS (`11771684`) |
| Multi-rank GPU | PASS, including fixed-BoxArray cross-rank physics (`11771696`) |
| Arena limit only in local harness | PASS |

## Reproduction notes

- Local repository HEAD: `cc22ff9ee08fe7330ce6467f8a59bde68a7e92a4`
- CUDA source: `ext/AMReX-Codes/amrex`, AMReX 26.06.
- NOVA isolated worktree:
  `/work/brunnels/jackplum/alamo-mlmg-nosync-20260727`
- NOVA runtime jobs: `11771684` and `11771696` (both exit zero).
- The user's `src/Integrator/Flame.{cpp,H}` edits were neither edited nor
  reverted; the exact files were copied into the isolated A100 worktree so the
  remote binary tests the same current source state.

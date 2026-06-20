# Result: Task B — Phase 1 elastic disposition

## Summary (D1 verdict: GPU | CPU | HYBRID, with confidence)

**D1 verdict: CPU** (elastic stays CPU-resident). **Confidence: high.**

GPU elastic (strict/no-fast-math binary) diverges (SIGABRT) above 1024^2 finest
grid for a reason unrelated to `--use_fast_math` (the strict binary diverges
too, with the exact same character: immediate blow-up from the first MLMG
V-cycle of the first Newton iteration). At every resolution where GPU *does*
converge (512^2, 1024^2), a fully-subscribed CPU `np8` node is faster — 3.40x
faster at 512^2, 2.27x faster at 1024^2 — using matched inputs and matched
MLMG iteration counts (so the comparison is apples-to-apples, not an artifact
of different convergence paths). CPU is also the only path that completes at
2048^2 at all. Per the task's decision tree ("GPU slower than CPU node ⇒
CPU-resident"), and the GPU's complete failure to scale past 1024^2, the
verdict is CPU, not hybrid — there is no resolution regime where GPU elastic
is currently competitive, let alone superior.

## 1.1 convergence: resolution achieved, memory ceiling, converge vs diverge, iters vs CPU

- Resolutions tested on GPU strict (`bin/alamo_gpu_phase1-nofast`, a private
  copy of `bin/alamo_gpu-2d-nofast-cuda86-g++`): 512^2, 1024^2, 2048^2.
- 512^2: converges, 22+22 MLMG iterations (Newton iter 1, 2).
- 1024^2: converges, 16+16 MLMG iterations.
- 2048^2: **diverges, SIGABRT.** Reproduced fresh (`gpu_strict_2048_retest.log`):
  the very first MLMG V-cycle of Newton iteration 1 already shows
  `Fine resid/bnorm = 1.089` (no improvement over the initial residual) and the
  residual grows monotonically and explosively from there, aborting after 126
  iterations at `resid/bnorm = 1.84e20`. (An earlier run aborted after 83
  iterations at `1.04e20` — same failure mode, different iteration count,
  consistent with genuine numerical instability rather than a single
  deterministic bug.)
- **Memory ceiling was NOT the limiting factor.** The 2048^2 run executed all
  126 (failing) MLMG iterations rather than aborting on allocation; a
  successful 1024^2 run reported `Free GPU global memory: 500 MB` remaining.
  The practical ceiling for GPU-resident elastic in this configuration is
  therefore a **numerical/conditioning ceiling at 1024^2**, not a memory
  ceiling — 2048^2 would likely fit in the 8 GB budget if it converged.
- CPU `bin/alamo-2d-clang++` `np8` converges at all three resolutions tested
  (512^2, 1024^2, 2048^2 — 2048^2 confirmed previously in
  `benchmark/phase1_elastic_2048/cpu_np8.log`, 13+13 iterations), with
  iteration counts matching GPU exactly at 512^2/1024^2 (22+22, 16+16) —
  confirming GPU and CPU follow the identical numerical path when GPU
  converges at all, and that the 2048^2 GPU failure is GPU-specific.

## 1.2 branch taken and rationale

Branch taken: **"still diverges ⇒ it is NOT the flag."** The no-fast-math/
strict binary (`--fmad=false`, no `--use_fast_math`) reproduces the same
divergence character as the original fast-math finding at 2048^2. This
confirms (does not merely repeat) the prior conclusion in
`benchmark/PHASE1_ELASTIC_DISPOSITION.md`: the fast-math compiler flag is
**not** the root cause. Per the task's instruction, this is a diagnosis, not a
fix — device-elastic performance work is flagged as frozen above 1024^2 for
the lead, and the candidate root-cause investigation (operator/conditioning at
thin features, `psi_floor`, masked-operator behavior at 4 AMR levels / 8 MG
levels) is documented but **not pursued** in this task (no source edits, no
solver debugging — out of scope per the task's non-goals).

## 1.3 occupancy/registers for Fapply / Diagonal / Newton

**Blocked — not obtainable on this machine.** `ncu` (v2022.4.1.0,
`.local/nsight-compute/usr/bin/ncu`) connects to the target process but every
attempt (both `--set default`/`--set basic` profile sections and an explicit
single metric, `sm__warps_active.avg.pct_of_peak_sustained_active`) returns
`No metrics to collect found in sections.` / `No kernels were profiled.` Root
cause: `ERR_NVGPUCTRPERM` (no permission to access NVIDIA GPU performance
counters). Confirmed via `/proc/sys/kernel/perf_event_paranoid` = `4` (most
restrictive) and no passwordless sudo available to relax it or set
`NVreg_RestrictProfilingToAdminUsers=0`. This matches the task's documented
stop condition for ncu unavailability. **No occupancy or registers/thread data
exists for Fapply, Diagonal, or Newton::prepareForSolve from this
investigation** — this is an environment/permissions gap, not a finding about
the kernels, and should not be used as evidence for or against the GPU path.

## 1.4 fair GPU-vs-CPU-node benchmark (solve regions only)

Reran CPU `np8` at 512^2 and 1024^2 with identical overrides to the existing
GPU strict logs (previously the only CPU log was at 2048^2, a resolution GPU
cannot reach — not a matched comparison). All timings are the MLMG `Solve`
timer (solve region only) summed over Newton iterations 1+2 of the first
elastic solve:

| Finest grid | GPU strict (np1) | CPU np8 | CPU/GPU (CPU faster) | Iters match? |
| --- | ---: | ---: | ---: | --- |
| 512^2  | 2.624 s | 0.773 s | **3.40x** | yes (22+22 both) |
| 1024^2 | 3.109 s | 1.371 s | **2.27x** | yes (16+16 both) |
| 2048^2 | SIGABRT | 2.091 s | N/A — GPU incomplete | n/a |

Source logs: `benchmark/phase1_elastic_2048/{gpu_strict_512,gpu_strict_1024,
cpu_np8_512_retest,cpu_np8_1024_retest,gpu_strict_2048_retest,cpu_np8}.log`.

## Does the coarse 1.65x win survive the multi-core baseline?

**No.** The original "~1.65x GPU faster" figure (coarse 64^2-base case, from
prior memory `gpu_perf_nsys_findings.md`) was measured against a
single-process/under-subscribed CPU baseline. Against a fully-subscribed
`np8` CPU node at matched resolution and matched MLMG iteration count, CPU is
faster at every resolution where both paths complete (3.40x at 512^2, 2.27x at
1024^2), and CPU is the only path that completes at 2048^2. The coarse GPU win
does not survive a fair multi-core comparison — it was an artifact of an
unfair baseline.

## Files changed/created

- Created (private binary, per hard rule, copied before any other action):
  `bin/alamo_gpu_phase1-nofast`.
- Created: `input_phase1_elastic` (copy of shared `input`; shared `input` was
  never modified).
- Created/extended under `benchmark/phase1_elastic_2048/`:
  `gpu_strict_2048_retest.log` (+ `gpu_strict_2048_retest/` plotfile),
  `cpu_np8_1024_retest.log` (+ plotfile), `cpu_np8_512_retest.log` (+
  plotfile). Pre-existing files (`cpu_np8.log`, `gpu_strict.log`,
  `gpu_strict_512.log`, `gpu_strict_1024.log`, and their `*_plot` dirs) were
  read, not modified.
- Extended: `benchmark/PHASE1_ELASTIC_DISPOSITION.md` (report R1) with a new
  "R1 update — 2026-06-20, Phase 1 task B" section covering all four steps and
  the final D1 verdict.
- This file: `docs/agent_plans/20260620-gpu-phase1-and-lever25/results/B-RESULT.md`.
- No edits to `src/`, `configure`, `Makefile`, `.alamo`, the shared `input`, or
  any file owned by Track A. No builds were run. No git commit was made.

## Issues found

- ncu cannot collect any GPU hardware counters in this environment
  (`ERR_NVGPUCTRPERM`, `perf_event_paranoid=4`, no passwordless sudo) — this
  blocks step 1.3 entirely and should be flagged to whoever controls this
  machine if occupancy/register data is needed in the future (requires either
  root to set `NVreg_RestrictProfilingToAdminUsers=0` or running as root).
- The 2048^2 GPU divergence is non-deterministic in iteration count to abort
  (83 vs 126 iterations across two runs with identical inputs) even though
  both are on the same private binary — worth noting for whoever picks up the
  root-cause investigation, as it may indicate uninitialized state, race-y
  reduction order, or atomic-add nondeterminism in the device kernels, not
  purely a static conditioning problem.

## Deviations from task

- Task asked for a "≥1024^2" fair benchmark; delivered 512^2 AND 1024^2 (both
  resolutions where GPU converges) for a more complete picture, plus confirmed
  the 2048^2 asymmetry (CPU converges, GPU does not) rather than stopping at
  the minimum required resolution.
- Step 1.3 (ncu) could not produce any data due to a hard permissions wall in
  this environment; per the task's explicit stop condition this is recorded
  and the investigation proceeded with 1.4 and the verdict using the data that
  was available. No occupancy-based bias was applied to the D1 decision tree
  (which would have required real occupancy numbers); the verdict rests
  entirely on the convergence result (1.1/1.2) and the timing comparison
  (1.4), both of which independently point to CPU.

## Follow-up needed (e.g. operator root-cause if diverging)

- Root-cause the 2048^2+ GPU elastic divergence (operator/conditioning, not a
  build flag — confirmed by this task). Recommended path is already recorded
  in `benchmark/PHASE1_ELASTIC_DISPOSITION.md` under "Next investigation":
  compare CPU vs GPU operator inputs (`model_mf`, `psi_mf`, BC masks, RHS,
  grid metadata) at 2048^2 before `solver.solve`; isolate
  `Operator::Elastic`/MLMG from Flame model construction with a minimal
  elastic-only case; sweep `elastic.use_psi`, `psi_floor`, AMR depth one knob
  at a time; then inspect `Fapply`/`Diagonal`/Newton `prepareForSolve` for
  nodal indexing, ghost-fill, mask, or precision-sensitive GPU-vs-CPU
  differences. The iteration-count nondeterminism noted above (83 vs 126
  iterations to abort) should be investigated as part of this — it suggests
  looking at reduction/atomic ordering in the device kernels, not only static
  operator conditioning.
- Re-attempt ncu profiling (step 1.3) if/when GPU performance-counter
  permissions are available (root access or driver param change) — needed
  before any occupancy-based engineering decision about Fapply/Diagonal/
  Newton kernel design can be made.
- Per the D1 verdict, do not invest further GPU elastic performance work until
  both the divergence is fixed and a GPU implementation is shown to beat a
  fully-subscribed CPU `np8` node — neither holds today.

# Phase 1 Elastic Disposition

> **✅ The 2048² GPU divergence reported here is SOLVED (2026-06-25).** Root cause:
> a GPU cross-stream use-after-free race on a per-box temporary in
> `Operator<Grid::Node>::interpolation()`; fixed with one line (`tmpfab.elixir()`).
> It was never a fast-math, conditioning, or CPU-vs-GPU-perf question — it was a
> multi-box-only transfer-temp lifetime bug. The "CPU-resident elastic" disposition
> below is **obsolete**; elastic runs correctly on GPU multi-box. See
> `benchmark/GPU_BRANCH_GUIDE.md` (D1 section) and
> `benchmark/elastic_sensitivity_20260621/GPU_ELASTIC_DEBUG_PLAN.md`. The perf
> numbers below (CPU np8 2.27–3.40× faster) remain valid data for the separate
> "is device-elastic worth it" question.

Date: 2026-06-20

Roadmap: `/home/jackplum/Desktop/GPU-OPT-ROADMAP.txt`

## Purpose

Execute D1 from the GPU optimization roadmap: test whether the strict/no-fast CUDA elastic solve converges at the 2048^2 finest-grid case and compare against the local CPU best-case baseline.

## Case

Base input: `input`

Overrides used for all runs unless noted:

```bash
max_step=1
stop_time=1e99_s
amr.max_level=3
elastic.interval=1
elastic.tstart=0.0
elastic.solver.verbose=4
elastic.print_model=1
amr.plot_int=-1
amr.thermo.plot_int=1
amr.thermo.int=1
```

Resolution mapping:

- 512^2 finest: `amr.n_cell="64 64 64"`, `max_level=3`
- 1024^2 finest: `amr.n_cell="128 128 128"`, `max_level=3`
- 2048^2 finest: `amr.n_cell="256 256 256"`, `max_level=3`

Logs and plotfiles are under `benchmark/phase1_elastic_2048/`.

## Results

| Run | Binary / ranks | Finest grid | Result | Newton | MLMG result |
| --- | --- | ---: | --- | ---: | --- |
| GPU strict | `bin/alamo_gpu-2d-nofast-cuda86-g++`, `-np 1` | 512^2 | converged | 2 | 22 + 22 iterations, final ratios `8.38e-9`, `8.22e-9` |
| GPU strict | `bin/alamo_gpu-2d-nofast-cuda86-g++`, `-np 1` | 1024^2 | converged | 2 | 16 + 16 iterations, final ratios `4.66e-9`, `4.60e-9` |
| GPU strict | `bin/alamo_gpu-2d-nofast-cuda86-g++`, `-np 1` | 2048^2 | failed | first Newton solve | MLMG diverged and aborted after 83 iterations; final `resid/bnorm = 1.044e20` |
| CPU | `bin/alamo-2d-clang++`, `-np 8` | 2048^2 | converged | 2 | 13 + 13 iterations, final ratios `4.04e-9`, `3.90e-9` |

## Interpretation

The 2048^2 divergence is not explained by CUDA fast-math. The strict/no-fast CUDA binary still diverges, while the CPU `-np 8` solve converges on the same input and resolution.

This points to a GPU elastic path bug or GPU-specific operator/MLMG behavior that appears only at the 2048^2 finest-grid case. The strict GPU path is healthy at 512^2 and 1024^2 for this same test setup, so the failure is scale-sensitive rather than universally broken.

## D1 verdict

Elastic disposition for now: **CPU-resident / freeze device-elastic performance work**.

Per the roadmap, do not proceed with GPU elastic performance optimization until the 2048^2 strict-GPU divergence is root-caused. Phase-field GPU optimization can proceed independently, but any coupled workflow should assume CPU elastic or hybrid CPU-elastic/GPU-phase-field until this is fixed.

## Next investigation

Recommended root-cause path:

1. Compare CPU vs GPU operator inputs at 2048^2: `model_mf`, `psi_mf`, BC masks, RHS, and geometry/grid metadata before `solver.solve`.
2. Run a minimal 2048^2 elastic-only/patch-style case to separate Flame model construction from `Operator::Elastic`/MLMG behavior.
3. Test GPU 2048^2 with simplified knobs one at a time: `elastic.use_psi=0`, higher `psi_floor`, uniform model, reduced AMR depth, and single-level 2048^2 if memory allows.
4. Once the first CPU/GPU divergence point is isolated, inspect `Operator::Elastic::Fapply`, `Diagonal`, and Newton `prepareForSolve` data paths for nodal indexing, ghost fill, mask, or precision-sensitive differences.

---

## R1 update — 2026-06-20, Phase 1 task B (D1 verdict finalized)

This section extends the report above with: (1) a 2048^2 reproduction using a
**private** copy of the strict binary (`bin/alamo_gpu_phase1-nofast`, isolated from
a concurrent Track A rebuild), (2) the D1 branch decision, (3) an ncu
occupancy/registers attempt, and (4) a fair np8-CPU-vs-GPU benchmark at matched
resolutions. Inputs used: `input_phase1_elastic` (a copy of the shared `input`,
never modified in place). All new artifacts are under
`benchmark/phase1_elastic_2048/*_retest*`.

### 1.1 — 2048^2 reproduction (strict/no-fast binary)

Reran the exact 2048^2 case (`amr.n_cell="256 256 256"`, `max_level=3`,
`elastic.tstart=0.0`) on the private nofast copy:

- Result: **SIGABRT, MLMG diverges**, confirming the original finding
  (`gpu_strict_2048_retest.log`). This run diverged even *worse* than the first
  attempt — 126 iterations to abort vs. 83 previously — consistent with a
  genuine numerical instability (sensitive to minor state/scheduling
  differences), not a deterministic single failure mode.
- Critically: **Newton Iteration 1's very first MLMG iteration already has
  `Fine resid/bnorm = 1.089`** (no improvement at all over the initial
  residual), and it is monotonically unstable from there — this is not "slow
  convergence that eventually gives up," it is immediate divergence at the
  first V-cycle. 4 AMR levels, 8 MG levels on the coarsest AMR level (vs. 7 MG
  levels at 1024^2).
- Not OOM: GPU memory was never the limiting factor (`Free GPU global memory:
  500 MB` reported at the end of a *successful* 1024^2 elastic run, and the
  2048^2 attempt also ran to completion of its (failing) iterations rather than
  aborting on allocation).
- Memory ceiling for this elastic configuration was not hit at 2048^2 — the
  failure is purely numerical, so the practical resolution ceiling for
  GPU-resident elastic in this configuration is **1024^2** (last resolution
  that both fits and converges).

### 1.2 — D1 branch taken

**Branch: "still diverges ⇒ NOT the fast-math flag."** The strict/no-fast binary
(`--fmad=false`, no `--use_fast_math`) diverges identically in character to the
original fast-math finding. This reconfirms and finalizes the original
conclusion in this file: fast-math is **not** the root cause of the 2048^2
failure. Per the task's decision tree, this **freezes device-elastic
performance work** above 1024^2 and flags an operator/conditioning
investigation for the lead (see "Next investigation" above, still applicable
and not pursued further here — out of scope for this task, which is
measurement/disposition only, not source debugging).

### 1.3 — ncu occupancy/registers attempt

Attempted via `.local/nsight-compute/usr/bin/ncu` (v2022.4.1.0) against the
private nofast binary, both with `--set default`/`--set basic` and with an
explicit single metric (`sm__warps_active.avg.pct_of_peak_sustained_active`),
targeting `.*Elastic.*Fapply.*`, `.*Elastic.*Diagonal.*`, and
`.*Newton.*prepareForSolve.*` kernel regexes at 1024^2 (the largest converging
resolution).

**Result: blocked.** `ncu` connects to the process (`==PROF== Connected to
process ...`) but every attempt returns `==WARNING== No metrics to collect
found in sections.` / `==WARNING== No kernels were profiled.` Root cause:
`ERR_NVGPUCTRPERM` — "The user does not have permission to access NVIDIA GPU
Performance Counters on the target device." Confirmed via
`/proc/sys/kernel/perf_event_paranoid` = `4` (most restrictive) and no
passwordless sudo (`sudo -n true` fails) to set
`NVreg_RestrictProfilingToAdminUsers=0` or relax `perf_event_paranoid`. This
matches the task's documented stop condition ("ncu unavailable/permission-
blocked → record and proceed with what you have"). **No occupancy or
register/thread numbers were obtainable on this machine for Fapply,
Diagonal, or Newton::prepareForSolve.** This is an environment limitation, not
evidence about the kernels themselves — do not infer occupancy from it either
way.

### 1.4 — Fair GPU-vs-CPU-node (np8) benchmark, matched resolutions

The existing logs compared GPU strict at 512^2/1024^2 against a CPU run that
was only captured at **2048^2** — not a matched-resolution comparison. Reran
CPU `np8` at 512^2 and 1024^2 (`cpu_np8_512_retest.log`,
`cpu_np8_1024_retest.log`) with identical overrides, so each resolution now has
a same-iteration-count (16+16 at 1024^2, 22+22 at 512^2), same-input GPU vs CPU
pair. Times are MLMG `Solve` timer sums (Newton iteration 1 + 2; first elastic
solve only, solve region only):

| Finest grid | GPU strict (np1, `Solve` sum) | CPU `np8` (`Solve` sum) | CPU/GPU speedup (CPU faster) |
| --- | ---: | ---: | ---: |
| 512^2  | 2.624 s | 0.773 s | **3.40×** |
| 1024^2 | 3.109 s | 1.371 s | **2.27×** |
| 2048^2 | SIGABRT (diverges) | 1.048+1.043 = 2.091 s, converges | N/A (GPU cannot complete) |

(MLMG iteration counts matched exactly between GPU and CPU at each resolution,
confirming this is an apples-to-apples timing comparison, not an artifact of
different convergence paths.)

### Does the coarse 1.65x win survive the multi-core baseline?

**No.** The originally reported "~1.65x GPU faster" result (memory:
`gpu_perf_nsys_findings.md`, coarse 64^2-base case) was against a
single-process/single-core CPU baseline. Once compared against a **fully
subscribed `np8` CPU node** at matched resolution and matched iteration count,
**CPU is faster at every resolution tested** — 3.40x faster at 512^2 and 2.27x
faster at 1024^2 — and CPU is the *only* path that completes at all at 2048^2.
The GPU elastic win is an artifact of comparing against an under-subscribed
CPU baseline; it does not survive a fair multi-core comparison.

### D1 verdict (final): **CPU**

- GPU is slower than the CPU `np8` node at every resolution that runs on both
  (2.27x–3.40x slower), and GPU elastic cannot complete at all above 1024^2 due
  to a non-fast-math-related divergence.
- Per the task's decision tree ("GPU slower than CPU node ⇒ CPU-resident"),
  elastic should remain CPU-resident. This is a high-confidence verdict: it is
  based on a matched-resolution, matched-iteration-count, solve-region-only
  benchmark, not an extrapolation.
- Device-elastic performance work should stay frozen until (a) the 2048^2+
  divergence is root-caused (operator/conditioning, not a build flag) and (b)
  a GPU implementation is shown to beat a fully-subscribed CPU node — neither
  condition currently holds.

**Relationship to `docs/gpu_elastic_device_port_plan.md`:** that plan
("Option B", full device-native elastic) is a controlled experimental branch
exploring whether de-virtualization fixes the device crash, not a
reopening of this D1 verdict. Progress there (e.g. completing a step-60 run)
shows the solve *runs*, not that it *wins* — it must still clear the three
hard gates that plan now states (no-fast-math convergence at high res, a
fair N-rank CPU comparison, `ncu` data on the actual target GPU) before D1
flips. Until then, CPU-resident elastic remains the standing, mainline
disposition.

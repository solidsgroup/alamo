# GPU Acceleration Roadmap v2 — Measurement-Driven Optimization (alpha-1.0 → beta)

Branch: `chamber-gpu` (never merged — see `benchmark/GPU_BRANCH_GUIDE.md`).
Supersedes-as-forward-plan: the original `~/Desktop/GPU-OPT-ROADMAP.txt` (Phases
0–5), which is **complete** — all 5 definition-of-done items are DONE
(`benchmark/PHASE5_BRANCH_DONE.md`). This v2 roadmap is the *next* arc:
alpha-1.0 → beta.

Created: 2026-06-25. Updated: 2026-06-28.
Status: **Phase A measured (A1+A2) · B1 gate PASSED · Phase C ACTIVE.**

> ### 2026-06-28 — Phase A is measured. Full analysis: `benchmark/PHASE_A_FINDINGS.md`
> The combined flame+elastic A100 profiling (nsys 11306663/11306722 + diag 11310480
> + ncu 11318197) landed. Headlines that drive everything below:
> - **`Operator::Elastic::Fapply` = 74.6% of GPU kernel time** (58% exclusive of the
>   elastic solve). The elastic MLMG solve = **95.7% of evolve time.** Flame = **0.20%.**
> - **GPU is 89.3% busy ⇒ compute/memory-bound**, not launch-bound.
> - The solve **converges cleanly in ~10 V-cycles** (not grinding); the Krylov bottom
>   solver is only ~7% of solve wall time. Cost is the **fine-level Fapply applies**
>   (255 registers/thread ⇒ ~12.5% occupancy; ~2.5 KB/node uncoalesced AoS stiffness).
> - **Combined GPU speedup measured: 22.9× per step (256³, 1 A100 vs 64-rank CPU)** —
>   *better* than flame-only (13.1×). The Amdahl-cap fear is refuted.
> - ⇒ **B1 gate PASSED (elastic ≈95% of wall).** Phase C kernel work is justified.
>   **100% of the available time savings are in the elastic solve;** flame needs no
>   GPU work. Levers + A/B harness: `PHASE_A_FINDINGS.md` §8,
>   `input_3d_centre_bore_256_a2_tuned`, `benchmark/phase_c_elastic_ab.sh`.

---

## 1. Why a v2 roadmap (the inversion)

The Phase 0–5 roadmap was **fix-and-disposition driven**: does the GPU build
run, does elastic diverge, is device-elastic worth it, does the CPU suite stay
green. It succeeded. But its recurring failure mode is one sentence repeated in
every report: *"launches/step, sync fraction, and achieved occupancy were not
collected — pending `ncu`."* Optimization and disposition calls were made
largely **blind to the device counters** (`ncu` was blocked the entire project
by `ERR_NVGPUCTRPERM` on the local A1000; `PHASE3_R3_crossover.md` "Standing
metric set" lists every counter row as *pending*).

v2 inverts that. The governing rule is:

> **No kernel-level optimization without a counter that justifies it.**

Phase A is instrumentation. Phases B–D spend only the levers the data proves
matter.

---

## 2. Current state — the alpha-1.0 baseline

What is **proven** (cite, don't re-litigate):

- **Phase-field wins ~10× on GPU in the saturating 3D regime.** Single A100 vs
  full 64-rank CPU node: **9.6× at 128³, 13.1× at 256³**, advantage growing
  with size (sublinear GPU cost-per-cell). HIGH confidence; 128³/1-GPU number
  reproduced within 6% across two NOVA sessions. `PHASE3_R3_crossover.md`.
- **Elastic now runs correctly on GPU multi-box.** The entire prior divergence
  history collapsed to one bug — a cross-stream use-after-free on the per-box
  `tmpfab` in `Operator<Grid::Node>::interpolation()`
  (`src/Operator/Operator.cpp:728`, one-line `elixir()` fix, commit
  `c00f69086`). Verified end-to-end: box sweep + real 3-material operator +
  300-step flame+elastic chamber sim all converge ~1e-9. `GPU_BRANCH_GUIDE.md`
  (D1 section). **NB (2026-06-26): that verification was all 2D.** The **3D**
  elastic device path was broken by a separate bug (`F.inverse().transpose()`
  faulting on device, `NeoHookean.H`) until fixed this date — it had never been
  run, since every 3D crossover used `elastic.type=disable`. See
  `docs/llm/changelog/2026-06-26-3d-elastic-gpu-fix.md`.
- **Branch DoD is 5/5 DONE**, CPU regression suite green (118 run, 92 verified,
  0 failed). `PHASE5_BRANCH_DONE.md`.

What is **unmeasured** (the gap this roadmap exists to close):

- **Combined flame+elastic performance on A100 is unknown.** The 9.6–13.1×
  crossover was measured with `elastic.type = disable`. The elastic fix is
  *correctness*-verified, never *performance*-benchmarked on A100. The only
  elastic perf data that exists (CPU np8 2.27–3.40× faster, `PHASE1_ELASTIC_DISPOSITION.md`)
  is **2D, on the 50W-capped shared local A1000, in the pre-saturation regime**
  — wrong hardware and wrong regime. It does not transfer to A100 3D.
- **No `ncu` counters anywhere** (launches/step, occupancy, sync fraction,
  registers/thread) on the real target hardware.
- **AMR depth is unvalidated on GPU** as a tunable — only the strong claim that
  deep subcycling AMR is launch-poison (see §4).

**Bottom line for stakeholders:** *phase-field GPU acceleration is proven (~10×);
combined flame+elastic GPU acceleration is unmeasured.* Any single
"combined speedup" number quoted before Phase A completes is a guess.

### Alpha-1.0 cut

Tag the current tip of `chamber-gpu` as the alpha-1.0 baseline-of-record
(annotated tag, **not** a merge — per the never-merge policy). This roadmap's
Phase A re-baselines *with elastic enabled*; alpha-1.0 is the frozen reference
those measurements are compared against. See task **A0**.

---

## 3. Guiding principles

1. **Measure before optimize.** A counter justifies every kernel change (§1).
2. **Two builds, two purposes.** Perf headlines use the fast-math binary;
   correctness claims use *only* the strict/no-fast-math binary
   (`GPU_BRANCH_GUIDE.md` build matrix). Never mix.
3. **Single-GPU is the supported shape** at 128³/256³. Multi-GPU is a measured
   negative result (halo-bound), parked until a domain large enough to amortize
   halos exists (§ Phase D).
4. **Capture the standing metric set every run** (§5) — no more one-off
   throughput-only logs.
5. **Combined physics is the unit of measurement** from here on. Phase-field-only
   and elastic-only numbers are diagnostic sub-measurements, not the headline.

---

## 4. On AMR (the standing hypothesis Phase B tests)

The current data says **deep subcycling AMR is a net negative on GPU**, and the
mechanism is understood, not guessed:

- Wide-shallow (`max_level=1`) beats deep AMR (`max_level=5`) at the same finest
  resolution by **~9.8× per step**; nsys shows deep AMR issues **646k vs 96k
  `cudaLaunchKernel`/step** (6.7×) with the *same* ~3µs avg kernel — the GPU is
  **launch/latency-bound, not compute-bound**. Deep AMR buys 6.7× more *tiny*
  kernels, not bigger ones. (`gpu_perf_nsys_findings.md`.)
- Launches scale as `nsubsteps^level`; nsubsteps 2→1 cut launches 5.7×.

So the GPU win regime is wide-shallow / uniform. The open question Phase B
answers is whether the chamber physics *needs* AMR's resolution concentration,
and whether on GPU a uniform-fine grid is simply cheaper than an adaptive
hierarchy (the CPU tradeoff **inverts** on GPU). Two unexplored levers:
**non-subcycling AMR** (removes the `nsubsteps^level` multiplier) and
**1-level + large blocking factor**. Prediction to be tested: **0–1 AMR levels
on GPU.**

---

## 5. The standing metric set (capture identically every run)

Aligned with `PHASE3_R3_crossover.md` and `benchmark/PERF_TRACKING.md`. Every
benchmark run in this roadmap emits these as a saved artifact, not a one-off log:

| Metric | Source | Tool |
| --- | --- | --- |
| wall/step, GPU vs CPU node | run log / `time` wrapper | always |
| `launches_per_step` (`cudaLaunchKernel` calls / steps) | `cuda_api_sum` | nsys |
| `kernel_avg_us` (Σ kernel time / Σ instances) | `cuda_gpu_kern_sum` | nsys |
| `sync_frac` (`cudaStreamSynchronize` / total CUDA-API time) | `cuda_api_sum` | nsys |
| achieved occupancy + registers/thread (per kernel family) | sections | **ncu** |
| **elastic vs phase-field wall-time split** | TinyProfiler regions | always |
| golden-compare residual vs CPU (no-fast-math) | `golden_compare_flame.sh` | strict build |

Existing harness: `benchmark/perf_regression_track.py` already computes
`launches_per_step`/`kernel_avg_us`/`sync_frac` from an nsys capture and appends
to `benchmark/perf_regression.csv` per commit (`PERF_TRACKING.md`). The
elastic/phase-field split and the ncu rows are the new additions.

---

## 6. Phases & tasks

Four phases. Each task has: **Goal · Steps · Done-when · Depends · Artifact.**
IDs are stable references (A0, A1, …).

### Phase A — Instrument & baseline

The gate. Nothing in B–D starts until A1+A2 produce data.

---

**A0 — Cut the alpha-1.0 baseline-of-record**
- **Goal:** freeze the current branch tip as the named reference point.
- **Steps:** annotated tag `chamber-gpu-alpha1` on the current tip (elastic-fix
  commit `c00f69086` or later); record the binary build commands and the exact
  input configs (`input_3d_flame_{128,256,512}`, `input_base` chamber) used as
  references. No merge.
- **Done-when:** tag exists; `benchmark/ALPHA1_BASELINE.md` records commit SHA,
  build flags, and the canonical input set.
- **Depends:** none.
- **Artifact:** git tag + `benchmark/ALPHA1_BASELINE.md`.

**A1 — Unblock device counters on NOVA (ncu + nsys)**
- **STATUS (2026-06-28): MOSTLY DONE.** `nsys` fully captured + analyzed
  (`profiles_alpha1/nsys_11306663`, kernel taxonomy + 89.3% busy fraction in
  `PHASE_A_FINDINGS.md`). `ncu` report `ncu_11318197/elastic_kernels.ncu-rep`
  (183 MB) was produced on NOVA but the local ncu (2022.4.1) can't open it
  (newer format) — the per-kernel Speed-of-Light (Fapply memory-vs-compute %)
  still needs a NOVA-side re-export or `ncu-ui` 2025.x. Registers/thread (255 for
  Fapply) came from the nsys kernel records instead. **Remaining:** the Fapply SoL.
- **Goal:** obtain the counters that were blocked the entire Phase 0–5 project.
- **Steps:** run `benchmark/g0_ncu_capture.sh` (prepared command) on a NOVA node
  where GPU performance counters are permitted; confirm `ncu` returns real
  section data (not `ERR_NVGPUCTRPERM`). Reuse the existing diagnostic harness
  (`benchmark/nova_flame_gpu_3d_diag.slurm` + `benchmark/phase3_nova_diag.sh`)
  as the submission vehicle; add an `ncu`/`nsys` profiling mode if the diag
  script only does `CUDA_LAUNCH_BLOCKING`/`compute-sanitizer` today.
- **Done-when:** one `ncu` report on NOVA with achieved occupancy +
  registers/thread for the Flame phase/thermal kernels and the elastic
  `Fapply`/`Diagonal`/`Newton::prepareForSolve` families; one `nsys` report with
  launches/step + sync fraction. Both archived.
- **Depends:** none (parallel with A0).
- **Artifact:** `benchmark/profiles_alpha1/{ncu,nsys}_*` + a short
  `PROFILE_NOTES.md` reading the numbers.

**A2 — Combined flame+elastic A100 baseline-of-record (THE key measurement)**
- **STATUS (2026-06-28): DONE (headline) — see `benchmark/PHASE_A_FINDINGS.md`.**
  Combined GPU 256³ = **22.9× per step** vs 64-rank CPU; elastic fraction
  `f ≈ 0.95` of GPU wall (Fapply 74.6% of kernel time). Stability clean (no
  NaN/diverge; physics matches CPU modulo grid-resolution). **Caveats for a
  publishable row:** (1) all four production a2 runs were SLURM-preempted/timed-out
  so the breakdown comes from the `interval=1` diag + nsys traces, not a completed
  run — a non-preemptible reservation is still wanted for a clean wall total;
  (2) the CPU baseline has a 32-idle-rank load-imbalance bug (`Operator.cpp:459`)
  that flatters the GPU ratio. Neither changes the optimization target.
- **Goal:** the single number that turns "unmeasured" into a real combined
  speedup, plus the elastic wall-time fraction `f`.
- **Steps:** re-run the Phase-3 crossover matrix **with elastic enabled**
  (`elastic.type` set to the real solver instead of `disable`) at 128³ and 256³,
  single A100, against the 64-rank CPU node baseline. Reuse
  `benchmark/phase3_scaling_sweep.sh` + `nova_flame_gpu_3d.slurm` /
  `nova_flame_cpu_3d.slurm`. Capture TinyProfiler region times so the
  elastic-solve vs phase-field-advance split is explicit. Confirm stability for
  ≥10 elastic solves (no NaN/diverge), consistent with the local 300-step
  validation.
- **Done-when:** a results table giving combined GPU/CPU ratio at 128³/256³ **and**
  the elastic share of GPU wall time; stability confirmed over the run (not just
  step 1).
- **Depends:** A0 (reference), benefits from A1 (counters) but does not block on
  it for the wall-clock headline.
- **Artifact:** `benchmark/PHASE_A2_combined_crossover.md` (parallels
  `PHASE3_R3_crossover.md`).
- **⚠️ Prerequisite hit + cleared (2026-06-26):** the very first attempt to run 3D
  elastic on GPU crashed immediately (`CUDA error 719` on the first elastic solve).
  Root cause: a chained Eigen expression `F.inverse().transpose()` faults on device;
  two sites in `src/Model/Solid/Finite/NeoHookean.H` (3D `DW`/`DDW`). Fixed by
  materializing the inverse first. This was the **first time the 3D combined device
  path was ever exercised** — exactly the "unmeasured" gap this task exists to close.
  Full writeup: `docs/llm/changelog/2026-06-26-3d-elastic-gpu-fix.md`. Before running
  A2 on NOVA, rebuild the **3D strict/nofast** binary with the fix (only the 3D fast
  binary was rebuilt locally) and re-confirm ≥10 stable solves.

**A3 — Lock the standing metric set + wire perf-regression into GPU CI**
- **Goal:** make §5 a per-run artifact, not a manual ritual; catch regressions
  automatically.
- **Steps:** extend `benchmark/perf_regression_track.py` (or a sibling) to also
  record the elastic/phase-field split; add a *gated, non-blocking* GPU job to
  `.github/workflows/chamber-gpu-correctness.yml` per the "Future CI wiring"
  block in `PERF_TRACKING.md` (`runs-on: [self-hosted, cuda]`), so it runs
  `--capture --compare` when a CUDA runner is attached and prints `skipped: no
  GPU` otherwise.
- **Done-when:** a perf-regression row is produced by CI on a GPU runner (or the
  job is in place and verified `skipped`-safe on CPU-only); the CSV schema
  includes the elastic split.
- **Depends:** A1 (defines what's captured), A2 (first real combined row).
- **Artifact:** updated `perf_regression_track.py`, workflow job, first CI row.

**A4 — Sibling GPU-only-fault audit (correctness insurance)**
- **Goal:** ensure the known device-only faults (the `interpolation()` elixir race
  **and** the `F.inverse().transpose()` expression-template fault) were the only
  instances of their patterns. Two distinct "works-on-CPU, faults-on-GPU" bug
  classes have now bitten elastic — treat the audit as covering both.
- **Steps:**
  1. *Elixir race:* grep ALAMO custom GPU `MFIter` loops for a temporary
     `FArrayBox`/`BaseFab` (`.resize(` / `FArrayBox ` inside `for (MFIter`) used in
     a device kernel without an `elixir()`; for each hit, decide elixir-needed or
     safe-by-construction (writes its own output directly, like `restriction()`).
     Add `elixir()` where the temp feeds an async device kernel.
  2. *Chained Eigen expressions on device:* `grep -rn "\.inverse()\.\|\.transpose()\.\|\.cross()\.\|\.eval()\." src --include=*.H --include=*.cpp`
     and inspect any in a `AMREX_GPU_*DEVICE` path. Nested Eigen expression templates
     (e.g. `Transpose<Inverse<Matrix3>>`) can fault on device even when each piece is
     device-safe alone. Fix by forcing evaluation to a concrete matrix first (two
     statements, or `.eval()`). As of 2026-06-26 only the two `NeoHookean.H` sites
     existed and are fixed.
- **Done-when:** an audit table of every in-loop device temp **and** every chained
  Eigen expression in a device path, each with a verdict; any fix verified by
  re-running the forced-multi-box elastic case (2D) and a 3D elastic smoke run.
- **Depends:** none.
- **Artifact:** `benchmark/elixir_race_audit.md` (rename/extend to cover both classes).
- **See:** `docs/llm/changelog/2026-06-26-3d-elastic-gpu-fix.md`.

### Phase B — Characterize the combined workload

Spends A1/A2 data to answer "where does the time go and what regime is optimal."

---

**B1 — Elastic time-fraction & kernel breakdown**
- **STATUS (2026-06-28): DONE — GATE PASSED.** `PHASE_A_FINDINGS.md` §2 is the
  breakdown: elastic = **≈95% of GPU wall** (Fapply 74.6% of kernel time), flame
  0.20%. Far above the 15% threshold ⇒ **Phase C is justified and unblocked.**
  The Amdahl sensitivity is moot — elastic is the workload, not a slice.
- **Goal:** decompose A2's combined number — is elastic a small slice or the
  bottleneck? This decides whether Phase C happens at all.
- **Steps:** from A2's TinyProfiler split and A1's `ncu`, attribute GPU wall time
  across phase-field advance, elastic MLMG (Fapply/Diagonal/smoother/bottom),
  FillBoundary/halo, and interp. Compute the Amdahl combined-speedup sensitivity
  to elastic's `f`.
- **Done-when:** a breakdown chart/table + a stated verdict: elastic is/ isn't a
  ≥X% wall-time fraction (X to be set, suggest 15%) → Phase C gate.
- **Depends:** A1, A2.
- **Artifact:** `benchmark/PHASE_B1_workload_breakdown.md`.

**B2 — AMR-depth sweep (the §4 hypothesis test)**
- **Goal:** find the GPU-optimal AMR configuration for the chamber problem at
  matched effective resolution **and matched accuracy**.
- **Steps:** sweep uniform-fine vs `max_level=1` vs deep AMR at fixed finest
  resolution; include the two unexplored levers — **non-subcycling AMR** (single
  dt across levels) and **1-level + large blocking factor**. Verify physical
  accuracy parity (burn-front position, mdot, pressure trajectory) so a faster
  config isn't just a coarser one. Capture launches/step for each.
- **Done-when:** a table of {config → wall/step, launches/step, accuracy-delta};
  a recommended AMR depth for the GPU chamber config.
- **Depends:** A1, A2.
- **Artifact:** `benchmark/PHASE_B2_amr_depth.md`.

**B3 — Resolution sweep & fair CPU baselines**
- **Goal:** map where combined speedup peaks; close the 512³ CPU-baseline gap if
  cheap.
- **Steps:** combined flame+elastic at 128³/256³/512³ single A100; produce a
  *valid* CPU baseline at each (fix the `--mem=32G` hardcode in
  `nova_flame_cpu_3d.slurm` that OOM'd 512³, only if a 512³ CPU point is wanted
  — otherwise leave 512³ GPU-only and note it). Plot speedup vs size.
- **Done-when:** speedup-vs-size curve for combined physics with the regime of
  peak advantage identified.
- **Depends:** A2; B2 (use the chosen AMR config).
- **Artifact:** `benchmark/PHASE_B3_resolution_sweep.md`.

### Phase C — Targeted, counter-justified kernel work

**Gate SATISFIED (B1, 2026-06-28): elastic ≈95% of wall.** All of Phase C is now
in-scope. Prioritized lever list + expected gains: `PHASE_A_FINDINGS.md` §8. Order
of attack: **C0 (cheap input levers, ship now) → C1 + C4 (kernel + reuse) →
C2/C3 from the re-measured counter.** Never land a kernel change without an A100
before/after.

---

**C0 — Elastic solver input tuning** *(no source change — DO FIRST)* ✅ staged
- **Goal:** cut the Fapply **call count** with zero risk before touching kernels.
- **Steps:** the a2 baseline left `tol_rel/abs` at the 1e-8 default and `bottom_*`
  unset (BiCGStab inherits 1e-8). Ship `input_3d_centre_bore_256_a2_tuned`
  (tol 1e-6, `bottom_tol_rel=1e-4`/`bottom_max_iter=50`, `interval=100`) and run
  `benchmark/phase_c_elastic_ab.sh` to attribute each lever vs the baseline at
  256³/A100, gated on stress-field parity (`compare_thermo.py`).
- **Done-when:** A/B table {baseline, tol, bottom, interval, all} → per-solve wall
  + V-cycle count + stress delta; expected ~3–5× on the elastic solve.
- **Depends:** none (executable now). **Artifact:** `benchmark/PHASE_C0_input_ab.md`.

**C1 — Elastic `Fapply` register pressure / occupancy** *(the structural win)*
- **STATUS (2026-06-28): SOURCE STARTED on branch `chamber-gpu-elastic-opt`** (isolated
  worktree `/home/jackplum/Projects/alamo-elastic-opt`; not in main, not committed). Two
  **bit-identical** edits to `Fapply` landed and compile GPU-3D/sm_86 clean: (b) grad(C)
  accumulated one direction at a time so only **one** `Matrix4` derivative temp is live
  (was `Cgrad1/2/3` = 3 live = 135 doubles in 3D); plus a boundary-only **`sig`-sink** that
  drops a `Matrix4·Matrix` product + a live `Set::Matrix` from every interior node (≈90% of
  Fapply time). (a) the per-kernel `__launch_bounds__` register cap is **deferred** (can
  regress; gate on the A1 ncu SoL). **Still owed: A100 before/after + strict golden compare.**
  Detail: `benchmark/PHASE_C1_fapply_occupancy.md`.
- **Goal:** raise the **~12.5%** occupancy of the kernel that is 74.6% of GPU time.
- **MEASURED (2026-06-28):** runtime register count is **255/thread** (the CUDA
  cap — the compiler is spilling), not the 87–101 in the G0 static dump ⇒ ~12.5%
  theoretical occupancy (1 block/SM). Root cause located: the `!m_uniform` branch
  of `Fapply` (`src/Operator/Elastic.cpp:615-624`) holds **4 live `Matrix4` temps**
  (DDW + Cgrad1/2/3 = 180 doubles).
- **Steps:** (a) `__launch_bounds__`/register cap on the Fapply launch; (b)
  restructure the grad(C) accumulation to not materialize 3 full `Cgrad` `Matrix4`
  temps — accumulate `(Cgrad_d·gradu).col(d)` one dimension at a time. Measure
  achieved occupancy (needs the A1 ncu SoL re-export) + wall/step each change.
- **Done-when:** occupancy up + wall/step down, strict-build golden compare clean,
  or a documented "no win" with the counter. **Depends:** A1 SoL (nice-to-have).
- **Artifact:** `benchmark/PHASE_C1_fapply_occupancy.md`.

**C2 — Fapply memory layout (AoS → SoA) + uniform-node fast-path** *(biggest single
kernel win; gate on the A1 ncu SoL)*
- **Goal:** kill the ~2.5 KB/node of **uncoalesced** stiffness traffic in the
  fine-level Fapply (`Array4<Matrix4>` is array-of-structs, 360 B/node stride).
- **Steps:** (a) precompute a thin material-interface mask so the >90% of
  interior-uniform nodes (3 piecewise-constant materials) skip the 6 neighbor
  `Matrix4` loads + the grad(C) branch entirely; (b) if the ncu SoL confirms
  memory-bound, store the 45 Matrix4 components as planes (SoA) so neighbor reads
  coalesce. Also revisit `max_coarsening_level` / bottom-solver to trim the 27,840
  tiny coarse Fapply launches (these are 66.8% of Fapply *calls* but only 5.3% of
  its *time* — a launch-count, not wall, lever).
- **Done-when:** measured wall/step reduction on the elastic region, correctness
  preserved. **Depends:** A1 SoL, C1. **Artifact:** `benchmark/PHASE_C2_fapply_memory.md`.

**C4 — Reuse the elastic operator across solves** *(setup-cost lever)*
- **Goal:** stop rebuilding the operator every solve. `Mechanics::TimeStepBegin`
  calls `Elastic::define()` (~2.4 GB alloc) + `prepareForSolve` per solve
  (`Mechanics.H:200-209`); only the coefficients change. Hoist construction; only
  re-`SetModel`. **Artifact:** `benchmark/PHASE_C4_operator_reuse.md`.

**C3 — Phase-field launch levers, only if A1 reveals a new hotspot** *(conditional)*
- **Goal:** revisit the previously-shelved 2D launch levers only if 3D ncu data
  changes the picture.
- **Steps:** the Phase-2 conclusion was that field-packing (lever 2.5) and
  Flame-loop fusion (2.4) are low-ROI — FillBoundary is 54% of launch *count* but
  15% of kernel *time*, and the packable path is architecturally blocked by the
  implicit-comp-0 convention (`gpu_perf_nsys_findings.md`). Re-open **only** if
  A1's 3D occupancy/launch data contradicts the 2D-derived conclusion.
- **Done-when:** explicit keep-shelved or re-open decision backed by 3D counters.
- **Depends:** A1, B1.
- **Artifact:** note appended to `PHASE_B1_workload_breakdown.md`.

### Phase D — Scaling envelope

Only after single-GPU is fully characterized and optimized.

---

**D1 — Multi-GPU revisit at large per-GPU domains** *(conditional)*
- **Goal:** find the per-GPU domain size where halo cost is amortized and 2 GPUs
  finally beat 1.
- **Steps:** multi-GPU is currently a measured loss (2 GPUs 3.5–13.3× *slower*
  than 1 at every tested size, halo-bound — `PHASE3_R3_crossover.md`). Revisit
  **only** with (a) much larger per-GPU subdomains and (b) GPU-aware MPI / NVLink
  confirmed in play; use A1's sync-fraction methodology to confirm halo
  staging is the bottleneck before/after.
- **Done-when:** either a crossover where multi-GPU wins, or a confirmed
  "single-GPU is the shape" with sync-fraction evidence.
- **Depends:** A1, B3.
- **Artifact:** `benchmark/PHASE_D1_multigpu.md`.

**D2 — H100/H200 rows**
- **Goal:** extend the crossover beyond A100-80 (the only device tested).
- **Steps:** re-run the A2/B3 combined matrix on H100/H200 (the NOVA build
  already targets `sm_90`, `build_alamo_nova_3d.sh`).
- **Done-when:** combined speedup rows on ≥1 newer device.
- **Depends:** A2, B3.
- **Artifact:** rows added to `PHASE_A2_combined_crossover.md` / `PHASE_B3_*`.

---

## 7. Task summary

| ID | Phase | Task | Depends | Gate/conditional | Artifact |
| --- | --- | --- | --- | --- | --- |
| A0 | A | Cut alpha-1.0 baseline | — | — | tag + `ALPHA1_BASELINE.md` |
| A1 | A | Unblock ncu/nsys on NOVA | — | — | `profiles_alpha1/` |
| A2 | A | Combined flame+elastic A100 baseline | A0 | — | `PHASE_A2_combined_crossover.md` |
| A3 | A | Lock metric set + GPU-CI perf regression | A1,A2 | — | workflow + CSV |
| A4 | A | Sibling elixir-race audit | — | — | `elixir_race_audit.md` |
| B1 | B | Elastic time-fraction breakdown | A1,A2 | **gates Phase C** | `PHASE_B1_workload_breakdown.md` |
| B2 | B | AMR-depth sweep | A1,A2 | — | `PHASE_B2_amr_depth.md` |
| B3 | B | Resolution sweep + fair CPU baselines | A2,B2 | — | `PHASE_B3_resolution_sweep.md` |
| C1 | C | Fapply register/occupancy | A1,B1 | if elastic ≥15% wall | `PHASE_C1_fapply_occupancy.md` |
| C2 | C | Elastic MLMG launch reduction | A1,B1 | if elastic ≥15% wall | `PHASE_C2_mlmg_launches.md` |
| C3 | C | Phase-field launch levers (re-open?) | A1,B1 | if 3D ncu contradicts 2D | note in B1 |
| D1 | D | Multi-GPU at large domains | A1,B3 | if domain large enough | `PHASE_D1_multigpu.md` |
| D2 | D | H100/H200 rows | A2,B3 | — | rows in A2/B3 |

Critical path: **A0/A1 → A2 → B1 → (C1/C2 if gated)**. A4 is independent
correctness insurance and can run anytime.

---

## 8. Exit criteria — alpha-1.0 → beta

The branch is "beta" when:

1. **A2 done:** a measured combined flame+elastic A100 speedup (with elastic
   *enabled*) and a known elastic wall-time fraction — the headline is no longer
   a guess.
2. **A1 done:** real `ncu`/`nsys` counters on NOVA exist for both solver
   families (the project's standing metric set is finally populated).
3. **B1+B2 done:** the workload is decomposed and the GPU-optimal AMR regime is
   chosen with accuracy parity verified.
4. **A3 done:** the standing metric set is captured automatically (CI perf
   regression on a GPU runner), so future changes can't silently regress.
5. Phase C is *either* executed (with measured wins) *or* explicitly skipped
   with B1 evidence that elastic isn't worth optimizing — both are valid beta
   states; what's not valid is optimizing blind.

Multi-GPU (D1) and H100/H200 (D2) are **beyond beta** — nice-to-have envelope
extension, not gating.

---

## 9. Cross-references

- Branch map / build matrix / IC-BC safety: `benchmark/GPU_BRANCH_GUIDE.md`
- Phase-field crossover (alpha baseline, elastic disabled): `benchmark/PHASE3_R3_crossover.md`
- Elastic fix root-cause: `GPU_BRANCH_GUIDE.md` D1 section + commit `c00f69086`
- Prior DoD checklist (Phases 0–5, all DONE): `benchmark/PHASE5_BRANCH_DONE.md`
- Static kernel resource dump (Fapply registers etc.): `benchmark/G0_BASELINE_OF_RECORD.md`
- Perf-regression harness: `benchmark/perf_regression_track.py` + `benchmark/PERF_TRACKING.md`
- AMR/launch-bound evidence: memory `gpu_perf_nsys_findings.md`
- NOVA run/diag harness: `benchmark/phase3_scaling_sweep.sh`,
  `benchmark/nova_flame_gpu_3d.slurm`, `benchmark/nova_flame_cpu_3d.slurm`,
  `benchmark/nova_flame_gpu_3d_diag.slurm`, `benchmark/phase3_nova_diag.sh`,
  `benchmark/g0_ncu_capture.sh`

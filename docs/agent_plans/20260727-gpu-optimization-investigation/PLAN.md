# TASK: gpu-optimization-investigation
# Folder: docs/agent_plans/20260727-gpu-optimization-investigation/

---

## Header

| Field         | Value                                                       |
|---------------|-------------------------------------------------------------|
| Risk tier     | 3 (Operator::Elastic / MLMG configuration)                   |
| Model         | opus                                                         |
| Verification  | partial-oracle (device lint + golden compare + strict A/B; convergence-count checks by judgment) |
| Est. scope    | investigation-wide; source edits confined to `src/Operator/Elastic.cpp`, `src/Operator/Operator.cpp` |
| Parallel-safe | no (exclusive GPU)                                           |

## Operating rules

Per `docs/llm/TASK_TEMPLATE.md`. Additional rule for this task: **no candidate is
adopted on a reasoned argument alone.** Each one is either measured by paired A/B
on this machine, or explicitly recorded as unmeasured with the reason.

## Context budget

Read first: `benchmark/status.sh` output, `docs/llm/PLAN.md`, this file.
Read: `src/Operator/Elastic.cpp`, `src/Operator/Operator.cpp`,
`src/Set/Matrix4_Major.H`, `src/Solver/Nonlocal/Linear.H`,
`src/Integrator/Base/Mechanics.H:141,208-212`,
`docs/agent_plans/20260726-gpu-optimization-sweep/results/RESULT.md`,
`docs/agent_plans/20260721-fapply-runtime-optimization/results/RESULT.md`,
`benchmark/MLMG_HIGH_CONTRAST_FINDINGS.md` (+ `mlmg_high_contrast_20260702/`).
Forbidden: `docs/archive/*`, unrelated task folders.

## Objective

Systematically evaluate the GPU optimization space on the elastic operator —
coalescing, shared-memory tiling, host/device transfers, kernel fusion and
fission, launch-count reduction, occupancy, divergence, precision, and the rest
of the requested list — against measured evidence rather than plausibility.
Retain what measures a win on this hardware; record the rest as measured-and-
rejected or as unmeasurable-here with the reason.

## Baseline (established before any change)

HEAD `19132c02c`, RTX A1000 (sm_86, 16 SM, 8 GB), gates green.

Static device resources, production `Elastic<Sym::Major>` (`Sym::Major == 1`,
`src/Set/Matrix4.H:12`), via `scripts/res_usage.sh` (`cuobjdump -res-usage`,
needs no root, unlike ncu which fails `ERR_NVGPUCTRPERM` here):

| kernel | 2D regs / stack | 2D occ ceiling | 3D regs / stack | 3D occ ceiling |
|---|---|---|---|---|
| `Fapply`   |  94 / 0 B   | 45.4% | 252 / 288 B  | 16.9% |
| `Diagonal` | 126 / 240 B | 33.9% | 255 / 1408 B | 16.7% |
| `Energy`   |  52 / 0 B   | 82.0% | 120 / 0 B    | 35.5% |

Dynamic, 2D-conservative case (`artifacts/nsys/` of the 20260726 task):
776 ms first-to-last-kernel window, 421.8 ms kernel busy, **354 ms (45.6%) GPU
idle** in 32,771 gaps of median 6.3 us; 91.2% of kernel time in `Fapply`
(76.7%) and fused `Fsmooth` (14.5%); 59% of `Fapply` launches run on <= 9 blocks
on a 16-SM device and cost ~25% of the window; explicit memcpy traffic 268
transfers / ~0 GB / 0.12 ms.

## Candidates closed by existing evidence (no experiment run)

| candidate | why closed |
|---|---|
| Minimizing host-device transfers; async transfers | Already nil: 268 transfers, ~0 GB, 0.12 ms total. Nothing to recover. |
| Loop unrolling / branch divergence in `Matrix4` algebra | Already done. `operator*(Matrix4,Matrix)`, `(Matrix4,Matrix3)` and `MulCol` are hand-unrolled to literal `data[]` indices (`Matrix4_Major.H:541-620`); the `uid` if-chain in `operator()` is not on the hot path. |
| Occupancy / register-pressure tuning | Measured dead end. Task 3.3 Step 4 cut `Fapply` 254->250 regs and 440->288 stack bytes for no material runtime gain. 2D `Fapply` has 94 regs and zero spill, so the 2D case cannot express a register win at all; and at <= 9 blocks on 16 SMs, per-SM occupancy is irrelevant because most SMs hold no block. |
| `__launch_bounds__` sweep | Already run as task 3.2 (`docs/agent_plans/20260713-launch-bounds-sweep/`): `Fapply` spills hard at `min_blocks >= 2` (128-reg cap vs 254 live). Knob exists, default off. |
| Fast math | Already the default build; `--cuda-fp strict` is the opt-out used for bit-identity checks. |
| FP64 tensor cores (DMMA) | Not present on sm_86 consumer silicon; the 9x9 contraction is too small and too irregular to map onto A100 DMMA. |
| Multi-stream overlap | No concurrency available on these decks: one box per level, so AMReX's per-MFIter stream rotation never advances (32,484 of 32,772 launches on one stream). |

## Experiments (each = one arm of a paired A/B, medians of 5)

- **E1 - launch-count reduction via MG coarsening depth.** Config-only
  (`<prefix>.max_coarsening_level`, `Mechanics.H:141`) plus bottom-solver
  variation (`elastic.solver.bottom_solver`, `Linear.H:280-282`). Targets the
  ~25% of window spent in <= 9-block launches and the 45.6% idle. **Hard
  constraint:** `MLMG_HIGH_CONTRAST_FINDINGS.md` forbids capping `mcl` without a
  smoother bottom solver; iteration counts must be reported alongside wall time,
  since a wall win bought with degraded convergence is not a win.
- **E2 - level-dependent `Fapply` block width.** The one occupancy-adjacent item
  not yet measured: 128 is currently global, chosen for the small-grid launches,
  while the grid=100 tier (32% of `Fapply` time) is SM-saturated.
- **E3 - face-flux kernel fission.** Each interior face flux is computed twice,
  once per adjacent node (`Elastic.cpp:272-289`). Halves `Matrix4` traffic at the
  cost of a 9-double/node round trip.
- **E4 - shared-memory tiling of `U`** and **E5 - AoS->SoA for `m_ddw_mf`**:
  scoped and costed only, unless E1-E3 leave time. Both need a custom launch or
  a cross-cutting storage change.

## Oracle

```bash
benchmark/status.sh                 # device-lint + golden-compare + a100-sanitizer
scripts/ab_timing.sh <case> ...     # paired wall-clock A/B (from the 20260726 task)
scripts/res_usage.sh <binary>       # static regs/stack/occupancy ceiling
```

Covers: bit-identity under `--cuda-fp strict`, device-lint bug classes, wall time.
Does NOT cover: A100 behaviour (FP64 rate differs by ~32x from sm_86, so any
FP64-throughput-bound result here may not transfer); the `use_psi=1` 3D branch,
which is not locally measurable (26 GB managed high-water on an 8 GB card).

## Closeout

- [ ] `status.sh` green; strict-build bit-identity or an adjudicated diff
- [ ] `results/RESULT.md` with retained / rejected / unmeasured for every candidate
- [ ] append-only `docs/llm/changelog/` entry; `results/DONE`
- [ ] `docs/llm/SESSION_LOG.tsv` line

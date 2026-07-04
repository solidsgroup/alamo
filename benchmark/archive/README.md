# GPU Port — History Archive

This folder holds the **closed-phase records** of the `chamber-gpu` GPU port:
superseded plans and point-in-time phase reports whose conclusions are settled.
They are kept for the record and for citation, not for active editing.

- **The live forward plan is `benchmark/GPU_ROADMAP_V3.md`.** Start there.
- Settled wins that should not be re-litigated are consolidated in
  `benchmark/SUCCESS_BOOK.md`.
- Operational/build/branch policy lives in `benchmark/GPU_BRANCH_GUIDE.md` (live).

> **Why these moved and others didn't.** The docs here are *phase narratives* —
> they describe work that is finished. Living references that are still consulted
> operationally (the build guide, the wins book, the metric harness, the current
> `ALPHA1_BASELINE.md`, the `LOCAL_A100_SPOOFING.md` gate policy, the
> `docs/gpu_*` technical references, the `docs/llm/` map) stayed at their original
> paths. Filenames here are unchanged, so any prose reference to e.g.
> `PHASE3_R3_crossover.md` still names the right file — now under `archive/`.

The bulky solver-output directories that the elastic debug work generated
(`benchmark/elastic_sensitivity_20260621/out_*`, ~5.6 G of run artifacts) were
**left in place**, not moved — only their markdown writeups are archived here.

---

## The three load-bearing evidence docs

If you read only three things from this archive, read these — the rest of the
project's forward plan rests on them:

1. **`PHASE3_R3_crossover.md`** — the **D3 = WIN @ single** crossover. Single A100
   beats a full 64-rank CPU node **9.6× @ 128³, 13.1× @ 256³** on 3D phase-field
   (elastic disabled). High confidence; the most-cited doc in the corpus.
2. **`PHASE_A_FINDINGS.md`** — the v2 Phase-A profiling that made elastic the
   project: combined flame+elastic **22.9× @ 256³**; `Operator::Elastic::Fapply`
   = **74.6% of GPU kernel time**, elastic solve **≈95% of wall**, flame 0.20%;
   Fapply at **255 regs/thread ⇒ ~12.5% occupancy**. *(Caveats now tracked in v3:
   unfair CPU baseline → task 2.A; runs were SLURM-preempted → task 2.B.)*
3. **`G0_BASELINE_OF_RECORD.md`** — the frozen Phase-0 correctness/timing/static-
   register baseline (Fapply 87–101 static regs → ~33% theoretical occupancy at
   sm_86; the runtime 255-reg spill discovered later in Phase A corrects this).

---

## Superseded plans

| Doc | What it was | Disposition |
|-----|-------------|-------------|
| `GPU_ROADMAP_V2.md` | The measurement-driven roadmap (Phases A–D). Got the first real A100 counters; B1 gate passed (elastic ≈95% of wall). | **Superseded by `GPU_ROADMAP_V3.md`.** Open items carried into v3 §2. |
| `PHASE_C0_input_ab.md` | First NOVA Phase-C sweep: tuned 256³ deck made ~2.5× more wall-progress than untuned — but no stress-field parity. | **Superseded;** folded into v3 task **2.B** (now budget-gated). |

The original P0–P5 master roadmap with decision trees D1–D4 lives at
`docs/llm/ROADMAP.md` (left in place; its forward arc is now `GPU_ROADMAP_V3.md`).

## Phase-by-phase records (Phases 0–5, the original arc)

| Doc | Phase | Settled conclusion |
|-----|-------|--------------------|
| `G0_BASELINE_OF_RECORD.md` | 0 | Baseline of record; golden compare passes; static register counts; `ncu` blocked locally (`ERR_NVGPUCTRPERM`). |
| `PHASE1_ELASTIC_DISPOSITION.md` | 1 | GPU-vs-CPU elastic at 512²–2048². **SUPERSEDED verdict** (D1=CPU-resident) — the 2048² divergence was later root-caused as the elixir UAF race and fixed; the perf numbers (CPU np8 2.27–3.40× at converging sizes) remain valid history. |
| `phase2_box_sweep/README.md` | 2 | Grid-strategy sweep; GPU/CPU parity best at `wide_512_bf32_mgs128` (~1.0×). Wide-shallow beats deep AMR on GPU. |
| `PHASE3_3D_READINESS.md` | 3 | 3D GPU readiness audit; standing input-authoring traps (BMP is 2D-only; `elastic.on=0` does NOT disable the solve — use `elastic.type=disable`; 3D needs six BC faces; ~255 regs/thread on sm_86). |
| `PHASE3_NOVA_TESTING_PROGRESS.md` | 3 | Pre-results NOVA workflow debug log (missing gres, async_out/MPI_THREAD_MULTIPLE). All resolved; **dead history.** |
| `PHASE3_R3_crossover.md` | 3 | **D3 = WIN @ single** (9.6×/13.1×). Multi-GPU is a confirmed regression (halo-bound). 512³ CPU baseline descoped. |
| `PHASE4_R4_dispatch.md` | 4 | **D4 = ISOLATE**: GPU device-model contained in the CUDA build only; shared CPU build untouched. |
| `PHASE5_BRANCH_DONE.md` | 5 | Branch definition-of-done 5/5; CPU suite 118 run / 92 verified / 0 failed. |

## v2 Phase-A / Phase-C records

| Doc | What it measured |
|-----|------------------|
| `PHASE_A_FINDINGS.md` | The combined flame+elastic A100 profiling — see "load-bearing" above. The single most important measurement in the project. |
| `PHASE_C0_input_ab.md` | First Phase-C input-lever sweep (see "Superseded plans"). |

## Correctness / debugging records

| Doc | What it recorded |
|-----|------------------|
| `GPU_TEST_SUITE_FIXES.md` | `tests/GPU/` 5P/4F → 9P/0F. Two **real source defects** fixed (`Integrator::Restart` node-fab OOB segfault; headerless restart `thermo.dat`), 2 stale decks, 1 stiffness-contrast divergence. Deferred: GPU traction-diagnostic race; C2 bit-exact restart. |
| `elixir_race_audit.md` | The v2 **A4** audit: swept all custom GPU `MFIter` loops for the local-`FArrayBox` UAF pattern — only `Operator.cpp:715 tmpfab` (already fixed via `elixir()`). v3 task **4.B** extends this to the broader async-lifetime pattern. |
| `elastic_sensitivity_20260621/REPORT.md` | Elastic sensitivity study: fast-math irrelevant to convergence; `psi_floor`/void only delay the stall. Conclusion (GPU-reduction nondeterminism) **superseded** by the elixir root cause. |
| `elastic_sensitivity_20260621/GPU_ELASTIC_DEBUG_PLAN.md` | The multi-phase divergence investigation; **SOLVED 2026-06-25** — cross-stream UAF on `tmpfab` in `interpolation()`, fixed via `elixir()`, verified end-to-end (2048² converges 7 iters → 8.5e-9). |
| `elastic_sensitivity_20260621/coarsening_sweep.md` | `max_coarsening_level` sweep on uniform 2048²; symptom characterization, now explained by the elixir fix. |
| `elastic_sensitivity_20260621/fix_notes.md` | ~700-line investigation log tracing the false leads (bicgstab acceptance → noise-floor residual → ghost-divide → coarse-operator) to the final UAF root cause. |

---

## See also (left in place — not archived)

- `docs/llm/{ROADMAP,VERSIONS,CONVENTIONS}.md` + `docs/llm/changelog/`,
  `docs/llm/perf/2026-06-20-gpu-port-report.md` — the LLM navigation map and its
  changelog/perf history. `CURRENT.md` and `ROADMAP.md` now carry a banner
  pointing forward to `GPU_ROADMAP_V3.md`.
- `docs/agent_plans/<date>-*/` — per-phase orchestration working folders
  (`PLAN.md` + `tasks/` + `results/`); the implementation-detail record behind
  these phase reports.
- `docs/gpu_device_capture_conventions.md`, `docs/gpu_elastic_device_port_plan.md`,
  `docs/gpu_safe_ic_bc_matrix.md` — standing technical references (still live).

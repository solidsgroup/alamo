# chamber-gpu version ledger

Semantic versions for the `chamber-gpu` branch, one row per era where the
GPU-vs-CPU performance picture materially changed. A "version" here is a real
git tag (`gpu-vX.Y.Z`) on a `chamber-gpu` commit, not just a doc — so
`git checkout gpu-v0.2.0` reproduces the numbers in that row's report.

**`gpu-v0.1.0` is tagged** (annotated, on `bde739bd0`, 2026-06-22). `gpu-v0.2.0`
is **not yet taggable**: `benchmark/PHASE3_R3_crossover.md`'s "FIRST CROSSOVER
MEASURED" content (the 39x/70x numbers) exists only in the **uncommitted
working tree** — `git show HEAD:benchmark/PHASE3_R3_crossover.md` still returns
the old "SKELETON" version. There is no commit yet that actually reproduces
those numbers on checkout, so tagging it now would point `gpu-v0.2.0` at a
commit that doesn't contain the thing the tag claims. Commit the
`PHASE3_R3_crossover.md` update (and the related in-progress
`Flame.{H,cpp}`/`Operator.cpp` edits noted in `CURRENT.md`, if they belong to
the same change) first, then tag.

Bump convention: **major** = a roadmap-phase gate (G0-G5) closes or a D-decision
flips; **minor** = a measured perf change worth citing; **patch** = doc/non-perf
housekeeping. Every bump gets a row here, a `changelog/` entry, and (if perf
moved) a `perf/` report.

| Version | Commit | Date | Headline | Decisions | Report |
|---|---|---|---|---|---|
| **gpu-v0.1.0** (tagged) | `bde739bd0` | 2026-06-20 | Elastic de-virtualized (Option B): GPU elastic 1.65x faster than CPU at coarse res **against a single-core baseline only — this does not survive a fair `np8` comparison: CPU is 2.27x-3.40x faster at matched resolution/iteration count, and GPU elastic fails to converge above 1024^2 (see `benchmark/PHASE1_ELASTIC_DISPOSITION.md` R1 update)**; phase-field/AMR path 1.50x-2.87x slower depending on AMR depth; wide-shallow grid halves the gap (1.74x) | Option B chosen over managed-memory/hybrid; **D1 = CPU-resident elastic (final, high confidence)** | `perf/2026-06-20-gpu-port-report.md`, `benchmark/PHASE1_ELASTIC_DISPOSITION.md` |
| gpu-v0.2.0 (blocked — see note above) | not yet committed | 2026-06-21 (16-rank, exploratory) + 2026-06-22 (64-rank, confirmed), uncommitted | CPU-node-vs-GPU crossover on NOVA 3D, elastic disabled: single A100 beats the full 64-rank CPU node **9.6x @ 128^3, 13.1x @ 256^3** (06-22 confirmation; supersedes the earlier 06-21 16-rank exploratory estimate of ~39x/~70x — do not cite the 39x/70x figures as current) | **D3 = WIN @ single** (HIGH confidence on 128^3/256^3; still open: `ncu` occupancy/registers). **512^3 CPU baseline and multi-GPU scaling descoped 2026-06-24 (user directive)** — 512^3 OOMs on a hardcoded `--mem=32G` and is not being pursued (128^3/256^3 are sufficient data); multi-GPU remains a confirmed regression (2 GPUs 3.5x-13.3x slower than 1 at every tested size) but is on the back burner, not a near-term goal. | `benchmark/PHASE3_R3_crossover.md` (working-tree version) |
| gpu-v0.3.0 (next) | TBD | TBD | Pending: Phase 4 dispatch decision (D4) folds into next version (512^3/multi-GPU no longer block this — see descope note above) | D4 pending | TBD |

## Other standing decisions (not separately versioned, carried forward)

- **D1 REVERSED 2026-06-22 (user override, not a re-measurement):** the user
  explicitly rejected CPU-resident elastic — "I don't want elastic to be CPU
  resident. I'd like for the elastic model to run on the GPU." D1 is back
  OPEN, not closed; Phase 1 / Gate G1 should be treated as re-opened, not
  passed. Original verdict and its evidence stand as history, not as the
  current decision — see `benchmark/PHASE1_ELASTIC_DISPOSITION.md` for that
  history and the new investigation memory below for current status.
  ~~D1 = CPU-RESIDENT elastic~~ (superseded; the line below is historical).
- D1 (historical, 2026-06-20) = CPU-RESIDENT elastic for the GPU phase-field
  path generally — see `benchmark/PHASE1_ELASTIC_DISPOSITION.md`.
- **D2 = regime problem, not launch-code-fixable** for residual 2D launch
  overhead — see `changelog/2026-06-20-session-handoff.md` section 1.
- **D4 = ISOLATE** (proposed, pending Runnels ratification) — see
  `benchmark/PHASE4_R4_dispatch.md`.

`gpu-v0.1.0` tag created 2026-06-22 (annotated, on `bde739bd0`). `gpu-v0.2.0`
awaits a commit — see the blocker note above.

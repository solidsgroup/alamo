# Phase 0 supplemental scout

Read-only inventory of sibling ALAMO clones/worktrees under `/home/jackplum/Projects`.
No source, refs, worktrees, or working trees were changed. The comparison tree is
`chamber-gpu`, HEAD `2e6a8f8f5430a58dd6a85b07b08d472e2bd1d5cd`, already dirty.
Verdicts use only the required vocabulary. “Committed” means Git tree/blob evidence;
dirty files are identified as working-tree observations, not branch evidence.

## Git roots and linked worktrees

| path | branch | full HEAD | repo state | relevant files/notes | verdict |
|---|---|---|---|---|---|
| `/home/jackplum/Projects/alamo-pf-gpu-base` | detached | `7e972f1e8a0ef8630fd99cb8b4899b1375bc297b` | clean | Historical GPU markdown listed below | new-info |
| `/home/jackplum/Projects/alamo-elastic-opt` | `chamber-gpu-elastic-opt` | `0ba677e37f10b0ac1b4cdcb1d7bb388924a0e00b` | dirty: `benchmark/GPU_ROADMAP_V2.md`, SLURM file, `input_copy` | Committed tree used; dirty files excluded | new-info |
| `/home/jackplum/Projects/alamo-pf-gpu-opt` | `codex/gpu-pf-structural-speedups` | `66f7fedf3824d4c53952c0bb8bdb440195c9ef25` | clean | Structural-speedup AB logs | new-info |
| `/home/jackplum/Projects/alamo-chamber-verify` | `chamber` | `8ae9ebaa6aba6e47cdc6bdf05a211a615a8fb749` | dirty: `src/Solver/Nonlocal/Newton.H`, `benchmark_rt_sweep/` | No unique relevant markdown | irrelevant |
| `/home/jackplum/Projects/alamo-fapply-322b` | `launch-bounds-sweep` | `23fc0f3b95f13929f122eff4bdf22a149beb8977` | clean | PTXAS, timing, gate logs | new-info |
| `/home/jackplum/Projects/alamo-pf-gpu-prev` | detached | `9c04ad45b9c975617c73633fdb0a88781bc362f9` | clean | Duplicate historical docs only | irrelevant |
| `/tmp/alamo-elastic-parent-d964cfab8-baseline` | detached | `d964cfab8c5ca063182dbef042a7bc1c129019be` | dirty: `baseline_logs/`, timing file | No unique relevant markdown; dirty artifacts excluded | irrelevant |
| `/tmp/alamo-elastic-void-baseline-a06c6f15d-clean` | detached | `a06c6f15d7f4825238395d7e1ac6d1dc55ec5877` | clean | No unique relevant markdown | irrelevant |
| `/tmp/alamo-elastic-void-continue` | `elastic-void-wip-20260713` | `5ae7dffa35ac5ec0bfeb9e2b9ced9200885f8227` | clean | No unique relevant markdown | irrelevant |
| `/home/jackplum/Projects/alamo-manual-build` | `manual-build` | `2e6a8f8f5430a58dd6a85b07b08d472e2bd1d5cd` | clean at scout start | Manual worktree, not evidence | irrelevant |

Other roots found (`simba`, `EM525`, `69ff9f85bd499974787aa962`, `makecirclescpp`,
`chamberutils`, `regressionrate`, and vendored dependencies) are not ALAMO GPU
implementation repositories and are irrelevant.

## Committed markdown findings absent by path from chamber-gpu

The rows below are committed tree blobs. Each row explicitly identifies the source
branch, full HEAD, and clean/dirty state. Blob checks show the final nine rows are
already present under a different `benchmark/archive/` path in `chamber-gpu`.

| path | bytes / blob | branch | full HEAD | repo state | relevant notes | verdict |
|---|---:|---|---|---|---|---|
| `benchmark/GPU_ROADMAP_V2.md` | 22,406 / `54bbdefbcd9cc63ec704fd7fb4fc0ae26acea2f3` | `chamber-gpu-elastic-opt` | `0ba677e37f10b0ac1b4cdcb1d7bb388924a0e00b` | dirty (listed dirty files above) | Historical GPU roadmap; committed blob | new-info |
| `benchmark/PHASE1_ELASTIC_DISPOSITION.md` | 10,998 / `9dca478099f4c86c0aa6c1658fc9b5151bd0dce9` | `chamber-gpu-elastic-opt` | `0ba677e37f10b0ac1b4cdcb1d7bb388924a0e00b` | dirty | Elastic failure/disposition history; possible dead ends | anti-pattern candidate |
| `benchmark/PHASE5_BRANCH_DONE.md` | 10,779 / `087cef479f20ad48f23751788967dec825bb64cd` | `chamber-gpu-elastic-opt` | `0ba677e37f10b0ac1b4cdcb1d7bb388924a0e00b` | dirty | Historical GPU branch completion record | new-info |
| `benchmark/PHASE_C1_cpu_golden_compare.md` | 3,592 / `65402aecaf9d9d5ffa8129c3455bf9a6a7e7a2c9` | `chamber-gpu-elastic-opt` | `0ba677e37f10b0ac1b4cdcb1d7bb388924a0e00b` | dirty | C1 CPU golden comparison | new-info |
| `benchmark/PHASE_C1_fapply_occupancy.md` | 9,668 / `326c4b60ab176ee907493794290c328c732fb7c4` | `chamber-gpu-elastic-opt` | `0ba677e37f10b0ac1b4cdcb1d7bb388924a0e00b` | dirty | C1 fapply occupancy study | new-info |
| `benchmark/PHASE_C1_nova_ab.md` | 5,192 / `ce28a18eaa1e7d50f5aa5fbd6a4c68d5c3046a2f` | `chamber-gpu-elastic-opt` | `0ba677e37f10b0ac1b4cdcb1d7bb388924a0e00b` | dirty | C1 NOVA A/B report | new-info |
| `docs/agent_plans/20260713-launch-bounds-sweep/results/ptxas_table.md` | 5,678 / `d4a0a540f00159a8ae06b49db2a533f464eac52e` | `launch-bounds-sweep` | `23fc0f3b95f13929f122eff4bdf22a149beb8977` | clean | PTXAS register/launch-bounds table | new-info |
| `benchmark/GPU_TEST_SUITE_FIXES.md` | 10,797 / `7b58673f2ebd7b9fc67df8fb7c7aa3a3e04110f9` | `chamber-gpu-elastic-opt` | `0ba677e37f10b0ac1b4cdcb1d7bb388924a0e00b` | dirty | Exact blob at `benchmark/archive/GPU_TEST_SUITE_FIXES.md` | irrelevant |
| `benchmark/PHASE3_3D_READINESS.md` | 3,767 / `67650d725b9426086b1e2160bbc314c726ce442c` | `chamber-gpu-elastic-opt` | `0ba677e37f10b0ac1b4cdcb1d7bb388924a0e00b` | dirty | Exact blob at archive path | irrelevant |
| `benchmark/PHASE3_NOVA_TESTING_PROGRESS.md` | 6,592 / `a2fef87859c8813e09ca934d850de8bd149858d4` | `chamber-gpu-elastic-opt` | `0ba677e37f10b0ac1b4cdcb1d7bb388924a0e00b` | dirty | Exact blob at archive path | irrelevant |
| `benchmark/PHASE3_R3_crossover.md` | 14,815 / `c3ca2ad6415d876efd7c6c7d9866dc46e224c81f` | `chamber-gpu-elastic-opt` | `0ba677e37f10b0ac1b4cdcb1d7bb388924a0e00b` | dirty | Exact blob at archive path | irrelevant |
| `benchmark/elastic_sensitivity_20260621/GPU_ELASTIC_DEBUG_PLAN.md` | 61,181 / `8da72860fb6dc78a3ea4ef05a125d436cea0de70` | `chamber-gpu-elastic-opt` | `0ba677e37f10b0ac1b4cdcb1d7bb388924a0e00b` | dirty | Exact blob at archive path; historical dead ends | irrelevant |
| `benchmark/elastic_sensitivity_20260621/fix_notes.md` | 60,693 / `fdaa4e4eac583f9b96d85bd83c3527589739ab1d` | `chamber-gpu-elastic-opt` | `0ba677e37f10b0ac1b4cdcb1d7bb388924a0e00b` | dirty | Exact blob at archive path | irrelevant |
| `benchmark/elixir_race_audit.md` | 5,047 / `2c875d088bbdf9eab61baa96f57b4426ba36f1a8` | `chamber-gpu-elastic-opt` | `0ba677e37f10b0ac1b4cdcb1d7bb388924a0e00b` | dirty | Exact blob at archive path | irrelevant |

The first seven rows are supplemental candidates. They are evidence/context only,
not Transform sources. Any contradiction is resolved in favor of `chamber-gpu`.

## Non-markdown evidence (working-tree artifacts)

These files are generated artifacts, not Git tree evidence. They are explicitly
tagged with the repository HEAD and state; they should be copied into Tier 2 only
after provenance is accepted. No same-path copies were found in the current tree.

| path | bytes | branch | full HEAD | repo state | relevant notes | verdict |
|---|---:|---|---|---|---|---|
| `docs/agent_plans/20260709-fapply-kernel-surgery/results/gate_golden_cpu.log` (under `/home/jackplum/Projects/alamo-fapply-322b`) | 262,072 | `launch-bounds-sweep` | `23fc0f3b95f13929f122eff4bdf22a149beb8977` | clean | Golden CPU gate log | new-info |
| `docs/agent_plans/20260709-fsmooth-launch-fusion/results/memcheck.log` (same repo) | 55,047 | `launch-bounds-sweep` | `23fc0f3b95f13929f122eff4bdf22a149beb8977` | clean | Memcheck result | new-info |
| `docs/agent_plans/20260709-fsmooth-launch-fusion/results/timing_summary.txt` (same repo) | 2,256 | `launch-bounds-sweep` | `23fc0f3b95f13929f122eff4bdf22a149beb8977` | clean | Timing summary | new-info |
| `docs/agent_plans/20260709-fapply-kernel-surgery/results/gate_lint.log` (same repo) | 2,341 | `launch-bounds-sweep` | `23fc0f3b95f13929f122eff4bdf22a149beb8977` | clean | Device-pattern lint gate | new-info |
| `docs/agent_plans/20260709-fapply-kernel-surgery/results/regs_baseline.txt` and `regs_after.txt` (same repo) | 960 each | `launch-bounds-sweep` | `23fc0f3b95f13929f122eff4bdf22a149beb8977` | clean | Register-count before/after pair | new-info |
| `benchmark/validate/references/gpu_alpha1_local/centre_bore_3d_128_a2/metrics.json` (same repo) | 5,670 | `launch-bounds-sweep` | `23fc0f3b95f13929f122eff4bdf22a149beb8977` | clean | GPU validation metrics | new-info |
| `benchmark/baseline_references/rod_and_tube_step2/gpu_strict.json` (same repo) | 2,098 | `launch-bounds-sweep` | `23fc0f3b95f13929f122eff4bdf22a149beb8977` | clean | Strict GPU baseline | new-info |
| `benchmark/local_ab_20260628_134838/base/run.log` (under `/home/jackplum/Projects/alamo-pf-gpu-opt`) | 30,980 | `codex/gpu-pf-structural-speedups` | `66f7fedf3824d4c53952c0bb8bdb440195c9ef25` | clean | Structural-speedup AB baseline log | new-info |
| `benchmark/local_ab_20260628_134838/opt/run.log` (same repo) | 29,478 | `codex/gpu-pf-structural-speedups` | `66f7fedf3824d4c53952c0bb8bdb440195c9ef25` | clean | Structural-speedup AB optimized log | new-info |

The pf-gpu-opt AB directory has 16 run logs; the two paths above are representative.
The current chamber-gpu worktree remains primary where newer gate logs overlap.

## Recommendation

Add the seven unique committed markdown rows and the named generated artifacts to a
separate supplemental inventory. Do not duplicate the nine archive-identical blobs,
and do not elevate any supplemental evidence to Tier 1 without direct
`chamber-gpu` diff support.

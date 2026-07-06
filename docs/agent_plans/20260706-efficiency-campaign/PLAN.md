# TASK: efficiency-campaign
# Folder: docs/agent_plans/20260706-efficiency-campaign/

---

## Header

| Field        | Value                                                        |
|--------------|--------------------------------------------------------------|
| Risk tier    | 1 (benchmark harness + docs; no src/)                        |
| Model        | fable (captain) + haiku (hygiene worker)                     |
| Verification | full-oracle (status.sh all-green)                            |
| Est. scope   | benchmark/{ci_golden_compare.sh,baseline_suite.py,status.sh}, .gitignore, docs/llm/PLAN.md, 4 stale doc headers, 6 task-folder closeouts, 3 root strays deleted |
| Parallel-safe| yes — hygiene worker file set disjoint from gate-fix file set |

## Context budget

Read first: benchmark/status.sh output, this PLAN.md
Read: benchmark/ci_golden_compare.sh, benchmark/baseline_suite.py:70-110,
      benchmark/status.sh, docs/llm/PLAN.md
Forbidden: docs/archive/*, unrelated task folders

## Objective

2026-07-06 efficiency review found: (a) golden-compare gate red because binary
selection is lexicographic (`sort | tail -1` picks stale Jun-29
`alamo-2d-perf-clang++` over the freshly built `alamo-2d-g++`; all three
"missing" solver params in fact exist at HEAD); (b) repo hygiene debt (root
strays, 954MB git garbage, stale-authority docs, 0/6 task folders closed);
(c) PLAN.md missing two newly found GPU wins. Fix all; leave gate green.

## Oracle

Command(s): benchmark/status.sh  (all three gates PASS; OPEN list accurate)
Covers: gate correctness end-to-end incl. rod_and_tube_step2 golden case
Does NOT cover: A100 A/B (needs NOVA login node); hygiene deletions (spot-check)

## Steps

### Step 1 - gate binary-selection fix (captain)
DO: ci_golden_compare.sh CPU/GPU legs select the deterministic just-built
    binary name (`bin/alamo-${DIM}d-${COMP}`, `bin/alamo_gpu-${DIM}d-nofast-
    cuda${ARCH}-${COMP}`); baseline_suite.py find_binary picks newest-mtime;
    status.sh run_gate tees gate output to benchmark/_gate_logs/<gate>.log and
    prints the path on FAIL; .gitignore covers _gate_logs.
CHECK: benchmark/ci_golden_compare.sh exits 0 (all 4 golden cases).

### Step 2 - hygiene sweep (haiku worker)
DO: delete root strays (TASK_TEMPLATE.md, claude1.md, chamber_gpu_changes.diff
    — verify duplicates/stale first); git gc to clear tmp_pack garbage;
    SUPERSEDED headers on remediation.md, docs/gpu_elastic_device_port_plan.md,
    benchmark/archive/README.md stale pointer, root tasks/ README note;
    retroactive results/RESULT.md + DONE for the 5 pre-template task folders;
    DONE for 20260705-meta-workflow-install (its only blocker was this gate).
CHECK: git count-objects -vH garbage=0; status.sh OPEN list shows only this folder.

### Step 3 - A100 A/B for C1 (PLAN.md task 3.1)
DO: per benchmark/PHASE_C1_nova_ab.md, from a NOVA login node. If NOVA
    unreachable from this machine, prepare handoff and record blocker.
CHECK: before/after wall-time + occupancy table in results/.

### Step 4 - PLAN.md backlog update (captain)
DO: extend task 3.2b to include Diagonal per-component DDW hoist
    (Elastic.cpp:860,868); add task 3.I Newton norm0 fusion (Newton.H:419,572,
    816 -> single ReduceOps pass per field, pattern Flame.cpp:616). Keep <100
    lines.
CHECK: docs/llm/PLAN.md diff; doc-budget guard passes on commit.

## Closeout

- [ ] Oracle passes; status.sh all green
- [ ] results/RESULT.md: what changed, evidence, deviations from plan
- [ ] changelog/ entry (append-only)
- [ ] touch results/DONE
- [ ] Session log line appended to docs/llm/SESSION_LOG.tsv

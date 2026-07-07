# TASK: a100-ab-c1
# Folder: docs/agent_plans/20260707-a100-ab-c1/

---

## Header

| Field        | Value                                                        |
|--------------|--------------------------------------------------------------|
| Risk tier    | 0 (docs + benchmark results; no src/ touched)                |
| Model        | sonnet                                                        |
| Verification | full-oracle (fcompare parity + ncu/TinyProfiler measurement)  |
| Est. scope   | benchmark/PHASE_C1_fapply_occupancy.md (elastic-opt branch), docs/llm/PLAN.md |
| Parallel-safe| yes — no file overlap with other open work                   |

## Context budget

Read first: benchmark/status.sh output, this PLAN.md
Read: docs/llm/PLAN.md, docs/agent_plans/20260706-efficiency-campaign/results/A100_AB_HANDOFF.md
Reference only if step names it: benchmark/PHASE_C1_nova_ab.md (chamber-gpu-elastic-opt branch)
Forbidden: docs/archive/*, unrelated task folders

## Objective

Close PLAN.md task 3.1: the A100 A/B for the already-committed C1 `Fapply`
register-pressure edits (`0bb893acc` on `chamber-gpu-elastic-opt`) was blocked
last session on NOVA ssh access from this machine. User set up SSH
ControlMaster multiplexing (`~/.ssh/config` host `nova`) and authenticated
interactively; this session drives the measurement to completion and records
it.

## Oracle

Command(s): fcompare.gnu.ex (parity), ncu --nvtx-include capture (registers/
  occupancy), TinyProfiler wall A/B (job 11433933 on NOVA)
Covers: registers/thread + achieved occupancy delta, Fapply/MLMG wall delta,
  GPU stress-field parity to tolerance
Does NOT cover: whether occupancy staying flat justifies the launch-bounds
  sweep effort/benefit tradeoff (judgment call, deferred to task 3.2)

## Steps

### Step 1 - confirm prior NOVA jobs and pull existing results
VERIFY: ssh nova reachable via multiplexed connection; jobs from 2026-07-06
  session (build 11433898, wall A/B 11433933, first ncu 11433932) completed
DO: query squeue/sacct for job state; identify retargeted ncu job 11448954
  (gates on the Fapply NVTX range directly, fixing the mis-scoped first ncu
  attempt) already running
CHECK: sacct shows COMPLETED for all jobs

### Step 2 - GPU stress-parity check
DO: run fcompare.gnu.ex on out_c1_baseline/plot vs out_c1_modified/plot at
  steps 00000/00300/00600 (cell + node)
CHECK: cell fields bit-identical; node fields within fast-math tolerance
  (~1e-7-1e-8 relative error)

### Step 3 - collect occupancy/register ncu result
VERIFY: job 11448954 COMPLETED
DO: read c1ab_ncu2.11448954.out, extract registers/thread + achieved
  occupancy for the Operator::Elastic::Fapply() NVTX-gated launches, both
  arms
CHECK: numbers recorded, delta computed

### Step 4 - record result and close out
DO: update benchmark/PHASE_C1_fapply_occupancy.md on chamber-gpu-elastic-opt
  with full before/after table; update docs/llm/PLAN.md to drop task 3.1 and
  promote 3.2/3.2b; update memory gpu_c1_fapply_source_opt.md
CHECK: PLAN.md diff, memory file diff

## Checkpoints (tier >= 2 only)

N/A — tier 0.

## Closeout

- [x] Oracle passes; status.sh all green
- [x] results/RESULT.md: what changed, evidence, deviations from plan
- [ ] changelog/ entry (append-only) — not release-worthy, measurement only
- [x] touch results/DONE
- [x] Session log line appended to docs/llm/SESSION_LOG.tsv

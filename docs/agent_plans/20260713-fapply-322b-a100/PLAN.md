# TASK: fapply-322b-a100
# Folder: docs/agent_plans/20260713-fapply-322b-a100/

---

## Header

| Field        | Value                                                        |
|--------------|--------------------------------------------------------------|
| Risk tier    | 1 (no src/ edits; NOVA runs + evidence collection)           |
| Model        | sonnet                                                       |
| Verification | partial-oracle (fcompare parity is scripted; perf verdict is judgment) |
| Est. scope   | 0 src files; NOVA job scripts + results docs                 |
| Parallel-safe| yes (disjoint from local golden-compare rebuild)             |

## Operating rules

1. Read ONLY the files listed in Context budget. docs/archive/ is forbidden.
2. Every step's VERIFY must pass before acting on that step. On failure: STOP,
   report discrepancy, wait.
3. One commit per step unless stated. Message: `<area>: <what> (<task-folder>)`.
4. No scope expansion. New ideas go to a NOTES.md in this folder, not into code.
5. If required knowledge is missing, ask. Do not wing it.

## Context budget

Read first: benchmark/status.sh output, this PLAN.md
Read: docs/agent_plans/20260707-a100-ab-c1/results/RESULT.md (procedure precedent),
      benchmark/NOVA_SLURM_RUNBOOK.md, docs/agent_plans/20260709-fapply-kernel-surgery/results/RESULT.md
Reference only if step names it: benchmark/select_nova_resources.sh
Forbidden: docs/archive/*, unrelated task folders

## Objective

Task 3.2b (PLAN.md item 2) local leg is done: commits 9470889b1 + dc02baf17 on
branch `fapply-322b` (kernel surgery: DDW hoist, column-restricted contraction,
unrolled Matrix4xMatrix3; plus Fsmooth 4D launch fusion). A1000 evidence is
indicative only (shared, 50W-capped). The remaining leg is the A100 judgment:
wall + ncu executed-instructions A/B on NOVA, plus stress parity. PASS unlocks
the merge to chamber-gpu (separate task).

## Oracle

Command(s): fcompare on plotfiles between arms (rel err <= 1e-6 on stress/strain/disp
  node fields, cell fields bit-identical, per a100-ab-c1 precedent);
  TinyProfiler Fapply exclusive wall MODIFIED <= BASELINE (no regression).
Covers: parity + perf non-regression on A100.
Does NOT cover: occupancy interpretation; verdict on whether the win justifies
  merge is orchestrator judgment.

## Steps

### Step 1 - Sync arms to NOVA
VERIFY: `ssh nova 'sinfo -p nova -h | head -1'` responds; local branch tips:
  BASELINE=d964cfab8, MODIFIED=dc02baf17 (`git -C ~/Projects/alamo-fapply-322b log --oneline -1`).
DO: git bundle both refs, scp to /work/brunnels/jackplum/alamo-322b-ab/, clone
  two arm dirs (baseline/, modified/). Reuse alamo-c1ab conventions.
CHECK: `git -C <arm> rev-parse HEAD` matches expected sha on both arms.

### Step 2 - Build both arms (A100, sm_80)
DO: sbatch build jobs per NOVA_SLURM_RUNBOOK.md (modules cuda/gcc/openmpi,
  GPU_TYPE=a100, arch 80). 3D CUDA build.
CHECK: bin/alamo_gpu-3d-* exists in both arms.

### Step 3 - Wall A/B
DO: run deck input_3d_centre_bore_256_a2 on 1x A100, identical steps both arms,
  TinyProfiler on. Same node class; plotfile writes minimized (Lustre noise
  precedent from a100-ab-c1).
CHECK: both jobs COMPLETED; TinyProfiler tables extracted for
  Operator::Elastic::Fapply(), MLMG::solve(), Fsmooth region.

### Step 4 - ncu capture
DO: ncu gated on NVTX range Operator::Elastic::Fapply() (precedent: job 11448954
  scoping fix — do NOT capture init kernels), 6 launches/arm across two grid
  sizes; collect registers/thread, achieved occupancy, executed instructions.
CHECK: ncu report contains Fapply kernels (name contains launch_global; match by
  NVTX range not kernel name).

### Step 5 - Parity + verdict evidence
DO: fcompare stress/strain/disp between arms at matched steps; assemble
  results/RESULT.md tables (wall, ncu, parity).
CHECK: oracle above.

## Checkpoints (tier >= 2 only)

n/a (tier 1) — but report to orchestrator before any src/ change (none expected).

## Closeout

- [ ] Oracle passes
- [ ] results/RESULT.md: wall table, ncu table, parity, verdict recommendation
- [ ] touch results/DONE
- [ ] Session log line appended by orchestrator

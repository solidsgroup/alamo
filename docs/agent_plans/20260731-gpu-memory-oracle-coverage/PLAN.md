# TASK: gpu-memory-oracle-coverage
# Folder: docs/agent_plans/20260731-gpu-memory-oracle-coverage/

---

## Header

| Field        | Value |
|--------------|-------|
| Risk tier    | 1 — benchmark/test harness and test input only; no `src/` change |
| Model        | opus |
| Verification | full-oracle — focused C2 checks plus `FULL=1 benchmark/status.sh` |
| Est. scope   | 6 files, baseline references, <180 lines |
| Parallel-safe| no — builds/runs shared GPU binaries and baseline output trees |

## Operating rules

1. Read ONLY the files listed in Context budget. `docs/archive/` is forbidden.
2. Every step's VERIFY must pass before acting. On failure: stop and report.
3. One commit per step. Message: `<area>: <what> (20260731-gpu-memory-oracle-coverage)`.
4. No scope expansion. New ideas go to this folder's `NOTES.md`.
5. Missing knowledge means stop and ask; do not weaken coverage to make it pass.
6. Tier 1: unattended, gates plus spot-check diff.

## Context budget

Read first: `benchmark/status.sh` output and this `PLAN.md`
Read: `benchmark/status.sh`, `benchmark/ci_golden_compare.sh`,
`benchmark/baseline_suite.py`, `tests/GPU/testlib_gpu.py`,
`tests/GPU/C2_restart_roundtrip/{input,test.py}`,
`tests/GPU/C4_amr_correctness/{input,test.py}`
Reference only if a step names it: `benchmark/local_a100_gate.sh`,
`benchmark/baseline_references/`
Forbidden: `docs/archive/*`, unrelated task folders, the Desktop sweep campaign

## Objective

The memory campaign's FULL oracle does not cover the four properties required
by PLAN v2 §5.2: a pressure history long enough to diverge, a runtime regrid,
two ranks, and a checkpoint/restart cycle. Extend the existing small C2 restart
case so one 20-step campaign oracle covers all four, compare the two-rank result
against a checked-in one-rank golden, and make the FULL gate invoke it
explicitly.

## Oracle

Commands:

- `PYTHONDONTWRITEBYTECODE=1 python3 benchmark/baseline_suite.py unit`
- `python3 tests/GPU/run_gpu_tests.py --dry-run --test C2_restart_roundtrip`
- focused one-rank record and two-rank `gpu_strict` check from Step 2
- direct C2 restart test with `ALAMO_GPU_STRICT_BIN`
- `FULL=1 bash benchmark/status.sh`

Covers: 20-step variable-pressure evolution, a level-1 regrid, two MPI ranks
versus the one-rank golden trajectory, and continuous-versus-restarted final
pressure/geometry.

Does NOT cover: running all four properties under compute-sanitizer. Tier 2
remains a single-rank direct sanitizer leg; the FULL oracle covers the four
properties collectively.

## Steps

### Step 1 — Make C2 a two-rank AMR restart case

VERIFY:

```bash
python3 tests/GPU/run_gpu_tests.py --dry-run --test C2_restart_roundtrip
```

DO: give C2 one refined level, multiple boxes, and a step-10 regrid using C4's
proven settings. Add an optional MPI-rank argument to `run_alamo`; use two ranks
for both C2 continuous and restart legs; prefer `ALAMO_GPU_STRICT_BIN` when set.

CHECK: dry-run clean; direct C2 run passes; step-10 cell/node outputs contain
`Level_1`; logs state two MPI processes.

### Step 2 — Add the 20-step golden trajectory

VERIFY: Step 1 green.

DO: add `campaign_pressure_regrid_step20` to `baseline_suite.py` using C2's
input. Require at least 20 chamber-pressure samples and a nonzero pressure span
before record/check can pass. Record CPU, GPU-fast, and GPU-strict references at
one rank; check GPU-strict at two ranks against that reference.

CHECK:

```bash
NP=1 python3 benchmark/baseline_suite.py record \
  --case campaign_pressure_regrid_step20 \
  --profiles=cpu,gpu_fast,gpu_strict
GPU_STRICT_NP=2 python3 benchmark/baseline_suite.py check \
  --case campaign_pressure_regrid_step20 --profiles=gpu_strict
```

### Step 3 — Wire the FULL gate

VERIFY: Step 2 green.

DO: make GPU mode in `ci_golden_compare.sh` run the campaign case with two
ranks and invoke C2 directly so a skipped aggregate GPU suite cannot pass it.
Update `status.sh` coverage text to name the now-enforced properties.

CHECK:

```bash
bash -n benchmark/status.sh benchmark/ci_golden_compare.sh
FULL=1 bash benchmark/status.sh
```

## Closeout

- [ ] Focused C2 and one-rank/two-rank golden checks pass
- [ ] `FULL=1 benchmark/status.sh` is green
- [ ] `results/RESULT.md` records coverage, commands, and limitations
- [ ] `touch results/DONE`
- [ ] Session log line appended to `docs/llm/SESSION_LOG.tsv`

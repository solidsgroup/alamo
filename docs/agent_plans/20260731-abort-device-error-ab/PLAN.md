# TASK: abort-device-error-ab
# Folder: docs/agent_plans/20260731-abort-device-error-ab/

---

## Header

| Field        | Value |
|--------------|-------|
| Risk tier    | 2 — scratch source variant disables a correctness tripwire |
| Model        | opus |
| Verification | partial-oracle — output equivalence plus repeated NOVA timing/profile |
| Est. scope   | 1 scratch-only header edit, 2 capture scripts, result record |
| Parallel-safe| no — reuses NOVA build and serialized timing resources |

## Operating rules

1. Read ONLY the files listed in Context budget. `docs/archive/` is forbidden.
2. Every step's VERIFY must pass before acting. On failure: stop and report.
3. The scratch source edit is never applied or committed in the working tree.
4. No scope expansion. New ideas go to this folder's `NOTES.md`.
5. Do not reinterpret a performance result as permission to remove the
   production tripwire.
6. Tier 2: stop after the plan checkpoint before creating the scratch variant.

## Context budget

Read first: `benchmark/status.sh` output and this `PLAN.md`
Read: `src/Util/Util.H`, `src/Integrator/Flame.cpp`,
`benchmark/phase0_capture.{sh,slurm}`, final Phase-0 capture artifacts
Reference only if Step 2 names it: `benchmark/build_alamo_nova.sh`,
`benchmark/phase0_analyze.py`, `benchmark/nsys_idle.py`
Forbidden: `docs/archive/*`, unrelated task folders

## Objective

Measure PLAN v2 §5.3's two unconditional `AbortIfDeviceError` stream
synchronizations before scoping Phase 2/3.  Build an identified scratch variant
that bypasses only the abort helper's synchronization and host landing, compare
it with the exact-source baseline using the same repeated protocol, and discard
the binary after recording evidence.

## Oracle

Commands:

- five randomized/interleaved repetitions of the exact and scratch variants on
  the Phase-0 production decks;
- nsys CUDA/NVTX capture of both variants;
- pressure/thermo output comparison over the capture horizon.

Covers: wall/step, CUDA synchronization count/time, GPU idle fraction, and
unchanged successful outputs for these inputs.

Does NOT cover: whether permanently removing the tripwire is safe, or paths
whose invalid values the tested decks do not trigger.

## Steps

### Step 1 — Freeze the comparison

VERIFY: final exact-source Phase-0 captures complete with zero required-leg
failures and the source key is `c88836ce414b44cc`.

DO: record compiler/toolkit, binary hashes, deck hashes, horizons, realized
timing order, and baseline pressure/thermo outputs.

CHECK: the baseline row is reproducible from its manifest and has five
successful repetitions.

### Step 2 — Build the scratch no-landing variant

VERIFY: Step 1 is green and the plan checkpoint is confirmed.

DO: in a separate NOVA scratch tree, make `Util::AbortIfDeviceError` return
without `streamSynchronizeAll()` or `flag.value()`.  Retain the flag allocation
and device writes so the comparison isolates the helper's synchronization and
host landing.  Save the one-file diff and its hash; do not alter the working
tree.

CHECK: plain and profile binaries build; a short run completes; binary and
scratch-diff hashes are recorded.

### Step 3 — Measure and report

VERIFY: scratch smoke is green.

DO: run the same five-repetition timing protocol and nsys horizon, serialized
on the same GPU type.  Compare medians under the two-standard-deviation rule,
CUDA sync count/time, idle fraction, and pressure/thermo outputs.

CHECK: `results/RESULT.md` reports a measured effect or an explicit
inside-uncertainty result.  Delete no baseline artifacts and make no production
source change.

## Checkpoints

- [ ] Plan checkpoint: human confirms the scratch-only isolation before the
      variant is created.
- [ ] Result checkpoint: human reviews hashes, output comparison, timing, and
      timeline evidence before any follow-on removal task is proposed.

## Closeout

- [ ] Five successful repetitions per arm
- [ ] CUDA synchronization and idle evidence captured
- [ ] Pressure/thermo outputs unchanged for the tested horizons
- [ ] `results/RESULT.md` records evidence and limitations
- [ ] `touch results/DONE`
- [ ] Session log line appended to `docs/llm/SESSION_LOG.tsv`

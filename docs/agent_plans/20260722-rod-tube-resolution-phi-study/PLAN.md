# TASK: rod-tube-resolution-phi-study
# Folder: docs/agent_plans/20260722-rod-tube-resolution-phi-study/

---

## Header

| Field        | Value                                                        |
|--------------|--------------------------------------------------------------|
| Risk tier    | 1 (input decks, run harness, post-processing; no src changes)|
| Model        | Codex root session                                           |
| Verification | partial-oracle                                               |
| Est. scope   | 3 input decks, 1 run script, task results and figures         |
| Parallel-safe| yes: two independent MPI simulations and disjoint outputs    |

## Operating rules

1. Read only the files listed in the context budget.
2. Preserve all pre-existing dirty-worktree changes.
3. Do not commit from this task; the repository already contains unrelated user work.
4. Stop a run only for a nonzero exit, non-finite field, or resource failure.
5. Keep all new outputs under uniquely named `output_rt1s_*` directories.

## Context budget

Read first: `benchmark/status.sh` output, this `PLAN.md`
Read: `input_rt1s_ideal`, `run_rt1s_sweep.sh`, `src/IC/Expression.H`,
      `src/IC/BMP.H`, `src/Integrator/Flame.cpp:115-142`,
      `src/Integrator/Integrator.cpp:885-1118`,
      `/home/jackplum/Projects/chamberutils/stress_validate/compare_alamo.py`,
      `/home/jackplum/Projects/chamberutils/stress_validate/quarter_rod_and_tube.json`
Reference outputs: `output_rt1s_ideal_ncell64_casingAl_void0.5_0.5/`
Forbidden: `docs/archive/*`, unrelated task folders

## Objective

Measure whether the rod-and-tube `phi`-boundary stress artifact changes with
base-grid resolution and AMR hierarchy. Run `n_cell=96, max_level=2` and
`n_cell=128, max_level=1` to one second concurrently at four MPI ranks each,
generate the existing analytical validation figures, and inspect nodal versus
face-consistent stress. Then run a baseline-resolution isolation arm whose
outer `phi` transition has width approximately `pf.eps`.

## Oracle

Commands: both resolution runs and the phi-width isolation run exit zero;
each `thermo.dat` reaches `t >= 0.999`; `compare_alamo.py` exits zero and emits
PNG/CSV results; diagnostic tables contain no non-finite stress values.
Covers: completion, output availability, resolution trend, radial/azimuthal
stress behavior, and face-versus-node reconstruction error.
Does NOT cover: full 6.5-second campaign behavior or production acceptance of
a changed material-interface model.

## Steps

### Step 1 - Derive resolution decks and parallel driver
VERIFY: baseline deck points to the completed aluminum-casing, 0.5/0.5 MPa
void recipe and stops at one second.
DO: create two decks differing only in output name, `amr.n_cell`, and
`amr.max_level`; create a two-run `np=4` concurrent driver.
CHECK: normalized diffs show only those intended deltas.

### Step 2 - Run the two resolution arms
VERIFY: the production CPU binary exists and no target output directory is
being overwritten.
DO: launch both simulations concurrently and monitor their logs.
CHECK: both exit zero, metadata records the intended hierarchy, and
`thermo.dat` reaches one second.

### Step 3 - Generate and inspect validation outputs
VERIFY: each run contains a one-second node plotfile.
DO: run `compare_alamo.py`, build a corrected spatial AMR composite, and
reconstruct the conservative face stress from displacement/model fields.
CHECK: report the `phi`-interface spike amplitude, face/node mismatch, and
comparison across the 64/96/128 base-grid cases.

### Step 4 - Epsilon-width phi isolation arm
VERIFY: determine the current `phi=0.5` radius and representable width relative
to bitmap pixels and mesh spacing.
DO: create a deterministic circular mask or analytic IC with exact plateaus
and transition width approximately `pf.eps`; run it at the 64/max-level-2
baseline hierarchy to one second at `np=4`.
CHECK: generate the same validation and stress diagnostics, explicitly noting
whether the transition is under-resolved and whether the aberration narrows or
changes amplitude.

## Closeout

- [x] `results/RESULT.md` records inputs, wall times, validation metrics, and conclusion
- [x] Figures and diagnostic CSVs saved under `results/`
- [x] `results/DONE` created
- [x] One task outcome line appended to `docs/llm/SESSION_LOG.tsv`

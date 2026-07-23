# TASK: ideal-cross-section-sweep

---

## Header

| Field        | Value                                                        |
|--------------|--------------------------------------------------------------|
| Risk tier    | 0 (decks + driver script, no src/)                           |
| Model        | fable (main session, background job)                         |
| Verification | partial-oracle                                               |
| Est. scope   | 5 new input decks, 1 driver script, 0 src files              |
| Parallel-safe| no (saturates local CPU np=8 sequentially)                   |

## Context budget

Read first: benchmark/status.sh output, this PLAN.md
Read: input_sweep_{anchor,star,multifin,centre_bore,rod_and_tube},
      input_rt1s_ideal, run_sweep.sh
Forbidden: docs/archive/*

## Objective

The 2026-07-18 psi-OFF sweep (5/5 geometries to t=6.5s) used soft casing
(5/2 GPa) and 1/1 MPa void at n_cell=96 (full) / 48 (quarter). Overnight
2026-07-18→19 the "ideal" recipe was validated on rod_and_tube quarter
(`input_rt1s_ideal`, output_rt1s_ideal_ncell64_casingAl_void0.5_0.5):
aluminum casing 70/26 GPa, void 0.5/0.5 MPa, n_cell 64 quarter
(= 128 full-domain equivalent), same psi-OFF solver recipe. Produce the
full five-geometry cross-section sweep at these ideal settings for
downstream chamberutils figure generation.

## Deck deltas (input_sweep_<g> -> input_ideal_<g>)

- plot_file: output_ideal_<g>
- amr.n_cell: 96 -> 128 (full four), 48 -> 64 (rod_and_tube quarter)
- model_void.kappa/mu: 1_MPa -> 0.5_MPa
- model_casing.kappa: 5_GPa -> 70_GPa; model_casing.mu: 2_GPa -> 26_GPa
- everything else IDENTICAL to the 2026-07-18 sweep decks (stop_time 6.5s,
  dt 2.0e-4, psi OFF recipe, At 9.0e-4 full / 5.5e-4 centre_bore /
  2.25e-4 rod_and_tube quarter)

dt=2.0e-4 at the finer dx already validated by last night's rt1s ideal run.

## Oracle

Command(s): each run exits 0 and output_ideal_<g>/thermo.dat reaches
t >= 6.49; driver log ideal_sweep_driver.log shows exit=0 x5.
Covers: completion + solver stability at ideal contrast (70 GPa / 0.5 MPa
= 1.4e5x, validated only on rod_and_tube quarter so far).
Does NOT cover: physics quality of traces; figures are the human check.

## Steps

### Step 1 - generate decks
DO: sed-derive input_ideal_<g> from input_sweep_<g> per deltas above.
CHECK: diff shows only intended lines changed.

### Step 2 - driver + launch
DO: run_ideal_sweep.sh (clone of run_sweep.sh, ideal names, np=8,
rod_and_tube first), launch background.
CHECK: driver log advancing; first plotfiles appearing.

### Step 3 - monitor to completion, then RESULT.md + DONE

## Closeout

- [ ] results/RESULT.md with per-geometry exit codes, wall times, disk
- [ ] touch results/DONE
- [ ] SESSION_LOG.tsv line

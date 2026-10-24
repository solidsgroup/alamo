# Resumed on user instruction

Resumed at 2026-09-14T17:41:04.579038+00:00. The prior pause marker is archived in pause_history/.
M03 and M24 startup100 completed. M21's interrupted attempt saved only t=0;
continuation must use the finalized initial preview, preserving the interrupted run.

Preparing bounded 2-microsecond startup extensions at approximately 34 atm,
with a 1000-step cap and snapshots every 0.25 microseconds. Cheap launch agents
run the jobs; root performs all scientific postprocessing. Current production
sweep: 0 of 12 accepted rates. Startup slopes are not Figure 10 results.

Latest initial-condition PNGs remain in initial_conditions/.

---

# Continuation record — compact, flame-initialized study

The user selected this local machine, one seed-101 disk pack per formulation,
and three pressures spanning Figure 10 (approximately 6.8, 34, 203–205 atm).
Current production manifest has 12 cases. The 84-case and original 12-case
planar manifests are archived. Do not restart the old watch_pilots.py watcher.

M03 resolves the 20 µm population and homogenizes only 0.7 µm AP. Others retain
≤20 µm homogenization. Frozen solid calibration SHA256 remains
2a7eb40ae12762fdf6074b6704b52b423bc5f93a094e685b674e18b284e5a497.
Reactive four-flame gas chemistry remains active; current decks have no
external heating and use freely evolving initial temperature/species fields.
The fixed gas coefficients differ from Gross's unavailable composition-specific
pseudo-binder/fine-particle correction inputs; retain that qualification.

Current immutable executable SHA256:
f2bc201ab63641fc1e2c95cb4076c9401f16dcb35c253075484c9a9be34dce6c.
It includes an optional clean stop-file hook, tested with one four-rank step.

## Completed

- All old laser pilots cancelled, receipts and fields preserved.
- Four planar full-height zero-step previews passed; archived images are in
  initial_conditions/planar_full_height/.
- Three 34-atm planar startup checks cleanly stopped at 100 steps, about
  0.54–0.60 µs. Their rates are not statistically meaningful. Measured maximum
  velocities 25–56 m/s and heat release 243–408 MW/m² indicate a strong startup
  adjustment. See analysis/warm_100step_startup.json.
- Prepared shortened domains with unchanged widths and finest spacings:
  M03 solid 150 µm; M17 350 µm; M21 1130.45 µm; M24 687.5 µm. Gas near 34 atm is
  150–187.5 µm, low-pressure gas 350–375 µm, high-pressure gas 100–125 µm.
- prepare_compact.py generates approximate local diffusion-flame seeds from
  surface disk intersections, conserved AP/binder stream mixing, mapped
  four-reaction progress, enthalpy-based temperature, EOS-consistent density,
  and a potential mass-flux field. These are initial guesses, not steady
  solutions. The AP seed's 830 K at 60 atm comes from Gross section 4.3;
  no calibrated parameter or surface-temperature boundary condition changed.
- Initial compact preview attempt failed before physics due to the native
  parser's 16-slot stack. Those failed receipts remain. The revised v2 inputs
  use a compressed enthalpy table (error <0.0001 K) and fewer intermediate
  variables. See preview_cases.txt for the latest attempts.
- assess_steady.py checks correlated blocks, <=5% drift/temporal CI, and at
  least one largest-particle diameter of fit-window recession. Synthetic
  constant/accelerating/short-recession verification passed. Ending a segment
  never establishes stationarity. The sweep aggregator now requires one
  accepted pack and reports temporal uncertainty, not three-pack SD.

## Work to continue

1. Four current compact PNGs are complete and visually checked. See STATUS.md,
   preview_cases.txt, and initial_conditions/initialization_audit.json. M03 uses
   v2, M17/M24 v3, M21 v4. v2 M17/M21/M24 serial attempts were blocked by
   sandbox OpenMPI sockets; v3 ran successfully outside the sandbox. M21 v3
   was superseded because cropping underrepresented its largest AP population.
   prepare_m21_compact_pack.py generated a complete smaller tile in packs_compact
   and corrected all three M21 production inputs. Its 400/200/50 µm mass
   fractions agree with Table 2 within 0.03 percentage points; width and grid
   spacing increase 0.48% to fit the whole tile on a square-cell grid.
2. M03 compact startup100 completed; M21 and M24 startup checks have been
   requested through agents (inspect actual run.json before assuming active).
   Compare compact/local-flame startup against the saved planar checks. Cheap
   agents launch; root handles scientific postprocessing. Three reusable agents
   launch_q500, launch_q200, launch_q1000 use gpt-5.6-luna low.
3. Extend only promising runs. Initial production segments end at 0.5 ms;
   there is no automatic continuation queue. Preserve restart ancestry/input
   and binary hashes. All production launches use launch_case.py.
4. Check the shortened boundaries, mesh/interface-width dependence, sustained
   gas burning, and sufficient particle-scale recession before final rates.
   One pack cannot quantify packing-to-packing variability.
5. Complete Figure 10 errors/plot after accepted rates exist. Currently zero
   accepted rates; do not present startup slopes or analytic estimates as
   Low-Mach validation results.

Python: /home/esandall/Software/anaconda3/bin/python. Writable repo and /tmp.
Existing MPI restart checks and prior solver fixes remain documented in
README_original_84case_plan.md. Preserve all unrelated working-tree changes.

## Latest diagnostic details

All initial species are nonnegative and T is 300–3126 K. Input expressions
satisfy the EOS, but initial AMR interpolation leaves a maximum 0.234% stored
volume residual in M24 (much smaller in other cases). Do not describe every
initialized cell as roundoff-consistent. M03 startup100 retained almost steady
gas heat (+2% versus +148% for the planar seed), but its mean timestep is still
5.5 ns and is not improved relative to the planar 6.0 ns mean. Wall timing
loads differ; do not claim a controlled speedup. compare_startups.py writes
analysis/startup_comparison.json. The first restart plot can contain zero
qdot before reference-pressure synchronization; use the checkpoint diagnostic
for this identical initial state (history now does so with explicit provenance).

The current binary remains f2bc201a... . No production case is launched.
All 12 current input hashes and the frozen calibration hash match. M03/M24
startup100 use preview v2/v3; M21 startup100 was re-prepared from corrected v4
before launch. Agent launch_q500 completed M03; launch_q200 owns M21 and
launch_q1000 owns M24. The 100-step runs are intentionally short diagnostics.

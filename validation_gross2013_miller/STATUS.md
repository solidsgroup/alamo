# Transient Figure 10 sweep with statistical stopping

The user cancelled the experimental runs and requested removal of the
quasi-steady implementation, then authorized a transient sweep that stops
as soon as each rate is statistically steady.

The implementation, tests, documentation and dedicated helpers were removed
from the active workspace. The rebuilt executable exactly matches the
pre-feature binary:
`f2bc201ab63641fc1e2c95cb4076c9401f16dcb35c253075484c9a9be34dce6c`.
All 12 original input/metadata pairs were restored byte for byte. Previous
calibration, mixing, restart and pressure-solver work is preserved.
Cancelled logs and feature diagnostics are archived at
`/tmp/alamo-removed-quasi-20260915`. Host process inspection confirmed both old
queues and MPI simulations were stopped.

Twelve new `_transient_steady` cases are prepared: one pack each for M03, M17,
M21 and M24 at three pressures spanning approximately 6.8–205 atm. The four
middle-pressure cases reuse finalized 2-microsecond transient checkpoints.
Other cases use the existing local-flame initial conditions, with no laser.
Calibrated material properties, gas chemistry, fine grids and packings are
unchanged. The calibrated solid parameter SHA256 remains
`2a7eb40ae12762fdf6074b6704b52b423bc5f93a094e685b674e18b284e5a497`.

The root-authored statistical monitor is active in session39810, and both
transient queues have launched successfully. M03 at34.3921 atm and M21 at34.0775
atm are advancing physical time beyond their 2-microsecond checkpoints; ten
cases are queued. The monitor has read both initial histories without errors.
There are no statistically accepted transient production rates yet.

Small queue: session89714 (M03/M17), four MPI ranks per active case.
Large queue: session17839 (M21/M24), four MPI ranks per active case.
Automatic early stopping and time recording are active.

The monitor checks completed plotfiles every 15 wall-clock seconds. Acceptance
retains the existing startup exclusion, at least one largest-particle diameter
of fit-window recession, four correlation blocks, <=5% relative temporal 95%
uncertainty and drift, sustained heat release, and a later confirmation window.
It writes STOP at acceptance; the existing solver hook stops cleanly with final
output. It records simulated detection time, first passing time, wall time,
rate, uncertainty and actual stopping time in each case's
`analysis/steady_stop.json` and the aggregate `analysis/time_to_steady.csv`.

If monitoring is unavailable for five minutes, the launch worker requests an
unaccepted stop and does not start another case. Ceiling time does not imply
acceptance. M03/M17 share one serial queue; M21/M24 share the other. Reference
rates set only output cadence and a ceiling, never statistical acceptance.
Detection precision is limited by the recorded output cadence and polling.

Validation: three control tests passed (short/drifting histories rejected,
later confirmation required, automatic STOP/timing recorded, incomplete
plotfiles excluded). All 12 input hashes and physical parameters were checked.
The executable restoration is exact by SHA256. No solver changes were added
for statistical stopping; the existing transient stop-file hook is used.

## Latest run check

Checked 2026-09-14T20:42:07.215584-05:00. Two runs are active, ten are queued, and none has completed. The statistical monitor is active and its heartbeat is fresh.

- M03 at 34.3921 atm: 41.9 minutes elapsed; physical time 7.982 microseconds; timestep 7.372 ns; step 1115.
- M21 at 34.0775 atm: 41.9 minutes elapsed; physical time 12.831 microseconds; timestep 20.925 ns; step 815.

Both solvers are advancing physical time. No new scheduled plotfile has been reached, so the monitor has only the initial restart histories so far. Neither run has a statistically accepted rate. Details: analysis/transient_run_progress.json.

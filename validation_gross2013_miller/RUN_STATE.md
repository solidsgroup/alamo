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

## Active processes and continuation

Latest user instruction authorizes these transient runs and statistical early
stops. No pause is pending. Production manifest: production_cases.txt. Original
restored manifest: production_cases_before_transient_steady.txt.

Root monitor session39810: monitor_transient_sweep.py --interval15.
Its PID in transient_monitor_state.json is inside a sandbox PID namespace;
use its unified session for control. It performs the scientific postprocessing.

Cheap launch-only workers:
- quasi_baseline_launch: launch_transient_queue.py queues/fig10_transient_small.json.
- quasi_sandwich_launch: launch_transient_queue.py queues/fig10_transient_large.json.
The small queue is session89714 and the large queue is session17839. Both
have real solver receipts and are advancing. Do not restart old queue
sessions91157/96249; those were cancelled. Do not start old pilot
watchers or launch duplicate cases while approval is pending.

Queue receipts: queues/fig10_transient_{small,large}.state.json.
Monitor heartbeat/results: analysis/transient_monitor_state.json.
Per-case run.json and stdout.log record actual execution and exit status.
Root monitor refreshes production_status.csv, time_to_steady.csv, and the
accepted Figure10 comparison. Plotfiles are read only after celloutput.visit
lists their completed writes. A postprocessing error repeated three times
requests an unaccepted stop, rather than allowing an unmonitored run to its
ceiling. Failed numerical cases remain failed and the queue proceeds to others.

Saved transient checkpoints remain under runs/:
- M03_p34.3921_seed101_compact_diffusion_relax2us_v2/output/00307cell
- M17_p34.3921_seed101_compact_diffusion_relax2us/output/00222cell
- M21_p34.0775_seed101_compact_diffusion_relax2us/output/00225cell
- M24_p33.6723_seed101_compact_diffusion_relax2us/output/00235cell

M03 resolves20-micron AP and homogenizes0.7-micron AP; others homogenize<=20
microns. No calibrated parameter changed. Gas coefficients remain fixed,
so this is a predictive comparison with Gross/Miller, not a reconstruction
of unavailable composition-specific gas/fine-particle corrections. One pack
per formulation supplies temporal uncertainty, not inter-pack variability.
Mesh/interface-width and shortened-domain checks remain separate from
statistical acceptance. Preserve all other preexisting dirty source files.

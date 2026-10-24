# Gross (2013), Figure 10: calibrated-solid 2D comparison

Work in progress. The current sweep is **one disk packing (seed 101) per
formulation, three pressures each, 12 cases total**, on this local machine.
The endpoints and middle selected Miller data are approximately 6.8, 34, and
203–205 atm. No accepted Figure 10 regression rate exists yet.

The user requested shortened solid/gas domains and preexisting flames instead
of laser ignition. Current decks contain no external heating. All four gas
reactions remain active. Cheap launch agents execute simulations; the root
model authors and runs scientific postprocessing.

## Transient sweep with statistical stopping

The transient executable is restored exactly to its pre-feature SHA256
`f2bc201ab63641fc1e2c95cb4076c9401f16dcb35c253075484c9a9be34dce6c`.
The removed experimental implementation, diagnostics and cancelled run logs
are archived at `/tmp/alamo-removed-quasi-20260915`. The original twelve input
and metadata pairs are restored byte for byte. Current production cases have
new `_transient_steady` directories and use the same calibrated physics and
packs. The four middle-pressure cases restart their completed 2-microsecond
transient checkpoints; other pressures use the existing flame initial fields.

`monitor_transient_sweep.py` is root-authored postprocessing running separately
from the cheap launch agents. It reads only plotfiles whose completed writes
appear in `celloutput.visit`, assesses every new snapshot, and requests a clean
stop as soon as `assess_steady.py` confirms statistical stationarity. The existing
criteria below are retained. Detection is limited by plot cadence and the
15-second monitor interval; the solver polls its STOP file every 100 steps.

Each accepted case records physical detection time, first passing time, elapsed
wall time, rate and temporal uncertainty in `analysis/steady_stop.json`. The
final exit and actual stop time are added after the solver exits. The aggregate
table is `analysis/time_to_steady.csv`. Saved-startup time is recorded separately
for restart cases. No startup slope or unfinished run is reported as an accepted
Figure 10 point.

Two four-rank serial queues (`queues/fig10_transient_small.json` for M03/M17 and
`queues/fig10_transient_large.json` for M21/M24) launch the cases automatically.
They pin binary/input hashes and require a live statistical monitor. If its
heartbeat disappears for five minutes, the active case is stopped without
accepting a rate and that queue ends. The worker scripts do no scientific
postprocessing. Queue state and per-case run receipts record actual execution.

Output spacing is at least 25 microseconds and approximately 1/40 of the
reference particle-traversal time for the larger packs. This limits output
volume while sampling the particle timescale. Reference rates determine only
cadence and a safety ceiling; acceptance depends entirely on computed fields.
A ceiling reached without stationarity is not an accepted result. No worker
should be started twice; the historical pilot watcher remains disabled.

## Calibrated solids

`parameters_solid_frozen.json` is an unchanged copy of the selected Chen
calibration, SHA256
`2a7eb40ae12762fdf6074b6704b52b423bc5f93a094e685b674e18b284e5a497`.

| Property | Pure binder | Pure AP |
|---|---:|---:|
| Density (kg/m³) | 920 | 1950 |
| Specific heat (J/kg/K) | 2130 | 1297.9 |
| Conductivity (W/m/K) | 0.213 | 0.4186 |
| Q (cal/g) | −66 | −100 |
| Speed prefactor A (m/s) | 24.6083329584 | 1179.584761920 |
| E/R (K) | 5568.850642399 | 10739.153221088 |
| E (kJ/mol) | 46.30200049 | 89.29028801 |
| E (kcal/mol) | 11.06644371 | 21.34089102 |

Density is volume-additive; cp and Q are mass-weighted; ln(A) and E/R are
volume-weighted. Conductivity uses the physical 2D Chen branch. Physical A is
converted to the phase-field multiplier using A/(1.5 × interface width).

Table 2 percentages are mass fractions of the whole propellant, with 12.63%
binder. M03 now resolves the 20 µm AP population, as explicitly requested;
only its 0.7 µm AP is homogenized. Other formulations homogenize AP ≤20 µm.

| Formulation | Resolved AP (µm) | AP mass fraction within matrix | Tile resolved area fraction | Matrix density (kg/m³) |
|---|---|---:|---:|---:|
| M03 | 20 | 0.7143180 | 0.4887852 | 1477.4523 |
| M17 | 90 | 0.8154049 | 0.276678 | 1616.0226 |
| M21 | 400, 200, 50 | 0.5199544 | 0.6456100 | 1268.3402 |
| M24 | 200, 50 | 0.5199544 | 0.6456100 | 1268.3402 |

`packs_reduced/` contains the seed-101 disk tiles. M03 has 256 disks of 20 µm.
The other seed-101 packs were preserved. M03, M17, and M24 compact domains crop these same
geometries; they do not rescale particles. M21 uses an entire smaller tile in
`packs_compact/`, with 5×400 µm, 20×200 µm, and 107×50 µm disks. Its width
changes by 0.48% to fit this tile on a square-cell grid. This avoids the
underrepresentation of 400 µm AP caused by cropping the old M21 tile. Tile fractions and actual cropped,
discretized bed fractions are different quantities, both recorded. Vertical
image disks implement the periodic packing tile; solver periodicity is only
horizontal. Old `packs/` and all old runs remain available.

## Compact domains and flame initialization

| Formulation | Width (µm) | Solid depth (µm) | Gas height near 34 atm (µm) | Finest spacing (µm) |
|---|---:|---:|---:|---:|
| M03 | 400 | 150 | 150 | 0.390625 |
| M17 | 800 | 350 | 150 | 0.78125 |
| M21 | 2009.69 | 1130.45 | 188.41 | 0.9812934 |
| M24 | 2000 | 687.5 | 187.5 | 0.9765625 |

Gas height increases to about 350–377 µm at the low pressure and decreases to
100–126 µm at the high pressure. M03/M17/M24 widths and finest spacings are
unchanged; M21 changes them by 0.48% for its complete compact tile. These dimensions remain subject to boundary checks: the reaction
zone must fit below the outlet, and the regressing front must retain a cold
solid buffer above the base. A steady enlarged-domain comparison is still
needed before claiming independence from domain height.

The first planar seeds produced substantial startup heat release: after
100 steps, about 0.54–0.60 µs, maximum velocities were 25–56 m/s. These checks
were cleanly stopped and are diagnostics, not steady rates.
`analysis/warm_100step_startup.json` records the measured fields.

`prepare_compact.py` builds approximate local diffusion-flame seeds above the
exposed AP/matrix boundaries. A constant-diffusivity convective mixing model
smooths the surface streams. Reaction progress from a saved developed flame
is mapped onto the local conserved AP/binder mixture. Six gas species remain
nonnegative and normalized. Temperature uses a local heat-release estimate;
density satisfies the EOS and diffuse gas occupancy. A potential mass-flux
field supplies an initial velocity consistent with varying surface streams.
This is an initial guess, **not a solved steady flame**. Surface temperature
and gas fields evolve freely after initialization. The AP seed uses Gross's
830 K at 60 atm reference with pressure scaling; this changes no calibrated
kinetic parameter and imposes no surface-temperature boundary condition.

The source flame and all assumptions are recorded in
`initial_conditions/flame_seed_source.json` and per-case metadata. Initial
fields are checked through actual zero-step AMReX output. Analytic input
densities satisfy the EOS exactly. AMR interpolation can introduce small
initial residuals; report the measured residuals from initialization_audit.json
rather than claiming all stored cells agree to roundoff. PNGs show exact input
disks, initialized gas temperature, reaction contours, and profiles. Cold
solids are coarsely represented away from the front in AMR and reconstructed
from the exact disks as the front approaches; the disk illustration is not a
claim that the entire cold bed is already finely resolved.

## Gas model and energy accounting

`gas_parameters.json` retains the existing six-species/four-reaction Rocfire
coefficients. Homogenized material emits AP_gas and HTPB_gas in mass proportion;
resolved AP emits AP_gas. Condensed Q is applied once during mass transfer;
legacy `chemistry.model.rocfire.qsolid` is zero. The pure-material phase-change
heat input defaults to zero for other decks. Expanded unit checks passed.

This remains a fixed-coefficient approximation to Gross. The paper's complete
composition-dependent pseudo-binder inputs and fine-particle correction table
are unavailable. They are not reconstructed by fitting Figure 10. Resolving
AP creates heterogeneous gas streams but does not guarantee that every stream
converts to the species named Final: Primary is also a product pool and Mono
may remain when its reaction partner is exhausted. Species inventories and
approximate outlet fluxes are recorded to assess this explicitly.

## Rates, uncertainty, and reference data

Figure 10 was digitized from PDF page 9 (printed page 990). Audited diamond
coordinates, solid-line knots, source hashes, and the overlay are in
`reference/`. Graphical rate uncertainty is approximately 2%, not an
experimental confidence interval. The caption's M02 conflicts with the panel
and Table 2; the study uses M03. M21's Gross solid curve is not extrapolated.

`analyze_runs.py` integrates condensed volume over valid AMR cells, without
counting covered coarse cells twice. Rate is recession per unit time. It also
records surface temperatures at η=0.5 crossings, gas heat release, extrema,
gas inventories, and cell-centered outlet-flux estimates. Restart ancestry is
included only before the selected checkpoint. The first restart plot precedes
reference-pressure synchronization and can have a zero qdot diagnostic; the
history uses the selected checkpoint's initialized heat-release diagnostic
for that identical physical instant, with explicit provenance. The original fit summary is a
diagnostic; it does not by itself accept a rate.

`fit_front_regression.py` plots the deepest eta=0.5 surface position as recession
from its initial position. It fits a line with a free intercept over the full
saved history and after a prescribed startup cutoff (by default, the first
half of elapsed time). `analysis/front_regression/` contains the PNG, source
histories with restart provenance, fitted slopes, and a sensitivity scan over
every cutoff leaving at least four samples. Run the script again to incorporate
new completed outputs; `--discard-fraction` controls the displayed cutoff.
These slopes describe the deepest front, whose horizontal location may change,
whereas integrated condensed volume describes mean regression. A linear fit
gives an interval-average rate. Truncation and a high R-squared alone do not
establish a statistically steady rate.

`assess_steady.py` excludes at least 0.25 ms of startup, requires recession of
at least one largest-particle diameter in the fit window, estimates temporal
correlation, and requires at least four blocks. Relative approximate 95% mean
uncertainty, drift between halves, block trend, and gas heat drift must each
be ≤5%. A passing result must persist for a further block duration. Statistical
stationarity does not replace mesh/interface-width or domain checks. One pack
provides no inter-pack uncertainty; no three-pack standard deviation may be
reported. Synthetic checks reject accelerating and insufficient-recession
histories and accept a sufficiently long constant-rate history.

The former fixed 0.5-ms segments are replaced by monitored transient runs.
The configured stop time is an upper ceiling; the normal stopping condition
is a confirmed statistically steady rate. The existing `run_control.stop_file`
hook supplies final output and a normal solver exit. Its four-rank one-step
check finalized normally. Runs that cool/extinguish or lose reliable monitoring
stop without accepting a burning rate.

## Provenance and existing verification

`launch_case.py` checks input hashes and runs an immutable content-addressed
binary. Receipts record actual MPI ranks, times, return codes, and binary/input
hashes. Launched cases cannot be overwritten or silently rerun. The current
executable hash is
`f2bc201ab63641fc1e2c95cb4076c9401f16dcb35c253075484c9a9be34dce6c`.

Pressure-solver options and restart-copy/initialization fixes are documented
in the archived `README_original_84case_plan.md`. Restart initialization agreed
exactly across 1/4/8 ranks; evolved temperature differed by at most 1.6e−6 K
in the recorded 100-step comparison. These are numerical checks, not validation
of a steady packed burning rate.

The old 84-case and 12-case planar manifests are archived. Current lists are
`production_cases.txt` and `preview_cases.txt`. Old laser pilots are stopped;
the old `watch_pilots.py` watcher must not be restarted for this reduced study.
Use `/home/esandall/Software/anaconda3/bin/python` for the Python dependencies.

## Files for transferring the study

Git includes the source scripts, frozen parameters, reference data, pack files,
initial-condition profiles, manifests, and case inputs/metadata. Full solver
outputs, archived executables, process receipts, and live analysis are local
run data. Rebuild the executable for the destination compiler/MPI environment
and record its new hash before creating destination launch manifests.

The four middle-pressure production inputs refer to saved transient
checkpoints. Transfer their complete checkpoint directories and the associated
restart ancestry separately, or regenerate the warmup runs before launching
those inputs. A Git checkout alone does not contain these saved fields.
The local queue manifests retain their original binary/input hashes and are
provenance for this machine, rather than ready-made launch manifests for a
different build. Use new case directories and receipts for a cluster run.

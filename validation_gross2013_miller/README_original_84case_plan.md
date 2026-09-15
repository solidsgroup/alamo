# Gross (2013), Figure 10: calibrated-solid 2D comparison

Work in progress: packing/reference preparation is complete; ignition and
resolution pilots are running. No production burning-rate comparison is valid
yet. Do not treat a cooling transient or a pilot rate as a Figure 10 result.

The requested study uses the previously selected Chen calibration, reactive
four-flame gas chemistry, external heating only during ignition, AP particles
of diameter ≤20 µm homogenized with binder, and three generated 2D disk packs
per formulation. The original Chen validation directory is retained.

## Materials and geometry

`parameters_solid_frozen.json` is an exact copy of the selected
`../validation_homogeneous_chen2002/parameters_current.json`, SHA256
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

Table 2 entries are mass fractions of the entire propellant, including 12.63%
binder. They are not normalized within the AP. `reference/formulations.csv`
contains the exact conversions and homogeneous matrix properties.

| Formulation | Resolved diameters (µm) | Fine AP in matrix (mass %) | Resolved area (%) | Matrix density (kg/m³) |
|---|---|---:|---:|---:|
| M03 | none | 87.3700 | 0 | 1708.4266 |
| M17 | 90 | 81.5405 | 27.6678 | 1616.0226 |
| M21 | 400, 200, 50 | 51.9954 | 64.5610 | 1268.3402 |
| M24 | 200, 50 | 51.9954 | 64.5610 | 1268.3402 |

Density is volume-additive; cp and Q are mass-weighted; ln(A) and E/R are
volume-weighted. Conductivity uses the physical 2D Chen branch. Physical A is
converted to the phase-field multiplier using A/(1.5 × interface width).

`generate_packs.py` generates nonoverlapping periodic disks by gradual
inflation and overlap-energy minimization. Seeds 101, 202, and 303 are
independent. Each `.json` records disk counts, realized size-bin mass
fractions, box size, and the minimum interparticle gap. Integer counts produce
small, explicitly recorded errors in the rounded Table 2 size-bin fractions.
Box height enforces the total resolved area fraction exactly. The underlying
tile is periodic in x and y; `.xyzr` files include vertical image disks so the
initial planar surface cuts an unbiased packing. The solver is periodic only
in x, with a solid base and gas outflow at the top.

M03 contains no resolved disks; its three realizations are identical planar
controls and cannot measure packing variability. `packs/pack_overview.png`
shows all generated geometries.

## Gas model and energy accounting

`gas_parameters.json` records the existing Alamo six-species/four-reaction
Rocfire coefficients. The study retains the existing homogeneous product
split: fine AP produces AP_gas and binder produces HTPB_gas, in mass
proportion. All four gas reactions remain active. Gas advection, conduction,
diffusion, and the low-Mach pressure projection are enabled.

This is a fixed-coefficient approximation to the paper. Gross uses a
composition-specific pseudo-binder gas, adiabatic flame-temperature inputs,
and a correction for unresolved AP particle size. The PDF supplies a reference
gas parameter set, but not the full composition-dependent input data or
correction table. Those missing functions are not silently reconstructed by
fitting Figure 10. Consequently, reproducing the plot is a predictive
comparison of the requested calibrated-solid model, not an exact reproduction
of the original Rocfire calculation.

Condensed Q is applied once, when solid mass becomes gas. The legacy
`chemistry.model.rocfire.qsolid` is explicitly zero. Resolved AP uses the new
generic `phase_change.heat_release` input; the homogeneous matrix uses its
existing mass-weighted constituent heat inputs. The default for existing
pure-material decks remains zero. The expanded phase-change unit test checks
the pure-material mass transfer and endothermic energy balance; the full unit
test executable passed (`analysis/unit_tests.log`).

Ignition uses a spatial Gaussian heating pulse with an exact time cutoff.
External heat is absent throughout any accepted rate-fit window. Initial
temperature is 300 K and no surface temperature is prescribed. The initial
0.3 ms ignition pilot cools after heating ends; additional ignition and mesh
pilots determine whether sustained burning can be established. Stronger pilots
use a 1 ms pulse at the same 30 MW/m² intensity.

The stronger ignition pilots encountered pressure-solver convergence failures
before the pulse ended. A semicoarsening trial and a conjugate-gradient bottom
solver trial failed and are excluded. The semicoarsening code was removed;
its trial binary is archived for reproducing the failed benchmark. Optional
`projection.bottom_solver` and `projection.bottom_max_iter` controls now allow
BiCGStab with 1000 bottom iterations (defaults remain unchanged for other
decks). An identical-checkpoint, one-step comparison had exactly equal
temperature, interface, solid mass, and heat release; relative pressure and
velocity differences were about 2.5e−10. See
`analysis/bottom_solver_comparison.json`. Longer restart pilots test robustness
through ignition; this one-step check alone does not validate a burning rate.

The coarse 1 ms ignition pilot has advanced beyond the heating cutoff with an
active flame. The two fine pilots have not yet reached the cutoff. There is
still no accepted steady burning rate or resolution-converged Figure 10 point.

MPI restart checks exposed two initialization problems in the generic
integrator. Restart now uses `ParallelCopy`, since the plotfile reader's rank
distribution can differ from the destination distribution. Cell allocations
also start from zero, as fresh and nodal initialization already do, before
loading valid saved cells and filling boundaries. This avoids uninitialized
scratch storage and outer ghost corners. The first restart correction loaded
all 48 saved fields exactly on 1, 4, and 8 ranks, but an additional four-rank
attempt still failed the check that includes ghost cells; the allocation fix
subsequently passed a 100-step four-rank run and two additional restart
initialization checks. Serial/eight-rank evolution over 100 steps
differed by at most 1.6e−6 K in temperature and 1.4e−9 in solid volume fraction
(`analysis/mpi_restart_100step_comparison.json`). No chemistry, material
parameters, or prescribed ignition history were changed by these restart fixes.
The final allocation fix also agrees with the original serial evolution to
1.6e−6 K in temperature and 1.1e−9 in solid fraction after 100 steps
(`analysis/mpi_restart_initialized_100step_comparison.json`). The full M24
seed-101 packed ignition pilot is now running on eight ranks.

## Reference data and postprocessing

Figure 10 is an embedded 1153×1068 JPEG in PDF page 9, printed page 990.
`prepare_reference.py` extracts manually audited diamond centers and solid-line
knots using logarithmic axes. The original image, coordinates, hashes,
digitization overlay, and graphical resolution are saved in `reference/`.
Graphical uncertainty is approximately 2% in rate (two source pixels), not an
experimental confidence interval. The caption says M02, but its panel,
Table 2, and text identify M03. The study uses M03.

Seven experimental pressure points per formulation and three realizations
give 84 requested production runs. Reference pressures are taken from the
digitized diamonds; their stored decimals do not represent measurement
precision. The Gross solid line is a separate comparator. The M21 solid curve
ends earlier than the experimental series and must not be extrapolated.

Cheap agents launch simulations and record actual terminal exit codes plus
input/binary hashes. The root/current model handles all scientific
postprocessing. `analyze_runs.py` integrates condensed volume over AMR cells
without counting covered coarse cells twice. Mean regression speed is the
slope of that volume divided by domain width. It also records surface
temperatures interpolated at η=0.5 crossings, integrated gas heat release,
pressure/temperature extrema, and multiple-crossing columns.
Restart histories include only the parent history before the selected
checkpoint, so total recession remains referenced to the original initial
surface. Surface interpolation allocates a finest-level strip around the
front; checks against a full-domain covering grid agree exactly
(`analysis/surface_strip_verification.json`).

Acceptance requires burning after ignition, adequate particle-layer recession,
stable rates across late fitting windows, and mesh/interface-width checks.
Time samples are correlated: regression-fit residuals must not be presented as
independent statistical uncertainty. Packing variation must use the three
independent geometries. A completed executable alone does not establish an
accurate mean burning rate.

## Local execution and cost

The user explicitly selected this machine. Production inputs are prepared but
remain unlaunched while ignition and resolution are tested. The 84 requested
durations sum to about 9.83 s of simulated time. The maximum allowed timestep
alone implies at least 19.7 million steps. A preliminary M24 100-step test
advanced 49.51 µs in 113.1 s on one rank and 63.7 s on eight ranks, including
restart and output. These are early ignition timings; the developed flame
can require smaller timesteps and different AMR coverage. Extrapolation
suggests months for the complete sweep on this machine at the current grid,
not a short overnight run. Local pilot execution continues before committing
the full batch to these provisional numerical settings.

`launch_case.py` checks each prepared input digest and snapshots the executable
under its content hash before launching. Receipts capture the actual command,
start/finish times, MPI rank count, return code, and immutable executable hash.
Cases cannot be silently overwritten or relaunched. One failed exploratory
fine v3 pilot had only a post-launch binary digest; its receipt explicitly
marks the executed binary hash unverified and it is excluded from validation.
Production decks now request a final plot even if a step limit is encountered,
and their step ceiling no longer truncates long low-pressure runs at three
million steps. Completed duration is still checked against physical time.

## Commands

Use `/home/esandall/Software/anaconda3/bin/python` for Python dependencies.

```bash
python validation_gross2013_miller/prepare_reference.py
python validation_gross2013_miller/generate_packs.py
python validation_gross2013_miller/prepare_runs.py --pilot
python validation_gross2013_miller/analyze_runs.py validation_gross2013_miller/runs/CASE
python validation_gross2013_miller/launch_case.py validation_gross2013_miller/runs/CASE --ranks 8
```

These preparation and analysis commands do not launch production simulations.
Each case stores its full standalone `input` and `case.json`. Preserve failed
pilots; subsequent attempts receive distinct case suffixes.

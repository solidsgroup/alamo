# Chamber Coupling Review Fixes - 2026-06-30

Scope: Flame, chamber ballistics, elastic chamber traction, mechanics coupling,
restart thermo state, and the focused GPU correctness/smoke tests that exercise
those paths.

## Bugs fixed

- `variable_pressure=1` no longer depends on `thermal.on=1` for the chamber
  mass-flow path. Flame now registers the temperature/mass-flux scratch fields
  unconditionally, registers chamber volume/area/mass-flux integrals when either
  thermal transport or variable pressure needs them, initializes those fields
  when thermal transport is disabled, and computes `mdot` for pressure evolution
  even in non-thermal runs.
- Chamber ballistics now aborts on invalid ODE inputs or pressure updates
  instead of silently propagating non-finite pressure.
- Restart now restores non-extensive thermo variables, including
  `chamber_pressure`, from the checkpoint-time `thermo.dat` row. The parser uses
  whitespace tokenization so tab-delimited thermo output is handled correctly.
- Time-based thermo integration now checks `amr.thermo.plot_dt` consistently
  instead of mixing `thermo.dt`, the global `plot_dt`, and `thermo.plot_dt`.
- `eta_old_mf` now carries the same three ghost cells as `eta_mf`; the two fabs
  are swapped during advance, so mismatched halos could leave the elastic
  eta-to-node blend reading past the old buffer.
- Flame's elastic psi path now keeps the base mechanics `psi_on` flag disabled
  while attaching Flame's eta-derived `psi_mf` directly to the solver when
  `elastic.use_psi=1`. That preserves the intended `elastic.use_psi=0`
  behavior and avoids reattaching the empty base mechanics psi field.
- GPU builds now reject `elastic.type=dynamic` early because that path captures
  host boundary-condition polymorphism inside device kernels.

## Test hardening

- `tests/GPU/F2_smoke_elastic` now exercises `elastic.traction_from_chamber=1`
  with the same static 1 MPa chamber pressure load.
- `tests/GPU/C1_correctness_elastic` now fails if required physics/elastic
  thermo columns are missing from either CPU or GPU output.
- `tests/GPU/C2_restart_roundtrip` now asserts that restarted
  `chamber_pressure` matches the continuous run at checkpoint time.
- `tests/GPU/C4_amr_correctness` now fails if required `eta` or `temp` plotfile
  fields are missing instead of warning and skipping the comparison.

## Verification

- `make -j2 bin/alamo_gpu-3d-nofast-cuda86-g++`
- `make -j2 POSTFIX=2d-nofast-cuda86-g++ AMREX_TARGET=ext/AMReX-Codes/amrex/2d-nofast-cuda86-g++-26.06-dirty bin/alamo_gpu-2d-nofast-cuda86-g++`
- `make -j2 POSTFIX=2d-g++ AMREX_TARGET=ext/AMReX-Codes/amrex/2d-g++-26.06-dirty CC=mpicxx COMP_CMD='mpicxx -c' LINK_CMD=mpicxx bin/alamo-2d-g++`
- `tests/GPU/F2_smoke_elastic/test.py` passed with the 2D nofast GPU binary.
- `tests/GPU/C1_correctness_elastic/test.py` passed with the 2D CPU and 2D
  nofast GPU binaries.
- `tests/GPU/C2_restart_roundtrip/test.py` passed with the 2D nofast GPU binary.
- `tests/GPU/F1_smoke_flame_only/input` passed through the GPU test harness with
  `max_step=20` using the 2D nofast GPU binary.

The first concurrent C1/F2/C2/F1 run exhausted GPU memory; each focused rerun
above passed when run alone.

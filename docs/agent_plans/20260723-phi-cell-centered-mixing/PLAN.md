# TASK: phi-cell-centered-mixing

## Header

| Field | Value |
|---|---|
| Risk tier | 3 (field centering changes AMR, mechanics coefficients, and restart layout) |
| Verification | paired visual experiment plus unit/build gates |
| Scope | `Flame.H`, `Flame.cpp`, task-owned run artifacts |

## Objective

Add an opt-in experiment that treats `phi` like `eta` for mechanics mixing:
store it cell-centered and use `CellToNodeAverage(phi)` when constructing the
nodal elastic model. Preserve the existing nodal behavior as the default.
Produce a run that can be compared in VisIt without simultaneously changing
the material interpolation law.

## Controls

- Do not change the arithmetic material weights.
- Do not change `eta`, pressure loading, elastic tolerances, or AMR criteria.
- Give cell-centered `phi` explicit physical boundary conditions; do not reuse
  the numerically different `eta` boundary values.
- Preserve all pre-existing dirty changes and record source hashes.
- Treat node- and cell-centered checkpoints as layout-incompatible.

## Steps

1. Add `phi.cell_centered_mixing=0` by default.
2. When enabled, register `phi` as a three-ghost cell field with `phi.bc`;
   otherwise retain the current two-ghost nodal field.
3. Branch only at centering-dependent consumers:
   - mechanics: direct nodal value versus `CellToNodeAverage`;
   - cell physics: `NodeToCellAverage` versus direct cell value;
   - nodal synchronization versus cell ghost fill.
4. Build `bin/test-2d-g++` and `bin/alamo-2d-g++`; run the 2-D unit suite.
5. Run a short paired case first. If stable, produce the requested
   VisIt-comparison output.
6. Record exact commands, hashes, and verdict. Retain the opt-in knob only if
   tests pass; do not make it the default based on visual evidence alone.

## Acceptance

- default-off behavior is unchanged;
- the experimental run completes without NaN/Inf or unused inputs;
- `phi` appears in the cell plot and its nodal mechanics reconstruction is
  finite on every AMR level;
- unit suite, `git diff --check`, device lint, golden compare, and A100
  sanitizer pass.


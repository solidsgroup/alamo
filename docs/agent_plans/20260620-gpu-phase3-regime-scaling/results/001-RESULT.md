# Result: Task 001 — 3D wide-shallow Flame input files

## Summary
Authored the canonical 3D GPU-resident Flame production input `input_3d_flame` plus
three crossover-sweep grid variants (`_128`, `_256`, `_512`) at the repo root. All four
mirror the proven, GPU-verified `input_3d_smoke` template's IC/BC/elastic structure
exactly (Expression ICs for `pf.eta` and `phi`, full six-face eta/temp Dirichlet BCs,
`elastic.type = disable` with NO elastic model/bc block), and carry the full
propellant/thermal/chamber physics blocks. Only the grid was scaled up to production
sizes following the wide-shallow box strategy (roadmap 3.3): `amr.max_level = 1`,
`amr.blocking_factor = 32`, `amr.max_grid_size = 128`, large base grid, thin z, and
`amr.node.all = 1`. `plot_file` is per-variant, `amr.plot_int = -1`, `amr.thermo.int = 1`
(sweeps override step counts on the command line).

## Files changed (created only)
- `input_3d_flame`     — canonical production input; `amr.n_cell = 256 256 128`,
  `geometry.prob_hi = 0.1754_m 0.1754_m 0.0877_m`, `plot_file = output_3d_flame`.
- `input_3d_flame_128` — `amr.n_cell = 128 128 64`,  `plot_file = output_3d_flame_128`.
- `input_3d_flame_256` — `amr.n_cell = 256 256 128`, `plot_file = output_3d_flame_256`
  (duplicates the canonical grid by design, per task step 5).
- `input_3d_flame_512` — `amr.n_cell = 512 512 256`, `plot_file = output_3d_flame_512`.

No existing/shared file was modified (`input`, `input_3d_smoke`, `src/`, `benchmark/`,
`bin/`, and task 002/003 files were left untouched).

## Tests run and results (no build, no GPU — structural validation only)
- `grep -n "bmp" input_3d_flame*` → no matches (exit 1). PASS — no BMP ICs.
- `grep -nE "elastic\.bc|^model_(prop|void|casing)" input_3d_flame*` → no matches
  (exit 1). PASS — no elastic model/bc blocks (strict-parser safe with `type=disable`).
- `elastic.type = disable` present in all four files (non-comment line). PASS.
- All field `*.ic.type` are `expression`, except `temp.ic.type = constant` — which is
  identical to the smoke template (temp uses a constant IC there too, not an Expression).
  PASS (matches template verbatim).
- Each file has exactly the six distinct BC faces (xlo, xhi, ylo, yhi, zlo, zhi) for
  both `pf.eta.bc.constant.type.*` and `thermal.temp.bc.constant.type.*`. PASS.
- Grids confirmed: `_128`=128 128 64, canonical & `_256`=256 256 128, `_512`=512 512 256;
  all with `max_level=1`, `blocking_factor=32`, `max_grid_size=128`. PASS.
- Propellant (`propellant.homogenize.*`, 16 keys), chamber (`chamber.ballistic.*`,
  5 keys), and `thermal.on = 1` blocks present and intact. PASS.
- Geometry wide-shallow: `prob_hi = 0.1754 0.1754 0.0877` (thin z). PASS.

## Issues found
None. The smoke template existed (stop condition not triggered) and was mirrored
cleanly. No new ParmParse keys were invented; `amr.max_grid_size` is the only key not
explicitly present in `input_3d_smoke`, but it is a standard AMReX/ALAMO key already
used in `input` and mandated by the task (step 2), so it is in-contract.

## Deviations from task
None. All four files created as specified; grid, box strategy, ICs, BCs, and elastic
disposition match the task and the template.

## Follow-up needed
- Lead-owned: actual 3D build + GPU smoke/run of `input_3d_flame` to confirm runtime
  acceptance (workers do not build/run CUDA). The `_512` variant (~33.6M base cells)
  will exceed the local A1000's 8 GB and is intended for NOVA A100/H200.
- Sweep driver (task 003) references these inputs by name; naming scheme is uniform.

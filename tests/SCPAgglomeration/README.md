# Singly degenerate Cahn–Hilliard agglomeration cases

These two proposed Flame cases put a nanopowder Al concentration field in the
HTPB-rich space around a random bimodal AP packing. **Neither Alamo case has
been run or validated yet.** Review the inputs, [2D IC sketch](ic_preview_2d.svg),
[3D periodic sections](ic_preview_3d.svg), and [3D axial section with gas](ic_preview_3d_axial.svg)
before using either one. The older
`input` and `2d.xyzr` are retained unchanged.

| Case | Input | AP geometry | Base mesh | Intended first run |
| --- | --- | --- | --- | --- |
| 2D slice | `input_2d` | `bimodal_2d.xyzr` | 168 × 96 | 20 steps at 1 µs |
| 3D parent | `input_3d` | `bimodal_3d.xyzr` | 96 × 96 × 168 | 20 steps at 1 µs |

## Composition and IC interpretation

The formulation in `~/code/convert_mass_to_volume/mass_to_vol.py` is 68 wt% AP,
17 wt% HTPB, and 15 wt% Al, with densities 1.95, 0.92, and 2.70 g/cm³.
The implied AP fraction of the entire solid is about 0.592. The Al fraction
*within the Al–HTPB mixture* is

```text
c_Al = (15/2.70) / ((15/2.70) + (17/0.92)) = 0.2311.
```

`agglomeration.alpha` is a Cahn–Hilliard order parameter on `[-1, 1]`, with
the local Al fraction represented by `(alpha + 1)/2`. `IC::Random` adds
`offset + mult*U`, where `U` is uniform on `[0, 1)`. The new inputs use
`offset = -0.60` and `mult = 0.1244`, so the *unrestricted* local Al fraction
is random on `[0.20, 0.2622)`, centered on 0.2311. This is a nanopowder
concentration field, not resolved Al particles. The random realization is
selected by AMReX at initialization and is not depicted in the SVG.

In `Flame::InitializeAgglomeration`, the seed is then restricted by
`eta * (1 - phi)`, where `eta` is unburned propellant and `phi` is AP.
Thus AP and the gas padding begin with essentially no Al, while HTPB-rich
cells retain the local random concentration. The plotted `agglom.alpha`
therefore also reflects the binder mask; converting its domain-wide mean
directly to an Al/HTPB ratio would be misleading. The intended 0.2311 value
applies in fully unburned binder away from diffuse AP boundaries.

## AP geometry

`generate_packing.py` uses the bimodal sizes and volume targets from
`~/code/packing_fr` and `~/code/convert_mass_to_volume`. It grows 36 coarse
and 72 fine spheres from random centers and relaxes all overlaps using
minimum-image distances in x/y. Their centers and radii stay inside the
z slab. The fixed seed is 8402; NumPy and SciPy are needed to regenerate
the checked-in XYZR and SVG files. This geometry generation has been run;
Alamo has not.

`bimodal_3d.xyzr` is a 512 µm cube with 36 coarse spheres (mean diameter
149.37 µm) and 72 fine spheres (mean diameter 75.05 µm). Their sharp-sphere
volume fractions are 0.4721 and 0.1199, totaling **0.5920 AP** with the
intended roughly 80/20 coarse/fine volume split. The packing is genuinely
random in position and repeats across x/y. The checked-in data have no
meaningful periodic pair overlaps; the largest numerical overlap is under
0.00002 µm. Every sphere lies inside `0 ≤ z ≤ 512 µm`.

`bimodal_2d.xyzr` is the `x = 179 µm` cross-section of that 3D packing,
mapped from parent `(z, y)` to 2D `(x, y)`. It has 20 circles, repeats across
y, and has sharp-circle AP area fraction **0.5922**. Sectioned diameters
differ from their parent sphere diameters. The diffuse `phi` interface will
change the measured AP fraction in both Alamo cases. At 0.592 AP loading,
the 0.2311 Al fraction of the remaining Al–HTPB mixture corresponds to about
0.094 of the whole solid by volume.

The 2D input has 192 µm of gas at the -x and +x faces, and its solid slab is
`0 ≤ x ≤ 512 µm`; y is periodic. The 3D input has 192 µm of gas at the -z
and +z faces, and its solid slab is `0 ≤ z ≤ 512 µm`; x and y are periodic.
The full computational boxes are periodic for the spectral Cahn–Hilliard
solver, with gas separating the two burn faces. The XYZR files use
micrometers, declared by `phi.ic.psread.file.unit = um`.

`IC::PSRead` now applies the minimum image in every periodic direction so
spheres crossing an x/y corner are complete. This source change and the
resulting IC fields still need compilation and plotfile inspection after
the cases are approved for a run.

## Model and review points

Both inputs use the current Flame full-feedback thermal/propellant settings
and the singly-degenerate spectral Cahn–Hilliard settings from the older
`tests/SCPAgglomeration/input`. The `eta_mobility_cutoff = 0.5` mask leaves
Al effectively immobile in unburned material and allows evolution as the
burn front passes. The 500 K initial temperature and 4 MPa pressure are
exploratory conditions. A 5 µs output interval would give an initial and
several early snapshots in a 20 µs run. No claim is made yet about numerical
stability, burning distance, Al conservation, or visible agglomeration over
that short interval.

Before running, review the IC sketches, the AP geometry,
the Al range, and the thermal/burn settings. When authorized to run, copy the
source tree to a temporary Alamo directory without simulation data, configure
with FFT and HDF5 in `alamo-llvm`, use the matching `*-hdf5-*` binary, and
inspect the actual initial plotfile. The SVG is a geometric sketch, not a
render of Alamo's initialized fields.

Run from the temporary repository root so the relative XYZR paths resolve.

# Singly degenerate Cahn–Hilliard agglomeration cases

The two bimodal-packing Flame cases put a nanopowder Al concentration field in
the HTPB-rich space around a random AP packing. The 2D case has been run on
Nova through 200 ms; the 3D case has not been run. The separate
`input` case uses a larger 1 × 1 mm 2D slab and a 96-circle packing based on
the original `2d.xyzr` size distribution. Review the inputs,
[2D IC sketch](ic_preview_2d.svg),
[3D periodic sections](ic_preview_3d.svg), and [3D axial section with gas](ic_preview_3d_axial.svg)
before using the bimodal cases. The original `2d.xyzr` remains unchanged.

| Case | Input | AP geometry | Base mesh | Intended first run |
| --- | --- | --- | --- | --- |
| 2D slice | `input_2d` | `bimodal_2d.xyzr` | 896 × 512 | 200 ms at 1 µs/step; plots every 2 ms |
| 3D parent | `input_3d` | `bimodal_3d.xyzr` | 96 × 96 × 168 | 2 ms at 1 µs/step; plots every 0.2 ms |
| Expanded legacy 2D | `input` | `2d_large.xyzr` | 320 × 256 | 75 ms at 1 µs/step; plots every 0.5 ms |
| Lower-left interface zoom, short | `input_interface_zoom_2026ses` | `2d_large.xyzr` | 320 × 256 | 10 ms at 1 µs/step; plots every 10 µs |
| Lower-left interface zoom, extended | `input_interface_zoom_extended_2026ses` | `2d_large.xyzr` | 320 × 256 | 50 ms at 1 µs/step; plots every 25 µs |

The expanded 2D run is retained in
`flagship_2d_result_random_alpha_2026ses/`, including the copied `input_used`
and matching `2d_large.xyzr`. The current `input` writes future runs to that
directory. Its 151 cell-centered HDF5 frames are rendered by
`../../scripts/render_scp_agglomeration_video_2026ses.py` at 30 fps, with a
tight crop, fixed 0–1 Jet scale and α legend, and no time overlay. It writes a
transparent ProRes 4444 MOV and looping transparent GIF, plus a tightly cropped
white-background H.264 MP4 for players that do not support alpha video. The uv
project declares h5py, Pillow, PyVista/VTK, Matplotlib, and ImageIO; FFmpeg and
gifski handle encoding. Re-render from the repository root with PyVista's
GUI-backed renderer, FFmpeg, and gifski available:

```sh
UV_PROJECT_ENVIRONMENT="$PWD/.venv" uv run --project . \
  python scripts/render_scp_agglomeration_video_2026ses.py
```

`input_interface_zoom_2026ses` retains the expanded case's geometry and
physical/numerical settings; only the output path, `amr.plot_dt` (10 µs), and
`stop_time` (10 ms) differ. Its 1,001 HDF5 snapshots, copied input, and matching
geometry are in `flagship_2d_interface_zoom_2026ses/`. The transparent GIF
zooms to x = −0.03…0.20 mm and y = 0…0.24 mm at the lower-left burn front. The
render samples every fourth saved state for a compact 30 fps loop:

```sh
UV_PROJECT_ENVIRONMENT="$PWD/.venv" uv run --project . \
  python scripts/render_scp_agglomeration_video_2026ses.py \
  --plot-dir tests/SCPAgglomeration/flagship_2d_interface_zoom_2026ses \
  --output tests/SCPAgglomeration/flagship_2d_interface_zoom_2026ses/flagship_2d_interface_zoom_2026ses.mp4 \
  --roi-mm -0.03 0.20 0 0.24 --frame-stride 4
```

The extended run retains the same geometry and model settings, with a 50 ms
stop time and 25 µs plot interval (2,001 states). Its GIF uses the slightly
wider x = −0.03…0.32 mm, y = 0…0.24 mm crop and every eighth saved state, so the
burn front travels visibly through the field while the loop stays compact:

```sh
UV_PROJECT_ENVIRONMENT="$PWD/.venv" uv run --project . \
  python scripts/render_scp_agglomeration_video_2026ses.py \
  --plot-dir tests/SCPAgglomeration/flagship_2d_interface_zoom_extended_2026ses \
  --output tests/SCPAgglomeration/flagship_2d_interface_zoom_extended_2026ses/flagship_2d_interface_zoom_extended_2026ses.mp4 \
  --roi-mm -0.03 0.32 0 0.24 --frame-stride 8
```

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
`offset + mult*U`, where `U` is uniform on `[0, 1)`. The bimodal `input_2d`
and `input_3d` cases use `offset = -0.60` and `mult = 0.1244`, giving an
unrestricted fraction on `[0.20, 0.2622)`, centered on 0.2311. The expanded
legacy `input` retains the original full-range field, with `offset = -1` and
`mult = 2`, so its unrestricted fraction spans `[0, 1)`. These are nanopowder
concentration fields, not resolved Al particles. The random realization is
selected by AMReX at initialization and is not depicted in the SVG.

In `Flame::InitializeAgglomeration`, the seed is then restricted by
`eta * (1 - phi)`, where `eta` is unburned propellant and `phi` is AP.
Thus AP and the gas padding begin with essentially no Al, while HTPB-rich
cells retain the local random concentration. The plotted `agglom.alpha`
therefore also reflects the binder mask; converting its domain-wide mean
directly to an Al/HTPB ratio would be misleading. The intended 0.2311 value
applies to the bimodal cases in fully unburned binder away from diffuse AP
boundaries. The expanded legacy input intentionally uses the original
full-range Al field instead.

## AP geometry

`generate_packing.py` uses the bimodal sizes and volume targets from
`~/code/packing_fr` and `~/code/convert_mass_to_volume`. It grows 36 coarse
and 72 fine spheres from random centers and relaxes all overlaps using
minimum-image distances in x/y. Their centers and radii stay inside the
z slab. The fixed seed is 8402; NumPy and SciPy are needed to regenerate
the checked-in XYZR and SVG files. This geometry generation is complete; the
3D Alamo case has not been run.

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
spheres crossing an x/y corner are complete. The 2D case compiled and wrote
HDF5 plotfiles; the 3D case still needs validation.

## Model and review points

All three inputs use the current Flame full-feedback thermal/propellant
settings and singly-degenerate spectral Cahn–Hilliard settings. The
`eta_mobility_cutoff = 0.5` mask leaves
Al effectively immobile in unburned material and allows evolution as the
burn front passes. The 500 K initial temperature and 4 MPa pressure are
exploratory conditions. The bimodal 2D case is configured for 101 plot times,
the expanded legacy 2D case for 151, and the 3D case for 11, including their
initial states. The expanded case retains the original `input` SDCH settings
(`L = 10`, clipping on, pointwise restoration off) and its full-range Al
initialization; its larger-domain behavior still needs inspection. The shorter
3D run limits the cost of a first full-volume pilot.
No claim is made yet about numerical stability, burning distance, Al
conservation, or visible agglomeration in the expanded 2D and 3D cases.

### Previous bimodal 2D comparison: kinetics rejected

Nova task `16854825_0` used the current `input_2d` parameters: 1 µm cells,
1 µs steps, `gamma=5e-7`, `L=0.5`, mobility floor `0.01`, eta cutoff `0.5`,
and clipping/restoration off. The geometry and IC definitions are unchanged.
Human review rejected its slow mobility: droplets form too long after the burn
front passes, and remain too small relative to their alpha interface thickness.
The current preset is retained for comparison while faster kinetics and better
interface resolution are investigated.
At 200 ms it has 88 connected alpha>0.9 cores and a median alpha>0 body
diameter of 15.4 µm. Overlap tracking identifies ten joining events; rendered
sequences confirm lasting joins. By 80 ms there are six tracked joins,
compared with none for the earlier same-grid/gamma/L run with floor `1e-6`.

The signed Al integral is conserved to about 6e-9 percent in saved frames;
bounding `(alpha+1)/2` to [0,1] changes the final measure by -0.315%.
Saved alpha extrema are -1.0069 and 1.0360, with no alpha>0.9 cells behind
eta>=0.5. Small changes near the cutoff remain: the spectral mask does not
freeze every protected cell exactly. Droplets remain larger than the original
1–5 µm target; this is a configuration for visual review, not a calibrated
physical model. Review videos use PyVista/ffmpeg, jet Al fraction and the
white eta=0.5 contour.

Before running, review the IC sketches, the AP geometry,
the Al range, and the thermal/burn settings. When authorized to run, copy the
source tree to a temporary Alamo directory without simulation data, configure
with FFT and HDF5 in `alamo-llvm`, use the matching `*-hdf5-*` binary, and
inspect the actual initial plotfile. The SVG is a geometric sketch, not a
render of Alamo's initialized fields.

Run from the temporary repository root so the relative XYZR paths resolve.

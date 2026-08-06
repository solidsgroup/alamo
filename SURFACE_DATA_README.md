# Cyclone burn-surface data

`extract_burn_surface.py` extracts reduced-order regression data from every
`*cell` plotfile in this directory. It samples each connected
`rigid_eta = 0.5` contour at approximately one finest-level cell spacing and
writes one table per plotfile to `surface_data/`.

Run it from the Alamo repository root with:

```bash
python3 output.lm.ap_htpb_cyclone_half_fast/extract_burn_surface.py
```

The script reads material and transport properties from `metadata`. It streams
the plotfiles and retains only neighboring snapshots in memory.

Render heat and mass flux along every extracted contour with:

```bash
python3 output.lm.ap_htpb_cyclone_half_fast/render_surface_flux.py --jobs 8
```

Render the same fluxes against signed distance from the nearest material
interface with:

```bash
python3 output.lm.ap_htpb_cyclone_half_fast/render_interface_distance_flux.py --jobs 8
```

Render heat flux divided by regression mass flux against the same signed
distance with:

```bash
python3 output.lm.ap_htpb_cyclone_half_fast/render_interface_heat_per_mass.py --jobs 8
```

Superimpose every output in one signed-distance figure with:

```bash
python3 output.lm.ap_htpb_cyclone_half_fast/render_all_interface_heat_per_mass.py
```

## Definitions

The unit normal points from condensed material into the gas:

```text
n = -grad(eta) / |grad(eta)|.
```

Except at the first and last output, eta time derivatives use centered
differences between the adjacent plotfiles. The positive regression mass flux
is

```text
j = -(rho_AP eta_AP_dot + rho_HTPB eta_HTPB_dot) / |grad(eta)|.
```

This uses the two species-resolved rigid eta fields, so it remains meaningful
at AP/HTPB transitions. The first and last outputs use one-sided differences.

The local conductivity reproduces the LowMach thermal mixture rule using the
Rocfire gas conductivity and AP/HTPB conductivities recorded in `metadata`.
The reported normal heat flux is positive from the gas into the solid:

```text
q_into_solid = k grad(T) dot n.
```

`heat_flux_outward_W_m2` stores the opposite sign.

The local condensed composition is
`AP_solid_fraction = eta_AP / (eta_AP + eta_HTPB)`. `species` is `AP` at and
above a fraction of 0.5 and `HTPB` below it. An AP/HTPB interface is a crossing
of this fraction through 0.5 along an ordered connected burn contour.
`distance_to_interface_m` is the shortest arc-length distance along that
contour, not the Euclidean distance. `signed_distance_to_interface_m` is
positive in AP and negative in HTPB.

## Outputs

- `surface_data/<plotfile>_surface.csv`: one row per contour sample.
- `surface_data_manifest.csv`: time, table name, contour/point counts, and
  mean fluxes for each plotfile.
- `surface_validation/*.png`: representative contour overlays. Blue samples
  are HTPB, red samples are AP, and white rings mark material interfaces.
- `surface_extraction.log`: extraction progress and any error message.
- `surface_flux_plots/*_flux.png`: heat flux into the solid and regression
  mass flux versus contour arc length for every plotfile. AP points are red,
  HTPB points are blue, and dotted vertical lines mark material transitions.
- `surface_flux_render.log`: flux-plot rendering progress and errors.
- `surface_interface_flux_plots/*_interface_flux.png`: scatter plots against
  signed interface distance. HTPB distance is negative and AP distance is
  positive. Scatter points are used because several distinct material
  interfaces can contribute values at the same distance.
- `surface_interface_flux_render.log`: signed-distance rendering progress and
  errors.
- `surface_interface_heat_per_mass_plots/*_heat_per_mass.png`: signed
  `heat_flux_into_solid / mass_flux` in MJ/kg against interface distance.
  HTPB distance is negative and AP distance is positive.
- `surface_interface_heat_per_mass_render.log`: heat-per-mass rendering
  progress and errors.
- `surface_interface_heat_per_mass_all.png`: all 667,110 heat-per-mass
  samples superimposed on linear axes. The default displayed range is 0 to
  5 MJ/kg; change it with `--ymin` and `--ymax`.
- `surface_interface_heat_per_mass_density.png`: the same linear-axis data,
  with each point colored by its local two-dimensional bin population using
  the `jet` colormap. Density resolution is configurable with
  `--density-xbins` and `--density-ybins`.
- `surface_interface_heat_per_mass_all.log`: aggregate-plot loading progress.

The per-plotfile tables include coordinates, contour and point identifiers,
arc length, normal, eta gradients and rates, total and species-resolved mass
fluxes, temperature, pressure, conductivity, normal temperature gradient,
both heat-flux sign conventions, condensed composition, species, reference
density, and unsigned/signed distance to the nearest material interface.

The present extraction covers 3,445 plotfiles from 0 to
0.17219977427749195 s and contains 667,110 surface samples.

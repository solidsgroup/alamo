#!/usr/bin/env python3
"""Build correlated 2-D low-Mach flame initial-condition rasters.

The saved gas and condensed fields remain registered to their developed
propellant surface, are extended with uniform products above the particle-free
sampled region, and are then combined with the requested cold aluminum branch
at its original location.
"""

import argparse

import numpy as np
from PIL import Image
from scipy.interpolate import RegularGridInterpolator
import yt


parser = argparse.ArgumentParser()
parser.add_argument("plotfile")
parser.add_argument("--prefix", default="ap_htpb_developed_flame")
parser.add_argument(
    "--particle-velocity-fraction", type=float, default=0.0,
    help="initial aluminum velocity divided by the local developed-flow "
         "velocity (default: 0.0)")
args = parser.parse_args()
if not np.isfinite(args.particle_velocity_fraction):
    parser.error("--particle-velocity-fraction must be finite")

ds = yt.load(args.plotfile)
# Avoid a yt floating-point boundary check that can place the last covering-
# grid face a few ulps beyond the declared nonperiodic upper y boundary.  The
# requested cells themselves remain exactly inside the saved domain.
ds.force_periodicity()
level = ds.max_level
dims = ds.domain_dimensions * ds.refine_by**level
grid = ds.covering_grid(level, ds.domain_left_edge, dims)

lo = np.asarray(ds.domain_left_edge[:2], dtype=float)
hi = np.asarray(ds.domain_right_edge[:2], dtype=float)
dx = (hi - lo) / dims[:2]
xc = lo[0] + (np.arange(dims[0]) + 0.5) * dx[0]
yc = lo[1] + (np.arange(dims[1]) + 0.5) * dx[1]
xx, yy = np.meshgrid(xc, yc, indexing="ij")

field_names = [
    "component_density_AP_gas",
    "component_density_HTPB_gas",
    "component_density_Mono",
    "component_density_Premixed",
    "component_density_Primary",
    "component_density_Final",
    "temperature",
    "velocityx",
    "velocityy",
    "pressure",
    "component_density_AP_solid",
    "component_density_HTPB_solid",
]
fields = {
    name: np.asarray(grid[("boxlib", name)][:, :, 0], dtype=float)
    for name in field_names
}
eta = (
    np.asarray(grid[("boxlib", "component_density_AP_solid")][:, :, 0],
               dtype=float) / 1950.0
    + np.asarray(grid[("boxlib", "component_density_HTPB_solid")][:, :, 0],
                 dtype=float) / 920.0
)

# Select the eta=0.5 crossing nearest the nominal propellant surface.  This
# rejects unrelated crossings through discrete AP particles below the surface.
surface = np.empty(dims[0])
for i in range(dims[0]):
    crossing = np.flatnonzero((eta[i, :-1] >= 0.5) & (eta[i, 1:] < 0.5))
    if not crossing.size:
        raise RuntimeError(f"no propellant surface crossing at x index {i}")
    j = crossing[np.argmin(np.abs(0.5 * (yc[crossing] + yc[crossing + 1])
                                       + 3.0e-4))]
    fraction = (0.5 - eta[i, j]) / (eta[i, j + 1] - eta[i, j])
    surface[i] = yc[j] + fraction * (yc[j + 1] - yc[j])

# Periodic interpolation in x and ordinary interpolation in y.
period = hi[0] - lo[0]
xext = np.concatenate(([xc[-1] - period], xc, [xc[0] + period]))
xwrapped = np.mod(xx - lo[0], period) + lo[0]
height = yy - np.mean(surface)
# Keep the developed flame and propellant state in their saved coordinates.
# Translating only the gas would break their diffuse-interface balance, while
# warping columns independently would alter the flame's spatial derivatives.
source_y = np.minimum(yy, np.mean(surface) + 5.0e-4)
source_y = np.clip(source_y,
                   yc[0], yc[-1])
points = np.stack((xwrapped, source_y), axis=-1)

sampled = {}
for name, data in fields.items():
    extended = np.concatenate((data[-1:, :], data, data[:1, :]), axis=0)
    sampled[name] = RegularGridInterpolator(
        (xext, yc), extended, bounds_error=False,
        fill_value=None)(points)

# The late plot is particle-free up to 0.5 mm above the propellant.  Blend
# from that correlated state to uniform fully burned products by 0.8 mm rather
# than copying the translated particle and its wake from the reference run.
blend = np.clip((height - 5.0e-4) / 3.0e-4, 0.0, 1.0)
blend = blend * blend * (3.0 - 2.0 * blend)
for name in field_names[:6]:
    sampled[name] *= 1.0 - blend
sampled["component_density_Final"] += blend
sampled["temperature"] = ((1.0 - blend) * sampled["temperature"]
                          + blend * 2553.85363)
sampled["velocityx"] *= 1.0 - blend
sampled["velocityy"] = ((1.0 - blend) * sampled["velocityy"]
                        + blend * 1.70013460)
sampled["pressure"] = ((1.0 - blend) * sampled["pressure"]
                       + blend * 3.0e6)
# Low-Mach pressure is a 3 MPa thermodynamic reference plus a zero-mean
# hydrodynamic correction.  Retain the developed gradients without shifting
# the thermodynamic pressure inferred from the initial field average.
sampled["pressure"] += 3.0e6 - np.mean(sampled["pressure"])

propellant = np.zeros((*xx.shape, 4), dtype=np.uint8)
ap_raw = np.rint(65535.0 * np.clip(
    sampled["component_density_AP_solid"] / 1950.0, 0.0, 1.0)).astype(np.uint16)
htpb_raw = np.rint(65535.0 * np.clip(
    sampled["component_density_HTPB_solid"] / 920.0, 0.0, 1.0)).astype(np.uint16)
propellant[..., 0] = ap_raw >> 8
propellant[..., 1] = ap_raw & 255
propellant[..., 2] = htpb_raw >> 8
propellant[..., 3] = htpb_raw & 255
Image.fromarray(np.transpose(propellant, (1, 0, 2)), "RGBA").save(
    f"{args.prefix}_propellant.png")

# Convert the six gas densities to composition seeds and use 16-bit integer
# largest-remainder quantization.  Thus all six encoded values sum exactly to
# 65535 at every node and bilinear sampling preserves a unit mass-fraction sum.
gas = np.stack([sampled[name] for name in field_names[:6]], axis=-1)
gas = np.maximum(gas, 0.0)
gas /= np.maximum(np.sum(gas, axis=-1, keepdims=True), 1.0e-300)
raw = 65535.0 * gas
quantized = np.floor(raw).astype(np.uint16)
remainder = 65535 - np.sum(quantized, axis=-1, dtype=np.int32)
order = np.argsort(-(raw - quantized), axis=-1)
for rank in range(5):
    mask = remainder > rank
    ii, jj = np.nonzero(mask)
    quantized[ii, jj, order[ii, jj, rank]] += 1

species0 = np.stack((quantized[..., 0] >> 8,
                     quantized[..., 0] & 255,
                     quantized[..., 1] >> 8,
                     quantized[..., 1] & 255), axis=-1).astype(np.uint8)
species1 = np.stack((quantized[..., 2] >> 8,
                     quantized[..., 2] & 255,
                     quantized[..., 3] >> 8,
                     quantized[..., 3] & 255), axis=-1).astype(np.uint8)
species2 = np.stack((quantized[..., 4] >> 8,
                     quantized[..., 4] & 255,
                     quantized[..., 5] >> 8,
                     quantized[..., 5] & 255), axis=-1).astype(np.uint8)
Image.fromarray(np.transpose(species0, (1, 0, 2)), "RGBA").save(
    f"{args.prefix}_species_0.png")
Image.fromarray(np.transpose(species1, (1, 0, 2)), "RGBA").save(
    f"{args.prefix}_species_1.png")
Image.fromarray(np.transpose(species2, (1, 0, 2)), "RGBA").save(
    f"{args.prefix}_species_2.png")

# Restore the current aluminum branch at its original location.  Define its
# smoothed outer signed distance once so the thermal and velocity fields cannot
# silently diverge from the initialized liquid geometry.
particle_lobes = [(-95.0e-6, 0.0, 97.0e-6),
                  (35.0e-6, 55.0e-6, 93.0e-6),
                  (82.0e-6, 177.0e-6, 85.0e-6)]
outer_distances = np.stack([
    np.hypot(xx - x0, yy - y0) - radius
    for x0, y0, radius in particle_lobes
])
union_smoothing = 40.0e-6
minimum_distance = np.min(outer_distances, axis=0)
particle_distance = minimum_distance - union_smoothing * np.log(np.sum(
    np.exp(-(outer_distances - minimum_distance) / union_smoothing), axis=0))

# Resolve the established thermal stand-off on the same raster as the
# developed flame, measured normally from the actual smoothed particle.
thermal_standoff = 60.0e-6
thermal_transition_width = 40.0e-6
thermal_mask = 0.5 + 0.5 * np.tanh(
    2.0 * (particle_distance - thermal_standoff) /
    thermal_transition_width)
sampled["temperature"] = (930.0 +
    (sampled["temperature"] - 930.0) * thermal_mask)

# The complete current aggregate, including its enlarged melt layer and
# rounded lobe junctions, starts at rest in the developed gas.  Reproduce the
# same smooth outer union and 20 um diffuse profile used by Al_liquid in the
# input.  Hold the resolved diffuse particle at the requested velocity, then
# use a compact smooth transition on the gas side so the initial slip has no
# grid-scale velocity jump or preferred direction.
interface_width = 20.0e-6
velocity_rest_standoff = 2.0 * interface_width
velocity_transition_width = 2.0 * interface_width
gas_transition = np.clip(
    (particle_distance - velocity_rest_standoff) /
    velocity_transition_width, 0.0, 1.0)
gas_transition = gas_transition * gas_transition * (3.0 - 2.0 * gas_transition)
velocity_scale = (args.particle_velocity_fraction +
                  (1.0 - args.particle_velocity_fraction) * gas_transition)
sampled["velocityx"] *= velocity_scale
sampled["velocityy"] *= velocity_scale

temperature = np.zeros((*xx.shape, 4), dtype=np.uint8)
temperature_raw = np.rint(65535.0 * np.clip(
    (sampled["temperature"] - 300.0) / 2700.0, 0.0, 1.0)).astype(np.uint16)
temperature[..., 0] = temperature_raw >> 8
temperature[..., 1] = temperature_raw & 255
temperature[..., 3] = 255
Image.fromarray(np.transpose(temperature, (1, 0, 2)), "RGBA").save(
    f"{args.prefix}_temperature.png")

velocity = np.zeros((*xx.shape, 4), dtype=np.uint8)
velocity_x_raw = np.rint(65535.0 * np.clip(
    (sampled["velocityx"] + 20.0) / 40.0, 0.0, 1.0)).astype(np.uint16)
velocity_y_raw = np.rint(65535.0 * np.clip(
    (sampled["velocityy"] + 20.0) / 50.0, 0.0, 1.0)).astype(np.uint16)
velocity[..., 0] = velocity_x_raw >> 8
velocity[..., 1] = velocity_x_raw & 255
velocity[..., 2] = velocity_y_raw >> 8
velocity[..., 3] = velocity_y_raw & 255
Image.fromarray(np.transpose(velocity, (1, 0, 2)), "RGBA").save(
    f"{args.prefix}_velocity.png")

pressure = np.zeros((*xx.shape, 4), dtype=np.uint8)
pressure_raw = np.rint(65535.0 * np.clip(
    (sampled["pressure"] - 2.5e6) / 1.5e6, 0.0, 1.0)).astype(np.uint16)
pressure[..., 0] = pressure_raw >> 8
pressure[..., 1] = pressure_raw & 255
pressure[..., 3] = 255
Image.fromarray(np.transpose(pressure, (1, 0, 2)), "RGBA").save(
    f"{args.prefix}_pressure.png")

print(f"surface range: {surface.min():.9e} to {surface.max():.9e} m")
print("mapped ranges:")
for name in ("temperature", "velocityx", "velocityy", "pressure"):
    print(f"  {name}: {sampled[name].min():.8g} to {sampled[name].max():.8g}")

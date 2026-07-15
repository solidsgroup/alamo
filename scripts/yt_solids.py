#!/usr/bin/env python3
"""Render 2-D LowMach/BoxLib solid diagnostics with yt and Matplotlib.

Each output frame contains one composite plot: z-vorticity in the fluid and
2-D von Mises stress overlaid where eta exceeds the interface threshold.

The plot includes the eta interface and masked reference-map (xix/xiy)
isolines.  Massless visualization tracers are advanced between plotfiles with
a midpoint rule and bilinear interpolation of the velocity FRB, and are hidden
inside the solid.

Examples
--------
Render every plotfile under output.lm using the finest AMR resolution::

    python scripts/yt_solids.py output.lm

Render every fifth plotfile from step 3000 onward and make a movie::

    python scripts/yt_solids.py output.lm --start 3000 --stride 5 \
        --movie output_lm.mp4

Render only the newest plotfile::

    python scripts/yt_solids.py output.lm --latest
"""

from __future__ import annotations

import argparse
from collections import deque
from dataclasses import dataclass
import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import warnings

warnings.filterwarnings("ignore", message="Unable to import Axes3D")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm, to_rgba
from matplotlib.collections import LineCollection
import numpy as np
import yt

# Loading through yt.load() probes unrelated HDF5 frontends first.  Importing
# the AMReX frontend directly is both faster and avoids an unusable h5py install
# from preventing otherwise valid BoxLib plotfiles from loading.
try:
    import yt.frontends.amrex.io as _amrex_io  # noqa: F401
    from yt.frontends.amrex.data_structures import AMReXDataset
except ImportError:  # yt < 4.4
    AMReXDataset = None


FIELD_TYPE = "boxlib"
DEFAULT_STRESS_PREFIXES = (
    "cell_weighted_solid_deviatoric_cauchy_stress",
    "cell_total_cauchy_stress",
    "cauchy_stress",
    "elastic_stress",
    "stress",
)


@dataclass
class FrameData:
    path: Path
    step: int
    time: float
    extent: tuple[float, float, float, float]
    x: np.ndarray
    y: np.ndarray
    velocity_x: np.ndarray
    velocity_y: np.ndarray
    cell_dx: np.ndarray
    cell_dy: np.ndarray
    eta: np.ndarray
    xix: np.ndarray
    xiy: np.ndarray
    stress_xx: np.ndarray
    stress_xy: np.ndarray
    stress_yx: np.ndarray
    stress_yy: np.ndarray

    @property
    def dx(self) -> float:
        return (self.extent[1] - self.extent[0]) / self.velocity_x.shape[1]

    @property
    def dy(self) -> float:
        return (self.extent[3] - self.extent[2]) / self.velocity_x.shape[0]

    @property
    def vorticity(self) -> np.ndarray:
        # A raw gradient of the finest-resolution image differentiates the
        # repeated pixels representing coarse AMR cells and produces stripes.
        # Use each pixel's native AMR cell width as the differencing distance.
        return amr_derivative(
            self.velocity_y, self.cell_dx, self.dx, axis=1
        ) - amr_derivative(self.velocity_x, self.cell_dy, self.dy, axis=0)

    @property
    def von_mises(self) -> np.ndarray:
        # Symmetrize shear in case the two stored off-diagonal values differ by
        # roundoff.  This is the standard plane-stress/2-D von Mises invariant.
        shear = 0.5 * (self.stress_xy + self.stress_yx)
        vm_squared = (
            self.stress_xx**2
            - self.stress_xx * self.stress_yy
            + self.stress_yy**2
            + 3.0 * shear**2
        )
        return np.sqrt(np.maximum(vm_squared, 0.0))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot vorticity, eta=0.5, masked xix/xiy isolines, advected "
            "tracers, and masked von Mises stress from 2-D BoxLib outputs."
        )
    )
    parser.add_argument(
        "plotfile_root",
        nargs="?",
        default="output.lm",
        type=Path,
        help="directory containing *cell plotfiles (default: output.lm)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("yt_plots"),
        help="PNG output directory (default: yt_plots)",
    )
    parser.add_argument("--start", type=int, help="minimum plotfile step")
    parser.add_argument("--stop", type=int, help="maximum plotfile step, inclusive")
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="render every Nth discovered plotfile (default: 1)",
    )
    parser.add_argument(
        "--latest", action="store_true", help="render only the newest plotfile"
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="do not rerender PNGs that already exist in --output-dir",
    )
    parser.add_argument(
        "--resolution",
        type=int,
        default=0,
        metavar="N",
        help=(
            "pixels along the longer domain direction; 0 uses the finest AMR "
            "cell resolution (default: 0)"
        ),
    )
    parser.add_argument(
        "--eta-threshold",
        type=float,
        default=0.5,
        help="solid mask and interface value (default: 0.5)",
    )
    parser.add_argument(
        "--xix-levels",
        type=int,
        default=6,
        help="number of evenly spaced xix levels (default: 6)",
    )
    parser.add_argument(
        "--xiy-levels",
        type=int,
        default=10,
        help="number of evenly spaced xiy levels (default: 10)",
    )
    parser.add_argument(
        "--particles",
        type=int,
        default=400,
        help="number of massless visualization tracers; 0 disables (default: 400)",
    )
    parser.add_argument(
        "--particle-color",
        default="white",
        help="Matplotlib color for tracer points and trails (default: white)",
    )
    parser.add_argument(
        "--trail-length",
        type=int,
        default=12,
        help="number of rendered positions retained in tracer trails (default: 12)",
    )
    parser.add_argument(
        "--tracer-cfl",
        type=float,
        default=0.75,
        help="maximum tracer displacement in grid cells per substep (default: 0.75)",
    )
    parser.add_argument(
        "--seed", type=int, default=1947, help="tracer jitter seed (default: 1947)"
    )
    parser.add_argument(
        "--scale-samples",
        type=int,
        default=12,
        help=(
            "number of time samples used for stable global color limits; 0 "
            "rescales every frame (default: 12)"
        ),
    )
    parser.add_argument(
        "--percentile",
        type=float,
        default=99.5,
        help="robust percentile used for automatic color maxima (default: 99.5)",
    )
    parser.add_argument(
        "--vorticity-max",
        type=float,
        help="fixed symmetric vorticity limit; overrides automatic scaling",
    )
    parser.add_argument(
        "--stress-max",
        type=float,
        help="fixed von Mises upper limit; overrides automatic scaling",
    )
    parser.add_argument(
        "--vorticity-cmap",
        default="jet",
        help="Matplotlib colormap for vorticity (default: jet)",
    )
    parser.add_argument(
        "--stress-cmap",
        default="viridis",
        help="Matplotlib colormap for von Mises stress (default: viridis)",
    )
    parser.add_argument(
        "--velocity-fields",
        nargs=2,
        metavar=("VELX", "VELY"),
        default=("velocityx", "velocityy"),
        help="on-disk x/y velocity field names",
    )
    parser.add_argument(
        "--stress-prefix",
        help=(
            "stress field prefix before _xx/_xy/_yx/_yy; by default the "
            "cell-weighted solid deviatoric stress is preferred"
        ),
    )
    parser.add_argument(
        "--movie",
        type=Path,
        help="optional MP4 path (a bare filename is placed in --output-dir)",
    )
    parser.add_argument("--fps", type=float, default=20.0, help="movie FPS (default: 20)")
    parser.add_argument(
        "--dpi", type=int, default=160, help="PNG resolution in dots per inch (default: 160)"
    )
    args = parser.parse_args()

    if args.stride < 1:
        parser.error("--stride must be at least 1")
    if args.resolution < 0:
        parser.error("--resolution cannot be negative")
    if args.xix_levels < 1:
        parser.error("--xix-levels must be at least 1")
    if args.xiy_levels < 1:
        parser.error("--xiy-levels must be at least 1")
    if args.particles < 0:
        parser.error("--particles cannot be negative")
    try:
        to_rgba(args.particle_color)
    except ValueError as error:
        parser.error(f"invalid --particle-color: {error}")
    if args.trail_length < 1:
        parser.error("--trail-length must be at least 1")
    if args.tracer_cfl <= 0.0:
        parser.error("--tracer-cfl must be positive")
    if args.scale_samples < 0:
        parser.error("--scale-samples cannot be negative")
    if not 0.0 < args.percentile <= 100.0:
        parser.error("--percentile must be in (0, 100]")
    if args.vorticity_max is not None and args.vorticity_max <= 0.0:
        parser.error("--vorticity-max must be positive")
    if args.stress_max is not None and args.stress_max <= 0.0:
        parser.error("--stress-max must be positive")
    if args.fps <= 0.0:
        parser.error("--fps must be positive")
    return args


def plotfile_step(path: Path) -> int:
    match = re.match(r"^(\d+)cell$", path.name)
    if match is None:
        raise ValueError(f"not a recognized cell plotfile name: {path.name}")
    return int(match.group(1))


def discover_plotfiles(args: argparse.Namespace) -> list[Path]:
    root = args.plotfile_root
    if not root.is_dir():
        raise FileNotFoundError(f"plotfile root does not exist: {root}")

    paths = []
    for path in root.iterdir():
        if not path.is_dir() or re.match(r"^\d+cell$", path.name) is None:
            continue
        # Cell_H is written after the level data and is a better completeness
        # check than the top-level Header alone for a live simulation.
        if (path / "Header").is_file() and (path / "Level_0" / "Cell_H").is_file():
            paths.append(path)
    paths.sort(key=plotfile_step)

    if args.start is not None:
        paths = [path for path in paths if plotfile_step(path) >= args.start]
    if args.stop is not None:
        paths = [path for path in paths if plotfile_step(path) <= args.stop]
    paths = paths[:: args.stride]
    if args.latest and paths:
        paths = paths[-1:]
    if not paths:
        raise FileNotFoundError(f"no complete *cell plotfiles selected under {root}")
    return paths


def load_dataset(path: Path):
    if AMReXDataset is not None:
        return AMReXDataset(str(path))
    return yt.load(str(path))


def choose_resolution(ds, requested: int) -> tuple[int, int]:
    left = np.asarray(ds.domain_left_edge, dtype=float)
    right = np.asarray(ds.domain_right_edge, dtype=float)
    width, height = right[0] - left[0], right[1] - left[1]
    if requested:
        if width >= height:
            return requested, max(2, round(requested * height / width))
        return max(2, round(requested * width / height)), requested

    refine_by = int(getattr(ds, "refine_by", 2))
    factor = refine_by ** int(ds.max_level)
    return (
        int(ds.domain_dimensions[0]) * factor,
        int(ds.domain_dimensions[1]) * factor,
    )


def resolve_field(ds, name: str) -> tuple[str, str]:
    preferred = (FIELD_TYPE, name)
    if preferred in ds.field_list:
        return preferred
    matches = [field for field in ds.field_list if field[1] == name]
    if matches:
        return matches[0]
    available = ", ".join(field[1] for field in ds.field_list)
    raise KeyError(f"field {name!r} is unavailable; plotfile fields: {available}")


def resolve_stress_fields(ds, prefix: str | None) -> dict[str, tuple[str, str]]:
    prefixes = (prefix,) if prefix else DEFAULT_STRESS_PREFIXES
    for candidate in prefixes:
        fields = {
            component: (FIELD_TYPE, f"{candidate}_{component}")
            for component in ("xx", "xy", "yx", "yy")
        }
        if all(field in ds.field_list for field in fields.values()):
            return fields
    requested = prefix or " or ".join(DEFAULT_STRESS_PREFIXES)
    raise KeyError(f"could not find _xx/_xy/_yx/_yy fields for stress prefix {requested!r}")


def read_frame(
    path: Path,
    resolution: tuple[int, int],
    velocity_names: tuple[str, str],
    stress_prefix: str | None,
) -> FrameData:
    ds = load_dataset(path)
    left = np.asarray(ds.domain_left_edge, dtype=float)
    right = np.asarray(ds.domain_right_edge, dtype=float)
    center = 0.5 * (left + right)
    extent = (float(left[0]), float(right[0]), float(left[1]), float(right[1]))
    nx, ny = resolution

    fields = {
        "velocity_x": resolve_field(ds, velocity_names[0]),
        "velocity_y": resolve_field(ds, velocity_names[1]),
        "eta": resolve_field(ds, "eta"),
        "xix": resolve_field(ds, "xix"),
        "xiy": resolve_field(ds, "xiy"),
    }
    fields.update(
        {f"stress_{key}": value for key, value in resolve_stress_fields(ds, stress_prefix).items()}
    )

    slice_2d = ds.slice(2, float(center[2]))
    frb = slice_2d.to_frb(
        float(right[0] - left[0]),
        resolution,
        center=center,
        height=float(right[1] - left[1]),
        periodic=False,
    )
    arrays = {key: np.asarray(frb[field], dtype=float) for key, field in fields.items()}
    arrays["cell_dx"] = np.asarray(frb[("index", "dx")], dtype=float)
    arrays["cell_dy"] = np.asarray(frb[("index", "dy")], dtype=float)
    dx = (extent[1] - extent[0]) / nx
    dy = (extent[3] - extent[2]) / ny
    x = np.linspace(extent[0] + 0.5 * dx, extent[1] - 0.5 * dx, nx)
    y = np.linspace(extent[2] + 0.5 * dy, extent[3] - 0.5 * dy, ny)
    return FrameData(
        path=path,
        step=plotfile_step(path),
        time=float(ds.current_time),
        extent=extent,
        x=x,
        y=y,
        **arrays,
    )


def amr_derivative(
    values: np.ndarray,
    cell_widths: np.ndarray,
    pixel_width: float,
    axis: int,
) -> np.ndarray:
    """Differentiate an FRB using the native AMR cell width at each pixel.

    yt repeats a coarse cell's value over several pixels when an FRB is made at
    the finest AMR resolution.  A one-pixel difference therefore alternates
    between zero and a jump at every coarse-cell edge.  This routine instead
    reaches one local cell width in each direction before differencing.
    """
    size = values.shape[axis]
    if size < 2:
        return np.zeros_like(values)
    strides = np.rint(cell_widths / pixel_width).astype(int)
    strides = np.clip(strides, 1, size - 1)
    derivative = np.empty_like(values, dtype=float)
    coordinates = np.arange(size)

    for stride in np.unique(strides):
        lower = np.clip(coordinates - stride, 0, size - 1)
        upper = np.clip(coordinates + stride, 0, size - 1)
        lower_values = np.take(values, lower, axis=axis)
        upper_values = np.take(values, upper, axis=axis)
        central = (coordinates >= stride) & (coordinates + stride < size)
        distance = np.where(central, 2.0 * stride * pixel_width, stride * pixel_width)
        reshape = [1] * values.ndim
        reshape[axis] = size
        candidate = (upper_values - lower_values) / distance.reshape(reshape)
        mask = strides == stride
        derivative[mask] = candidate[mask]
    return derivative


def robust_max(values: np.ndarray, percentile: float) -> float:
    finite = np.asarray(values)[np.isfinite(values)]
    if finite.size == 0:
        return 1.0
    value = float(np.percentile(finite, percentile))
    if value <= np.finfo(float).eps:
        value = float(np.max(finite))
    return value if value > np.finfo(float).eps else 1.0


def estimate_color_limits(
    paths: list[Path],
    resolution: tuple[int, int],
    args: argparse.Namespace,
) -> tuple[float | None, float | None]:
    vorticity_max = args.vorticity_max
    stress_max = args.stress_max
    if args.scale_samples == 0 or (vorticity_max is not None and stress_max is not None):
        return vorticity_max, stress_max

    count = min(args.scale_samples, len(paths))
    indices = np.unique(np.linspace(0, len(paths) - 1, count, dtype=int))
    vorticity_samples = []
    stress_samples = []
    print(f"Estimating global color limits from {len(indices)} time sample(s) ...")
    for index in indices:
        frame = read_frame(
            paths[index], resolution, tuple(args.velocity_fields), args.stress_prefix
        )
        if vorticity_max is None:
            vorticity_samples.append(np.abs(frame.vorticity).ravel())
        if stress_max is None:
            inside = frame.eta > args.eta_threshold
            if np.any(inside):
                stress_samples.append(frame.von_mises[inside])

    if vorticity_max is None:
        vorticity_max = robust_max(np.concatenate(vorticity_samples), args.percentile)
    if stress_max is None:
        stress_max = (
            robust_max(np.concatenate(stress_samples), args.percentile)
            if stress_samples
            else 1.0
        )
    print(f"Color limits: |vorticity| <= {vorticity_max:.6g}, von Mises <= {stress_max:.6g}")
    return vorticity_max, stress_max


def seed_particles(
    count: int, extent: tuple[float, float, float, float], seed: int
) -> np.ndarray:
    if count == 0:
        return np.empty((0, 2), dtype=float)
    xlo, xhi, ylo, yhi = extent
    aspect = (xhi - xlo) / (yhi - ylo)
    columns = max(1, math.ceil(math.sqrt(count * aspect)))
    rows = max(1, math.ceil(count / columns))
    xx, yy = np.meshgrid(
        (np.arange(columns) + 0.5) / columns,
        (np.arange(rows) + 0.5) / rows,
    )
    particles = np.column_stack((xx.ravel(), yy.ravel()))[:count]
    rng = np.random.default_rng(seed)
    jitter = rng.uniform(-0.28, 0.28, size=particles.shape)
    particles[:, 0] = np.clip(particles[:, 0] + jitter[:, 0] / columns, 0.0, 1.0)
    particles[:, 1] = np.clip(particles[:, 1] + jitter[:, 1] / rows, 0.0, 1.0)
    particles[:, 0] = xlo + particles[:, 0] * (xhi - xlo)
    particles[:, 1] = ylo + particles[:, 1] * (yhi - ylo)
    return particles


def bilinear_sample(field: np.ndarray, points: np.ndarray, frame: FrameData) -> np.ndarray:
    if points.size == 0:
        return np.empty(0, dtype=float)
    ny, nx = field.shape
    qx = np.clip((points[:, 0] - frame.extent[0]) / frame.dx - 0.5, 0.0, nx - 1.0)
    qy = np.clip((points[:, 1] - frame.extent[2]) / frame.dy - 0.5, 0.0, ny - 1.0)
    i0 = np.floor(qx).astype(int)
    j0 = np.floor(qy).astype(int)
    i1 = np.minimum(i0 + 1, nx - 1)
    j1 = np.minimum(j0 + 1, ny - 1)
    fx = qx - i0
    fy = qy - j0
    return (
        (1.0 - fx) * (1.0 - fy) * field[j0, i0]
        + fx * (1.0 - fy) * field[j0, i1]
        + (1.0 - fx) * fy * field[j1, i0]
        + fx * fy * field[j1, i1]
    )


def velocity_at(points: np.ndarray, frame: FrameData) -> np.ndarray:
    return np.column_stack(
        (
            bilinear_sample(frame.velocity_x, points, frame),
            bilinear_sample(frame.velocity_y, points, frame),
        )
    )


def advect_particles(
    particles: np.ndarray, frame: FrameData, dt: float, tracer_cfl: float
) -> np.ndarray:
    if particles.size == 0 or dt <= 0.0:
        return particles.copy()
    max_speed = float(
        np.nanmax(np.hypot(frame.velocity_x, frame.velocity_y), initial=0.0)
    )
    cell_size = min(frame.dx, frame.dy)
    substeps = max(1, math.ceil(dt * max_speed / (tracer_cfl * cell_size)))
    h = dt / substeps
    result = particles.copy()
    lower = np.array([frame.extent[0], frame.extent[2]])
    upper = np.array([frame.extent[1], frame.extent[3]])
    for _ in range(substeps):
        velocity_0 = velocity_at(result, frame)
        midpoint = np.clip(result + 0.5 * h * velocity_0, lower, upper)
        result = np.clip(result + h * velocity_at(midpoint, frame), lower, upper)
    return result


def visible_levels(field: np.ma.MaskedArray, levels: np.ndarray) -> np.ndarray:
    compressed = field.compressed()
    if compressed.size == 0:
        return np.empty(0)
    low, high = float(np.nanmin(compressed)), float(np.nanmax(compressed))
    return levels[(levels >= low) & (levels <= high)]


def solid_reference_levels(
    frame: FrameData, eta_threshold: float, count: int
) -> tuple[np.ndarray, np.ndarray]:
    """Choose fixed, evenly spaced reference-map levels inside the solid."""
    inside = frame.eta > eta_threshold
    if not np.any(inside):
        return np.empty(0), np.empty(0)

    def levels(values: np.ndarray) -> np.ndarray:
        finite = values[inside & np.isfinite(values)]
        if finite.size == 0 or float(np.max(finite)) <= float(np.min(finite)):
            return np.empty(0)
        return np.linspace(float(np.min(finite)), float(np.max(finite)), count + 2)[1:-1]

    return levels(frame.xix), levels(frame.xiy)


def add_solid_overlays(
    ax,
    frame: FrameData,
    eta_threshold: float,
    xix_levels: np.ndarray,
    xiy_levels: np.ndarray,
) -> None:
    ax.contour(
        frame.x,
        frame.y,
        frame.eta,
        levels=[eta_threshold],
        colors="black",
        linewidths=1.35,
        zorder=6,
    )
    outside = frame.eta <= eta_threshold
    masked_xix = np.ma.masked_where(outside, frame.xix)
    masked_xiy = np.ma.masked_where(outside, frame.xiy)
    shown_xix = visible_levels(masked_xix, xix_levels)
    shown_xiy = visible_levels(masked_xiy, xiy_levels)
    if shown_xix.size:
        ax.contour(
            frame.x,
            frame.y,
            masked_xix,
            levels=shown_xix,
            colors="#666666",
            linewidths=0.8,
            linestyles="solid",
            corner_mask=False,
            zorder=5,
        )
    if shown_xiy.size:
        ax.contour(
            frame.x,
            frame.y,
            masked_xiy,
            levels=shown_xiy,
            colors="#888888",
            linewidths=0.8,
            linestyles="dashed",
            corner_mask=False,
            zorder=5,
        )


def add_particle_trails(
    ax,
    history: deque[np.ndarray],
    frame: FrameData,
    eta_threshold: float,
    particle_color: str,
) -> None:
    if not history or history[-1].size == 0:
        return
    if len(history) > 1:
        snapshots = list(history)
        segments = []
        colors = []
        red, green, blue, base_alpha = to_rgba(particle_color)
        for age, (start, end) in enumerate(zip(snapshots[:-1], snapshots[1:]), start=1):
            segment_array = np.stack((start, end), axis=1)
            midpoint = 0.5 * (start + end)
            fluid = (
                (bilinear_sample(frame.eta, start, frame) <= eta_threshold)
                & (bilinear_sample(frame.eta, midpoint, frame) <= eta_threshold)
                & (bilinear_sample(frame.eta, end, frame) <= eta_threshold)
            )
            segment_array = segment_array[fluid]
            segments.extend(segment_array)
            alpha = 0.08 + 0.48 * age / (len(snapshots) - 1)
            colors.extend(
                [(red, green, blue, base_alpha * alpha)] * len(segment_array)
            )
        if segments:
            ax.add_collection(
                LineCollection(segments, colors=colors, linewidths=0.45, zorder=7)
            )
    particles = history[-1]
    particles = particles[
        bilinear_sample(frame.eta, particles, frame) <= eta_threshold
    ]
    if particles.size == 0:
        return
    ax.scatter(
        particles[:, 0],
        particles[:, 1],
        s=4.0,
        c=particle_color,
        edgecolors="black",
        linewidths=0.18,
        alpha=0.9,
        zorder=8,
    )


def render_frame(
    frame: FrameData,
    history: deque[np.ndarray],
    output_path: Path,
    args: argparse.Namespace,
    vorticity_max: float | None,
    stress_max: float | None,
    xix_levels: np.ndarray,
    xiy_levels: np.ndarray,
) -> None:
    vorticity = frame.vorticity
    von_mises = frame.von_mises
    if vorticity_max is None:
        frame_vorticity_max = robust_max(np.abs(vorticity), args.percentile)
    else:
        frame_vorticity_max = vorticity_max
    solid = frame.eta > args.eta_threshold
    if stress_max is None:
        frame_stress_max = (
            robust_max(von_mises[solid], args.percentile) if np.any(solid) else 1.0
        )
    else:
        frame_stress_max = stress_max

    fig, ax = plt.subplots(figsize=(7.7, 6.25), constrained_layout=True)

    vort_image = ax.imshow(
        vorticity,
        origin="lower",
        extent=frame.extent,
        cmap=args.vorticity_cmap,
        norm=TwoSlopeNorm(
            vmin=-frame_vorticity_max, vcenter=0.0, vmax=frame_vorticity_max
        ),
        interpolation="bicubic",
    )

    stress_cmap = matplotlib.colormaps[args.stress_cmap].copy()
    stress_cmap.set_bad((0.0, 0.0, 0.0, 0.0))
    masked_stress = np.ma.masked_where(~solid, von_mises)
    stress_image = ax.imshow(
        masked_stress,
        origin="lower",
        extent=frame.extent,
        cmap=stress_cmap,
        norm=Normalize(vmin=0.0, vmax=frame_stress_max),
        interpolation="bicubic",
    )

    vorticity_bar = fig.colorbar(
        vort_image,
        ax=ax,
        location="left",
        shrink=0.86,
        pad=0.08,
        label="vorticity [code units]",
    )
    vorticity_bar.ax.yaxis.set_label_position("left")
    fig.colorbar(
        stress_image,
        ax=ax,
        location="right",
        shrink=0.86,
        pad=0.04,
        label="von Mises stress [code units]",
    )

    add_solid_overlays(
        ax, frame, args.eta_threshold, xix_levels=xix_levels, xiy_levels=xiy_levels
    )
    ax.set_xlim(frame.extent[0], frame.extent[1])
    ax.set_ylim(frame.extent[2], frame.extent[3])
    ax.set_aspect("equal")
    ax.set_xlabel("x [code length]")
    ax.set_ylabel("y [code length]")
    ax.set_title(
        r"Vorticity $\omega_z=\partial_xv_y-\partial_yv_x$; "
        + fr"von Mises stress for $\eta>{args.eta_threshold:g}$"
    )

    add_particle_trails(
        ax, history, frame, args.eta_threshold, args.particle_color
    )
    fig.suptitle(
        f"{frame.path.parent.name}/{frame.path.name}    step={frame.step}    t={frame.time:.6g}",
        fontsize=11,
    )
    fig.savefig(output_path, dpi=args.dpi)
    plt.close(fig)


def make_movie(frame_paths: list[Path], movie: Path, fps: float) -> None:
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg is None:
        raise RuntimeError("--movie requires ffmpeg, but ffmpeg is not on PATH")
    movie.parent.mkdir(parents=True, exist_ok=True)
    frame_duration = 1.0 / fps
    manifest_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".ffconcat", dir=frame_paths[0].parent, delete=False
        ) as manifest:
            manifest_path = Path(manifest.name)
            manifest.write("ffconcat version 1.0\n")
            for frame_path in frame_paths:
                escaped = frame_path.resolve().as_posix().replace("'", "'\\''")
                manifest.write(f"file '{escaped}'\n")
                manifest.write(f"duration {frame_duration:.12g}\n")
            escaped = frame_paths[-1].resolve().as_posix().replace("'", "'\\''")
            manifest.write(f"file '{escaped}'\n")
        subprocess.run(
            [
                ffmpeg,
                "-y",
                "-loglevel",
                "error",
                "-f",
                "concat",
                "-safe",
                "0",
                "-i",
                str(manifest_path),
                "-vf",
                "pad=ceil(iw/2)*2:ceil(ih/2)*2,format=yuv420p",
                "-movflags",
                "+faststart",
                str(movie),
            ],
            check=True,
        )
    finally:
        if manifest_path is not None:
            manifest_path.unlink(missing_ok=True)


def main() -> None:
    args = parse_args()
    yt.set_log_level(40)
    paths = discover_plotfiles(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_paths = [
        args.output_dir / f"frame_{index:06d}_step_{plotfile_step(path):09d}.png"
        for index, path in enumerate(paths)
    ]

    if args.skip_existing and all(path.is_file() for path in output_paths):
        print(f"All {len(output_paths)} selected frame(s) already exist; nothing to render")
        if args.movie is not None:
            movie = args.movie
            if not movie.is_absolute() and movie.parent == Path("."):
                movie = args.output_dir / movie
            print(f"Writing movie -> {movie}")
            make_movie(output_paths, movie, args.fps)
        return

    first_ds = load_dataset(paths[0])
    resolution = choose_resolution(first_ds, args.resolution)
    print(
        f"Selected {len(paths)} plotfile(s), steps {plotfile_step(paths[0])} "
        f"through {plotfile_step(paths[-1])}; FRB resolution {resolution[0]}x{resolution[1]}"
    )
    vorticity_max, stress_max = estimate_color_limits(paths, resolution, args)

    first_frame = read_frame(
        paths[0], resolution, tuple(args.velocity_fields), args.stress_prefix
    )
    xix_levels, _ = solid_reference_levels(
        first_frame, args.eta_threshold, args.xix_levels
    )
    _, xiy_levels = solid_reference_levels(
        first_frame, args.eta_threshold, args.xiy_levels
    )
    particles = seed_particles(args.particles, first_frame.extent, args.seed)
    history: deque[np.ndarray] = deque(maxlen=args.trail_length)
    history.append(particles.copy())
    rendered_paths = []

    for index, path in enumerate(paths):
        frame = (
            first_frame
            if index == 0
            else read_frame(path, resolution, tuple(args.velocity_fields), args.stress_prefix)
        )
        output_path = output_paths[index]
        progress = f"[{index + 1:>{len(str(len(paths)))}}/{len(paths)}]"
        if args.skip_existing and output_path.is_file():
            print(f"{progress} skipping existing {output_path}")
        else:
            print(f"{progress} {path} -> {output_path}")
            render_frame(
                frame,
                history,
                output_path,
                args,
                vorticity_max=vorticity_max,
                stress_max=stress_max,
                xix_levels=xix_levels,
                xiy_levels=xiy_levels,
            )
        rendered_paths.append(output_path)

        if index + 1 < len(paths):
            next_ds = load_dataset(paths[index + 1])
            dt = float(next_ds.current_time) - frame.time
            particles = advect_particles(particles, frame, dt, args.tracer_cfl)
            history.append(particles.copy())

    if args.movie is not None:
        movie = args.movie
        if not movie.is_absolute() and movie.parent == Path("."):
            movie = args.output_dir / movie
        print(f"Writing movie -> {movie}")
        make_movie(rendered_paths, movie, args.fps)

    print(f"Done: wrote {len(rendered_paths)} frame(s) to {args.output_dir}")


if __name__ == "__main__":
    main()

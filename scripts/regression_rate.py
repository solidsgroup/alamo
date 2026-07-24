#!/usr/bin/env python3
"""Compute solid regression rate from LowMach HDF5 (Chombo/AMReX) plotfiles.

Reads every ``*cell.h5`` plotfile in a directory, reconstructs a
finest-resolution profile of a phase field (``rigid_eta`` by default) merged
across AMR levels, locates the y-position where it crosses a threshold
(0.5 by default), and reports the regression rate as the time derivative of
that front position.

Each HDF5 plotfile stores, per AMR level, a set of rectangular boxes. Every
box's data is a flat array laid out component-major, and within a component,
Fortran order (i fastest) over the box's cells -- this matches the AMReX FAB
layout used when the HDF5 plotfile writer flattens each box. Levels finer
than 0 only cover the (adaptively refined) sub-region tracking the interface,
so coarser-level data is upsampled by repetition to fill gaps outside the
finer level's boxes.

Examples
--------
Regression rate of the AP surface (rigid_eta = 0.5 contour) over the whole run::

    python scripts/regression_rate.py output.lm.ap_monopropellant_fullfeedback

Track a different field/threshold and dump the raw front-position history::

    python scripts/regression_rate.py output.lm.ap_htpb_fullfeedback \\
        --field rigid_eta --threshold 0.5 --csv front_position.csv
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

PLOTFILE_RE = re.compile(r"^(\d+)cell\.h5$")


@dataclass
class Snapshot:
    step: int
    time: float
    front_y: float


def discover_plotfiles(root: Path) -> list[Path]:
    if not root.is_dir():
        raise FileNotFoundError(f"plotfile directory does not exist: {root}")
    paths = [p for p in root.iterdir() if PLOTFILE_RE.match(p.name)]
    if not paths:
        raise FileNotFoundError(f"no *cell.h5 plotfiles found under {root}")
    paths.sort(key=lambda p: int(PLOTFILE_RE.match(p.name).group(1)))
    return paths


def component_index(f: h5py.File, field: str) -> int:
    names = {v.decode() if isinstance(v, bytes) else v: int(k.split("_")[1])
             for k, v in f.attrs.items() if k.startswith("component_")}
    if field not in names:
        available = ", ".join(sorted(names))
        raise KeyError(f"field {field!r} not found; available fields: {available}")
    return names[field]


def level_shape(f: h5py.File, level: int) -> tuple[int, int]:
    prob_domain = f[f"level_{level}"].attrs["prob_domain"]
    lo_i, lo_j, hi_i, hi_j = (int(x) for x in prob_domain)
    return hi_i - lo_i + 1, hi_j - lo_j + 1


def level_field(f: h5py.File, level: int, comp: int, ncomp: int) -> np.ndarray:
    """Return this level's field as a (ni, nj) array, NaN where not covered."""
    grp = f[f"level_{level}"]
    ni, nj = level_shape(f, level)
    grid = np.full((ni, nj), np.nan)
    boxes = grp["boxes"][:]
    offsets = grp["data:offsets=0"][:]
    data = grp["data:datatype=0"][:]
    for b, (lo_i, lo_j, hi_i, hi_j) in enumerate(boxes):
        bni, bnj = hi_i - lo_i + 1, hi_j - lo_j + 1
        ncells = bni * bnj
        chunk = data[offsets[b]:offsets[b + 1]]
        comp_chunk = chunk[comp * ncells:(comp + 1) * ncells]
        grid[lo_i:hi_i + 1, lo_j:hi_j + 1] = comp_chunk.reshape((bni, bnj), order="F")
    return grid


def merged_field(f: h5py.File, comp: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Merge all AMR levels onto the finest level's resolution.

    Returns (x, y, field) where field has shape (nx_finest, ny_finest).
    """
    finest_level = int(f.attrs["finest_level"][0])
    ncomp = int(f.attrs["num_components"][0])
    ref_ratio = {lvl: int(f[f"level_{lvl}"].attrs["ref_ratio"][0])
                 for lvl in range(finest_level)}

    def upsample_factor(level: int) -> int:
        factor = 1
        for lvl in range(level, finest_level):
            factor *= ref_ratio[lvl]
        return factor

    ni_f, nj_f = level_shape(f, finest_level)
    merged = np.full((ni_f, nj_f), np.nan)
    for level in range(finest_level + 1):
        grid = level_field(f, level, comp, ncomp)
        factor = upsample_factor(level)
        if factor > 1:
            grid = np.repeat(np.repeat(grid, factor, axis=0), factor, axis=1)
        valid = ~np.isnan(grid)
        merged[valid] = grid[valid]

    dx = float(f[f"level_{finest_level}"].attrs["dx"][0])
    prob_lo = f[f"level_{finest_level}"].attrs["prob_lo"]
    x = prob_lo[0] + (np.arange(ni_f) + 0.5) * dx
    y = prob_lo[1] + (np.arange(nj_f) + 0.5) * dx
    return x, y, merged


def front_position(y: np.ndarray, field_1d: np.ndarray, threshold: float) -> float | None:
    """Linearly interpolate the y where field_1d crosses threshold.

    If multiple crossings exist (e.g. transient noise), the one nearest the
    domain center is returned, since the tracked interface starts there.
    """
    values = field_1d - threshold
    sign_changes = np.where(np.diff(np.sign(values)) != 0)[0]
    if sign_changes.size == 0:
        return None
    if sign_changes.size > 1:
        center = 0.5 * (y[0] + y[-1])
        crossing_y = y[sign_changes] + 0.5 * (y[1] - y[0])
        sign_changes = sign_changes[[np.argmin(np.abs(crossing_y - center))]]
    i0 = sign_changes[0]
    i1 = i0 + 1
    frac = values[i0] / (values[i0] - values[i1])
    return float(y[i0] + frac * (y[i1] - y[i0]))


def read_snapshot(path: Path, field: str, threshold: float) -> Snapshot:
    with h5py.File(path, "r") as f:
        time = float(f.attrs["time"][0])
        comp = component_index(f, field)
        _, y, grid = merged_field(f, comp)
        # Average the front position across all x-columns; the geometry here
        # is periodic/uniform in x, so this also smooths out per-column noise.
        crossings = [front_position(y, grid[i, :], threshold) for i in range(grid.shape[0])]
        crossings = [c for c in crossings if c is not None]
        if not crossings:
            raise RuntimeError(f"no {field}={threshold} crossing found in {path}")
        front_y = float(np.mean(crossings))
    step = int(PLOTFILE_RE.match(path.name).group(1))
    return Snapshot(step=step, time=time, front_y=front_y)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plotfile_root", type=Path,
                         help="directory containing *cell.h5 plotfiles")
    parser.add_argument("--field", default="rigid_eta",
                         help="field to track (default: rigid_eta)")
    parser.add_argument("--threshold", type=float, default=0.5,
                         help="contour value defining the front (default: 0.5)")
    parser.add_argument("--csv", type=Path,
                         help="optional path to write step,time,front_y,rate as CSV")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = discover_plotfiles(args.plotfile_root)
    print(f"Found {len(paths)} plotfile(s) under {args.plotfile_root}")

    snapshots = [read_snapshot(p, args.field, args.threshold) for p in paths]

    rows = []
    print(f"{'step':>10} {'time [s]':>14} {'front_y [m]':>14} {'rate [m/s]':>14}")
    prev = None
    for snap in snapshots:
        rate = None
        if prev is not None and snap.time > prev.time:
            rate = (snap.front_y - prev.front_y) / (snap.time - prev.time)
        rate_str = f"{rate:14.6e}" if rate is not None else " " * 14
        print(f"{snap.step:10d} {snap.time:14.6e} {snap.front_y:14.6e} {rate_str}")
        rows.append((snap.step, snap.time, snap.front_y, rate))
        prev = snap

    times = np.array([s.time for s in snapshots])
    fronts = np.array([s.front_y for s in snapshots])
    if len(snapshots) >= 2:
        slope, intercept = np.polyfit(times, fronts, 1)
        print(f"\nLinear fit over full run: d(front_y)/dt = {slope:.6e} m/s "
              f"=> regression rate magnitude = {abs(slope):.6e} m/s "
              f"({abs(slope) * 1000.0:.6g} mm/s)")

    if args.csv is not None:
        with open(args.csv, "w") as fh:
            fh.write("step,time,front_y,rate\n")
            for step, time, front_y, rate in rows:
                fh.write(f"{step},{time},{front_y},{'' if rate is None else rate}\n")
        print(f"Wrote {args.csv}")


if __name__ == "__main__":
    main()

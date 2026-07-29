#!/usr/bin/env python3
"""Compute solid regression rate from LowMach AMReX plotfiles.

Reads every ``*cell`` plotfile directory in a directory (plain AMReX/BoxLib
plotfiles -- this build does not write the HDF5 plotfile variant, so plotfiles
are read via yt rather than h5py), extracts a finest-resolution profile of a
phase field (``rigid_eta`` by default), locates the y-position where it
crosses a threshold (0.5 by default), and reports the regression rate as the
time derivative of that front position.

Each plotfile is loaded with ``yt.load`` and resampled onto a uniform grid at
the finest AMR level present via ``ds.covering_grid``, which yt fills from
whichever level actually covers each cell -- equivalent to the coarse-to-fine
merge this script used to do by hand against the HDF5 writer's per-level boxes.

Examples
--------
Regression rate of the AP surface (rigid_eta = 0.5 contour) over the whole run::

    python scripts/regression_rate.py output.lm.ap_monopropellant

Track a different field/threshold and dump the raw front-position history::

    python scripts/regression_rate.py output.lm.ap_monopropellant \\
        --field rigid_eta --threshold 0.5 --csv front_position.csv
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yt

yt.set_log_level(50)  # suppress per-plotfile INFO spam; errors still surface

PLOTFILE_RE = re.compile(r"^(\d+)cell$")


@dataclass
class Snapshot:
    step: int
    time: float
    front_y: float


def discover_plotfiles(root: Path) -> list[Path]:
    if not root.is_dir():
        raise FileNotFoundError(f"plotfile directory does not exist: {root}")
    paths = [p for p in root.iterdir() if p.is_dir() and PLOTFILE_RE.match(p.name)]
    if not paths:
        raise FileNotFoundError(f"no *cell plotfiles found under {root}")
    paths.sort(key=lambda p: int(PLOTFILE_RE.match(p.name).group(1)))
    return paths


def merged_field(ds, field: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Resample ``field`` onto the finest AMR level's resolution.

    Returns (x, y, field) where field has shape (nx_finest, ny_finest); the
    domain here is 2D (or a thin 3D slab), so the trailing/third axis is
    dropped.
    """
    finest_level = ds.index.max_level
    dims = ds.domain_dimensions * ds.refine_by ** finest_level
    dims = np.array([dims[0], dims[1], 1])
    cg = ds.covering_grid(level=finest_level, left_edge=ds.domain_left_edge,
                          dims=dims)
    data = cg["boxlib", field].to_ndarray()[:, :, 0]

    dx = (ds.domain_width / dims)[:2].to("m").ndarray_view()
    lo = ds.domain_left_edge[:2].to("m").ndarray_view()
    x = lo[0] + (np.arange(dims[0]) + 0.5) * dx[0]
    y = lo[1] + (np.arange(dims[1]) + 0.5) * dx[1]
    return x, y, data


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


def read_snapshot(path: Path, field: str, threshold: float) -> Snapshot | None:
    """Read one plotfile's front position, or None if the interface is gone.

    A fast/high-pressure run can fully consume the solid before ``stop_time``
    -- later plotfiles then have no ``field=threshold`` crossing anywhere in
    the domain. That is not a data error, just the run continuing past
    burnout, so it is reported as None (skip this snapshot) rather than
    raised as an exception.
    """
    ds = yt.load(str(path))
    # covering_grid can otherwise refuse to sample flush against a
    # non-periodic domain boundary (it wants ghost cells past the edge);
    # the field being tracked here never needs those ghost cells since
    # front_position only interpolates strictly inside the domain.
    ds.force_periodicity()
    time = float(ds.current_time.to("s").value)
    _, y, grid = merged_field(ds, field)
    # Average the front position across all x-columns; the geometry here
    # is periodic/uniform in x, so this also smooths out per-column noise.
    crossings = [front_position(y, grid[i, :], threshold) for i in range(grid.shape[0])]
    crossings = [c for c in crossings if c is not None]
    if not crossings:
        return None
    front_y = float(np.mean(crossings))
    step = int(PLOTFILE_RE.match(path.name).group(1))
    return Snapshot(step=step, time=time, front_y=front_y)


def find_steady_state(
    times: np.ndarray,
    fronts: np.ndarray,
    burnout_frac: float,
    steady_tol: float,
    min_points: int,
) -> tuple[np.ndarray, int, int]:
    """Classify per-interval rates and locate the steady-burning window.

    Returns (phases, steady_start, steady_end) where ``phases`` has one entry
    per interval (length ``len(times) - 1``) labeled "transient", "steady",
    or "extinguished", and ``steady_start``/``steady_end`` are snapshot
    indices (inclusive) bounding the steady-state interval run used for the
    reported regression rate. If no steady window is found, both are -1.
    """
    n = len(times)
    rates = np.diff(fronts) / np.diff(times)
    phases = np.full(n - 1, "transient", dtype=object)
    if n < 3:
        return phases, -1, -1

    max_rate = float(np.max(np.abs(rates)))
    if max_rate <= 0.0:
        phases[:] = "extinguished"
        return phases, -1, -1

    # A cell keeps regressing once ignited; a sustained drop to a small
    # fraction of the peak rate near the end of the run means the AP burned
    # through (or the run stopped advancing) and later data should be
    # dropped, not just a single noisy interval mid-run.
    burning = np.abs(rates) > burnout_frac * max_rate
    burning_indices = np.where(burning)[0]
    active_end = int(burning_indices[-1]) + 1 if burning_indices.size else 0
    phases[active_end:] = "extinguished"
    if active_end < min_points:
        return phases, -1, -1

    active_rates = rates[:active_end]
    median_rate = float(np.median(active_rates))
    steady = np.abs(active_rates - median_rate) <= steady_tol * abs(median_rate)

    # Keep the longest contiguous run of steady intervals; ties favor the
    # latest run, since transients are expected at the start of a burn.
    best_start, best_len = -1, 0
    run_start = None
    for i in range(active_end + 1):
        is_steady = i < active_end and steady[i]
        if is_steady and run_start is None:
            run_start = i
        if not is_steady and run_start is not None:
            run_len = i - run_start
            if run_len >= best_len:
                best_start, best_len = run_start, run_len
            run_start = None
    if run_start is not None:
        run_len = active_end - run_start
        if run_len >= best_len:
            best_start, best_len = run_start, run_len

    phases[:active_end][steady] = "steady"
    if best_len < min_points:
        return phases, -1, -1
    return phases, best_start, best_start + best_len


#: Conversion factor from m/s to the unit named by --unit.
UNIT_SCALE = {"m/s": 1.0, "mm/s": 1.0e3, "cm/s": 1.0e2}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("plotfile_root", type=Path,
                         help="directory containing *cell plotfile directories")
    parser.add_argument("--field", default="rigid_eta",
                         help="field to track (default: rigid_eta)")
    parser.add_argument("--threshold", type=float, default=0.5,
                         help="contour value defining the front (default: 0.5)")
    parser.add_argument("--burnout-frac", type=float, default=0.1,
                         help="interval rate below this fraction of the peak rate, "
                              "sustained through the end of the run, is treated as "
                              "extinguished and excluded (default: 0.1)")
    parser.add_argument("--steady-tol", type=float, default=0.1,
                         help="relative tolerance to the median active rate used to "
                              "classify an interval as steady state (default: 0.1)")
    parser.add_argument("--min-steady-points", type=int, default=3,
                         help="minimum number of intervals required to report a "
                              "steady-state rate (default: 3)")
    parser.add_argument("--min-time", type=float, default=0.0,
                         help="exclude plotfiles before this time [s] from the "
                              "regression-rate calculation entirely (e.g. to skip "
                              "a known startup transient before steady burning is "
                              "established; default: 0.0, no exclusion)")
    parser.add_argument("--csv", type=Path,
                         help="optional path to write step,time,front_y,rate,phase as CSV")
    parser.add_argument("--rate-only", action="store_true",
                         help="print only the steady-state regression-rate magnitude "
                              "(in --unit) to stdout and exit; exit status is nonzero "
                              "and 'nan' is printed if no steady window is found -- "
                              "for use in scripted parameter sweeps")
    parser.add_argument("--unit", choices=sorted(UNIT_SCALE), default="m/s",
                         help="unit for the --rate-only value (default: m/s)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = discover_plotfiles(args.plotfile_root)
    raw_snapshots = [read_snapshot(p, args.field, args.threshold) for p in paths]
    n_missing = sum(1 for s in raw_snapshots if s is None)
    if n_missing:
        print(f"note: {n_missing}/{len(raw_snapshots)} plotfile(s) had no "
              f"{args.field}={args.threshold} crossing (interface fully "
              f"consumed/left the domain); excluded from the fit",
              file=sys.stderr)
    snapshots = [s for s in raw_snapshots if s is not None]
    if len(snapshots) < 2:
        raise RuntimeError(
            f"fewer than 2 usable plotfiles under {args.plotfile_root} "
            f"after excluding missing-crossing snapshots")
    if args.min_time > 0.0:
        snapshots = [s for s in snapshots if s.time >= args.min_time]
        if len(snapshots) < 2:
            raise RuntimeError(
                f"--min-time {args.min_time} leaves fewer than 2 plotfiles "
                f"under {args.plotfile_root}")
    times = np.array([s.time for s in snapshots])
    fronts = np.array([s.front_y for s in snapshots])

    phases, steady_start, steady_end = find_steady_state(
        times, fronts, args.burnout_frac, args.steady_tol, args.min_steady_points)

    rate_magnitude = None
    if steady_start >= 0:
        steady_times = times[steady_start:steady_end + 1]
        steady_fronts = fronts[steady_start:steady_end + 1]
        slope, _ = np.polyfit(steady_times, steady_fronts, 1)
        rate_magnitude = abs(slope)

    if args.rate_only:
        if rate_magnitude is None:
            print("nan")
            raise SystemExit(1)
        print(f"{rate_magnitude * UNIT_SCALE[args.unit]:.6e}")
        return

    print(f"Found {len(paths)} plotfile(s) under {args.plotfile_root}")
    rows = []
    print(f"{'step':>10} {'time [s]':>14} {'front_y [m]':>14} "
          f"{'rate [m/s]':>14} {'phase':>13}")
    for i, snap in enumerate(snapshots):
        rate = None
        phase = ""
        if i > 0:
            rate = (fronts[i] - fronts[i - 1]) / (times[i] - times[i - 1])
            phase = phases[i - 1]
        rate_str = f"{rate:14.6e}" if rate is not None else " " * 14
        print(f"{snap.step:10d} {snap.time:14.6e} {snap.front_y:14.6e} "
              f"{rate_str} {phase:>13}")
        rows.append((snap.step, snap.time, snap.front_y, rate, phase))

    if steady_start < 0:
        print("\nNo sustained steady-burning window found; "
              "no regression rate reported.")
    else:
        print(f"\nSteady-state window: steps {snapshots[steady_start].step}-"
              f"{snapshots[steady_end].step} "
              f"(t = {times[steady_start]:.6g}-{times[steady_end]:.6g} s)")
        print(f"Steady-state regression rate: "
              f"regression rate magnitude = {rate_magnitude:.6e} m/s "
              f"({rate_magnitude * 1000.0:.6g} mm/s)")
        if np.any(phases == "extinguished"):
            first_ext = int(np.argmax(phases == "extinguished"))
            print(f"AP stopped burning after step {snapshots[first_ext].step} "
                  f"(t = {times[first_ext]:.6g} s); later data excluded.")

    if args.csv is not None:
        with open(args.csv, "w") as fh:
            fh.write("step,time,front_y,rate,phase\n")
            for step, time, front_y, rate, phase in rows:
                fh.write(f"{step},{time},{front_y},"
                         f"{'' if rate is None else rate},{phase}\n")
        print(f"Wrote {args.csv}")


if __name__ == "__main__":
    main()

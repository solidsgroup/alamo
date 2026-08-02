#!/usr/bin/env python3
"""Measure a planar condensed/gas interface without lateral averaging."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yt

yt.set_log_level(50)


def _crossing_position(profile, coordinates, level):
    crossings = np.flatnonzero(
        (profile[:-1] >= level) & (profile[1:] < level))
    if not len(crossings):
        return np.nan

    # The first downward crossing is the surface connected to the solid below.
    index = crossings[0]
    difference = profile[index + 1] - profile[index]
    if difference == 0.0:
        return np.nan
    return coordinates[index] + (
        (level - profile[index])
        * (coordinates[index + 1] - coordinates[index])
        / difference)


def read_interface_history(output_directory):
    plotfiles = list(Path(output_directory).glob("*cell"))
    if len(plotfiles) < 3:
        raise RuntimeError(
            f"{output_directory} must contain at least three plotfiles")

    history = []
    for plotfile in plotfiles:
        dataset = yt.load(str(plotfile))
        # Avoid yt rejecting an exact covering grid when floating-point
        # reconstruction places its upper edge a few ulps outside the domain.
        dataset.force_periodicity()
        level = dataset.index.max_level
        refinement = int(dataset.refine_by) ** level
        dimensions = np.asarray(
            dataset.domain_dimensions, dtype=int) * refinement
        grid = dataset.covering_grid(
            level=level,
            left_edge=dataset.domain_left_edge,
            dims=dimensions)

        eta = np.asarray(grid[("boxlib", "rigid_eta")]).squeeze()
        temperature = np.asarray(
            grid[("boxlib", "temperature")]).squeeze()
        if eta.ndim != 2:
            raise RuntimeError(
                f"Expected a two-dimensional eta field in {plotfile}")
        if not np.all(np.isfinite(eta)):
            raise AssertionError(f"Non-finite rigid_eta in {plotfile}")
        if not np.all(np.isfinite(temperature)):
            raise AssertionError(f"Non-finite temperature in {plotfile}")

        nx, ny = eta.shape
        lower = np.asarray(dataset.domain_left_edge, dtype=float)
        upper = np.asarray(dataset.domain_right_edge, dtype=float)
        x = lower[0] + (np.arange(nx) + 0.5) * (
            upper[0] - lower[0]) / nx
        y = lower[1] + (np.arange(ny) + 0.5) * (
            upper[1] - lower[1]) / ny

        positions = np.asarray([
            _crossing_position(eta[i, :], y, 0.5) for i in range(nx)])
        widths = np.asarray([
            _crossing_position(eta[i, :], y, 0.1)
            - _crossing_position(eta[i, :], y, 0.9)
            for i in range(nx)])
        if not np.all(np.isfinite(positions)):
            raise RuntimeError(
                f"Not every x column has a connected eta=0.5 crossing "
                f"in {plotfile}")

        extreme_index = np.argmin(positions)
        history.append({
            "time_s": float(dataset.current_time),
            "extreme_position_m": positions[extreme_index],
            "mean_position_m": np.mean(positions),
            "maximum_position_m": np.max(positions),
            "extreme_x_m": x[extreme_index],
            "median_width_m": np.nanmedian(widths),
        })

    return sorted(history, key=lambda row: row["time_s"])


def fit_interface_rates(history, fit_start_fraction=0.5):
    if not 0.0 <= fit_start_fraction < 1.0:
        raise ValueError("fit_start_fraction must lie in [0, 1)")
    times = np.asarray([row["time_s"] for row in history])
    fit_start = fit_start_fraction * times[-1]
    fit = times >= fit_start
    if np.count_nonzero(fit) < 3:
        raise RuntimeError("The interface-rate fit requires at least three points")

    rates = {}
    for label in ("extreme", "mean", "maximum"):
        positions = np.asarray([
            row[f"{label}_position_m"] for row in history])
        slope, intercept = np.polyfit(times[fit], positions[fit], 1)
        rates[f"{label}_burn_rate_mm_per_s"] = -1000.0 * slope
        rates[f"{label}_fit_intercept_m"] = intercept

    rates["fit_start_s"] = fit_start
    rates["fit_end_s"] = times[-1]
    rates["fit_start_fraction"] = fit_start_fraction
    return rates


def write_results(output_directory, history, rates):
    output_directory = Path(output_directory)
    csv_path = output_directory / "interface-position.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=history[0].keys())
        writer.writeheader()
        writer.writerows(history)

    mean_rate = rates["mean_burn_rate_mm_per_s"]
    summary = {
        "calibration_metric":
            "minimum y of the solid-connected eta=0.5 crossing",
        **rates,
        "extreme_to_mean_rate_ratio":
            rates["extreme_burn_rate_mm_per_s"]
            / mean_rate if mean_rate != 0.0 else None,
    }
    with (output_directory / "interface-rates.json").open("w") as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")

    times = np.asarray([row["time_s"] for row in history])
    extreme = np.asarray([
        row["extreme_position_m"] for row in history])
    mean = np.asarray([row["mean_position_m"] for row in history])
    maximum = np.asarray([
        row["maximum_position_m"] for row in history])
    fit = times >= rates["fit_start_s"]
    fitted_extreme = (
        -rates["extreme_burn_rate_mm_per_s"] / 1000.0 * times[fit]
        + rates["extreme_fit_intercept_m"])

    plt.clf()
    plt.plot(1.0e3 * times, 1.0e6 * extreme, "o-",
             markerfacecolor="none", label="Extreme (calibration)")
    plt.plot(1.0e3 * times, 1.0e6 * mean, "--",
             label="Lateral mean (diagnostic only)")
    plt.plot(1.0e3 * times, 1.0e6 * maximum, ":",
             label="Maximum")
    plt.plot(1.0e3 * times[fit], 1.0e6 * fitted_extreme,
             color="black", linewidth=1.0, label="Extreme fit")
    plt.xlabel("Time [ms]")
    plt.ylabel("eta=0.5 position [um]")
    plt.title(
        "Extreme burn rate: "
        f"{rates['extreme_burn_rate_mm_per_s']:.3f} mm/s")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_directory / "interface-position.png", dpi=160)


def analyze(output_directory, fit_start_fraction=0.5, write=True):
    history = read_interface_history(output_directory)
    rates = fit_interface_rates(history, fit_start_fraction)
    if write:
        write_results(output_directory, history, rates)
    return history, rates


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Fit the deepest solid-connected eta=0.5 interface position."))
    parser.add_argument("output_directory")
    parser.add_argument(
        "--fit-start-fraction", type=float, default=0.5,
        help="fraction of the final time at which the linear fit begins")
    args = parser.parse_args()

    _, rates = analyze(
        args.output_directory,
        fit_start_fraction=args.fit_start_fraction)
    print(
        "Calibration burn rate (extreme eta=0.5 position): "
        f"{rates['extreme_burn_rate_mm_per_s']:.8g} mm/s")
    print(
        "Lateral-mean diagnostic (do not calibrate): "
        f"{rates['mean_burn_rate_mm_per_s']:.8g} mm/s")


if __name__ == "__main__":
    main()

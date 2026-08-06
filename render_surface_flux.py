#!/usr/bin/env python3
"""Render heat and mass flux versus burn-surface arc length."""

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def render(arguments):
    csv_path, output_directory, dpi = arguments
    with csv_path.open() as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        return csv_path.name, 0

    figure, axes = plt.subplots(
        2, 1, figsize=(9.0, 6.5), sharex=True,
        gridspec_kw={"hspace": 0.08})
    contour_ids = sorted({int(row["contour_id"]) for row in rows})
    for contour_id in contour_ids:
        contour = [row for row in rows if int(row["contour_id"]) == contour_id]
        arc = 1.0e3 * np.asarray([
            float(row["arc_length_m"]) for row in contour])
        heat_flux = 1.0e-6 * np.asarray([
            float(row["heat_flux_into_solid_W_m2"]) for row in contour])
        mass_flux = np.asarray([
            float(row["mass_flux_kg_m2_s"]) for row in contour])
        ap = np.asarray([row["species"] == "AP" for row in contour])

        axes[0].plot(arc, heat_flux, color="0.25", linewidth=0.8)
        axes[1].plot(arc, mass_flux, color="0.25", linewidth=0.8)
        for axis, values in zip(axes, (heat_flux, mass_flux)):
            axis.scatter(
                arc[ap], values[ap], s=9.0, color="#c43c39",
                linewidths=0.0, label="AP" if contour_id == contour_ids[0] else None)
            axis.scatter(
                arc[~ap], values[~ap], s=9.0, color="#3569b7",
                linewidths=0.0,
                label="HTPB" if contour_id == contour_ids[0] else None)
            for index in np.flatnonzero(ap[1:] != ap[:-1]):
                axis.axvline(
                    0.5 * (arc[index] + arc[index + 1]), color="0.55",
                    linewidth=0.6, linestyle=":")

    axes[0].axhline(0.0, color="0.5", linewidth=0.5)
    axes[1].axhline(0.0, color="0.5", linewidth=0.5)
    axes[0].set_ylabel("Heat flux into solid [MW/m$^2$]")
    axes[1].set_ylabel("Regression mass flux [kg/m$^2$/s]")
    axes[1].set_xlabel("Distance along eta=0.5 contour [mm]")
    axes[0].legend(loc="best", frameon=False, ncol=2)
    for axis in axes:
        axis.grid(True, linewidth=0.35, alpha=0.35)
        axis.margins(x=0.01)

    time_ms = 1.0e3 * float(rows[0]["time_s"])
    figure.suptitle(f"{rows[0]['plotfile']}   t = {time_ms:.4f} ms")
    figure.subplots_adjust(left=0.12, right=0.98, top=0.92, bottom=0.10)
    output_path = output_directory / f"{csv_path.stem}_flux.png"
    figure.savefig(output_path, dpi=dpi)
    plt.close(figure)
    return csv_path.name, len(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Plot heat and mass flux along each extracted burn surface")
    parser.add_argument(
        "output_directory", nargs="?", default=Path(__file__).resolve().parent,
        type=Path)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--dpi", type=int, default=130)
    parser.add_argument("--max-outputs", type=int, default=None)
    args = parser.parse_args()

    root = args.output_directory.resolve()
    csv_paths = sorted((root / "surface_data").glob("*_surface.csv"))
    if args.max_outputs is not None:
        csv_paths = csv_paths[:args.max_outputs]
    if not csv_paths:
        raise RuntimeError(f"No surface CSV files found beneath {root}")

    output_directory = root / "surface_flux_plots"
    output_directory.mkdir(exist_ok=True)
    work = [(path, output_directory, args.dpi) for path in csv_paths]
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        for index, (name, points) in enumerate(executor.map(render, work), start=1):
            print(
                f"[{index}/{len(work)}] {name}: {points} points",
                flush=True)


if __name__ == "__main__":
    main()

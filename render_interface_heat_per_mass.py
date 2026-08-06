#!/usr/bin/env python3
"""Render heat flux per regression mass flux versus interface distance."""

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

    distance = 1.0e3 * np.asarray([
        float(row["signed_distance_to_interface_m"]) for row in rows])
    heat_flux = np.asarray([
        float(row["heat_flux_into_solid_W_m2"]) for row in rows])
    mass_flux = np.asarray([
        float(row["mass_flux_kg_m2_s"]) for row in rows])
    heat_per_mass = 1.0e-6 * heat_flux / mass_flux
    ap = np.asarray([row["species"] == "AP" for row in rows])

    figure, axis = plt.subplots(figsize=(9.0, 5.0))
    axis.scatter(
        distance[ap], heat_per_mass[ap], s=12.0, color="#c43c39",
        linewidths=0.0, label="AP")
    axis.scatter(
        distance[~ap], heat_per_mass[~ap], s=12.0, color="#3569b7",
        linewidths=0.0, label="HTPB")
    axis.axvline(0.0, color="0.2", linewidth=0.8)
    axis.axhline(0.0, color="0.5", linewidth=0.5)
    axis.grid(True, linewidth=0.35, alpha=0.35)
    axis.margins(x=0.02)
    axis.set_ylabel("Heat flux / mass flux [MJ/kg]")
    axis.set_xlabel(
        "Signed distance from AP/HTPB interface [mm]\n"
        "HTPB $\\leftarrow$ 0 $\\rightarrow$ AP")
    axis.legend(loc="best", frameon=False, ncol=2)

    time_ms = 1.0e3 * float(rows[0]["time_s"])
    axis.set_title(f"{rows[0]['plotfile']}   t = {time_ms:.4f} ms")
    figure.tight_layout()
    output_path = output_directory / f"{csv_path.stem}_heat_per_mass.png"
    figure.savefig(output_path, dpi=dpi)
    plt.close(figure)
    return csv_path.name, len(rows)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot heat flux per mass flux against signed interface distance"))
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

    output_directory = root / "surface_interface_heat_per_mass_plots"
    output_directory.mkdir(exist_ok=True)
    work = [(path, output_directory, args.dpi) for path in csv_paths]
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        for index, (name, points) in enumerate(executor.map(render, work), start=1):
            print(
                f"[{index}/{len(work)}] {name}: {points} points",
                flush=True)


if __name__ == "__main__":
    main()

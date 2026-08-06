#!/usr/bin/env python3
"""Superimpose all heat-flux/mass-flux surface data in one figure."""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description="Superimpose all heat-per-mass interface data")
    parser.add_argument(
        "output_directory", nargs="?", default=Path(__file__).resolve().parent,
        type=Path)
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--ymin", type=float, default=0.0)
    parser.add_argument("--ymax", type=float, default=5.0)
    parser.add_argument("--density-xbins", type=int, default=320)
    parser.add_argument("--density-ybins", type=int, default=260)
    args = parser.parse_args()

    root = args.output_directory.resolve()
    csv_paths = sorted((root / "surface_data").glob("*_surface.csv"))
    if not csv_paths:
        raise RuntimeError(f"No surface CSV files found beneath {root}")

    distance = []
    heat_per_mass = []
    ap = []
    times = []
    for index, csv_path in enumerate(csv_paths, start=1):
        with csv_path.open() as stream:
            for row in csv.DictReader(stream):
                distance.append(
                    1.0e3 * float(row["signed_distance_to_interface_m"]))
                heat_per_mass.append(
                    1.0e-6 * float(row["heat_flux_into_solid_W_m2"])
                    / float(row["mass_flux_kg_m2_s"]))
                ap.append(row["species"] == "AP")
                times.append(float(row["time_s"]))
        if index % 250 == 0 or index == len(csv_paths):
            print(f"[{index}/{len(csv_paths)}] {csv_path.name}", flush=True)

    distance = np.asarray(distance)
    heat_per_mass = np.asarray(heat_per_mass)
    ap = np.asarray(ap)
    times = np.asarray(times)
    finite = np.isfinite(distance) & np.isfinite(heat_per_mass)
    distance = distance[finite]
    heat_per_mass = heat_per_mass[finite]
    ap = ap[finite]
    times = times[finite]

    figure, axis = plt.subplots(figsize=(10.0, 6.0))
    axis.scatter(
        distance[~ap], heat_per_mass[~ap], s=1.5, color="#3569b7",
        alpha=0.025, linewidths=0.0, rasterized=True, label="HTPB")
    axis.scatter(
        distance[ap], heat_per_mass[ap], s=1.5, color="#c43c39",
        alpha=0.025, linewidths=0.0, rasterized=True, label="AP")
    axis.axvline(0.0, color="0.2", linewidth=0.8)
    axis.axhline(0.0, color="0.5", linewidth=0.5)
    axis.set_ylim(args.ymin, args.ymax)
    axis.grid(True, linewidth=0.35, alpha=0.35)
    axis.set_ylabel("Heat flux / mass flux [MJ/kg]")
    axis.set_xlabel(
        "Signed distance from AP/HTPB interface [mm]\n"
        "HTPB $\\leftarrow$ 0 $\\rightarrow$ AP")
    axis.legend(loc="best", frameon=False, ncol=2, markerscale=6.0)
    axis.set_title(
        f"All outputs: {len(heat_per_mass):,} samples, "
        f"t = {1.0e3 * times.min():.4f}--{1.0e3 * times.max():.4f} ms")
    figure.tight_layout()
    figure.savefig(root / "surface_interface_heat_per_mass_all.png", dpi=args.dpi)
    plt.close(figure)

    visible = ((heat_per_mass >= args.ymin) &
               (heat_per_mass <= args.ymax))
    x_edges = np.linspace(
        distance.min(), distance.max(), args.density_xbins + 1)
    y_edges = np.linspace(args.ymin, args.ymax, args.density_ybins + 1)
    counts, _, _ = np.histogram2d(
        distance[visible], heat_per_mass[visible], bins=(x_edges, y_edges))
    x_bin = np.clip(
        np.searchsorted(x_edges, distance[visible], side="right") - 1,
        0, args.density_xbins - 1)
    y_bin = np.clip(
        np.searchsorted(y_edges, heat_per_mass[visible], side="right") - 1,
        0, args.density_ybins - 1)
    density = counts[x_bin, y_bin]
    order = np.argsort(density)

    figure, axis = plt.subplots(figsize=(10.0, 6.0))
    points = axis.scatter(
        distance[visible][order], heat_per_mass[visible][order],
        c=density[order], s=2.0, cmap="jet", vmin=1.0,
        vmax=density.max(), alpha=0.65, linewidths=0.0, rasterized=True)
    axis.axvline(0.0, color="0.2", linewidth=0.8)
    axis.axhline(0.0, color="0.5", linewidth=0.5)
    axis.set_ylim(args.ymin, args.ymax)
    axis.grid(True, linewidth=0.35, alpha=0.35)
    axis.set_ylabel("Heat flux / mass flux [MJ/kg]")
    axis.set_xlabel(
        "Signed distance from AP/HTPB interface [mm]\n"
        "HTPB $\\leftarrow$ 0 $\\rightarrow$ AP")
    axis.set_title(
        f"All outputs: {visible.sum():,} displayed samples, "
        f"t = {1.0e3 * times.min():.4f}--{1.0e3 * times.max():.4f} ms")
    colorbar = figure.colorbar(points, ax=axis)
    colorbar.set_label(
        f"Local sample count ({args.density_xbins} x "
        f"{args.density_ybins} bins)")
    figure.tight_layout()
    figure.savefig(
        root / "surface_interface_heat_per_mass_density.png", dpi=args.dpi)
    plt.close(figure)


if __name__ == "__main__":
    main()

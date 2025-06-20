#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt


def load_thermo(folder: str) -> dict[str, np.ndarray]:
    file_path = os.path.join(folder, "thermo.dat")
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"{file_path} not found")

    with open(file_path, "r") as f:
        headers = f.readline().strip().split()
    data = np.loadtxt(file_path, skiprows=1)
    return {h: data[:, i] for i, h in enumerate(headers)}


def save_plots_per_timestep(folder: str, out_dir: str, fmt: str, dpi: int):
    data = load_thermo(folder)
    strain = data["disp_yhi_y"] / 4
    stress = data["trac_yhi_y"] / 2
    label = os.path.basename(folder.rstrip('/'))

    # Create output subfolder for this case
    save_dir = os.path.join(out_dir, label)
    os.makedirs(save_dir, exist_ok=True)

    for t in range(1, len(strain) + 1):
        plt.figure(figsize=(6, 4))
        plt.plot(-strain[:t], -stress[:t], color="blue", lw=2)
        plt.xlabel("Displacement")
        plt.ylabel("Traction")
        plt.title(f"{label} — Step {t}/{len(strain)}")
        plt.tight_layout()

        fig_path = os.path.join(save_dir, f"{label}_step_{t:04d}.{fmt}")
        plt.savefig(fig_path, dpi=dpi)
        plt.close()


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Save traction–displacement plots at each timestep.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("folders", nargs="+", help="Simulation output folders")
    parser.add_argument("--out_dir", default="stepwise_figures", help="Directory to save all plots")
    parser.add_argument("--fmt", default="png", choices=["png", "pdf", "svg"], help="Image format")
    parser.add_argument("--dpi", type=int, default=150, help="Resolution for raster formats")

    args = parser.parse_args(argv)

    for folder in args.folders:
        try:
            save_plots_per_timestep(folder, args.out_dir, args.fmt, args.dpi)
            print(f"Saved timestep plots for {folder}")
        except Exception as e:
            print(f"Error in {folder}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()


#!/usr/bin/env python3
# ────────────────────────────────────────────────────────────────
#  traction disp curve

# ────────────────────────────────────────────────────────────────
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


def plot_stress_strain_combined(folders: list[str],
                                 out_dir: str = "figures",
                                 fmt: str = "png",
                                 dpi: int = 300) -> str:

    plt.figure(figsize=(6, 4))

    for folder in folders:
        try:
            data = load_thermo(folder)
            strain = data["disp_yhi_y"] / 4
            stress = data["trac_yhi_y"] / 2
            label = os.path.basename(folder.rstrip('/'))
            plt.plot(-strain, -stress, lw=1.5, label=label)
        except Exception as e:
            print(f"Failed on {folder}: {e}", file=sys.stderr)

    plt.xlabel("Displacement")
    plt.ylabel("Traction")
    plt.title("Traction vs Displacement")
    plt.legend()
    plt.tight_layout()

    os.makedirs(out_dir, exist_ok=True)
    fig_name = "combined_trac_disp." + fmt
    fig_path = os.path.join(out_dir, fig_name)
    plt.savefig(fig_path, dpi=dpi, bbox_inches="tight")
    plt.close()

    return fig_path


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Plot combined disp–trac curve from multiple thermo.dat folders.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("folders", nargs="+",
                        help="Simulation output folders")
    parser.add_argument("--out_dir", default="figures",
                        help="Directory where the figure is written")
    parser.add_argument("--fmt", default="png", choices=["png", "pdf", "svg"],
                        help="Image format for saved figure")
    parser.add_argument("--dpi", type=int, default=300,
                        help="Resolution for raster formats")

    args = parser.parse_args(argv)

    fig = plot_stress_strain_combined(args.folders, args.out_dir, args.fmt, args.dpi)
    print(f"saved combined plot: {fig}")

#main
if __name__ == "__main__":
    main()



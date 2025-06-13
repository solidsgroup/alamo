#!/usr/bin/env python3
# ────────────────────────────────────────────────────────────────
#  stress_strain_plot.py
#
#  Read thermo.dat from one or more simulation folders, plot the
#  (disp_yhi_x / 4 , trac_yhi_x / 2) “engineering” stress–strain
#  curve(s), and save the figure(s) to disk.
#
#  • Run as a script  →  python stress_strain_plot.py output_*
#  • Use in a notebook →  just execute the bottom-half cells
# ────────────────────────────────────────────────────────────────
import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

# ────────────────────────────────────────────────────────────────
#  I/O helper
# ────────────────────────────────────────────────────────────────
def load_thermo(folder: str) -> dict[str, np.ndarray]:
    """
    Read `thermo.dat` in *folder* and return a dict of column vectors.

    Assumes first row contains whitespace-delimited headers.
    """
    file_path = os.path.join(folder, "thermo.dat")
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"{file_path} not found")

    with open(file_path, "r") as f:
        headers = f.readline().strip().split()
    data = np.loadtxt(file_path, skiprows=1)

    return {h: data[:, i] for i, h in enumerate(headers)}


# ────────────────────────────────────────────────────────────────
#  core plotting routine
# ────────────────────────────────────────────────────────────────
def plot_stress_strain(folder: str,
                       out_dir: str = "figures",
                       fmt: str = "png",
                       dpi: int = 300) -> str:
    """
    Plot stress–strain curve for *folder* and save to *out_dir*.

    Returns the saved figure path.
    """
    data = load_thermo(folder)

    # engineering strain & stress (scale factors are problem-specific)
    strain = data["disp_yhi_y"] /4
    stress = data["trac_yhi_y"] /2

    # prepare output path
    os.makedirs(out_dir, exist_ok=True)
    fig_name = f"{os.path.basename(folder.rstrip('/'))}_stress_strain.{fmt}"
    fig_path = os.path.join(out_dir, fig_name)

    # create the plot
    plt.figure(figsize=(4, 3))
    plt.plot(-(strain), -(stress), lw=1.6, label=folder)
    plt.xlabel(r"Displacement")
    plt.ylabel(r"Traction")
    plt.title("Traction vs displacement")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=dpi, bbox_inches="tight")
    plt.close()

    return fig_path


# ────────────────────────────────────────────────────────────────
#  command-line entry point
# ────────────────────────────────────────────────────────────────
def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Plot stress–strain curves from thermo.dat files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("folders", nargs="+",
                        help="Simulation output folders (each must contain thermo.dat)")
    parser.add_argument("--out_dir", default="figures",
                        help="Directory where figures are written")
    parser.add_argument("--fmt", default="png", choices=["png", "pdf", "svg"],
                        help="Image format for saved figure")
    parser.add_argument("--dpi", type=int, default=300,
                        help="Resolution (DPI) for raster formats")

    args = parser.parse_args(argv)

    for f in args.folders:
        try:
            fig = plot_stress_strain(f, args.out_dir, args.fmt, args.dpi)
            print(f"[✓] saved {fig}")
        except Exception as exc:
            print(f"[✗] {f}: {exc}", file=sys.stderr)


# ────────────────────────────────────────────────────────────────
#  notebook-friendly usage
# ────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    main()

"""
### Minimal notebook cell -----------------------------------------
import stress_strain_plot as ssp
ssp.plot_stress_strain("output")      # writes figures/output_stress_strain.png
"""


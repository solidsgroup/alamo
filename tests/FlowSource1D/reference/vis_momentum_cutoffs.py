#!/usr/bin/env python3
# Plot momentumx vs x for each momentum-source-term reference file, to
# compare how the cutoff parameter affects the simulation result. A second
# row zooms into the x-region where the cutoff value changes the result most.
import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

script_dir = os.path.dirname(os.path.abspath(__file__))


def cutoff_label(path):
    name = os.path.basename(path)
    m = re.search(r"momentum-source-term-([0-9-]+)\.dat", name)
    if not m:
        return "baseline"
    return "cutoff=" + m.group(1).replace("-", ".", 1).replace("-", "")


def zoom_region(dfs, column="momentumx", spread_frac=0.05, pad_frac=0.15):
    """Find the x-range where the spread across files is largest."""
    x_lo = max(df["x"].min() for df in dfs)
    x_hi = min(df["x"].max() for df in dfs)
    grid = np.linspace(x_lo, x_hi, 500)
    values = np.array([np.interp(grid, df["x"], df[column]) for df in dfs])
    spread = values.max(axis=0) - values.min(axis=0)

    threshold = spread_frac * spread.max()
    mask = spread > threshold
    if not mask.any():
        return x_lo, x_hi

    region_x = grid[mask]
    lo, hi = region_x.min(), region_x.max()
    pad = pad_frac * (hi - lo)
    return max(x_lo, lo - pad), min(x_hi, hi + pad)


files = sorted(glob.glob(os.path.join(script_dir, "momentum-source-term*.dat")))
dfs = [pd.read_csv(path) for path in files]
labels = [cutoff_label(path) for path in files]

zoom_lo, zoom_hi = zoom_region(dfs, column="momentumx")

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
columns = ["momentumx", "density", "energy"]

for df, label in zip(dfs, labels):
    for col, ax in zip(columns, axes[0]):
        ax.plot(df["x"], df[col], marker="o", markersize=3, label=label)
    for col, ax in zip(columns, axes[1]):
        ax.plot(df["x"], df[col], marker="o", markersize=3, label=label)

for col, ax in zip(columns, axes[0]):
    ax.set_ylabel(col)
    ax.set_xlabel("x")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title(col)

for col, ax in zip(columns, axes[1]):
    ax.set_ylabel(col)
    ax.set_xlabel("x")
    ax.set_xlim(zoom_lo, zoom_hi)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title(f"{col} (zoomed: x in [{zoom_lo:.3f}, {zoom_hi:.3f}])")

fig.suptitle("FlowSource1D momentum source term: effect of cutoff")
fig.tight_layout()

outpath = os.path.join(script_dir, "momentum_cutoffs.png")
fig.savefig(outpath, dpi=150)
print(f"Wrote {outpath}")
plt.show()

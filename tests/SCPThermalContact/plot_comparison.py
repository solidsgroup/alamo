#!/usr/bin/env python3
"""
Plot the simulation's final temperature profile against the analytical
Robin/convective-boundary solution (see generate_reference.py and
Readme.rst). Run this after `./bin/alamo-... tests/SCPThermalContact/input`
has produced tests/SCPThermalContact/output/; writes reference/comparison.png.
"""
import sys
import glob
import numpy as np
import yt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from generate_reference import T_robin, stop_time, y_lo, y_hi

outdir = sys.argv[1] if len(sys.argv) > 1 else "output"
path = sorted(glob.glob("{}/*cell/".format(outdir)))[-1]

ds = yt.load(path)
t = float(ds.current_time)
ad = ds.all_data()
y = np.array(ad[("index", "y")].to("code_length"))
T = np.array(ad[("boxlib", "temperature")])
order = np.argsort(y)
y_sorted, T_sorted = y[order], T[order]
uniq_y = np.unique(y_sorted)
T_avg = np.array([T_sorted[y_sorted == yy].mean() for yy in uniq_y])

y_fine = np.linspace(y_lo, y_hi, 400)
depth_fine = np.where(y_fine < 0, -y_fine, 0.0)
T_fine = np.where(y_fine < 0, T_robin(depth_fine, t), np.nan)

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(uniq_y * 1e9, T_avg, 'o', color='C1', markerfacecolor='none',
        label='Alamo (Hydro::AdvanceSolidEnergy)')
ax.plot(y_fine * 1e9, T_fine, '-', color='C0', linewidth=2,
        label='Analytic (Carslaw & Jaeger)')
ax.axvline(0.0, color='gray', linestyle=':', linewidth=1)
ax.set_xlabel("y [nm]  (solid: y<0, gas: y>0)")
ax.set_ylabel("Temperature [K]")
ax.set_title(f"Solid/gas interfacial conduction, t={t:.2e} s")
ax.legend()
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("reference/comparison.png", dpi=150)
print("wrote reference/comparison.png")

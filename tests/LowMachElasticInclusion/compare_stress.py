#!/usr/bin/env python3
"""
Plot simulated vs. analytic displacement and stress for
LowMachElasticInclusion.

Samples the same y~=0 ray as `test` from a finished run's final plotfile,
overlays it against the closed-form solution (far-field confined-compression
state plus the circular-inhomogeneity correction near the AP inclusion --
see Readme.rst / generate_reference.py for the derivation), and writes a
2x2-panel PNG: disp_x comparison (top-left) and stress_xx/stress_yy
comparison (top-right), each with its absolute error directly below it
(with the excluded interface band shaded).

Usage:
    python3 compare_stress.py <outdir> [output.png]

<outdir> is a run directory containing timestamped `*cell/` plotfiles (i.e.
what you point `test` at) -- e.g. the directory scripts/runtests.py or a
manual `lowmach input ... plot_file=<outdir>/...` run writes to.
"""
import glob
import sys

sys.path.insert(0, "../../scripts")
import numpy as np
import matplotlib.pyplot as plt
import testlib
from generate_reference import analytic, a, w, sigma_xx_in, sigma_yy_in, p0, s

if len(sys.argv) < 2:
    raise SystemExit(f"usage: {sys.argv[0]} <outdir> [output.png]")
outdir = sys.argv[1]
outpng = sys.argv[2] if len(sys.argv) > 2 else "comparison.png"

plotfiles = sorted(glob.glob(f"{outdir}/*cell/"))
if len(plotfiles) < 1:
    raise RuntimeError(f"no plotfiles found under {outdir}")
final = plotfiles[-1]

# --- Sample simulation along y~=0 (same ray as `test`) ---------------------
x_lo, x_hi = -1.3e-3, 1.3e-3
y_ray = 9.375e-6

df = testlib.readContours(
    path=final,
    start=[x_lo, y_ray, 0.0],
    end=[x_hi, y_ray, 0.0],
    vars=["disp_x", "disp_y", "stress_xx", "stress_yy"],
)
x = df["x"].to_numpy()
r = np.abs(x)
disp_x = df["disp_x"].to_numpy()
stress_xx = df["stress_xx"].to_numpy()
stress_yy = df["stress_yy"].to_numpy()

band = 3.0 * w

# --- Analytic curve on a fine grid (for smooth plotting) -------------------
x_fine = np.linspace(x_lo, x_hi, 2000)
ux_fine, _, sxx_fine, syy_fine = analytic(x_fine)

# --- Analytic curve on the sim's own sample points (for error panels) ------
ux_exact, _, sxx_exact, syy_exact = analytic(x)
err_ux = disp_x - ux_exact
err_xx = stress_xx - sxx_exact
err_yy = stress_yy - syy_exact

# --- Plot --------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True,
                          gridspec_kw={"height_ratios": [2.2, 1]})

title = (
    f"LowMachElasticInclusion: simulated vs. analytic displacement/stress (y~=0 ray)\n"
    f"AP in HTPB block, top-loaded confined compression, R/w={a/w:g}, "
    f"p0={p0:.4f}, s={s:.4f}, sigma_xx_in={sigma_xx_in:.4f}, sigma_yy_in={sigma_yy_in:.4f}"
)
fig.suptitle(title, fontsize=10)

# -- disp_x --
ax = axes[0, 0]
ax.plot(x * 1e3, disp_x, '-', color='tab:green', lw=1.5, label='sim disp_x')
ax.plot(x_fine * 1e3, ux_fine, '--', color='k', lw=1.2, label='analytic disp_x')
ax.axvspan(-(a + band) * 1e3, -(a - band) * 1e3, color='red', alpha=0.12)
ax.axvspan((a - band) * 1e3, (a + band) * 1e3, color='red', alpha=0.12)
ax.set_ylabel("disp_x (m)")
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.18), ncol=2, fontsize=8, frameon=False)

ax = axes[1, 0]
ax.plot(x * 1e3, np.abs(err_ux), '-', color='tab:green', lw=1.2, label='|err| disp_x')
ax.axvspan(-(a + band) * 1e3, -(a - band) * 1e3, color='red', alpha=0.12, label='excluded interface band')
ax.axvspan((a - band) * 1e3, (a + band) * 1e3, color='red', alpha=0.12)
ax.set_xlabel("x (mm)")
ax.set_ylabel("|absolute error| (m)")
ax.legend(loc='upper right', fontsize=8)

# -- stress_xx / stress_yy --
ax = axes[0, 1]
ax.plot(x * 1e3, stress_xx, '-', color='tab:blue', lw=1.5, label='sim stress_xx')
ax.plot(x * 1e3, stress_yy, '-', color='tab:orange', lw=1.5, label='sim stress_yy')
ax.plot(x_fine * 1e3, sxx_fine, '--', color='k', lw=1.2, label='analytic stress_xx')
ax.plot(x_fine * 1e3, syy_fine, ':', color='k', lw=1.6, label='analytic stress_yy')
ax.axvspan(-(a + band) * 1e3, -(a - band) * 1e3, color='red', alpha=0.12)
ax.axvspan((a - band) * 1e3, (a + band) * 1e3, color='red', alpha=0.12)
ax.set_ylabel("stress / P")
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.28), ncol=2, fontsize=8, frameon=False)

ax = axes[1, 1]
ax.plot(x * 1e3, np.abs(err_xx), '-', color='tab:blue', lw=1.2, label='|err| stress_xx')
ax.plot(x * 1e3, np.abs(err_yy), '-', color='tab:orange', lw=1.2, label='|err| stress_yy')
ax.axvspan(-(a + band) * 1e3, -(a - band) * 1e3, color='red', alpha=0.12)
ax.axvspan((a - band) * 1e3, (a + band) * 1e3, color='red', alpha=0.12)
ax.set_xlabel("x (mm)")
ax.set_ylabel("|absolute error|")
ax.legend(loc='upper right', fontsize=8)

fig.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig(outpng, dpi=150)
print(f"Wrote {outpng}")

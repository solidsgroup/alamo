#!/usr/bin/env python3
"""
Plot simulated vs. analytic stress for LowMachElasticInclusion.

Samples the same y~=0 ray as `test` from a finished run's final plotfile,
overlays it against the closed-form solution (far-field confined-compression
state plus the circular-inhomogeneity correction near the AP inclusion --
see Readme.rst / generate_reference.py for the derivation), and writes a
two-panel PNG: stress_xx/stress_yy comparison on top, absolute error on the
bottom (with the excluded interface band shaded).

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
stress_xx = df["stress_xx"].to_numpy()
stress_yy = df["stress_yy"].to_numpy()

band = 3.0 * w

# --- Analytic curve on a fine grid (for smooth plotting) -------------------
x_fine = np.linspace(x_lo, x_hi, 2000)
_, _, sxx_fine, syy_fine = analytic(x_fine)

# --- Analytic curve on the sim's own sample points (for error panel) -------
_, _, sxx_exact, syy_exact = analytic(x)
err_xx = stress_xx - sxx_exact
err_yy = stress_yy - syy_exact

# --- Plot --------------------------------------------------------------
fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True,
                          gridspec_kw={"height_ratios": [2.2, 1]})

ax = axes[0]
ax.plot(x * 1e3, stress_xx, '-', color='tab:blue', lw=1.5, label='sim stress_xx')
ax.plot(x * 1e3, stress_yy, '-', color='tab:orange', lw=1.5, label='sim stress_yy')
ax.plot(x_fine * 1e3, sxx_fine, '--', color='k', lw=1.2, label='analytic stress_xx')
ax.plot(x_fine * 1e3, syy_fine, ':', color='k', lw=1.6, label='analytic stress_yy')
ax.axvspan(-(a + band) * 1e3, -(a - band) * 1e3, color='red', alpha=0.12)
ax.axvspan((a - band) * 1e3, (a + band) * 1e3, color='red', alpha=0.12)
ax.set_ylabel("stress / P")
ax.set_title(
    f"LowMachElasticInclusion: simulated vs. analytic stress (y~=0 ray)\n"
    f"AP in HTPB block, top-loaded confined compression, R/w={a/w:g}, "
    f"p0={p0:.4f}, s={s:.4f}, sigma_xx_in={sigma_xx_in:.4f}, sigma_yy_in={sigma_yy_in:.4f}",
    fontsize=10,
)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.32), ncol=4, fontsize=8, frameon=False)

ax2 = axes[1]
ax2.plot(x * 1e3, np.abs(err_xx), '-', color='tab:blue', lw=1.2, label='|err| stress_xx')
ax2.plot(x * 1e3, np.abs(err_yy), '-', color='tab:orange', lw=1.2, label='|err| stress_yy')
ax2.axvspan(-(a + band) * 1e3, -(a - band) * 1e3, color='red', alpha=0.12, label='excluded interface band')
ax2.axvspan((a - band) * 1e3, (a + band) * 1e3, color='red', alpha=0.12)
ax2.set_xlabel("x (mm)")
ax2.set_ylabel("|absolute error|")
ax2.legend(loc='upper right', fontsize=8)

fig.tight_layout()
fig.savefig(outpng, dpi=150)
print(f"Wrote {outpng}")

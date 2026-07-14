import math
from pathlib import Path

import alamo
import numpy
import pylab

alamo.include("Util/Util.H")
alamo.include("Numeric/Advect/Limiter/Koren.H")
alamo.include("Numeric/Advect/Limiter/MinMod.H")
alamo.include("Numeric/Advect/Limiter/MC.H")
alamo.include("Numeric/Advect/Limiter/Superbee.H")
alamo.include("Numeric/Advect/Limiter/UMIST.H")
alamo.include("Numeric/Advect/Limiter/VanAlbada.H")
alamo.include("Numeric/Advect/Limiter/VanLeer.H")


alamo.Util.Initialize()

limiters = [
    ("minmod", alamo.Numeric.Advect.Limiter.MinMod(),
     lambda r: numpy.maximum(0.0, numpy.minimum(r, 1.0))),
    ("MC", alamo.Numeric.Advect.Limiter.MC(),
     lambda r: numpy.maximum(0.0, numpy.minimum.reduce([2.0 * r, 0.5 * (1.0 + r), 2.0 + 0.0 * r]))),
    ("superbee", alamo.Numeric.Advect.Limiter.Superbee(),
     lambda r: numpy.maximum(0.0, numpy.maximum(numpy.minimum(2.0 * r, 1.0), numpy.minimum(r, 2.0)))),
    ("van Leer", alamo.Numeric.Advect.Limiter.VanLeer(),
     lambda r: (r + numpy.abs(r)) / (1.0 + numpy.abs(r))),
    ("van Albada", alamo.Numeric.Advect.Limiter.VanAlbada(),
     lambda r: numpy.where(r > 0.0, (r * r + r) / (r * r + 1.0), 0.0)),
    ("Koren", alamo.Numeric.Advect.Limiter.Koren(),
     lambda r: numpy.maximum(0.0, numpy.minimum.reduce([2.0 * r, (1.0 + 2.0 * r) / 3.0, 2.0 + 0.0 * r]))),
    ("UMIST", alamo.Numeric.Advect.Limiter.UMIST(),
     lambda r: numpy.maximum(
         0.0,
         numpy.minimum.reduce([2.0 * r, (1.0 + 3.0 * r) / 4.0, (3.0 + r) / 4.0, 2.0 + 0.0 * r]))),
]

r = numpy.linspace(0.0, 3.0, 301)
tvd_lo = next(expected_fn for name, _, expected_fn in limiters if name == "minmod")(r)
tvd_hi = next(expected_fn for name, _, expected_fn in limiters if name == "superbee")(r)
outdir = Path(__file__).with_name("output")
outdir.mkdir(exist_ok=True)

ncols = math.ceil(math.sqrt(len(limiters)))
nrows = math.ceil(len(limiters) / ncols)
fig, axes = pylab.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.5 * nrows), squeeze=False)
axes = axes.flat

for ax, (name, limiter, expected_fn) in zip(axes, limiters):
    # The C++ limiter returns a limited slope from (lo, mid, hi).
    # Choose lo=0, mid=r, hi=r+1 so the downwind slope is one and
    # the returned slope is directly the Sweby limiter function phi(r).
    phi = numpy.array([limiter(0.0, float(_r), float(_r + 1.0)) for _r in r])
    phi_decreasing = numpy.array([limiter(0.0, -float(_r), -float(_r + 1.0)) for _r in r])
    expected = expected_fn(r)

    assert numpy.max(numpy.abs(phi - expected)) < 1.0e-13
    assert numpy.max(numpy.abs(phi_decreasing + expected)) < 1.0e-13
    assert abs(limiter(0.0, -1.0, 0.0)) < 1.0e-13
    assert numpy.all(phi >= -1.0e-13)
    assert numpy.all(phi <= numpy.minimum(2.0 * r, 2.0) + 1.0e-13)

    ax.fill_between(r, tvd_lo, tvd_hi, color="gray", alpha=0.18, label="TVD region")
    ax.plot(r, tvd_lo, color="gray", linestyle="--", linewidth=0.9)
    ax.plot(r, tvd_hi, color="gray", linestyle="--", linewidth=0.9)
    ax.plot(r, expected, color="black", linewidth=1.8, label="analytic")
    ax.plot(
        r,
        phi,
        linestyle="None",
        marker="o",
        markevery=10,
        markersize=4.0,
        markerfacecolor="none",
        markeredgecolor="tab:blue",
        markeredgewidth=1.0,
        label="C++",
    )
    ax.set_title(name)
    ax.set_xlabel("r")
    ax.set_ylabel(r"$\phi(r)$")
    ax.set_xlim(0.0, 3.0)
    ax.set_ylim(-0.05, 2.1)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left", fontsize=8)

for ax in axes[len(limiters):]:
    ax.axis("off")

fig.suptitle("Slope limiters")
fig.tight_layout()
fig.savefig(outdir / "limiters.png", dpi=160)
alamo.Util.Finalize()

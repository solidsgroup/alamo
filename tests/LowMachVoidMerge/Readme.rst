LowMach void-merge regression test
==================================

This is a small, fast version of ``input.lm.ap_void``: a circular void carved
into a solid AP domain, sized so the real regressing solid/gas front reaches
and merges with it within a short ``stop_time``. It exists to catch two
things in one run:

1. The pressure-projection MLMG failure that ``input.lm.ap_void`` originally
   hit -- an unreachably tight ``projection.tol_rel`` given this problem's
   RHS dynamic range (see ``input.lm.ap_void``'s header comment for the full
   derivation and the empirical evidence). No existing LowMach test exercised
   ``rigid_solid``/``mechanisms`` before this one, so this path was untested.
2. The "front reaches a pre-existing void" merge path, end to end.

Run it with::

    ./configure --dim=2 --yaml
    make -j bin/lowmach
    PATH=/path/to/python-with-yt/bin:$PATH \
        ./scripts/runtests.py tests/LowMachVoidMerge --serial --dim=2

How to define a void
---------------------

A void is not a distinct IC type -- it's a smooth multiplicative mask applied
consistently across every field that needs to "know" about it: the solid
density, the gas density, and the temperature. All three use the same
``xc``/``yc``/``Rvoid``/``wvoid`` constants and the same mask shape::

    voidmask_out = 0.5 + 0.5*tanh((r - Rvoid) / wvoid)     # ~0 inside, ~1 outside
    r = sqrt((x-xc)^2 + (y-yc)^2)

- The **solid** density IC multiplies its normal (non-void) profile by
  ``voidmask_out`` (zero inside the void).
- The **gas** density IC adds a ``(1 - voidmask_out)`` term so the void
  interior is filled with product gas at the ambient pressure.
- The **temperature** IC multiplies its preheat term by ``voidmask_out`` so
  the void interior sits at ambient temperature instead of the preheated
  profile used elsewhere.

See ``AP_solid.density.ic``, ``Final.density.ic``, and ``temperature.ic`` in
this test's ``input`` file (or ``input.lm.ap_void``) for the literal
expressions.

Sizing rule (read this before changing Rvoid/wvoid)
-----------------------------------------------------

``wvoid`` is **not** a free/arbitrary parameter -- it must match the
Allen-Cahn phase-field model's own natural interface width, or the phase
field starts far out of equilibrium at the void boundary and the
``AP_decomposition`` mechanism fires an artificial relaxation rate large
enough to destabilize the pressure projection (this was the original root
cause of ``input.lm.ap_void``'s MLMG failure, on top of the tolerance issue
above). For the quartic double-well built from ``(w0, w12, w1)``::

    a2 = -11*w0 + 16*w12 - 5*w1
    delta = 2*sqrt(kappa / (lambda * 2*a2))

For the parameters used here and in ``input.lm.ap_void``
(``w0=0, w12=2, w1=1, kappa=1.0e-8, lambda=8.333333``): ``a2=27``,
``delta =~ 9.4 um``, rounded to ``wvoid = 10 um``. This also happens to match
the width already used for the ordinary (non-void) solid/gas interface in
these inputs, ``tanh(2*y/20e-6)`` = ``tanh(y/10e-6)`` -- that agreement is the
cross-check that the formula is right.

Given ``wvoid``, size ``Rvoid`` so the void actually opens up (gas fraction at
the center is ``0.5*(1+tanh(Rvoid/wvoid))`` -- at ``Rvoid=2*wvoid`` that's only
98.2%, at ``Rvoid=4*wvoid`` it's 99.97%) and so curvature-driven relaxation of
the circular interface (a real physical effect, not a bug -- small closed
curves in an Allen-Cahn field have an intrinsic surface-tension-like
relaxation) stays gentle: ``Rvoid >= 3*wvoid`` is the minimum, ``4*wvoid`` is
comfortable. Resolve ``wvoid`` with at least 2 level-0 cells and 4 finest-level
cells. Finally, keep the void's full smoothing skirt, ``Rvoid + 3*wvoid``,
clear of both the burning surface and the domain edges by a healthy margin.

This durable design rule is shape-independent (see below): a sharp corner is
locally near-infinite curvature and strictly worse than a smooth curve of the
same size, so it needs the same minimum-radius treatment.

Non-circular voids
-------------------

The circular mask above is the simplest case, not the only one:

- **Any smooth analytic shape** (ellipse, superellipse/rounded square, star,
  etc.) works today with zero code changes -- just substitute a different
  signed-distance-like formula ``f(x,y)`` for
  ``sqrt((x-xc)^2+(y-yc)^2) - Rvoid`` in the ``voidmask_out`` expression
  above. The same curvature-vs-``delta`` sizing rule applies at every point
  of the boundary, not just a single radius.
- **Polygons from a points file** or **image masks** are supported via
  ``IC::PointList`` (``pointlist``) and ``IC::PNG`` (``png``), which are wired
  into ``AP_solid.density.ic.type`` / ``Final.density.ic.type`` /
  ``temperature.ic.type`` alongside ``expression`` and ``constant``.
  ``IC::PointList`` smooths with ``tanh(d/(sqrt(2)*eps))`` rather than
  ``tanh(d/wvoid)``, so match widths with ``eps = wvoid/sqrt(2)``. Give it a
  points file (``x y z [ObjNum]`` per line, one polygon per ``ObjNum`` run)
  and mark the void polygon(s) with ``invert``::

      AP_solid.density.ic.type = pointlist
      AP_solid.density.ic.pointlist.file.name = my_shape.dat
      AP_solid.density.ic.pointlist.file.unit = um
      AP_solid.density.ic.pointlist.eps = 7.07_um
      AP_solid.density.ic.pointlist.value = 1950.0_kg/m^3
      AP_solid.density.ic.pointlist.invert = 0 1   # polygon 0 = solid block, polygon 1 = void

Permanently-inert voids (embedded filler/manufacturing defects)
-------------------------------------------------------------------

The void above is a *physical* gas pocket: nothing prevents the real front
from decomposing right up to its wall and merging with it (and nothing
should -- ``rigid_eta`` is recomputed from density every step with no memory
of "this cell used to be a void", so a merge is not a special event; it's
just what the existing physics does once density on both sides reaches
zero). If instead you want a region that must **never** react even once the
front reaches it -- an inert filler particle, say -- that needs the separate
``psi_mechanism`` gate:

::

    psi.ic.type = expression
    psi.ic.expression.constant.xcd = 70.0_um
    psi.ic.expression.constant.ycd = -70.0_um
    psi.ic.expression.constant.Rdisc = 12.0_um
    psi.ic.expression.constant.wdisc = 5.0_um
    psi.ic.expression.region0 = "0.5+0.5*tanh((sqrt((x-xcd)^2+(y-ycd)^2)-Rdisc)/wdisc)"
    # psi.release_temperature = 500.0_K   # omit for a permanently-inert filler

``psi_mechanism`` is 1 (mechanism fully active) outside the region and 0
(fully inert) inside it, multiplying the ``AP_decomposition`` mechanism's rate
functions directly. Left at its default (``psi.release_temperature`` = a huge
value), the masked region never reacts, no matter how long the front sits
next to it. If ``psi.release_temperature`` is set to something reachable, the
region releases (permanently -- the update is monotone, so it cannot re-freeze
on cooling) once **both** of the following hold, checked once per timestep:

- local temperature exceeds ``psi.release_temperature``, and
- a neighboring cell's real solid (``rigid_eta``) has measurably dropped
  (``psi.release_eta_threshold``, default 0.01).

Both signals are required together deliberately: temperature alone is not
enough, since conduction can heat a region well before any front has
physically consumed neighboring mass, which would release the gate before
anything has actually "arrived".

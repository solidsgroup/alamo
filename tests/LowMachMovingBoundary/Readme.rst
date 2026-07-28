LowMach elastic fluid-structure interaction with a moving boundary
====================================================================

Extension of ``tests/LowMachElastic`` that adds a moving solid/fluid
boundary: instead of a frozen composition, the two rigid solid species
(``AP_solid``, ``HTPB_solid``) regress directly into the single fluid
species (``Gas``) at a fixed, prescribed interface speed. No chemistry
kinetics are involved -- ``chemistry.model.type = frozen`` throughout; the
mass transfer comes entirely from ``Model::Mechanism::PhaseChange``'s
``prescribed_speed`` option, which converts solid to gas at a constant
normal rate regardless of temperature or pressure (see
``Model/Mechanism/PhaseChange.H``: ``MassSource`` computes
``eta_dot = -(input_eta/eta) * prescribed_speed * |grad(rigid_eta)|``
when ``prescribed_speed`` is set, bypassing the Arrhenius/Allen-Cahn rate
path entirely).

This is a regression/smoke check, not a comparison against an analytic
solution: it verifies that ``LowMach``'s elastic FSI solve (added for
``tests/LowMachElastic``) continues to mix the model field and drive the
elastic solver correctly as ``rigid_eta`` itself evolves in time, not just
for a static configuration.

Setup
-----

Same domain and initial condition as ``tests/LowMachElastic``: AP_solid and
HTPB_solid split by a diffuse interface at ``x = x0``, occupying the bottom
half of the domain, with Gas above at a uniform reference pressure (and the
same ideal-gas-consistent Gas density IC needed for LowMach's mixture-volume
pressure-projection constraint -- see that test's Readme for why).

Two ``mechanisms.names`` entries regress each solid species directly to Gas:

- ``AP_regression``: ``AP_solid -> Gas`` at ``prescribed_speed = 4.0 m/s``
- ``HTPB_regression``: ``HTPB_solid -> Gas`` at ``prescribed_speed = 2.0 m/s``

Because the two species regress at different rates, and the local
regression rate at any point is blended by that point's AP/HTPB volume
fraction, the initially-flat solid/gas interface tilts as it retreats
instead of moving uniformly -- a non-planar moving boundary rather than a
simple 1-D translation. The elastic solve uses ``elastic.interval = 5`` (see
``Base::Mechanics``'s ``elastic.interval``/``m_interval`` throttle) so the
300-fixed-iteration MLMG solve does not have to run every timestep to track
the (comparatively slowly) moving boundary.

The ``test`` script checks, using the initial and final plotfiles:

- Total solid volume (``sum(rigid_eta)``) decreases -- mass is genuinely
  being converted to gas, not held constant.
- AP_solid's fractional volume loss (``rigid_species_eta_AP_solid``) exceeds
  HTPB_solid's (``rigid_species_eta_HTPB_solid``), confirming the per-species
  prescribed speed is actually being applied (not e.g. some averaged rate).
  Per-species volume fraction is used instead of an interface-height/grid
  threshold since the run is deliberately short (bounded by the elastic MLMG
  caveat below) and the regression distance is sub-cell.
- On the final, moved plotfile: ``model_mu``/``model_kappa`` equal the
  (shared) AP/HTPB elastic constants everywhere, and
  ``disp``/``stress``/``strain`` are finite and show a nonzero displacement
  response to the (nonzero) interfacial traction -- i.e. the elastic solve
  is still healthy after the boundary has moved.

Numerical caveats
------------------

The Gas density IC must be volume-fraction-consistent, as in
``tests/LowMachElastic`` -- see that test's Readme.rst.

Unlike ``LowMachElastic``, AP_solid and HTPB_solid are given **the same**
elastic constants here rather than a modest contrast. During development, a
modest (~7%) contrast that converges reliably for the *static* interface in
``LowMachElastic`` still triggered the same pre-existing
``Operator::Elastic`` MLMG divergence-for-heterogeneous-media gap partway
through this test, once the interface had retreated and tilted enough to
become more geometrically irregular than the static case. Since this test's
purpose is specifically the moving boundary/elastic-interval tracking (the
arbitrary-species mixing itself is already covered by ``LowMachElastic``), a
spatially uniform elastic model sidesteps that unrelated, already-documented
robustness gap rather than working around it here too.

Even with a spatially uniform elastic model, the elastic MLMG solve
eventually diverged (around ``t ~ 1.0e-5 s``, roughly the eighth periodic
solve at ``elastic.interval = 5``) once the retreating/tilting ``rigid_eta``
interface became irregular enough -- i.e. this is not purely a
material-contrast problem but a broader ``Operator::Elastic``
robustness-to-interface-shape gap. ``stop_time`` is set short enough
(``9.0e-6 s``) to stay well inside the region that reliably converges while
still producing several periodic elastic solves and measurable (if sub-cell)
regression. Extending this test to run longer, or root-causing why an
evolving/thinning psi mask destabilizes the MLMG solve, is a worthwhile
follow-up (same category of issue as the contrast gap noted above and in
``LowMachElastic``).

Run the case with::

    ./configure --dim=2
    make -j bin/lowmach
    ./scripts/runtests.py tests/LowMachMovingBoundary --serial --dim=2

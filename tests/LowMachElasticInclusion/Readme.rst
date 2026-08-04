LowMach elastic circular inclusion (Kolosov-Muskhelishvili validation)
========================================================================

Analytic validation of the LowMach elastic fluid-structure-interaction solve
(``elastic.apply_fluid_pressure``, see ``LowMach::UpdateModel`` and
``tests/LowMachElasticPressure``) against the classical two-phase circular
inhomogeneity solution, obtained via Kolosov-Muskhelishvili complex
potentials.

This is a simplified, single-particle descendant of
``input.lm.ap_htpb_packed_elastic`` (the production AP/HTPB packed-bed
case): one circular AP particle at the center of a square domain, uniform
frozen pressure, uniform 300K temperature, frozen chemistry -- no HTPB, no
regression, no reaction. Where that input drives the elastic solve with the
real local gas pressure field on a complex multi-particle geometry (and has
no closed-form answer to check against), this test isolates the single
piece of physics that *does* have one: a stiff circular inclusion pressure-
loaded at its own boundary, embedded in a softer "matrix."
``AP_solid.elastic.model.mu/kappa`` below are the real (rescaled) AP crystal
moduli from that production file (E ~ 21.6 GPa, nu ~ 0.14).

Setup
-----

A single ``rigid_solid`` species (``AP_solid``) occupies a disk of radius
``a = 1.5e-4`` m at the center of a ``3.0mm x 3.0mm`` square domain; a
single fluid species (``Gas``) fills the rest, at a uniform, time-frozen
pressure ``P`` (``projection.update_pressure = 0``, see
``tests/LowMachElasticPressure/Readme.rst``). Chemistry is frozen and the
temperature field is a uniform 300K everywhere, so the AP/Gas interface
never moves.

Unlike ``tests/LowMachElasticPressure`` (where the fluid region is masked
out of the elastic solve and only the solid modulus matters),
``elastic.void.model.mu/kappa`` here is a **first-class physical
parameter**: the Gas region is given a real (if soft) elastic modulus so
that the whole domain -- solid and fluid alike -- solves as a genuine
two-phase elastic continuum, which is what makes this a non-trivial
Kolosov-Muskhelishvili inhomogeneity problem in the first place. It is a
soft ersatz value, not HTPB's real modulus, chosen purely so MLMG converges
at all: AP's real modulus is ~1,894x (mu) / ~250x (kappa) stiffer, and even
this contrast requires capping multigrid coarsening
(``elastic.max_coarsening_level``, see "Grid resolution" below) to avoid
divergence through the diffuse interface.

**Loading and boundary conditions.** All physical loading enters through
``elastic.apply_fluid_pressure``'s interfacial body force at the AP/Gas
diffuse interface (see "How the loading works" below) -- there is no
separately hand-derived domain-boundary traction. The domain's outer
boundary carries **rollers** on all four faces (normal displacement
clamped, tangential traction-free; the corners get both components
clamped), the same convention ``input.lm.ap_htpb_packed_elastic`` already
uses on its x-faces, applied symmetrically here since this problem has no
privileged direction. Rollers, not a fully traction-free boundary, are
required for a subtle reason: an all-traction boundary is well-posed only
up to rigid-body motion (2 translations + 1 rotation in 2D). The
interfacial load is self-equilibrated (net force/torque zero by symmetry),
so a fully free boundary wouldn't make MLMG blow up outright, but the
discrete system is singular/rank-deficient, and empirically the "smoother"
bottom solver handles that badly: with only one face clamped instead of
rollers on all four, the residual bottomed out around ``3e-6`` by iteration
~100, then crept back up geometrically, diverging by iteration ~10,000.
Rollers distribute a much weaker (purely normal) constraint along the
*entire* boundary instead, which converges cleanly (161 iterations, ~1.6s).
Since the domain half-width is 10x the inclusion radius, the correction the
rollers introduce relative to a "fully free at infinity" solution is
``O((a/L)^2) ~ 1e-4`` relative -- negligible next to the ~0.3-2% error this
test actually measures (see "Comparison").

How the loading works
----------------------

``elastic.apply_fluid_pressure`` adds RHS ``-P*grad(eta)`` to
``div(sigma) = rhs`` (``eta`` = the AP/Gas phase field, 1 in AP / 0 in Gas;
``elastic.traction = 0`` so this pressure term is the only source). Since
``P`` is spatially uniform, ``-P*grad(eta) = -grad(P*eta)`` exactly, so
``div(sigma + P*eta*I) = 0``: the combination ``sigma + P*eta*I`` is
continuous across the interface, which means ``sigma`` itself jumps there::

    sigma_rr(a+) = sigma_rr(a-) + P

This is **not** ordinary matched-material traction continuity -- it is the
correct statement of a fluid at pressure P pushing on a solid surface (a
delta-function body force localized at the, here diffuse, interface). This
was confirmed directly against raw simulation output during test
development: assuming plain continuity under-predicts the interior stress
magnitude by roughly a factor of 2.

Analytic solution
------------------

For purely radial (no angular dependence) loading, a circular inhomogeneity
problem has no angular (``cos 2*theta``, ``sin 2*theta``, ...) harmonic
content -- the general Kolosov-Muskhelishvili complex potentials
``phi(z) = A*z + B/z``, ``psi(z) = C*z + D/z`` collapse to their elementary
``n = +-1`` (monopole/dipole) terms, which is exactly the axisymmetric
Lame solution ``u_r(r) = A*r + B/r``, with:

- ``B = 0`` required for the (bounded) inclusion, and
- ``A = 0`` required in the matrix, since the domain boundary is (up to
  the negligible roller correction above) traction-free at large r: only
  the localized ``B/r`` dipole term can decay to zero stress as
  ``r -> infinity``.

Plane-strain constitutive law (matching
``Model::Solid::Finite::NeoHookeanPredeformed`` linearized about ``F = I``
-- see ``tests/LowMachElasticPressure/Readme.rst``): with
``lambda = kappa - 2*mu/3``,

    sigma_rr = 2*(lambda+mu)*A - 2*mu*B/r^2
    sigma_tt = 2*(lambda+mu)*A + 2*mu*B/r^2

Define the "2-D areal modulus" ``k = 2*(lambda+mu) = 2*kappa + (2/3)*mu``
for each phase (1 = matrix/Gas, 2 = inclusion/AP). Matching ``u_r``
(continuous) and the ``sigma_rr`` jump above at ``r = a``, with ``A1 = 0``
and ``B2 = 0``, gives::

    A2 = -P / (k2 + 2*mu1)
    B1 = a^2 * A2

Uniform stress inside the inclusion (``r < a``)::

    sigma_xx = sigma_yy = sigma_in = k2*A2

Outside (``r >= a``), pure dipole decay, no remote offset::

    sigma_rr(r) = -2*mu1*B1/r^2
    sigma_tt(r) = +2*mu1*B1/r^2

Along the ``y = 0`` sampling ray (``theta = 0`` or ``pi``),
``sigma_xx = sigma_rr``, ``sigma_yy = sigma_tt`` directly, and ``u_r`` maps
to ``disp_x = u_r * sign(x)``, ``disp_y = 0``.

With ``mu1 = 5.0e3``, ``kappa1 = 4.0e5`` (Gas/void ersatz), ``mu2 =
9.47e6``, ``kappa2 = 1.00e8`` (AP, real rescaled moduli), ``P = 1``:
``k1 = 803333``, ``k2 = 2.06313e8``, ``A2 = -4.8468e-9``,
``B1 = -1.0905e-16``, giving ``sigma_in = -0.999952`` -- i.e. this
configuration sits extremely close to the **rigid-inclusion limit**
(``k2 >> k1``): AP is so much stiffer than the void that it barely deforms
further under the interfacial load and simply carries almost exactly
``-P`` throughout, with a correspondingly tiny exterior dipole field
(``B1`` is ~16 orders of magnitude smaller than ``a^2``). See
``generate_reference.py`` for the code that produces these numbers and the
reference CSV.

Comparison
----------

``test`` samples a ray along ``y ~= 0`` (offset by half a cell -- see
`test` for why an exact ``y=0.0`` ray degenerates) from ``x = -1.3mm`` to
``x = +1.3mm`` (``disp_x``, ``disp_y``, ``stress_xx``, ``stress_yy``),
excludes a band ``|r - a| < 3*w`` around the AP/Gas interface where the
analytic sharp jump and the simulation's smoothed ``eta`` transition
necessarily disagree, and checks:

1. Interior stress: ``stress_xx``/``stress_yy`` uniform and equal to
   ``sigma_in`` well inside the inclusion (8% tolerance; achieved is
   ~0.3-0.4%).
2. Full-profile **absolute** (not relative) error of ``disp_x``,
   ``stress_xx``, ``stress_yy`` against the closed-form solution, via
   ``testlib.validate`` against the reference CSV when one is supplied by
   the run's ``check-file``, and again inline as a fallback. Absolute
   tolerances are used deliberately: because this configuration sits so
   close to the rigid-inclusion limit, the analytic exterior field decays
   to ``O(1e-5)*P`` within a few interface widths of ``r=a`` -- a
   relative-L2 metric over the whole ray would be dominated by that
   near-zero exterior region, where even the diffuse interface's own
   residual smoothing (numerically ``~1e-4*P``) reads as enormous relative
   error despite being physically negligible.
3. ``disp_y`` negligible (absolute) along the sampling ray, a consequence
   of the purely radial (no angular dependence) loading.

The ``[2d-pressure-off]`` sub-case is a wiring control:
``elastic.apply_fluid_pressure = 0`` makes the interfacial RHS identically
zero, and the boundary carries no loading of its own (rollers are
zero-traction tangentially, zero-displacement normally -- there's no
"pressure" leaking in through the boundary either way), so the solve must
return zero displacement everywhere.

LowMach elastic inclusion (LARGE-DOMAIN variant): AP in a wider HTPB block, top-loaded
=======================================================================================

**Domain-size sensitivity study.** This is a variant of
``tests/LowMachElasticInclusion`` that widens the HTPB block (and hence the
whole domain) 3x, from ``a/L = 1/10`` (inclusion radius to block half-width)
to ``a/L = 1/30``, while holding the mesh spacing ``dx`` fixed at the
baseline test's value (so ``amr.n_cell`` scales up 3x with it -- see
``input``). The analytic reference (see ``generate_reference.py``)
superposes the classical circular-inhomogeneity-in-*infinite*-matrix
solution onto the finite block's confined-compression far field; the error
this neglects is the Saint-Venant/image correction, ``O((a/L)^2)`` -- so
widening the domain 3x should shrink that neglected term by roughly ``9x``
and correspondingly reduce this test's RMS error against the baseline
test's. Everything else (materials, loading, mesh resolution, tolerances)
is unchanged from the baseline test; see that test's Readme for the full
derivation, reproduced below with this variant's numbers substituted in.

Analytic validation of the LowMach elastic fluid-structure-interaction solve
(``elastic.apply_fluid_pressure``, see ``LowMach::UpdateModel`` and
``tests/LowMachElasticPressure``) with **three** coexisting species in the
domain -- a circular AP inclusion embedded inside a square HTPB block,
itself sitting in a slightly taller rectangular domain with a Gas layer
above it that delivers the loading -- exercising the elastic solver's
species-mixing path (``LowMach::UpdateModel``'s composition-weighted model,
see below) rather than the simple two-phase (solid-in-void) setup this test
started as.

This is a simplified, single-particle descendant of
``input.lm.ap_htpb_packed_elastic`` (the production AP/HTPB packed-bed
case): one circular AP particle centered in a square HTPB block, uniform
frozen pressure, uniform 300K temperature, frozen chemistry -- no
regression, no reacting chemistry. Where that input drives the elastic
solve with the real local gas pressure field on a complex multi-particle
geometry (and has no closed-form answer to check against), this test
isolates a geometry simple enough to have one: a circular AP inhomogeneity
embedded in an "infinite" HTPB matrix, itself under a uniform confined-
compression remote load. ``AP_solid.elastic.model.mu/kappa`` and
``HTPB_solid.elastic.model.mu/kappa`` below are the real (rescaled) moduli
from that production file (AP: E ~ 21.6 GPa, nu ~ 0.14; HTPB: nu ~ 0.499, a
nearly incompressible rubbery binder).

Setup
-----

A circular ``AP_solid`` inclusion (``rigid_solid``) of radius ``a =
1.5e-4`` m sits at the origin, embedded in a square ``HTPB_solid``
(``rigid_solid``) block spanning ``x,y in [-4.5e-3,4.5e-3]`` m -- i.e. the
HTPB block fills the domain's full width, 3x wider than the baseline test's
``[-1.5e-3,1.5e-3]``. A ``Gas`` layer (``fluid``, 0.6mm thick, same
physical thickness as the baseline test) sits above the block, from
``y=4.5e-3`` to the domain top at ``y=5.1e-3``; the domain is therefore a
**slight rectangle** (``9.0mm x 9.6mm``), not the square the original
circular-inclusion-in-void version of this test used. All three species
carry real elastic moduli
(``AP_solid``, ``HTPB_solid`` genuinely; ``Gas`` a soft ersatz value, see
"Numerical stability" below) so the whole domain solves as a genuine
three-phase elastic continuum.

Chemistry is frozen and the temperature field is a uniform 300K everywhere,
so none of the AP/HTPB or HTPB/Gas interfaces ever move.

**Loading and boundary conditions.** All physical loading enters through
``elastic.apply_fluid_pressure``'s interfacial body force at the HTPB/Gas
diffuse interface -- Gas's uniform, time-frozen pressure ``P``
(``projection.update_pressure = 0``, see
``tests/LowMachElasticPressure/Readme.rst``) pushes straight down on the
top face of the HTPB block, the same "fluid pressure IS the interfacial
traction" mechanism as the original version of this test (see "How the
loading works" below), just now delivered onto a flat solid/fluid boundary
instead of a circular one. AP/HTPB is an **ordinary bonded interface**
(plain continuity of traction and displacement) -- the special fluid-
pressure jump only applies where an actual fluid (Gas) touches a solid.

The domain's side and bottom faces carry **rollers** (normal displacement
clamped, tangential traction-free) -- these coincide exactly with the HTPB
block's own side/bottom faces, since the block spans the full domain width
and its bottom sits at the domain bottom. The top face (Gas's own free
surface) is **fully traction-free**, not a roller: nothing supports the
domain from above, and the only load anywhere is the interfacial pressure
term. Unlike the original version of this test (which needed rollers on
*every* face purely to remove rigid-body null modes from a self-
equilibrated interfacial load), this system is **statically determinate**:
the net downward force from the pressure interface is genuinely reacted by
the bottom roller, not merely constrained against drift.

How the loading works
----------------------

Unchanged from the original version of this test (see
``tests/LowMachElasticPressure/Readme.rst`` and ``LowMach::UpdateModel``):
``elastic.apply_fluid_pressure`` adds RHS ``-P*grad(eta)`` to
``div(sigma) = rhs``, where ``eta`` is the total rigid-solid phase fraction
(1 in AP or HTPB, 0 in Gas). Since ``P`` is spatially uniform,
``div(sigma + P*eta*I) = 0``, so ``sigma`` jumps across the diffuse
interface: ``sigma_nn(fluid side) = sigma_nn(solid side) + P`` (crossing in
the direction of the outward normal from solid into fluid). Applied at the
flat HTPB/Gas interface (outward normal ``+y``), with the Gas layer itself
carrying zero stress (see next section):
``0 = sigma_yy(HTPB top) + P``, i.e. ``sigma_yy(HTPB top) = -P``.

Step 1: far-field confined-compression state
----------------------------------------------

Ignoring the AP inclusion, the block/roller/top-load system is **exactly**
1-D: the HTPB block spans the full domain width with rollers on both sides
(``u_x = 0`` at ``x = +-4.5e-3`` for every ``y``), so by symmetry
``u_x = 0`` and all fields are uniform in ``x`` -- this is confined
("oedometer") compression, not free uniaxial stress: ``eps_xx = eps_zz =
0`` (plane strain). This is an **exact** solution of the no-inclusion
problem (not merely a far-field approximation), since nothing breaks the
x-uniformity anywhere in the block or the Gas layer above it. The Gas
layer, by the same 1-D argument (traction-free top + roller sides + a
uniform interfacial load at its bottom), is itself statically determinate
with ``sigma_yy = 0`` throughout, **independent of Gas's own modulus**
(this is why the "Numerical stability" ersatz value below does not need to
match anything physical -- see also "Comparison").

Plane-strain constitutive law (``lambda = kappa - 2*mu/3``) then gives,
uniformly through the HTPB block::

    sigma_yy_inf = -P
    sigma_xx_inf = -P * lambda1 / (lambda1 + 2*mu1)

-- **not** equibiaxial (that only holds in the incompressible limit
``lambda1 -> infinity``; HTPB's ``nu ~ 0.499`` puts this test close to, but
not exactly at, that limit: with the real moduli below,
``sigma_xx_inf ~= -0.9960*P``, ``sigma_yy_inf = -P``).

Step 2: circular-inhomogeneity correction near the AP inclusion
-------------------------------------------------------------------

Superpose the classical two-phase circular-inhomogeneity-under-remote-
stress solution (AP embedded in "infinite" HTPB -- valid since the block is
30x the inclusion radius on every side here (3x wider than the baseline
test's 10x), so Saint-Venant/image corrections are ``O((a/L)^2)`` --
~9x smaller here than in the baseline test, negligible next to the ~1-6%
error the baseline test tolerates). Decompose the remote state
``(sigma_xx_inf, sigma_yy_inf)`` into:

- an **isotropic** part ``p0 = (sigma_xx_inf+sigma_yy_inf)/2``, handled by
  the same axisymmetric Lame solution (``u_r = A*r + B/r``) the original
  version of this test used, but now sourced by a REMOTE stress rather than
  an interfacial jump, so ``A1 != 0`` in the matrix::

      A1 = p0/k1
      A2 = p0*(k1+2*mu1) / (k1*(k2+2*mu1))
      B1 = a^2*(A2 - p0/k1)
      sigma_in_iso = k2*A2                                              (r<a)
      sigma_rr(r) = p0 - 2*mu1*B1/r^2, sigma_tt(r) = p0 + 2*mu1*B1/r^2   (r>=a)

  where ``k = 2*kappa + (2/3)*mu`` (1 = HTPB matrix, 2 = AP inclusion),
  same definition as before.

- a **deviatoric** part ``s = (sigma_xx_inf-sigma_yy_inf)/2`` (equivalent
  to a remote pure-shear state at 45 degrees), handled by the classical
  circular-inhomogeneity-under-remote-shear solution via Kolosov-
  Muskhelishvili complex potentials. Solving the bonded-interface matching
  problem gives, with ``kM1 = 3-4*nu1`` (Muskhelishvili's plane-strain
  material constant, for the **matrix only** -- remarkably, for a circular
  inhomogeneity under remote shear, the inclusion's own ``kM2`` drops out
  of the solution entirely)::

      B  = s*a^2*(mu1-mu2) / (mu1 + kM1*mu2)
      D3 = a^2*B
      gamma2p = -s*mu2*(1+kM1) / (mu1 + kM1*mu2)

  Interior (``r<a``, **uniform** Cartesian stress -- the classic 2-D
  "Eshelby" result that a circular inhomogeneity's interior field under
  remote uniform stress is itself uniform)::

      sigma_xx_dev_in = -gamma2p,  sigma_yy_dev_in = +gamma2p

  Along the ``y=0`` ray (``theta=0/pi``), exterior (``r=|x|>=a``)::

      sigma_xx_dev(x) = -4*B/x^2 + s + 3*D3/x^4
      sigma_yy_dev(x) =        -s - 3*D3/x^4

  (even in ``x``, so this holds for ``x<0`` too without extra sign
  handling). Displacement, evaluated on the real axis
  (``uy_dev = 0`` there by symmetry, matching the isotropic part)::

      ux_dev(x) = [(kM1+1)*B/x + s*x - D3/x^3] / (2*mu1)   (r>=a, odd in x)
      ux_dev(x) = -gamma2p*x / (2*mu2)                     (r<a, odd in x)

**Total**: ``sigma_xx = sigma_xx_iso + sigma_xx_dev``, etc.; ``ux = ux_iso
+ ux_dev`` with ``ux_iso(x) = A1*x + B1/x`` (``r>=a``), ``A2*x`` (``r<a``).
See ``generate_reference.py`` for the full derivation notes and the code
that evaluates these formulas.

With ``mu1=1.67e6, kappa1=8.33e8`` (HTPB), ``mu2=9.47e6, kappa2=1.00e8``
(AP), ``a=1.5e-4``, ``P=1.0e6`` (1 MPa, physical units): ``p0 = -998001 Pa``,
``s = 1999.46 Pa`` (the deviatoric correction is small -- about 0.2% of P --
because HTPB is nearly incompressible, so the confined-compression far field
sits close to hydrostatic already), giving ``sigma_xx_in = -980674 Pa``,
``sigma_yy_in = -987464 Pa`` -- notably **not equal** to each other, unlike
the original equibiaxial-loading version of this test.

``disp_y`` along the sampling ray is **not** zero here (unlike the original
version): the bottom roller pins ``u_y=0`` at ``y=ylo``, not at ``y=0``, so
the block's own rigid vertical compaction contributes a uniform offset
``u_y_base(y) = eps_yy*(y-ylo)``, ``eps_yy = sigma_yy_inf/(lambda1+2*mu1)``,
on top of the (still-zero, by symmetry) inclusion perturbation to
``disp_y`` along ``y~=0``. See ``test``'s disp_y check.

Numerical stability
--------------------

``elastic.void.model.mu/kappa`` (Gas) is a soft ersatz modulus, needed
purely so MLMG stays well-conditioned -- its exact value does not enter the
analytic reference at all, since the Gas layer's stress state is
statically determinate given its free top and roller sides (see Step 1
above), independent of Gas's own modulus. This value is **stiffer** than
the original circular-inclusion-in-void version of this test used
(``5.0e3/4.0e5``): with HTPB's real, nearly-incompressible ``kappa=8.33e8``
now a first-class matrix material (rather than AP alone sitting in void),
the HTPB/void kappa contrast at the old ersatz value is ~2,082x, well above
the ~[250x,2100x] window ``input.lm.ap_htpb_packed_elastic`` verifies
stable (see its WARNING) -- confirmed directly: at
``elastic.max_coarsening_level=3`` MLMG genuinely diverges (blows up within
1-2 V-cycles) at the old ersatz value, but converges cleanly, if very
slowly, at ``elastic.max_coarsening_level=0`` (no coarsening at all),
which rules out a physically singular setup and points squarely at a
multigrid-coarsening instability specific to this contrast. Raising the
ersatz Gas modulus 10x (``5.0e4/4.0e6``) brings the HTPB/void kappa
contrast down to ~208x, back in the verified window, and
``elastic.max_coarsening_level=3`` converges cleanly again.

Comparison
----------

``test`` samples a ray along ``y ~= 0`` (offset by half a cell -- see
``test`` for why an exact ``y=0.0`` ray degenerates) from ``x = -1.3mm`` to
``x = +1.3mm`` (``disp_x``, ``disp_y``, ``stress_xx``, ``stress_yy``),
excludes a band ``|r - a| < 3*w`` around the AP/HTPB interface where the
analytic sharp jump and the simulation's smoothed ``eta`` transition
necessarily disagree, and checks:

1. Interior stress: ``stress_xx`` uniform and equal to ``sigma_xx_in``,
   ``stress_yy`` uniform and equal to ``sigma_yy_in``, well inside the
   inclusion (8% tolerance on each; baseline test achieves ~1.1-1.4% at
   ``a/L=1/10`` -- see below for this variant's achieved error at
   ``a/L=1/30``).
2. Full-profile **absolute** RMS error of ``disp_x``, ``stress_xx``,
   ``stress_yy`` against the closed-form solution (interior + exterior
   branches, interface band excluded), via ``testlib.validate`` against the
   reference CSV when one is supplied by the run's ``check-file``, and
   again inline as a fallback (baseline test achieves disp_x ~6.4e-8 vs
   5e-7 tolerance, stress ~2.5e3-4.1e3 Pa vs 6e4 Pa tolerance at
   ``a/L=1/10``; this is the error this domain-widening variant is checking
   for improvement in).
3. ``disp_y`` matches the uniform rigid compaction offset
   ``eps_yy*(y_ray-ylo)`` along the sampling ray (absolute tolerance) --
   **not** zero, since (unlike the original version) the bottom roller
   pins ``u_y=0`` at the domain bottom, not at ``y=0`` -- see "Step 2"
   above.

**Achieved error at a/L=1/30 vs. the baseline's a/L=1/10** (both runs at
``elastic.tol_rel=1e-8``, this variant needing
``elastic.max_coarsening_level=3``/``elastic.solver.max_iter=20000`` -- see
"MLMG coarsening instability" above): the ``disp_y`` relative error (the
motivation for this variant -- see "Step 2") drops by roughly 4-5x when the
domain is widened from 10x to 30x the inclusion radius, confirming that the
baseline's ``disp_y`` mismatch is a finite-domain confinement effect
(coupling through HTPB's near-incompressible nu~0.499, not a discretization
or solver defect -- an inclusion-induced lateral strain that can't relax
against the roller side walls gets converted into extra vertical
compaction, integrated cumulatively over the column height into disp_y):

================================  ============  ==============
metric                            baseline      LargeDomain
                                   (a/L=1/10)    (a/L=1/30)
================================  ============  ==============
disp_y relative error, far-field  3.7% (max)    0.87% (max)
                                   3.7% (mean)   0.80% (mean)
disp_y relative error, overall    5.4% (max)    1.4% (max)
                                   3.9% (mean)   0.87% (mean)
stress_xx relative error (mean)   0.23%         0.25%
stress_yy relative error (mean)   0.38%         0.36%
stress_xy relative error (mean)   0.005%        0.006%
================================  ============  ==============

(relative errors as a percentage of the applied pressure ``P`` for
stresses, and of the exact rigid-offset value for ``disp_y``; interface
band excluded as above.) Stress errors are essentially unchanged between
the two domain sizes, as expected -- stress equilibrates locally, so it
was never sensitive to the confinement effect that ``disp_y`` (a
column-integrated quantity) picks up.

The ``[2d-pressure-off]`` sub-case is a wiring control:
``elastic.apply_fluid_pressure = 0`` makes the interfacial RHS identically
zero, and none of the boundaries carry any loading of their own (rollers
are zero-traction tangentially/zero-displacement normally; the top is
simply traction-free), so the solve must return zero displacement
everywhere.

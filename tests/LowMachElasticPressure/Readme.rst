LowMach elastic fluid-structure interaction driven by fluid pressure
======================================================================

Analytic validation of ``elastic.apply_fluid_pressure``, the option added to
``LowMach::UpdateModel`` that drives the elastic interfacial body force from
the local fluid pressure field::

    rhs = elastic.traction * grad(rigid_eta) - pressure * grad(rigid_eta)

Unlike ``tests/LowMachElastic`` (an explicit smoke/regression test -- see
that test's Readme.rst) and ``tests/LowMachMovingBoundary``, this test
compares the solver's output against a closed-form elasticity solution, so
it pins down the *sign* and *magnitude* of the fluid-pressure coupling, not
just that it produces finite output.

Setup
-----

A single ``rigid_solid`` species (``Solid``) occupies the bottom half of the
domain (``y < y0 = 3.2e-3``); a single fluid species (``Gas``) occupies the
top half at a uniform pressure ``P``. Only one rigid solid species is used
(rather than the two-species AP/HTPB split in ``LowMachElastic``) so the
elastic model field is spatially uniform -- the arbitrary-species mixing is
already covered by that test, and a uniform model sidesteps its documented
``Operator::Elastic`` heterogeneous-media MLMG-divergence gap.

Elastic boundary conditions are **rollers** on the x-faces
(``elastic.bc.constant.type.xlo/xhi = disp trac`` -- displacement in x,
traction in y) and clamped at the base (``ylo = disp disp``), with the top
(fluid side) traction-free. This forces ``u_x = 0`` everywhere, reducing the
problem to exact 1-D confined (oedometric) compression.

Freezing the pressure field
----------------------------

LowMach's pressure projection is non-incremental: every step it recomputes
``pressure_mf = max(pressure_reference + phi/pressure_scale, pressure_floor)``,
where ``phi`` is the projection potential -- a full overwrite of whatever the
IC prescribed. ``projection.update_pressure = 0`` disables only that final
write; the projection (and the rigid-solid velocity relaxation it contains)
still runs in full, but ``pressure_mf`` retains the constant IC value for the
entire run. ``projection.enabled = 0`` is not an option here: LowMach aborts
at parse time when rigid solid species are present with projection disabled.

Analytic solution
------------------

``Model::Solid::Finite::NeoHookeanPredeformed`` is a decoupled
(deviatoric/volumetric-split) compressible neo-Hookean model. Expanding its
strain energy to second order about ``F = I`` gives the standard linear
plane-strain form (2D is plane strain here -- ``F33`` is pinned to 1, and
there is no plane-stress option for the finite models)::

    mu_eff     = mu
    lambda_eff = kappa - 2*mu/3
    K          = kappa                       (true 3D bulk modulus)
    M          = lambda_eff + 2*mu           (confined/oedometric modulus)
              = kappa + 4*mu/3

With ``mu = 140``, ``kappa = 150``: ``lambda = 56.667``, ``M = 336.667``.

The elastic solve is a single-iteration Newton step
(``elastic.solver.nriters = 1``) from a zeroed initial guess
(``elastic.zero_out_displacement = 1``), which linearizes exactly about
``F = I`` -- so this is not an approximate small-strain limit, it is the
exact discrete target for this configuration (the only error sources are
the diffuse solid/fluid interface and the mesh).

For a solid half-space ``y >= 0`` confined at the base and loaded by a
uniform pressure ``P`` at ``y = y0`` with ``u_x = 0``, force balance
``d(sigma_yy)/dy = 0`` and the constitutive law give::

    sigma_yy(y) = -P                          (uniform in the solid)
    disp_y(y)   = -P*y/M                      (linear, disp_y(0) = 0)
    disp_x(y)   = 0
    sigma_xx    = -P * lambda/M               (= -0.1683*P here)

**Sign check.** The elastic operator solves ``div(sigma) = rhs``. Since
``rhs = -P * grad(rigid_eta)`` and ``rigid_eta`` goes from 1 (solid) to 0
(fluid) across the diffuse interface, integrating the balance across the
band gives ``[sigma_nn] = +P`` (i.e. ``sigma_nn`` jumps up by ``P`` moving
from fluid to solid), so ``sigma_nn = -P`` in the solid -- compression, as
expected for a compressive fluid pressure pushing on the solid.

Why the pressure is O(1) Pa
----------------------------

As in ``tests/LowMachElastic`` (see that input's comments), the elastic
solver's multigrid V-cycle is only well-behaved when the modulus/``dx^2``
ratio stays moderate; literal SI moduli (~1e8 Pa) combined with this test's
SI grid spacing (~1e-4 m) push that ratio far out of the solver's converged
regime. Keeping ``mu``/``kappa`` at the same nondimensional-style magnitudes
used by the other LowMach elastic tests, and scaling the applied pressure
down to O(1) Pa accordingly, keeps strains at a physically sane ~3e-3
without touching the (separately validated) solver conditioning question.
This is a mechanics validation, not a simulation at realistic chamber
pressure.

What the test checks
---------------------

Three sub-cases, distinguished by ``elastic.apply_fluid_pressure`` and the
applied pressure ``P`` (1 Pa and 4 Pa):

- ``2d-pressure-off``: the fluid-pressure branch disabled and
  ``elastic.traction = 0`` -- the interfacial RHS is identically zero, so
  displacement must be (numerically) zero. Control case.
- ``2d-pressure-1Pa`` / ``2d-pressure-4Pa``: sampling a vertical ray through
  the domain interior (``x = 1.6e-3``), restricted to ``y in [0.5e-3, 2.4e-3]``
  -- away from the clamped base and more than 10 interface widths below the
  diffuse solid/fluid band at ``y0 = 3.2e-3`` -- the test checks:

  1. ``disp``/``stress`` are finite everywhere.
  2. ``sigma_yy ~= -P`` (asserted within 2%; achieved ~0.2-0.8%) -- the
     headline sign/magnitude check.
  3. ``disp_y(y)`` matches ``-P*y/M`` (asserted within 1% relative L2 error;
     achieved ~1e-9, confirming this is genuinely the exact linear solve
     described above, not an approximation).
  4. ``disp_x`` is negligible relative to ``disp_y`` (confinement holds).
  5. ``sigma_xx/sigma_yy ~= lambda/M`` (asserted within 5%; achieved
     ~1-4%, the loosest of the checks -- likely residual diffuse-interface/
     discretization effects on the shear-coupled component).

Boundary-node caveat
---------------------

``Base::Mechanics`` overwrites the interfacial RHS at every domain-boundary
node with the prescribed elastic BC value (this is how ``elastic.bc.type``
is enforced). The analytic comparison therefore samples the domain interior
only, well clear of the ``x = 0``/``x = 3.2e-3`` boundary columns.

Run the case with::

    ./configure --dim=2
    make -j bin/lowmach
    ./scripts/runtests.py tests/LowMachElasticPressure --serial --dim=2

Known follow-on
----------------

This remains one-way coupling (fluid pressure -> solid traction only); the
resulting displacement/stress is not fed back into the LowMach momentum
equation. A curved-interface (Lame pressurized-cavity) validation and a
mesh-convergence-rate study were both considered and set aside for now.

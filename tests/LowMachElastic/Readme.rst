LowMach elastic fluid-structure interaction (smoke test)
=========================================================

Smoke test for the elastic solve added to ``Integrator::LowMach`` to support
fluid-structure interaction. ``LowMach``'s parse system already supports an
arbitrary number of ``rigid_solid`` species; this feature gives each such
species elastic constants (``<species>.elastic.model.mu/kappa``), mixes them
into a nodal model field by their volume fraction (``rigid_species_eta_mf``),
and drives ``Base::Mechanics``'s existing Newton/``Operator::Elastic`` solve
from it, following the pattern already used by ``Integrator::Flame``.

This test is **not** a comparison against an analytic solution -- it is a
regression/smoke check that the pipeline (arbitrary-species mixing, ghost
cells, the elastic solve itself, and a simple constant-traction interfacial
body force) runs to completion and produces physically sane results.

Setup
-----

Two ``rigid_solid`` species (``AP_solid``, ``HTPB_solid``) occupy the bottom
half of the domain (``y < y0``), split by a diffuse interface at ``x = x0``
(AP to the left, HTPB to the right). A single fluid species (``Gas``)
occupies the top half at a uniform reference pressure. Chemistry is frozen
and there are no phase change mechanisms, so the composition -- and hence
``rigid_eta`` and the AP/HTPB split -- never evolves; this isolates the
elastic solve and model-mixing from the rest of LowMach's physics.

The two sub-cases (``2d-traction-zero``, ``2d-traction-applied``) only
differ in ``elastic.traction``, the constant interfacial body force
``rhs = elastic.traction * grad(rigid_eta)`` applied at the solid/fluid
interface (the same form used by ``Flame::UpdateModel``). Replacing this
constant with the actual local fluid pressure, and feeding the resulting
displacement/stress back into the momentum equation, is the natural next
step and is explicitly out of scope here.

The ``test`` script checks:

- ``model_mu``/``model_kappa`` form a convex combination of the two
  species' constants everywhere (proving the arbitrary-species mixing is
  correct), and that both species' values are actually reached (proving
  mixing occurred, not just one species winning everywhere).
- ``disp``/``stress``/``strain`` are finite everywhere.
- The zero-traction case produces (numerically) zero displacement.
- The traction-applied case produces a nonzero displacement response.

Numerical caveats found during development
-------------------------------------------

Two pre-existing LowMach/solver behaviors had to be worked around to get a
stable case, neither of which is specific to the new elastic code:

1. **Gas density IC must be volume-fraction-consistent.** LowMach's pressure
   projection includes a mixture-volume-constraint term (the ``mixed_phase``
   branch of ``RHS()``: ``rhs += (volume_fraction - 1) / dt``) that requires
   ``gas_volume_fraction + rigid_eta == 1`` everywhere, not just far from the
   interface. Choosing the gas density profile independently of ``rigid_eta``
   (e.g. a spatially uniform gas density) leaves a residual that gets divided
   by ``dt``, injecting an ``O(1/dt)`` pressure-Poisson source that produces
   unbounded velocity as ``dt -> 0`` -- indistinguishable from a solver
   instability unless you know to look for it. The fix (matching the
   ``Final`` species IC in ``tests/LMRFMonoAP/input``) is to define the gas
   density as the exact ideal-gas-consistent complement of the solid volume
   fraction: ``rho_gas = (P/(R*T)) * (1 - rigid_eta(x,y))``.

2. **The elastic multigrid solve is sensitive to material contrast at this
   grid scale.** With literal SI moduli (~1e8 Pa) at this test's SI grid
   spacing (~1e-4 m), or with a large modulus contrast between the two rigid
   solid species (e.g. AP/HTPB values of 140/8, as used in
   ``tests/SCPSpheresElastic``), the ``MLMG`` V-cycle for the elastic solve
   diverges on this mesh -- independent of the void fraction, the traction
   magnitude, x-periodicity, or the diffuse-interface profile shape used for
   the species split. A spatially *uniform* model (no contrast at all)
   converges reliably in a few tens of iterations regardless of magnitude;
   the mixed model field itself was verified correct by direct inspection
   (finite, correctly weighted, no NaN/Inf) in every case, including the
   diverging ones -- so this looks like a pre-existing robustness gap in
   ``Operator::Elastic``'s multigrid for heterogeneous media rather than a
   bug in the new mixing/model code. This test therefore uses non-SI-scale
   moduli (matching the convention used by every other elastic test in this
   suite) and a deliberately modest AP/HTPB contrast (~7%). Reproducing the
   divergence with a genuinely realistic (e.g. 140/8) contrast, and either
   fixing or better characterizing it, is worth a follow-up investigation.

Run the case with::

    ./configure --dim=2
    make -j bin/lowmach
    ./scripts/runtests.py tests/LowMachElastic --serial --dim=2

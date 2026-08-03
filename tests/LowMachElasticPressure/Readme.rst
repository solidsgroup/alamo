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

Four sub-cases, distinguished by ``elastic.apply_fluid_pressure``, the
applied pressure ``P`` (1 Pa and 4 Pa), and whether a soft void model is set:

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
- ``2d-pressure-1Pa-softvoid``: repeats ``2d-pressure-1Pa`` with
  ``elastic.void.model.mu/kappa`` set to 1e-4 of the solid's (which also
  disables the ``psi`` mask -- see ``LowMach::Parse``) and
  ``elastic.max_coarsening_level=2``, checked against the *same* reference
  CSV as ``2d-pressure-1Pa``. This exercises the soft "void modulus" path
  added for ``input.lm.ap_htpb_packed_elastic``, where masking the void
  purely via ``psi``'s ``1e-8`` floor reintroduced the full solid modulus
  there and made the MLMG solve diverge.

  The coarsening cap is not optional here: at this test's grid spacing and
  diffuse interface width (~2 cells), MLMG genuinely diverges -- not just
  converges slowly -- once the solid/void contrast passes roughly ``1e2``
  to ``1e3``, unless coarsening is capped (coarsening through a couple of
  cells spanning a factor-``1e4`` coefficient jump produces garbage
  coarse-level operators). This is the same heterogeneous-media
  MLMG-divergence gap noted above, and is why
  ``input.lm.ap_htpb_packed_elastic`` also sets
  ``elastic.max_coarsening_level``. Aside from the coarsening cap, the void
  modulus here is soft enough to leave the solid's analytic solution
  unaffected, confirming
  that path reduces to the same answer as the default psi-masked one.

Void modulus vs. solver stability and accuracy
------------------------------------------------

``2d-pressure-1Pa-softvoid`` above pins down one contrast (``1e4``, capped)
against the analytic solution, but doesn't show how that tradeoff moves as
the void modulus changes. To characterize it, the 1 Pa case was re-run
directly (not through ``runtests.py``) at eleven solid/void contrasts from
``1x`` to ``1e5x`` (``elastic.void.model.mu/kappa`` scaled proportionally
to the solid's 140/150), each both with and without
``elastic.max_coarsening_level=2``, and compared against the same
``disp_y = -P*y/M`` / ``sigma_yy = -P`` analytic targets used by ``test``:

.. list-table::
   :header-rows: 1

   * - contrast
     - no cap
     - cap = 2
     - disp_y err
     - sigma_yy err
   * - 1x
     - OK (17 it.)
     - OK (12 it.)
     - 0.0003%
     - 0.19%
   * - 3x
     - OK (16 it.)
     - OK (13 it.)
     - 0.53%
     - 0.72%
   * - 10x
     - OK (16 it.)
     - OK (14 it.)
     - 0.90%
     - 1.10%
   * - 33x
     - OK (17 it.)
     - OK (14 it.)
     - 0.83%
     - 1.02%
   * - 100x
     - OK (27 it.)
     - OK (14 it.)
     - 0.61%
     - 0.81%
   * - 333x
     - OK (24 it.)
     - OK (14 it.)
     - 0.44%
     - 0.63%
   * - 1,000x
     - **diverges**
     - OK (15 it.)
     - 0.36%
     - 0.56%
   * - 3,333x
     - **diverges**
     - OK (15 it.)
     - 0.33%
     - 0.52%
   * - 10,000x
     - **diverges**
     - OK (18 it.)
     - 0.32%
     - 0.51%
   * - 33,333x
     - **diverges**
     - OK (52 it.)
     - 0.31%
     - 0.50%
   * - 100,000x
     - **diverges**
     - **diverges**
     - n/a
     - n/a

Two things fall out of this that aren't obvious from a single data point:

- **Stability is governed by the coarsening cap, not raw contrast.**
  Uncapped, MLMG diverges outright (not just slowly) once contrast passes
  roughly 300-1,000x -- consistent with the coarsening argument in the
  sub-case description above. Capped at 2 levels, it stays stable out to
  ~30,000x with a flat iteration count (12-18), and only starts costing
  more iterations near the edge of that range (52 at 33,000x) before
  failing at 100,000x.
- **Accuracy vs. contrast is not monotonic**, and the smaller-contrast end
  is not actually the more physically faithful one. At contrast ~1 the void
  is essentially as stiff as the solid, so there is no real interface to
  get wrong and the ~1e-3% error is not evidence the interface physics is
  being captured -- it is not really modeling a void at all. Error is
  *worst* in the middle, around contrast 3-30 (~0.9-1.1%), then improves
  monotonically as contrast grows, settling near a ~0.3%/0.5% floor by
  ~1,000-3,000x. That floor is set by the diffuse-interface discretization
  itself (dx and interface width -- see "Why the pressure is O(1) Pa"
  above), not by how soft the void is, so pushing contrast higher than
  that buys no more accuracy.

Net implication: once ``elastic.max_coarsening_level`` is capped, there is
no accuracy reason to run the void softer than roughly 1,000-3,000x below
the solid -- doing so only spends iteration budget for no benefit, and
going much softer than that risks the stability cliff above.
``input.lm.ap_htpb_packed_elastic`` currently sits at ~1,670x (HTPB) to
~1e5x (AP kappa); the AP side is past where this sweep's accuracy plateau
sets in and close to where even the capped solver starts needing
meaningfully more iterations -- a candidate to revisit (e.g. raising
``elastic.void.model.mu/kappa`` toward AP-contrast ~1,000-3,000x) if that
input file's own MLMG iteration count runs high in practice.

Scripted follow-up: coarsening depth and resolution
-----------------------------------------------------

The sweep above was run by hand and never varied resolution, so it could not
show whether iteration count is h-independent -- the actual acceptance test
for any coarse-grid fix, and the one axis ``elastic.max_coarsening_level``
capping can never recover on its own. ``scripts/solver_benchmark.py`` (added
alongside ``scripts/solverlib.py`` for the MLMG diagnostics-first pass, see
``~/.claude/plans/improve-MLMG-solver-ideas.md``) scripts this case by
invoking ``bin/lowmach-<dim>d-<comp>`` directly with the same
``elastic.void.model.mu/kappa`` scaling used above. Its ``--preset contrast``
reproduces every cell of the table above (iteration counts and errors both),
confirming the harness before trusting its other presets. One cell needed
manual attribution: at cap=2, 33,333x sits close enough to the stability edge
that a ~1e-5 relative change in the void modulus (e.g. rounding ``140/33333``
to 4 vs. 17 significant figures) flips it between converging in ~51
iterations and diverging in 2 -- a real property of the operator there, not
a harness bug.

Note: this test does not set ``elastic.solver.conservative_face_flux = 1``
(unlike ``input.lm.ap_htpb_packed_elastic``), so everything below runs on the
*non-conservative* ``Fapply`` path, not the conservative flux-difference one.

**A2 -- bisecting the coarsening depth** (``--preset cap-bisect``), at three
contrasts, ``amr.n_cell = 32 64``:

.. list-table::
   :header-rows: 1

   * - contrast
     - cap=0
     - cap=1
     - cap=2
     - cap=3
     - cap=4
     - cap=5
     - cap=6
     - uncapped
   * - 1,000x
     - OK (3 it.)
     - OK (11 it.)
     - OK (15 it.)
     - OK (17 it.)
     - **diverges**
     - **diverges**
     - **diverges**
     - **diverges**
   * - 10,000x
     - OK (3 it.)
     - OK (11 it.)
     - OK (18 it.)
     - **diverges**
     - **diverges**
     - **diverges**
     - **diverges**
     - **diverges**
   * - 33,333x
     - OK (2 it.)
     - OK (11 it.)
     - **diverges**
     - **diverges**
     - **diverges**
     - **diverges**
     - **diverges**
     - **diverges**

Two things fall out of this that the hand-run sweep couldn't show:

- **``max_coarsening_level = 0`` always converges, in 2-3 iterations, at
  every contrast tested -- including 33,333x.** With zero coarsening the
  fine-level operator is well-behaved regardless of contrast; every failure
  mode in this table is introduced by coarsening, not present without it.
  Read against the plan's root-cause table, this rules the fine-level
  averaging story (R8) *out* as the dominant issue for this
  (non-conservative-path) configuration and points at the coarse-grid
  coefficient averaging (R1) instead.
- **The highest coarsening depth that still converges falls as contrast
  rises** (3 -> 2 -> 1, for 1,000x -> 10,000x -> 33,333x), rather than
  cliffing at the same fixed depth regardless of contrast. That is more
  consistent with coefficient-averaging error accumulating per level and
  crossing MLMG's divergence threshold sooner at higher contrast (R1) than
  with a fixed structural/layout defect (R6) that would bite at a constant
  depth. A4 (offline two-grid spectral-radius analysis, deferred from this
  pass) is the instrument that can confirm this directly.
- Accuracy is flat across every coarsening depth that converges (e.g.
  ``disp_y_err`` at 1,000x is 0.0036 at cap=0 through cap=3, to 4
  significant figures) -- confirming again that ``max_coarsening_level``
  trades stability, not accuracy.

**H-refine -- the missing axis** (``--preset h-refine``), ``amr.n_cell`` in
{32x64, 64x128, 128x256}, at cap=2 and uncapped:

At 32x64 and 64x128 the qualitative picture above holds (uncapped diverges,
cap=2 converges except at the 33,333x edge case noted above), with iteration
counts if anything falling slightly with resolution rather than growing --
no sign of the iteration count blowing up under refinement in the range
tested. **128x256 could not be evaluated: every configuration at that
resolution crashes** with a bus error (invalid address alignment) inside
``Operator::Elastic::averageDownCoeffsSameAmrLevel``, not a convergence
failure. ``amr.max_grid_size = 64`` (set in this test's ``input``) splits a
128x256 domain into an 8-box (2x4) layout, vs. 2 boxes (64x128, split along
one direction only) or a single box (32x64). Bisecting by hand (outside the
harness): ``amr.n_cell = 128 64`` (2 boxes along x, 1 along y) crashes the
same way; ``amr.n_cell = 64 256`` (1 box along x, 4 along y) does not --
and the crash reproduces even at contrast ~1x, so it is unrelated to modulus
contrast entirely. This looks like a real, previously-undiscovered defect
specific to having more than one box along the x-direction in
coefficient coarsening, distinct from the convergence/divergence behavior
documented everywhere else on this page, and worth its own investigation
before drawing any h-independence conclusion past 64x128. Not fixed here
(this page documents the operator's numerics as found, per the
diagnostics-first pass's own
scope; see the plan for where this would fit against D1b/R1).

Boundary-node caveat
---------------------

``Base::Mechanics`` overwrites the interfacial RHS at every domain-boundary
node with the prescribed elastic BC value (this is how ``elastic.bc.type``
is enforced). The analytic comparison therefore samples the domain interior
only, well clear of the ``x = 0``/``x = 3.2e-3`` boundary columns.

Reference data
---------------

``disp_y`` and ``stress_yy`` are checked against ``reference/pressure-{off,1Pa,4Pa}.csv``
via ``testlib.validate`` (the run's ``check-file``, wired up per sub-case in
``input``), rather than an inline formula in ``test`` -- consistent with
``tests/SCPThermalContact`` and ``tests/LowMachConduction``. Each reference
CSV is produced directly from the closed-form solution above by
``generate_reference.py`` (not from a prior simulation run), so the
comparison is against a known-correct answer, not just a snapshot of past
behavior. ``mu``/``kappa`` in that script must match ``input`` exactly; if
either changes, regenerate::

    cd tests/LowMachElasticPressure
    python3 generate_reference.py

``disp_x`` (confinement) and ``sigma_xx/sigma_yy`` (constitutive ratio) are
still checked inline in ``test`` against their analytic targets, since a
zero-valued or ratio-based reference doesn't fit the ``testlib.validate``
CSV-comparison idiom cleanly.

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

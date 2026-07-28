LowMach conduction (tensorial diffuse-interface heat conduction)
==================================================================

Reproduces the 1D semi-infinite two-phase contact problem of Ettrich et al.
2014, *Modelling Simul. Mater. Sci. Eng.* 22 085006, section 8.1 (eqs
40-42): two semi-infinite bodies of different thermal conductivity and
volumetric heat capacity, initially at different uniform temperatures, are
brought into contact at ``t=0``. The interfacial temperature and the
subsequent ``erfc``/``erf`` temperature profiles have a closed-form
solution, used here as the reference.

The test uses a single fluid ("air") species and a single
``deformable_solid`` species (aluminium-like) to form the diffuse interface.
It is configured as a *pure conduction* problem so it matches the analytic
assumptions exactly:

- ``projection.enabled = 0`` -- no pressure solve, so no thermal-expansion
  velocity is ever generated; velocity stays at its zero initial condition
  for the whole run (checked directly by the ``test`` script).
- ``advect_temperature = 0`` -- temperature evolves only through the
  implicit conduction solve.
- ``chemistry.model.type = frozen`` -- no reaction, so the density field
  (and hence the interface) never moves.

This exercises ``LowMach::ComputeThermalState``'s harmonic (series)
conductivity and heat-capacity mixing rules and the tensorial mobility path
in ``Operator::Diffusion``, which is used because heat flux is normal to a
diffuse interface that is not aligned with a single grid direction only by
coincidence of this 1D setup -- the same code path applies for curved
interfaces in 2D/3D problems.

The fluid's conductivity and reference pressure are chosen to give a large
but numerically tractable conductivity/heat-capacity contrast with the
aluminium-like solid (rather than using literal air properties, which -- at
a several-cell-wide diffuse interface -- push the discretization into a
regime dominated by interface-resolution error rather than the
conductivity-mixing rule being tested).

Run the case with::

    ./configure --dim=2
    make -j bin/lowmach
    PATH=/path/to/python-with-yt/bin:$PATH \
        ./scripts/runtests.py tests/LowMachConduction --serial --dim=2

The check compares the final ``temperature`` profile along the centerline
against the analytic reference and verifies velocity remains exactly zero
throughout the run.

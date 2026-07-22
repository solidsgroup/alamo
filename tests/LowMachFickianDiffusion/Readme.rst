LowMach Fickian diffusion
=========================

This is the LowMach counterpart of
``multicomponent_master:tests/FlowDiffusion1D``.  It uses the same three-species
mixture, smoothed step initial condition, Lennard-Jones transport parameters,
timestep, final time, and analytic error-function reference profile.  The
transverse mesh is reduced because this effectively one-dimensional test does
not need the production case's resolution.

Run the case with::

    ./configure --dim=2
    make -j bin/lowmach
    PATH=/path/to/python-with-yt/bin:$PATH \
        ./scripts/runtests.py tests/LowMachFickianDiffusion --serial --dim=2

The check compares ``mole_fraction_left`` along the centerline with the hydro
test's analytic reference.  It also verifies mass-fraction closure and that the
occupancy-weighted mole-fraction sum remains unchanged.

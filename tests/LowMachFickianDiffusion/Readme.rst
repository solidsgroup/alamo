LowMach Fickian diffusion
=========================

This is the LowMach counterpart of
``multicomponent_master:tests/FlowDiffusion1D``.  It uses the same three-species
mixture, smoothed step initial condition, Lennard-Jones transport parameters,
mesh, timestep, final time, and analytic error-function reference profile.

Run the case with::

    ./configure --dim=2
    make -j bin/lowmach
    PATH=/path/to/python-with-yt/bin:$PATH \
        ./scripts/runtests.py tests/LowMachFickianDiffusion --serial --dim=2

The check compares ``mole_fraction_left`` along the centerline with the hydro
test's analytic reference and independently verifies mole-fraction closure in
the initial and final states.

LowMach mixture property consistency
====================================

These four uniform LowMach cases reproduce the mixtures in
``multicomponent_master:tests/FlowDiffusion0D``.  They compare viscosity,
thermal conductivity, and the first mixture-averaged diffusion coefficient
against that branch's hydro reference values for both constant and
Lennard-Jones transport models.

LowMach temperature initial conditions are set to the temperature reconstructed
by hydro from each case's partial densities and 1 atm pressure.  This matters
for the two species-diffusion cases, whose density data imply about 297.10 K
rather than the nominal 293 K used by the other mixture.

Build a two-dimensional executable and run only these cases with::

    ./configure --dim=2
    make -j bin/lowmach
    PATH=/path/to/python-with-yt/bin:$PATH \
        ./scripts/runtests.py tests/LowMachMixtureProperties --serial --dim=2

The checks use fields emitted by ``diagnostics.extended_fields = 1``.  The
mixtures are uniform, so both initial and final plotfiles must match the same
golden values.

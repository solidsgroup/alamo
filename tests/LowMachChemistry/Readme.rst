LowMach H2/O2 chemistry consistency
===================================

This test uses the ten-species ``h2o2.yaml`` mechanism and stoichiometric
H2/O2 mixture from ``multicomponent_master:tests/Chemistry``.  The
``2d-source`` case checks the instantaneous species production and heat
release rates at 1000 K and 1 atm against Cantera 3.2 values for that
mechanism.  The existing explicit and backward-Euler cases then exercise the
coupled LowMach chemistry, temperature, dilatation, projection, and AMR path.

Configure with YAML support and run this test family with::

    ./configure --dim=2 --yaml
    make -j bin/lowmach
    PATH=/path/to/python-with-yt/bin:$PATH \
        ./scripts/runtests.py tests/LowMachChemistry --dim=2

Use ``--serial --sections 2d-source`` for the inexpensive cell-local source
check alone.

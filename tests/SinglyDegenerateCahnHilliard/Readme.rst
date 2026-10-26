This problem tests singly degenerate mobility in :ref:`Integrator::CahnHilliard`.

Use :code:`mobility = singly_degenerate` with :code:`method = realspace` or
:code:`method = spectral` (requires a build configured with :code:`--fft`).
The default :code:`mobility_floor = 0` makes the pure phases immobile;
:code:`spectral_stabilization` defaults to :math:`\gamma(L/16+M_0)`.
Both methods use face fluxes with arithmetic face mobility and a centered
Laplacian for the chemical potential. The spectral method uses FFTs to apply
a discrete biharmonic stabilizer to the increment, including Nyquist modes.
It does not guarantee energy decrease for arbitrary timesteps.
Constant mobility remains the default. Periodic boundary conditions are required.
The realspace timestep must satisfy the explicit fourth-order stability constraint.
Neither method clips the solution by default; optional clipping can alter the
mean, so it is disabled in the conservation checks.

The 12 regression sections keep a compact set of distinct checks. In both two
and three dimensions, realspace and spectral updates are checked against the
same linearized Fourier-mode solution, which also checks their agreement and
mass conservation. A spectral case checks a nonzero mobility floor; realspace
and spectral cases keep a pure phase stationary, and a further realspace case
checks an initial value outside the pure-phase interval. Spectral checkerboards
must decay in two and three dimensions, a strongly varying mobility case must
dissipate energy at a small timestep, and a regridded spectral AMR case checks
partial fine-level coverage. The 2D spectral Fourier-mode case also uses custom
input and field names and enables diagnostics, checking the namespaced fields.
Realspace AMR does not reflux coarse/fine fluxes; use a uniform grid when exact
mass conservation is required. Spectral AMR assembles a uniform finest-level
grid, which requires memory proportional to that entire grid.

The checker reads both yt native plotfiles and Alamo HDF5 plotfiles. For example,
run the two-dimensional serial sections with:

.. code-block:: bash

    scripts/runtests.py --dim=2 --serial --fft tests/SinglyDegenerateCahnHilliard

For an HDF5-enabled build, add :code:`--hdf5` to select its matching executable.

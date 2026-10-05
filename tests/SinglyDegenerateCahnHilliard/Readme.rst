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
The realspace timestep must satisfy the explicit fourth-order stability constraint;
neither method clips the solution.

The regression sections check a linearized Fourier mode on anisotropic domains,
nonzero mobility floors, stationary phases, conservation, and energy dissipation
in two and three dimensions. Grid-scale checkerboards must decay, and a
strongly varying mobility case must dissipate energy at small timesteps. MPI
sections exercise patch/rank boundaries, and the spectral AMR section exercises
a partially refined grid with regridding.
Realspace AMR does not reflux coarse/fine fluxes; use a uniform grid when exact
mass conservation is required. Spectral AMR assembles a uniform finest-level
grid, which requires memory proportional to that entire grid.

For example, run the two-dimensional serial sections with:

.. code-block:: bash

    scripts/runtests.py --dim=2 --serial --fft tests/SinglyDegenerateCahnHilliard

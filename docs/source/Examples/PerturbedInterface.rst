.. role:: cpp(code)
   :language: c++
.. role:: bash(code)
   :language: bash

Perturbed Interface
===================

This example simluates the evolution of a perturbed interface with strongly anisotropic grain boundary energy.
For best results (currently) use 2D.
The following image is the output in 2D.

.. image:: perturbedinterface.png
   :scale: 50%
   :align: center

Notice that the equilibrium shape has a sharp corner - this is the result of the anisotropy.


Test problem
------------

Run the test problem (2D or 3D) via

.. code-block:: bash

      ./bin/alamo... tests/PerturbedInterface/input

This will generate an output directory :code:`tests/PerturbedInterface/output/`, and will re-name the old directory if it exists already.
Use VisIt to open :code:`tests/Flame/output/output.visit`.
:code:`tests/Flame/output/metadata` contains all of the input parameters.
The input file has the following form:

.. code-block:: bash

   timestep = 0.1                                 # INITIAL timestep (for isotropic)
   stop_time = 50                                 # Simulation should converge by this time
   
   plot_file = tests/PerturbedInterface/output
   
   amr.plot_int = 1000                            # Timestep decreases with anisotropy
   amr.max_level = 2                              # 2 AMR levels of refinement
   amr.n_cell = 64 64 64                          # 64x64(x64) grid
   amr.blocking_factor = 8                        # Prevent overly small patches
   amr.regrid_int = 10                            # You can probably increase this...                            
   amr.grid_eff = 1.0                            
   amr.max_grid_size = 8
   
   ic.type=perturbed_interface                    # Sets up the initial perturbation
   ic.wave_numbers=2                              # Play with this to add additional modes
   ic.wave_amplitudes=1.0                         # ... and additional amplitudes
   
   geometry.prob_lo = 0 -4 0                      # Geometry - periodic in the x (and z)
   geometry.prob_hi = 8 4 8                       # directions only
   geometry.is_periodic= 1 0 1                    #
   
   bc.hi = INT_DIR EXT_DIR INT_DIR                # Boundary conditions - Dirichlet on top
   bc.lo = INT_DIR EXT_DIR INT_DIR                # and bottom, periodic elsewhere
   bc.lo_2 = 1.0 0.0 0.0 0.0                      #
   bc.hi_2 = 0.0 1.0 0.0 0.0                      #
   
   pf.number_of_grains = 2                        # Only a top and bottom grain here
   pf.M = 1.0                                     # Mobility
   pf.mu = 10.0                                   # Numerical parameter
   pf.gamma = 1.0                                 # Numerical parameters
   pf.l_gb=0.1                                    # GB width
   pf.sigma0=0.075                                # Initial isotropic GB energy
   
   anisotropy.on=1                                # This turns on anisotropy
   anisotropy.tstart= 1.                          # Wait one second before activating
                                                  # (this gives GB time to diffuse)
   anisotropy.timestep=0.001                      # We need a smaller timestep with anisotropy
   anisotropy.theta0= 45                          # Equilibrium angles at pm 45 degrees
   anisotropy.sigma0=0.075                        # GBA parameter 1
   anisotropy.sigma1=0.07                         # GBA parameter 2
   anisotropy.beta= 0.00001                       # Curvature coefficient

See also:
---------

- :bash:`./src/alamo.cc`: Entry point for the solver
- :ref:`API-Integrator-PhaseFieldMicrostructure`: Integrator that explicitly evolves the order parameter
- :ref:`API-Operator-Elastic`: Operator to do elastic solves
- :ref:`API-IC-PerturbedInterface`: Parameterization of initial perturbation

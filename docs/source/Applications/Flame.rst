.. role:: cpp(code)
   :language: c++
.. role:: bash(code)
   :language: bash

Flame
=====

The flame simulation uses a phase field type model to simulate the burn front of a solid material undergoing burn.

Theory
------

The burned / unburned region is tracked using an order parameter :math:`\eta`.
:math:`\eta` is evolved using the kinetic equation

:math:`\frac{\partial\eta}{\partial t} = -L\frac{\delta F}{\delta \eta}`

where the free energy F is

:math:`F[\eta] = \int_\Omega f(\eta,\nabla\eta)\,dx`

and

:math:`f(\eta,\nabla\eta) = \underbrace{w(\eta)}_{\text{Chemical}} + \underbrace{\frac{1}{2}\kappa|\nabla\eta|^2}_{\text{Interface regularization}} + \underbrace{c_p\Delta T}_{\text{Thermal}} + \underbrace{\frac{1}{2}\sigma\cdot\varepsilon}_{\text{Elastic - in progress}}`


Test problem
------------

Run the test problem (2D or 3D) via

.. code-block:: bash

      ./bin/flame tests/Flame/input

This will generate an output directory :code:`tests/Flame/output/`, and will re-name the old directory if it exists already.
Use VisIt to open :code:`tests/Flame/output/output.visit`.
:code:`tests/Flame/output/metadata` contains all of the input parameters.
The input file has the following form:

.. code-block:: bash

   # * = you probably don't need to change this

   amr.plot_file        = tests/Flame/output # The name of the output file 
   amr.plot_int		= 100                # How frequently to dump output
   amr.max_level        = 3                  # How many levels of grid refinement are allowed
   amr.n_cell		= 30 10 2            # Number of cells in the x,y,z directions
   #amr.max_grid_size	= 50                 # * Maximum allowable grid size
   amr.blocking_factor	= 1                  # * Minimum size of each block
   amr.regrid_int	= 10                 # * How frequently (in number of timesteps) to regrid
   amr.grid_eff		= 1.0                # * Closeness of regrid patch to tagged region (default = 0.7) 
   geometry.prob_lo	= -1.5 -0.5 -0.1     # Problem domain: furthest point to the lower left
   geometry.prob_hi	= 1.5 0.5 0.1        # Problem domain: furthest point to the upper right
   geometry.is_periodic	= 0 1 1              # * Boundary conditions - currently under construction 

   timestep		= 0.001       # Timestep
   stop_time            = 50          # Simulation stop time
   
   physics.M		= 1.0         # Mobility - controls flame speed
   physics.kappa	= 0.01	      # Numerical width of burn region
   physics.w1		= 1.0	      # Chemical energy before burn
   physics.w12		= 2.0	      # Chemical energy during burn
   physics.w0		= 0.0	      # Chemical energy after burn
   physics.rho1		= 1.0	      # Density before burn
   physics.rho0		= 0.0	      # Density after burn
   physics.k1		= 0.1         # Conductivity before burn
   physics.k0		= 0.0	      # Conductivity after burn
   physics.cp1		= 1.0	      # Specific heat before burn
   physics.cp0		= 1.0	      # Specific heat after burn
   physics.qdotburn	= 0.0         # Rate of heat loss due to flame
   physics.fs_number	= 100         # Number of flame speed regions (1 for homogeneous)
   physics.fs_min	= -0.5        # Minimum flame speed 
   physics.fs_max	= 0.5         # Maximum flame speed
   
   # Boundary conditions for temperature
   TempBC.hi	        = EXT_DIR INT_DIR INT_DIR 
   TempBC.lo	        = EXT_DIR INT_DIR INT_DIR 
   TempBC.lo_1	        = 0.0 
   TempBC.hi_1	        = 0.0 
   
   # Boundary conditions for order parameter
   EtaBC.hi	        = EXT_DIR INT_DIR INT_DIR
   EtaBC.lo	        = EXT_DIR INT_DIR INT_DIR
   EtaBC.lo_1	        = 0.0 
   EtaBC.hi_1	        = 1.0 


See also:
---------

- :bash:`./src/flame.cc`: Entry point for the Flame solver
- :ref:`API-Integrator-Flame`: Integrator that explicitly evolves the order parameter
- :ref:`API-Operator-Elastic`: Operator to do elastic solves

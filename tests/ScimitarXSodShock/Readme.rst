Sod Shock Tube Simulation
=========================

This demonstrates the evolution of a 1D Sod shock tube problem in a 2D domain, solving for pressure, velocity, internal energy, and density. The numerical results can be extracted and compared with the exact solution.

See :ref:`Integrator::ScimitarX` for additional information.

Extracting Data Using VisIt
--------------------------
The simulation data (pressure, velocity, internal energy, density) can be extracted using the script `pyscript_visit_extract.py`.

Run the following command from the terminal inside the output folder:

visit -cli -s visit_extract.py ./

*********************************************************************************************

Plot the results and Compare it with Exact Solution using the script `plot_sodshock.py`.

Run the following command to execute the script.

python3 pyscript_sodshock_plots.py


*********************************************************************************************



 

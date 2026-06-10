
This test simulates the deflagration of a solid composite propellant (SCP) consisting of packed AP in an HTPB matrix with voids.
It combines the phase field method for interface tracking with thermal transport.

The following *full-length* simulation results from running the input file directly with no modification.

.. figure:: ../../../tests/SCPThermalVoid/reference/initial_eta.png
   :align: center
	   The inital void is placed to take up most of the domain. Where eta=0 is void and eta=1 is solid
.. figure:: ../../../tests/SCPThermalVoid/reference/eta_mesh.gif
   :align: center
	   This is how the eta field evolves over time. AMR is used to only refine in areas of interest
.. figure:: ../../../tests/SCPThermalVoid/reference/temp.gif
   :align: center
	   This is the temperature field over time. The temperature is hottest in the HTPB binder (not shown)
**Notes**

-  The methods for this model and updates to the full-feedback model can be found here: https://doi.org/10.31224/7215
   



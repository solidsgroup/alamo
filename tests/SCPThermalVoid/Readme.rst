This test simulates the regression of a solid composite propellant (SCP) with a void in the domain.
An AP sphere is placed in each corner of the domain with a void in the middle,

.. figure:: ../../../tests/SCPThermalVoid/reference/initial_phi.png
   :align: center
.. figure:: ../../../tests/SCPThermalVoid/reference/initial_eta.png
   :align: center

as the domain burns regression rate increases around the void due to an increase in surface area,
but the maximum temperature is in the center of the AP particles becuase AP burns hotter than HTPB.

.. figure:: ../../../tests/SCPThermalVoid/reference/temp.gif
   :align: center

This test also uses the new remeshing to increase preformance, where when regridding if T<Tcutoff,
the inital eta field is used to prevent squares of eta appearing instead of circles

.. figure:: ../../../tests/SCPThermalVoid/reference/eta_mesh.gif
   :align: center

**Notes**

-  The entire simulation runs in ~20 minutes locally on 8 cores, the test case doesn't run to completation.
-  Discussion of the model and the methods can be found in the paper below

**References**

.. bibliography::
   :list: bullet
   :filter: False

   meier2024finite
   meier2024diffuse



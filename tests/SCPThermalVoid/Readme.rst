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

-  The domain is :code:`0.0005 x 0.0005`; however, the sphere packing file is
   defined over a larger domain (:code:`0.001 x 0.001`).
-  The parameters are set to be as course as possible in both space and time without causing instability or solver divergence.
   This means that the results may or may not be completely accurate.
   It also may mean that you see some jitter or spurious fluctuation in the elastic solution (especially in the stress and strain fields.)
   This should go away with increased resolution, so it is recommended that you increase the AMR level and decrease the timestep in order to get properly converged results.
-  The parameter :code:`elastic.solver.fixed_iter` replaces :code:`elastic.solver.max_iter`.
   This allows the solution to continue even if full convergence is not attained.
   However, as long as it is inside the Newton solver, the NR convergence check will prevent the solver from continuing unless it is within tolerance.

**References**

.. bibliography::
   :list: bullet
   :filter: False

   meier2024finite
   meier2024diffuse



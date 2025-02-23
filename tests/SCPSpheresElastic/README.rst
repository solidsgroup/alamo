
This test simulates the deflagration of a solid composite propellant (SCP) consisting of packed AP in an HTPB matrix.
It combines the phase field method for interface tracking with thermal transport and finite-deformation elasticity.

The following *full-length* simulation results from running the input file directly with no modification.
On eight processors, it should take about one hour to run.

.. figure:: ../../../tests/SCPSpheresElastic/reference/combine.gif
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



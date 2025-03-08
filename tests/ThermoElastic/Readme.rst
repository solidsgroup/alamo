
This test demonstrates multiphysics coupling using the :ref:`Integrator::ThermoElastic` 
integrator that inherits both from :ref:`Integrator::HeatConduction` and :ref:`Integrator::Mechanics`.
In the test, the domain is initially stress-free and at T=0.
Thermal evolution is driven by the lower and left boundaries that are at T=1.
The circular inclusion has a coefficient of thermal expansion, while the matrix 
does not.
So, as the inclusion heats up, it expands, inducing a stress distribution.

.. figure:: ../../../tests/ThermoElastic/reference/movie.gif
   :scale: 100%
   :align: center

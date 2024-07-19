This test problem illustrates the use of the phase field method when initialized with a Voronoi tesselation.
This example tests the following capabilities

- :ref:`Integrator::PhaseFieldMicrostructure`
- :ref:`IC::Voronoi`

The following is the result in 2D; the problem can be run in 3D as well.

.. image:: ../../../tests/Voronoi/reference/voronoi.gif
   :scale: 75%
   :align: center

Notice that there are 100 initial grains in the tesselation, but only 10 actual grains in the simulation.
Therefore, some grains will eventually merge together when the boundary between them disappears.

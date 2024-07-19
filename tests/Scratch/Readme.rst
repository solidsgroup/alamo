
This input file simulates the effect of a moving point load applied over the surface of the domain. 
The J2 plastic material accumulates residual stress left by the application of the point load.
This test exercises the following capabilities:

* :ref:`Integrator::Mechanics`
* :ref:`BC::Operator::Elastic::Expression`
* :ref:`Model::Solid::Affine::J2`

.. figure:: ../../../tests/Scratch/reference/movie.gif
   :scale: 50%
   :align: center

To generate the movie here, run in 3D.
A 3D test is not included in the regression test suite due to the extensive time required.
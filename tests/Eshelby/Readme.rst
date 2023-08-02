This example demonstrates the Mechanics solver using the `Eshelby inclusion problem <https://en.wikipedia.org/wiki/Eshelby%27s_inclusion>`_.
The example primarily tests the following Alamo capabilities

* :ref:`Integrator::Mechanics`
* :ref:`Model::Solid::Affine::Isotropic`

In 2D, the :code:`sigma_xx` solution is:

.. figure:: ../../../tests/Eshelby/reference/stress_xx_2d.png
   :scale: 50%
   :align: center

In 3D, the solution is tested against the exact solution derived using Eshelby inclusion theory.
For additional background see Runnels *et al*: https://doi.org/10.1016/j.jcp.2020.110065

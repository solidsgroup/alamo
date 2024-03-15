

This example simluates the evolution of a perturbed interface with strongly anisotropic grain boundary energy.
For background see Ribot *et al*: https://doi.org/10.1088/1361-651X/ab47a0

This problem tests the following

* :ref:`Integrator::PhaseFieldMicrostructure`
* :ref:`Model::Interface::GB`

Example output:

.. figure:: ../../../tests/PerturbedInterface/reference/perturbedinterface.png
   :scale: 50%
   :align: center

Notice that the sharp corner in the equilibrium shape is the result of strong anisotropy. :cite:t:`ribot2019new`

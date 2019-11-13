.. role:: cpp(code)
   :language: c++
.. role:: bash(code)
   :language: bash

Perturbed Interface
===================

For background see Ribot *et al*: https://doi.org/10.1088/1361-651X/ab47a0

This example simluates the evolution of a perturbed interface with strongly anisotropic grain boundary energy.
For best results (currently) use 2D.
The following image is the output in 2D.

.. image:: perturbedinterface.png
   :scale: 50%
   :align: center

Notice that the equilibrium shape has a sharp corner - this is the result of the anisotropy.

.. literalinclude:: ../../../tests/PerturbedInterface/input
   :caption: PerturbedInterface input file (tests/PerturbedInterface/input)
   :language: makefile


Test problem
------------

Run the test problem (2D or 3D) via

.. code-block:: bash

      ./bin/alamo... tests/PerturbedInterface/input

This will generate an output directory :code:`tests/PerturbedInterface/output/`, and will re-name the old directory if it exists already.
Use VisIt to open :code:`tests/Flame/output/output.visit`.
:code:`tests/Flame/output/metadata` contains all of the input parameters.
The input file has the following form:


See also:
---------

- :bash:`./src/alamo.cc`: Entry point for the solver
- :ref:`API-Integrator-PhaseFieldMicrostructure`: Integrator that explicitly evolves the order parameter
- :ref:`API-Operator-Elastic`: Operator to do elastic solves
- :ref:`API-IC-PerturbedInterface`: Parameterization of initial perturbation

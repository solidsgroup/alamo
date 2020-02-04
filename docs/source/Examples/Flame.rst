.. role:: cpp(code)
   :language: c++
.. role:: bash(code)
   :language: bash

Flame
=====

The flame simulation uses a phase field type model to simulate the burn front of a solid material undergoing burn.

.. literalinclude:: ../../../tests/Flame/input
   :caption: Flame input file (tests/Flame/input)
   :language: makefile


Run the test problem (2D or 3D) via

.. code-block:: bash

      ./bin/flame tests/Flame/input

This will generate an output directory :code:`tests/Flame/output/`, and will re-name the old directory if it exists already.
Use VisIt to open :code:`tests/Flame/output/output.visit`.

See also:

- :bash:`./src/flame.cc`: Entry point for the Flame solver
- :ref:`API-Integrator-Flame`: Integrator that explicitly evolves the order parameter
- :ref:`API-Operator-Elastic`: Operator to do elastic solves

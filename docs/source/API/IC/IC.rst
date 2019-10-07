IC
---

This namespace contains a set of Initial Condition (`IC`) objects.
These specify the inital geometry that the area of interest for each problem assumes.
It should be noted that this namespace includes a purely abstract IC object that the other IC objects inherit from.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Affine
   Constant
   Cylinder
   PerturbedInterface
   Random
   Sphere
   TabulatedInterface
   Trig
   Trig2
   Voronoi
   Wedge

.. doxygenclass:: IC::IC
   :project: alamo
   :members:
   :protected-members:
   :private-members:
   
.. todolist::


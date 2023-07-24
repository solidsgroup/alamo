This example demonstrates the Finite Kinematics newton solver for a composite consisting of a compliant rubber-like material and a (near) rigid metal or ceramic-like material.
The difference in elastic moduli is 10x.
The low faces (xlo, ylo, zlo) are homogeneous dirichlet in the normal direction and naumann in the lateral directions.
The three-dimensional test yields the following:

.. figure:: ../../../tests/RubberWithInclusion/reference/movie.gif
   :scale: 50%
   :align: center

For regression testing purposes, only the two-dimensional case is considered for efficiency.
In two dimensions, the distortion of the mesh is more severe, which usually means that additional Newton iterations are required.
           
.. figure:: ../../../tests/RubberWithInclusion/reference/movie-2d.gif
   :scale: 50%
   :align: center




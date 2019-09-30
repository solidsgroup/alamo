Sphere
------
Define a area of interest in the shape of either a sphere or cylinder

There are four methods of applying this IC:

   * XYZ: This defines the area as a 3D sphere, with the center and radius provided.

   * XY, YZ, XZ: These three options define the area as a cylinder,
     with the circular faces parallel with the plane chosen.
     In these cases, only the radius input is actually used

Both methods assign the field values assigned to the relevant nodes.


.. doxygenclass:: IC::Sphere
   :project: alamo
   :members:
   :protected-members:
   :private-members:
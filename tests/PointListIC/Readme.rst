
This test is for the point list initial condition. It tests the IC by placing a rotated NACA 9612 airfoil on in the domain.
The point list IC works with one or more closed polygons that can be rotated counter-clockwise about a center of rotation. See below for example figure

.. figure:: ../../../tests/PointListIC/alpha_ic_img.png
   :align: center

Each line of the input file gives a point as ``x y z ObjNum``. The ``ObjNum`` column is optional; if it is
omitted, every point in the file is treated as a single polygon. When it is present, a new polygon starts
wherever ``ObjNum`` changes by more than ``obj_num_threshold`` (default ``0.5``) from the previous point.

By default every polygon is unioned together (via a pointwise max) to form the solid field. Individual
polygons can instead be marked as void/negative regions using the ``invert`` parameter, which takes one
entry per polygon (nonzero = void). Void polygons are unioned separately and then multiplied into the
solid field, cutting holes wherever they overlap it. ``invert_all`` is a separate, single flag that flips
the whole finished field (solid field with voids cut out) at the end, e.g. to make the IC 1 outside the
solid region instead of inside it.


.. _API-Numeric-Stencil:

.. role:: cpp(code)
   :language: c++

Stencil
=======

Very high order derivatives (greater than second) are needed for some operations
in Alamo. These functions compute a variety of derivatives using different kinds
of stencils.

All functions are compiled with the :cpp:`AMREX_FORCE_INLINE` directive. It is
*extremely important* to ensure that there are no inlining errors when compiling,
as any failure to inline these functions will result in a *huge* performance hit!

To use: if you have an :cpp:`amrex::Array4` called :cpp:`phi`, and you need to compute

.. math::
   :nowrap:

      \begin{gather*}
      \frac{\partial^4}{\partial x_1^3 \partial x_2}\text{(phi)}
      \end{gather*}

You can do it using the following code:

.. code-block:: c++
   :linenos:

   #include "Numeric/Stencil.H"

   ...

   Set::Scalar derivative = Numeric::Stencil::D<3,1,0>(phi,i,j,k,DX[0],DX[1],DX[2])

where :cpp:`DX` is an array of grid dimensions, and :cpp:`i,j,k` are the lattice coordinates.
You can optionally specify stencil types using 
.. code-block:: c++

   Numeric::Stencil::Type::Lo;
   Numeric::Stencil::Type::Hi;
   Numeric::Stencil::Type::Central;

where :cpp:`Central` is always the default.



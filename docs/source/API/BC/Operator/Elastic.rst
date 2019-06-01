.. _API-BC-Operator-Elastic:

.. role:: cpp(code)
   :language: c++

Elastic
=======
There are four boundary conditions that can be used with :cpp:`Operator::Elastic`:

Dirchlet/Displacement
---------------------
For a Dirichlet point the Dirichlet boundary operator is used, :math:`\f$D_{\Omega_1}\f$`.
This is nothing other than the identity operator, i.e.

:math:`D_{\Omega_1}(\mathbf{u}) = \mathbf{u},`

that is, it is the identity.
The value of the displacement is thus determined by the value of the corresponding point
in the right hand side during the solve.

Natural/Traction
----------------
This operator is similar to (but not the same as) a typical Neumann boundary condition.
The operator, denoted :math:`D_{\Omega_2}` is defined as

:math:`D_{\Omega_2}(\mathbf{u}) = \mathbb{C}\nabla\mathbf{u} \cdot \mathbf{n} = \sigma\cdot\mathbf{n}`

The output of this operator is a surface traction. The value of the surface traction is
specified in the right hand side during the solve.

Periodic
--------

Neumann
-------


.. doxygenclass:: BC::Operator::Elastic
   :project: alamo
   :members:
   :protected-members:
   :private-members:

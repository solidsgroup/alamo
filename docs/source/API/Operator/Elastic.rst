.. _API-Operator-Elastic:

.. role:: cpp(code)
   :language: c++

Elastic
=======

Stress divergence operator definition
-------------------------------------

The stress divergence equation is

:math:`f_i = D_{\Omega}(\mathbf{u})_i = \frac{\partial}{\partial x_j}\mathbb{C}_{ijkl}\frac{\partial u_k}{\partial x_l} = \mathbb{C}_{ijkl}u_{k,lj} + \mathbb{C}_{ijkl,j}u_{k,l}`

Where :math:`\mathbf{f}` is the vector of nodal forces, :math:`\mathbf{u}` the vector of nodal displacements, and :math:`\mathbb{C}_{ijkl}` the fourth-order elasticity tensor.
 
Boundary operator
-----------------

If a point is located on the boundary then a boundary operator is used rather than the stress divergence operator.
This can be found on the :cpp:`::BC::Operator::Elastic` page of the documentation.

Solid Model Template Class
--------------------------

This operator is templated with class T.
T should be of type :cpp:`Model::Solid::Solid` or derived.
These models are data structures containing all necessary elastic constants to compute the
stress tensor :math:`\sigma` given :math:`\nabla\mathbf{u}`.
The models are stored in an :cpp:`amrex::FabArray` of :cpp:`amrex::BaseFab<T>`.
Models have basic arithmetic operators defined, and use the parentheses operator to do the calculation.
 




.. doxygenclass:: Operator::Elastic
   :project: alamo
   :members:
   :protected-members:
   :private-members:

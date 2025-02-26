
This is a 1D test to determine Riemann solver performance by comparing against the known
solution for the shock tube problem.
This implements Example 5.2 from *Computational Gasdynamics* by Culbert B. Laney (page 91):

.. admonition:: Example 5.2 
   :class: reference

   Find the solution to Roe’s approximate Riemann problem at :math:`t=0.01s`
   if :math:`p_L=100,000 N/m^2`, :math:`\rho_L=1kg/m^3`, :math:`u_L=100m/s` and
   :math:`p_R=10,000N/m^2`, :math:`\rho_R=0.125kg/m^3`, :math:`u_R=-50m/s`



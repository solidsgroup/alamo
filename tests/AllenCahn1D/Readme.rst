**Derivation of the exact solution**:
The Allen Cahn equation admits a 1D solution that can be used to test the implementation.

* The energy functional in 1D is

  .. math ::
  
     F = \int_\Omega \Big(\lambda W(\alpha) + \frac{1}{2}\kappa \Big(\frac{d\alpha}{dx}\Big)^2 \Big)d\,\mathbf{x}

   
* At steady-state, the Euler-Lagrange equation is:

  .. math ::
  
     -\frac{1}{L}\frac{\partial\alpha}{\partial t} = \lambda \frac{dW}{\partial \alpha} - \kappa \frac{d^2\alpha}{dx^2} = 0

  Or, simply

  .. math ::
    
     \lambda \frac{dW}{\partial \alpha} = \kappa \frac{d^2\alpha}{dx^2} 

* This is a nonlinar ODE, however, it can be reduced by multiplying both sides by :math:`\frac{d\alpha}{dx}`  (prime indicates derivative with respect to argument)

  .. math :: 
    
     \lambda\frac{dW}{d\alpha}\frac{d\alpha}{dx} &= \lambda\frac{dW}{dx}
     
     \kappa\frac{d^2\alpha}{dx^2}\frac{d\alpha}{dx^2} &= \frac{1}{2}\kappa\frac{d}{dx}\Big(\frac{d\alpha}{dx}\Big)^2 
     
  Substituting into the above expression we obtain:

  .. math :: 
    
     \lambda\frac{dW}{dx} = \frac{1}{2}\kappa\frac{d}{dx}\Big(\frac{d\alpha}{dx}\Big)^2 
  
  which, when integrated, becomes

  .. math :: 
    
     \frac{d\alpha}{dx}  = \sqrt{\frac{2\lambda}{\kappa} W(\alpha) + C_0} 

  We will let :math:`C_0` be zero for the moment. 
  
* Next, we substitute the expression for :math:`W`

  .. math :: 
    
     W(\alpha) = \alpha^2(1-\alpha)^2
  
  which has a very convenient square root.
  
* Now we have a nonlinar, but almost solvable differential equation:

  .. math :: 

     \frac{d\alpha}{dx} = \gamma \alpha(1-\alpha)

  where we define :math:`\gamma=\pm\sqrt{\frac{2\lambda}{\kappa}}`.
  This differential equation can be solved by the substitution

  .. math :: 

     u = \alpha^{-1} \implies \alpha=u^{-1}, \frac{d\alpha}{dx} = -u^{-2}\frac{du}{dx}
     
  which, when substituting, gives

  .. math :: 

     \frac{du}{dx} + \gamma u = \gamma
     
  This has the straightforward solution

  .. math :: 

     u(x) = C_1e^{-\gamma x} + 1

  resulting in the following solution for :math:`\alpha`

  .. math :: 

     \alpha(x) = \frac{1}{C_1e^{-\gamma x} + 1}
  
  
  We desire a solution centered at :math:`x=0`, i.e. such that :math:`y(-x)=1-y(x)`.
  One can see that this is satisfied only when :math:`C_1=1`.
  So the interfacial solution is:

  .. math :: 

     \boxed{\alpha(x) = \frac{1}{e^{-\gamma x} + 1}}
  
**Determination of interface width** :math:`\varepsilon`: we estimate the interface width based on the slope of :math:`\alpha` at :math:`x=0`.
The slope at the center point is

.. math :: 

   \frac{d+\alpha}{dx}\Bigg|_{x=0}
   = \frac{\gamma e^{-\gamma x}}{(e^{-\gamma x} + 1)^2}\Bigg|_{x=0}
   = \frac{\gamma}{4} \approx \frac{1}{\varepsilon} \implies \varepsilon = \boxed{\sqrt{\frac{8\kappa}{\lambda}}}
     
**Determination of interface energy** :math:`\sigma`: this is given by substituting the solution into the free energy functional:

.. math ::

   \sigma = F\Big[\frac{1}{e^{-\lambda x} + 1}\Big]
  

The gradient portion is:

.. math ::

   \int_{-\infty}^\infty \frac{1}{2}\kappa \Big(\frac{d\alpha}{dx}\Big)^2\,dx = \frac{1}{12}\kappa\gamma
   = \frac{1}{6}\sqrt{\frac{\kappa\lambda}{2}}
 
And the chemical potential portion is

.. math ::

   \int_{-\infty}^\infty \lambda W[\alpha] dx = \lambda \frac{5}{6\gamma}
   = \frac{5}{6}\sqrt{\frac{\kappa\lambda}{2}}

Giving the final relationship

.. math ::

   \boxed{\sigma = \sqrt{\frac{\kappa\lambda}{2}}}


**Inverting the relationship**: It may be more desirable to define behavior in terms of interface energy and width.
One can then show that

.. math ::

   \boxed{\lambda = \frac{4 \,\sigma}{\varepsilon}}, \, \, \,
   \boxed{\kappa = \frac{\sigma\varepsilon}{2}}

Writing the free energy in terms of these quantities gives

  .. math ::
  
     F = \int_\Omega \Big(\frac{4\sigma}{\varepsilon} W(\alpha) + \frac{\sigma\varepsilon}{4} \Big(\frac{d\alpha}{dx}\Big)^2 \Big)d\,\mathbf{x}


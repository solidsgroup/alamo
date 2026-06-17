**Derivation of the exact solution**

The Allen Cahn equation admits a 1D solution that can be used to test the implementation. The energy functional in 1D is

.. math ::
  
   F[\alpha]=\int_\Omega\left(\lambda W(\alpha)+\frac{1}{2}\kappa\left(\frac{\partial\alpha}{\partial x}\right)^2\right)dx

The Allen-Cahn equation is

.. math ::

   \frac{\partial\alpha}{\partial t}=-L\frac{\delta F}{\delta\alpha}

where the variational derivative can be derived from the Euler-Lagrange equation resulting in

.. math ::

   \frac{\partial\alpha}{\partial t}=-L\left(\lambda\frac{dW}{d\alpha}-\kappa\frac{\partial^2\alpha}{\partial x^2}\right).
   
At steady-state, (:math:`\frac{\partial\alpha}{\partial t}=0`) the Allen-Cahn equation simplifies to

.. math ::

   \lambda\frac{dW}{d\alpha}=\kappa\frac{d^2\alpha}{dx^2}.

This is a nonlinear ODE, but it can be reduced by multiplying both sides by :math:`\frac{d\alpha}{dx}`:

.. math ::

   \lambda\frac{dW}{d\alpha}\frac{d\alpha}{dx}=\kappa\frac{d^2\alpha}{dx^2}\frac{d\alpha}{dx}

Using the chain rule, we can further simplify:

.. math ::
  
   \lambda\frac{dW}{dx}=\frac{\kappa}{2}\frac{d}{dx}\left(\frac{d\alpha}{dx}\right)^2

When integrated, we obtain

.. math ::

   \int\lambda\frac{dW}{dx}dx&=\int\frac{\kappa}{2}\frac{d}{dx}\left(\frac{d\alpha}{dx}\right)^2dx

   \lambda W(\alpha)+C_1&=\frac{\kappa}{2}\left(\frac{d\alpha}{dx}\right)^2+C_2

Combining the constants of integration into a single constant :math:`C_0=C_2-C_1`:

.. math ::

   \lambda W(\alpha)=\frac{\kappa}{2}\left(\frac{d\alpha}{dx}\right)^2+C_0

Solving for :math:`\frac{d\alpha}{dx}`:

.. math ::

   \frac{d\alpha}{dx}=\pm\sqrt{\frac{2}{\kappa}(\lambda W(\alpha)-C_0)}

For a monotonic interface profile connecting two phases, we require :math:`C_0=0` to ensure that :math:`\frac{d\alpha}{dx}=0` at the minima of :math:`W(\alpha)`:

.. math ::

   \frac{d\alpha}{dx}=\pm\sqrt{\frac{2\lambda}{\kappa}W(\alpha)}

As is often the convention, we let :math:`W` represent a double-well chemical potential

.. math :: 
    
   W(\alpha)=\alpha^2(1-\alpha)^2
  
which has a convenient square root. After substituting, we now have a nonlinear, but trivially solved differential equation:

.. math :: 

   \frac{d\alpha}{dx}=\gamma\alpha(1-\alpha)

where :math:`\gamma=\pm\sqrt{\frac{2\lambda}{\kappa}}`. This differential equation can be solved by separation of variables or by substitution as shown here:

.. math :: 

   u=\alpha^{-1}\implies\alpha=u^{-1},\frac{d\alpha}{dx}=-u^{-2}\frac{du}{dx}

Substituting :math:`u` for :math:`\alpha` we find

.. math :: 

   \frac{du}{dx}+\gamma u=\gamma,
     
which has the solution:

.. math :: 

   u(x)=Ce^{-\gamma x}+1

where :math:`C` is a constant. Reversing the substitution, we get the following solution for :math:`\alpha`:

.. math :: 

   \alpha(x)=\frac{1}{Ce^{-\gamma x}+1}
  
We desire a solution centered at :math:`x=0` such that :math:`\alpha(-x)=1-\alpha(x)`. This condition is only satisfied when :math:`C=1`; so, the interfacial solution is:

.. math ::

   \boxed{\alpha(x)=\frac{1}{e^{-\gamma x}+1}}
  
**Determination of interface width** :math:`\varepsilon`

The slope at the center of the interface (:math:`x=0`) is

.. math :: 

   \left.\frac{d\alpha}{dx}\right|_{x=0}=\left.\frac{\gamma e^{-\gamma x}}{(e^{-\gamma x} + 1)^2}\right|_{x=0}=\frac{\gamma}{4}

The function changes from approximately 0 to 1 across the interface. Using a linear approximation slope—*i.e.*, slope :math:`\approx` (change in :math:`\alpha`)/(change in :math:`x`)—the interface width :math:`\varepsilon` is as follows:

.. math ::
  
   \frac{\gamma}{4}\approx\frac{1}{\varepsilon}\implies\varepsilon=\frac{4}{\gamma}=\boxed{\sqrt{\frac{8\kappa}{\lambda}}}
     
**Determination of interface energy density** :math:`\sigma`

The interface energy :math:`\sigma` is given by substituting the solution for :math:`\alpha` into the free energy functional:

.. math ::

   \sigma&=F\left[\frac{1}{e^{-\gamma x}+1}\right]
   
   &=\int_{-\infty}^\infty\left(\lambda W(\alpha)+\frac{1}{2}\kappa\left(\frac{d\alpha}{dx}\right)^2\right)dx

Using the chain rule, the derivative of :math:`\alpha` is

.. math ::

   \frac{d\alpha}{dx}=-\frac{1}{\left(e^{-\gamma x}+1\right)^2}\left(-\gamma e^{-\gamma x}\right)=\frac{\gamma e^{-\gamma x}}{\left(e^{-\gamma x}+1\right)^2}.

Conveniently, this can be written as

.. math ::

   \frac{d\alpha}{dx}=\gamma\alpha(1-\alpha).

Squaring this derivative, we can make another convenient substitution:

.. math ::

   \left(\frac{d\alpha}{dx}\right)^2=\gamma^2\alpha^2(1-\alpha)^2=\gamma^2W(\alpha)

Plugging this back into the integrand:

.. math ::

   \sigma=\int_{-\infty}^\infty\left(\lambda W(\alpha)+\frac{1}{2}\kappa\gamma^2W(\alpha)\right)dx

Recalling that :math:`\gamma^2=\frac{2\lambda}{\kappa}`:

.. math ::

   \sigma&=2\lambda\int_{-\infty}^{\infty}W(\alpha)dx

   &=2\lambda\int_{-\infty}^\infty\alpha^2(1-\alpha)^2dx

We can trivially solve this integral with a simple substitution:

.. math ::

   u=\alpha&\implies\frac{du}{dx}=\frac{d\alpha}{dx}=\gamma\alpha(1-\alpha)

   &\implies dx=\frac{du}{\gamma\alpha(1-\alpha)}=\frac{du}{\gamma u(1-u)}

To adjust the bounds of integration, we note:

.. math ::

   x&:-\infty\rightarrow\infty
   
   u&:0\rightarrow1

Substituting these variables into the integrand:

.. math ::

   \sigma&=2\lambda\int_0^1u^2(1-u)^2\frac{du}{\gamma u(1-u)}

   &=\frac{2\lambda}{\gamma}\int_0^1u(1-u)du

   &=\frac{\lambda}{3\gamma}

Applying one final substitution, we find the following expression for the interface energy :math:`\sigma`:

.. math ::

   \boxed{\sigma=\frac{1}{3}\sqrt{\frac{\lambda\kappa}{2}}}

**Inverting the relationship**

It may be more desirable to define behavior in terms of the interface energy :math:`\sigma` and interface width :math:`\varepsilon`. From the previously derived equation for :math:`\sigma`, we can show:

.. math ::

   \lambda&=\frac{18\sigma^2}{\kappa}

   \kappa&=\frac{18\sigma^2}{\lambda}

From the definition of :math:`\gamma`, we can use the following substitutions:

.. math ::

   \kappa&=\frac{2\lambda}{\gamma^2}

   \lambda&=\frac{\kappa\gamma^2}{2}

These relationships then give us the following expressions:

.. math ::

   \lambda&=\frac{9\sigma^2\gamma^2}{\lambda}=3\sigma\gamma

   \kappa&=\frac{36\sigma^2}{\kappa\gamma^2}=\frac{6\sigma}{\gamma}


Recalling that :math:`\gamma=\frac{4}{\varepsilon}`, we can finalize the inversion:

.. math ::

   \lambda&=\frac{12\sigma}{\varepsilon}

   \kappa&=\frac{3\sigma\varepsilon}{2}

Writing the free energy in terms of these quantities gives:

.. math ::
  
   F\left[\alpha\right]=\int_\Omega\left(\frac{12\sigma}{\varepsilon}W(\alpha)+\frac{3\sigma\varepsilon}{4}\left(\frac{d\alpha}{dx}\right)^2\right)dx

**Dimensional analysis**

Given that :math:`\sigma` has units of :math:`\left[\mathrm{E}\ \mathrm{L}^{-2}\right]` and :math:`\varepsilon` has units of length :math:`\left[\mathrm{L}\right]`, it becomes clear from the inverted relationships that :math:`\lambda` and :math:`\kappa` must have the following units:

.. math ::

   \left[\lambda\right]&=\left[\mathrm{E}\ \mathrm{L}^{-3}\right]

   \left[\kappa\right]&=\left[\mathrm{E}\ \mathrm{L}^{-1}\right]


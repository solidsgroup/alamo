SH (Spherical Harmonic) 
--------------------------


Form of the function in terms of the angles :math:`\theta,\phi`

.. math::

   W(\phi,\theta) &= \sigma_0 - \sigma_1 \cos^2(2 \phi) \sin^2(2\theta)

   \frac{\partial W}{\partial \phi} &= + 4\sigma_1 \sin(2\phi) \cos(2 \phi) \sin^2(2\theta)

   \frac{\partial W}{\partial \theta} &= - 4\sigma_1 \cos^2(2 \phi) \cos(2\theta)\sin(2\theta)
   
Alternatively we have :math:`W(\mathbf{n})` where :math:`\mathbf{n}` is a normal vector, such that

.. math::

   \begin{align}
   \mathbf{n} &= 
   \begin{bmatrix}
   \sin\theta\cos\phi \\ \sin\theta\sin\phi \\ \cos\theta
   \end{bmatrix} &
   \theta &= \cos^{-1}(n_2) &
   \phi &= \tan^{-1}\Big(\frac{n_1}{n_0}\Big)
   \end{align}

Energy as a function of normal vector is implemented as the simple pushforward :math:`W(\theta(\mathbf{n}),\phi(mathbf{n}))`.
However the derivative is more complex.
The derivative operation is :math:`DS_2 \to TS_2`, the tangent bundle on the unit 2-sphere.
The most useful form, however, is the embedding of the derivative within :math:`\mathbb{R}^3`.
The following form is readily implemented:

.. math::

   DW = \frac{\partial W}{\partial\mathbf{n}} - \Big(\frac{\partial W}{\partial \mathbf{n}}\cdot \mathbf{n}\Big)\cdot\mathbf{n}

ensuring that :math:`DW\in T S_2` (i.e. the derivative is tangent to the sphere at all points.)
Computing for all of the components:

.. math::

   \frac{\partial W}{\partial\mathbf{n}} =
   \begin{bmatrix}
   W_\phi \, \partial\phi / \partial n_0 \\ W_\phi \, \partial\phi / \partial n_1 \\ W_\theta \, \partial \theta / \partial n_2
   \end{bmatrix}
   =
   \begin{bmatrix}
   - W_\phi \, n_1 / \sqrt{n_0^2 + n_1^2} \\ W_\phi \, n_0 / \sqrt{n_0^2 + n_1^2} \\ - W_\theta / \sqrt{1-n_2^2}
   \end{bmatrix}
   =
   \begin{bmatrix}
     - W_\phi \, \sin\theta \\ W_\phi \, \cos\theta \\ - W_\theta / \sin\phi
   \end{bmatrix}

Then computing the normal projection

.. math::
   
   \frac{\partial W}{\partial\mathbf{n}} \cdot\mathbf{n} =  - W_\theta \cos\phi/ \sin\phi

   \Big(\frac{\partial W}{\partial\mathbf{n}} \cdot\mathbf{n}\Big)\mathbf{n}
   =
   \begin{bmatrix}
   - W_\theta \cos\phi\cos\theta\\
   - W_\theta \cos\phi\sin\theta \\
   - W_\theta \cos^2\phi / \sin\phi
   \end{bmatrix}

Therefore the final expression is

.. math::

   DW = 
   \begin{bmatrix}
     - W_\phi \, \sin\theta + W_\theta \cos\phi\cos\theta\\
       W_\phi \, \cos\theta + W_\theta \cos\phi\sin\theta\\
     - W_\theta \sin\phi
   \end{bmatrix}

.. danger::

   This does not actually work.
   Instead, a numerical approximation is implemented for :code:`DW` and :code:`DDW`


.. doxygenclass:: Model::Interface::GB::SH
   :project: alamo
   :members: 
   :protected-members:
   :private-members:



.. role:: cpp(code)
   :language: c++

CG
--

This documents the function :code:`Solver::Local::CG` which is a basic solver for 
linear systems.
The following example illustrates one possible usage case:

.. code-block:: cpp

    // Create a linear elastic cubic model
    Model::Solid::Linear::Cubic model(2.0,3.0,1.0);
    // Create a strain tensor
    Set::Matrix eps = Set::Matrix::Zero();
    // Create a stress tensor
    Set::Matrix sig = Set::Matrix::Zero();
    // Solve a CONSTRAINED problem: we want the solution 
    // to satisfy eps(0,0) = 0.1. So set eps equal to that value,
    // and set the value of mask to 1 to indicate to the solver 
    // not to change it.
    Set::iMatrix mask = Set::iMatrix::Zero();
    eps(0,0) = 0.1;    mask(0,0) = 1;
    // Run the CG solver and get the output
    eps = Solver::Local::CG(model.DDW(eps),sig,eps,mask,true);

.. NOTE::
    
    This is a rudimentary implementation with very few features.
    It will eventually be upgraded to a fully-featured C++ class, so 
    keep in mind that you may need to modify your code if you use it.

.. doxygenfunction:: Solver::Local::CG
   :project: alamo



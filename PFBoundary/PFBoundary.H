#ifndef PFBoundary_H
#define PFBoundary_H

#include "AMReX.H"

class PFBoundary
{
public:
  PFBoundary() {};
  virtual amrex::Real operator()(amrex::Real theta) = 0;
};


#endif

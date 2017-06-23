#ifndef PFBoundary_H
#define PFBoundary_H

#include "AMReX.H"

class PFBoundary
{
public:
  PFBoundary() {};
  virtual amrex::Real W(amrex::Real theta) = 0;
  virtual amrex::Real DW(amrex::Real theta) = 0;
};


#endif

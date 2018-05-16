#ifndef PFAmrBC_H_
#define PFAmrBC_H_

#include <AMReX_PhysBCFunct.H>

// In this tutorial, we assume the boundaries are all periodic.  
// Thus there are no physical boundaries. And FillBoundary does nothing.

class PFAmrPhysBC
  : public amrex::PhysBCFunctBase
{
public:
  PFAmrPhysBC () {}
  virtual ~PFAmrPhysBC () {}
  virtual void FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real time) override { }
  using amrex::PhysBCFunctBase::FillBoundary;
};

#endif

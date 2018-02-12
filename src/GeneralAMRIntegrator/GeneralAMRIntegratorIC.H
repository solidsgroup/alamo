#ifndef GeneralAMRIntegratorIC_H_
#define GeneralAMRIntegratorIC_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>

class GeneralAMRIntegratorIC
{
public:
  GeneralAMRIntegratorIC (amrex::Array<amrex::Geometry> &_geom)
    : geom(_geom) {} ;
  
  virtual void Initialize(const int lev, amrex::Array<std::unique_ptr<amrex::MultiFab> > &fab) = 0;
  
protected:
  amrex::Array<amrex::Geometry> &geom;
};

#endif

#ifndef IC_IC_H_
#define IC_IC_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>

namespace IC
{
class IC
{
public:
  IC (amrex::Array<amrex::Geometry> &_geom)
    : geom(_geom) {} ;
  
  virtual void Initialize(const int lev, amrex::Array<std::unique_ptr<amrex::MultiFab> > &fab) = 0;
  
protected:
  amrex::Array<amrex::Geometry> &geom;
};
}
#endif

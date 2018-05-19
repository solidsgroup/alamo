#ifndef BC_BC_H_
#define BC_BC_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>

#define NEUMANN 5

namespace BC
{
class BC
  : public amrex::PhysBCFunctBase
{

public:
  //BC (amrex::Vector<amrex::Geometry> &_geom) : geom(_geom) {};
  BC (amrex::Vector<amrex::Geometry> &_geom,
      amrex::Vector<std::string> bc_hi_str,
      amrex::Vector<std::string> bc_lo_str,
      amrex::Vector<amrex::Real> _bc_lo_1,
      amrex::Vector<amrex::Real> _bc_hi_1,
      amrex::Vector<amrex::Real> _bc_lo_2,
      amrex::Vector<amrex::Real> _bc_hi_2,
      amrex::Vector<amrex::Real> _bc_lo_3,
      amrex::Vector<amrex::Real> _bc_hi_3);
  virtual void FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real /*time*/) override;
  void SetLevel(int _lev);
  amrex::BCRec GetBCRec() ;

protected:
  std::vector<amrex::Geometry> &geom;
  int bc_lo[BL_SPACEDIM];
  int bc_hi[BL_SPACEDIM];
  int bc_hi_flux[BL_SPACEDIM];
  int bc_lo_flux[BL_SPACEDIM];
  amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
  amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
  amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;

  int lev=0;

};
}
#endif

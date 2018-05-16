#ifndef IC_VORONOI_H_
#define IC_VORONOI_H_

#include "IC/IC.H"

namespace IC
{
class Voronoi : public IC
{
public:
  Voronoi (amrex::Vector<amrex::Geometry> &_geom, int _number_of_grains);
  
  void Initialize(const int lev,
		  amrex::Vector<std::unique_ptr<amrex::MultiFab> > &field);
  
private:
  int number_of_grains;
  amrex::Vector<amrex::Real> voronoi_x;
  amrex::Vector<amrex::Real> voronoi_y;
#if BL_SPACEDIM==3
  amrex::Vector<amrex::Real> voronoi_z;
#endif

};
}
#endif

#ifndef IC_VORONOI_H_
#define IC_VORONOI_H_

#include "Set/Set.H"
#include "IC/IC.H"

namespace IC
{
class Voronoi : public IC
{
public:
  enum Type {Partition, Values};

  Voronoi (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {}
  Voronoi (amrex::Vector<amrex::Geometry> &_geom, int _number_of_grains, Set::Scalar a_alpha = 1.0) : IC(_geom) 
  {
    Define(_number_of_grains,a_alpha);
  }

  void Define (int _number_of_grains,
               Set::Scalar a_alpha = 1.0)
  {
    Define(_number_of_grains, std::vector<Set::Scalar>(_number_of_grains,a_alpha),Type::Partition);
  }

  void Define (int _number_of_grains,
               std::vector<Set::Scalar> a_alpha,
               Type t = Type::Values);
  
  void Add(const int lev,
		  amrex::Vector<amrex::MultiFab * > &field);
  
private:
  int number_of_grains;
  std::vector<Set::Scalar> alpha;
  std::vector<Set::Vector> voronoi;
  Type type;
  amrex::Vector<amrex::Real> voronoi_x;
  amrex::Vector<amrex::Real> voronoi_y;
#if BL_SPACEDIM==3
  amrex::Vector<amrex::Real> voronoi_z;
#endif

};
}
#endif

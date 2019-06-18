#ifndef IC_FLAME_WEDGE_H_
#define IC_FLAME_WEDGE_H_

#include "IC/IC.H"

namespace IC
{
class Wedge : public IC
{
public:
	Wedge (amrex::Vector<amrex::Geometry> &_geom);
  
	void Add(const int lev,
			amrex::Vector<amrex::MultiFab * > &field);
  
private:
};
}
#endif

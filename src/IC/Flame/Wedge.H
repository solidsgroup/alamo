#ifndef IC_FLAME_WEDGE_H_
#define IC_FLAME_WEDGE_H_

#include "IC/IC.H"

namespace IC
{
namespace Flame
{
class Wedge : public IC
{
public:
	Wedge (amrex::Vector<amrex::Geometry> &_geom);
  
	void Initialize(const int lev,
			amrex::Vector<std::unique_ptr<amrex::MultiFab> > &field);
  
private:
};
}
}
#endif

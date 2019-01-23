#ifndef IC_CONSTANT_H_
#define IC_CONSTANT_H_

#include "IC/IC.H"

/// \class Constant
/// \brief 
///
namespace IC
{
class Constant : public IC
{
public:
	Constant (amrex::Vector<amrex::Geometry> &_geom, amrex::Vector<amrex::Real> _value);  

	void Add(const int lev,amrex::Vector<amrex::MultiFab * > &field);

private:
	amrex::Vector<amrex::Real> value;
};
}
#endif

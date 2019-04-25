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
	Constant (amrex::Vector<amrex::Geometry> &_geom, std::vector<amrex::Real> _value);  
	~Constant() { };
	void Add(const int lev,amrex::Vector<amrex::MultiFab * > &field);

private:
	std::vector<amrex::Real> value;
};
}
#endif

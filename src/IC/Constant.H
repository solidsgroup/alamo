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
  	Constant (amrex::Vector<amrex::Geometry> &_geom) :
    IC::IC(_geom)
	{}
	Constant(amrex::Vector<amrex::Geometry> &_geom, std::vector<amrex::Real> _value)
	: IC(_geom), value(_value)
	{}
	~Constant() { };
	void Add(const int lev, amrex::Vector<amrex::MultiFab * > &field)
	{
		for (amrex::MFIter mfi(*field[lev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.tilebox();
	
			amrex::FArrayBox &field_box = (*field[lev])[mfi];
	
			for (int i = box.loVect()[0]-field[lev]->nGrow(); i<=box.hiVect()[0]+field[lev]->nGrow(); i++) 
				for (int j = box.loVect()[1]-field[lev]->nGrow(); j<=box.hiVect()[1]+field[lev]->nGrow(); j++)
	#if BL_SPACEDIM==3
					for (int k = box.loVect()[2]-field[lev]->nGrow(); k<=box.hiVect()[2]+field[lev]->nGrow(); k++)
	#endif
				{
					for (unsigned int m = 0; m < value.size(); m++)
						field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),m) = value[m];     
				}
		}
	
	}
private:
	std::vector<amrex::Real> value;
};
}
#endif

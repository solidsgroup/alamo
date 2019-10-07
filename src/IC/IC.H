#ifndef IC_IC_H_
#define IC_IC_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Util/Util.H"

namespace IC
{
/// \brief Pure abstract IC object from which all other IC objects inherit.
class IC
{
public:
	IC (amrex::Vector<amrex::Geometry> &_geom)
		: geom(_geom) {} ;
	virtual ~IC() {}

	virtual void Add(const int lev, amrex::Vector<amrex::MultiFab *> &field) = 0;
/// Function that IC's inherit to initialize the geometry's mesh
	void Initialize(const int lev,
			amrex::Vector<amrex::MultiFab *> &field)
	{
		field[lev]->setVal(0.0);
		Add(lev,field);
	};
	void Initialize(const int lev,
			amrex::Vector<std::unique_ptr<amrex::MultiFab> > &field)
	{
		amrex::Vector<amrex::MultiFab *> pfield = amrex::GetVecOfPtrs(field);
		Initialize(lev,pfield);
	};
	void Initialize(const int lev,
			amrex::Vector<amrex::MultiFab> &field)
	{
		amrex::Vector<amrex::MultiFab *> pfield = amrex::GetVecOfPtrs(field);
		Initialize(lev,pfield);
	};
	void Add(const int lev,
			amrex::Vector<amrex::MultiFab> &field)
	{
		amrex::Vector<amrex::MultiFab *> pfield = amrex::GetVecOfPtrs(field);
		Add(lev,pfield);
	};

	void Add(const int lev,
			amrex::Vector<std::unique_ptr<amrex::MultiFab> > &field)
	{
		amrex::Vector<amrex::MultiFab *> pfield = amrex::GetVecOfPtrs(field);
		Add(lev,pfield);
	};

	virtual void SetComp(int a_comp) final {comp = a_comp;}
protected:
	amrex::Vector<amrex::Geometry> &geom;
	int comp = 0;
};
}
#endif

#ifndef BC_BC_H_
#define BC_BC_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>

namespace BC
{
class BC
	: public amrex::PhysBCFunctBase
{

public:
	BC (amrex::Vector<amrex::Geometry> &_geom) : geom(_geom) {}

	virtual void FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real /*time*/) override = 0;

	void SetLevel(int _lev) { lev = _lev; };
	virtual amrex::BCRec GetBCRec() = 0;

	amrex::Array<int,AMREX_SPACEDIM> IsPeriodic()
	{return {{AMREX_D_DECL(geom[0].isPeriodic(0),geom[0].isPeriodic(1),geom[0].isPeriodic(2))}}; }

	/// \todo Not sure how to implement this with a pure virtual class...
	// template<class T>
	// virtual const amrex::Array<amrex::Array<T,AMREX_SPACEDIM>,2> GetBCTypes() = 0;

private:
	std::vector<amrex::Geometry> &geom;
	int lev=0;
};
}
#endif

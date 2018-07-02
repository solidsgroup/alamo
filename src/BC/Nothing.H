#ifndef BC_NOTHING_H_
#define BC_NOTHING_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>

namespace BC
{
class Nothing
	: public BC
{

public:
	Nothing (amrex::Vector<amrex::Geometry> &_geom) : BC(_geom) {};
	virtual void FillBoundary (amrex::MultiFab&, int, int, amrex::Real /*time*/) override {};
	virtual amrex::BCRec GetBCRec() {return amrex::BCRec();}
};
}
#endif

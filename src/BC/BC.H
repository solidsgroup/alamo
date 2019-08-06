#ifndef BC_BC_H_
#define BC_BC_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_Mask.H>
#include <AMReX_Periodicity.H>

#include "Util/Util.H"

/// \brief Collection of boundary condition (BC) objects
namespace BC
{
enum Orientation {
	All = -1,
	AMREX_D_DECL(xlo = 0,
		     ylo = 1,
		     zlo = 2),
	AMREX_D_DECL(xhi = 0 + AMREX_SPACEDIM,
		     yhi = 1 + AMREX_SPACEDIM,
		     zhi = 2 + AMREX_SPACEDIM)};

class BC
	: public amrex::PhysBCFunctBase
{

public:
	virtual ~BC() {};
	
	void define (const amrex::Geometry& a_geom) {m_geom = a_geom;};

	virtual void FillBoundary (amrex::FArrayBox &in,
				   const amrex::Box &box,
				   int ngrow, int dcomp, int ncomp,
				   amrex::Real time,
				   Orientation face = Orientation::All,
				   const amrex::Mask *mask = nullptr) = 0;

	virtual void FillBoundary (amrex::MultiFab& mf,
				   int dcomp, int ncomp,
				   amrex::Real time, int /*bccomp*/) override
	{
		mf.FillBoundary(m_geom.periodicity());
		for (amrex::MFIter mfi(mf,true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.tilebox();
			amrex::FArrayBox &in = mf[mfi];
			FillBoundary(in, box, mf.nGrow(), dcomp, ncomp, time);
		}
	}



	virtual amrex::BCRec GetBCRec() = 0;

	virtual amrex::Array<int,AMREX_SPACEDIM> IsPeriodic()
	{return {{AMREX_D_DECL(m_geom.isPeriodic(0),m_geom.isPeriodic(1),m_geom.isPeriodic(2))}}; }

	virtual amrex::Periodicity Periodicity () const {
		return m_geom.periodicity();
	}
	
	virtual amrex::Periodicity Periodicity (const amrex::Box& b) {
		return m_geom.periodicity(b);
	}


protected:
	amrex::Geometry m_geom;
	//std::vector<amrex::Geometry> &geom;
	//int lev=0;
};

namespace BCUtil
{
int ReadString(std::string bcstring);
bool IsPeriodic(int bctype);
bool IsNeumann(int bctype);
bool IsReflectEven(int bctype);
bool IsReflectOdd(int bctype);
bool IsDirichlet(int bctype);
}

}
#endif

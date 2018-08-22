#ifndef BC_CONSTANT_H_
#define BC_CONSTANT_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Array.H>

#include "BC/BC.H"

namespace BC
{
class Constant
	: public BC
{

public:
	//Constant (amrex::Vector<amrex::Geometry> &_geom) : geom(_geom) {};
	Constant (amrex::Vector<std::string> bc_hi_str,
		  amrex::Vector<std::string> bc_lo_str,
		  AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_lo_1,
			       amrex::Vector<amrex::Real> _bc_lo_2,
			       amrex::Vector<amrex::Real> _bc_lo_3),
		  AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_hi_1,
			       amrex::Vector<amrex::Real> _bc_hi_2,
			       amrex::Vector<amrex::Real> _bc_hi_3));
	virtual void FillBoundary (amrex::FArrayBox &in, const amrex::Box &box,
				   int ngrow, int dcomp, int ncomp, amrex::Real time,
				   amrex::MFIter &mfi,
				   Orientation face = Orientation::All,
				   const amrex::Mask *mask = nullptr) override;

	using BC::FillBoundary;
	
	amrex::BCRec GetBCRec() ;
	virtual amrex::Array<int,AMREX_SPACEDIM> IsPeriodic() override;
	virtual amrex::Periodicity Periodicity () const override;
	virtual amrex::Periodicity Periodicity (const amrex::Box& b) override;



	template<class T>
	const amrex::Array<amrex::Array<T,AMREX_SPACEDIM>,2> GetBCTypes()
	{
		return {{{AMREX_D_DECL((T)bc_lo[0],(T)bc_lo[1],(T)bc_lo[2])},
		 			{AMREX_D_DECL((T)bc_lo[0],(T)bc_lo[1],(T)bc_lo[2])}}};
	}


private:
	int bc_lo[BL_SPACEDIM];
	int bc_hi[BL_SPACEDIM];
	int bc_hi_flux[BL_SPACEDIM];
	int bc_lo_flux[BL_SPACEDIM];
	amrex::Vector<amrex::Real> AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3);
	amrex::Vector<amrex::Real> AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3);
};
}
#endif

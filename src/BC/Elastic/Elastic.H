#ifndef BC_CONSTANT_H_
#define BC_CONSTANT_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Array.H>

#include "BC/BC.H"
#include "Set/Set.H"

namespace BC
{
class Elastic
	: public BC<Set::Scalar>
{

public:
	Elastic (amrex::Vector<std::string> bc_hi_str
		 ,amrex::Vector<std::string> bc_lo_str
		 ,amrex::Vector<amrex::Real> _bc_lo_1
		 ,amrex::Vector<amrex::Real> _bc_hi_1
		 ,amrex::Vector<amrex::Real> _bc_lo_2
		 ,amrex::Vector<amrex::Real> _bc_hi_2
#if AMREX_SPACEDIM > 2
		 ,amrex::Vector<amrex::Real> _bc_lo_3
		 ,amrex::Vector<amrex::Real> _bc_hi_3
#endif
		 );
	virtual void FillBoundary (amrex::FArrayBox &in, const amrex::Box &box,
				   int ngrow, int dcomp, int ncomp, amrex::Real time,
				   Orientation face = Orientation::All,
				   const amrex::Mask *mask = nullptr) override;

	using BC<Set::Scalar>::FillBoundary;
	
	amrex::BCRec GetBCRec() override;
	amrex::Array<int,AMREX_SPACEDIM> IsPeriodic() override;

	template<class T>
	const amrex::Array<amrex::Array<T,AMREX_SPACEDIM>,2> GetBCTypes()
	{
		return {{{AMREX_D_DECL((T)bc_lo[0],(T)bc_lo[1],(T)bc_lo[2])},
		 			{AMREX_D_DECL((T)bc_lo[0],(T)bc_lo[1],(T)bc_lo[2])}}};
	}


private:
	int bc_lo[BL_SPACEDIM];
	int bc_hi[BL_SPACEDIM];
	amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
	amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
#if AMREX_SPACEDIM > 2
	amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;
#endif
};
}
#endif

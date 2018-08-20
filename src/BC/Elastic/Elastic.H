#ifndef BC_ELASTIC_H_
#define BC_ELASTIC_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Array.H>
#include <AMReX.H>

#include "BC/BC.H"
#include "Operator/Operator.H"
#include "Operator/Elastic/Elastic.H"
#include "Set/Set.H"

namespace Operator
{
	namespace Elastic
	{
		class Elastic;
	}
}

namespace BC
{
class Elastic
	: public BC
{

public:
	Elastic ( amrex::Vector<std::string> bc_hi_str,
		  amrex::Vector<std::string> bc_lo_str,
		  AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_lo_1,
			       amrex::Vector<amrex::Real> _bc_lo_2,
			       amrex::Vector<amrex::Real> _bc_lo_3),
		  AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_hi_1,
			       amrex::Vector<amrex::Real> _bc_hi_2,
			       amrex::Vector<amrex::Real> _bc_hi_3));

	virtual void FillBoundary (amrex::MultiFab& mf,
				   amrex::Real time) override;

	using BC::FillBoundary;
	
	amrex::BCRec GetBCRec() ;
	amrex::Array<int,AMREX_SPACEDIM> IsPeriodic();
	void SetElasticOperator(Operator::Elastic::Elastic* a_operator);

	template<class T>
	const amrex::Array<amrex::Array<T,AMREX_SPACEDIM>,2> GetBCTypes()
	{
		return {{{AMREX_D_DECL((T)bc_lo[0],(T)bc_lo[1],(T)bc_lo[2])},
		 			{AMREX_D_DECL((T)bc_lo[0],(T)bc_lo[1],(T)bc_lo[2])}}};
	}
protected:
	void StencilFill(	amrex::Vector<Set::Vector> &stencil,
			const amrex::Vector<Set::Vector> &traction,
			const amrex::Vector<int> &points,
			const amrex::IntVect &m,
			const int amrlev,
			const int mglev,
			const amrex::MFIter &mfi);

private:
	int bc_lo[BL_SPACEDIM];
	int bc_hi[BL_SPACEDIM];
	int bc_hi_flux[BL_SPACEDIM];
	int bc_lo_flux[BL_SPACEDIM];
	amrex::Vector<amrex::Real> AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3);
	amrex::Vector<amrex::Real> AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3);
	Operator::Elastic::Elastic* m_operator;
};
}
#endif

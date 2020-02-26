#ifndef BC_STEP_H_
#define BC_STEP_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Array.H>

#include "Set/Set.H"
#include "BC/BC.H"
#include "BC/Constant.H"
namespace BC
{
class Step
	: public BC
{
public:
Step ()
{
}

void FillBoundary (amrex::FArrayBox &a_in,
			const amrex::Box &a_box,
			int ngrow, int /*dcomp*/, int /*ncomp*/, amrex::Real /*time*/,
			Orientation face, const amrex::Mask * /*mask*/) override
{
	
	const amrex::Real* DX = m_geom.CellSize();

	amrex::Box box = a_box;
	box.grow(ngrow);
	const amrex::Dim3 lo= amrex::lbound(m_geom.Domain()), hi = amrex::ubound(m_geom.Domain());

	amrex::Array4<amrex::Real> const& in = a_in.array();

	amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k)
	{
		amrex::IntVect glevel;
		AMREX_D_TERM(glevel[0] = std::max(std::min(0,i-lo.x),i-hi.x); ,
					 glevel[1] = std::max(std::min(0,j-lo.y),j-hi.y); ,
					 glevel[2] = std::max(std::min(0,k-lo.z),k-hi.z); );

	    Set::Vector x;
	    AMREX_D_TERM(x(0) = m_geom.ProbLo()[0] + ((Set::Scalar)(i) + 0.5) * DX[0];,
		    		 x(1) = m_geom.ProbLo()[1] + ((Set::Scalar)(j) + 0.5) * DX[1];,
			    	 x(2) = m_geom.ProbLo()[2] + ((Set::Scalar)(k) + 0.5) * DX[2];);

		
		if (glevel[0]<0 && (face == Orientation::xlo || face == Orientation::All)) 
		{
	        if (x(1) > m_h1) 
			{
				in(i,j,k,0) = 0.;
				in(i,j,k,1) = 1.;
			}
	        else 
			{
				in(i,j,k,0) = 1.;
				in(i,j,k,1) = 0.;
			}
		}
		else if (glevel[0]>0)
		{
	        if (x(1) > m_h2) 
			{
				in(i,j,k,0) = 0.;
				in(i,j,k,1) = 1.;
			}
	        else 
			{
				in(i,j,k,0) = 1.;
				in(i,j,k,1) = 0.;
			}
		}
		else if (glevel[1]<0) // Bottom boundary
		{
				in(i,j,k,0) = 1.;
				in(i,j,k,1) = 0.;
		}
		else if (glevel[1]>0 && (face == Orientation::yhi || face == Orientation::All)) // Top boundary
		{
				in(i,j,k,0) = 0.;
				in(i,j,k,1) = 1.;
		}
	});

	
}
	amrex::BCRec GetBCRec() override
	{
		int bc_lo[AMREX_SPACEDIM] = {AMREX_D_DECL(INT_DIR,INT_DIR,INT_DIR)};
		int bc_hi[AMREX_SPACEDIM] = {AMREX_D_DECL(INT_DIR,INT_DIR,INT_DIR)};

		return amrex::BCRec(bc_lo,bc_hi);
	}

private:
	Set::Scalar m_h1 = NAN, m_h2 = NAN, m_val1 = NAN, m_val2=NAN;

public:
    static void Parse(Step & value, amrex::ParmParse & pp)
    {
		pp.query("h1",value.m_h1);
		pp.query("h2",value.m_h2);
    }

};
}
#endif

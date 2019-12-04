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
	: public Constant
{

public:
Step (amrex::Vector<std::string> bc_hi_str,
		    amrex::Vector<std::string> bc_lo_str,
		    AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_lo_1,
				 amrex::Vector<amrex::Real> _bc_lo_2,
				 amrex::Vector<amrex::Real> _bc_lo_3),
		    AMREX_D_DECL(amrex::Vector<amrex::Real> _bc_hi_1,
				 amrex::Vector<amrex::Real> _bc_hi_2,
				 amrex::Vector<amrex::Real> _bc_hi_3))
	: Constant(bc_hi_str,bc_lo_str,_bc_lo_1,_bc_lo_2,_bc_hi_1,_bc_hi_2)
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



	for (int n = 0; n < a_in.nComp(); n++)
	amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k)
	{
		amrex::IntVect glevel;
		AMREX_D_TERM(glevel[0] = std::max(std::min(0,i-lo.x),i-hi.x); ,
					 glevel[1] = std::max(std::min(0,j-lo.y),j-hi.y); ,
					 glevel[2] = std::max(std::min(0,k-lo.z),k-hi.z); );
    
        Set::Vector x;
        AMREX_D_TERM(x(0) = m_geom.ProbLo()[0] + ((Set::Scalar)(i) + 0.5) * m_geom.CellSize()[0];,
		    		 x(1) = m_geom.ProbLo()[1] + ((Set::Scalar)(j) + 0.5) * m_geom.CellSize()[1];,
			    	 x(2) = m_geom.ProbLo()[2] + ((Set::Scalar)(k) + 0.5) * m_geom.CellSize()[2];);

		
		if (glevel[0]<0 && (face == Orientation::xlo || face == Orientation::All)) // Left boundary
		{
            if (x(1) > 0.3) in(i,j,k,n) = 0.5;
            else in(i,j,k,n) = -0.5;
		}
		else if (glevel[0]>0)// && (face == Orientation::xhi || face == Orientation::All)) // Right boundary
		{
            if (x(1) > 0.7) in(i,j,k,n) = 0.5;
            else in(i,j,k,n) = -0.5;
		}
		
		else if (glevel[1]<0)// && (face == Orientation::ylo || face == Orientation::All)) // Bottom boundary
		{
            in(i,j,k,n) = -0.5;
		}
		else if (glevel[1]>0 && (face == Orientation::yhi || face == Orientation::All)) // Top boundary
		{
            in(i,j,k,n) = 0.5;
		}
	});
}
//Step (amrex::Vector<amrex::Geometry> &_geom) : geom(_geom) {};
//	Step () {};

//	//virtual ~Step() {};
//
//	virtual void FillBoundary (amrex::FArrayBox &in, const amrex::Box &box,
//				   int ngrow, int dcomp, int ncomp, amrex::Real time,
//				   Orientation face = Orientation::All,
//				   const amrex::Mask *mask = nullptr) override
//    {
//        Util::Message(INFO,"time = ", time);
//    }
//	
//	amrex::BCRec GetBCRec() override;
//	virtual amrex::Array<int,AMREX_SPACEDIM> IsPeriodic() override;
//	virtual amrex::Periodicity Periodicity () const override;
//	virtual amrex::Periodicity Periodicity (const amrex::Box& b) override;
//
//
//
//	template<class T>
//	const amrex::Array<amrex::Array<T,AMREX_SPACEDIM>,2> GetBCTypes()
//	{
//		return {{{AMREX_D_DECL((T)bc_lo[0],(T)bc_lo[1],(T)bc_lo[2])},
//		 			{AMREX_D_DECL((T)bc_lo[0],(T)bc_lo[1],(T)bc_lo[2])}}};
//	}
//
//
//private:
//	int bc_lo[BL_SPACEDIM];
//	int bc_hi[BL_SPACEDIM];
//	amrex::Vector<amrex::Real> AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3);
//	amrex::Vector<amrex::Real> AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3);
};
}
#endif

#include <AMReX_MultiFabUtil.H>
#include "eigen3/Eigen/Core"
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>

#include <AMReX_ArrayLim.H>

#include "Elastic.H"

/// \fn Operator::Elastic::MLPFStiffnessMatrix
///
/// Relay to the define function
/// Also define elastic constants here.
Operator::Elastic::Elastic::Elastic (const Vector<Geometry>& a_geom,
					  const Vector<BoxArray>& a_grids,
					  const Vector<DistributionMapping>& a_dmap,
					  const LPInfo& a_info)
{
  define(a_geom, a_grids, a_dmap, a_info);
}

Operator::Elastic::Elastic::~Elastic ()
{}

void
Operator::Elastic::Elastic::Fapply (int amrlev, ///<[in] AMR Level
			     int mglev,  ///<[in]
			     MultiFab& f,///<[out] The force vector
			     const MultiFab& u ///<[in] The displacements vector
			     ) const
{
  const Real* dx = m_geom[amrlev][mglev].CellSize();

  for (MFIter mfi(f, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const amrex::BaseFab<amrex::Real> &ufab  = u[mfi];
      amrex::BaseFab<amrex::Real>       &ffab  = f[mfi];
      
      for (int m = bx.loVect()[0]; m<=bx.hiVect()[0]; m++)
	for (int n = bx.loVect()[1]; n<=bx.hiVect()[1]; n++)
	  {
	    for (int i=0; i<AMREX_SPACEDIM; i++)
	      {
		ffab(amrex::IntVect(m,n),i) = 0.0;
		for (int k=0; k<AMREX_SPACEDIM; k++)
		  {
		    // C_{ijkl} u_{k,jl}

		    ffab(amrex::IntVect(m,n),i) -= 
		      C(i,0,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi)
		      *(ufab(amrex::IntVect(m+1,n),k) - 2.0*ufab(amrex::IntVect(m,n),k) + ufab(amrex::IntVect(m-1,n),k))/dx[0]/dx[0] 
		      + 
		      (C(i,0,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi) + C(i,1,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi))
		      *(ufab(amrex::IntVect(m+1,n+1),k) + ufab(amrex::IntVect(m-1,n-1),k) - ufab(amrex::IntVect(m+1,n-1),k) - ufab(amrex::IntVect(m-1,n+1),k))/(2.0*dx[0])/(2.0*dx[1])
		      +
		      C(i,1,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi)
		      *(ufab(amrex::IntVect(m,n+1),k) - 2.0*ufab(amrex::IntVect(m,n),k) + ufab(amrex::IntVect(m,n-1),k))/dx[1]/dx[1];

		    // C_{ijkl,j} u_{k,l}

		    ffab(amrex::IntVect(m,n),i) -= 
		      ((C(i,0,k,0,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,0,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(2.0*dx[0]) +
		       (C(i,1,k,0,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,0,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(2.0*dx[1])) *
		      ((ufab(amrex::IntVect(m+1,n),k) - ufab(amrex::IntVect(m-1,n),k))/(2.0*dx[0]))
		      +
		      ((C(i,0,k,1,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,1,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(2.0*dx[0]) +
		       (C(i,1,k,1,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,1,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(2.0*dx[1])) *
		      ((ufab(amrex::IntVect(m,n+1),k) - ufab(amrex::IntVect(m,n-1),k))/(2.0*dx[1]));
		      
		  }
	      }
	  }
    }
}


/// \fn Operator::Elastic::Fsmooth
/// 
/// Perform one half Gauss-Seidel iteration corresponding to the operator specified
/// in Operator::Elastic::Fapply.
/// The variable redblack corresponds to whether to smooth "red" nodes or "black"
/// nodes, where red and black nodes are distributed in a checkerboard pattern.
/// 
/// \todo Extend to 3D
///
void
Operator::Elastic::Elastic::Fsmooth (int amrlev,          ///<[in] AMR level
			      int mglev,           ///<[in] 
			      MultiFab& u,       ///<[inout] Solution (displacement field)
			      const MultiFab& rhs, ///<[in] Body force vectors (rhs=right hand side)
			      int redblack         ///<[in] Smooth even vs. odd modes
			      ) const
{
  BL_PROFILE("Operator::Elastic::Fsmooth()");

  const Real* dx = m_geom[amrlev][mglev].CellSize();

  for (MFIter mfi(u,MFItInfo().EnableTiling().SetDynamic(true));
       mfi.isValid(); ++mfi)
    {
      const Box&       tbx     = mfi.tilebox();
      FArrayBox&       ufab    = u[mfi];
      const FArrayBox& rhsfab  = rhs[mfi];

      for (int n = tbx.loVect()[1]; n<=tbx.hiVect()[1]; n++)
	{
	  int noffset = (tbx.loVect()[0] + n + redblack)%2;
	  for (int m = tbx.loVect()[0] + noffset; m <= tbx.hiVect()[0]; m+= 2)
	    {
	      for (int i=0; i<AMREX_SPACEDIM; i++)
		{
		  amrex::Real rho = 0.0, aa = 0.0;
		      
		  for (int k=0; k<AMREX_SPACEDIM; k++)
		    {
		      
		      rho -= 
			C(i,0,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi)
			*(ufab(amrex::IntVect(m+1,n),k) - (i==k ? 0.0 : 2.0*ufab(amrex::IntVect(m,n),k)) + ufab(amrex::IntVect(m-1,n),k))/dx[0]/dx[0]
			+ 
			(C(i,0,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi) + C(i,1,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi))
			*(ufab(amrex::IntVect(m+1,n+1),k) + ufab(amrex::IntVect(m-1,n-1),k) - ufab(amrex::IntVect(m+1,n-1),k) - ufab(amrex::IntVect(m-1,n+1),k))/(2.0*dx[0])/(2.0*dx[1])
			+
			C(i,1,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi)
			*(ufab(amrex::IntVect(m,n+1),k) - (i==k ? 0.0 : 2.0*ufab(amrex::IntVect(m,n),k)) + ufab(amrex::IntVect(m,n-1),k))/dx[1]/dx[1];

		      // C_{ijkl,j} u_{k,l}

		      rho -= 
		        ((C(i,0,k,0,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,0,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(2.0*dx[0]) +
		         (C(i,1,k,0,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,0,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(2.0*dx[1])) *
		        ((ufab(amrex::IntVect(m+1,n),k) - ufab(amrex::IntVect(m-1,n),k))/(2.0*dx[0]))
		        +
		        ((C(i,0,k,1,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,1,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(2.0*dx[0]) +
		         (C(i,1,k,1,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,1,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(2.0*dx[1])) *
		        ((ufab(amrex::IntVect(m,n+1),k) - ufab(amrex::IntVect(m,n-1),k))/(2.0*dx[1]));
		      
		    }

		  aa -=
		    -2.0*C(i,0,i,0,amrex::IntVect(m,n),amrlev,mglev,mfi)/dx[0]/dx[0]
		    -2.0*C(i,1,i,1,amrex::IntVect(m,n),amrlev,mglev,mfi)/dx[1]/dx[1];

		  //std::cout << "nans not detetected, rho=" << rho << ", aa=" << aa << std::endl;
		  if (rho != rho) std::cout << "nans detetected, rho=" << rho << ", aa=" << aa << std::endl;
		  if (rho != rho) amrex::Abort("nans detected");
		  
		  ufab(amrex::IntVect(m,n),i) = (rhsfab(amrex::IntVect(m,n),i) - rho) / aa;
		}
	    }
	}
    }
}

/// \fn Operator::Elastic::FFlux
/// 
/// Compute the "flux" corresponding to the operator in Operator::Elastic::Fapply.
/// Because the operator is self-adjoint and positive-definite, the flux is not
/// required for adequate convergence (?)
/// Therefore, the fluxes are simply set to zero and returned.
/// 
/// \todo Extend to 3D
///
void
Operator::Elastic::Elastic::FFlux (int /*amrlev*/, const MFIter& /*mfi*/,
			  const std::array<FArrayBox*,AMREX_SPACEDIM>& sigmafab,
			  const FArrayBox& /*ufab*/, const int /*face_only*/) const
{

  amrex::BaseFab<amrex::Real> &fxfab = *sigmafab[0];
  amrex::BaseFab<amrex::Real> &fyfab = *sigmafab[1];
  fxfab.setVal(0.0);
  fyfab.setVal(0.0);
}




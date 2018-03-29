#include <AMReX_MultiFabUtil.H>
#include "MLStiffnessMatrix.H"
#include "eigen3/Eigen/Core"
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>

#include <AMReX_ArrayLim.H>

/// \fn MLStiffnessMatrix::MLStiffnessMatrix
///
/// Relay to the define function
/// Also define elastic constants here.
MLStiffnessMatrix::MLStiffnessMatrix (const Vector<Geometry>& a_geom,
				      const Vector<BoxArray>& a_grids,
				      const Vector<DistributionMapping>& a_dmap,
				      const LPInfo& a_info)
{
  define(a_geom, a_grids, a_dmap, a_info);

  // amrex::Real E1 = 1.0, E2=2.0;
  // amrex::Real nu1 = 0.33, nu2=0.33;//0.33;
  
  // mu1 = E1*nu1 / (1.0 - 2.0*nu1) / (1.0 + nu1);
  // lambda1 = E1 / 2.0 / (1.0 + nu1);

  // mu2 = E2*nu2 / (1.0 - 2.0*nu2) / (1.0 + nu2);
  // lambda2 = E2 / 2.0 / (1.0 + nu2);

  mu1 = 1.0; mu2 = 2.0;
  lambda1 = 1.0; lambda2 = 1.0;

}

/// \fn define
///
/// Currently this is just a relay to call the parent class' define function
/// However if there are operator-specific Fabs to be created, this is the place
/// to do it.
void
MLStiffnessMatrix::define (const Vector<Geometry>& a_geom,
			   const Vector<BoxArray>& a_grids,
			   const Vector<DistributionMapping>& a_dmap,
			   const LPInfo& a_info)
{
  MLCellLinOp::define(a_geom, a_grids, a_dmap, a_info);
}

MLStiffnessMatrix::~MLStiffnessMatrix ()
{}

/// \fn prepareForSolve
///
/// Currently just a relay function. It exists because it overrides a pure
/// virtual function.
void
MLStiffnessMatrix::prepareForSolve ()
{
  MLCellLinOp::prepareForSolve();
}

/// \fn MLStiffnessMatrix::Fapply
/// 
/// Numerically implement the operator
///
/// \f[ f_i = \mathbb{C}_{ijkl} u_{k,jl}\f]
///
/// For each node, a strain gradient variable is constructed containing all second
/// partial derivatives of the displacement.
///
/// \f[\text{gradepsilon[i](j,k)} = \frac{\partial^2 u_i}{\partial x_j\partial x_k}\f]
///
/// The force vector is then computed by
///
/// \f[f_p = 2\mu\varepsilon_{pq,q} + \lambda\delta_{pq}\varepsilon_{kk,q}\f]
///
/// where \f$\mu\f$ = \link #mu \endlink member variable and 
///
/// \f$\lambda\f$ = \link #lambda \endlink member variable.
///
/// \note This implementation is for an isotropic material with uniform properties.
///       and will not work for materials with properties that vary in space
///
void
MLStiffnessMatrix::Fapply (int amrlev, ///<[in] AMR Level
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
		    //ffab(amrex::IntVect(m,n),i) = 0.0;

		    // C_{ijkl} u_{k,jl}

		    ffab(amrex::IntVect(m,n),i) -= 
		      C(i,0,k,0,amrex::IntVect(m,n),amrlev,mglev)
		      *(ufab(amrex::IntVect(m+1,n),k) - 2.0*ufab(amrex::IntVect(m,n),k) + ufab(amrex::IntVect(m-1,n),k))/dx[0]/dx[0] 
		      + 
		      (C(i,0,k,1,amrex::IntVect(m,n),amrlev,mglev) + C(i,1,k,0,amrex::IntVect(m,n),amrlev,mglev))
		      *(ufab(amrex::IntVect(m+1,n+1),k) + ufab(amrex::IntVect(m-1,n-1),k) - ufab(amrex::IntVect(m+1,n-1),k) - ufab(amrex::IntVect(m-1,n+1),k))/(2.0*dx[0])/(2.0*dx[1])
		      +
		      C(i,1,k,1,amrex::IntVect(m,n),amrlev,mglev)
		      *(ufab(amrex::IntVect(m,n+1),k) - 2.0*ufab(amrex::IntVect(m,n),k) + ufab(amrex::IntVect(m,n-1),k))/dx[1]/dx[1];

		    // C_{ijkl,j} u_{k,l}

		    ffab(amrex::IntVect(m,n),i) -= 
		      ((C(i,0,k,0,amrex::IntVect(m+1,n),amrlev,mglev) - C(i,0,k,0,amrex::IntVect(m-1,n),amrlev,mglev))/(2.0*dx[0]) +
		       (C(i,1,k,0,amrex::IntVect(m,n+1),amrlev,mglev) - C(i,1,k,0,amrex::IntVect(m,n-1),amrlev,mglev))/(2.0*dx[1])) *
		      ((ufab(amrex::IntVect(m+1,n),k) - ufab(amrex::IntVect(m-1,n),k))/(2.0*dx[0]))
		      +
		      ((C(i,0,k,1,amrex::IntVect(m+1,n),amrlev,mglev) - C(i,0,k,1,amrex::IntVect(m-1,n),amrlev,mglev))/(2.0*dx[0]) +
		       (C(i,1,k,1,amrex::IntVect(m,n+1),amrlev,mglev) - C(i,1,k,1,amrex::IntVect(m,n-1),amrlev,mglev))/(2.0*dx[1])) *
		      ((ufab(amrex::IntVect(m,n+1),k) - ufab(amrex::IntVect(m,n-1),k))/(2.0*dx[1]));
		      
		  }
	      }
	  }
    }
}


/// \fn MLStiffnessMatrix::Fsmooth
/// 
/// Perform one half Gauss-Seidel iteration corresponding to the operator specified
/// in MLStiffnessMatrix::Fapply.
/// The variable redblack corresponds to whether to smooth "red" nodes or "black"
/// nodes, where red and black nodes are distributed in a checkerboard pattern.
/// 
/// \todo Extend to 3D
///
void
MLStiffnessMatrix::Fsmooth (int amrlev,  ///<[in] AMR level
			    int mglev, ///<[in] 
			    MultiFab& u, ///<[inout] Solution (displacement field)
			    const MultiFab& rhs, ///<[in] Body force vectors (rhs=right hand side)
			    int redblack ///<[in] Variable to determine whether to smooth even or odd modes
			    ) const
{
  BL_PROFILE("MLStiffnessMatrix::Fsmooth()");

  const int nComp = AMREX_SPACEDIM;
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
			C(i,0,k,0,amrex::IntVect(m,n),amrlev,mglev)
			*(ufab(amrex::IntVect(m+1,n),k) - (i==k ? 0.0 : 2.0*ufab(amrex::IntVect(m,n),k)) + ufab(amrex::IntVect(m-1,n),k))/dx[0]/dx[0]
			+ 
			(C(i,0,k,1,amrex::IntVect(m,n),amrlev,mglev) + C(i,1,k,0,amrex::IntVect(m,n),amrlev,mglev))
			*(ufab(amrex::IntVect(m+1,n+1),k) + ufab(amrex::IntVect(m-1,n-1),k) - ufab(amrex::IntVect(m+1,n-1),k) - ufab(amrex::IntVect(m-1,n+1),k))/(2.0*dx[0])/(2.0*dx[1])
			+
			C(i,1,k,1,amrex::IntVect(m,n),amrlev,mglev)
			*(ufab(amrex::IntVect(m,n+1),k) - (i==k ? 0.0 : 2.0*ufab(amrex::IntVect(m,n),k)) + ufab(amrex::IntVect(m,n-1),k))/dx[1]/dx[1];

		      // C_{ijkl,j} u_{k,l}

		      rho -= 
		        ((C(i,0,k,0,amrex::IntVect(m+1,n),amrlev,mglev) - C(i,0,k,0,amrex::IntVect(m-1,n),amrlev,mglev))/(2.0*dx[0]) +
		         (C(i,1,k,0,amrex::IntVect(m,n+1),amrlev,mglev) - C(i,1,k,0,amrex::IntVect(m,n-1),amrlev,mglev))/(2.0*dx[1])) *
		        ((ufab(amrex::IntVect(m+1,n),k) - ufab(amrex::IntVect(m-1,n),k))/(2.0*dx[0]))
		        +
		        ((C(i,0,k,1,amrex::IntVect(m+1,n),amrlev,mglev) - C(i,0,k,1,amrex::IntVect(m-1,n),amrlev,mglev))/(2.0*dx[0]) +
		         (C(i,1,k,1,amrex::IntVect(m,n+1),amrlev,mglev) - C(i,1,k,1,amrex::IntVect(m,n-1),amrlev,mglev))/(2.0*dx[1])) *
		        ((ufab(amrex::IntVect(m,n+1),k) - ufab(amrex::IntVect(m,n-1),k))/(2.0*dx[1]));
		      
		    }

		  aa -= 
		    -2.0*C(i,0,i,0,amrex::IntVect(m,n),amrlev,mglev)/dx[0]/dx[0]
		    -2.0*C(i,1,i,1,amrex::IntVect(m,n),amrlev,mglev)/dx[1]/dx[1];

		  if (rho != rho) amrex::Abort("nans detected");
		  
		  ufab(amrex::IntVect(m,n),i) = (rhsfab(amrex::IntVect(m,n),i) - rho) / aa;
		}
	    }
	}

    }
}

/// \fn MLStiffnessMatrix::Fsmooth
/// 
/// Compute the "flux" corresponding to the operator in MLStiffnessMatrix::Fapply.
/// Because the operator is self-adjoint and positive-definite, the flux is not
/// required for adequate convergence.
/// Therefore, the fluxes are simply set to zero and returned.
/// 
/// \todo Extend to 3D
///
void
MLStiffnessMatrix::FFlux (int /*amrlev*/, const MFIter& /*mfi*/,
 			  const std::array<FArrayBox*,AMREX_SPACEDIM>& sigma,
 			  const FArrayBox& /*sol*/, const int /*face_only*/) const
{

  amrex::BaseFab<amrex::Real> &fxfab = *sigma[0];
  amrex::BaseFab<amrex::Real> &fyfab = *sigma[1];

  fxfab.setVal(0.0);
  fyfab.setVal(0.0);
}

int
MLStiffnessMatrix::getNComp() const
{
  return AMREX_SPACEDIM;
}


amrex::Real
MLStiffnessMatrix::C(const int i, const int j, const int k, const int l,
		     const amrex::IntVect loc,
		     int amrlev, int mglev) const
{
  amrex::Real mu, lambda;
  amrex::Real x = m_geom[amrlev][mglev].ProbLo()[0] + ((amrex::Real)(loc[0]) + 0.5) * m_geom[amrlev][mglev].CellSize()[0];
  amrex::Real y = m_geom[amrlev][mglev].ProbLo()[1] + ((amrex::Real)(loc[1]) + 0.5) * m_geom[amrlev][mglev].CellSize()[1];
  //mu = mu1 + y*(mu2-mu1);
  //lambda = lambda1 + y*(lambda2-lambda1);
  if (y>0.5) {mu=mu1; lambda=lambda1; }
  else {mu=mu2; lambda=lambda2;}

  amrex::Real ret = 0.0;
  if (i==k && j==l) ret += mu;
  if (i==l && j==k) ret += mu;
  if (i==j && k==l) ret += lambda;
};



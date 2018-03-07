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

  amrex::Real E = 68;
  amrex::Real nu = 0.45;//0.33;
  
  mu = E*nu / (1.0 - 2.0*nu) / (1.0 + nu);
  lambda = E / 2.0 / (1.0 + nu);
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

      for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
	for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
	  {

	    amrex::Array<Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,AMREX_SPACEDIM> > gradepsilon(AMREX_SPACEDIM);
	    
	    for (int n = 0; n < AMREX_SPACEDIM; n++)
	      gradepsilon[n] <<
		(ufab(amrex::IntVect(i+1,j),n) + ufab(amrex::IntVect(i-1,j),n) - 2.*ufab(amrex::IntVect(i,j),n))/dx[0]/dx[0],
		(ufab(amrex::IntVect(i+1,j+1),n) + ufab(amrex::IntVect(i-1,j-1),n) - ufab(amrex::IntVect(i+1,j-1),n) - ufab(amrex::IntVect(i-1,j+1),n))/(2.*dx[0])/(2.*dx[1]),
		(ufab(amrex::IntVect(i+1,j+1),n) + ufab(amrex::IntVect(i-1,j-1),n) - ufab(amrex::IntVect(i+1,j-1),n) - ufab(amrex::IntVect(i-1,j+1),n))/(2.*dx[0])/(2.*dx[1]),
		(ufab(amrex::IntVect(i,j+1),n) + ufab(amrex::IntVect(i,j-1),n) - 2.*ufab(amrex::IntVect(i,j),n))/dx[1]/dx[1];

	    for (int p = 0; p < AMREX_SPACEDIM; p++)
	      {
		ffab(amrex::IntVect(i,j),p) = 0.0;
		for (int q = 0; q < AMREX_SPACEDIM; q++)
		  {
		    ffab(amrex::IntVect(i,j),p) -= 2.0*mu*gradepsilon[p](q,q);
		    if (p==q)
		      for (int k=0; k<AMREX_SPACEDIM; k++)
		      	ffab(amrex::IntVect(i,j),p) -= lambda * gradepsilon[k](k,q);
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
			    MultiFab& sol, ///<[inout] Solution (displacement field)
			    const MultiFab& rhs, ///<[in] Body force vectors (rhs=right hand side)
			    int redblack ///<[in] Variable to determine whether to smooth even or odd modes
			    ) const
{
  BL_PROFILE("MLStiffnessMatrix::Fsmooth()");

  const int nComp = AMREX_SPACEDIM;
  const Real* dx = m_geom[amrlev][mglev].CellSize();

  for (MFIter mfi(sol,MFItInfo().EnableTiling().SetDynamic(true));
       mfi.isValid(); ++mfi)
    {
      const Box&       tbx     = mfi.tilebox();
      FArrayBox&       solnfab = sol[mfi];
      const FArrayBox& rhsfab  = rhs[mfi];

      for (int n = 0; n < nComp; n++)
      	for (int j = tbx.loVect()[1]; j<=tbx.hiVect()[1]; j++)
       	  {
       	    int ioffset = (tbx.loVect()[0] + j + redblack)%2;
       	    for (int i = tbx.loVect()[0] + ioffset; i <= tbx.hiVect()[0]; i+= 2)
       	      {
		amrex::Array<Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,AMREX_SPACEDIM> > gradepsilonD(AMREX_SPACEDIM);
		amrex::Array<Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,AMREX_SPACEDIM> > gradepsilonR(AMREX_SPACEDIM);
	    
		for (int p = 0; p < AMREX_SPACEDIM; p++)
		  if (p==n)
		    {
		      gradepsilonR[p] <<
			(solnfab(amrex::IntVect(i+1,j),p) + solnfab(amrex::IntVect(i-1,j),p))/dx[0]/dx[0],
			(solnfab(amrex::IntVect(i+1,j+1),p) + solnfab(amrex::IntVect(i-1,j-1),p) - solnfab(amrex::IntVect(i+1,j-1),p) - solnfab(amrex::IntVect(i-1,j+1),p))/(2.*dx[0])/(2.*dx[1]),
			(solnfab(amrex::IntVect(i+1,j+1),p) + solnfab(amrex::IntVect(i-1,j-1),p) - solnfab(amrex::IntVect(i+1,j-1),p) - solnfab(amrex::IntVect(i-1,j+1),p))/(2.*dx[0])/(2.*dx[1]),
			(solnfab(amrex::IntVect(i,j+1),p) + solnfab(amrex::IntVect(i,j-1),p))/dx[1]/dx[1];
		
		      gradepsilonD[p] <<
			-2.0/dx[0]/dx[0],
			0.0,
			0.0,
			-2.0/dx[1]/dx[1];
		    }
		  else
		    {
		      gradepsilonR[p] <<
			(solnfab(amrex::IntVect(i+1,j),p) + solnfab(amrex::IntVect(i-1,j),p) - 2.*solnfab(amrex::IntVect(i,j),p))/dx[0]/dx[0],
			(solnfab(amrex::IntVect(i+1,j+1),p) + solnfab(amrex::IntVect(i-1,j-1),p) - solnfab(amrex::IntVect(i+1,j-1),p) - solnfab(amrex::IntVect(i-1,j+1),p))/(2.*dx[0])/(2.*dx[1]),
			(solnfab(amrex::IntVect(i+1,j+1),p) + solnfab(amrex::IntVect(i-1,j-1),p) - solnfab(amrex::IntVect(i+1,j-1),p) - solnfab(amrex::IntVect(i-1,j+1),p))/(2.*dx[0])/(2.*dx[1]),
			(solnfab(amrex::IntVect(i,j+1),p) + solnfab(amrex::IntVect(i,j-1),p) - 2.*solnfab(amrex::IntVect(i,j),p))/dx[1]/dx[1];
		
		      gradepsilonD[p] <<
			0.0,
			0.0,
			0.0,
			0.0;
		    }
		
		amrex::Real rho = 0.0;

		for (int q = 0; q < AMREX_SPACEDIM; q++)
		  {
		    rho -= 2.0*mu*gradepsilonR[n](q,q);
		    if (n==q)
		      for (int k=0; k<AMREX_SPACEDIM; k++)
			rho -= lambda * gradepsilonR[k](k,q);
		  }

		amrex::Real aa = 0.0;
		for (int q = 0; q < AMREX_SPACEDIM; q++)
		  {
		    aa -= 2.0*mu*gradepsilonD[n](q,q);
		    if (n==q)
		      for (int k=0; k<AMREX_SPACEDIM; k++)
			aa -= lambda * gradepsilonD[k](k,q);
		  }

		solnfab(amrex::IntVect(i,j),n) = (rhsfab(amrex::IntVect(i,j),n) - rho) / aa;

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



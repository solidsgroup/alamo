#include <AMReX_MultiFabUtil.H>
#include "MLPFStiffnessMatrix.H"
#include "eigen3/Eigen/Core"
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>

#include <AMReX_ArrayLim.H>

/// \fn MLPFStiffnessMatrix::MLPFStiffnessMatrix
///
/// Relay to the define function
/// Also define elastic constants here.
MLPFStiffnessMatrix::MLPFStiffnessMatrix (const Vector<Geometry>& a_geom,
					  const Vector<BoxArray>& a_grids,
					  const Vector<DistributionMapping>& a_dmap,
					  const LPInfo& a_info)
{
  define(a_geom, a_grids, a_dmap, a_info);

  // amrex::Real E1 = 1.0, E2 = 1.0;
  // amrex::Real nu1 = 0.1, nu2 = 0.33;
  // mu1 = E1*nu1 / (1.0 - 2.0*nu1) / (1.0 + nu1);
  // lambda1 = E1 / 2.0 / (1.0 + nu1);
  // mu2 = E2*nu2 / (1.0 - 2.0*nu2) / (1.0 + nu2);
  // lambda2 = E2 / 2.0 / (1.0 + nu2);

  mu1 = 1.0; mu2=2.0;
  lambda1 = 2.0; lambda2=1.0;
   
  // amrex::Real G1 = 2.0, G2 = 1.0;
  // amrex::Real nu1 = 0.33, nu2 = 0.33;
  // mu1 = G1;
  // lambda1 = 2.0*G1*nu1/(1.0 - 2.0*nu1);
  // mu2 = G2;
  // lambda2 = 2.0*G2*nu2/(1.0 - 2.0*nu2);
}


/// \fn define
///
/// (This documentation is depricated)
/// However if there are operator-specific Fabs to be created, this is the place
/// to do it.
void
MLPFStiffnessMatrix::define (const Vector<Geometry>& a_geom,
			   const Vector<BoxArray>& a_grids,
			   const Vector<DistributionMapping>& a_dmap,
			   const LPInfo& a_info)
{
  MLCellLinOp::define(a_geom, a_grids, a_dmap, a_info);
  m_a_coeffs.resize(m_num_amr_levels);
  for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
      m_a_coeffs[amrlev].resize(m_num_mg_levels[amrlev]);
      for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
	/// \todo MAKE THIS MORE VERSATILE FOR MULTI-GRAIN
	m_a_coeffs[amrlev][mglev].define(m_grids[amrlev][mglev],m_dmap[amrlev][mglev],2,0);
    }
}

void
MLPFStiffnessMatrix::setACoeffs(int amrlev,
				const MultiFab& alpha)
{
  MultiFab::Copy(m_a_coeffs[amrlev][0], alpha, 0, 0, alpha.nComp(), 0);
}
void
MLPFStiffnessMatrix::averageDownCoeffs ()
{

  for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
    {
      auto& fine_a_coeffs = m_a_coeffs[amrlev];

      averageDownCoeffsSameAmrLevel(fine_a_coeffs);
    }

  averageDownCoeffsSameAmrLevel(m_a_coeffs[0]);
}
void
MLPFStiffnessMatrix::averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a)
{
  int nmglevs = a.size();
  for (int mglev = 1; mglev < nmglevs; ++mglev)
    {
      amrex::average_down(a[mglev-1], a[mglev], 0, a[0].nComp(), mg_coarsen_ratio);
    }
}
void
MLPFStiffnessMatrix::applyMetricTermsCoeffs ()
{
#if (AMREX_SPACEDIM != 3)
    for (int alev = 0; alev < m_num_amr_levels; ++alev)
    {
        const int mglev = 0;
        applyMetricTerm(alev, mglev, m_a_coeffs[alev][mglev]);
    }
#endif
}

MLPFStiffnessMatrix::~MLPFStiffnessMatrix ()
{}

/// \fn prepareForSolve
///
/// Relay function and distribute coefficients to all MG levels (?)
void
MLPFStiffnessMatrix::prepareForSolve ()
{
  MLCellLinOp::prepareForSolve();
  applyMetricTermsCoeffs();
  averageDownCoeffs();
}

/// \fn MLPFStiffnessMatrix::Fapply
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
MLPFStiffnessMatrix::Fapply (int amrlev, ///<[in] AMR Level
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
      const amrex::BaseFab<amrex::Real> &etafab = m_a_coeffs[amrlev][mglev][mfi];

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

	    // amrex::Real x = m_geom[amrlev][mglev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * m_geom[amrlev][mglev].CellSize()[0];
	    // amrex::Real y = m_geom[amrlev][mglev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * m_geom[amrlev][mglev].CellSize()[1];
	    // amrex::Real mu,lambda;
	    // if (x>2.5) {mu = mu1; lambda=lambda1;}
	    // else {mu = mu2; lambda=lambda2;}
	    amrex::Real mu = (mu1*etafab(amrex::IntVect(i,j),0) + mu2*etafab(amrex::IntVect(i,j),1))/(etafab(amrex::IntVect(i,j),0) + etafab(amrex::IntVect(i,j),1));
	    amrex::Real lambda = (lambda1*etafab(amrex::IntVect(i,j),0) + lambda2*etafab(amrex::IntVect(i,j),1))/(etafab(amrex::IntVect(i,j),0) + etafab(amrex::IntVect(i,j),1));

	    //if (mglev>2) std::cout << "mu = " << mu << " lambda = " << lambda << std::endl;
	    // ffab(amrex::IntVect(i,j),0) = mu;
	    // ffab(amrex::IntVect(i,j),1) = lambda;

	    for (int p = 0; p < AMREX_SPACEDIM; p++)
	      {
	    	ffab(amrex::IntVect(i,j),p) = 0.0;
	    	for (int q = 0; q < AMREX_SPACEDIM; q++)
	    	  {
	    	    ffab(amrex::IntVect(i,j),p) -= 2.0*mu*gradepsilon[p](q,q);
	    	    for (int k=0; k<AMREX_SPACEDIM; k++)
	    	      ffab(amrex::IntVect(i,j),p) -= lambda * gradepsilon[k](k,p);
	    	  }
	      }	    
	  }
    }
}


/// \fn MLPFStiffnessMatrix::Fsmooth
/// 
/// Perform one half Gauss-Seidel iteration corresponding to the operator specified
/// in MLPFStiffnessMatrix::Fapply.
/// The variable redblack corresponds to whether to smooth "red" nodes or "black"
/// nodes, where red and black nodes are distributed in a checkerboard pattern.
/// 
/// \todo Extend to 3D
///
void
MLPFStiffnessMatrix::Fsmooth (int amrlev,          ///<[in] AMR level
			      int mglev,           ///<[in] 
			      MultiFab& sol,       ///<[inout] Solution (displacement field)
			      const MultiFab& rhs, ///<[in] Body force vectors (rhs=right hand side)
			      int redblack         ///<[in] Variable to determine whether to smooth even or odd modes
			      ) const
{
  BL_PROFILE("MLPFStiffnessMatrix::Fsmooth()");

  const int nComp = AMREX_SPACEDIM;
  const Real* dx = m_geom[amrlev][mglev].CellSize();

  for (MFIter mfi(sol,MFItInfo().EnableTiling().SetDynamic(true));
       mfi.isValid(); ++mfi)
    {
      const Box&       tbx     = mfi.tilebox();
      FArrayBox&       solnfab = sol[mfi];
      const FArrayBox& rhsfab  = rhs[mfi];
      const amrex::BaseFab<amrex::Real> &etafab = m_a_coeffs[amrlev][mglev][mfi];

      	for (int j = tbx.loVect()[1]; j<=tbx.hiVect()[1]; j++)
       	  {
       	    int ioffset = (tbx.loVect()[0] + j + redblack)%2;
       	    for (int i = tbx.loVect()[0] + ioffset; i <= tbx.hiVect()[0]; i+= 2)
       	      {

		amrex::Array<Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,AMREX_SPACEDIM> > gradepsilonD(AMREX_SPACEDIM);
		amrex::Array<Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,AMREX_SPACEDIM> > gradepsilonR(AMREX_SPACEDIM);

		// amrex::Real x = m_geom[amrlev][mglev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * m_geom[amrlev][mglev].CellSize()[0];
		// amrex::Real y = m_geom[amrlev][mglev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * m_geom[amrlev][mglev].CellSize()[1];
		// amrex::Real mu,lambda;
		// if (x>2.5) {mu = mu1; lambda=lambda1;}
		// else {mu = mu2; lambda=lambda2;}

		amrex::Real mu = (mu1*etafab(amrex::IntVect(i,j),0) + mu2*etafab(amrex::IntVect(i,j),1))/(etafab(amrex::IntVect(i,j),0) + etafab(amrex::IntVect(i,j),1));
		amrex::Real lambda = (lambda1*etafab(amrex::IntVect(i,j),0) + lambda2*etafab(amrex::IntVect(i,j),1))/(etafab(amrex::IntVect(i,j),0) + etafab(amrex::IntVect(i,j),1));

		for (int p = 0; p < AMREX_SPACEDIM; p++)
		  {
		    for (int q = 0; q < AMREX_SPACEDIM; q++)
		      if (q==p)
			{
			  gradepsilonR[q] <<
			    (solnfab(amrex::IntVect(i+1,j),q) + solnfab(amrex::IntVect(i-1,j),q))/dx[0]/dx[0],
			    (solnfab(amrex::IntVect(i+1,j+1),q) + solnfab(amrex::IntVect(i-1,j-1),q) - solnfab(amrex::IntVect(i+1,j-1),q) - solnfab(amrex::IntVect(i-1,j+1),q))/(2.*dx[0])/(2.*dx[1]),
			    (solnfab(amrex::IntVect(i+1,j+1),q) + solnfab(amrex::IntVect(i-1,j-1),q) - solnfab(amrex::IntVect(i+1,j-1),q) - solnfab(amrex::IntVect(i-1,j+1),q))/(2.*dx[0])/(2.*dx[1]),
			    (solnfab(amrex::IntVect(i,j+1),q) + solnfab(amrex::IntVect(i,j-1),q))/dx[1]/dx[1];
		
			  gradepsilonD[q] <<
			    -2.0/dx[0]/dx[0],
			    0.0,
			    0.0,
			    -2.0/dx[1]/dx[1];
			}
		      else
			{
			  gradepsilonR[q] <<
			    (solnfab(amrex::IntVect(i+1,j),q) + solnfab(amrex::IntVect(i-1,j),q) - 2.*solnfab(amrex::IntVect(i,j),q))/dx[0]/dx[0],
			    (solnfab(amrex::IntVect(i+1,j+1),q) + solnfab(amrex::IntVect(i-1,j-1),q) - solnfab(amrex::IntVect(i+1,j-1),q) - solnfab(amrex::IntVect(i-1,j+1),q))/(2.*dx[0])/(2.*dx[1]),
			    (solnfab(amrex::IntVect(i+1,j+1),q) + solnfab(amrex::IntVect(i-1,j-1),q) - solnfab(amrex::IntVect(i+1,j-1),q) - solnfab(amrex::IntVect(i-1,j+1),q))/(2.*dx[0])/(2.*dx[1]),
			    (solnfab(amrex::IntVect(i,j+1),q) + solnfab(amrex::IntVect(i,j-1),q) - 2.*solnfab(amrex::IntVect(i,j),q))/dx[1]/dx[1];
		
			  gradepsilonD[q] <<
			    0.0,
			    0.0,
			    0.0,
			    0.0;
			}
		
		    amrex::Real rho = 0.0;

		    for (int q = 0; q < AMREX_SPACEDIM; q++)
		      {
			rho -= 2.0*mu*gradepsilonR[p](q,q);
			for (int k=0; k<AMREX_SPACEDIM; k++)
			  rho -= lambda * gradepsilonR[k](k,p);
		      }

		    amrex::Real aa = 0.0;
		    for (int q = 0; q < AMREX_SPACEDIM; q++)
		      {
			aa -= 2.0*mu*gradepsilonD[p](q,q);
			for (int k=0; k<AMREX_SPACEDIM; k++)
			  aa -= lambda * gradepsilonD[k](k,p);
		      }

		    solnfab(amrex::IntVect(i,j),p) = (rhsfab(amrex::IntVect(i,j),p) - rho) / aa;
		  }
       	      }
       	  }
    }
}

/// \fn MLPFStiffnessMatrix::FFlux
/// 
/// Compute the "flux" corresponding to the operator in MLPFStiffnessMatrix::Fapply.
/// Because the operator is self-adjoint and positive-definite, the flux is not
/// required for adequate convergence (?)
/// Therefore, the fluxes are simply set to zero and returned.
/// 
/// \todo Extend to 3D
///
void
MLPFStiffnessMatrix::FFlux (int /*amrlev*/, const MFIter& /*mfi*/,
 			  const std::array<FArrayBox*,AMREX_SPACEDIM>& sigma,
 			  const FArrayBox& /*sol*/, const int /*face_only*/) const
{

  amrex::BaseFab<amrex::Real> &fxfab = *sigma[0];
  amrex::BaseFab<amrex::Real> &fyfab = *sigma[1];

  fxfab.setVal(0.0);
  fyfab.setVal(0.0);
}

int
MLPFStiffnessMatrix::getNComp() const
{
  return AMREX_SPACEDIM;
}



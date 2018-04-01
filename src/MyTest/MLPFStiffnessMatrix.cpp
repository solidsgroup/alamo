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

  mu1 = 1.0; mu2 = 2.0;
  lambda1 = 1.0; lambda2 = 1.0;

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
	m_a_coeffs[amrlev][mglev].define(m_grids[amrlev][mglev],m_dmap[amrlev][mglev],2,1);
    }
}

void
MLPFStiffnessMatrix::setACoeffs(int amrlev,
				const MultiFab& alpha)
{
  MultiFab::Copy(m_a_coeffs[amrlev][0], alpha, 0, 0, alpha.nComp(), 1);
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
		      C(i,0,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi)
		      *(ufab(amrex::IntVect(m+1,n),k) - 2.0*ufab(amrex::IntVect(m,n),k) + ufab(amrex::IntVect(m-1,n),k))/dx[0]/dx[0]
		      +
		      (C(i,0,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi) + C(i,1,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi))
		      *(ufab(amrex::IntVect(m+1,n+1),k) + ufab(amrex::IntVect(m-1,n-1),k) - ufab(amrex::IntVect(m+1,n-1),k) - ufab(amrex::IntVect(m-1,n+1),k))/(2.0*dx[0])/(2.0*dx[1])
		      +
		      C(i,1,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi)
		      *(ufab(amrex::IntVect(m,n+1),k) - 2.0*ufab(amrex::IntVect(m,n),k) + ufab(amrex::IntVect(m,n-1),k))/dx[1]/dx[1];

		    // C_{ijkl,j} u_{k,l}
				amrex::Real temp1, temp2, temp3, temp4;
				if(m == bx.loVect()[0])
				{
					temp1 = (C(i,0,k,0,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi))/(dx[0]);
					temp3 = (C(i,0,k,1,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi))/(dx[0]);
				}
				else if(m == bx.hiVect()[0])
				{
					temp1 = (C(i,0,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi) - C(i,0,k,0,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(dx[0]);
					temp3 = (C(i,0,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi) - C(i,0,k,1,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(dx[0]);
				}
				else
				{
					temp1 = (C(i,0,k,0,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,0,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(2.0*dx[0]);
					temp3 = (C(i,0,k,1,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,1,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(2.0*dx[0]);
				}
				if(n == bx.loVect()[1])
				{
					temp2 = (C(i,1,k,0,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi))/(dx[1]);
					temp4 = (C(i,1,k,1,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi))/(dx[1]);
				}
				else if(n == bx.hiVect()[1])
				{
					temp2 = (C(i,1,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi) - C(i,1,k,0,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(dx[1]);
					temp4 = (C(i,1,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi) - C(i,1,k,1,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(dx[1]);
				}
				else
				{
					temp2 = (C(i,1,k,0,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,0,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(2.0*dx[1]);
					temp4 = (C(i,1,k,1,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,1,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(2.0*dx[1]);
				}

				ffab(amrex::IntVect(m,n),i) -= (temp1 + temp2)*((ufab(amrex::IntVect(m+1,n),k) - ufab(amrex::IntVect(m-1,n),k))/(2.0*dx[0])) +
					(temp3 + temp4)*((ufab(amrex::IntVect(m,n+1),k) - ufab(amrex::IntVect(m,n-1),k))/(2.0*dx[1]));
				//ffab(amrex::IntVect(m,n),i) -=
		      //((C(i,0,k,0,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,0,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(2.0*dx[0]) +
		       //(C(i,1,k,0,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,0,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(2.0*dx[1])) *
		      //((ufab(amrex::IntVect(m+1,n),k) - ufab(amrex::IntVect(m-1,n),k))/(2.0*dx[0]))
		      //+
		      //((C(i,0,k,1,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,1,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(2.0*dx[0]) +
		       //(C(i,1,k,1,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,1,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(2.0*dx[1])) *
		      //((ufab(amrex::IntVect(m,n+1),k) - ufab(amrex::IntVect(m,n-1),k))/(2.0*dx[1]));

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
			      MultiFab& u,       ///<[inout] Solution (displacement field)
			      const MultiFab& rhs, ///<[in] Body force vectors (rhs=right hand side)
			      int redblack         ///<[in] Variable to determine whether to smooth even or odd modes
			      ) const
{
  BL_PROFILE("MLPFStiffnessMatrix::Fsmooth()");

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
			C(i,0,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi)
			*(ufab(amrex::IntVect(m+1,n),k) - (i==k ? 0.0 : 2.0*ufab(amrex::IntVect(m,n),k)) + ufab(amrex::IntVect(m-1,n),k))/dx[0]/dx[0]
			+
			(C(i,0,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi) + C(i,1,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi))
			*(ufab(amrex::IntVect(m+1,n+1),k) + ufab(amrex::IntVect(m-1,n-1),k) - ufab(amrex::IntVect(m+1,n-1),k) - ufab(amrex::IntVect(m-1,n+1),k))/(2.0*dx[0])/(2.0*dx[1])
			+
			C(i,1,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi)
			*(ufab(amrex::IntVect(m,n+1),k) - (i==k ? 0.0 : 2.0*ufab(amrex::IntVect(m,n),k)) + ufab(amrex::IntVect(m,n-1),k))/dx[1]/dx[1];

		      // C_{ijkl,j} u_{k,l}
					amrex::Real temp1, temp2, temp3, temp4;
					if(m == tbx.loVect()[0])
					{
						temp1 = (C(i,0,k,0,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi))/(dx[0]);
						temp3 = (C(i,0,k,1,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi))/(dx[0]);
					}
					else if(m == tbx.hiVect()[0])
					{
						temp1 = (C(i,0,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi) - C(i,0,k,0,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(dx[0]);
						temp3 = (C(i,0,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi) - C(i,0,k,1,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(dx[0]);
					}
					else
					{
						temp1 = (C(i,0,k,0,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,0,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(2.0*dx[0]);
						temp3 = (C(i,0,k,1,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,1,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(2.0*dx[0]);
					}
					if(n == tbx.loVect()[1])
					{
						temp2 = (C(i,1,k,0,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi))/(dx[1]);
						temp4 = (C(i,1,k,1,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi))/(dx[1]);
					}
					else if(n == tbx.hiVect()[1])
					{
						temp2 = (C(i,1,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi) - C(i,1,k,0,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(dx[1]);
						temp4 = (C(i,1,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi) - C(i,1,k,1,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(dx[1]);
					}
					else
					{
						temp2 = (C(i,1,k,0,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,0,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(2.0*dx[1]);
						temp4 = (C(i,1,k,1,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,1,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(2.0*dx[1]);
					}

					rho -= (temp1+temp2)*((ufab(amrex::IntVect(m+1,n),k) - ufab(amrex::IntVect(m-1,n),k))/(2.0*dx[0])) +
						(temp3+temp4)*((ufab(amrex::IntVect(m,n+1),k) - ufab(amrex::IntVect(m,n-1),k))/(2.0*dx[1]));
		      //rho -=
		        //((C(i,0,k,0,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,0,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(2.0*dx[0]) +
		         //(C(i,1,k,0,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,0,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(2.0*dx[1])) *
		        //((ufab(amrex::IntVect(m+1,n),k) - ufab(amrex::IntVect(m-1,n),k))/(2.0*dx[0]))
		        //+
		        //((C(i,0,k,1,amrex::IntVect(m+1,n),amrlev,mglev,mfi) - C(i,0,k,1,amrex::IntVect(m-1,n),amrlev,mglev,mfi))/(2.0*dx[0]) +
		         //(C(i,1,k,1,amrex::IntVect(m,n+1),amrlev,mglev,mfi) - C(i,1,k,1,amrex::IntVect(m,n-1),amrlev,mglev,mfi))/(2.0*dx[1])) *
		        //((ufab(amrex::IntVect(m,n+1),k) - ufab(amrex::IntVect(m,n-1),k))/(2.0*dx[1]));



		      if (rho != rho){
		      std::cout << __LINE__ << ": " << C(i,0,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi) << std::endl;
		      std::cout << __LINE__ << ": " << ufab(amrex::IntVect(m+1,n),k)<< std::endl;
		      std::cout << __LINE__ << ": " << (i==k ? 0.0 : 2.0*ufab(amrex::IntVect(m,n),k))<< std::endl;
		      std::cout << __LINE__ << ": " << ufab(amrex::IntVect(m-1,n),k)<< std::endl;
		      std::cout << __LINE__ << ": " << C(i,0,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi) << std::endl;
		      std::cout << __LINE__ << ": " << C(i,1,k,0,amrex::IntVect(m,n),amrlev,mglev,mfi)<< std::endl;
		      std::cout << __LINE__ << ": " << ufab(amrex::IntVect(m+1,n+1),k)<< std::endl;
		      std::cout << __LINE__ << ": " << ufab(amrex::IntVect(m-1,n-1),k)<< std::endl;
		      std::cout << __LINE__ << ": " << ufab(amrex::IntVect(m+1,n-1),k)<< std::endl;
		      std::cout << __LINE__ << ": " << ufab(amrex::IntVect(m-1,n+1),k)<< std::endl;
		      std::cout << __LINE__ << ": " << C(i,1,k,1,amrex::IntVect(m,n),amrlev,mglev,mfi)<< std::endl;
		      std::cout << __LINE__ << ": " << ufab(amrex::IntVect(m,n+1),k)<< std::endl;
		      std::cout << __LINE__ << ": " << (i==k ? 0.0 : 2.0*ufab(amrex::IntVect(m,n),k))<< std::endl;
		      std::cout << __LINE__ << ": " << ufab(amrex::IntVect(m,n-1),k)<< std::endl;
		      std::cout << __LINE__ << ": " << C(i,0,k,0,amrex::IntVect(m+1,n),amrlev,mglev,mfi)<< std::endl;
		      std::cout << __LINE__ << ": " << C(i,0,k,0,amrex::IntVect(m-1,n),amrlev,mglev,mfi)<< std::endl; // NAN NAN
		      std::cout << __LINE__ << ": " << C(i,1,k,0,amrex::IntVect(m,n+1),amrlev,mglev,mfi)<< std::endl;
		      std::cout << __LINE__ << ": " << C(i,1,k,0,amrex::IntVect(m,n-1),amrlev,mglev,mfi)<< std::endl;
		      std::cout << __LINE__ << ": " << ufab(amrex::IntVect(m+1,n),k)<< std::endl;
		      std::cout << __LINE__ << ": " << ufab(amrex::IntVect(m-1,n),k)<< std::endl;
		      std::cout << __LINE__ << ": " << C(i,0,k,1,amrex::IntVect(m+1,n),amrlev,mglev,mfi)<< std::endl;
		      std::cout << __LINE__ << ": " << C(i,0,k,1,amrex::IntVect(m-1,n),amrlev,mglev,mfi)<< std::endl;
		      std::cout << __LINE__ << ": " << C(i,1,k,1,amrex::IntVect(m,n+1),amrlev,mglev,mfi)<< std::endl;
		      std::cout << __LINE__ << ": " << C(i,1,k,1,amrex::IntVect(m,n-1),amrlev,mglev,mfi)<< std::endl; // NAN NAN
		      std::cout << __LINE__ << ": " << ufab(amrex::IntVect(m,n+1),k)<< std::endl;
		      std::cout << __LINE__ << ": " << ufab(amrex::IntVect(m,n-1),k)<< std::endl;
		      std::cout << std::endl;
		      amrex::Abort("nans detected");}




		    }

		  aa -=
		    -2.0*C(i,0,i,0,amrex::IntVect(m,n),amrlev,mglev,mfi)/dx[0]/dx[0]
		    -2.0*C(i,1,i,1,amrex::IntVect(m,n),amrlev,mglev,mfi)/dx[1]/dx[1];

		  //if (rho != rho) amrex::Abort("nans detected");
		  //if (rho != rho) std::cout << "nans detetected" << std::endl;

		  ufab(amrex::IntVect(m,n),i) = (rhsfab(amrex::IntVect(m,n),i) - rho) / aa;
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


amrex::Real
MLPFStiffnessMatrix::C(const int i, const int j, const int k, const int l,
		       const amrex::IntVect loc,
		       int amrlev, int mglev, MFIter &mfi) const
{
  amrex::Real mu, lambda;

  // amrex::Real x = m_geom[amrlev][mglev].ProbLo()[0] + ((amrex::Real)(loc[0]) + 0.5) * m_geom[amrlev][mglev].CellSize()[0];
  //  amrex::Real y = m_geom[amrlev][mglev].ProbLo()[1] + ((amrex::Real)(loc[1]) + 0.5) * m_geom[amrlev][mglev].CellSize()[1];
  //  if (y>0.0) {mu=mu1; lambda=lambda1; }
  //  else {mu=mu2; lambda=lambda2;}

  const amrex::BaseFab<amrex::Real> &etafab = m_a_coeffs[amrlev][mglev][mfi];
  mu = (mu1*etafab(loc,0) + mu2*etafab(loc,1)) / (etafab(loc,0) + etafab(loc,1));
  lambda = (lambda1*etafab(loc,0) + lambda2*etafab(loc,1)) / (etafab(loc,0) + etafab(loc,1));

  // TODO: This is a hack, needs to be fixed.
  if (mu != mu)
	{
		std::cout << "Nans detected i = " << i << " j = " << j << " k = " << k << " l = " << l << std:: endl;
		std::cout << "etafab = " << etafab(loc,0) << " , " << etafab(loc,1) << std::endl;
		mu = 0.5*(mu1+mu2);
	}
  if (lambda != lambda) lambda = 0.5*(lambda1+lambda2);

  //std::cout << "(" << mu << ","<<lambda << ") ";

  amrex::Real ret = 0.0;
  if (i==k && j==l) ret += mu;
  if (i==l && j==k) ret += mu;
  if (i==j && k==l) ret += lambda;

  if (ret != ret)
    {
      //return 0.
	//std::cout << "eta1 = " << etafab(loc,0) << " eta2 = " << etafab(loc,1) << std::endl;
      //std::cout << "mu="<<mu<<" lambda="<<lambda<<std::endl;
      amrex::Abort("NAN DETECTED");
    }

  return ret;
};

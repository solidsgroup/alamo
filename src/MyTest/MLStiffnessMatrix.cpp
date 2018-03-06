#include <AMReX_MultiFabUtil.H>
#include "MLStiffnessMatrix.H"
#include "eigen3/Eigen/Core"
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>

#include <AMReX_ArrayLim.H>

namespace amrex {

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

  // mu = 1.5;
  // lambda = 1.0;

}

void
MLStiffnessMatrix::define (const Vector<Geometry>& a_geom,
			   const Vector<BoxArray>& a_grids,
			   const Vector<DistributionMapping>& a_dmap,
			   const LPInfo& a_info)
{
  MLCellLinOp::define(a_geom, a_grids, a_dmap, a_info);

  m_a_coeffs.resize(m_num_amr_levels);
  m_b_coeffs.resize(m_num_amr_levels);
  for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
      m_a_coeffs[amrlev].resize(m_num_mg_levels[amrlev]);
      m_b_coeffs[amrlev].resize(m_num_mg_levels[amrlev]);
      for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
	  m_a_coeffs[amrlev][mglev].define(m_grids[amrlev][mglev],
					   m_dmap[amrlev][mglev],
					   1, 0);
	  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
	      const BoxArray& ba = amrex::convert(m_grids[amrlev][mglev], IntVect::TheDimensionVector(idim));
	      m_b_coeffs[amrlev][mglev][idim].define(ba,
						     m_dmap[amrlev][mglev],
						     1, 0);
            }
        }
    }
}

MLStiffnessMatrix::~MLStiffnessMatrix ()
{}

void
MLStiffnessMatrix::setScalars (Real a, Real b)
{
  m_a_scalar = a;
  m_b_scalar = b;
  if (a == 0.0)
    {
      for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
        {
	  m_a_coeffs[amrlev][0].setVal(0.0);
        }
    }
}

void
MLStiffnessMatrix::setACoeffs (int amrlev, const MultiFab& alpha)
{
  MultiFab::Copy(m_a_coeffs[amrlev][0], alpha, 0, 0, 1, 0);
}

void
MLStiffnessMatrix::setBCoeffs (int amrlev,
			       const std::array<MultiFab const*,AMREX_SPACEDIM>& beta)
{
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    MultiFab::Copy(m_b_coeffs[amrlev][0][idim], *beta[idim], 0, 0, 1, 0);
  }
}

void
MLStiffnessMatrix::averageDownCoeffs ()
{

  for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
    {
      auto& fine_a_coeffs = m_a_coeffs[amrlev];
      auto& fine_b_coeffs = m_b_coeffs[amrlev];

      averageDownCoeffsSameAmrLevel(fine_a_coeffs, fine_b_coeffs);
      averageDownCoeffsToCoarseAmrLevel(amrlev);
    }

  averageDownCoeffsSameAmrLevel(m_a_coeffs[0], m_b_coeffs[0]);
}

void
MLStiffnessMatrix::averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a,
						  Vector<std::array<MultiFab,AMREX_SPACEDIM> >& b)
{
  int nmglevs = a.size();
  for (int mglev = 1; mglev < nmglevs; ++mglev)
    {
      if (m_a_scalar == 0.0)
        {
	  a[mglev].setVal(0.0);
        }
      else
        {
	  amrex::average_down(a[mglev-1], a[mglev], 0, 1, mg_coarsen_ratio);
        }
        
      Vector<const MultiFab*> fine {AMREX_D_DECL(&(b[mglev-1][0]),
						 &(b[mglev-1][1]),
						 &(b[mglev-1][2]))};
      Vector<MultiFab*> crse {AMREX_D_DECL(&(b[mglev][0]),
					   &(b[mglev][1]),
					   &(b[mglev][2]))};
      IntVect ratio {mg_coarsen_ratio};
      amrex::average_down_faces(fine, crse, ratio, 0);
    }
}

void
MLStiffnessMatrix::averageDownCoeffsToCoarseAmrLevel (int flev)
{
  auto& fine_a_coeffs = m_a_coeffs[flev  ].back();
  auto& fine_b_coeffs = m_b_coeffs[flev  ].back();
  auto& crse_a_coeffs = m_a_coeffs[flev-1].front();
  auto& crse_b_coeffs = m_b_coeffs[flev-1].front();
  auto& crse_geom     = m_geom    [flev-1][0];

  if (m_a_scalar != 0.0) {
    amrex::average_down(fine_a_coeffs, crse_a_coeffs, 0, 1, mg_coarsen_ratio);
  }
     
  std::array<MultiFab,AMREX_SPACEDIM> bb;
  Vector<MultiFab*> crse(AMREX_SPACEDIM);
  Vector<MultiFab const*> fine(AMREX_SPACEDIM);
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    BoxArray ba = fine_b_coeffs[idim].boxArray();
    ba.coarsen(mg_coarsen_ratio);
    bb[idim].define(ba, fine_b_coeffs[idim].DistributionMap(), 1, 0);
    crse[idim] = &bb[idim];
    fine[idim] = &fine_b_coeffs[idim];
  }
  IntVect ratio {mg_coarsen_ratio};
  amrex::average_down_faces(fine, crse, ratio, 0);

  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    crse_b_coeffs[idim].ParallelCopy(bb[idim], crse_geom.periodicity());
  }
}

void
MLStiffnessMatrix::applyMetricTermsCoeffs ()
{
#if (AMREX_SPACEDIM != 3)
    for (int alev = 0; alev < m_num_amr_levels; ++alev)
    {
        const int mglev = 0;
        applyMetricTerm(alev, mglev, m_a_coeffs[alev][mglev]);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            applyMetricTerm(alev, mglev, m_b_coeffs[alev][mglev][idim]);
        }
    }
#endif
}

void
MLStiffnessMatrix::prepareForSolve ()
{
  BL_PROFILE("MLStiffnessMatrix::prepareForSolve()");

  MLCellLinOp::prepareForSolve();

  applyMetricTermsCoeffs();

  averageDownCoeffs();

  m_is_singular.clear();
  m_is_singular.resize(m_num_amr_levels, false);
  auto itlo = std::find(m_lobc.begin(), m_lobc.end(), BCType::Dirichlet);
  auto ithi = std::find(m_hibc.begin(), m_hibc.end(), BCType::Dirichlet);
  if (itlo == m_lobc.end() && ithi == m_hibc.end())
    {  // No Dirichlet
      for (int alev = 0; alev < m_num_amr_levels; ++alev)
        {
	  if (m_domain_covered[alev])
            {
	      if (m_a_scalar == 0.0)
                {
		  m_is_singular[alev] = true;
                }
	      else
                {
		  Real asum = m_a_coeffs[alev].back().sum();
		  Real amax = m_a_coeffs[alev].back().norm0();
		  m_is_singular[alev] = (asum <= amax * 1.e-12);
                }
            }
        }
    }
}

void
MLStiffnessMatrix::Fapply (int amrlev, int mglev, MultiFab& f, const MultiFab& u) const
{
  BL_PROFILE("MLStiffnessMatrix::Fapply()");

  // const MultiFab& acoef = m_a_coeffs[amrlev][mglev];
  // const MultiFab& bxcoef = m_b_coeffs[amrlev][mglev][0];
  // const MultiFab& bycoef = m_b_coeffs[amrlev][mglev][1];

  // const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();
  const Real* dx = m_geom[amrlev][mglev].CellSize();

  // if (f.nComp() != AMREX_SPACEDIM)
  //   amrex::Abort("ncomp not correct");

  for (MFIter mfi(f, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      //amrex::Real alpha=m_a_scalar, beta=m_b_scalar;
      //amrex::Real dhx = beta*dxinv[0]*dxinv[0], dhy = beta*dxinv[1]*dxinv[1];

      const amrex::BaseFab<amrex::Real> &ufab  = u[mfi];
      amrex::BaseFab<amrex::Real>       &ffab  = f[mfi];
      // const amrex::BaseFab<amrex::Real> &afab  = acoef[mfi];
      // const amrex::BaseFab<amrex::Real> &bxfab = bxcoef[mfi];
      // const amrex::BaseFab<amrex::Real> &byfab = bycoef[mfi];

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

void
MLStiffnessMatrix::Fsmooth (int amrlev,
			    int mglev,
			    MultiFab& sol,
			    const MultiFab& rhs,
			    int redblack) const
{
  BL_PROFILE("MLStiffnessMatrix::Fsmooth()");

  const MultiFab& acoef = m_a_coeffs[amrlev][mglev];
  AMREX_D_TERM(const MultiFab& bxcoef = m_b_coeffs[amrlev][mglev][0];,
	       const MultiFab& bycoef = m_b_coeffs[amrlev][mglev][1];,
	       const MultiFab& bzcoef = m_b_coeffs[amrlev][mglev][2];);

  const auto& undrrelxr = m_undrrelxr[amrlev][mglev];
  const auto& maskvals  = m_maskvals [amrlev][mglev];

  OrientationIter oitr;

  const FabSet& f0 = undrrelxr[oitr()]; ++oitr;
  const FabSet& f1 = undrrelxr[oitr()]; ++oitr;
  const FabSet& f2 = undrrelxr[oitr()]; ++oitr;
  const FabSet& f3 = undrrelxr[oitr()]; ++oitr;

  const MultiMask& mm0 = maskvals[0];
  const MultiMask& mm1 = maskvals[1];
  const MultiMask& mm2 = maskvals[2];
  const MultiMask& mm3 = maskvals[3];

  const int nComp = AMREX_SPACEDIM;
  const Real* h = m_geom[amrlev][mglev].CellSize();
  const Real* dx = m_geom[amrlev][mglev].CellSize();

  for (MFIter mfi(sol,MFItInfo().EnableTiling().SetDynamic(true));
       mfi.isValid(); ++mfi)
    {
      const Mask& m0 = mm0[mfi];
      const Mask& m1 = mm1[mfi];
      const Mask& m2 = mm2[mfi];
      const Mask& m3 = mm3[mfi];

      const Box&       tbx     = mfi.tilebox();
      const Box&       vbx     = mfi.validbox();
      FArrayBox&       solnfab = sol[mfi];
      const FArrayBox& rhsfab  = rhs[mfi];
      const FArrayBox& afab    = acoef[mfi];

      AMREX_D_TERM(const FArrayBox& bxfab = bxcoef[mfi];,
		   const FArrayBox& byfab = bycoef[mfi];,
		   const FArrayBox& bzfab = bzcoef[mfi];);

      
      amrex::Real alpha = m_a_scalar, beta=m_b_scalar;
      amrex::Real dhx = beta/h[0]/h[0];
      amrex::Real dhy = beta/h[1]/h[1];

      const FArrayBox& f0fab = f0[mfi];
      const FArrayBox& f1fab = f1[mfi];
      const FArrayBox& f2fab = f2[mfi];
      const FArrayBox& f3fab = f3[mfi];

      // phi, DIMS(phi)    ---> solnfab
      // rhs, DIMS(rhs)    ---> rhsfab
      // f0, DIMS(f0)      ---> f0fab
      // m0, DIMS(m0)      ---> m0
      // f1, DIMS(f1)      ---> f1fab
      // m1, DIMS(m1)      ---> m1
      // f2, DIMS(f2)      ---> f2fab
      // m2, DIMS(m2)      ---> f2
      // f3, DIMS(f3)      ---> f3fab
      // m3, DIMS(m3)      ---> f3
      // lo, hi,           ---> tbx.loVect, tbx.hiVect
      // blo, bhi          ---> vbx.loVect, vbx.hiVect
      // nc, h, redblack   ---> nc, h, redblack

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

		// for (int p = 0; p < AMREX_SPACEDIM; p++)
		//   {
		//     ffab(amrex::IntVect(i,j),p) = 0.0; 
		
		for (int q = 0; q < AMREX_SPACEDIM; q++)
		  {
		    rho -= 2.0*mu*gradepsilonR[n](q,q);
		    if (n==q)
		      for (int k=0; k<AMREX_SPACEDIM; k++)
			rho -= lambda * gradepsilonR[k](k,q);
		  }
		  // }	    

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


/// \todo Currently planning to abandon FEM...but may want to hang onto this just in case
amrex::Real
MLStiffnessMatrix::K (amrex::IntVect m, amrex::IntVect n,int i, int k)
{

  //  ___________________________
  // |             |             |
  // |		   |	         |
  // |    Phi2	   |	Phi1     |
  // |		   |	         |
  // |		   |	         |
  // |_____________|_____________|
  // |		   |	         |
  // |		   |	         |
  // |	  Phi3	   |	Phi4     |
  // |		   |	         |
  // |		   |	         |
  // |_____________|_____________|
  //
  //

  amrex::Real DPhi1DPhi1[2][2] =  {{  1./3.,   1./4.  },{  1./4.,  1./3.  }};
  amrex::Real DPhi1DPhi2[2][2] =  {{ -1./3.,   1./4.  },{ -1./4.,  1./6.  }};
  amrex::Real DPhi1DPhi3[2][2] =  {{ -1./6.,  -1./4.  },{ -1./4., -1./6.  }};
  amrex::Real DPhi1DPhi4[2][2] =  {{  1./6.,  -1./4.  },{  1./4., -1./3.  }};
  
  amrex::Real DPhi2DPhi1[2][2] =  {{ -1./3.,  -1./4.  },{  1./4.,  1./6.  }};
  amrex::Real DPhi2DPhi2[2][2] =  {{  1./3.,  -1./4.  },{ -1./4.,  1./3.  }};
  amrex::Real DPhi2DPhi3[2][2] =  {{  1./6.,   1./4.  },{ -1./4., -1./3.  }};
  amrex::Real DPhi2DPhi4[2][2] =  {{ -1./6.,   1./4.  },{  1./4., -1./6.  }};

  amrex::Real DPhi3DPhi1[2][2] =  {{ -1./6.,  -1./4.  },{ -1./4., -1./6.  }};
  amrex::Real DPhi3DPhi2[2][2] =  {{  1./6.,  -1./4.  },{  1./4., -1./3.  }};
  amrex::Real DPhi3DPhi3[2][2] =  {{  1./3.,   1./4.  },{  1./4.,  1./3.  }};
  amrex::Real DPhi3DPhi4[2][2] =  {{ -1./3.,   1./4.  },{ -1./4.,  1./6.  }};

  amrex::Real DPhi4DPhi1[2][2] =  {{  1./6.,   1./4.  },{ -1./4., -1./3.  }};
  amrex::Real DPhi4DPhi2[2][2] =  {{ -1./6.,   1./4.  },{  1./4., -1./6.  }};
  amrex::Real DPhi4DPhi3[2][2] =  {{ -1./3.,  -1./4.  },{  1./4.,  1./6.  }};
  amrex::Real DPhi4DPhi4[2][2] =  {{  1./3.,  -1./4.  },{ -1./4.,  1./3.  }};


  amrex::Real DPhiDPhi[2][2] = {{0,0},{0,0}};


  if (m - n == amrex::IntVect(1,1))
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) DPhiDPhi[i][j] = DPhi1DPhi3[i][j];

  if (m - n == amrex::IntVect(-1,1))
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) DPhiDPhi[i][j] = DPhi2DPhi4[i][j];

  if (m - n == amrex::IntVect(-1,-1))
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) DPhiDPhi[i][j] = DPhi3DPhi1[i][j];

  if (m - n == amrex::IntVect(1,-1))
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) DPhiDPhi[i][j] = DPhi4DPhi2[i][j];

  if (m - n == amrex::IntVect(1,0))
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) DPhiDPhi[i][j] = DPhi1DPhi2[i][j] + DPhi4DPhi3[i][j];

  if (m - n == amrex::IntVect(-1,0))
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) DPhiDPhi[i][j] = DPhi2DPhi1[i][j] + DPhi3DPhi4[i][j];

  if (m - n == amrex::IntVect(0,1))
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) DPhiDPhi[i][j] = DPhi1DPhi4[i][j] + DPhi2DPhi3[i][j];

  if (m - n == amrex::IntVect(0,-1))
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) DPhiDPhi[i][j] = DPhi4DPhi1[i][j] + DPhi3DPhi2[i][j];

  if (m == n)
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) DPhiDPhi[i][j] = DPhi1DPhi1[i][j] + DPhi2DPhi2[i][j] + DPhi3DPhi3[i][j] + DPhi4DPhi4[i][j];

  return mu * (i==k ? 1. : 0.) * (DPhiDPhi[0][0] + DPhiDPhi[1][1]) + (mu + lambda) * DPhiDPhi[i][k];

}

}


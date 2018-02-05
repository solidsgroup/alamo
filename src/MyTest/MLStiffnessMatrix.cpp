#if AMREX_SPACEDIM==2
hello
#include <AMReX_MultiFabUtil.H>
#include "MLStiffnessMatrix.H"
#include "AMReX_MLABecLap_F.H"
#include "AMReX_ABec_F.H"

namespace amrex {

MLStiffnessMatrix::MLStiffnessMatrix (const Vector<Geometry>& a_geom,
				      const Vector<BoxArray>& a_grids,
				      const Vector<DistributionMapping>& a_dmap,
				      const LPInfo& a_info)
{
  define(a_geom, a_grids, a_dmap, a_info);
}

void
MLStiffnessMatrix::define (const Vector<Geometry>& a_geom,
			   const Vector<BoxArray>& a_grids,
			   const Vector<DistributionMapping>& a_dmap,
			   const LPInfo& a_info)
{
  BL_PROFILE("MLStiffnessMatrix::define()");

  MLLinOp::define(a_geom, a_grids, a_dmap, a_info);

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
  BL_PROFILE("MLStiffnessMatrix::averageDownCoeffs()");

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
  for (int alev = 0; alev < m_num_amr_levels; ++alev)
    {
      for (int mglev = 0; mglev < m_num_mg_levels[alev]; ++mglev)
        {
	  applyMetricTerm(alev, mglev, m_a_coeffs[alev][mglev]);
	  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
	      applyMetricTerm(alev, mglev, m_b_coeffs[alev][mglev][idim]);
            }
        }
    }
}

void
MLStiffnessMatrix::prepareForSolve ()
{
  BL_PROFILE("MLStiffnessMatrix::prepareForSolve()");

  MLLinOp::prepareForSolve();

  applyMetricTermsCoeffs();

  averageDownCoeffs();

  m_Anorm.resize(m_num_amr_levels);
  for (int alev = 0; alev < m_num_amr_levels; ++alev)
    {
      m_Anorm[alev].assign(m_num_mg_levels[alev], -1.0);
    }

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
MLStiffnessMatrix::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
  BL_PROFILE("MLStiffnessMatrix::Fapply()");

  const MultiFab& acoef = m_a_coeffs[amrlev][mglev];
  const MultiFab& bxcoef = m_b_coeffs[amrlev][mglev][0];
  const MultiFab& bycoef = m_b_coeffs[amrlev][mglev][1];

  const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

  for (MFIter mfi(out, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      amrex::Real alpha=m_a_scalar, beta=m_b_scalar;
      amrex::Real dhx = beta*dxinv[0]*dxinv[0], dhy = beta*dxinv[1]*dxinv[1];

      const amrex::BaseFab<amrex::Real> &xfab  = in[mfi];
      amrex::BaseFab<amrex::Real>       &yfab  = out[mfi];
      const amrex::BaseFab<amrex::Real> &afab  = acoef[mfi];
      const amrex::BaseFab<amrex::Real> &bxfab = bxcoef[mfi];
      const amrex::BaseFab<amrex::Real> &byfab = bycoef[mfi];

      for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
	for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
	  {
	    yfab(amrex::IntVect(i,j)) = alpha*afab(amrex::IntVect(i,j))*xfab(amrex::IntVect(i,j))
	      - dhx * (bxfab(amrex::IntVect(i+1,j))*(xfab(amrex::IntVect(i+1,j)) - xfab(amrex::IntVect(i,j)))
		       - bxfab(amrex::IntVect(i,j))*(xfab(amrex::IntVect(i,j)) - xfab(amrex::IntVect(i-1,j))))
	      - dhy * (byfab(amrex::IntVect(i,j+1))*(xfab(amrex::IntVect(i,j+1)) - xfab(amrex::IntVect(i,j)))
		       - byfab(amrex::IntVect(i,j))*(xfab(amrex::IntVect(i,j)) - xfab(amrex::IntVect(i,j-1))));
	  }
    }
}

void
MLStiffnessMatrix::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const
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

  const int nc = 1;
  const Real* h = m_geom[amrlev][mglev].CellSize();

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

      const FArrayBox& f0fab = f0[mfi];
      const FArrayBox& f1fab = f1[mfi];
      const FArrayBox& f2fab = f2[mfi];
      const FArrayBox& f3fab = f3[mfi];

      FORT_GSRB(solnfab.dataPtr(), ARLIM(solnfab.loVect()),ARLIM(solnfab.hiVect()),
		rhsfab.dataPtr(), ARLIM(rhsfab.loVect()), ARLIM(rhsfab.hiVect()),
		&m_a_scalar, &m_b_scalar,
		afab.dataPtr(), ARLIM(afab.loVect()),    ARLIM(afab.hiVect()),
		bxfab.dataPtr(), ARLIM(bxfab.loVect()),   ARLIM(bxfab.hiVect()),
		byfab.dataPtr(), ARLIM(byfab.loVect()),   ARLIM(byfab.hiVect()),
		f0fab.dataPtr(), ARLIM(f0fab.loVect()),   ARLIM(f0fab.hiVect()),
		m0.dataPtr(), ARLIM(m0.loVect()),   ARLIM(m0.hiVect()),
		f1fab.dataPtr(), ARLIM(f1fab.loVect()),   ARLIM(f1fab.hiVect()),
		m1.dataPtr(), ARLIM(m1.loVect()),   ARLIM(m1.hiVect()),
		f2fab.dataPtr(), ARLIM(f2fab.loVect()),   ARLIM(f2fab.hiVect()),
		m2.dataPtr(), ARLIM(m2.loVect()),   ARLIM(m2.hiVect()),
		f3fab.dataPtr(), ARLIM(f3fab.loVect()),   ARLIM(f3fab.hiVect()),
		m3.dataPtr(), ARLIM(m3.loVect()),   ARLIM(m3.hiVect()),
		tbx.loVect(), tbx.hiVect(), vbx.loVect(), vbx.hiVect(),
		&nc, h, &redblack);

    }
}

void
MLStiffnessMatrix::FFlux (int amrlev, const MFIter& mfi,
 			  const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
 			  const FArrayBox& sol, const int face_only) const
{
  BL_PROFILE("MLStiffnessMatrix::FFlux()");

  const int mglev = 0;
  const auto& bx = m_b_coeffs[amrlev][mglev][0][mfi];
  const auto& by = m_b_coeffs[amrlev][mglev][1][mfi];

  const Box& box = mfi.tilebox();
  const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

  amrex::Real beta=m_b_scalar;
  
  amrex::Real dhx = beta*dxinv[0], dhy = beta*dxinv[1];

  amrex::BaseFab<amrex::Real> &fxfab = *flux[0];
  amrex::BaseFab<amrex::Real> &fyfab = *flux[1];
  const amrex::BaseFab<amrex::Real> &solfab = sol;

  if (face_only)
    {
      for (int i = box.loVect()[0]; i<=box.hiVect()[0]+1; i+= 1 + (box.hiVect()[0]-box.loVect()[0]))
	for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++)
	  {
	    fxfab(amrex::IntVect(i,j)) = -dhx *
	      bx(amrex::IntVect(i,j)) *
	      (solfab(amrex::IntVect(i,j))-solfab(amrex::IntVect(i,j)));
	  }
      for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++)
	for (int j = box.loVect()[1]; j<=box.hiVect()[1]+1; j+= 1 + (box.hiVect()[1]-box.loVect()[1]))
	  {
	    fyfab(amrex::IntVect(i,j)) = -dhy *
	      by(amrex::IntVect(i,j)) *
	      (solfab(amrex::IntVect(i,j))-solfab(amrex::IntVect(i,j-1)));
	  }
    }
  else
    {
      for (int i = box.loVect()[0]; i<=box.hiVect()[0]+1; i++)
	for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++)
	  {
	    fxfab(amrex::IntVect(i,j)) = -dhx *
	      bx(amrex::IntVect(i,j)) *
	      (solfab(amrex::IntVect(i,j))-solfab(amrex::IntVect(i,j)));
	  }
       
      for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++)
	for (int j = box.loVect()[1]; j<=box.hiVect()[1]+1; j++)
	  {
	    fyfab(amrex::IntVect(i,j)) = -dhy *
	      by(amrex::IntVect(i,j)) *
	      (solfab(amrex::IntVect(i,j))-solfab(amrex::IntVect(i,j-1)));
	  }
    }
}

Real
MLStiffnessMatrix::Anorm (int amrlev, int mglev) const
{
  BL_PROFILE("MLStiffnessMatrix::Anorm()");

  if (m_Anorm[amrlev][mglev] < 0.0)
    {
      const MultiFab& acoef = m_a_coeffs[amrlev][mglev];
      AMREX_D_TERM(const MultiFab& bxcoef = m_b_coeffs[amrlev][mglev][0];,
		   const MultiFab& bycoef = m_b_coeffs[amrlev][mglev][1];,
		   const MultiFab& bzcoef = m_b_coeffs[amrlev][mglev][2];);
      const Real* dx = m_geom[amrlev][mglev].CellSize();

      const int nc = 1;
      Real res = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(max:res)
#endif
      {
	for (MFIter mfi(acoef,true); mfi.isValid(); ++mfi)
	  {
	    Real tres;
	    
	    const Box&       tbx  = mfi.tilebox();
	    const FArrayBox& afab = acoef[mfi];
	    AMREX_D_TERM(const FArrayBox& bxfab = bxcoef[mfi];,
			 const FArrayBox& byfab = bycoef[mfi];,
			 const FArrayBox& bzfab = bzcoef[mfi];);

	    FORT_NORMA(&tres,
		       &m_a_scalar, &m_b_scalar,
		       afab.dataPtr(),  ARLIM(afab.loVect()), ARLIM(afab.hiVect()),
		       bxfab.dataPtr(), ARLIM(bxfab.loVect()), ARLIM(bxfab.hiVect()),
		       byfab.dataPtr(), ARLIM(byfab.loVect()), ARLIM(byfab.hiVect()),
		       tbx.loVect(), tbx.hiVect(), &nc, dx);
                
	    res = std::max(res, tres);
	  }
      }
        
      ParallelAllReduce::Max(res, Communicator(amrlev,mglev));
      m_Anorm[amrlev][mglev] = res;
    }

  return m_Anorm[amrlev][mglev];
}


#if (AMREX_SPACEDIM == 2)
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
#endif


}

#endif

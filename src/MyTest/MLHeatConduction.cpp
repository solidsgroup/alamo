#include <AMReX_MultiFabUtil.H>
#include "MLHeatConduction.H"
#include "eigen3/Eigen/Core"
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>

#include <AMReX_ArrayLim.H>

namespace amrex {

MLHeatConduction::MLHeatConduction (const Vector<Geometry>& a_geom,
				      const Vector<BoxArray>& a_grids,
				      const Vector<DistributionMapping>& a_dmap,
				      const LPInfo& a_info)
{
  define(a_geom, a_grids, a_dmap, a_info);

  mu = 0.1;
  lambda = 1.;

}

void
MLHeatConduction::define (const Vector<Geometry>& a_geom,
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

MLHeatConduction::~MLHeatConduction ()
{}

void
MLHeatConduction::setScalars (Real a, Real b)
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
MLHeatConduction::setACoeffs (int amrlev, const MultiFab& alpha)
{
  MultiFab::Copy(m_a_coeffs[amrlev][0], alpha, 0, 0, 1, 0);
}

void
MLHeatConduction::setBCoeffs (int amrlev,
			       const std::array<MultiFab const*,AMREX_SPACEDIM>& beta)
{
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    MultiFab::Copy(m_b_coeffs[amrlev][0][idim], *beta[idim], 0, 0, 1, 0);
  }
}

void
MLHeatConduction::averageDownCoeffs ()
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
MLHeatConduction::averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a,
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
MLHeatConduction::averageDownCoeffsToCoarseAmrLevel (int flev)
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
MLHeatConduction::applyMetricTermsCoeffs ()
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
MLHeatConduction::prepareForSolve ()
{
  BL_PROFILE("MLHeatConduction::prepareForSolve()");

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
MLHeatConduction::Fapply (int amrlev, int mglev, MultiFab& f, const MultiFab& u) const
{
  BL_PROFILE("MLHeatConduction::Fapply()");

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
	    /// \todo This is temporary! Just checking to see if I can get it to work with a Laplacian...

	    // for (int n = 0; n<u.nComp() ; n++)
	    //   ffab(amrex::IntVect(i,j),n) = 
	    //  	- dhx * ((ufab(amrex::IntVect(i+1,j),n) - ufab(amrex::IntVect(i,j),n))
	    // 		 - (ufab(amrex::IntVect(i,j),n) - ufab(amrex::IntVect(i-1,j),n)))
	    //  	- dhy * ((ufab(amrex::IntVect(i,j+1),n) - ufab(amrex::IntVect(i,j),n))
	    //  		 - (ufab(amrex::IntVect(i,j),n) - ufab(amrex::IntVect(i,j-1),n)));
	    // continue;

	    /// \todo This does **NOT** generalize to materials with non-constant \f$\mathbb{C}\f$! A more general finite difference scheme is needed.
	    /// \note `gradepsilon[n](i,j) gives the ij component of epsilon differentiated with respect to \f$x_n\f$
	    /// \note Using `AMREX_SPACEDIM` here--that does NOT mean this generalizes to 3D (yet)

	    amrex::Array<Eigen::Matrix<amrex::Real,AMREX_SPACEDIM,AMREX_SPACEDIM> > gradepsilon(AMREX_SPACEDIM);
	    
	    for (int n = 0; n < AMREX_SPACEDIM; n++)
	    gradepsilon[n] <<
	      (ufab(amrex::IntVect(i+1,j),n) + ufab(amrex::IntVect(i-1,j),n) - 2.*ufab(amrex::IntVect(i,j),n))/dx[0]/dx[0],
	      (ufab(amrex::IntVect(i+1,j+1),n) + ufab(amrex::IntVect(i-1,j-1),n) - ufab(amrex::IntVect(i+1,j-1),n) - ufab(amrex::IntVect(i-1,j+1),n))/(2.*dx[0])/(2.*dx[1]),
	      (ufab(amrex::IntVect(i+1,j+1),n) + ufab(amrex::IntVect(i-1,j-1),n) - ufab(amrex::IntVect(i+1,j-1),n) - ufab(amrex::IntVect(i-1,j+1),n))/(2.*dx[0])/(2.*dx[1]),
	      (ufab(amrex::IntVect(i,j+1),n) + ufab(amrex::IntVect(i,j-1),n) - 2.*ufab(amrex::IntVect(i,j),n))/dx[1]/dx[1];
	    
	    ffab(amrex::IntVect(i,j),0) = - (gradepsilon[0](0,0) + gradepsilon[0](1,1));
	    ffab(amrex::IntVect(i,j),1) = - (gradepsilon[1](0,0) + gradepsilon[1](1,1));

	    continue;


	    for (int p = 0; p < AMREX_SPACEDIM; p++)
	      {
		ffab(amrex::IntVect(i,j),p) = 0.0;
		for (int q = 0; q < AMREX_SPACEDIM; q++)
		  {
		    ffab(amrex::IntVect(i,j),p) += 2.0*mu*gradepsilon[0](p,0);
		    if (p==q)
		      for (int k=0; k<AMREX_SPACEDIM; k++)
			ffab(amrex::IntVect(i,j)) += lambda * gradepsilon[q](k,k);
		  }
		  }	    
	  }
	  }
}

void
MLHeatConduction::Fsmooth (int amrlev,
			    int mglev,
			    MultiFab& sol,
			    const MultiFab& rhs,
			    int redblack) const
{
  BL_PROFILE("MLHeatConduction::Fsmooth()");

  OrientationIter oitr;

  const FabSet& f0 = m_undrrelxr[amrlev][mglev][oitr()]; ++oitr;
  const FabSet& f1 = m_undrrelxr[amrlev][mglev][oitr()]; ++oitr;
  const FabSet& f2 = m_undrrelxr[amrlev][mglev][oitr()]; ++oitr;
  const FabSet& f3 = m_undrrelxr[amrlev][mglev][oitr()]; ++oitr;

  const MultiMask& mm0 = m_maskvals [amrlev][mglev][0];
  const MultiMask& mm1 = m_maskvals [amrlev][mglev][1];
  const MultiMask& mm2 = m_maskvals [amrlev][mglev][2];
  const MultiMask& mm3 = m_maskvals [amrlev][mglev][3];

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
      const FArrayBox& afab    = m_a_coeffs[amrlev][mglev][mfi];

      AMREX_D_TERM(const FArrayBox& bxfab = m_b_coeffs[amrlev][mglev][0][mfi];,
		   const FArrayBox& byfab = m_b_coeffs[amrlev][mglev][1][mfi];,
		   const FArrayBox& bzfab = m_b_coeffs[amrlev][mglev][2][mfi];);

      
      amrex::Real alpha = m_a_scalar, beta=m_b_scalar;
      amrex::Real dhx = beta/h[0]/h[0];
      amrex::Real dhy = beta/h[1]/h[1];

      const FArrayBox& f0fab = f0[mfi];
      const FArrayBox& f1fab = f1[mfi];
      const FArrayBox& f2fab = f2[mfi];
      const FArrayBox& f3fab = f3[mfi];

      for (int n = 0; n < getNComp(); n++)
      	for (int j = tbx.loVect()[1]; j<=tbx.hiVect()[1]; j++)
       	  {
       	    int ioffset = (tbx.loVect()[0] + j + redblack)%2;
       	    for (int i = tbx.loVect()[0] + ioffset; i <= tbx.hiVect()[0]; i+= 2)
       	      {

       		amrex::Real cf0 = 0.0, cf1 = 0.0, cf2 = 0.0, cf3=0.0;

       		if ( (i == vbx.loVect()[0])  &&  (m0(amrex::IntVect(vbx.loVect()[0]-1,j)) > 0))
       		  cf0 = f0fab(amrex::IntVect(vbx.loVect()[0],j));

      		if ( (j == vbx.loVect()[1])  &&  (m1(amrex::IntVect(i,vbx.loVect()[1]-1)) > 0))
       		  cf1 = f1fab(amrex::IntVect(i,vbx.loVect()[1]));

       		if ( (i == vbx.hiVect()[0])  &&  (m2(amrex::IntVect(vbx.hiVect()[0]+1,j)) > 0))
       		  cf2 = f2fab(amrex::IntVect(vbx.hiVect()[0],j));

       		if ( (j == vbx.hiVect()[1])  &&  (m3(amrex::IntVect(i,vbx.hiVect()[1]+1)) > 0))
       		  cf3 = f3fab(amrex::IntVect(i,vbx.hiVect()[1]));

       		amrex::Real delta =
		  dhx*(bxfab(amrex::IntVect(i,j))*cf0 + bxfab(amrex::IntVect(i+1,j))*cf2) +
		  dhy*(byfab(amrex::IntVect(i,j))*cf1 + byfab(amrex::IntVect(i,j+1))*cf3);

		amrex::Real gamma =
		  alpha * afab(amrex::IntVect(i,j)) +
		  dhx*(bxfab(amrex::IntVect(i,j)) + bxfab(amrex::IntVect(i+1,j))) +
		  dhy*(byfab(amrex::IntVect(i,j)) + byfab(amrex::IntVect(i,j+1)));
		  

       		amrex::Real rho = 
       		  dhx*(bxfab(amrex::IntVect(i,j))  *solnfab(amrex::IntVect(i-1,j),n) +
		       bxfab(amrex::IntVect(i+1,j))*solnfab(amrex::IntVect(i+1,j),n)) + 
       		  dhy*(byfab(amrex::IntVect(i,j))  *solnfab(amrex::IntVect(i,j-1),n) +
		       byfab(amrex::IntVect(i,j+1))*solnfab(amrex::IntVect(i,j+1),n));

       		solnfab(amrex::IntVect(i,j),n) =
       		  (rhsfab(amrex::IntVect(i,j),n) + rho - solnfab(amrex::IntVect(i,j),n)*delta)
		  / (gamma - delta);

       	      }
       	  }
    }
}

void
MLHeatConduction::FFlux (int /*amrlev*/, const MFIter& /*mfi*/,
 			  const std::array<FArrayBox*,AMREX_SPACEDIM>& sigma,
 			  const FArrayBox& /*sol*/, const int /*face_only*/) const
{

  amrex::BaseFab<amrex::Real> &fxfab = *sigma[0];
  amrex::BaseFab<amrex::Real> &fyfab = *sigma[1];

  fxfab.setVal(0.0);
  fyfab.setVal(0.0);
}

int
MLHeatConduction::getNComp() const
{
  return AMREX_SPACEDIM;
}



}

#include <AMReX_MultiFabUtil.H>
#include "MLHeatConduction.H"
#include "AMReX_ABec_F.H"

namespace amrex {

MLHeatConduction::MLHeatConduction (const Vector<Geometry>& a_geom,
				      const Vector<BoxArray>& a_grids,
				      const Vector<DistributionMapping>& a_dmap,
				      const LPInfo& a_info)
{
  define(a_geom, a_grids, a_dmap, a_info);
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
MLHeatConduction::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
  BL_PROFILE("MLHeatConduction::Fapply()");

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
MLHeatConduction::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const
{
  BL_PROFILE("MLHeatConduction::Fsmooth()");

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

      /// \note GSRB = Gauss-Seidel Red Black
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

      // amrex::Real alpha = m_a_scalar, beta = m_b_scalar;
      // int lo[2] = {tbx.loVect()[0], tbx.loVect()[1]};
      // int hi[2] = {tbx.hiVect()[0], tbx.hiVect()[1]};
      // int blo[2] = {vbx.loVect()[0], vbx.loVect()[1]};
      // int bhi[2] = {vbx.hiVect()[0], vbx.hiVect()[1]};
      // int do_line;
      // int ilen,jlen;
      // static int LSDIM = 127;
      // amrex::Real a_ls[LSDIM];// REAL_T a_ls(0:LSDIM)
      // amrex::Real b_ls[LSDIM];// REAL_T b_ls(0:LSDIM)
      // amrex::Real c_ls[LSDIM];// REAL_T c_ls(0:LSDIM)
      // amrex::Real r_ls[LSDIM];// REAL_T r_ls(0:LSDIM)
      // amrex::Real u_ls[LSDIM];// REAL_T u_ls(0:LSDIM)
      // if (h[1] > 1.5*h[0]) {        // if (h(2). gt. 1.5D0*h(1)) then 
      // 	  do_line = 1;              // do_line = 1
      // 	  ilen = hi[0] - lo[0] + 1; // hi[0]-lo[0]+1;
      // 	  if (ilen > LSDIM)         // if (ilen .gt. LSDIM) then
      // 	    amrex::Abort("Too big for Line Solve in GSRB: ilen="+ilen);  // print *,'TOO BIG FOR LINE SOLVE IN GSRB: ilen = ',ilen
      //                                                                    // call bl_error("stop")
      // 	                                                                 // end if
      // }
      // else if (h[0] > 1.5*h[1]) {   // else if (h(1) .gt. 1.5D0*h(2)) then
      // 	  do_line = 2;              // do_line = 2
      // 	  jlen = hi[1] - lo[1] + 1; // jlen = hi(2)-lo(2)+1
      // 	  if (jlen > LSDIM)         // if (jlen .gt. LSDIM) then
      // 	    amrex::Abort("Too big for Line Solve in GSRB: ilen="+ilen);  // print *,'TOO BIG FOR LINE SOLVE IN GSRB: jlen = ',jlen
      //                                                                    // call bl_error("stop")
      //                                                                    // end if
      // }
      // else           //  else 
      // 	do_line = 0; // do_line = 0
      //                // end if
      // amrex::Real dhx, dhy;
      // dhx = beta/h[0]/h[0]; // dhx = beta/h(1)**2
      // dhy = beta/h[1]/h[1]; // dhy = beta/h(2)**2
      // for (int n = 0; n < nc; n++) // do n = 1, nc
      // 	if (do_line == 0)  {    //if (do_line .eq. 0) then
      // 	  for (int j = lo[1]; j <= hi[1]; j++){ //do j = lo(2), hi(2)
      // 	    int ioff = (lo[0] + j + redblack)%2; // ioff = MOD(lo(1) + j + redblack, 2)
      // 	    for (int i = lo[0]+ioff; i <= hi[0]; i += 2) {//do i = lo(1) + ioff,hi(1),2 
      // 	      amrex::Real cf0 = (i==blo[0])&&(m0(amrex::IntVect(blo[0]-1,j))>0) ? f0fab(amrex::IntVect(blo[0],j)) : 0.0; // cf0 = merge(f0(blo(1),j), 0.0D0,  (i .eq. blo(1)) .and. (m0(blo(1)-1,j).gt.0))
      // 	      amrex::Real cf1 = (j==blo[1])&&(m1(amrex::IntVect(i,blo[1]-1))>0) ? f1fab(amrex::IntVect(i,blo[1])) : 0.0; // cf1 = merge(f1(i,blo(2)), 0.0D0,  (j .eq. blo(2)) .and. (m1(i,blo(2)-1).gt.0))
      // 	      amrex::Real cf2 = (i==bhi[0])&&(m2(amrex::IntVect(bhi[0]+1,j))>0) ? f2fab(amrex::IntVect(bhi[0],j)) : 0.0; // cf2 = merge(f2(bhi(1),j), 0.0D0, (i .eq. bhi(1)) .and. (m2(bhi(1)+1,j).gt.0))
      // 	      amrex::Real cf3 = (j==bhi[1])&&(m3(amrex::IntVect(i,bhi[1]+1))>0) ? f3fab(amrex::IntVect(i,bhi[1])) : 0.0; // cf3 = merge(f3(i,bhi(2)), 0.0D0, (j .eq. bhi(2)) .and. (m3(i,bhi(2)+1).gt.0))
      // 	      amrex::Real delta =
      // 		dhx * (bxfab(amrex::IntVect(i,j))*cf0 + bxfab(amrex::IntVect(i+1,j))*cf2) +  // delta = dhx*(bX(i,j)*cf0 + bX(i+1,j)*cf2)  +  dhy*(bY(i,j)*cf1 + bY(i,j+1)*cf3)
      // 		dhy * (byfab(amrex::IntVect(i,j))*cf1 + byfab(amrex::IntVect(i,j+1))*cf3);
      // 	      amrex::Real gamma =   // gamma = alpha*a(i,j) +   dhx*( bX(i,j) + bX(i+1,j) ) +   dhy*( bY(i,j) + bY(i,j+1) )
      // 		alpha*afab(amrex::IntVect(i,j)) +
      // 		dhx*(bxfab(amrex::IntVect(i,j)) + bxfab(amrex::IntVect(i+1,j))) +
      // 		dhy*(byfab(amrex::IntVect(i,j)) + byfab(amrex::IntVect(i,j+1)));
      // 	      amrex::Real rho =  //rho = dhx*(bX(i,j)*phi(i-1,j,n) + bX(i+1,j)*phi(i+1,j,n))+dhy*(bY(i,j)*phi(i,j-1,n) + bY(i,j+1)*phi(i,j+1,n))
      // 		dhx*(bxfab(amrex::IntVect(i,j))*solnfab(amrex::IntVect(i-1,j),n) + bxfab(amrex::IntVect(i+1,j))*solnfab(amrex::IntVect(i+1,j),n)) +
      // 		dhy*(byfab(amrex::IntVect(i,j))*solnfab(amrex::IntVect(i,j-1),n) + byfab(amrex::IntVect(i,j+1))*solnfab(amrex::IntVect(i,j+1),n));

      // 	      solnfab(amrex::IntVect(i,j),n) = (rhsfab(amrex::IntVect(i,j),n) + rho - rhsfab(amrex::IntVect(i,j),n)*delta) / (gamma-delta); //phi(i,j,n) = (rhs(i,j,n) + rho - phi(i,j,n)*delta)/(gamma - delta)
      // 	    }   //  end do
      // 	  }  // end do
      // 	}
      // 	else if (do_line == 2) {  // else if (do_line .eq. 2) then
      // 	  int ioff = (lo[0]+redblack)%2; //          ioff = MOD(lo(1) + redblack, 2)
      // 	  for (int i = lo[0]+ioff; i <= hi[0]; i+=2) {//          do i = lo(1) + ioff,hi(1),2
      // 	    for (int j = lo[1]; j<=hi[1]; j++) {//              do j = lo(2), hi(2)
      // 	      amrex::Real cf0 = (i==blo[0])&&(m0(amrex::IntVect(blo[0]-1,j))>0) ? f0fab(amrex::IntVect(blo[0],j)) : 0.0;//cf0 = merge(f0(blo(1),j), 0.0D0, (i .eq. blo(1)) .and. (m0(blo(1)-1,j).gt.0))
      // 	      amrex::Real cf1 = (j==blo[1])&&(m1(amrex::IntVect(i,blo[1]-1))>0) ? f1fab(amrex::IntVect(i,blo[1])) : 0.0;//cf1 = merge(f1(i,blo(2)), 0.0D0, (j .eq. blo(2)) .and. (m1(i,blo(2)-1).gt.0))
      // 	      amrex::Real cf2 = (i==bhi[0])&&(m2(amrex::IntVect(bhi[0]+1,j))>0) ? f2fab(amrex::IntVect(bhi[0],j)) : 0.0;//cf2 = merge(f2(bhi(1),j), 0.0D0, (i .eq. bhi(1)) .and. (m2(bhi(1)+1,j).gt.0))
      // 	      amrex::Real cf3 = (j==bhi[1])&&(m3(amrex::IntVect(i,bhi[1]+1))>0) ? f3fab(amrex::IntVect(i,bhi[1])) : 0.0;//cf3 = merge(f3(i,bhi(2)), 0.0D0, (j .eq. bhi(2)) .and. (m3(i,bhi(2)+1).gt.0))
      // 	      amrex::Real delta =
      // 		dhx * (bxfab(amrex::IntVect(i,j))*cf0 + bxfab(amrex::IntVect(i+1,j))*cf2) +  // delta = dhx*(bX(i,j)*cf0 + bX(i+1,j)*cf2)  +  dhy*(bY(i,j)*cf1 + bY(i,j+1)*cf3)
      // 		dhy * (byfab(amrex::IntVect(i,j))*cf1 + byfab(amrex::IntVect(i,j+1))*cf3);
      // 	      amrex::Real gamma =   // gamma = alpha*a(i,j) +   dhx*( bX(i,j) + bX(i+1,j) ) +   dhy*( bY(i,j) + bY(i,j+1) )
      // 		alpha*afab(amrex::IntVect(i,j)) +
      // 		dhx*(bxfab(amrex::IntVect(i,j)) + bxfab(amrex::IntVect(i+1,j))) +
      // 		dhy*(byfab(amrex::IntVect(i,j)) + byfab(amrex::IntVect(i,j+1)));
      // 	      amrex::Real rho_x =  // rho_x = dhx*(bX(i,j)*phi(i-1,j,n) + bX(i+1,j)*phi(i+1,j,n)=))
      // 		dhx*(bxfab(amrex::IntVect(i,j))*solnfab(amrex::IntVect(i-1,j),n) + bxfab(amrex::IntVect(i+1,j))*solnfab(amrex::IntVect(i+1,j),n));
      // 	      a_ls[j-lo[1]] = -dhy*byfab(amrex::IntVect(i,j)); //  a_ls(j-lo(2)) = -dhy*bY(i,j)
      // 	      b_ls[j-lo[1]] = gamma-delta; // b_ls(j-lo(2)) = gamma - delta
      // 	      c_ls[j-lo[1]] = -dhy*byfab(amrex::IntVect(i,j+1)); // c_ls(j-lo(2)) = -dhy*bY(i,j+1)
      // 	      r_ls[j-lo[1]] = rhsfab(amrex::IntVect(i,j),n) + rho_x - solnfab(amrex::IntVect(i,j),n)*delta; // r_ls(j-lo(2)) = rhs(i,j,n) + rho_x - phi(i,j,n)*delta
      // 	      if (j==lo[1])      //  if (j .eq. lo(2)) 
      // 		r_ls[j-lo[1]] =  r_ls[j-lo[1]] + dhy*byfab(amrex::IntVect(i,j))*solnfab(amrex::IntVect(i,j-1),n); //  r_ls(j-lo(2)) = r_ls(j-lo(2)) + dhy*bY(i,j)*phi(i,j-1,n)
      // 	      if (j==hi[1])  //if (j .eq. hi(2)) 
      // 		r_ls[j-lo[1]] = r_ls[j-lo[1]] + dhy*byfab(amrex::IntVect(i,j+1))*solnfab(amrex::IntVect(i,j+1),n); // r_ls(j-lo(2)) = r_ls(j-lo(2)) + dhy*bY(i,j+1)*phi(i,j+1,n)
      // 	    } //end do
      // 	    // TODO TODO TODO TODO TODO // call tridiag(a_ls,b_ls,c_ls,r_ls,u_ls,jlen)
      // 	    for (int j = lo[1]; j <= hi[1]; j++) {  // do j = lo(2), hi(2)
      // 	      solnfab(amrex::IntVect(i,j),n) = u_ls[j-lo[1]]; //  phi(i,j,n) = u_ls(j-lo(2))
      // 	    }// end do
      // 	  } // end do
      // 	}
      // 	else if (do_line == 1) {// else if (do_line .eq. 1) then
      // 	  int joff = (lo[1] + redblack)%2; // joff = MOD(lo(2) + redblack, 2)
      // 	  for (int j=lo[1]+joff; j<= hi[1]; j+= 2) { // do j = lo(2) + joff,hi(2),2
      // 	    for (int i=lo[0]; i<= hi[0]; i++) {// do i = lo(1), hi(1)
      // 	      amrex::Real cf0 = (i==blo[0])&&(m0(amrex::IntVect(blo[0]-1,j))>0) ? f0fab(amrex::IntVect(blo[0],j)) : 0.0;//cf0 = merge(f0(blo(1),j), 0.0D0, (i .eq. blo(1)) .and. (m0(blo(1)-1,j).gt.0))
      // 	      amrex::Real cf1 = (j==blo[1])&&(m1(amrex::IntVect(i,blo[1]-1))>0) ? f1fab(amrex::IntVect(i,blo[1])) : 0.0;//cf1 = merge(f1(i,blo(2)), 0.0D0, (j .eq. blo(2)) .and. (m1(i,blo(2)-1).gt.0))
      // 	      amrex::Real cf2 = (i==bhi[0])&&(m2(amrex::IntVect(bhi[0]+1,j))>0) ? f2fab(amrex::IntVect(bhi[0],j)) : 0.0;//cf2 = merge(f2(bhi(1),j), 0.0D0, (i .eq. bhi(1)) .and. (m2(bhi(1)+1,j).gt.0))
      // 	      amrex::Real cf3 = (j==bhi[1])&&(m3(amrex::IntVect(i,bhi[1]+1))>0) ? f3fab(amrex::IntVect(i,bhi[1])) : 0.0;//cf3 = merge(f3(i,bhi(2)), 0.0D0, (j .eq. bhi(2)) .and. (m3(i,bhi(2)+1).gt.0))
      // 	      amrex::Real delta =
      // 		dhx * (bxfab(amrex::IntVect(i,j))*cf0 + bxfab(amrex::IntVect(i+1,j))*cf2) +  // delta = dhx*(bX(i,j)*cf0 + bX(i+1,j)*cf2)  +  dhy*(bY(i,j)*cf1 + bY(i,j+1)*cf3)
      // 		dhy * (byfab(amrex::IntVect(i,j))*cf1 + byfab(amrex::IntVect(i,j+1))*cf3);
      // 	      amrex::Real gamma =   // gamma = alpha*a(i,j) +   dhx*( bX(i,j) + bX(i+1,j) ) +   dhy*( bY(i,j) + bY(i,j+1) )
      // 		alpha*afab(amrex::IntVect(i,j)) +
      // 		dhx*(bxfab(amrex::IntVect(i,j)) + bxfab(amrex::IntVect(i+1,j))) +
      // 		dhy*(byfab(amrex::IntVect(i,j)) + byfab(amrex::IntVect(i,j+1)));
      // 	      amrex::Real rho_y =   // rho_y = dhy*(bY(i,j)*phi(i,j-1,n) + bY(i,j+1)*phi(i,j+1,n))
      // 		dhy * (byfab(amrex::IntVect(i,j))*solnfab(amrex::IntVect(i,j-1),n) + byfab(amrex::IntVect(i,j+1))*solnfab(amrex::IntVect(i,j+1),n));
      // 	      a_ls[i-lo[0]] = -dhx*bxfab(amrex::IntVect(i,j)); // a_ls(i-lo(1)) = -dhx*bX(i,j)
      // 	      b_ls[i-lo[0]] = gamma-delta; // b_ls(i-lo(1)) = gamma - delta
      // 	      c_ls[i-lo[0]] = -dhx*byfab(amrex::IntVect(i+1,j)); // c_ls(i-lo(1)) = -dhx*bX(i+1,j)
      // 	      r_ls[i-lo[0]] = rhsfab(amrex::IntVect(i,j),n) + rho_y - solnfab(amrex::IntVect(i,j),n)*delta; // r_ls(i-lo(1)) = rhs(i,j,n) + rho_y - phi(i,j,n)*delta
      // 	      if (i==lo[0]) // if (i .eq. lo(1)) 
      // 		r_ls[i-lo[0]] = r_ls[i-lo[0]] + dhx*bxfab(amrex::IntVect(i,j))*solnfab(amrex::IntVect(i-1,j),n); // r_ls(i-lo(1)) = r_ls(i-lo(1)) + dhx*bX(i,j)*phi(i-1,j,n)
      // 	      if (i==hi[0]) // if (i .eq. hi(1)) 
      // 		r_ls[i-lo[0]] = r_ls[i-lo[0]] + dhx*bxfab(amrex::IntVect(i+1,j))*solnfab(amrex::IntVect(i+1,j),n); //  r_ls(i-lo(1)) = r_ls(i-lo(1)) + dhx*bX(i+1,j)*phi(i+1,j,n)
      // 	    }//  end do
      // 	    // TODO TODO TODO TODO TODO // call tridiag(a_ls,b_ls,c_ls,r_ls,u_ls,ilen)
      // 	    for (int i=lo[0]; i<= hi[0]; i++) { // do i = lo(1), hi(1)
      // 	      solnfab(amrex::IntVect(i,j),n) = u_ls[i-lo[0]]; // phi(i,j,n) = u_ls(i-lo(1))
      // 	    } //  end do
      // 	  } // end do
      // 	}
      // 	else {//        else
      // 	  amrex::Abort("Bogus Do Line"); //          print *,'BOGUS DO_LINE '
      // 	                                 //          call bl_error("stop")
      // 	} //        end if
    }
}

void
MLHeatConduction::FFlux (int amrlev, const MFIter& mfi,
 			  const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
 			  const FArrayBox& sol, const int face_only) const
{
  // std::cout << "in FFlux, face_only=" << face_only << std::endl;
  // BL_PROFILE("MLHeatConduction::FFlux()");

  // const int mglev = 0;
  // const auto& bx = m_b_coeffs[amrlev][mglev][0][mfi];
  // const auto& by = m_b_coeffs[amrlev][mglev][1][mfi];

  // const Box& box = mfi.tilebox();
  // const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

  // amrex::Real beta=m_b_scalar;
  
  // amrex::Real dhx = beta*dxinv[0], dhy = beta*dxinv[1];

  // amrex::BaseFab<amrex::Real> &fxfab = *flux[0];
  // amrex::BaseFab<amrex::Real> &fyfab = *flux[1];
  // const amrex::BaseFab<amrex::Real> &solfab = sol;



  // if (face_only)
  //   {
  //     for (int i = box.loVect()[0]; i<=box.hiVect()[0]+1; i+= 1 + (box.hiVect()[0]-box.loVect()[0]))
  // 	for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++)
  // 	  {
  // 	    //fxfab(amrex::IntVect(i,j)) = 0.0;
  // 	      // -dhx *
  // 	      // bx(amrex::IntVect(i,j)) *
  // 	      // (solfab(amrex::IntVect(i,j))-solfab(amrex::IntVect(i,j)));
  // 	  }
  //     for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++)
  // 	for (int j = box.loVect()[1]; j<=box.hiVect()[1]+1; j+= 1 + (box.hiVect()[1]-box.loVect()[1]))
  // 	  {
  // 	    //fyfab(amrex::IntVect(i,j)) = 0.0;
  // 	      // -dhy *
  // 	      // by(amrex::IntVect(i,j)) *
  // 	      // (solfab(amrex::IntVect(i,j))-solfab(amrex::IntVect(i,j-1)));
  // 	  }
  //   }
  // else
  //   {
  //     for (int i = box.loVect()[0]; i<=box.hiVect()[0]+1; i++)
  // 	for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++)
  // 	  {
  // 	    //fxfab(amrex::IntVect(i,j)) = 0.0;
  // 	      // -dhx *
  // 	      // bx(amrex::IntVect(i,j)) *
  // 	      // (solfab(amrex::IntVect(i,j))-solfab(amrex::IntVect(i,j)));
  // 	  }
       
  //     for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++)
  // 	for (int j = box.loVect()[1]; j<=box.hiVect()[1]+1; j++)
  // 	  {
  // 	    //fyfab(amrex::IntVect(i,j)) = 0.0;
  // 	      // -dhy *
  // 	      // by(amrex::IntVect(i,j)) *
  // 	      // (solfab(amrex::IntVect(i,j))-solfab(amrex::IntVect(i,j-1)));
  // 	  }
  //   }
}

}


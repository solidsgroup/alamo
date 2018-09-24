
#include "Operator.H"
#include <AMReX_MLLinOp_F.H>
#include <AMReX_MultiFabUtil.H>
#include "Util/Color.H"


#define TRACER	//std::cout << Color::FG::Yellow << __FILE__ << ":" << __LINE__ << Color::FG::Default << " " << __func__ << std::endl;
#define PROBE	//std::cout << Color::FG::Red << __FILE__ << ":" << __LINE__ << Color::FG::Default << " " << __func__ << std::endl;


using namespace amrex;
namespace Operator {


void Operator::Diagonal (int amrlev,
			 int mglev,
			 amrex::MultiFab& diag) const
{
	int ncomp = diag.nComp();
	int nghost = diag.nGrow();
	amrex::MultiFab x(diag.boxArray(), diag.DistributionMap(), ncomp, nghost);
	amrex::MultiFab Ax(diag.boxArray(), diag.DistributionMap(), ncomp, nghost);

	x.setVal(0.0);
	Ax.setVal(0.0);
	diag.setVal(0.0);

	for (MFIter mfi(x, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		amrex::FArrayBox &diagfab = diag[mfi];
		amrex::FArrayBox       &xfab    = x[mfi];
		amrex::FArrayBox       &Axfab   = Ax[mfi];

		AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
			     for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
			     for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
		{
			amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
			for (int i = 0; i < ncomp; ++i)
			{
				xfab(m,i) = 1.0;
				Fapply(amrlev,mglev,Ax,x);
				diagfab(m,i) = amrex::MultiFab::Dot(x, 0, Ax, 0, ncomp, nghost);
				xfab.setVal(0.0);
				Axfab.setVal(0.0);
			}
		}
	}

}


Operator::Operator (const Vector<Geometry>& a_geom,
		    const Vector<BoxArray>& a_grids,
		    const Vector<DistributionMapping>& a_dmap,
		    const LPInfo& a_info,
		    const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
	TRACER;
	define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

 Operator::~Operator ()
 {}

 void
	 Operator::define (const Vector<Geometry>& a_geom,
		const Vector<BoxArray>& a_grids,
		const Vector<DistributionMapping>& a_dmap,
		const LPInfo& a_info,
		const Vector<FabFactory<FArrayBox> const*>& a_factory)
 {
	 TRACER;
	 BL_PROFILE("MLNodeLaplacian::define()");

	 // This makes sure grids are cell-centered;
	 Vector<BoxArray> cc_grids = a_grids;
	 for (auto& ba : cc_grids)
		 ba.enclosedCells();
	
	 MLNodeLinOp::define(a_geom, cc_grids, a_dmap, a_info, a_factory);

	 // m_sigma.resize(m_num_amr_levels);
	 // for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	 // {
	 // 	 m_sigma[amrlev].resize(m_num_mg_levels[amrlev]);
	 // 	 const int mglev = 0;
	 // 	 const int idim = 0;
	 // 	 m_sigma[amrlev][mglev][idim].reset
	 // 		 (new MultiFab(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 1));
	 // 	 m_sigma[amrlev][mglev][idim]->setVal(0.0);
	 // }

#if (AMREX_SPACEDIM == 2)
	 m_is_rz = Geometry::IsRZ();
#endif
 }

void
Operator::setSigma (int amrlev, const MultiFab& a_sigma)
{
	TRACER;
	MultiFab::Copy(*m_sigma[amrlev][0][0], a_sigma, 0, 0, 1, 0);
}

void
Operator::compRHS (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel,
		   const Vector<const MultiFab*>& rhnd,
		   const Vector<MultiFab*>& a_rhcc)
{
	Util::Abort("here in compRHS");
}

void
Operator::updateVelocity (const Vector<MultiFab*>& vel, const Vector<MultiFab const*>& sol) const
{
	TRACER;
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		const auto& sigma = *m_sigma[amrlev][0][0];
		const Real* dxinv = m_geom[amrlev][0].InvCellSize();
		for (MFIter mfi(*vel[amrlev], true); mfi.isValid(); ++mfi)
		{
			const Box& bx = mfi.tilebox();
			auto& vfab = (*vel[amrlev])[mfi];
			{
				amrex_mlndlap_mknewu(BL_TO_FORTRAN_BOX(bx),
						     BL_TO_FORTRAN_ANYD(vfab),
						     BL_TO_FORTRAN_ANYD((*sol[amrlev])[mfi]),
						     BL_TO_FORTRAN_ANYD(sigma[mfi]),
						     dxinv);
			}
		}
	}
}

void
Operator::averageDownCoeffs ()
{
	TRACER;
	// BL_PROFILE("Operator::averageDownCoeffs()");

	// if (m_coarsening_strategy == CoarseningStrategy::Sigma)
	// {
	// 	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	// 	{
	// 		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
	// 		{
	// 			int ndims = (m_use_harmonic_average) ? AMREX_SPACEDIM : 1;
	// 			for (int idim = 0; idim < ndims; ++idim)
	// 			{
	// 				if (m_sigma[amrlev][mglev][idim] == nullptr) {
	// 					if (mglev == 0) {
	// 						m_sigma[amrlev][mglev][idim].reset
	// 							(new MultiFab(*m_sigma[amrlev][mglev][0], amrex::make_alias, 0, 1));
	// 					} else {
	// 						m_sigma[amrlev][mglev][idim].reset
	// 							(new MultiFab(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 1));
	// 						m_sigma[amrlev][mglev][idim]->setVal(0.0);
	// 					}
	// 				}
	// 			}
	// 		}
	// 	}
	// }

	// for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
	// {
	// 	averageDownCoeffsSameAmrLevel(amrlev);
	// 	averageDownCoeffsToCoarseAmrLevel(amrlev);
	// }

	// averageDownCoeffsSameAmrLevel(0);

	// for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	// {
	//     if (m_use_harmonic_average) {
	//         int mglev = 0;
	//         FillBoundaryCoeff(*m_sigma[amrlev][mglev][0], m_geom[amrlev][mglev]);
	//         for (mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev)
	//         {
	//             for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
	//                 if (m_sigma[amrlev][mglev][idim]) {
	//                     FillBoundaryCoeff(*m_sigma[amrlev][mglev][idim], m_geom[amrlev][mglev]);
	//                 }
	//             }
	//         }
	//     } else {
	//         int idim = 0;
	//         for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
	//         {
	//             if (m_sigma[amrlev][mglev][idim]) {
	//                 FillBoundaryCoeff(*m_sigma[amrlev][mglev][idim], m_geom[amrlev][mglev]);
	//             }
	//         }
	//     }
	// }
}

void
Operator::averageDownCoeffsToCoarseAmrLevel (int flev)
{
	TRACER;
	const int mglev = 0;
	const int idim = 0;  // other dimensions are just aliases
	amrex::average_down(*m_sigma[flev][mglev][idim], *m_sigma[flev-1][mglev][idim], 0, 1,
			    m_amr_ref_ratio[flev-1]);
}

void
Operator::averageDownCoeffsSameAmrLevel (int amrlev)
{
	TRACER;
	if (m_coarsening_strategy != CoarseningStrategy::Sigma) return;

	const int nsigma = (m_use_harmonic_average) ? AMREX_SPACEDIM : 1;

	for (int mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev)
	{
		for (int idim = 0; idim < nsigma; ++idim)
		{
			const MultiFab& fine = *m_sigma[amrlev][mglev-1][idim];
			MultiFab& crse = *m_sigma[amrlev][mglev][idim];
			bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
			MultiFab cfine;
			if (need_parallel_copy) {
				const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
				cfine.define(ba, fine.DistributionMap(), 1, 0);
			}

			MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;

#ifdef _OPENMP
#pragma omp parallel
#endif
			for (MFIter mfi(*pcrse, true); mfi.isValid(); ++mfi)
			{
				const Box& bx = mfi.tilebox();
				amrex_mlndlap_avgdown_coeff(BL_TO_FORTRAN_BOX(bx),
							    BL_TO_FORTRAN_ANYD((*pcrse)[mfi]),
							    BL_TO_FORTRAN_ANYD(fine[mfi]),
							    &idim);
			}

			if (need_parallel_copy) {
				crse.ParallelCopy(cfine);
			}
		}
	}
}

void
Operator::FillBoundaryCoeff (MultiFab& sigma, const Geometry& geom)
{
	TRACER;
	//     BL_PROFILE("Operator::FillBoundaryCoeff()");

	//     sigma.FillBoundary(geom.periodicity());

	//     if (m_coarsening_strategy == CoarseningStrategy::Sigma)
	//     {
	//         const Box& domain = geom.Domain();

	// #ifdef _OPENMP
	// #pragma omp parallel
	// #endif
	//         for (MFIter mfi(sigma, MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
	//         {
	//             if (!domain.contains(mfi.fabbox()))
	//             {
	//                 amrex_mlndlap_fillbc_cc(BL_TO_FORTRAN_ANYD(sigma[mfi]),
	//                                         BL_TO_FORTRAN_BOX(domain),
	//                                         m_lobc.data(), m_hibc.data());
	//             }
	//         }
	//     }
}

void
Operator::buildMasks ()
{
	TRACER;
	if (m_masks_built) return;
	BL_PROFILE("Operator::buildMasks()");

	m_masks_built = true;

	m_is_bottom_singular = false;
	auto itlo = std::find(m_lobc.begin(), m_lobc.end(), BCType::Dirichlet);
	auto ithi = std::find(m_hibc.begin(), m_hibc.end(), BCType::Dirichlet);

	if (itlo == m_lobc.end() && ithi == m_hibc.end())
	{  // No Dirichlet
		PROBE;
		/// \todo need to work out BCs a little more rigorously...
		//m_is_bottom_singular = m_domain_covered[0];
	}

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		std::vector< std::pair<int,Box> > isects;
		IArrayBox ccfab;

		for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
		{
			for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
			{
				const Geometry& geom = m_geom[amrlev][mglev];
				const auto& period = geom.periodicity();
				const Box& ccdomain = geom.Domain();
				const Box& nddomain = amrex::surroundingNodes(ccdomain);
				const std::vector<IntVect>& pshifts = period.shiftIntVect();

				Box ccdomain_p = ccdomain;
				for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
					if (Geometry::isPeriodic(idim)) {
						ccdomain_p.grow(idim, 1);
					}
				}

				{
					auto& dmask = *m_dirichlet_mask[amrlev][mglev];
					const BoxArray& ccba = m_grids[amrlev][mglev];

					for (MFIter mfi(dmask, MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
					{
						const Box& ndbx = mfi.validbox();
						const Box& ccbx = amrex::enclosedCells(ndbx);
						const Box& ccbxg1 = amrex::grow(ccbx,1);
						IArrayBox& mskfab = dmask[mfi];
                        
						ccfab.resize(ccbxg1);
						ccfab.setVal(1);
						ccfab.setComplement(2,ccdomain_p,0,1);

						for (const auto& iv : pshifts)
						{
							ccba.intersections(ccbxg1+iv, isects);
							for (const auto& is : isects)
							{
								ccfab.setVal(0, is.second-iv, 0, 1);
							}
						}
                        
						amrex_mlndlap_set_dirichlet_mask(BL_TO_FORTRAN_ANYD(mskfab),
										 BL_TO_FORTRAN_ANYD(ccfab),
										 BL_TO_FORTRAN_BOX(nddomain),
										 m_lobc.data(), m_hibc.data());
					}
				}
			}
		}
	}

	for (int amrlev = 0; amrlev < m_num_amr_levels-1; ++amrlev)
	{
		iMultiFab& cc_mask = *m_cc_fine_mask[amrlev];
		iMultiFab& nd_mask = *m_nd_fine_mask[amrlev];
		LayoutData<int>& has_cf = *m_has_fine_bndry[amrlev];
		const BoxArray& fba = m_grids[amrlev+1][0];
		const BoxArray& cfba = amrex::coarsen(fba, AMRRefRatio(amrlev));

		const Box& ccdom = m_geom[amrlev][0].Domain();

		AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMRRefRatio(amrlev) == 2, "ref_ratio != 0 not supported");

		cc_mask.setVal(0);  // coarse by default

		const std::vector<IntVect>& pshifts = m_geom[amrlev][0].periodicity().shiftIntVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
		{
			std::vector< std::pair<int,Box> > isects;

			for (MFIter mfi(cc_mask); mfi.isValid(); ++mfi)
			{
				has_cf[mfi] = 0;
				IArrayBox& fab = cc_mask[mfi];
				const Box& bx = fab.box();
				for (const auto& iv : pshifts)
				{
					cfba.intersections(bx+iv, isects);
					for (const auto& is : isects)
					{
						fab.setVal(1, is.second-iv, 0, 1);
					}
					if (!isects.empty()) has_cf[mfi] = 1;
				}

				amrex_mlndlap_fillbc_cc_i(BL_TO_FORTRAN_ANYD(fab),
							  BL_TO_FORTRAN_BOX(ccdom),
							  m_lobc.data(), m_hibc.data());
			}
		}

#ifdef _OPENMP
#pragma omp parallel
#endif
		for (MFIter mfi(nd_mask,true); mfi.isValid(); ++mfi)
		{
			const Box& bx = mfi.tilebox();
			amrex_mlndlap_set_nodal_mask(BL_TO_FORTRAN_BOX(bx),
						     BL_TO_FORTRAN_ANYD(nd_mask[mfi]),
						     BL_TO_FORTRAN_ANYD(cc_mask[mfi]));
		}
	}

	auto& has_cf = *m_has_fine_bndry[m_num_amr_levels-1];
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(has_cf); mfi.isValid(); ++mfi)
	{
		has_cf[mfi] = 0;
	}

	{
		int amrlev = 0;
		int mglev = m_num_mg_levels[amrlev]-1;
		const iMultiFab& omask = *m_owner_mask[amrlev][mglev];
		m_bottom_dot_mask.define(omask.boxArray(), omask.DistributionMap(), 1, 0);

		const Geometry& geom = m_geom[amrlev][mglev];
		Box nddomain = amrex::surroundingNodes(geom.Domain());

		if (m_coarsening_strategy != CoarseningStrategy::Sigma) {
			nddomain.grow(1000); // hack to avoid masks being modified at Neuman boundary
		}

#ifdef _OPENMP
#pragma omp parallel
#endif
		for (MFIter mfi(m_bottom_dot_mask,true); mfi.isValid(); ++mfi)
		{
			const Box& bx = mfi.tilebox();
			auto& dfab = m_bottom_dot_mask[mfi];
			const auto& sfab = omask[mfi];
			amrex_mlndlap_set_dot_mask(BL_TO_FORTRAN_BOX(bx),
						   BL_TO_FORTRAN_ANYD(dfab),
						   BL_TO_FORTRAN_ANYD(sfab),
						   BL_TO_FORTRAN_BOX(nddomain),
						   m_lobc.data(), m_hibc.data());
		}
	}

	if (m_is_bottom_singular)
	{
		int amrlev = 0;
		int mglev = 0;
		const iMultiFab& omask = *m_owner_mask[amrlev][mglev];
		TRACER;TRACER;TRACER;TRACER;TRACER;TRACER;TRACER;
		m_coarse_dot_mask.define(omask.boxArray(), omask.DistributionMap(), 1, 0);

		const Geometry& geom = m_geom[amrlev][mglev];
		Box nddomain = amrex::surroundingNodes(geom.Domain());

		if (m_coarsening_strategy != CoarseningStrategy::Sigma) {
			nddomain.grow(1000); // hack to avoid masks being modified at Neuman boundary
		}

#ifdef _OPENMP
#pragma omp parallel
#endif
		for (MFIter mfi(m_coarse_dot_mask,true); mfi.isValid(); ++mfi)
		{
			const Box& bx = mfi.tilebox();
			auto& dfab = m_coarse_dot_mask[mfi];
			const auto& sfab = omask[mfi];
			amrex_mlndlap_set_dot_mask(BL_TO_FORTRAN_BOX(bx),
						   BL_TO_FORTRAN_ANYD(dfab),
						   BL_TO_FORTRAN_ANYD(sfab),
						   BL_TO_FORTRAN_BOX(nddomain),
						   m_lobc.data(), m_hibc.data());
		}
	}
}

void
Operator::buildStencil ()
{
	TRACER;
	// currently using Sigma coarsening strategy
	m_stencil.resize(m_num_amr_levels);
	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		m_stencil[amrlev].resize(m_num_mg_levels[amrlev]);
	}
    
	if (m_coarsening_strategy != CoarseningStrategy::RAP) return;

	const int ncomp_s = (AMREX_SPACEDIM == 2) ? 5 : 9;
	const int ncomp_c = (AMREX_SPACEDIM == 2) ? 6 : 12;
	AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMREX_SPACEDIM != 1,
					 "Operator::buildStencil: 1d not supported");
	AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!Geometry::IsRZ(),
					 "Operator::buildStencil: cylindrical not supported for RAP");

	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
		{
			m_stencil[amrlev][mglev].reset
				(new MultiFab(amrex::convert(m_grids[amrlev][mglev],
							     IntVect::TheNodeVector()),
					      m_dmap[amrlev][mglev], ncomp_s, 4));
			m_stencil[amrlev][mglev]->setVal(0.0);
		}

		{
			const Geometry& geom = m_geom[amrlev][0];
			const Real* dxinv = geom.InvCellSize();


#ifdef _OPENMP
#pragma omp parallel
#endif
			{
				FArrayBox sgfab;
				FArrayBox cnfab;
				for (MFIter mfi(*m_stencil[amrlev][0], MFItInfo().EnableTiling().SetDynamic(true));
				     mfi.isValid(); ++mfi)
				{
					Box vbx = mfi.validbox();
					AMREX_D_TERM(vbx.growLo(0,1);, vbx.growLo(1,1);, vbx.growLo(2,1));
					Box bx = mfi.growntilebox(1);
					bx &= vbx;
					const Box& ccbx = amrex::enclosedCells(bx);
					FArrayBox& stfab = (*m_stencil[amrlev][0])[mfi];
					const FArrayBox& sgfab_orig = (*m_sigma[amrlev][0][0])[mfi];
                    
					{
						Box bx2 = amrex::grow(ccbx,1);
						sgfab.resize(bx2);
						sgfab.setVal(0.0);
						bx2 &= sgfab_orig.box();
						sgfab.copy(sgfab_orig, bx2, 0, bx2, 0, 1);
						amrex_mlndlap_set_stencil(BL_TO_FORTRAN_BOX(bx),
									  BL_TO_FORTRAN_ANYD(stfab),
									  BL_TO_FORTRAN_ANYD(sgfab),
									  dxinv);
					}
				}

				for (MFIter mfi(*m_stencil[amrlev][0],true); mfi.isValid(); ++mfi)
				{
					const Box& bx = mfi.tilebox();
					FArrayBox& stfab = (*m_stencil[amrlev][0])[mfi];
                    
					amrex_mlndlap_set_stencil_s0(BL_TO_FORTRAN_BOX(bx),
								     BL_TO_FORTRAN_ANYD(stfab));
				}
			}

			m_stencil[amrlev][0]->FillBoundary(geom.periodicity());
		}

		for (int mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev)
		{
			const MultiFab& fine = *m_stencil[amrlev][mglev-1];
			MultiFab& crse = *m_stencil[amrlev][mglev];
			bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
			MultiFab cfine;
			if (need_parallel_copy) {
				const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
				cfine.define(ba, fine.DistributionMap(), fine.nComp(), 1);
				cfine.setVal(0.0);
			}

			MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;

#ifdef _OPENMP
#pragma omp parallel
#endif
			{
				for (MFIter mfi(*pcrse, true); mfi.isValid(); ++mfi)
				{
					Box vbx = mfi.validbox();
					AMREX_D_TERM(vbx.growLo(0,1);, vbx.growLo(1,1);, vbx.growLo(2,1));
					Box bx = mfi.growntilebox(1);
					bx &= vbx;
					amrex_mlndlap_stencil_rap(BL_TO_FORTRAN_BOX(bx),
								  BL_TO_FORTRAN_ANYD((*pcrse)[mfi]),
								  BL_TO_FORTRAN_ANYD(fine[mfi]));
				}

				for (MFIter mfi(*pcrse,true); mfi.isValid(); ++mfi)
				{
					const Box& bx = mfi.tilebox();
					FArrayBox& stfab = (*pcrse)[mfi];
                    
					amrex_mlndlap_set_stencil_s0(BL_TO_FORTRAN_BOX(bx),
								     BL_TO_FORTRAN_ANYD(stfab));
				}
			}

			if (need_parallel_copy) {
				crse.ParallelCopy(cfine);
			}

			m_stencil[amrlev][mglev]->FillBoundary(m_geom[amrlev][mglev].periodicity());
		}
	}
}

void
Operator::fixUpResidualMask (int amrlev, iMultiFab& resmsk)
{
	TRACER;
	if (!m_masks_built) buildMasks();

	const iMultiFab& cfmask = *m_nd_fine_mask[amrlev];

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(resmsk,true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		amrex_mlndlap_fixup_res_mask(BL_TO_FORTRAN_BOX(bx),
					     BL_TO_FORTRAN_ANYD(resmsk[mfi]),
					     BL_TO_FORTRAN_ANYD(cfmask[mfi]));
	}
}

void
Operator::prepareForSolve ()
{
	TRACER;
	BL_PROFILE("Operator::prepareForSolve()");

	AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_num_amr_levels == 1 ||
					 m_coarsening_strategy != CoarseningStrategy::RAP,
					 "Operator::prepareForSolve RAP TODO");

	MLNodeLinOp::prepareForSolve();

	buildMasks();

	averageDownCoeffs();

#if (AMREX_SPACEDIM == 2)
	amrex_mlndlap_set_rz(&m_is_rz);
#endif

	buildStencil();
}

void
Operator::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
	TRACER;
	Util::Abort("multigrid not yet allowed");
}

void
Operator::interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const
{
	TRACER;
	Util::Abort("multigrid not yet supported");
}

void
Operator::averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
				  const MultiFab& fine_sol, const MultiFab& fine_rhs)
{
	TRACER;
	Util::Abort("LINE 836");
}

void
Operator::applyBC (int amrlev, int mglev, MultiFab& phi, BCMode/* bc_mode*/,
		   amrex::MLLinOp::StateMode /**/, bool skip_fillboundary) const
{
	TRACER;
	BL_PROFILE("Operator::applyBC()");

	const Geometry& geom = m_geom[amrlev][mglev];
	const Box& nd_domain = amrex::surroundingNodes(geom.Domain());

	if (!skip_fillboundary) {
		phi.FillBoundary(geom.periodicity());
	}

	if (m_coarsening_strategy == CoarseningStrategy::Sigma)
	{
#ifdef _OPENMP
#pragma omp parallel
#endif
		for (MFIter mfi(phi); mfi.isValid(); ++mfi)
		{
			if (!nd_domain.strictly_contains(mfi.fabbox()))
			{
				amrex_mlndlap_applybc(BL_TO_FORTRAN_ANYD(phi[mfi]),
						      BL_TO_FORTRAN_BOX(nd_domain),
						      m_lobc.data(), m_hibc.data());
			}
		}
	}
}

// void
// Operator::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
// {
//     BL_PROFILE("Operator::Fapply()");

//     const auto& sigma = m_sigma[amrlev][mglev];
//     const auto& stencil = m_stencil[amrlev][mglev];
//     const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

//     const Box& domain_box = amrex::surroundingNodes(m_geom[amrlev][mglev].Domain());
//     const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][mglev];

// #ifdef _OPENMP
// #pragma omp parallel
// #endif
//     for (MFIter mfi(out,true); mfi.isValid(); ++mfi)
//     {
//         const Box& bx = mfi.tilebox();
//         const FArrayBox& xfab = in[mfi];
//         FArrayBox& yfab = out[mfi];

//         if (m_coarsening_strategy == CoarseningStrategy::RAP)
//         {
//             amrex_mlndlap_adotx_sten(BL_TO_FORTRAN_BOX(bx),
//                                      BL_TO_FORTRAN_ANYD(yfab),
//                                      BL_TO_FORTRAN_ANYD(xfab),
//                                      BL_TO_FORTRAN_ANYD((*stencil)[mfi]),
//                                      BL_TO_FORTRAN_ANYD(dmsk[mfi]));
//         }
//         else if (m_use_harmonic_average && mglev > 0)
//         {
//             AMREX_D_TERM(const FArrayBox& sxfab = (*sigma[0])[mfi];,
//                          const FArrayBox& syfab = (*sigma[1])[mfi];,
//                          const FArrayBox& szfab = (*sigma[2])[mfi];);

//             amrex_mlndlap_adotx_ha(BL_TO_FORTRAN_BOX(bx),
//                                    BL_TO_FORTRAN_ANYD(yfab),
//                                    BL_TO_FORTRAN_ANYD(xfab),
//                                    AMREX_D_DECL(BL_TO_FORTRAN_ANYD(sxfab),
//                                                 BL_TO_FORTRAN_ANYD(syfab),
//                                                 BL_TO_FORTRAN_ANYD(szfab)),
//                                    BL_TO_FORTRAN_ANYD(dmsk[mfi]),
//                                    dxinv, BL_TO_FORTRAN_BOX(domain_box),
//                                    m_lobc.data(), m_hibc.data());
//         }
//         else
//         {
//             const FArrayBox& sfab = (*sigma[0])[mfi];

//             amrex_mlndlap_adotx_aa(BL_TO_FORTRAN_BOX(bx),
//                                    BL_TO_FORTRAN_ANYD(yfab),
//                                    BL_TO_FORTRAN_ANYD(xfab),
//                                    BL_TO_FORTRAN_ANYD(sfab),
//                                    BL_TO_FORTRAN_ANYD(dmsk[mfi]),
//                                    dxinv, BL_TO_FORTRAN_BOX(domain_box),
//                                    m_lobc.data(), m_hibc.data());
//         }
//     }
// }

// void
// Operator::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const
// {
//     BL_PROFILE("Operator::Fsmooth()");

//     const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][mglev];

//     if (m_use_gauss_seidel)
//     {
//         const auto& sigma = m_sigma[amrlev][mglev];
//         const auto& stencil = m_stencil[amrlev][mglev];
//         const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

//         const Box& domain_box = amrex::surroundingNodes(m_geom[amrlev][mglev].Domain());

//         if (m_coarsening_strategy == CoarseningStrategy::RAP)
//         {
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
//             for (MFIter mfi(sol); mfi.isValid(); ++mfi)
//             {
//                 const Box& bx = mfi.validbox();
//                 amrex_mlndlap_gauss_seidel_sten(BL_TO_FORTRAN_BOX(bx),
//                                                 BL_TO_FORTRAN_ANYD(sol[mfi]),
//                                                 BL_TO_FORTRAN_ANYD(rhs[mfi]),
//                                                 BL_TO_FORTRAN_ANYD((*stencil)[mfi]),
//                                                 BL_TO_FORTRAN_ANYD(dmsk[mfi]));
//             }
//         }
//         else if (m_use_harmonic_average && mglev > 0)
//         {
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
//             for (MFIter mfi(sol); mfi.isValid(); ++mfi)
//             {
//                 const Box& bx = mfi.validbox();
//                 AMREX_D_TERM(const FArrayBox& sxfab = (*sigma[0])[mfi];,
//                              const FArrayBox& syfab = (*sigma[1])[mfi];,
//                              const FArrayBox& szfab = (*sigma[2])[mfi];);

//                 amrex_mlndlap_gauss_seidel_ha(BL_TO_FORTRAN_BOX(bx),
//                                               BL_TO_FORTRAN_ANYD(sol[mfi]),
//                                               BL_TO_FORTRAN_ANYD(rhs[mfi]),
//                                               AMREX_D_DECL(BL_TO_FORTRAN_ANYD(sxfab),
//                                                            BL_TO_FORTRAN_ANYD(syfab),
//                                                            BL_TO_FORTRAN_ANYD(szfab)),
//                                               BL_TO_FORTRAN_ANYD(dmsk[mfi]),
//                                               dxinv, BL_TO_FORTRAN_BOX(domain_box),
//                                               m_lobc.data(), m_hibc.data());
//             }
//         }
//         else
//         {
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
//             for (MFIter mfi(sol); mfi.isValid(); ++mfi)
//             {
//                 const Box& bx = mfi.validbox();
//                 const FArrayBox& sfab = (*sigma[0])[mfi];

//                 amrex_mlndlap_gauss_seidel_aa(BL_TO_FORTRAN_BOX(bx),
//                                               BL_TO_FORTRAN_ANYD(sol[mfi]),
//                                               BL_TO_FORTRAN_ANYD(rhs[mfi]),
//                                               BL_TO_FORTRAN_ANYD(sfab),
//                                               BL_TO_FORTRAN_ANYD(dmsk[mfi]),
//                                               dxinv, BL_TO_FORTRAN_BOX(domain_box),
//                                               m_lobc.data(), m_hibc.data());
//             }
//         }

//         nodalSync(amrlev, mglev, sol);
//     }
//     else
//     {
//         MultiFab Ax(sol.boxArray(), sol.DistributionMap(), 1, 0);
//         Fapply(amrlev, mglev, Ax, sol);

//         const auto& sigma = m_sigma[amrlev][mglev];
//         const auto& stencil = m_stencil[amrlev][mglev];
//         const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

//         const Box& domain_box = amrex::surroundingNodes(m_geom[amrlev][mglev].Domain());

//         if (m_coarsening_strategy == CoarseningStrategy::RAP)
//         {
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
//             for (MFIter mfi(sol); mfi.isValid(); ++mfi)
//             {
//                 const Box& bx = mfi.validbox();
//                 amrex_mlndlap_jacobi_sten(BL_TO_FORTRAN_BOX(bx),
//                                           BL_TO_FORTRAN_ANYD(sol[mfi]),
//                                           BL_TO_FORTRAN_ANYD(Ax[mfi]),
//                                           BL_TO_FORTRAN_ANYD(rhs[mfi]),
//                                           BL_TO_FORTRAN_ANYD((*stencil)[mfi]),
//                                           BL_TO_FORTRAN_ANYD(dmsk[mfi]));
//             }
//         }
//         else if (m_use_harmonic_average && mglev > 0)
//         {
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
//             for (MFIter mfi(sol,true); mfi.isValid(); ++mfi)
//             {
//                 const Box& bx = mfi.tilebox();
//                 AMREX_D_TERM(const FArrayBox& sxfab = (*sigma[0])[mfi];,
//                              const FArrayBox& syfab = (*sigma[1])[mfi];,
//                              const FArrayBox& szfab = (*sigma[2])[mfi];);

//                 amrex_mlndlap_jacobi_ha(BL_TO_FORTRAN_BOX(bx),
//                                         BL_TO_FORTRAN_ANYD(sol[mfi]),
//                                         BL_TO_FORTRAN_ANYD(Ax[mfi]),
//                                         BL_TO_FORTRAN_ANYD(rhs[mfi]),
//                                         AMREX_D_DECL(BL_TO_FORTRAN_ANYD(sxfab),
//                                                      BL_TO_FORTRAN_ANYD(syfab),
//                                                      BL_TO_FORTRAN_ANYD(szfab)),
//                                         BL_TO_FORTRAN_ANYD(dmsk[mfi]),
//                                         dxinv, BL_TO_FORTRAN_BOX(domain_box),
//                                         m_lobc.data(), m_hibc.data());
//             }
//         }
//         else
//         {
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
//             for (MFIter mfi(sol,true); mfi.isValid(); ++mfi)
//             {
//                 const Box& bx = mfi.tilebox();
//                 const FArrayBox& sfab = (*sigma[0])[mfi];

//                 amrex_mlndlap_jacobi_aa(BL_TO_FORTRAN_BOX(bx),
//                                         BL_TO_FORTRAN_ANYD(sol[mfi]),
//                                         BL_TO_FORTRAN_ANYD(Ax[mfi]),
//                                         BL_TO_FORTRAN_ANYD(rhs[mfi]),
//                                         BL_TO_FORTRAN_ANYD(sfab),
//                                         BL_TO_FORTRAN_ANYD(dmsk[mfi]),
//                                         dxinv, BL_TO_FORTRAN_BOX(domain_box),
//                                         m_lobc.data(), m_hibc.data());
//             }
//         }
//     }
// }


void
Operator::compSyncResidualCoarse (MultiFab& sync_resid, const MultiFab& a_phi,
				  const MultiFab& vold, const MultiFab* rhcc,
				  const BoxArray& fine_grids, const IntVect& ref_ratio)
{
	TRACER;
	Util::Abort("LINE 1091");
}

void
Operator::compSyncResidualFine (MultiFab& sync_resid, const MultiFab& phi, const MultiFab& vold,
				const MultiFab* rhcc)
{
	TRACER;
	Util::Abort("LINE 1098");
}

void
Operator::reflux (int crse_amrlev,
		  MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
		  MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const
{
	TRACER;
	Util::Abort("LINE 1436");
}




// Operator::Operator ()
// {
// 	m_ixtype = amrex::IntVect::TheNodeVector();
// }

// Operator::~Operator () {}

// void
// Operator::define (amrex::Vector<amrex::Geometry> a_geom,
// 		   const amrex::Vector<amrex::BoxArray>& a_grids,
// 		   const amrex::Vector<amrex::DistributionMapping>& a_dmap,
// 		   BC::BC& a_bc,
// 		   const amrex::LPInfo& a_info,
// 		   const amrex::Vector<amrex::FabFactory<amrex::FArrayBox> const*>& a_factory)
// {
// 	m_bc = &a_bc;

// 	std::array<int,AMREX_SPACEDIM> is_periodic = m_bc->IsPeriodic();
// 	// for (int ilev=0; ilev < a_geom.size(); ilev++)
// 	// 	a_geom[ilev].SetPeriodicity(is_periodic);
	
// 	MLLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory);
// 	defineAuxData();
// 	defineBC();


// 	//setDomainBC

// 	m_lobc = {AMREX_D_DECL(is_periodic[0] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet,
// 			       is_periodic[1] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet,
// 			       is_periodic[2] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet)};
// 	m_hibc = {AMREX_D_DECL(is_periodic[0] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet,
// 			       is_periodic[1] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet,
// 			       is_periodic[2] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet)};
				  
// 	for (int ilev = 0; ilev < a_geom.size(); ++ilev)
// 		setLevelBC(ilev,nullptr);

// }

// void
// Operator::defineAuxData ()
// {
// 	BL_PROFILE("Operator::defineAuxData()");

// 	m_undrrelxr.resize(m_num_amr_levels);
// 	m_maskvals.resize(m_num_amr_levels);
// 	//m_fluxreg.resize(m_num_amr_levels-1);

// 	const int ncomp = getNComp();

// 	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
// 	{
// 		m_undrrelxr[amrlev].resize(m_num_mg_levels[amrlev]);
// 		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
// 		{
// 			m_undrrelxr[amrlev][mglev].define(m_grids[amrlev][mglev],
// 							  m_dmap[amrlev][mglev],
// 							  1, 0, 0, ncomp);
// 		}
// 	}
    
// 	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
// 	{
// 		m_maskvals[amrlev].resize(m_num_mg_levels[amrlev]);
// 		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
// 		{
// 			for (amrex::OrientationIter oitr; oitr; ++oitr)
// 			{
// 				const amrex::Orientation face = oitr();
// 				const int ngrow = 1;
// 				const int extent = 1;
// 				m_maskvals[amrlev][mglev][face].define(m_grids[amrlev][mglev],
// 								       m_dmap[amrlev][mglev],
// 								       m_geom[amrlev][mglev],
// 								       face, 0, ngrow, extent, ncomp, true);
// 			}
// 		}
// 	}
// }

// void
// Operator::defineBC ()
// {
// 	BL_PROFILE("Operator::defineBC()");

// 	const int ncomp = getNComp();

// 	m_bndry_sol.resize(m_num_amr_levels);
// 	m_crse_sol_br.resize(m_num_amr_levels);

// 	m_bndry_cor.resize(m_num_amr_levels);
// 	m_crse_cor_br.resize(m_num_amr_levels);

// 	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
// 	{
// 		m_bndry_sol[amrlev].reset(new amrex::MLMGBndry(m_grids[amrlev][0], m_dmap[amrlev][0],
// 							       ncomp, m_geom[amrlev][0]));
// 	}

// 	for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
// 	{
// 		const int in_rad = 0;
// 		const int out_rad = 1;
// 		const int extent_rad = 2;
// 		const int crse_ratio = m_amr_ref_ratio[amrlev-1];
// 		amrex::BoxArray cba = m_grids[amrlev][0];
// 		cba.coarsen(crse_ratio);
// 		m_crse_sol_br[amrlev].reset(new amrex::BndryRegister(cba, m_dmap[amrlev][0],
// 								     in_rad, out_rad, extent_rad, ncomp));
// 	}

// 	for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
// 	{
// 		const int in_rad = 0;
// 		const int out_rad = 1;
// 		const int extent_rad = 2;
// 		const int crse_ratio = m_amr_ref_ratio[amrlev-1];
// 		amrex::BoxArray cba = m_grids[amrlev][0];
// 		cba.coarsen(crse_ratio);
// 		m_crse_cor_br[amrlev].reset(new amrex::BndryRegister(cba, m_dmap[amrlev][0],
// 								     in_rad, out_rad, extent_rad, ncomp));
// 		m_crse_cor_br[amrlev]->setVal(0.0);
// 	}

// 	// This has be to done after m_crse_cor_br is defined.
// 	for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
// 	{
// 		m_bndry_cor[amrlev].reset(new amrex::MLMGBndry(m_grids[amrlev][0], m_dmap[amrlev][0],
// 							       ncomp, m_geom[amrlev][0]));
// 		amrex::MultiFab bc_data(m_grids[amrlev][0], m_dmap[amrlev][0], ncomp, 1);
// 		bc_data.setVal(0.0);

// 		m_bndry_cor[amrlev]->setBndryValues(*m_crse_cor_br[amrlev], 0, bc_data, 0, 0, ncomp,
// 						    m_amr_ref_ratio[amrlev-1], amrex::BCRec());

// 		m_bndry_cor[amrlev]->setLOBndryConds({AMREX_D_DECL(BCType::Dirichlet,
// 								   BCType::Dirichlet,
// 								   BCType::Dirichlet)},
// 			{AMREX_D_DECL(BCType::Dirichlet,
// 				      BCType::Dirichlet,
// 				      BCType::Dirichlet)},
// 			m_amr_ref_ratio[amrlev-1], amrex::RealVect{});
// 	}

// 	m_bcondloc.resize(m_num_amr_levels);
// 	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
// 	{
// 		m_bcondloc[amrlev].resize(m_num_mg_levels[amrlev]);
// 		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
// 		{
// 			m_bcondloc[amrlev][mglev].reset(new BndryCondLoc(m_grids[amrlev][mglev],
// 										m_dmap[amrlev][mglev]));
// 		} 
// 	}
// }

// // void
// // Operator::setLevelBC (int amrlev, const amrex::MultiFab* a_levelbcdata)
// // {
// // 	BL_PROFILE("Operator::setLevelBC()");

// // 	AMREX_ALWAYS_ASSERT(amrlev >= 0 && amrlev < m_num_amr_levels);

// // 	const int ncomp = getNComp();

// // 	amrex::MultiFab zero;
// // 	if (a_levelbcdata == nullptr) {
// // 		zero.define(m_grids[amrlev][0], m_dmap[amrlev][0], ncomp, 1);
// // 		zero.setVal(0.0);
// // 	} else {
// // 		AMREX_ALWAYS_ASSERT(a_levelbcdata->nGrow() >= 1);
// // 	}
// // 	const amrex::MultiFab& bcdata = (a_levelbcdata == nullptr) ? zero : *a_levelbcdata;

// // 	int br_ref_ratio = -1;

// // 	if (amrlev == 0)
// // 	{
// // 		if (needsCoarseDataForBC())
// // 		{
// // 			br_ref_ratio = m_coarse_data_crse_ratio > 0 ? m_coarse_data_crse_ratio : 2;
// // 			if (m_crse_sol_br[amrlev] == nullptr && br_ref_ratio > 0)
// // 			{
// // 				const int in_rad = 0;
// // 				const int out_rad = 1;
// // 				const int extent_rad = 2;
// // 				const int crse_ratio = br_ref_ratio;
// // 				amrex::BoxArray cba = m_grids[amrlev][0];
// // 				cba.coarsen(crse_ratio);
// // 				m_crse_sol_br[amrlev].reset(new amrex::BndryRegister(cba, m_dmap[amrlev][0],
// // 										     in_rad, out_rad,
// // 										     extent_rad, ncomp));
// // 			}
// // 			if (m_coarse_data_for_bc != nullptr) {
// // 				AMREX_ALWAYS_ASSERT(m_coarse_data_crse_ratio > 0);
// // 				const amrex::Box& cbx = amrex::coarsen(m_geom[0][0].Domain(), m_coarse_data_crse_ratio);
// // 				m_crse_sol_br[amrlev]->copyFrom(*m_coarse_data_for_bc, 0, 0, 0, ncomp,
// // 								m_bc->Periodicity(cbx));
// // 			} else {
// // 				m_crse_sol_br[amrlev]->setVal(0.0);
// // 			}
// // 			m_bndry_sol[amrlev]->setBndryValues(*m_crse_sol_br[amrlev], 0,
// // 							    bcdata, 0, 0, ncomp,
// // 							    br_ref_ratio, amrex::BCRec());
// // 			br_ref_ratio = m_coarse_data_crse_ratio;
// // 		}
// // 		else
// // 		{
// // 			m_bndry_sol[amrlev]->setBndryValues(bcdata,0,0,ncomp,amrex::BCRec());
// // 			br_ref_ratio = 1;
// // 		}
// // 	}
// // 	else
// // 	{
// // 		m_bndry_sol[amrlev]->setBndryValues(bcdata,0,0,ncomp, m_amr_ref_ratio[amrlev-1], amrex::BCRec());
// // 		br_ref_ratio = m_amr_ref_ratio[amrlev-1];
// // 	}

// // 	m_bndry_sol[amrlev]->setLOBndryConds(m_lobc, m_hibc, br_ref_ratio, m_coarse_bc_loc);

// // 	const amrex::Real* dx = m_geom[amrlev][0].CellSize();
// // 	for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
// // 	{
// // 		m_bcondloc[amrlev][mglev]->setLOBndryConds(m_geom[amrlev][mglev], dx,
// // 							   m_lobc, m_hibc,
// // 							   br_ref_ratio, m_coarse_bc_loc);
// // 	}
// // }

// amrex::BoxArray
// Operator::makeNGrids (int grid_size) const
// {
// 	const amrex::Box& dombx = m_geom[0].back().Domain();

// 	const amrex::BoxArray& old_ba = m_grids[0].back();
// 	const int N = old_ba.size();
// 	amrex::Vector<amrex::Box> bv;
// 	bv.reserve(N);
// 	for (int i = 0; i < N; ++i)
// 	{
// 		amrex::Box b = old_ba[i];
// 		b.coarsen(grid_size);
// 		b.refine(grid_size);
// 		amrex::IntVect sz = b.size();
// 		const amrex::IntVect nblks {AMREX_D_DECL(sz[0]/grid_size, sz[1]/grid_size, sz[2]/grid_size)};
        
// 		amrex::IntVect big = b.smallEnd() + grid_size - 1;
// 		b.setBig(big);

// #if (AMREX_SPACEDIM == 3)
// 		for (int kk = 0; kk < nblks[2]; ++kk) {
// #endif
// #if (AMREX_SPACEDIM >= 2)
// 			for (int jj = 0; jj < nblks[1]; ++jj) {
// #endif
// 				for (int ii = 0; ii < nblks[0]; ++ii)
// 				{
// 					amrex::IntVect shft{AMREX_D_DECL(ii*grid_size,jj*grid_size,kk*grid_size)};
// 					amrex::Box bb = amrex::shift(b,shft);
// 					bb &= dombx;
// 					bv.push_back(bb);
// 				}
// #if (AMREX_SPACEDIM >= 2)
// 			}
// #endif
// #if (AMREX_SPACEDIM == 3)
// 		}
// #endif
// 	}

// 	std::sort(bv.begin(), bv.end());
// 	bv.erase(std::unique(bv.begin(), bv.end()), bv.end());

// 	amrex::BoxList bl(std::move(bv));

// 	return amrex::BoxArray{std::move(bl)};
// }

// void
// Operator::restriction (int, int, amrex::MultiFab& crse, amrex::MultiFab& fine) const
// {
// 	const int ncomp = getNComp();
// 	amrex::average_down(fine, crse, 0, ncomp, 2);
// }

// void
// Operator::interpolation (int /*amrlev*/, int /*fmglev*/,
// 			  amrex::MultiFab& fine, const amrex::MultiFab& crse) const
// {
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
// 	for (amrex::MFIter mfi(crse,true); mfi.isValid(); ++mfi)
// 	{
// 		const amrex::Box&         bx    = mfi.tilebox();
// 		const int          ncomp = getNComp();
// 		const amrex::FArrayBox& cfab    = crse[mfi];
// 		amrex::FArrayBox&       ffab    = fine[mfi];

// 		amrex_mg_interp(ffab.dataPtr(),
// 				AMREX_ARLIM(ffab.loVect()), AMREX_ARLIM(ffab.hiVect()),
// 				cfab.dataPtr(),
// 				AMREX_ARLIM(cfab.loVect()), AMREX_ARLIM(cfab.hiVect()),
// 				bx.loVect(), bx.hiVect(), &ncomp);
// 	}    
// }

// void
// Operator::averageDownSolutionRHS (int camrlev, amrex::MultiFab& crse_sol, amrex::MultiFab& crse_rhs,
// 				   const amrex::MultiFab& fine_sol, const amrex::MultiFab& fine_rhs)
// {
// 	const auto amrrr = AMRRefRatio(camrlev);
// 	const int ncomp = getNComp();
// 	amrex::average_down(fine_sol, crse_sol, 0, ncomp, amrrr);
// 	amrex::average_down(fine_rhs, crse_rhs, 0, ncomp, amrrr);
// }

// // void
// // Operator::apply (int amrlev, int mglev, amrex::MultiFab& out, amrex::MultiFab& in, BCMode bc_mode,
// // 		  const amrex::MLMGBndry* bndry) const
// // {
// // 	BL_PROFILE("Operator::apply()");
// // 	applyBC(amrlev, mglev, in, bc_mode, bndry);
// // 	Fapply(amrlev, mglev, out, in);
// // }

// // void
// // Operator::smooth (int amrlev, int mglev, amrex::MultiFab& sol, const amrex::MultiFab& rhs,
// // 		   bool skip_fillboundary) const
// // {
// // 	BL_PROFILE("Operator::smooth()");
// // 	for (int redblack = 0; redblack < 2; ++redblack)
// // 	{
// // 		applyBC(amrlev, mglev, sol, BCMode::Homogeneous, nullptr, skip_fillboundary);
// // 		Fsmooth(amrlev, mglev, sol, rhs, redblack);
// // 		skip_fillboundary = false;
// // 	}
// // }

// // void
// // Operator::updateSolBC (int amrlev, const amrex::MultiFab& crse_bcdata) const
// // {
// // 	BL_PROFILE("Operator::updateSolBC()");

// // 	AMREX_ALWAYS_ASSERT(amrlev > 0);
// // 	const int ncomp = getNComp();
// // 	m_bc->define(m_geom[amrlev-1][0]);
// // 	m_crse_sol_br[amrlev]->copyFrom(crse_bcdata, 0, 0, 0, ncomp, m_bc->Periodicity());
// // 	m_bndry_sol[amrlev]->updateBndryValues(*m_crse_sol_br[amrlev], 0, 0, ncomp, m_amr_ref_ratio[amrlev-1]);
// // }

// // void
// // Operator::updateCorBC (int amrlev, const amrex::MultiFab& crse_bcdata) const
// // {
// // 	BL_PROFILE("Operator::updateCorBC()");
// // 	AMREX_ALWAYS_ASSERT(amrlev > 0);
// // 	const int ncomp = getNComp();
// // 	m_bc->define(m_geom[amrlev-1][0]);
// // 	m_crse_cor_br[amrlev]->copyFrom(crse_bcdata, 0, 0, 0, ncomp, m_bc->Periodicity());
// // 	m_bndry_cor[amrlev]->updateBndryValues(*m_crse_cor_br[amrlev], 0, 0, ncomp, m_amr_ref_ratio[amrlev-1]);
// // }

// // void
// // Operator::solutionResidual (int amrlev, amrex::MultiFab& resid, amrex::MultiFab& x,
// // 			     const amrex::MultiFab& b, const amrex::MultiFab* crse_bcdata)
// // {
// // 	BL_PROFILE("Operator::solutionResidual()");
// // 	const int ncomp = getNComp();
// // 	const int mglev = 0;
// // 	if (crse_bcdata != nullptr) {
// // 		updateSolBC(amrlev, *crse_bcdata);
// // 	}
// // 	apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, m_bndry_sol[amrlev].get());

// // 	AMREX_ALWAYS_ASSERT(resid.nComp() == b.nComp());
// // 	amrex::MultiFab::Xpay(resid, -1.0, b, 0, 0, ncomp, 0);
// // }

// // void
// // Operator::fillSolutionBC (int amrlev, amrex::MultiFab& sol, const amrex::MultiFab* crse_bcdata)
// // {
// // 	BL_PROFILE("Operator::fillSolutionBC()");
// // 	const int mglev = 0;

// // 	if (crse_bcdata != nullptr) {
// // 		updateSolBC(amrlev, *crse_bcdata);
// // 	}
// // 	applyBC(amrlev, mglev, sol, BCMode::Inhomogeneous, m_bndry_sol[amrlev].get());    
// // }

// // void
// // Operator::correctionResidual (int amrlev, int mglev, amrex::MultiFab& resid, amrex::MultiFab& x,
// // 			       const amrex::MultiFab& b,
// // 			       BCMode bc_mode, const amrex::MultiFab* crse_bcdata)
// // {
// // 	BL_PROFILE("Operator::correctionResidual()");
// // 	const int ncomp = getNComp();
// // 	if (bc_mode == BCMode::Inhomogeneous)
// // 	{
// // 		if (crse_bcdata)
// // 		{
// // 			AMREX_ALWAYS_ASSERT(mglev == 0);
// // 			AMREX_ALWAYS_ASSERT(amrlev > 0);
// // 			updateCorBC(amrlev, *crse_bcdata);
// // 		}
// // 		apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, m_bndry_cor[amrlev].get());
// // 	}
// // 	else
// // 	{
// // 		AMREX_ALWAYS_ASSERT(crse_bcdata == nullptr);
// // 		apply(amrlev, mglev, resid, x, BCMode::Homogeneous, nullptr);
// // 	}

// // 	amrex::MultiFab::Xpay(resid, -1.0, b, 0, 0, ncomp, 0);
// // }

// // void
// // MLNodeLaplacian::applyBC (int amrlev, int mglev, MultiFab& phi, BCMode/* bc_mode*/,
// //                           bool skip_fillboundary) const

// void
// Operator::applyBC (int amrlev, int mglev, amrex::MultiFab& in, amrex::MLLinOp::BCMode bc_mode, bool skip_fillboundary) const
// {
//     BL_PROFILE("MLNodeLaplacian::applyBC()");

//     const amrex::Geometry& geom = m_geom[amrlev][mglev];
//     const amrex::Box& nd_domain = amrex::surroundingNodes(geom.Domain());

//     if (!skip_fillboundary) {
//         in.FillBoundary(geom.periodicity());
//     }

// // //    int inhom = (bc_mode == BCMode::Inhomogeneous);

// //     if (m_coarsening_strategy == amrex::CoarseningStrategy::Sigma)
// //     {
// // #ifdef _OPENMP
// // #pragma omp parallel
// // #endif
// //         for (MFIter mfi(in); mfi.isValid(); ++mfi)
// //         {
// //             if (!nd_domain.strictly_contains(mfi.fabbox()))
// //             {
// //                 amrex_mlndlap_applybc(BL_TO_FORTRAN_ANYD(in[mfi]),
// //                                       BL_TO_FORTRAN_BOX(nd_domain),
// //                                       m_lobc.data(), m_hibc.data());
// //             }
// //         }
// //     }
// }

// void
// Operator::reflux (int,
// 		   amrex::MultiFab&, const amrex::MultiFab&, const amrex::MultiFab&,
// 		   amrex::MultiFab&, amrex::MultiFab&, const amrex::MultiFab&) const
// {
// }

// // /// \note This function was stripped and replaced with setVal.
// // void
// // Operator::compFlux (int /*amrlev*/, const std::array<amrex::MultiFab*,AMREX_SPACEDIM>& fluxes,
// // 		     amrex::MultiFab& /*sol*/) const
// // {
// // 	BL_PROFILE("Operator::compFlux()");
// // 	for (int idim=0; idim < AMREX_SPACEDIM; idim++)
// // 		fluxes[idim]->setVal(0.0);
// // 	return;
// // }

// // void
// // Operator::compGrad (int amrlev, const std::array<amrex::MultiFab*,AMREX_SPACEDIM>& grad,
// // 		    amrex::MultiFab& sol) const
// // {
// // 	BL_PROFILE("Operator::compGrad()");

// // 	if (sol.nComp() > 1)
// // 		amrex::Abort("Operator::compGrad called, but only works for single-component solves");

// // 	const int mglev = 0;

// // 	// m_bc->define(m_geom[amrlev][mglev]);
// // 	// m_bc->FillBoundary(sol,0,0,0.0);
// // 	// auto *nonconst_this = const_cast<Operator*>(this);
// // 	// nonconst_this->setLevelBC(amrlev,&sol);

// // 	applyBC(amrlev, mglev, sol, BCMode::Inhomogeneous, m_bndry_sol[amrlev].get());

// // 	const amrex::Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

// // #ifdef _OPENMP
// // #pragma omp parallel
// // #endif
// // 	for (amrex::MFIter mfi(sol, amrex::MFItInfo().EnableTiling().SetDynamic(true));  mfi.isValid(); ++mfi)
// // 	{
// // 		AMREX_D_TERM(const amrex::Box& xbx = mfi.nodaltilebox(0);,
// // 			     const amrex::Box& ybx = mfi.nodaltilebox(1);,
// // 			     const amrex::Box& zbx = mfi.nodaltilebox(2););
// // 		amrex_mllinop_grad(AMREX_D_DECL(BL_TO_FORTRAN_BOX(xbx),
// // 						BL_TO_FORTRAN_BOX(ybx),
// // 						BL_TO_FORTRAN_BOX(zbx)),
// // 				   BL_TO_FORTRAN_ANYD(sol[mfi]),
// // 				   AMREX_D_DECL(BL_TO_FORTRAN_ANYD((*grad[0])[mfi]),
// // 						BL_TO_FORTRAN_ANYD((*grad[1])[mfi]),
// // 						BL_TO_FORTRAN_ANYD((*grad[2])[mfi])),
// // 				   dxinv);
// // 	}
// // }

// void
// Operator::prepareForSolve ()
// {
// 	BL_PROFILE("Operator::prepareForSolve()");

// 	const int ncomp = getNComp();
// 	for (int amrlev = 0;  amrlev < m_num_amr_levels; ++amrlev)
// 	{
// 		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
// 		{
// 			const auto& bcondloc = *m_bcondloc[amrlev][mglev];
// 			const auto& maskvals = m_maskvals[amrlev][mglev];
// 			const amrex::Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

// 			amrex::BndryRegister& undrrelxr = m_undrrelxr[amrlev][mglev];
// 			amrex::MultiFab foo(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], ncomp, 0, amrex::MFInfo().SetAlloc(false));
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
// 			for (amrex::MFIter mfi(foo, amrex::MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
// 			{
// 				const amrex::Box& vbx = mfi.validbox();

// 				const RealTuple & bdl = bcondloc.bndryLocs(mfi);
// 				const BCTuple   & bdc = bcondloc.bndryConds(mfi);

// 				for (amrex::OrientationIter oitr; oitr; ++oitr)
// 				{
// 					const amrex::Orientation ori = oitr();
                    
// 					int  cdr = ori;
// 					amrex::Real bcl = bdl[ori];
// 					int  bct = bdc[ori];
                    
// 					amrex::FArrayBox& ffab = undrrelxr[ori][mfi];
// 					const amrex::Mask& m   =  maskvals[ori][mfi];

// 					amrex_mllinop_comp_interp_coef0(BL_TO_FORTRAN_BOX(vbx),
// 									BL_TO_FORTRAN_ANYD(ffab),
// 									BL_TO_FORTRAN_ANYD(m),
// 									cdr, bct, bcl, maxorder, dxinv, ncomp);
// 				}
// 			}
// 		}
// 	}
// 	averageDownCoeffs();
// }

// // amrex::Real
// // Operator::xdoty (int amrlev, int mglev,
// // 		  const amrex::MultiFab& x, const amrex::MultiFab& y, bool local) const
// // {
// // 	const int ncomp = getNComp();
// // 	const int nghost = 0;
// // 	amrex::Real result = amrex::MultiFab::Dot(x,0,y,0,ncomp,nghost,true);
// // 	if (!local) {
// // 		amrex::ParallelAllReduce::Sum(result, Communicator(amrlev, mglev));
// // 	}
// // 	return result;
// // }

// Operator::BndryCondLoc::BndryCondLoc (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
// 	: bcond(ba, dm),
// 	  bcloc(ba, dm)
// {
// }

// void
// Operator::BndryCondLoc::setLOBndryConds (const amrex::Geometry& geom, const amrex::Real* dx,
// 					  const std::array<BCType,AMREX_SPACEDIM>& lobc,
// 					  const std::array<BCType,AMREX_SPACEDIM>& hibc,
// 					  int ratio, const amrex::RealVect& a_loc)
// {
// 	const amrex::Box&  domain = geom.Domain();

// #ifdef _OPENMP
// #pragma omp parallel
// #endif
// 	for (amrex::MFIter mfi(bcloc); mfi.isValid(); ++mfi)
// 	{
// 		const amrex::Box& bx = mfi.validbox();
// 		RealTuple & bloc  = bcloc[mfi];
// 		BCTuple   & bctag = bcond[mfi];

// 		amrex::MLMGBndry::setBoxBC(bloc, bctag, bx, domain, lobc, hibc, dx, ratio, a_loc);
// 	}
// }

// // void
// // Operator::applyMetricTerm (int, int, amrex::MultiFab&) const
// // {
// // 	// This is a required method needed only if the geometry is 
// // 	// non-Cartesian. This operator is used for Cartesian coordinates
// // 	// only.
// // }

// // void
// // Operator::unapplyMetricTerm (int, int, amrex::MultiFab&) const
// // {
// // 	// This is a required method needed only if the geometry is 
// // 	// non-Cartesian. This operator is used for Cartesian coordinates
// // 	// only.
// // }
// void
// Operator::averageDownCoeffs ()
// {
// 	for (int i = 0; i < m_num_a_fabs; i++)
// 	{
// 		for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
// 		{
// 			auto& fine_a_coeffs = m_a_coeffs[i][amrlev];
// 			averageDownCoeffsSameAmrLevel(fine_a_coeffs);
// 		}
// 		averageDownCoeffsSameAmrLevel(m_a_coeffs[i][0]);
// 	}
// }

// void
// Operator::averageDownCoeffsSameAmrLevel (amrex::Vector<amrex::MultiFab>& a)
// {
// 	int nmglevs = a.size();
// 	for (int mglev = 1; mglev < nmglevs; ++mglev)
// 	{
// 		amrex::average_down(a[mglev-1], a[mglev], 0, a[0].nComp(), mg_coarsen_ratio);
// 	}
// }


const amrex::FArrayBox &
Operator::GetFab(const int num, const int amrlev, const int mglev, const amrex::MFIter &mfi) const
{
	TRACER;
 	return m_a_coeffs[num][amrlev][mglev][mfi];
}


void
Operator::RegisterNewFab(amrex::Vector<amrex::MultiFab> &input)
{
	TRACER;
	/// \todo assertions here
	m_a_coeffs.resize(m_a_coeffs.size() + 1);
	m_a_coeffs[m_num_a_fabs].resize(m_num_amr_levels);
	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		m_a_coeffs[m_num_a_fabs][amrlev].resize(m_num_mg_levels[amrlev]);
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
			m_a_coeffs[m_num_a_fabs][amrlev][mglev].define(m_grids[amrlev][mglev],
								       m_dmap[amrlev][mglev],
								       input[amrlev].nComp(),
								       input[amrlev].nGrow());

		amrex::MultiFab::Copy(m_a_coeffs[m_num_a_fabs][amrlev][0],
				      input[amrlev], 0, 0,
				      input[amrlev].nComp(),
				      input[amrlev].nGrow());
	}
	m_num_a_fabs++;
}


void
Operator::RegisterNewFab(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &input)
{
	TRACER;
	/// \todo assertions here
	m_a_coeffs.resize(m_a_coeffs.size() + 1);
	m_a_coeffs[m_num_a_fabs].resize(m_num_amr_levels);
	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		m_a_coeffs[m_num_a_fabs][amrlev].resize(m_num_mg_levels[amrlev]);
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
			m_a_coeffs[m_num_a_fabs][amrlev][mglev].define(m_grids[amrlev][mglev],
								       m_dmap[amrlev][mglev],
								       input[amrlev]->nComp(),
								       input[amrlev]->nGrow());

		amrex::MultiFab::Copy(m_a_coeffs[m_num_a_fabs][amrlev][0],
				      *input[amrlev], 0, 0,
				      input[amrlev]->nComp(),
				      input[amrlev]->nGrow());
	}
	m_num_a_fabs++;
}





}

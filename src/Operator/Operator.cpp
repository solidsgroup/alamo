
#include "Operator.H"
#include <AMReX_MLLinOp_F.H>
#include <AMReX_MultiFabUtil.H>
#include "Util/Color.H"


#define TRACER	//std::cout << Color::FG::Yellow << __FILE__ << ":" << __LINE__ << Color::FG::Default << " " << __func__ << std::endl;
#define PROBE	//std::cout << Color::FG::Red << __FILE__ << ":" << __LINE__ << Color::FG::Default << " " << __func__ << std::endl;


using namespace amrex;
namespace Operator {



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

	 Vector<BoxArray> cc_grids = a_grids;
	 for (auto& ba : cc_grids)
		 ba.enclosedCells();
	
	 MLNodeLinOp::define(a_geom, cc_grids, a_dmap, a_info, a_factory);
 }

void
Operator::compRHS (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel,
		   const Vector<const MultiFab*>& rhnd,
		   const Vector<MultiFab*>& a_rhcc) // TODO
{
	TRACER;
	Util::Abort("do not call compRHS");
}

void
Operator::averageDownCoeffs () // TODO
{
	TRACER;
	Util::Abort("Do not call averageDownCoeffs");
}

void
Operator::averageDownCoeffsToCoarseAmrLevel (int flev) // TODO
{
	TRACER;
	Util::Abort("Do not call averageDownCoeffsToCoarseAmrLevel");
}

void
Operator::averageDownCoeffsSameAmrLevel (int amrlev) // TODO
{
	TRACER;
	Util::Abort("Do not call averageDownCoeffsSameAmrLevel");
}

void
Operator::FillBoundaryCoeff (MultiFab& sigma, const Geometry& geom) // TODO
{
	TRACER;
	Util::Abort("Do not call FillBoundaryCoeff");
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

	//averageDownCoeffs(); // TODO

	buildStencil();
}

void
Operator::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const // TODO
{
	TRACER;
	Util::Abort("multigrid not yet implemented");
}

void
Operator::interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const // TODO
{
	TRACER;
	Util::Abort("multigrid not yet implemented");
}

void
Operator::averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
				  const MultiFab& fine_sol, const MultiFab& fine_rhs) 
{
	TRACER;
	const auto& amrrr = AMRRefRatio(camrlev);
	const int ncomp = fine_sol.nComp();
	amrex::average_down(fine_sol, crse_sol, 0, ncomp, amrrr);

	if (isSingular(0))
	{
		Util::Abort("should not be using with singular operator");
		
		MultiFab frhs(fine_rhs.boxArray(), fine_rhs.DistributionMap(), ncomp, 1);
		MultiFab::Copy(frhs, fine_rhs, 0, 0, ncomp, 0);
		restrictInteriorNodes(camrlev, crse_rhs, frhs);
	}

}
void
Operator::restrictInteriorNodes (int camrlev, MultiFab& crhs, MultiFab& a_frhs) const // TODO
{
	TRACER;
	Util::Abort("multigrid not supported");
}


void
Operator::applyBC (int amrlev, int mglev, MultiFab& phi, BCMode/* bc_mode*/,
		   amrex::MLLinOp::StateMode /**/, bool skip_fillboundary) const // TODO
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

void
Operator::normalize (int amrlev, int mglev, MultiFab& mf) const
{
	TRACER;
	BL_PROFILE("MLNodeLaplacian::normalize()");
}

void
Operator::compSyncResidualCoarse (MultiFab& sync_resid, const MultiFab& a_phi,
				  const MultiFab& vold, const MultiFab* rhcc,
				  const BoxArray& fine_grids, const IntVect& ref_ratio)
{
	TRACER;
	Util::Abort("This function not implemented yet");
}

void
Operator::compSyncResidualFine (MultiFab& sync_resid, const MultiFab& phi, const MultiFab& vold,
				const MultiFab* rhcc)
{
	TRACER;
	Util::Abort("This function not implemented yet");
}

void
Operator::reflux (int crse_amrlev,
		  MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
		  MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const
{
	TRACER;
	BL_PROFILE("MLNodeLaplacian::reflux()");
	const Geometry& cgeom = m_geom[crse_amrlev  ][0];
	const Geometry& fgeom = m_geom[crse_amrlev+1][0];
	const Real* cdxinv = cgeom.InvCellSize();
	const Real* fdxinv = fgeom.InvCellSize();
	const Box& c_cc_domain = cgeom.Domain();
	const Box& c_nd_domain = amrex::surroundingNodes(c_cc_domain);

	const BoxArray& fba = fine_sol.boxArray();
	const DistributionMapping& fdm = fine_sol.DistributionMap();

	const iMultiFab& fdmsk = *m_dirichlet_mask[crse_amrlev+1][0];

	MultiFab fine_res_for_coarse(amrex::coarsen(fba, 2), fdm, 1, 0);

	applyBC(crse_amrlev+1, 0, fine_res, BCMode::Inhomogeneous, StateMode::Solution);

	for (MFIter mfi(fine_res_for_coarse, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		amrex_mlndlap_restriction(BL_TO_FORTRAN_BOX(bx),
					  BL_TO_FORTRAN_ANYD(fine_res_for_coarse[mfi]), // <<<<< OUTPUT
					  BL_TO_FORTRAN_ANYD(fine_res[mfi]),
					  BL_TO_FORTRAN_ANYD(fdmsk[mfi]),
					  BL_TO_FORTRAN_BOX(c_nd_domain),
					  m_lobc.data(), m_hibc.data());
	}
	res.ParallelCopy(fine_res_for_coarse, cgeom.periodicity());

	MultiFab fine_contrib(amrex::coarsen(fba, 2), fdm, 1, 0);
	fine_contrib.setVal(0.0);

	//const auto& fsigma = *m_sigma[crse_amrlev+1][0][0];

	{
		//FArrayBox sigfab;
		FArrayBox Axfab;
		for (MFIter mfi(fine_contrib, MFItInfo().EnableTiling().SetDynamic(true));
		     mfi.isValid(); ++mfi)
		{
			const Box& cvbx = mfi.validbox();
			const Box& fvbx = amrex::refine(cvbx,2);
			const Box& cbx = mfi.tilebox();
			const Box& fbx = amrex::refine(cbx,2);

			const Box& cc_fbx = amrex::enclosedCells(fbx);
			const Box& cc_fvbx = amrex::enclosedCells(fvbx);
			//const Box& bx_sig = amrex::grow(cc_fbx,2) & amrex::grow(cc_fvbx,1);
			//const Box& b = bx_sig & cc_fvbx;
			// sigfab.resize(bx_sig, 1);
			// sigfab.setVal(0.0);
			// sigfab.copy(fsigma[mfi], b, 0, b, 0, 1);

			const Box& bx_Ax = amrex::grow(fbx,1);
			const Box& b2 = bx_Ax & amrex::grow(fvbx,-1);
			Axfab.resize(bx_Ax);
			Axfab.setVal(0.0);
			Axfab.copy(fine_rhs[mfi], b2, 0, b2, 0, 1);
			Axfab.minus(fine_res[mfi], b2, 0, 0, 1);

			FineResidualContribution(crse_amrlev + 1,
						 0,
						 cbx, // box
						 cvbx, // cbbox
						 fine_contrib[mfi], //f (OUT)
						 fine_sol[mfi],//phi
						 Axfab, /// <<<< Ax (OUT)
						 fdmsk[mfi], // dmsk
						 fdxinv); // dxinv
		}
	}
	/*

	MultiFab fine_contrib_on_crse(crse_sol.boxArray(), crse_sol.DistributionMap(), 1, 0);
	fine_contrib_on_crse.setVal(0.0);
	fine_contrib_on_crse.ParallelAdd(fine_contrib, cgeom.periodicity());

	const iMultiFab& cdmsk = *m_dirichlet_mask[crse_amrlev][0];
	const auto& nd_mask     = m_nd_fine_mask[crse_amrlev];
	const auto& cc_mask     = m_cc_fine_mask[crse_amrlev];
	const auto& has_fine_bndry = m_has_fine_bndry[crse_amrlev];

	const auto& csigma = *m_sigma[crse_amrlev][0][0];

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(res, MFItInfo().EnableTiling().SetDynamic(true)); mfi.isValid(); ++mfi)
	{
		if ((*has_fine_bndry)[mfi])
		{
			const Box& bx = mfi.tilebox();
			amrex_mlndlap_res_cf_contrib(BL_TO_FORTRAN_BOX(bx),
						     BL_TO_FORTRAN_ANYD(res[mfi]), /// <<<< output
						     BL_TO_FORTRAN_ANYD(crse_sol[mfi]),
						     BL_TO_FORTRAN_ANYD(crse_rhs[mfi]),
						     BL_TO_FORTRAN_ANYD(csigma[mfi]),
						     BL_TO_FORTRAN_ANYD(cdmsk[mfi]),
						     BL_TO_FORTRAN_ANYD((*nd_mask)[mfi]),
						     BL_TO_FORTRAN_ANYD((*cc_mask)[mfi]),
						     BL_TO_FORTRAN_ANYD(fine_contrib_on_crse[mfi]),
						     cdxinv, BL_TO_FORTRAN_BOX(c_nd_domain),
						     m_lobc.data(), m_hibc.data());
		}
	}
	*/
}


//
// NON-INHERITED FUNCTIONS FOR MANAGING AUXILIARY DATA
// (Currently don't do anything)
//

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

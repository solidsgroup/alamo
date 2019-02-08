
#include "Operator.H"
#include <AMReX_MLLinOp_F.H>
#include <AMReX_MultiFabUtil.H>
#include "Util/Color.H"

namespace OperatorCell {

Operator::Operator ()
{
	m_ixtype = amrex::IntVect::TheCellVector();
}

Operator::~Operator () {}

void
Operator::define (amrex::Vector<amrex::Geometry> a_geom,
		   const amrex::Vector<amrex::BoxArray>& a_grids,
		   const amrex::Vector<amrex::DistributionMapping>& a_dmap,
		   BC::BC& a_bc,
		   const amrex::LPInfo& a_info,
		   const amrex::Vector<amrex::FabFactory<amrex::FArrayBox> const*>& a_factory)
{
	m_bc = &a_bc;

	std::array<int,AMREX_SPACEDIM> is_periodic = m_bc->IsPeriodic();
	// for (int ilev=0; ilev < a_geom.size(); ilev++)
	// 	a_geom[ilev].SetPeriodicity(is_periodic);
	
	MLCellLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory);
	defineAuxData();
	defineBC();


	//setDomainBC

	m_lobc = {AMREX_D_DECL(is_periodic[0] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet,
			       is_periodic[1] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet,
			       is_periodic[2] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet)};
	m_hibc = {AMREX_D_DECL(is_periodic[0] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet,
			       is_periodic[1] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet,
			       is_periodic[2] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet)};
				  
	for (int ilev = 0; ilev < a_geom.size(); ++ilev)
		setLevelBC(ilev,nullptr);

}

void
Operator::defineAuxData ()
{
	BL_PROFILE("Operator::defineAuxData()");

	m_undrrelxr.resize(m_num_amr_levels);
	m_maskvals.resize(m_num_amr_levels);
	//m_fluxreg.resize(m_num_amr_levels-1);

	const int ncomp = getNComp();

	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		m_undrrelxr[amrlev].resize(m_num_mg_levels[amrlev]);
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
		{
			m_undrrelxr[amrlev][mglev].define(m_grids[amrlev][mglev],
							  m_dmap[amrlev][mglev],
							  1, 0, 0, ncomp);
		}
	}
    
	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		m_maskvals[amrlev].resize(m_num_mg_levels[amrlev]);
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
		{
			for (amrex::OrientationIter oitr; oitr; ++oitr)
			{
				const amrex::Orientation face = oitr();
				const int ngrow = 1;
				const int extent = 1;
				m_maskvals[amrlev][mglev][face].define(m_grids[amrlev][mglev],
								       m_dmap[amrlev][mglev],
								       m_geom[amrlev][mglev],
								       face, 0, ngrow, extent, ncomp, true);
			}
		}
	}
}

void
Operator::defineBC ()
{
	BL_PROFILE("Operator::defineBC()");

	const int ncomp = getNComp();

	m_bndry_sol.resize(m_num_amr_levels);
	m_crse_sol_br.resize(m_num_amr_levels);

	m_bndry_cor.resize(m_num_amr_levels);
	m_crse_cor_br.resize(m_num_amr_levels);

	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		m_bndry_sol[amrlev].reset(new amrex::MLMGBndry(m_grids[amrlev][0], m_dmap[amrlev][0],
							       ncomp, m_geom[amrlev][0]));
	}

	for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
	{
		const int in_rad = 0;
		const int out_rad = 1;
		const int extent_rad = 2;
		const int crse_ratio = m_amr_ref_ratio[amrlev-1];
		amrex::BoxArray cba = m_grids[amrlev][0];
		cba.coarsen(crse_ratio);
		m_crse_sol_br[amrlev].reset(new amrex::BndryRegister(cba, m_dmap[amrlev][0],
								     in_rad, out_rad, extent_rad, ncomp));
	}

	for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
	{
		const int in_rad = 0;
		const int out_rad = 1;
		const int extent_rad = 2;
		const int crse_ratio = m_amr_ref_ratio[amrlev-1];
		amrex::BoxArray cba = m_grids[amrlev][0];
		cba.coarsen(crse_ratio);
		m_crse_cor_br[amrlev].reset(new amrex::BndryRegister(cba, m_dmap[amrlev][0],
								     in_rad, out_rad, extent_rad, ncomp));
		m_crse_cor_br[amrlev]->setVal(0.0);
	}

	// This has be to done after m_crse_cor_br is defined.
	for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
	{
		m_bndry_cor[amrlev].reset(new amrex::MLMGBndry(m_grids[amrlev][0], m_dmap[amrlev][0],
							       ncomp, m_geom[amrlev][0]));
		amrex::MultiFab bc_data(m_grids[amrlev][0], m_dmap[amrlev][0], ncomp, 1);
		bc_data.setVal(0.0);

		m_bndry_cor[amrlev]->setBndryValues(*m_crse_cor_br[amrlev], 0, bc_data, 0, 0, ncomp,
						    m_amr_ref_ratio[amrlev-1], amrex::BCRec());

		m_bndry_cor[amrlev]->setLOBndryConds({AMREX_D_DECL(BCType::Dirichlet,
								   BCType::Dirichlet,
								   BCType::Dirichlet)},
			{AMREX_D_DECL(BCType::Dirichlet,
				      BCType::Dirichlet,
				      BCType::Dirichlet)},
			m_amr_ref_ratio[amrlev-1], amrex::RealVect{});
	}

	m_bcondloc.resize(m_num_amr_levels);
	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		m_bcondloc[amrlev].resize(m_num_mg_levels[amrlev]);
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
		{
			m_bcondloc[amrlev][mglev].reset(new BndryCondLoc(m_grids[amrlev][mglev],
										m_dmap[amrlev][mglev]));
		} 
	}
}





void
Operator::interpolation (int /*amrlev*/, int /*fmglev*/,
			  amrex::MultiFab& fine, const amrex::MultiFab& crse) const
{
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (amrex::MFIter mfi(crse,true); mfi.isValid(); ++mfi)
	{
		const amrex::Box&         bx    = mfi.tilebox();
		const int          ncomp = getNComp();
		const amrex::FArrayBox& cfab    = crse[mfi];
		amrex::FArrayBox&       ffab    = fine[mfi];

		amrex_mg_interp(ffab.dataPtr(),
				AMREX_ARLIM(ffab.loVect()), AMREX_ARLIM(ffab.hiVect()),
				cfab.dataPtr(),
				AMREX_ARLIM(cfab.loVect()), AMREX_ARLIM(cfab.hiVect()),
				bx.loVect(), bx.hiVect(), &ncomp);
	}    
}

void
Operator::averageDownSolutionRHS (int camrlev, amrex::MultiFab& crse_sol, amrex::MultiFab& crse_rhs,
				   const amrex::MultiFab& fine_sol, const amrex::MultiFab& fine_rhs)
{
	const auto amrrr = AMRRefRatio(camrlev);
	const int ncomp = getNComp();
	amrex::average_down(fine_sol, crse_sol, 0, ncomp, amrrr);
	amrex::average_down(fine_rhs, crse_rhs, 0, ncomp, amrrr);
}

void
Operator::apply (int amrlev, int mglev,
		 amrex::MultiFab& out,
		 amrex::MultiFab& in,
		 BCMode bc_mode,
		 StateMode /*s_mode*/,
		 const amrex::MLMGBndry* bndry) const
{
	BL_PROFILE("Operator::apply()");
	applyBC(amrlev, mglev, in, bc_mode, bndry);
	Fapply(amrlev, mglev, out, in);
}

void
Operator::updateSolBC (int amrlev, const amrex::MultiFab& crse_bcdata) const
{
	BL_PROFILE("Operator::updateSolBC()");

	AMREX_ALWAYS_ASSERT(amrlev > 0);
	const int ncomp = getNComp();
	m_bc->define(m_geom[amrlev-1][0]);
	m_crse_sol_br[amrlev]->copyFrom(crse_bcdata, 0, 0, 0, ncomp, m_bc->Periodicity());
	m_bndry_sol[amrlev]->updateBndryValues(*m_crse_sol_br[amrlev], 0, 0, ncomp, m_amr_ref_ratio[amrlev-1]);
}

void
Operator::updateCorBC (int amrlev, const amrex::MultiFab& crse_bcdata) const
{
	BL_PROFILE("Operator::updateCorBC()");
	AMREX_ALWAYS_ASSERT(amrlev > 0);
	const int ncomp = getNComp();
	m_bc->define(m_geom[amrlev-1][0]);
	m_crse_cor_br[amrlev]->copyFrom(crse_bcdata, 0, 0, 0, ncomp, m_bc->Periodicity());
	m_bndry_cor[amrlev]->updateBndryValues(*m_crse_cor_br[amrlev], 0, 0, ncomp, m_amr_ref_ratio[amrlev-1]);
}


void
Operator::applyBC (int amrlev, int mglev, amrex::MultiFab& in, BCMode bc_mode,
		    const amrex::MLMGBndry* bndry, bool skip_fillboundary) const
{
	BL_PROFILE("Operator::applyBC()");
	// No coarsened boundary values, cannot apply inhomog at mglev>0.
	BL_ASSERT(mglev == 0 || bc_mode == BCMode::Homogeneous);
	BL_ASSERT(bndry != nullptr || bc_mode == BCMode::Homogeneous);

	m_bc->define(m_geom[amrlev][mglev]);

	const int ncomp = getNComp();
	const int cross = false;

	if (!skip_fillboundary) {
		in.FillBoundary(0, ncomp, m_bc->Periodicity(),cross); 
	}

	int flagbc = (bc_mode == BCMode::Homogeneous) ? 0 : 1;

	const amrex::Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

	const auto& maskvals = m_maskvals[amrlev][mglev];
	const auto& bcondloc = *m_bcondloc[amrlev][mglev];

	amrex::FArrayBox foo(amrex::Box::TheUnitBox(),ncomp);

#ifdef _OPENMP
#pragma omp parallel
#endif
	
	if (mglev==0 && bc_mode == BCMode::Inhomogeneous) 
	{
		m_bc->FillBoundary(in,0,0,0.0,0);
		auto *nonconst_this = const_cast<Operator*>(this);
		nonconst_this->setLevelBC(amrlev,&in);
	}

	//m_bc->define(m_geom[amrlev][mglev]);
	//m_bc->FillBoundary(bndry->bndryValues)
	for (amrex::MFIter mfi(in, amrex::MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
	{
		const amrex::Box& vbx   = mfi.validbox();
		amrex::FArrayBox& iofab = in[mfi];
		const RealTuple & bdl = bcondloc.bndryLocs(mfi);
		const BCTuple   & bdc = bcondloc.bndryConds(mfi);
		for (amrex::OrientationIter oitr; oitr; ++oitr)
		{
			const amrex::Orientation ori = oitr();
			int  cdr = ori;
			amrex::Real bcl = bdl[ori];
			int  bct = bdc[ori];
			foo.setVal(10.0);
			const amrex::FArrayBox& fsfab = (bndry != nullptr) ? bndry->bndryValues(ori)[mfi] : foo;
			const amrex::Mask& m = maskvals[ori][mfi];
			
			amrex_mllinop_apply_bc(BL_TO_FORTRAN_BOX(vbx),
			  		       BL_TO_FORTRAN_ANYD(iofab),
			 		       BL_TO_FORTRAN_ANYD(m),
			  		       cdr, bct, bcl,
			  		       BL_TO_FORTRAN_ANYD(fsfab),
			  		       maxorder, dxinv, flagbc, ncomp, cross);
		}
	}
}

void
Operator::compGrad (int amrlev, const std::array<amrex::MultiFab*,AMREX_SPACEDIM>& grad,
		    amrex::MultiFab& sol, amrex::MLLinOp::Location /*loc*/) const
{
	BL_PROFILE("Operator::compGrad()");

	if (sol.nComp() > 1)
		amrex::Abort("Operator::compGrad called, but only works for single-component solves");

	const int mglev = 0;

	// m_bc->define(m_geom[amrlev][mglev]);
	// m_bc->FillBoundary(sol,0,0,0.0);
	// auto *nonconst_this = const_cast<Operator*>(this);
	// nonconst_this->setLevelBC(amrlev,&sol);

	applyBC(amrlev, mglev, sol, BCMode::Inhomogeneous, m_bndry_sol[amrlev].get());

	const amrex::Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (amrex::MFIter mfi(sol, amrex::MFItInfo().EnableTiling().SetDynamic(true));  mfi.isValid(); ++mfi)
	{
		AMREX_D_TERM(const amrex::Box& xbx = mfi.nodaltilebox(0);,
			     const amrex::Box& ybx = mfi.nodaltilebox(1);,
			     const amrex::Box& zbx = mfi.nodaltilebox(2););
		amrex_mllinop_grad(AMREX_D_DECL(BL_TO_FORTRAN_BOX(xbx),
						BL_TO_FORTRAN_BOX(ybx),
						BL_TO_FORTRAN_BOX(zbx)),
				   BL_TO_FORTRAN_ANYD(sol[mfi]),
				   AMREX_D_DECL(BL_TO_FORTRAN_ANYD((*grad[0])[mfi]),
						BL_TO_FORTRAN_ANYD((*grad[1])[mfi]),
						BL_TO_FORTRAN_ANYD((*grad[2])[mfi])),
				   dxinv);
	}
}

void
Operator::prepareForSolve ()
{
	BL_PROFILE("Operator::prepareForSolve()");

	const int ncomp = getNComp();
	for (int amrlev = 0;  amrlev < m_num_amr_levels; ++amrlev)
	{
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
		{
			const auto& bcondloc = *m_bcondloc[amrlev][mglev];
			const auto& maskvals = m_maskvals[amrlev][mglev];
			const amrex::Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

			amrex::BndryRegister& undrrelxr = m_undrrelxr[amrlev][mglev];
			amrex::MultiFab foo(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], ncomp, 0, amrex::MFInfo().SetAlloc(false));
#ifdef _OPENMP
#pragma omp parallel
#endif
			for (amrex::MFIter mfi(foo, amrex::MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
			{
				const amrex::Box& vbx = mfi.validbox();

				const RealTuple & bdl = bcondloc.bndryLocs(mfi);
				const BCTuple   & bdc = bcondloc.bndryConds(mfi);

				for (amrex::OrientationIter oitr; oitr; ++oitr)
				{
					const amrex::Orientation ori = oitr();
                    
					int  cdr = ori;
					amrex::Real bcl = bdl[ori];
					int  bct = bdc[ori];
                    
					amrex::FArrayBox& ffab = undrrelxr[ori][mfi];
					const amrex::Mask& m   =  maskvals[ori][mfi];

					amrex_mllinop_comp_interp_coef0(BL_TO_FORTRAN_BOX(vbx),
									BL_TO_FORTRAN_ANYD(ffab),
									BL_TO_FORTRAN_ANYD(m),
									cdr, bct, bcl, maxorder, dxinv, ncomp);
				}
			}
		}
	}
	averageDownCoeffs();
}

Operator::BndryCondLoc::BndryCondLoc (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
	: bcond(ba, dm),
	  bcloc(ba, dm)
{
}

void
Operator::BndryCondLoc::setLOBndryConds (const amrex::Geometry& geom, const amrex::Real* dx,
					  const std::array<BCType,AMREX_SPACEDIM>& lobc,
					  const std::array<BCType,AMREX_SPACEDIM>& hibc,
					  int ratio, const amrex::RealVect& a_loc)
{
	const amrex::Box&  domain = geom.Domain();

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (amrex::MFIter mfi(bcloc); mfi.isValid(); ++mfi)
	{
		const amrex::Box& bx = mfi.validbox();
		RealTuple & bloc  = bcloc[mfi];
		BCTuple   & bctag = bcond[mfi];

		amrex::MLMGBndry::setBoxBC(bloc, bctag, bx, domain, lobc, hibc, dx, ratio, a_loc);
	}
}


void
Operator::averageDownCoeffs ()
{
	for (int i = 0; i < m_num_a_fabs; i++)
	{
		for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
		{
			auto& fine_a_coeffs = m_a_coeffs[i][amrlev];
			averageDownCoeffsSameAmrLevel(fine_a_coeffs);
		}
		averageDownCoeffsSameAmrLevel(m_a_coeffs[i][0]);
	}
}

void
Operator::averageDownCoeffsSameAmrLevel (amrex::Vector<amrex::MultiFab>& a)
{
	int nmglevs = a.size();
	for (int mglev = 1; mglev < nmglevs; ++mglev)
	{
		amrex::average_down(a[mglev-1], a[mglev], 0, a[0].nComp(), mg_coarsen_ratio);
	}
}



const amrex::FArrayBox &
Operator::GetFab(const int num, const int amrlev, const int mglev, const amrex::MFIter &mfi) const
{
	return m_a_coeffs[num][amrlev][mglev][mfi];
}


void
Operator::RegisterNewFab(amrex::Vector<amrex::MultiFab> &input)
{
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

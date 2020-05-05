
#include <AMReX_MLNodeLinOp.H>
#include <AMReX_MLCellLinOp.H>
#include <AMReX_MLNodeLap_K.H>
#include <AMReX_MultiFabUtil.H>
#include "Util/Color.H"
#include "Set/Set.H"
#include "Operator.H"

using namespace amrex;
namespace Operator {

// constexpr amrex::IntVect AMREX_D_DECL(Operator<Grid::Node>::dx,Operator<Grid::Node>::dy,Operator<Grid::Node>::dz);
constexpr amrex::IntVect AMREX_D_DECL(Operator<Grid::Cell>::dx,Operator<Grid::Cell>::dy,Operator<Grid::Cell>::dz);

void Operator<Grid::Node>::Diagonal (bool recompute)
{
	BL_PROFILE(Color::FG::Yellow + "Operator::Diagonal()" + Color::Reset);
	//Util::Message(INFO);
	
	if ( !recompute && m_diagonal_computed ) return;
	m_diagonal_computed = true;

	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
		{
			Diagonal(amrlev,mglev,*m_diag[amrlev][mglev]);
		}
	}
}

void Operator<Grid::Node>::Diagonal (int amrlev, int mglev, amrex::MultiFab &diag)
{
	BL_PROFILE("Operator::Diagonal()");
	//Util::Message(INFO);

	int ncomp = diag.nComp();
	int nghost = 0;

	int sep = 2;
	int num = AMREX_D_TERM(sep,*sep,*sep);
	int cntr = 0;

	amrex::MultiFab x(m_diag[amrlev][mglev]->boxArray(), m_diag[amrlev][mglev]->DistributionMap(), ncomp, nghost);
	amrex::MultiFab Ax(m_diag[amrlev][mglev]->boxArray(), m_diag[amrlev][mglev]->DistributionMap(), ncomp, nghost);

	for (MFIter mfi(x, false); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.validbox();
		amrex::FArrayBox       &diagfab = diag[mfi];
		amrex::FArrayBox       &xfab    = x[mfi];
		amrex::FArrayBox       &Axfab   = Ax[mfi];

		diagfab.setVal(0.0);
				
		for (int i = 0; i < num; i++)
		{
			for (int n = 0; n < ncomp; n++)
			{
				xfab.setVal(0.0);
				Axfab.setVal(0.0);
				
				//BL_PROFILE_VAR("Operator::Part1", part1); 
				AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
					     for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
					     for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
				{
					amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
				
					if ( m1%sep == i/sep   &&   m2%sep == i%sep ) xfab(m,n) = 1.0;
					else xfab(m,n) = 0.0;
				}
				//BL_PROFILE_VAR_STOP(part1);

				BL_PROFILE_VAR("Operator::Part2", part2); 
				Util::Message(INFO,"Calling fapply...",cntr++);
				Fapply(amrlev,mglev,Ax,x);
				BL_PROFILE_VAR_STOP(part2);
						
				//BL_PROFILE_VAR("Operator::Part3", part3); 
				Axfab.mult(xfab,n,n,1);
				diagfab.plus(Axfab,n,n,1);
				//BL_PROFILE_VAR_STOP(part3);
			}
		}
	}
}

void Operator<Grid::Node>::Fsmooth (int amrlev, int mglev, amrex::MultiFab& x, const amrex::MultiFab& b) const
{
	BL_PROFILE("Operator::Fsmooth()");

	amrex::Box domain(m_geom[amrlev][mglev].Domain());

	int ncomp = b.nComp();
	int nghost = 2; //b.nGrow();
	
	Set::Scalar omega = 2./3.; // Damping factor (very important!)

	amrex::MultiFab Ax(x.boxArray(), x.DistributionMap(), ncomp, nghost);
	amrex::MultiFab Dx(x.boxArray(), x.DistributionMap(), ncomp, nghost);
	amrex::MultiFab Rx(x.boxArray(), x.DistributionMap(), ncomp, nghost);
	
	if (!m_diagonal_computed) Util::Abort(INFO,"Operator::Diagonal() must be called before using Fsmooth");

	// This is a JACOBI iteration, not Gauss-Seidel.
	// So we need to do twice the number of iterations to get the same behavior as GS.
	for (int ctr = 0; ctr < 2; ctr++)
	{
		Fapply(amrlev,mglev,Ax,x); // find Ax

		amrex::MultiFab::Copy(Dx,x,0,0,ncomp,nghost); // Dx = x
		amrex::MultiFab::Multiply(Dx,*m_diag[amrlev][mglev],0,0,ncomp,nghost); // Dx *= diag  (Dx = x*diag)

		amrex::MultiFab::Copy(Rx,Ax,0,0,ncomp,nghost); // Rx = Ax
		amrex::MultiFab::Subtract(Rx,Dx,0,0,ncomp,nghost); // Rx -= Dx  (Rx = Ax - Dx)

		for (MFIter mfi(x, false); mfi.isValid(); ++mfi)
		{
			const Box& bx = mfi.validbox();
			amrex::FArrayBox       &xfab    = x[mfi];
			const amrex::FArrayBox &bfab    = b[mfi];
			//const amrex::FArrayBox &Axfab   = Ax[mfi];
			const amrex::FArrayBox &Rxfab   = Rx[mfi];
			const amrex::FArrayBox &diagfab = (*m_diag[amrlev][mglev])[mfi];

			for (int n = 0; n < ncomp; n++)
			{
				AMREX_D_TERM(for (int m1 = bx.loVect()[0] - 2; m1<=bx.hiVect()[0] + 2; m1++),
					     for (int m2 = bx.loVect()[1] - 2; m2<=bx.hiVect()[1] + 2; m2++),
					     for (int m3 = bx.loVect()[2] - 2; m3<=bx.hiVect()[2] + 2; m3++))
				{

					amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));

					// Skip ghost cells outside problem domain
					if (AMREX_D_TERM(m[0] < domain.loVect()[0], ||
							 m[1] < domain.loVect()[1], ||
							 m[2] < domain.loVect()[2])) continue;
					if (AMREX_D_TERM(m[0] > domain.hiVect()[0] + 1, ||
							 m[1] > domain.hiVect()[1] + 1, ||
							 m[2] > domain.hiVect()[2] + 1)) continue;

					if (AMREX_D_TERM(m[0] == bx.loVect()[0] - nghost || m[0] == bx.hiVect()[0] + nghost, ||
							 m[1] == bx.loVect()[1] - nghost || m[1] == bx.hiVect()[1] + nghost, ||
							 m[2] == bx.loVect()[2] - nghost || m[2] == bx.hiVect()[2] + nghost))
					{
						xfab(m,n) = 0.0;
						continue;
					}

					//xfab(m,n) = xfab(m,n) + omega*(bfab(m,n) - Axfab(m,n))/diagfab(m,n);
					xfab(m,n) = (1.-omega)*xfab(m,n) + omega*(bfab(m,n) - Rxfab(m,n))/diagfab(m,n);
				}
			}
		}
	}
	amrex::Geometry geom = m_geom[amrlev][mglev];
	realFillBoundary(x,geom);
	nodalSync(amrlev, mglev, x);
}

void Operator<Grid::Node>::normalize (int amrlev, int mglev, MultiFab& a_x) const
{
	BL_PROFILE("Operator::normalize()");
	amrex::Box domain(m_geom[amrlev][mglev].Domain());
	domain.convert(amrex::IntVect::TheNodeVector());

	int ncomp = getNComp();
	int nghost = 1; //x.nGrow();

	if (!m_diagonal_computed)
		Util::Abort(INFO,"Operator::Diagonal() must be called before using normalize");

	for (MFIter mfi(a_x, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{

		Box bx = mfi.tilebox();
		bx.grow(nghost);
		bx = bx & domain;

		amrex::Array4<amrex::Real> const& x = a_x.array(mfi);
		amrex::Array4<const amrex::Real> const& diag = m_diag[amrlev][mglev]->array(mfi);

		for (int n = 0; n < ncomp; n++)
		{
			amrex::ParallelFor (bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
					
					x(i,j,k) /= diag(i,j,k);

				} );
		}
	}
}

Operator<Grid::Node>::Operator (const Vector<Geometry>& a_geom,
		    const Vector<BoxArray>& a_grids,
		    const Vector<DistributionMapping>& a_dmap,
		    const LPInfo& a_info,
		    const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
	BL_PROFILE("Operator::Operator()");
	Util::Message(INFO);

	if (!(a_grids[0].ixType()  == amrex::IndexType::TheNodeType()))
		Util::Abort(INFO,"Operator must be defined using CELL CENTERED boxarrays.");

	define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

 Operator<Grid::Node>::~Operator ()
 {}

void Operator<Grid::Node>::define (const Vector<Geometry>& a_geom,
		       const Vector<BoxArray>& a_grids,
		       const Vector<DistributionMapping>& a_dmap,
		       const LPInfo& a_info,
		       const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
	BL_PROFILE("Operator::~Operator()");

	 // Make sure we're not trying to parallelize in vain.
	 if (amrex::ParallelDescriptor::NProcs() > a_grids[0].size())
	 {
		 Util::Warning(INFO,"There are more processors than there are boxes in the amrlev=0 boxarray!!\n",
			       "(NProcs = ",amrex::ParallelDescriptor::NProcs(),", a_grids[0].size() = ",a_grids[0].size(),")\n",
			       "You should decrease max_grid_size or you will not get proper scaling!");
	 }

	 // This makes sure grids are node-centered;
	 Vector<BoxArray> cc_grids = a_grids;
	 for (auto& ba : cc_grids) {
		 ba.enclosedCells();
	 }

	 MLNodeLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory);

	 int nghost = 2;
	 // Resize the multifab containing the operator diagonal
	 m_diag.resize(m_num_amr_levels);
	 for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	 {
		 m_diag[amrlev].resize(m_num_mg_levels[amrlev]);

		 for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
		 {
		 	 m_diag[amrlev][mglev].reset(new MultiFab(amrex::convert(m_grids[amrlev][mglev], amrex::IntVect::TheNodeVector()),
		 						  m_dmap[amrlev][mglev], getNComp(), nghost));
		 }
	 }

	// We need to instantiate the m_lobc objects.
	// WE DO NOT USE THEM - our BCs are implemented differently.
	// But they need to be the right size or the code will segfault.
	m_lobc.resize(getNComp(),{{AMREX_D_DECL(BCType::bogus,BCType::bogus,BCType::bogus)}});
	m_hibc.resize(getNComp(),{{AMREX_D_DECL(BCType::bogus,BCType::bogus,BCType::bogus)}});
}


void Operator<Grid::Node>::buildMasks ()
{
	BL_PROFILE("Operator::buildMasks()");
	if (m_masks_built) return;


	m_masks_built = true;

	m_is_bottom_singular = false;

	const auto lobc = LoBC();
	const auto hibc = HiBC();

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
				//const Box& nddomain = amrex::surroundingNodes(ccdomain);
				const std::vector<IntVect>& pshifts = period.shiftIntVect();

				Box ccdomain_p = ccdomain;
				for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
					if (geom.isPeriodic(idim)) {
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
				const Box& bx = mfi.fabbox();
				Array4<int> const& fab = cc_mask.array(mfi);
				for (const auto& iv : pshifts)
				{
					cfba.intersections(bx+iv, isects);
					for (const auto& is : isects)
					{
						Box const& b = is.second-iv;
                        			AMREX_HOST_DEVICE_PARALLEL_FOR_3D(b,i,j,k,
                        			{
                            				fab(i,j,k) = 1;
                        			});
                    			}
					if (!isects.empty()) has_cf[mfi] = 1;
				}

				mlndlap_fillbc_cc<int>(mfi.validbox(),fab,ccdom,lobc,hibc);
			}
		}

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
		for (MFIter mfi(nd_mask,TilingIfNotGPU()); mfi.isValid(); ++mfi)
		{
			const Box& bx = mfi.tilebox();
			Array4<int> const& nmsk = nd_mask.array(mfi);
			Array4<int const> const& cmsk = cc_mask.const_array(mfi);
			AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
			{
				mlndlap_set_nodal_mask(i,j,k,nmsk,cmsk);
			});
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


#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
		for (MFIter mfi(m_bottom_dot_mask,true); mfi.isValid(); ++mfi)
		{
			const Box& bx = mfi.tilebox();
			Array4<Real> const& dfab = m_bottom_dot_mask.array(mfi);
			Array4<int const> const& sfab = omask.const_array(mfi);
			AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
			{
				mlndlap_set_dot_mask(tbx, dfab, sfab, nddomain, lobc, hibc);
			});
		}
	}
}


void Operator<Grid::Node>::fixUpResidualMask (int amrlev, iMultiFab& resmsk)
{
	BL_PROFILE("Operator::fixUpResidualMask()");

	if (!m_masks_built) buildMasks();

	const iMultiFab& cfmask = *m_nd_fine_mask[amrlev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	for (MFIter mfi(resmsk,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		Array4<int> const& rmsk = resmsk.array(mfi);
		Array4<int const> const& fmsk = cfmask.const_array(mfi);
		AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
		{
			if (fmsk(i,j,k) == crse_fine_node) rmsk(i,j,k) = 1;
		});
	}
}

void Operator<Grid::Node>::prepareForSolve ()
{
	BL_PROFILE("Operator::prepareForSolve()");
	MLNodeLinOp::prepareForSolve();
	buildMasks();
	averageDownCoeffs();
	Diagonal(true);
}

void Operator<Grid::Node>::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
	BL_PROFILE("Operator::restriction()");

	applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous, StateMode::Solution);

	amrex::Box cdomain = m_geom[amrlev][cmglev].Domain();
	cdomain.convert(amrex::IntVect::TheNodeVector());

	bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
	MultiFab cfine;
	if (need_parallel_copy) {
		const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
		cfine.define(ba, fine.DistributionMap(), fine.nComp(), fine.nGrow());
	}

	MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;

	for (MFIter mfi(*pcrse, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();

		amrex::Array4<const amrex::Real> const& fdata = fine.array(mfi);
		amrex::Array4<amrex::Real> const& cdata       = pcrse->array(mfi);

		const Dim3 lo= amrex::lbound(cdomain), hi = amrex::ubound(cdomain);


		for (int n = 0; n < crse.nComp(); n++)
		{
			// I,J,K == coarse coordinates
			// i,j,k == fine coordinates
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int I, int J, int K) {
					int i=2*I, j=2*J, k=2*K;

					if ((I == lo.x || I == hi.x) &&
					    (J == lo.y || J == hi.y) &&
					    (K == lo.z || K == hi.z)) // Corner
						{
							cdata(I,J,K,n) = fdata(i,j,k,n);
						}
					else if ((J == lo.y || J == hi.y) &&
						 (K == lo.z || K == hi.z)) // X edge
						 {
							 cdata(I,J,K,n) = 0.25*fdata(i-1,j,k,n) + 0.5*fdata(i,j,k,n) + 0.25*fdata(i+1,j,k,n);
						 }
					else if ((K == lo.z || K == hi.z) &&
						 (I == lo.x || I == hi.x)) // Y edge
						{
							cdata(I,J,K,n) = 0.25*fdata(i,j-1,k,n) + 0.5*fdata(i,j,k,n) + 0.25*fdata(i,j+1,k,n);
						}
					else if ((I == lo.x || I == hi.x) &&
						 (J == lo.y || J == hi.y)) // Z edge
						{
							cdata(I,J,K,n) = 0.25*fdata(i,j,k-1,n) + 0.5*fdata(i,j,k,n) + 0.25*fdata(i,j,k+1,n);
						}
					else if (I == lo.x || I == hi.x) // X face
						{
							cdata(I,J,K,n) =
							(+     fdata(i,j-1,k-1,n) + 2.0*fdata(i,j,k-1,n) +     fdata(i,j+1,k-1,n)
							 + 2.0*fdata(i,j-1,k  ,n) + 4.0*fdata(i,j,k  ,n) + 2.0*fdata(i,j+1,k  ,n) 
							 +     fdata(i,j-1,k+1,n) + 2.0*fdata(i,j,k+1,n) +     fdata(i,j+1,k+1,n))/16.0;
						}
					else if (J == lo.y || J == hi.y) // Y face
						{
							cdata(I,J,K,n) =
							(+     fdata(i-1,j,k-1,n) + 2.0*fdata(i-1,j,k,n) +     fdata(i-1,j,k+1,n)
							 + 2.0*fdata(i  ,j,k-1,n) + 4.0*fdata(i  ,j,k,n) + 2.0*fdata(i  ,j,k+1,n) 
							 +     fdata(i+1,j,k-1,n) + 2.0*fdata(i+1,j,k,n) +     fdata(i+1,j,k+1,n))/16.0;
						}
					else if (K == lo.z || K == hi.z) // Z face
						{cdata(I,J,K,n) =
							(+     fdata(i-1,j-1,k,n) + 2.0*fdata(i,j-1,k,n) +     fdata(i+1,j-1,k,n)
							 + 2.0*fdata(i-1,j  ,k,n) + 4.0*fdata(i,j  ,k,n) + 2.0*fdata(i+1,j  ,k,n) 
							 +     fdata(i-1,j+1,k,n) + 2.0*fdata(i,j+1,k,n) +     fdata(i+1,j+1,k,n))/16.0;
							 }
					else // Interior
						cdata(I,J,K,n) =
							(fdata(i-1,j-1,k-1,n) + fdata(i-1,j-1,k+1,n) + fdata(i-1,j+1,k-1,n) + fdata(i-1,j+1,k+1,n) +
							 fdata(i+1,j-1,k-1,n) + fdata(i+1,j-1,k+1,n) + fdata(i+1,j+1,k-1,n) + fdata(i+1,j+1,k+1,n)) / 64.0
							+
							(fdata(i,j-1,k-1,n) + fdata(i,j-1,k+1,n) + fdata(i,j+1,k-1,n) + fdata(i,j+1,k+1,n) +
							 fdata(i-1,j,k-1,n) + fdata(i+1,j,k-1,n) + fdata(i-1,j,k+1,n) + fdata(i+1,j,k+1,n) +
							 fdata(i-1,j-1,k,n) + fdata(i-1,j+1,k,n) + fdata(i+1,j-1,k,n) + fdata(i+1,j+1,k,n)) / 32.0
							+
							(fdata(i-1,j,k,n) + fdata(i,j-1,k,n) + fdata(i,j,k-1,n) +
							 fdata(i+1,j,k,n) + fdata(i,j+1,k,n) + fdata(i,j,k+1,n)) / 16.0
							+
							fdata(i,j,k,n) / 8.0;
				});
		}
	}

	if (need_parallel_copy) {
		crse.ParallelCopy(cfine);
	}

	amrex::Geometry geom = m_geom[amrlev][cmglev];
	realFillBoundary(crse,geom);
	nodalSync(amrlev, cmglev, crse);
}

void Operator<Grid::Node>::interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const
{
	BL_PROFILE("Operator::interpolation()");
	//int nghost = getNGrow();
	amrex::Box fdomain = m_geom[amrlev][fmglev].Domain(); fdomain.convert(amrex::IntVect::TheNodeVector());
	
	bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
	MultiFab cfine;
	const MultiFab* cmf = &crse;
	if (need_parallel_copy) {
		const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
		cfine.define(ba, fine.DistributionMap(), crse.nComp(), crse.nGrow());
		cfine.ParallelCopy(crse);
		cmf = &cfine;
	}

	for (MFIter mfi(fine, false); mfi.isValid(); ++mfi)
	{
		const Box& fine_bx = mfi.validbox() & fdomain;
		const Box& course_bx = amrex::coarsen(fine_bx,2);
		const Box& tmpbx = amrex::refine(course_bx,2);
		FArrayBox tmpfab;
		tmpfab.resize(tmpbx,fine.nComp());
		tmpfab.setVal(0.0);
		const amrex::FArrayBox &crsefab = (*cmf)[mfi];

		amrex::Array4<const amrex::Real> const& cdata = crsefab.const_array();
		amrex::Array4<amrex::Real> const& fdata       = tmpfab.array();

		for (int n = 0; n < crse.nComp(); n++)
		{
			// I,J,K == coarse coordinates
			// i,j,k == fine coordinates
			amrex::ParallelFor (fine_bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
					
					int I=i/2, J=j/2, K=k/2;

					if (i%2 == 0 && j%2 == 0 && k%2 ==0) // Coincident
						fdata(i,j,k,n) = cdata(I,J,K,n);
					else if (j%2 == 0 && k%2 == 0) // X Edge
						fdata(i,j,k,n) = 0.5 * (cdata(I,J,K,n) + cdata(I+1,J,K,n));
					else if (k%2 == 0 && i%2 == 0) // Y Edge
						fdata(i,j,k,n) = 0.5 * (cdata(I,J,K,n) + cdata(I,J+1,K,n));
					else if (i%2 == 0 && j%2 == 0) // Z Edge
						fdata(i,j,k,n) = 0.5 * (cdata(I,J,K,n) + cdata(I,J,K+1,n)); 
					else if (i%2 == 0) // X Face
						fdata(i,j,k,n) = 0.25 * (cdata(I,J,K,n)   + cdata(I,J+1,K,n) +
									 cdata(I,J,K+1,n) + cdata(I,J+1,K+1,n));
					else if (j%2 == 0) // Y Face
						fdata(i,j,k,n) = 0.25 * (cdata(I,J,K,n)   + cdata(I,J,K+1,n) +
									 cdata(I+1,J,K,n) + cdata(I+1,J,K+1,n));
					else if (k%2 == 0) // Z Face
						fdata(i,j,k,n) = 0.25 * (cdata(I,J,K,n)   + cdata(I+1,J,K,n) +
									 cdata(I,J+1,K,n) + cdata(I+1,J+1,K,n));
					else // Center
						fdata(i,j,k,n) = 0.125 * (cdata(I,J,K,n) +
									  cdata(I+1,J,K,n)   + cdata(I,J+1,K,n)   + cdata(I,J,K+1,n) +
									  cdata(I,J+1,K+1,n) + cdata(I+1,J,K+1,n) + cdata(I+1,J+1,K,n) +
									  cdata(I+1,J+1,K+1,n));

				});
		}
		fine[mfi].plus(tmpfab,fine_bx,fine_bx,0,0,fine.nComp());
	}

	amrex::Geometry geom = m_geom[amrlev][fmglev];
	realFillBoundary(fine,geom);
	nodalSync(amrlev, fmglev, fine);
}
  
void Operator<Grid::Node>::averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& /*crse_rhs*/,
				                                  const MultiFab& fine_sol, const MultiFab& /*fine_rhs*/)
{
	BL_PROFILE("Operator::averageDownSolutionRHS()");
	const auto& amrrr = AMRRefRatio(camrlev);
	amrex::average_down(fine_sol, crse_sol, 0, crse_sol.nComp(), amrrr);
	
	if (isSingular(0))
	{
		Util::Abort(INFO,"Singular operators not supported!");
		// MultiFab frhs(fine_rhs.boxArray(), fine_rhs.DistributionMap(), 1, 1);
		// MultiFab::Copy(frhs, fine_rhs, 0, 0, 1, 0);
		// restrictInteriorNodes(camrlev, crse_rhs, frhs);
	}

}

void Operator<Grid::Node>::realFillBoundary(MultiFab &phi, const Geometry &geom) 
{
	for (int i = 0; i < 2; i++)
	{
		MultiFab & mf = phi;
		mf.FillBoundary(geom.periodicity());
		const int ncomp = mf.nComp();
		const int ng1 = 1;
		const int ng2 = 2;
		MultiFab tmpmf(mf.boxArray(), mf.DistributionMap(), ncomp, ng1);
		MultiFab::Copy(tmpmf, mf, 0, 0, ncomp, ng1); 
		mf.ParallelCopy   (tmpmf, 0, 0, ncomp, ng1, ng2, geom.periodicity());
	}
}

void Operator<Grid::Node>::applyBC (int amrlev, int mglev, MultiFab& phi, BCMode/* bc_mode*/,
		   amrex::MLLinOp::StateMode /**/, bool skip_fillboundary) const
{
	BL_PROFILE("Operator::applyBC()");

	const Geometry& geom = m_geom[amrlev][mglev];

	if (!skip_fillboundary) {

		realFillBoundary(phi,geom);
	}
}

const amrex::FArrayBox &
Operator<Grid::Node>::GetFab(const int num, const int amrlev, const int mglev, const amrex::MFIter &mfi) const
{
	BL_PROFILE("Operator::GetFab()");
	Util::Message(INFO);
 	return m_a_coeffs[num][amrlev][mglev][mfi];
}

void Operator<Grid::Node>::RegisterNewFab(amrex::Vector<amrex::MultiFab> &input)
{
	BL_PROFILE("Operator::RegisterNewFab()");
	Util::Message(INFO);
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


void Operator<Grid::Node>::RegisterNewFab(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &input)
{
	BL_PROFILE("Operator::RegisterNewFab()");
	Util::Message(INFO);
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

void Operator<Grid::Node>::reflux (int crse_amrlev,
		  MultiFab& res, const MultiFab& /*crse_sol*/, const MultiFab& /*crse_rhs*/,
		  MultiFab& fine_res, MultiFab& /*fine_sol*/, const MultiFab& /*fine_rhs*/) const
{
	BL_PROFILE("Operator::Elastic::reflux()");

	int ncomp = AMREX_SPACEDIM;

	amrex::Box cdomain(m_geom[crse_amrlev][0].Domain());
	cdomain.convert(amrex::IntVect::TheNodeVector());

	const Geometry& cgeom = m_geom[crse_amrlev  ][0];

 	const BoxArray&            fba = fine_res.boxArray();
 	const DistributionMapping& fdm = fine_res.DistributionMap();

 	MultiFab fine_res_for_coarse(amrex::coarsen(fba, 2), fdm, ncomp, 2);
	fine_res_for_coarse.ParallelCopy(res,0,0,ncomp,0,0,cgeom.periodicity());

 	applyBC(crse_amrlev+1, 0, fine_res, BCMode::Inhomogeneous, StateMode::Solution);

	/// \todo Replace with Enum
	// const int coarse_coarse_node = 0;
	const int coarse_fine_node = 1;
	const int fine_fine_node = 2;

	amrex::iMultiFab nodemask(amrex::coarsen(fba,2), fdm, 1, 2);
	nodemask.ParallelCopy(*m_nd_fine_mask[crse_amrlev],0,0,1,0,0,cgeom.periodicity());

	amrex::iMultiFab cellmask(amrex::convert(amrex::coarsen(fba,2),amrex::IntVect::TheCellVector()), fdm, 1, 2);
	cellmask.ParallelCopy(*m_cc_fine_mask[crse_amrlev],0,0,1,1,1,cgeom.periodicity());
	
	for (MFIter mfi(fine_res_for_coarse, false); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.validbox();

		amrex::Array4<const int> const& nmask = nodemask.array(mfi);
		//amrex::Array4<const int> const& cmask = cellmask.array(mfi);

		amrex::Array4<amrex::Real> const& cdata = fine_res_for_coarse.array(mfi);
		amrex::Array4<const amrex::Real> const& fdata       = fine_res.array(mfi);

		const Dim3 lo= amrex::lbound(cdomain), hi = amrex::ubound(cdomain);

		for (int n = 0; n < fine_res.nComp(); n++)
		{
			// I,J,K == coarse coordinates
			// i,j,k == fine coordinates
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int I, int J, int K) {
					int i=I*2, j=J*2, k=K*2;
					
					if (nmask(I,J,K) == fine_fine_node || nmask(I,J,K) == coarse_fine_node)
						{
							if ((I == lo.x || I == hi.x) &&
							    (J == lo.y || J == hi.y) &&
							    (K == lo.z || K == hi.z)) // Corner
								cdata(I,J,K,n) = fdata(i,j,k,n);
							else if ((J == lo.y || J == hi.y) &&
								 (K == lo.z || K == hi.z)) // X edge
								cdata(I,J,K,n) = 0.25*fdata(i-1,j,k,n) + 0.5*fdata(i,j,k,n) + 0.25*fdata(i+1,j,k,n);
							else if ((K == lo.z || K == hi.z) &&
								 (I == lo.x || I == hi.x)) // Y edge
								cdata(I,J,K,n) = 0.25*fdata(i,j-1,k,n) + 0.5*fdata(i,j,k,n) + 0.25*fdata(i,j+1,k,n);
							else if ((I == lo.x || I == hi.x) &&
								 (J == lo.y || J == hi.y)) // Z edge
								cdata(I,J,K,n) = 0.25*fdata(i,j,k-1,n) + 0.5*fdata(i,j,k,n) + 0.25*fdata(i,j,k+1,n);
							else if (I == lo.x || I == hi.x) // X face
								cdata(I,J,K,n) =
									(+     fdata(i,j-1,k-1,n) + 2.0*fdata(i,j,k-1,n) +     fdata(i,j+1,k-1,n)
									 + 2.0*fdata(i,j-1,k  ,n) + 4.0*fdata(i,j,k  ,n) + 2.0*fdata(i,j+1,k  ,n) 
									 +     fdata(i,j-1,k+1,n) + 2.0*fdata(i,j,k+1,n) +     fdata(i,j+1,k+1,n))/16.0;
							else if (J == lo.y || J == hi.y) // Y face
								cdata(I,J,K,n) =
									(+     fdata(i-1,j,k-1,n) + 2.0*fdata(i-1,j,k,n) +     fdata(i-1,j,k+1,n)
									 + 2.0*fdata(i  ,j,k-1,n) + 4.0*fdata(i  ,j,k,n) + 2.0*fdata(i  ,j,k+1,n) 
									 +     fdata(i+1,j,k-1,n) + 2.0*fdata(i+1,j,k,n) +     fdata(i+1,j,k+1,n))/16.0;
							else if (K == lo.z || K == hi.z) // Z face
								cdata(I,J,K,n) =
									(+     fdata(i-1,j-1,k,n) + 2.0*fdata(i,j-1,k,n) +     fdata(i+1,j-1,k,n)
									 + 2.0*fdata(i-1,j  ,k,n) + 4.0*fdata(i,j  ,k,n) + 2.0*fdata(i+1,j  ,k,n) 
									 +     fdata(i-1,j+1,k,n) + 2.0*fdata(i,j+1,k,n) +     fdata(i+1,j+1,k,n))/16.0;
							else // Interior
								cdata(I,J,K,n) =
									(fdata(i-1,j-1,k-1,n) + fdata(i-1,j-1,k+1,n) + fdata(i-1,j+1,k-1,n) + fdata(i-1,j+1,k+1,n) +
									 fdata(i+1,j-1,k-1,n) + fdata(i+1,j-1,k+1,n) + fdata(i+1,j+1,k-1,n) + fdata(i+1,j+1,k+1,n)) / 64.0
									+
									(fdata(i,j-1,k-1,n) + fdata(i,j-1,k+1,n) + fdata(i,j+1,k-1,n) + fdata(i,j+1,k+1,n) +
									 fdata(i-1,j,k-1,n) + fdata(i+1,j,k-1,n) + fdata(i-1,j,k+1,n) + fdata(i+1,j,k+1,n) +
									 fdata(i-1,j-1,k,n) + fdata(i-1,j+1,k,n) + fdata(i+1,j-1,k,n) + fdata(i+1,j+1,k,n)) / 32.0
									+
									(fdata(i-1,j,k,n) + fdata(i,j-1,k,n) + fdata(i,j,k-1,n) +
									 fdata(i+1,j,k,n) + fdata(i,j+1,k,n) + fdata(i,j,k+1,n)) / 16.0
									+
									fdata(i,j,k,n) / 8.0;
						}

				});
		}
	}

	// Copy the fine residual restricted onto the coarse grid
	// into the final residual.
	res.ParallelCopy(fine_res_for_coarse,0,0,ncomp,0,0,cgeom.periodicity());

	const int mglev = 0;

	// Sync up ghost nodes
	amrex::Geometry geom = m_geom[crse_amrlev][mglev];
	realFillBoundary(res,geom);
	nodalSync(crse_amrlev,mglev, res);
	return;
}

void
Operator<Grid::Node>::solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
			    const MultiFab* /*crse_bcdata*/)
{
	const int mglev = 0;
	const int ncomp = b.nComp();
	apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, StateMode::Solution);
	MultiFab::Xpay(resid, -1.0, b, 0, 0, ncomp, 2);
	amrex::Geometry geom = m_geom[amrlev][mglev];
	realFillBoundary(resid,geom);
}

void
Operator<Grid::Node>::correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
			      BCMode /*bc_mode*/, const MultiFab* /*crse_bcdata*/)
{
	resid.setVal(0.0);
	apply(amrlev, mglev, resid, x, BCMode::Homogeneous, StateMode::Correction);
	int ncomp = b.nComp();
	MultiFab::Xpay(resid, -1.0, b, 0, 0, ncomp, resid.nGrow());
	amrex::Geometry geom = m_geom[amrlev][mglev];
	realFillBoundary(resid,geom);
}




Operator<Grid::Cell>::Operator ()
{
	m_ixtype = amrex::IntVect::TheCellVector();
}

void
Operator<Grid::Cell>::define (amrex::Vector<amrex::Geometry> a_geom,
		   const amrex::Vector<amrex::BoxArray>& a_grids,
		   const amrex::Vector<amrex::DistributionMapping>& a_dmap,
		   BC::BC& a_bc,
		   const amrex::LPInfo& a_info,
		   const amrex::Vector<amrex::FabFactory<amrex::FArrayBox> const*>& a_factory)
{
	m_bc = &a_bc;

	std::array<int,AMREX_SPACEDIM> is_periodic = m_bc->IsPeriodic();
	
	MLCellLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory);

	Util::Warning(INFO,"This section of code has not been tested.");
	for (int n = 0; n < getNComp(); n++)
	{
		m_lobc.push_back( {AMREX_D_DECL(is_periodic[0] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet,
				       is_periodic[1] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet,
				       is_periodic[2] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet)});
		m_hibc.push_back( {AMREX_D_DECL(is_periodic[0] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet,
				       is_periodic[1] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet,
				       is_periodic[2] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet)});
	}

	for (int ilev = 0; ilev < a_geom.size(); ++ilev)
		setLevelBC(ilev,nullptr);

}


void
Operator<Grid::Cell>::prepareForSolve ()
{
	MLCellLinOp::prepareForSolve();
}

Operator<Grid::Cell>::BndryCondLoc::BndryCondLoc (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
	: bcond(ba, dm),
	  bcloc(ba, dm)
{
}

void
Operator<Grid::Cell>::BndryCondLoc::setLOBndryConds (const amrex::Geometry& /*geom*/, const amrex::Real* /*dx*/,
					  const amrex::Array<BCType,AMREX_SPACEDIM>& /*lobc*/,
					  const amrex::Array<BCType,AMREX_SPACEDIM>& /*hibc*/,
					  int /*ratio*/, const amrex::RealVect& /*a_loc*/)
{
//	const amrex::Box&  domain = geom.Domain();
	Util::Warning(INFO,"This code has not been properlyt tested");
//#ifdef _OPENMP
//#pragma omp parallel
//#endif
//	for (amrex::MFIter mfi(bcloc); mfi.isValid(); ++mfi)
//	{
//		const amrex::Box& bx = mfi.validbox();
//		RealTuple & bloc  = bcloc[mfi];
//		BCTuple   & bctag = bcond[mfi];
//		//amrex::MLMGBndry::setBoxBC(bloc, bctag, bx, domain, lobc, hibc, dx, ratio, a_loc,geom.isPeriodicArray());
//	}
}


void
Operator<Grid::Cell>::averageDownCoeffs ()
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
Operator<Grid::Cell>::averageDownCoeffsSameAmrLevel (amrex::Vector<amrex::MultiFab>& a)
{
	int nmglevs = a.size();
	for (int mglev = 1; mglev < nmglevs; ++mglev)
	{
		amrex::average_down(a[mglev-1], a[mglev], 0, a[0].nComp(), mg_coarsen_ratio);
	}
}



const amrex::FArrayBox &
Operator<Grid::Cell>::GetFab(const int num, const int amrlev, const int mglev, const amrex::MFIter &mfi) const
{
	return m_a_coeffs[num][amrlev][mglev][mfi];
}


void
Operator<Grid::Cell>::RegisterNewFab(amrex::Vector<amrex::MultiFab> &input)
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
Operator<Grid::Cell>::RegisterNewFab(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &input)
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

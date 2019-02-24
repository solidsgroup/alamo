
#include <AMReX_MLNodeLinOp.H>
#include <AMReX_MLCellLinOp.H>
#include <AMReX_MLLinOp_F.H>
#include <AMReX_MultiFabUtil.H>
#include "Util/Color.H"
#include "Set/Set.H"
#include "Operator.H"

using namespace amrex;
namespace Operator {

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

	for (MFIter mfi(x, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
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
	BL_PROFILE(Color::FG::Yellow + "Operator::Fsmooth()" + Color::Reset);

	amrex::Box domain(m_geom[amrlev][mglev].Domain());

	int ncomp = b.nComp();
	int nghost = b.nGrow();
	
	amrex::MultiFab Ax(x.boxArray(), x.DistributionMap(), ncomp, nghost);
	amrex::MultiFab Dx(x.boxArray(), x.DistributionMap(), ncomp, nghost);
	amrex::MultiFab Rx(x.boxArray(), x.DistributionMap(), ncomp, nghost);

	if (!m_diagonal_computed) Util::Abort(INFO,"Operator::Diagonal() must be called before using Fsmooth");

	Set::Scalar residual = 0.0;
	for (int redblack = 0; redblack < 2; redblack++)
	{
		Fapply(amrlev,mglev,Ax,x); // find Ax

		amrex::MultiFab::Copy(Dx,x,0,0,ncomp,2); // Dx = x
		amrex::MultiFab::Multiply(Dx,*m_diag[amrlev][mglev],0,0,ncomp,2); // Dx *= diag  (Dx = x*diag)

		amrex::MultiFab::Copy(Rx,Ax,0,0,ncomp,2); // Rx = Ax
		amrex::MultiFab::Subtract(Rx,Dx,0,0,ncomp,2); // Rx -= Dx  (Rx = Ax - Dx)

		for (MFIter mfi(x, true); mfi.isValid(); ++mfi)
		{
			const Box& bx = mfi.tilebox();
			amrex::FArrayBox       &xfab    = x[mfi];
			const amrex::FArrayBox &bfab    = b[mfi];
			const amrex::FArrayBox &Rxfab   = Rx[mfi];
			const amrex::FArrayBox &diagfab = (*m_diag[amrlev][mglev])[mfi];

			for (int n = 0; n < ncomp; n++)
			{
				AMREX_D_TERM(for (int m1 = bx.loVect()[0] - 2; m1<=bx.hiVect()[0] + 2; m1++),
					     for (int m2 = bx.loVect()[1] - 2; m2<=bx.hiVect()[1] + 2; m2++),
					     for (int m3 = bx.loVect()[2] - 2; m3<=bx.hiVect()[2] + 2; m3++))
				{

					if ((AMREX_D_TERM(m1, + m2, + m3))%2 == redblack) continue;
					amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));

					// Skip ghost cells outside problem domain
					if (AMREX_D_TERM(m[0] < domain.loVect()[0], ||
							 m[1] < domain.loVect()[1], ||
							 m[2] < domain.loVect()[2])) continue;
					if (AMREX_D_TERM(m[0] > domain.hiVect()[0] + 1, ||
							 m[1] > domain.hiVect()[1] + 1, ||
							 m[2] > domain.hiVect()[2] + 1)) continue;

					if (AMREX_D_TERM(m[0] == bx.loVect()[0] - 2 || m[0] == bx.hiVect()[0] + 2, ||
							 m[1] == bx.loVect()[1] - 2 || m[1] == bx.hiVect()[1] + 2, ||
							 m[2] == bx.loVect()[2] - 2 || m[2] == bx.hiVect()[2] + 2))
					{
						xfab(m,n) = 0.0;
						continue;
					}

					Set::Scalar xold = xfab(m,n);
					xfab(m,n) = (bfab(m,n) - Rxfab(m,n))/diagfab(m,n);
					residual += fabs(xold - xfab(m,n));

				}
			}
		}
	}
	nodalSync(amrlev, mglev, x);
}

void Operator<Grid::Node>::normalize (int amrlev, int mglev, MultiFab& x) const
{
	BL_PROFILE("Operator::normalize()");
	//Util::Message(INFO);
	bool debug = false;		// Enable this for inverse approximation
 
	amrex::Box domain(m_geom[amrlev][mglev].Domain());

	int ncomp = getNComp();
	int nghost = 1; //x.nGrow();

	//const Real* DX = m_geom[amrlev][mglev].CellSize();

	if (!m_diagonal_computed)
		Util::Abort(INFO,"Operator::Diagonal() must be called before using normalize");
	
	if(debug)
	{
		// We are trying to do a first order inverse correction here.
		amrex::MultiFab xtemp(x.boxArray(), x.DistributionMap(), ncomp, nghost);
		xtemp.setVal(0.0);
		amrex::MultiFab::Copy(xtemp,x,0,0,ncomp,0); // xtemp = x
		amrex::MultiFab R0x(x.boxArray(), x.DistributionMap(), ncomp, nghost);
		R0x.setVal(0,0);
		Error0x(amrlev,mglev,R0x,xtemp); 	// R0x = R0 * x = (I - A D0) * x
		amrex::MultiFab::Add(x,R0x,0,0,ncomp,0); // x = (I + R0)*x
	}

	for (MFIter mfi(x, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		amrex::FArrayBox       &xfab    = x[mfi];
		const amrex::FArrayBox &diagfab = (*m_diag[amrlev][mglev])[mfi];

		for (int n = 0; n < ncomp; n++)
		{
			AMREX_D_TERM(for (int m1 = bx.loVect()[0]-1; m1<=bx.hiVect()[0]+1; m1++),
				     for (int m2 = bx.loVect()[1]-1; m2<=bx.hiVect()[1]+1; m2++),
				     for (int m3 = bx.loVect()[2]-1; m3<=bx.hiVect()[2]+1; m3++))
			{
				amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));

				if (m[0] < domain.loVect()[0]) continue;
				if (m[1] < domain.loVect()[1]) continue;
				if (m[0] > domain.hiVect()[0]+1) continue;
				if (m[1] > domain.hiVect()[1]+1) continue;

				xfab(m,n) /= diagfab(m,n);
			}
		}
	}
}



bool Operator<Grid::Node>::VerificationCheck (int amrlev,
				  int mglev,
				  amrex::MultiFab& test) const
{
	BL_PROFILE("Operator::VerificationCheck()");
	Util::Message(INFO);
	bool result = false;
	int ncomp = test.nComp();
	int nghost = test.nGrow();
	amrex::MultiFab x(test.boxArray(), test.DistributionMap(), ncomp, nghost);
	amrex::MultiFab Ax(test.boxArray(), test.DistributionMap(), ncomp, nghost);

	x.setVal(0.0);
	Ax.setVal(0.0);
	test.setVal(0.0);

	for (MFIter mfi(x, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		amrex::FArrayBox	&testfab = test[mfi];
		amrex::FArrayBox	&xfab    = x[mfi];
		amrex::FArrayBox	&Axfab   = Ax[mfi];

		// Random point on inside
		int pointX = bx.loVect()[0] + rand() % (bx.hiVect()[0]-bx.loVect()[0]);
		int pointY = bx.loVect()[1] + rand() % (bx.hiVect()[1]-bx.loVect()[1]);
		int pointZ = bx.loVect()[2] + rand() % (bx.hiVect()[2]-bx.loVect()[2]);
		std::cout << "Point inside = (" << pointX << ", " << pointY << ", " << pointZ << ")" << std::endl;
		for (int i = 0; i < ncomp; i++)
		{
			amrex::IntVect m(AMREX_D_DECL(pointX,pointY,pointZ));
			xfab(m,i) = 1.0;
			Fapply(amrlev,mglev,Ax,x);
			testfab(m,i) = amrex::MultiFab::Dot(x,0,Ax,0,ncomp,nghost);
			std::cout << "test value = " << testfab(m,i) << std::endl;
			normalize(amrlev,mglev,x);
			std::cout << "Normalized = " << xfab(m,i) << ". Inverse = " << 1.0/xfab(m,i) << std::endl;
			xfab.setVal(0.0);
			Axfab.setVal(0.0);
		}

		// Random point on left face
		pointX = bx.loVect()[0];
		pointY = bx.loVect()[1] + rand() % (bx.hiVect()[1]-bx.loVect()[1]);
		pointZ = bx.loVect()[2] + rand() % (bx.hiVect()[2]-bx.loVect()[2]);
		std::cout << "Point left = (" << pointX << ", " << pointY << ", " << pointZ << ")" << std::endl;
		for (int i = 0; i < ncomp; i++)
		{
			amrex::IntVect m(AMREX_D_DECL(pointX,pointY,pointZ));
			xfab(m,i) = 1.0;
			Fapply(amrlev,mglev,Ax,x);
			testfab(m,i) = amrex::MultiFab::Dot(x,0,Ax,0,ncomp,nghost);
			std::cout << "test value = " << testfab(m,i) << std::endl;
			normalize(amrlev,mglev,x);
			std::cout << "Normalized = " << xfab(m,i) << ". Inverse = " << 1.0/xfab(m,i) << std::endl;
			xfab.setVal(0.0);
			Axfab.setVal(0.0);
		}

		// Random point on the right face
		pointX = bx.hiVect()[0];
		pointY = bx.loVect()[1] + rand() % (bx.hiVect()[1]-bx.loVect()[1]);
		pointZ = bx.loVect()[2] + rand() % (bx.hiVect()[2]-bx.loVect()[2]);
		std::cout << "Point right = (" << pointX << ", " << pointY << ", " << pointZ << ")" << std::endl;
		for (int i = 0; i < ncomp; i++)
		{
			amrex::IntVect m(AMREX_D_DECL(pointX,pointY,pointZ));
			xfab(m,i) = 1.0;
			Fapply(amrlev,mglev,Ax,x);
			testfab(m,i) = amrex::MultiFab::Dot(x,0,Ax,0,ncomp,nghost);
			std::cout << "test value = " << testfab(m,i) << std::endl;
			xfab.setVal(0.0);
			Axfab.setVal(0.0);
		}
		// Random point on bottom face
		pointY = bx.loVect()[1];
		pointX = bx.loVect()[0] + rand() % (bx.hiVect()[0]-bx.loVect()[0]);
		pointZ = bx.loVect()[2] + rand() % (bx.hiVect()[2]-bx.loVect()[2]);
		std::cout << "Point bottom = (" << pointX << ", " << pointY << ", " << pointZ << ")" << std::endl;
		for (int i = 0; i < ncomp; i++)
		{
			amrex::IntVect m(AMREX_D_DECL(pointX,pointY,pointZ));
			xfab(m,i) = 1.0;
			Fapply(amrlev,mglev,Ax,x);
			testfab(m,i) = amrex::MultiFab::Dot(x,0,Ax,0,ncomp,nghost);
			std::cout << "test value = " << testfab(m,i) << std::endl;
			normalize(amrlev,mglev,x);
			std::cout << "Normalized = " << xfab(m,i) << ". Inverse = " << 1.0/xfab(m,i) << std::endl;
			xfab.setVal(0.0);
			Axfab.setVal(0.0);
		}

		// Random point on the top face
		pointY = bx.hiVect()[1];
		pointX = bx.loVect()[0] + rand() % (bx.hiVect()[0]-bx.loVect()[0]);
		pointZ = bx.loVect()[2] + rand() % (bx.hiVect()[2]-bx.loVect()[2]);
		std::cout << "Point top = (" << pointX << ", " << pointY << ", " << pointZ << ")" << std::endl;
		for (int i = 0; i < ncomp; i++)
		{
			amrex::IntVect m(AMREX_D_DECL(pointX,pointY,pointZ));
			xfab(m,i) = 1.0;
			Fapply(amrlev,mglev,Ax,x);
			testfab(m,i) = amrex::MultiFab::Dot(x,0,Ax,0,ncomp,nghost);
			std::cout << "test value = " << testfab(m,i) << std::endl;
			normalize(amrlev,mglev,x);
			std::cout << "Normalized = " << xfab(m,i) << ". Inverse = " << 1.0/xfab(m,i) << std::endl;
			xfab.setVal(0.0);
			Axfab.setVal(0.0);
		}
		// Random point on back face
		pointZ = bx.loVect()[2];
		pointY = bx.loVect()[1] + rand() % (bx.hiVect()[1]-bx.loVect()[1]);
		pointX = bx.loVect()[0] + rand() % (bx.hiVect()[0]-bx.loVect()[0]);
		std::cout << "Point back = (" << pointX << ", " << pointY << ", " << pointZ << ")" << std::endl;
		for (int i = 0; i < ncomp; i++)
		{
			amrex::IntVect m(AMREX_D_DECL(pointX,pointY,pointZ));
			xfab(m,i) = 1.0;
			Fapply(amrlev,mglev,Ax,x);
			testfab(m,i) = amrex::MultiFab::Dot(x,0,Ax,0,ncomp,nghost);
			std::cout << "test value = " << testfab(m,i) << std::endl;
			normalize(amrlev,mglev,x);
			std::cout << "Normalized = " << xfab(m,i) << ". Inverse = " << 1.0/xfab(m,i) << std::endl;
			xfab.setVal(0.0);
			Axfab.setVal(0.0);
		}

		// Random point on the front face
		pointZ = bx.hiVect()[2];
		pointY = bx.loVect()[1] + rand() % (bx.hiVect()[1]-bx.loVect()[1]);
		pointX = bx.loVect()[0] + rand() % (bx.hiVect()[0]-bx.loVect()[0]);
		std::cout << "Point back = (" << pointX << ", " << pointY << ", " << pointZ << ")" << std::endl;
		for (int i = 0; i < ncomp; i++)
		{
			amrex::IntVect m(AMREX_D_DECL(pointX,pointY,pointZ));
			xfab(m,i) = 1.0;
			Fapply(amrlev,mglev,Ax,x);
			testfab(m,i) = amrex::MultiFab::Dot(x,0,Ax,0,ncomp,nghost);
			std::cout << "test value = " << testfab(m,i) << std::endl;
			normalize(amrlev,mglev,x);
			std::cout << "Normalized = " << xfab(m,i) << ". Inverse = " << 1.0/xfab(m,i) << std::endl;
			xfab.setVal(0.0);
			Axfab.setVal(0.0);
		}
	}
	return result;
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

	 // This makes sure grids are cell-centered;
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
 }


void Operator<Grid::Node>::compRHS (const Vector<MultiFab*>& /*rhs*/, const Vector<MultiFab*>& /*vel*/,
		   const Vector<const MultiFab*>& /*rhnd*/,
		   const Vector<MultiFab*>& /*a_rhcc*/)
{
	//Util::Message(INFO);
	Util::Abort(INFO, "compRHS not implemented");
}



void Operator<Grid::Node>::buildMasks ()
{
	BL_PROFILE("Operator::buildMasks()");
	if (m_masks_built) return;


	m_masks_built = true;

	m_is_bottom_singular = false;
	// auto itlo = std::find(m_lobc.begin(), m_lobc.end(), BCType::Dirichlet);
	// auto ithi = std::find(m_hibc.begin(), m_hibc.end(), BCType::Dirichlet);

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
                        
						// amrex_mlndlap_set_dirichlet_mask(BL_TO_FORTRAN_ANYD(mskfab),
						// 				 BL_TO_FORTRAN_ANYD(ccfab),
						// 				 BL_TO_FORTRAN_BOX(nddomain),
						// 				 m_lobc.data(), m_hibc.data());
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
}


void Operator<Grid::Node>::fixUpResidualMask (int amrlev, iMultiFab& resmsk)
{
	BL_PROFILE("Operator::fixUpResidualMask()");

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

	//Util::Message(INFO, "Not implemented (and shouldn't need to be!)");
}

void Operator<Grid::Node>::prepareForSolve ()
{
	BL_PROFILE("Operator::prepareForSolve()");

	AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_num_amr_levels == 1 ||
					 m_coarsening_strategy != CoarseningStrategy::RAP,
					 "Operator::prepareForSolve RAP TODO");

	MLNodeLinOp::prepareForSolve();

	buildMasks();


	averageDownCoeffs();
	Diagonal(true);

}

void Operator<Grid::Node>::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
	BL_PROFILE("Operator::restriction()");

	applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous, StateMode::Solution);

	//const Box& nd_domain = amrex::surroundingNodes(m_geom[amrlev][cmglev].Domain());

	bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
	MultiFab cfine;
	if (need_parallel_copy) {
		const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
		cfine.define(ba, fine.DistributionMap(), fine.nComp(), 0);
	}

	MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;

	// C++ version of Fortran (amrex_mlndlap_restriction)
	// Set::Scalar fac1 = 1.0/64.0;
	// Set::Scalar fac2 = 1.0/32.0;
	// Set::Scalar fac3 = 1.0/16.0; 
	// Set::Scalar fac4 = 1.0/8.0;

	static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
					   dy(AMREX_D_DECL(0,1,0)),
					   dz(AMREX_D_DECL(0,0,1)));

	for (MFIter mfi(*pcrse, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		const amrex::FArrayBox &finefab = fine[mfi];
		amrex::FArrayBox       &crsefab = (*pcrse)[mfi];
		AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
			     for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
			     for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
		{
			amrex::IntVect m_crse(AMREX_D_DECL(m1,m2,m3));
			amrex::IntVect m_fine(AMREX_D_DECL(2*m1, 2*m2, 2*m3));
			for (int i=0; i<crse.nComp(); i++)
			{
#if AMREX_SPACEDIM == 2
		
				// // amrex_mlndlap_restriction
				// for (int n = 0 ; n < ncomp; n++)
				// {
				// 	for (int m2 = bx.loVect()[1] +1; m2<=bx.hiVect()[1] -1; m2++)
				// 		for (int m1 = bx.loVect()[0] +1; m1<=bx.hiVect()[0] -1; m1++)
				// 	{
				// 		amrex::IntVect m_crse(AMREX_D_DECL(m1,m2,m3));
				// 		amrex::IntVect m_fine(AMREX_D_DECL(m1*2,m2*2,m3*2));

				crsefab(m_crse,i) =
					(+     finefab(m_fine-dx[0]-dx[1],i) + 2.0*finefab(m_fine-dx[1],i) +     finefab(m_fine+dx[0]-dx[1],i)
					 + 2.0*finefab(m_fine-dx[0]      ,i) + 4.0*finefab(m_fine      ,i) + 2.0*finefab(m_fine+dx[0]      ,i) 
					 +     finefab(m_fine-dx[0]+dx[1],i) + 2.0*finefab(m_fine+dx[1],i) +     finefab(m_fine+dx[0]+dx[1],i))/16.0;
				//}
#endif
#if AMREX_SPACEDIM == 3
				crsefab(m_crse,i) =
					(finefab(m_fine-dx-dy-dz,i) +
					 finefab(m_fine-dx-dy+dz,i) +
					 finefab(m_fine-dx+dy-dz,i) +
					 finefab(m_fine-dx+dy+dz,i) +
					 finefab(m_fine+dx-dy-dz,i) +
					 finefab(m_fine+dx-dy+dz,i) +
					 finefab(m_fine+dx+dy-dz,i) +
					 finefab(m_fine+dx+dy+dz,i)) / 64.0
					+
					(finefab(m_fine-dy-dz,i) +
					 finefab(m_fine-dy+dz,i) +
					 finefab(m_fine+dy-dz,i) +
					 finefab(m_fine+dy+dz,i) +
					 finefab(m_fine-dz-dx,i) +
					 finefab(m_fine-dz+dx,i) +
					 finefab(m_fine+dz-dx,i) +
					 finefab(m_fine+dz+dx,i) +
					 finefab(m_fine-dx-dy,i) +
					 finefab(m_fine-dx+dy,i) +
					 finefab(m_fine+dx-dy,i) +
					 finefab(m_fine+dx+dy,i)) / 32.0
					+
					(finefab(m_fine-dx,i) +
					 finefab(m_fine-dy,i) +
					 finefab(m_fine-dz,i) +
					 finefab(m_fine+dx,i) +
					 finefab(m_fine+dy,i) +
					 finefab(m_fine+dz,i)) / 16.0
					+
					finefab(m_fine,i) / 8.0;
#endif
			}
		}
	}

	if (need_parallel_copy) {
		crse.ParallelCopy(cfine);
	}

	// if (fine.contains_nan() || fine.contains_inf()) Util::Abort(INFO, "restriction (end) - nan or inf detected in fine");
	// if (crse.contains_nan() || crse.contains_inf()) Util::Abort(INFO, "restriction (end) - nan or inf detected in crse");
}

void Operator<Grid::Node>::interpolation (int /*amrlev*/, int /*fmglev*/, MultiFab& fine, const MultiFab& crse) const
{
	BL_PROFILE("Operator::interpolation()");
	// if (fine.contains_nan() || fine.contains_inf()) Util::Abort(INFO, "interpolation (beginning) - nan or inf detected in fine");
	// if (crse.contains_nan() || crse.contains_inf()) Util::Abort(INFO, "interpolation (beginning) - nan or inf detected in crse");
	bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
	MultiFab cfine;
	const MultiFab* cmf = &crse;
	if (need_parallel_copy) {
		const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
		cfine.define(ba, fine.DistributionMap(), crse.nComp(), 0);
		cfine.ParallelCopy(crse);
		cmf = &cfine;
	}
	static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
					   dy(AMREX_D_DECL(0,1,0)),
					   dz(AMREX_D_DECL(0,0,1)));
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		for (MFIter mfi(fine, true); mfi.isValid(); ++mfi)
		{
			const Box& fine_bx = mfi.tilebox();
			const Box& course_bx = amrex::coarsen(fine_bx,2);
			const Box& tmpbx = amrex::refine(course_bx,2);
			FArrayBox tmpfab;
			tmpfab.resize(tmpbx,fine.nComp());
			tmpfab.setVal(0.0);
			
			const amrex::FArrayBox &crsefab = (*cmf)[mfi];
			
			for (int i=0; i<crse.nComp(); i++)
			{
				AMREX_D_TERM(for (int m1 = fine_bx.loVect()[0]; m1<=fine_bx.hiVect()[0]; m1++),
					     for (int m2 = fine_bx.loVect()[1]; m2<=fine_bx.hiVect()[1]; m2++),
					     for (int m3 = fine_bx.loVect()[2]; m3<=fine_bx.hiVect()[2]; m3++))
				{
					amrex::IntVect m(AMREX_D_DECL(m1, m2, m3));
					amrex::IntVect M(AMREX_D_DECL(m1/2, m2/2, m3/2));

#if AMREX_SPACEDIM == 2
					if (m[0]==2*M[0] && m[1]==2*M[1]) // Coincident
						tmpfab(m,i) = crsefab(M,i);
					else if (m[1]==2*M[1]) // X Edge
						tmpfab(m,i) = 0.5 * (crsefab(M,i) + crsefab(M+dx,i));
					else if (m[0]==2*M[0]) // Y Edge
						tmpfab(m,i) = 0.5 * (crsefab(M,i) + crsefab(M+dy,i));
					else // Center
						tmpfab(m,i) = 0.25 * (crsefab(M,i) + crsefab(M+dx,i) +
								      crsefab(M+dy,i) + crsefab(M+dx+dy,i));
#endif
#if AMREX_SPACEDIM == 3
					if (m[0]==2*M[0] && m[1]==2*M[1] && m[2]==2*M[2]) // Coincident
						tmpfab(m,i) = crsefab(M,i);
					else if (m[1]==2*M[1] && m[2]==2*M[2]) // X Edge
						tmpfab(m,i) = 0.5 * (crsefab(M,i) + crsefab(M+dx,i));
					else if (m[2]==2*M[2] && m[0]==2*M[0]) // Y Edge
						tmpfab(m,i) = 0.5 * (crsefab(M,i) + crsefab(M+dy,i));
					else if (m[0]==2*M[0] && m[1]==2*M[1]) // Z Edge
						tmpfab(m,i) = 0.5 * (crsefab(M,i) + crsefab(M+dz,i));
					else if (m[0]==2*M[0]) // X Face
						tmpfab(m,i) = 0.25 * (crsefab(M,i) + crsefab(M+dy,i) +
								      crsefab(M+dz,i) + crsefab(M+dy+dz,i));
					else if (m[1]==2*M[1]) // Y Face
						tmpfab(m,i) = 0.25 * (crsefab(M,i) + crsefab(M+dz,i) +
								      crsefab(M+dx,i) + crsefab(M+dz+dx,i));
					else if (m[2]==2*M[2]) // Z Face
						tmpfab(m,i) = 0.25 * (crsefab(M,i) + crsefab(M+dx,i) +
								      crsefab(M+dy,i) + crsefab(M+dx+dy,i));
					else // Center
					{
						tmpfab(m,i) = 0.125 * (crsefab(M,i) +
								       crsefab(M+dx,i) + crsefab(M+dy,i) + crsefab(M+dz,i) +
								       crsefab(M+dy+dz,i) + crsefab(M+dz+dx,i) + crsefab(M+dx+dy,i) +
								       crsefab(M+dx+dy+dz,i));
					}

					if (std::isinf(tmpfab(m,i)))
					{
						std::cout << m << M << std::endl;
						Util::Abort(INFO, "Is Infinity");
					}
					if (std::isnan(tmpfab(m,i)))
					{
						std::cout << m << M << std::endl;
						Util::Abort(INFO, "Is NaN");
					}
#endif
				}
			}
			fine[mfi].plus(tmpfab,fine_bx,fine_bx,0,0,fine.nComp());
		}
	}
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

void Operator<Grid::Node>::applyBC (int amrlev, int mglev, MultiFab& phi, BCMode/* bc_mode*/,
		   amrex::MLLinOp::StateMode /**/, bool skip_fillboundary) const
{
	BL_PROFILE("Operator::applyBC()");

	const Geometry& geom = m_geom[amrlev][mglev];

	if (!skip_fillboundary) {

		//
		// This is special code to fill a boundary when there are TWO ghost nodes
		//

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

	const Geometry& cgeom = m_geom[crse_amrlev  ][0];
 	//const Geometry& fgeom = m_geom[crse_amrlev+1][0];
 	const Box& c_cc_domain = cgeom.Domain();
 	//const Box& c_nd_domain = amrex::surroundingNodes(c_cc_domain);
	//const Real* cDX = cgeom.CellSize();
	//const Real* fDX = fgeom.CellSize();
	

 	const BoxArray&            fba = fine_res.boxArray();
 	const DistributionMapping& fdm = fine_res.DistributionMap();

 	MultiFab fine_res_for_coarse(amrex::coarsen(fba, 2), fdm, ncomp, 2);
	fine_res_for_coarse.ParallelCopy(res,0,0,ncomp,0,0,cgeom.periodicity());

 	applyBC(crse_amrlev+1, 0, fine_res, BCMode::Inhomogeneous, StateMode::Solution);

	const int coarse_coarse_node = 0;
	const int coarse_fine_node = 1;
	const int fine_fine_node = 2;

	amrex::iMultiFab nodemask(amrex::coarsen(fba,2), fdm, 1, 2);
	nodemask.ParallelCopy(*m_nd_fine_mask[crse_amrlev],0,0,1,0,0,cgeom.periodicity());

	amrex::iMultiFab cellmask(amrex::convert(amrex::coarsen(fba,2),amrex::IntVect::TheCellVector()), fdm, 1, 2);
	cellmask.ParallelCopy(*m_cc_fine_mask[crse_amrlev],0,0,1,1,1,cgeom.periodicity());
	


	for (MFIter mfi(fine_res_for_coarse, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		FArrayBox &fine = fine_res[mfi];
		FArrayBox &crse = fine_res_for_coarse[mfi];
		//const FArrayBox &crserhs = crse_rhs[mfi];
		//const FArrayBox &finerhs = fine_rhs[mfi];
		//const FArrayBox &ufab = fine_sol[mfi];
		//TArrayBox &C = (*(model[crse_amrlev+1][0]))[mfi];
		
#if AMREX_SPACEDIM == 1
		Util::Abort(INFO, "reflux not implemented in 1D. Turn AMR off or switch to 2D.");
#elif AMREX_SPACEDIM == 2
		for (int n = 0 ; n < ncomp; n++)
			for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++)
				for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++)
				{
					amrex::IntVect m_crse(m1,  m2);
					amrex::IntVect m_fine(m1*2,m2*2);
					
					// Domain boundaries
					if (m1 == c_cc_domain.loVect()[0] || m1 == c_cc_domain.hiVect()[0] + 1 ||
					    m2 == c_cc_domain.loVect()[1] || m2 == c_cc_domain.hiVect()[1] + 1)
					{
						// Domain corner
						if ((m1 == c_cc_domain.loVect()[0]     && m2 == c_cc_domain.loVect()[1]     ) ||
						    (m1 == c_cc_domain.loVect()[0]     && m2 == c_cc_domain.hiVect()[1] + 1 ) ||
						    (m1 == c_cc_domain.hiVect()[0] + 1 && m2 == c_cc_domain.loVect()[1]     ) ||
						    (m1 == c_cc_domain.hiVect()[0] + 1 && m2 == c_cc_domain.hiVect()[1] + 1 ))
							crse(m_crse,n) = fine(m_fine,n);
						// Domain xlo or xhi edge
						else if (m1 == c_cc_domain.loVect()[0] || m1 == c_cc_domain.hiVect()[0] +1)
							crse(m_crse,n) = 0.25*fine(m_fine-dx[1],n) + 0.5*fine(m_fine,n) + 0.25*fine(m_fine+dx[1],n);
						// Domain ylo or yhi edge
						else if (m2 == c_cc_domain.loVect()[1] || m2 == c_cc_domain.hiVect()[1] +1)
							crse(m_crse,n) = 0.25*fine(m_fine-dx[0],n) + 0.5*fine(m_fine,n) + 0.25*fine(m_fine+dx[0],n);
						continue;
					}
					// Interior boundaries
					if ((nodemask[mfi])(m_crse) == fine_fine_node ||
					    (nodemask[mfi])(m_crse) == coarse_fine_node )
					{
						// if (m_crse == amrex::IntVect(7,16) ||
						//     m_crse == amrex::IntVect(8,16) ||
						//     m_crse == amrex::IntVect(24,16) ||
						//     m_crse == amrex::IntVect(25,16) 
						//     ) {
						// 	Util::Message(INFO,"m = ", m_crse, " val = ", fine(m_fine      ,n));
						// 	crse(m_crse,n) = 0.0;
						// 	continue;}
						crse(m_crse,n) = 
						 	((+     fine(m_fine-dx[0]-dx[1],n) + 2.0*fine(m_fine-dx[1],n) +     fine(m_fine+dx[0]-dx[1],n)
						  	  + 2.0*fine(m_fine-dx[0]      ,n) + 4.0*fine(m_fine      ,n) + 2.0*fine(m_fine+dx[0]      ,n) 
						  	  +     fine(m_fine-dx[0]+dx[1],n) + 2.0*fine(m_fine+dx[1],n) +     fine(m_fine+dx[0]+dx[1],n))/16.0);
					}
				}
#elif AMREX_SPACEDIM == 3
		for (int n = 0 ; n < ncomp; n++)
			for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++)
				for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++)
					for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++)
					{
						amrex::IntVect m_crse(m1,  m2,  m3  );
						amrex::IntVect m_fine(m1*2,m2*2,m3*2);
				
						bool xmin = (m1 == bx.loVect()[0]) || (m1 == c_cc_domain.loVect()[0]);
						bool xmax = (m1 == bx.hiVect()[0]) || (m1 == c_cc_domain.hiVect()[0] +1);
						bool ymin = (m2 == bx.loVect()[1]) || (m2 == c_cc_domain.loVect()[1]);
						bool ymax = (m2 == bx.hiVect()[1]) || (m2 == c_cc_domain.hiVect()[1] +1);
						bool zmin = (m3 == bx.loVect()[2]) || (m3 == c_cc_domain.loVect()[2]);
						bool zmax = (m3 == bx.hiVect()[2]) || (m3 == c_cc_domain.hiVect()[2] +1);
						
						if ((nodemask[mfi])(m_crse) == fine_fine_node ||
						    (nodemask[mfi])(m_crse) == coarse_fine_node)
						{
							crse(m_crse,n) =
								0.25*
								((+ 0.25 *fine(m_fine-dx[0]-dx[1]-dx[2],n) + 0.5 *fine(m_fine-dx[1]-dx[2],n) + 0.25 *fine(m_fine+dx[0]-dx[1]-dx[2],n)
								  + 0.5  *fine(m_fine-dx[0]      -dx[2],n) + 1.0 *fine(m_fine      -dx[2],n) + 0.5  *fine(m_fine+dx[0]      -dx[2],n) 
								  + 0.25 *fine(m_fine-dx[0]+dx[1]-dx[2],n) + 0.5 *fine(m_fine+dx[1]-dx[2],n) + 0.25 *fine(m_fine+dx[0]+dx[1]-dx[2],n)) / 4.0)
								+
								0.5*
								((+ 0.25 *fine(m_fine-dx[0]-dx[1]      ,n) + 0.5 *fine(m_fine-dx[1]      ,n) + 0.25 *fine(m_fine+dx[0]-dx[1]      ,n)
								  + 0.5  *fine(m_fine-dx[0]            ,n) + 1.0 *fine(m_fine            ,n) + 0.5  *fine(m_fine+dx[0]            ,n) 
								  + 0.25 *fine(m_fine-dx[0]+dx[1]      ,n) + 0.5 *fine(m_fine+dx[1]      ,n) + 0.25 *fine(m_fine+dx[0]+dx[1]      ,n)) / 4.0)
								+
								0.25*
								((+ 0.25 *fine(m_fine-dx[0]-dx[1]+dx[2],n) + 0.5 *fine(m_fine-dx[1]+dx[2],n) + 0.25 *fine(m_fine+dx[0]-dx[1]+dx[2],n)
								  + 0.5  *fine(m_fine-dx[0]      +dx[2],n) + 1.0 *fine(m_fine      +dx[2],n) + 0.5  *fine(m_fine+dx[0]      +dx[2],n) 
								  + 0.25 *fine(m_fine-dx[0]+dx[1]+dx[2],n) + 0.5 *fine(m_fine+dx[1]+dx[2],n) + 0.25 *fine(m_fine+dx[0]+dx[1]+dx[2],n)) / 4.0);
						}
						else if ((nodemask[mfi])(m_crse) == coarse_fine_node)
						{

							Util::Abort(INFO,"you should not be here!");

							// Corners
							if ((xmin && ymin && zmin) ||
							    (xmin && ymin && zmax) ||
							    (xmin && ymax && zmin) ||
							    (xmin && ymax && zmax) ||
							    (xmax && ymin && zmin) ||
							    (xmax && ymin && zmax) ||
							    (xmax && ymax && zmin) ||
							    (xmax && ymax && zmax))
								crse(m_crse,n) = fine(m_fine,n);
							// X edges
							else if ( (ymin && zmin) || (ymin && zmax) || (ymax && zmin) || (ymax && zmax) )
								crse(m_crse,n) = (0.25*fine(m_fine-dx[0],n) + 0.5*fine(m_fine,n) + 0.25*fine(m_fine+dx[0],n));
							// Y edges
							else if ( (zmin && zmin) || (zmin && xmax) || (zmax && xmin) || (zmax && xmax) )
								crse(m_crse,n) = (0.25*fine(m_fine-dx[1],n) + 0.5*fine(m_fine,n) + 0.25*fine(m_fine+dx[1],n));
							// Z edges
							else if ( (xmin && ymin) || (xmin && ymax) || (xmax && ymin) || (xmax && ymax) )
								crse(m_crse,n) = (0.25*fine(m_fine-dx[2],n) + 0.5*fine(m_fine,n) + 0.25*fine(m_fine+dx[2],n));
							// YZ face (X=const)
							else if ( xmin || xmax )
								crse(m_crse,n) =
									((+     fine(m_fine-dx[1]-dx[2],n) + 2.0*fine(m_fine-dx[2],n) +     fine(m_fine+dx[1]-dx[2],n)
									  + 2.0*fine(m_fine-dx[1]      ,n) + 4.0*fine(m_fine      ,n) + 2.0*fine(m_fine+dx[1]      ,n) 
									  +     fine(m_fine-dx[1]+dx[2],n) + 2.0*fine(m_fine+dx[2],n) +     fine(m_fine+dx[1]+dx[2],n))/16.0);
							// ZX face (Y=const)
							else if ( xmin || xmax )
								crse(m_crse,n) =
									((+     fine(m_fine-dx[2]-dx[0],n) + 2.0*fine(m_fine-dx[0],n) +     fine(m_fine+dx[2]-dx[2],n)
									  + 2.0*fine(m_fine-dx[2]      ,n) + 4.0*fine(m_fine      ,n) + 2.0*fine(m_fine+dx[2]      ,n) 
									  +     fine(m_fine-dx[2]+dx[0],n) + 2.0*fine(m_fine+dx[0],n) +     fine(m_fine+dx[2]+dx[2],n))/16.0);
							// XY face (Z=const)
							else if ( xmin || xmax )
								crse(m_crse,n) =
									((+     fine(m_fine-dx[0]-dx[1],n) + 2.0*fine(m_fine-dx[1],n) +     fine(m_fine+dx[0]-dx[1],n)
									  + 2.0*fine(m_fine-dx[0]      ,n) + 4.0*fine(m_fine      ,n) + 2.0*fine(m_fine+dx[0]      ,n) 
									  +     fine(m_fine-dx[0]+dx[1],n) + 2.0*fine(m_fine+dx[1],n) +     fine(m_fine+dx[0]+dx[1],n))/16.0);
							// Internal
							else
							{
								crse(m_crse,n) =
									0.25*
									((+ 0.25 *fine(m_fine-dx[0]-dx[1]-dx[2],n) + 0.5 *fine(m_fine-dx[1]-dx[2],n) + 0.25 *fine(m_fine+dx[0]-dx[1]-dx[2],n)
									  + 0.5  *fine(m_fine-dx[0]      -dx[2],n) + 1.0 *fine(m_fine      -dx[2],n) + 0.5  *fine(m_fine+dx[0]      -dx[2],n) 
									  + 0.25 *fine(m_fine-dx[0]+dx[1]-dx[2],n) + 0.5 *fine(m_fine+dx[1]-dx[2],n) + 0.25 *fine(m_fine+dx[0]+dx[1]-dx[2],n)) / 4.0)
									+
									0.5*
									((+ 0.25 *fine(m_fine-dx[0]-dx[1]      ,n) + 0.5 *fine(m_fine-dx[1]      ,n) + 0.25 *fine(m_fine+dx[0]-dx[1]      ,n)
									  + 0.5  *fine(m_fine-dx[0]            ,n) + 1.0 *fine(m_fine            ,n) + 0.5  *fine(m_fine+dx[0]            ,n) 
									  + 0.25 *fine(m_fine-dx[0]+dx[1]      ,n) + 0.5 *fine(m_fine+dx[1]      ,n) + 0.25 *fine(m_fine+dx[0]+dx[1]      ,n)) / 4.0)
									+
									0.25*
									((+ 0.25 *fine(m_fine-dx[0]-dx[1]+dx[2],n) + 0.5 *fine(m_fine-dx[1]+dx[2],n) + 0.25 *fine(m_fine+dx[0]-dx[1]+dx[2],n)
									  + 0.5  *fine(m_fine-dx[0]      +dx[2],n) + 1.0 *fine(m_fine      +dx[2],n) + 0.5  *fine(m_fine+dx[0]      +dx[2],n) 
									  + 0.25 *fine(m_fine-dx[0]+dx[1]+dx[2],n) + 0.5 *fine(m_fine+dx[1]+dx[2],n) + 0.25 *fine(m_fine+dx[0]+dx[1]+dx[2],n)) / 4.0);
							}
						}
						else if ((nodemask[mfi])(m_crse) == coarse_coarse_node)
						{
							//Util::Message(INFO,"Discovered coarse node while on fine fab: crse amrlev = ", crse_amrlev,", m_crse = ", m_crse, " box = ", bx);
						}

					}
#endif		
	}

	// Copy the fine residual restricted onto the coarse grid
	// into the final residual.
	res.ParallelCopy(fine_res_for_coarse,0,0,ncomp,0,0,cgeom.periodicity());

	const int mglev = 0;

	// Sync up ghost nodes
	amrex::Geometry geom = m_geom[crse_amrlev][mglev];
	for (int i = 0; i < 2; i++)
	{
		MultiFab & mf = res;
		mf.FillBoundary(geom.periodicity());
		const int ncomp = mf.nComp();
		const int ng1 = 1;
		const int ng2 = 2;
		MultiFab tmpmf(mf.boxArray(), mf.DistributionMap(), ncomp, ng1);
		MultiFab::Copy(tmpmf, mf, 0, 0, ncomp, ng1); 
		mf.ParallelCopy   (tmpmf, 0, 0, ncomp, ng1, ng2, geom.periodicity());
	}

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
}

void
Operator<Grid::Node>::correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
			      BCMode /*bc_mode*/, const MultiFab* /*crse_bcdata*/)
{
    apply(amrlev, mglev, resid, x, BCMode::Homogeneous, StateMode::Correction);
    int ncomp = b.nComp();
    MultiFab::Xpay(resid, -1.0, b, 0, 0, ncomp, resid.nGrow());
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
	// for (int ilev=0; ilev < a_geom.size(); ilev++)
	// 	a_geom[ilev].SetPeriodicity(is_periodic);
	
	MLCellLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory);

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
Operator<Grid::Cell>::prepareForSolve ()
{
	BL_PROFILE("Operator<Grid::Cell>::prepareForSolve()");

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

Operator<Grid::Cell>::BndryCondLoc::BndryCondLoc (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
	: bcond(ba, dm),
	  bcloc(ba, dm)
{
}

void
Operator<Grid::Cell>::BndryCondLoc::setLOBndryConds (const amrex::Geometry& geom, const amrex::Real* dx,
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

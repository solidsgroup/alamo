
#include "Operator.H"
#include <AMReX_MLLinOp_F.H>
#include <AMReX_MultiFabUtil.H>
#include "Util/Color.H"
#include "Set/Set.H"

using namespace amrex;
namespace Operator {


void Operator::Diagonal (bool recompute)
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

void Operator::Diagonal (int amrlev, int mglev, amrex::MultiFab &diag)
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



void Operator::Fsmooth(int amrlev, int mglev, amrex::MultiFab& x, const amrex::MultiFab& b) const
{
	BL_PROFILE(Color::FG::Yellow + "Operator::Fsmooth()" + Color::Reset);
	//Util::Message(INFO);

	int ncomp = b.nComp();
	int nghost = b.nGrow();
	//Util::Message(INFO,"ncomp = ", b.nComp(), " ngrow = ",b.nGrow());
	
	//amrex::MultiFab diag(x.boxArray(), x.DistributionMap(), ncomp, nghost);
	amrex::MultiFab Ax(x.boxArray(), x.DistributionMap(), ncomp, nghost);
	amrex::MultiFab Dx(x.boxArray(), x.DistributionMap(), ncomp, nghost);
	amrex::MultiFab Rx(x.boxArray(), x.DistributionMap(), ncomp, nghost);

	if (!m_diagonal_computed) Util::Abort(INFO,"Operator::Diagonal() must be called before using Fsmooth");

	// if ((b.contains_nan() || b.contains_inf()) && !NanInfOnBoundary(b))
	// 	Util::Abort(INFO,"amrlev = ",amrlev," mglev = ",mglev, " b fab contains nan/inf on interior!");

	Set::Scalar residual = 0.0;
	for (int redblack = 0; redblack < 2; redblack++)
	{
		Fapply(amrlev,mglev,Ax,x); // find Ax

		amrex::MultiFab::Copy(Dx,x,0,0,ncomp,0); // Dx = x
		amrex::MultiFab::Multiply(Dx,*m_diag[amrlev][mglev],0,0,ncomp,0); // Dx *= diag  (Dx = x*diag)

		amrex::MultiFab::Copy(Rx,Ax,0,0,ncomp,0); // Rx = Ax
		amrex::MultiFab::Subtract(Rx,Dx,0,0,ncomp,0); // Rx -= Dx  (Rx = Ax - Dx)

		for (MFIter mfi(x, true); mfi.isValid(); ++mfi)
		{
			const Box& bx = mfi.tilebox();
			amrex::FArrayBox       &xfab    = x[mfi];
			const amrex::FArrayBox &bfab    = b[mfi];
			const amrex::FArrayBox &Rxfab   = Rx[mfi];
			const amrex::FArrayBox &diagfab = (*m_diag[amrlev][mglev])[mfi];

			for (int n = 0; n < ncomp; n++)
			{
				AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
					     for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
					     for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
				{

					if ((AMREX_D_TERM(m1, + m2, + m3))%2 == redblack) continue;
					amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));


					Set::Scalar xold = xfab(m,n);
					xfab(m,n) = (bfab(m,n) - Rxfab(m,n))/diagfab(m,n);
					residual += fabs(xold - xfab(m,n));
					// Util::Message(INFO,m);
					// Util::Message(INFO,xfab(m,n));
					// Util::Message(INFO,bfab(m,n));
					// Util::Message(INFO,Rxfab(m,n));
					// Util::Message(INFO,diagfab(m,n));

				}
			}
		}
	}
	//Util::Message(INFO,"residual = ", residual);
}

void Operator::normalize (int amrlev, int mglev, MultiFab& x) const
{
	BL_PROFILE("Operator::normalize()");
	//Util::Message(INFO);
	bool debug = false;		// Enable this for inverse approximation

	int ncomp = getNComp();
	int nghost = 1; //x.nGrow();

	amrex::Box domain(m_geom[amrlev][mglev].Domain());
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
			AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
				     for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
				     for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
			{
				amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
				xfab(m,n) /= diagfab(m,n);
			}
		}
	}
}



bool Operator::VerificationCheck (int amrlev,
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


Operator::Operator (const Vector<Geometry>& a_geom,
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

 Operator::~Operator ()
 {}

 void
 Operator::define (const Vector<Geometry>& a_geom,
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

	 int nghost = 0;
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


void
Operator::compRHS (const Vector<MultiFab*>& /*rhs*/, const Vector<MultiFab*>& /*vel*/,
		   const Vector<const MultiFab*>& /*rhnd*/,
		   const Vector<MultiFab*>& /*a_rhcc*/)
{
	//Util::Message(INFO);
	Util::Abort(INFO, "compRHS not implemented");
}



void
Operator::buildMasks ()
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


void
Operator::fixUpResidualMask (int /*amrlev*/, iMultiFab& /*resmsk*/)
{
	BL_PROFILE("Operator::fixUpResidualMask()");
	Util::Message(INFO, "Not implemented (and shouldn't need to be!)");
}

void
Operator::prepareForSolve ()
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

void
Operator::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
	BL_PROFILE("Operator::restriction()");
	Util::Message(INFO);

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
	amrex::Box domain(m_geom[amrlev][cmglev].Domain());

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

void
Operator::interpolation (int /*amrlev*/, int /*fmglev*/, MultiFab& fine, const MultiFab& crse) const
{
	BL_PROFILE("Operator::interpolation()");
	Util::Message(INFO,"Not yet implemented");
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
				AMREX_D_TERM(for (int m1 = fine_bx.loVect()[0]-1; m1<=fine_bx.hiVect()[0]+1; m1++),
					     for (int m2 = fine_bx.loVect()[1]-1; m2<=fine_bx.hiVect()[1]+1; m2++),
					     for (int m3 = fine_bx.loVect()[2]-1; m3<=fine_bx.hiVect()[2]+1; m3++))
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

void
Operator::averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& /*crse_rhs*/,
				  const MultiFab& fine_sol, const MultiFab& fine_rhs)
{
	BL_PROFILE("Operator::averageDownSolutionRHS()");
	Util::Message(INFO,"Suspect implementation!");
	const auto& amrrr = AMRRefRatio(camrlev);
	amrex::average_down(fine_sol, crse_sol, 0, fine_rhs.nComp(), amrrr);
}

void
Operator::applyBC (int amrlev, int mglev, MultiFab& phi, BCMode/* bc_mode*/,
		   amrex::MLLinOp::StateMode /**/, bool skip_fillboundary) const
{
	BL_PROFILE("Operator::applyBC()");

	const Geometry& geom = m_geom[amrlev][mglev];
	const Box& nd_domain = amrex::surroundingNodes(geom.Domain());

	if (!skip_fillboundary) {
	 	phi.FillBoundary(geom.periodicity());
	}

	if (m_coarsening_strategy == CoarseningStrategy::Sigma)
	{
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

const amrex::FArrayBox &
Operator::GetFab(const int num, const int amrlev, const int mglev, const amrex::MFIter &mfi) const
{
	BL_PROFILE("Operator::GetFab()");
	Util::Message(INFO);
 	return m_a_coeffs[num][amrlev][mglev][mfi];
}

void
Operator::RegisterNewFab(amrex::Vector<amrex::MultiFab> &input)
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


void
Operator::RegisterNewFab(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &input)
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

}

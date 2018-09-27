
#include "Operator.H"
#include <AMReX_MLLinOp_F.H>
#include <AMReX_MultiFabUtil.H>
#include "Util/Color.H"
#include "Set/Set.H"

#define TRACER	std::cout << Color::FG::Yellow << __FILE__ << ":" << __LINE__ << Color::FG::Default << " " << __func__ << std::endl;
#define PROBE	std::cout << Color::FG::Red << __FILE__ << ":" << __LINE__ << Color::FG::Default << " " << __func__ << std::endl;


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

bool Operator::VerificationCheck (int amrlev,
				  int mglev,
				  amrex::MultiFab& test) const
{
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
 }

void
Operator::reflux (int crse_amrlev,
		     MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
		     MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const
{
	Util::Abort("reflux not yet implemented");
}

void
Operator::compRHS (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel,
		   const Vector<const MultiFab*>& rhnd,
		   const Vector<MultiFab*>& a_rhcc)
{
	Util::Abort("compRHS not implemented");
}


void
Operator::averageDownCoeffs ()
{
	Util::Abort("averageDownCoeffs not implemented");
}

void
Operator::averageDownCoeffsToCoarseAmrLevel (int flev)
{
	Util::Abort("averageDownCoeffsToCoarseAmrLevel not implemented");
	const int mglev = 0;
	const int idim = 0;  // other dimensions are just aliases
	// amrex::average_down(*m_sigma[flev][mglev][idim], *m_sigma[flev-1][mglev][idim], 0, 1,
	// 		    m_amr_ref_ratio[flev-1]);
}

void
Operator::averageDownCoeffsSameAmrLevel (int amrlev)
{
	Util::Abort("averageDownCoeffsToSameAmrLevel not implemented");
}

void
Operator::FillBoundaryCoeff (MultiFab& sigma, const Geometry& geom)
{
	Util::Abort("FillBoundaryCoeff not implemented");
}

void
Operator::buildMasks ()
{
	if (m_masks_built) return;
	BL_PROFILE("Operator::buildMasks()");
	m_masks_built = true;
	m_is_bottom_singular = false;
	int amrlev = 0;
	int mglev = m_num_mg_levels[amrlev]-1;
	const iMultiFab& omask = *m_owner_mask[amrlev][mglev];
	m_bottom_dot_mask.define(omask.boxArray(), omask.DistributionMap(), 1, 0);
	m_bottom_dot_mask.setVal(1);
}


void
Operator::fixUpResidualMask (int amrlev, iMultiFab& resmsk)
{
	Util::Abort("fixUpResidualMask not implemented");
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

	//averageDownCoeffs();
}

void
Operator::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
	TRACER;
	if (AMREX_SPACEDIM != 3) Util::Abort("restriction implemented in 3D only!");

	if (fine.contains_nan() || fine.contains_inf()) Util::Abort("restriction (beginning) - nan or inf detected in fine");
	if (crse.contains_nan() || crse.contains_inf()) Util::Abort("restriction (beginning) - nan or inf detected in crse");

	BL_PROFILE("MLNodeLaplacian::restriction()");

	applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous, StateMode::Solution);

	const Box& nd_domain = amrex::surroundingNodes(m_geom[amrlev][cmglev].Domain());

	bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
	MultiFab cfine;
	if (need_parallel_copy) {
		const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
		cfine.define(ba, fine.DistributionMap(), fine.nComp(), 0);
	}

	MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;

	// C++ version of Fortran (amrex_mlndlap_restriction)
	Set::Scalar fac1 = 1.0/64.0;
	Set::Scalar fac2 = 1.0/32.0;
	Set::Scalar fac3 = 1.0/16.0; 
	Set::Scalar fac4 = 1.0/8.0;
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
			amrex::IntVect m_fine(2*m1, 2*m2, 2*m3);
			for (int i=0; i<crse.nComp(); i++)
			{
				crsefab(m_crse,i) =
					fac1*(finefab(m_fine-dx-dy-dz,i) +
					      finefab(m_fine-dx-dy+dz,i) +
					      finefab(m_fine-dx+dy-dz,i) +
					      finefab(m_fine-dx+dy+dz,i) +
					      finefab(m_fine+dx-dy-dz,i) +
					      finefab(m_fine+dx-dy+dz,i) +
					      finefab(m_fine+dx+dy-dz,i) +
					      finefab(m_fine+dx+dy+dz,i))
					+
					fac2*(finefab(m_fine-dy-dz,i) +
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
					      finefab(m_fine+dx+dy,i))
					+
					fac3*(finefab(m_fine-dx,i) +
					      finefab(m_fine-dy,i) +
					      finefab(m_fine-dz,i) +
					      finefab(m_fine+dx,i) +
					      finefab(m_fine+dy,i) +
					      finefab(m_fine+dz,i))
					+
					fac4*finefab(m_fine,i);
			}
		}
	}

	if (need_parallel_copy) {
		crse.ParallelCopy(cfine);
	}

	if (fine.contains_nan() || fine.contains_inf()) Util::Abort("restriction (end) - nan or inf detected in fine");
	if (crse.contains_nan() || crse.contains_inf()) Util::Abort("restriction (end) - nan or inf detected in crse");

}

void
Operator::interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const
{
	TRACER;
	BL_PROFILE("MLNodeLaplacian::interpolation()");

	if (fine.contains_nan() || fine.contains_inf()) Util::Abort("interpolation (beginning) - nan or inf detected in fine");
	if (crse.contains_nan() || crse.contains_inf()) Util::Abort("interpolation (beginning) - nan or inf detected in crse");

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
			
			AMREX_D_TERM(for (int m1 = fine_bx.loVect()[0]; m1<=fine_bx.hiVect()[0]; m1++),
				     for (int m2 = fine_bx.loVect()[1]; m2<=fine_bx.hiVect()[1]; m2++),
				     for (int m3 = fine_bx.loVect()[2]; m3<=fine_bx.hiVect()[2]; m3++))
			{
				amrex::IntVect m(m1, m2, m3);
				amrex::IntVect M(m1/2, m2/2, m3/2);

				for (int i=0; i<crse.nComp(); i++)
				{
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
					else // Centroid
					{
						tmpfab(m,i) = 0.125 * (crsefab(M,i) +
								       crsefab(M+dx,i) + crsefab(M+dy,i) + crsefab(M+dz,i) +
								       crsefab(M+dy+dz,i) + crsefab(M+dz+dx,i) + crsefab(M+dx+dy,i) +
								       crsefab(M+dx+dy+dz,i));
					}

					if (std::isinf(tmpfab(m,i)))
					{
						std::cout << m << M << std::endl;
						Util::Abort("Is Infinity");
					}
					if (std::isnan(tmpfab(m,i)))
					{
						std::cout << m << M << std::endl;
						Util::Abort("Is NaN");
					}
				}
			}

			fine[mfi].plus(tmpfab,fine_bx,fine_bx,0,0,fine.nComp());

			if (fine[mfi].contains_nan()) std::cout << __LINE__ << " fine[mfi] contains nan" << std::endl;
			if (fine[mfi].contains_inf()) std::cout << __LINE__ << " fine[mfi] contains nan" << std::endl;

		}
	}

	if (fine.contains_nan()) Util::Abort("interpolation (end) - nan detected in fine");
	if (fine.contains_inf()) Util::Abort("interpolation (end) - inf detected in fine");
	if (crse.contains_nan()) Util::Abort("interpolation (end) - nan detected in crse");
	if (crse.contains_inf()) Util::Abort("interpolation (end) - inf detected in crse");
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

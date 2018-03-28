
#include "MyTest.H"
#include "MyTest_F.H"

using namespace amrex;

void
MyTest::initProbPoisson ()
{
    for (int ilev = 0; ilev <= max_level; ++ilev)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(rhs[ilev], true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            actual_init_poisson(BL_TO_FORTRAN_BOX(bx),
                                BL_TO_FORTRAN_ANYD(rhs[ilev][mfi]),
                                BL_TO_FORTRAN_ANYD(exact_solution[ilev][mfi]),
                                geom[ilev].ProbLo(), geom[ilev].ProbHi(),
                                geom[ilev].CellSize());
        }

        solution[ilev].setVal(0.0);
    }
}

void
MyTest::initProbABecLaplacian ()
{
    for (int ilev = 0; ilev <= max_level; ++ilev)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(rhs[ilev], true); mfi.isValid(); ++mfi)
        {
	  const amrex::Box& box = mfi.tilebox();
	  amrex::BaseFab<amrex::Real> &sol_box = solution[ilev][mfi];
	  amrex::Box domain(geom[ilev].Domain());

	  // for (int i = box.loVect()[0]-bcoef[ilev].nGrow(); i<=box.hiVect()[0]+bcoef[ilev].nGrow(); i++)
	  //   for (int j = box.loVect()[1]-bcoef[ilev].nGrow(); j<=box.hiVect()[1]+bcoef[ilev].nGrow(); j++)
	  //     { 
	  // 	amrex::Real y = geom[ilev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[ilev].CellSize()[1];
	  // 	if (y > 0.5)
	  // 	  bcoef[ilev][mfi](amrex::IntVect(i,j)) = 2.0;
	  // 	else
	  // 	  bcoef[ilev][mfi](amrex::IntVect(i,j)) = 1.0;
	  //     }

	   acoef[ilev].setVal(0.0);
	   bcoef[ilev].setVal(1.0);
	   // bcoef[ilev].setVal(1.0);
        }
	ascalar = 1.0;
	bscalar = 1.0;
        solution[ilev].setVal(0.0);
        rhs[ilev].setVal(0.0);
        exact_solution[ilev].setVal(0.0);
    }
}

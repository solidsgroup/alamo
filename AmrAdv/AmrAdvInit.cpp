
#include <AmrAdv.H>
#include <AmrAdv_F.H>

using namespace amrex;

void
AmrAdv::InitData ()
{
    if (restart_chkfile.empty())
    {
	const Real time = 0.0;
	InitFromScratch(time);
	AverageDown();

	if (plot_int > 0) {
	    WritePlotFile();
	}
    }
    else
    {
	InitFromCheckpoint();
    }
}

void AmrAdv::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
				      const DistributionMapping& dm)
{
    const int ncomp = 1;
    const int nghost = 0;

    phi_new[lev].reset(new MultiFab(ba, dm, ncomp, nghost));
    phi_old[lev].reset(new MultiFab(ba, dm, ncomp, nghost));

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    Real cur_time = t_new[lev];

    MultiFab& state = *phi_new[lev];

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();

	amrex::BaseFab<Real> &phi_box = state[mfi];

	//double offset_x = 0.25, offset_y=0.25;
	double offset_x = 0., offset_y=0.;
	
	for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++)
	  for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++)
	    {
	      amrex::Real x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];
	      amrex::Real y = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];
	      phi_box(amrex::IntVect(i,j)) =  1. + exp(-( (x-offset_x)*(x-offset_x) + (y-offset_y)*(y-offset_y))/0.01);
	    }
	// initdata(lev, cur_time, ARLIM_3D(lo), ARLIM_3D(hi),
	// 	 BL_TO_FORTRAN_3D(state[mfi]), ZFILL(dx),
	// 	 ZFILL(prob_lo));
    }
}

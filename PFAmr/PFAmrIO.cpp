
#include <AMReX_PlotFileUtil.H>

#include <PFAmr.H>

using namespace amrex;

std::string
PFAmr::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

Array<const MultiFab*>
PFAmr::PlotFileMF () const
{
    Array<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
	r.push_back(phi_new[i].get());
    }
    return r;
}

Array<std::string>
PFAmr::PlotFileVarNames () const
{
    return {"phi"};
}

void
PFAmr::WritePlotFile () const
{
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();
    
    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
				    Geom(), t_new[0], istep, refRatio());
}

void
PFAmr::InitFromCheckpoint ()
{
    amrex::Abort("PFAmr::InitFromCheckpoint: todo");
}


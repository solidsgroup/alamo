#include <AMReX_Utility.H>
#include <AMReX_PlotFileUtil.H>

#include "PFAmr.H"

#include <fstream>


using namespace amrex;

void
PFAmr::CreateCleanDirectory () const
{
  amrex::UtilCreateCleanDirectory(plot_file, false);
}

std::string
PFAmr::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file+"/", lev, 5);
}

Array<const MultiFab*>
PFAmr::PlotFileMF (int fab) const
{
    Array<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
	r.push_back(phi_new[fab][i].get());
    }
    return r;
}

Array<std::string>
PFAmr::PlotFileVarNames () const
{
  Array<std::string> names;
  for (int n = 0; n < number_of_grains; n++)
    names.push_back(amrex::Concatenate("phi",n));
  names.push_back("phi_sum");
  names.push_back("boundaries");
  return names;
}

void
PFAmr::WritePlotFile () const
{
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF(0);
    const auto& varnames = PlotFileVarNames();
    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
				   Geom(), t_new[0], istep, refRatio());

    if (ParallelDescriptor::IOProcessor())
      {
	std::ofstream outfile;
	if (istep[0]==0) outfile.open(plot_file+"/output.visit",std::ios_base::out);
	else outfile.open(plot_file+".visit",std::ios_base::app);
	outfile << plotfilename + "/Header" << std::endl;
      }
}

void
PFAmr::InitFromCheckpoint ()
{
    amrex::Abort("PFAmr::InitFromCheckpoint: todo");
}


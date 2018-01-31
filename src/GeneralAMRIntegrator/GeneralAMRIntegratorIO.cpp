#include <AMReX_Utility.H>
#include <AMReX_PlotFileUtil.H>

#include "GeneralAMRIntegrator.H"

#include <fstream>


using namespace amrex;

void
GeneralAMRIntegrator::CreateCleanDirectory () const
{
  amrex::UtilCreateCleanDirectory(plot_file, false);
}

std::vector<std::string>
GeneralAMRIntegrator::PlotFileName (int lev) const
{
  std::vector<std::string> name;
  name.push_back(plot_file+"/");
  name.push_back(amrex::Concatenate("", lev, 5));
  return name;
}

Array<const MultiFab*>
GeneralAMRIntegrator::PlotFileMF (int fab) const
{
  Array<const MultiFab*> r;
  // TODO
  for (int i = 0; i <= finest_level; ++i) {
    //r.push_back(phi_new[fab][i].get());
  }
  return r;
}

Array<std::string>
GeneralAMRIntegrator::PlotFileVarNames () const
{
  //Array<std::string> names;
  // TODO
  // for (int n = 0; n < number_of_grains; n++)
  //   names.push_back(amrex::Concatenate("phi",n));
  // names.push_back("phi_sum");
  // names.push_back("boundaries");
  //return names;
}

void
GeneralAMRIntegrator::WritePlotFile () const
{
  // TODO
  // const std::vector<std::string>& plotfilename = PlotFileName(istep[0]);
  // const auto& mf = PlotFileMF(0);

 
    // /**
    // * \brief Copy from src to dst including nghost ghost cells.
    // * The two MultiFabs MUST have the same underlying BoxArray.
    // * The copy is local.  The parallel copy function is in FabArray.H
    // */




  const int nlevels = finest_level+1;

  int total_components = 0;
  for (int i = 0; i < number_of_fabs; i++)
    total_components += ncomp_array[i];

  Vector<MultiFab> plotmf(nlevels);
  
  for (int ilev = 0; ilev < nlevels; ++ilev)
    {
      plotmf[ilev].define(grids[ilev], dmap[ilev], total_components, 0);

      int n = 0;
      for (int i = 0; i < number_of_fabs; i++)
	{
	  // REFERENCE
	  // static void Copy (MultiFab&       dst,
	  //                   const MultiFab& src,
	  //                   int             srccomp,
	  //                   int             dstcomp,
	  //                   int             numcomp,
	  //                   int             nghost);

	  MultiFab::Copy(plotmf[ilev], *(*fab_array[i])[ilev], 0, n, ncomp_array[i], 0);
	  n += ncomp_array[i];
	}
    }

  const std::vector<std::string>& plotfilename = PlotFileName(istep[0]);
  
  WriteMultiLevelPlotfile(plotfilename[0]+plotfilename[1], nlevels, amrex::GetVecOfConstPtrs(plotmf), name_array,
			  Geom(), t_new[0],istep, refRatio());


  if (ParallelDescriptor::IOProcessor())
    {
      std::ofstream outfile;
      if (istep[0]==0) outfile.open(plot_file+"/output.visit",std::ios_base::out);
      else outfile.open(plot_file+"/output.visit",std::ios_base::app);
      outfile << plotfilename[1] + "/Header" << std::endl;
    }
}

void
GeneralAMRIntegrator::InitFromCheckpoint ()
{
  amrex::Abort("GeneralAMRIntegrator::InitFromCheckpoint: todo");
}


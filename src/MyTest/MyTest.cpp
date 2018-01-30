#include "MyTest.H"

#include "MLStiffnessMatrix.H"

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_IndexType.H>

using namespace amrex;

#define DEBUG std::cout << "\033[34m" << "DEBUG:     " << "\033[0m" << __FILE__ << ":" << __LINE__ << std::endl;


MyTest::MyTest ()
{
  readParameters();

  int nlevels = max_level+1;
  geom.resize(nlevels);
  grids.resize(nlevels);
  dmap.resize(nlevels);

  solution.resize(nlevels);
  bcdata.resize(nlevels);
  rhs.resize(nlevels);
  acoef.resize(nlevels);
  bcoef.resize(nlevels);

  // define simulation domain
  RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});

  // set periodicity
  std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
  Geometry::Setup(&rb, 0, is_periodic.data());

  Box domain0(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
  // Box domain0(IntVect{AMREX_D_DECL(0,0,0)},
  // 	      IntVect{AMREX_D_DECL(n_cell,n_cell,n_cell)},
  // 	      IntVect(AMREX_D_DECL(IndexType::NODE,IndexType::NODE,IndexType::NODE)));
  Box domain = domain0;

  for (int ilev = 0; ilev < nlevels; ++ilev)
    {
      geom[ilev].define(domain);
      domain.refine(ref_ratio);
    }

  domain = domain0;
  for (int ilev = 0; ilev < nlevels; ++ilev)
    {
      grids[ilev].define(domain);
      grids[ilev].maxSize(max_grid_size);
      domain.grow(-n_cell/4);   // fine level cover the middle of the coarse domain
      domain.refine(ref_ratio); 
    }

  for (int ilev = 0; ilev < nlevels; ++ilev)
    {
      dmap     [ilev].define(grids[ilev]);
      solution [ilev].define(grids[ilev], dmap[ilev], BL_SPACEDIM, 1); 
      bcdata   [ilev].define(grids[ilev], dmap[ilev], BL_SPACEDIM, 1);
      rhs      [ilev].define(grids[ilev], dmap[ilev], BL_SPACEDIM, 0);
      acoef    [ilev].define(grids[ilev], dmap[ilev], 1, 0);
      bcoef    [ilev].define(grids[ilev], dmap[ilev], 1, 1);
    }

  for (int ilev = 0; ilev < nlevels; ++ilev)
    {
      //acoef[ilev].setVal(0.0,domain,0); // problem here when using nodes
      //bcoef[ilev].setVal(1.0,domain,0);
      //rhs[ilev].setVal(0.0,domain,0);
      //solution[ilev].setVal(0.0,domain,0);
      acoef[ilev].setVal(0.0);
      bcoef[ilev].setVal(1.0);
      rhs[ilev].setVal(0.0);
      solution[ilev].setVal(0.0);
    }
}

void
MyTest::solve ()
{
  LPInfo info;
  info.setAgglomeration(agglomeration);
  info.setConsolidation(consolidation);

  const Real tol_rel = 1.e-10;
  const Real tol_abs = 0.0;

  const int nlevels = geom.size();

  MLStiffnessMatrix mlabec(geom, grids, dmap, info);

  mlabec.setMaxOrder(linop_maxorder);

  //
  // SET BOUNDARY CONDITIONS
  //

  mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,LinOpBCType::Dirichlet,LinOpBCType::Dirichlet)},
		     {AMREX_D_DECL(LinOpBCType::Dirichlet,LinOpBCType::Dirichlet,LinOpBCType::Dirichlet)});

  for (int ilev = 0; ilev < nlevels; ++ilev)
    {
      amrex::Box domain(geom[ilev].Domain());
      for (MFIter mfi(bcdata[ilev], true); mfi.isValid(); ++mfi)
	{
          const Box& box = mfi.tilebox();
          //const Box& box = mfi.nodaltilebox();
	  amrex::BaseFab<amrex::Real> &bcdata_box = bcdata[ilev][mfi];
	  for (int i = box.loVect()[0] - bcdata[ilev].nGrow(); i<=box.hiVect()[0] + bcdata[ilev].nGrow(); i++)
	    for (int j = box.loVect()[1] - bcdata[ilev].nGrow(); j<=box.hiVect()[1] + bcdata[ilev].nGrow(); j++)
	      {
		if (j > domain.hiVect()[1]) // Top boundary
		  {
		    bcdata_box(amrex::IntVect(i,j)) = 1.;
		  }
		else
		  bcdata_box(amrex::IntVect(i,j)) = 0.;
	      }
	}
      solution[ilev].setVal(0.0);
      mlabec.setLevelBC(ilev,&bcdata[ilev]);
    }

  //
  // SET COEFFICIENTS
  //

  mlabec.setScalars(ascalar, bscalar);
  for (int ilev = 0; ilev < nlevels; ++ilev)
    {
      mlabec.setACoeffs(ilev, acoef[ilev]);
            
      std::array<MultiFab,AMREX_SPACEDIM> face_bcoef;
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
	{
	  const BoxArray& ba = amrex::convert(bcoef[ilev].boxArray(),
					      IntVect::TheDimensionVector(idim)); 
	  face_bcoef[idim].define(ba, bcoef[ilev].DistributionMap(), 1, 0);
	}

      amrex::average_cellcenter_to_face({AMREX_D_DECL(&face_bcoef[0],&face_bcoef[1],&face_bcoef[2])},
					bcoef[ilev],
					geom[ilev]);

      mlabec.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(face_bcoef));
    }

  //
  // CONFIGURE SOLVER
  //

  MLMG mlmg(mlabec);
  mlmg.setMaxIter(max_iter);
  mlmg.setMaxFmgIter(max_fmg_iter);
  mlmg.setVerbose(verbose);
  mlmg.setCGVerbose(cg_verbose);
  mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);
}

void
MyTest::readParameters ()
{
  ParmParse pp;
  pp.query("max_level", max_level);
  pp.query("ref_ratio", ref_ratio);
  pp.query("n_cell", n_cell);
  pp.query("max_grid_size", max_grid_size);

  pp.query("composite_solve", composite_solve);

  pp.query("verbose", verbose);
  pp.query("cg_verbose", cg_verbose);
  pp.query("max_iter", max_iter);
  pp.query("max_fmg_iter", max_fmg_iter);
  pp.query("linop_maxorder", linop_maxorder);
  pp.query("agglomeration", agglomeration);
  pp.query("consolidation", consolidation);
}

void
MyTest::writePlotfile () const
{
    const int ncomp = (acoef.empty()) ? 4 : 6;
    Vector<std::string> varname = {"solution", "boundary", "rhs", "error"};
    if (!acoef.empty()) {
        varname.emplace_back("acoef");
        varname.emplace_back("bcoef");
    }

    const int nlevels = max_level+1;

    Vector<MultiFab> plotmf(nlevels);
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        plotmf[ilev].define(grids[ilev], dmap[ilev], ncomp, 0);
        MultiFab::Copy(plotmf[ilev], solution      [ilev], 0, 0, 1, 0);
        MultiFab::Copy(plotmf[ilev], bcdata        [ilev], 0, 1, 1, 0);
        MultiFab::Copy(plotmf[ilev], rhs           [ilev], 0, 2, 1, 0);
        MultiFab::Copy(plotmf[ilev], solution      [ilev], 0, 3, 1, 0);
        if (!acoef.empty()) {
            MultiFab::Copy(plotmf[ilev], acoef[ilev], 0, 4, 1, 0);
            MultiFab::Copy(plotmf[ilev], bcoef[ilev], 0, 5, 1, 0);
        }
    }

    WriteMultiLevelPlotfile("plot", nlevels, amrex::GetVecOfConstPtrs(plotmf),
                            varname, geom, 0.0, Vector<int>(nlevels, 0),
                            Vector<IntVect>(nlevels, IntVect{ref_ratio}));
}

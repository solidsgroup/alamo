
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

#include "Operator/Elastic/CubicRotation/CubicRotation.H"
//#include "MyTest/MyTest.H"
//#include "MyTest/MLStiffnessMatrix.H"
//#include "MyTest/MLHeatConduction.H"
//#include "MyTest/MyMLMG.H"
//#include "MyTest/MyMLCGSolver.H"

using namespace amrex;

int main (int argc, char* argv[])
{
  amrex::Initialize(argc, argv);


  int max_level = 1;//0;
  int ref_ratio = 2;
  int n_cell = 128;
  int max_grid_size = 64;
    
  bool composite_solve = true;

  int verbose = 2;
  int cg_verbose = 0;
  int max_iter = 100;//100;
  int max_fmg_iter = 0;
  int linop_maxorder = 2;
  bool agglomeration = true;
  bool consolidation = true;

  amrex::Vector<amrex::Geometry> geom;
  amrex::Vector<amrex::BoxArray> grids;
  amrex::Vector<amrex::DistributionMapping> dmap;
  amrex::Vector<amrex::MultiFab> solution;
  amrex::Vector<amrex::MultiFab> bcdata;
  amrex::Vector<amrex::MultiFab> rhs;
  amrex::Vector<amrex::MultiFab> acoef;
  amrex::Vector<amrex::MultiFab> bcoef;

  //
  // READ PARAMETERS
  //
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



  //
  // CONSTRUCTOR
  //
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
  std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,0,0)};
  Geometry::Setup(&rb, 0, is_periodic.data());
  Box domain0(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
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

  int number_of_components = AMREX_SPACEDIM;
  int number_of_ghost_cells = 2;
  for (int ilev = 0; ilev < nlevels; ++ilev)
    {
      dmap[ilev].define(grids[ilev]);
      solution      [ilev].define(grids[ilev], dmap[ilev], number_of_components, number_of_ghost_cells); 
      bcdata        [ilev].define(grids[ilev], dmap[ilev], number_of_components, number_of_ghost_cells);
      rhs           [ilev].define(grids[ilev], dmap[ilev], number_of_components, number_of_ghost_cells);
      acoef[ilev].define(grids[ilev], dmap[ilev], 1, 0);
      bcoef[ilev].define(grids[ilev], dmap[ilev], 1, 1);
    }

  for (int ilev = 0; ilev < nlevels; ++ilev)
    {
      acoef[ilev].setVal(1.0);
      bcoef[ilev].setVal(1.0);
      rhs[ilev].setVal(0.0,0,1);
      rhs[ilev].setVal(0.0,1,1);
      solution[ilev].setVal(0.0);
    }


  //
  // SOLVE
  //
  
  LPInfo info;
  info.setAgglomeration(agglomeration);
  info.setConsolidation(consolidation);
  //const Real tol_rel = 1.e-10;
  const Real tol_rel = 1e-6;
  const Real tol_abs = 0.0;
  nlevels = geom.size();
  //info.setMaxCoarseningLevel(0);
  Eigen::Matrix<amrex::Real, 3, 3> R;
  R << 0.8755949, -0.3817528,  0.2959702,
       0.4200312,  0.9043038, -0.0762129,
       -0.2385525,  0.1910484, 0.9521519;  //3d rotation matrix for 30 deg rotation angle about [1,2,3]
  //R << 0.866025404, -0.5, 0,
  //     0.5, 0.866025404, 0,
  //     0.0, 0.0, 1.0;  //3d rotation matrix for 30 deg rotation angle
  //R << 1.0, 0.0, 0.0,
  //     0.0, 1.0, 0.0,
  //     0.0, 0.0, 1.0;  //3d rotation matrix with 0 deg rotation angle
  amrex::Real E,nu,mu;
  E = 1.0; nu = 0.25; mu = 2.0;
  amrex::Real C11, C12, C44;
  C11 = E*(1-nu)/(1-nu-2.0*nu*nu);
  C12 = E*nu/(1-nu-2.0*nu*nu);
  C44 = mu;
  
  Operator::Elastic::CubicRotation mlabec(R, C11, C12, C44);
  mlabec.define(geom, grids, dmap, info);
  mlabec.setMaxOrder(linop_maxorder);
  
  // set boundary conditions

  mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Periodic,
				   LinOpBCType::Dirichlet,
				   LinOpBCType::Dirichlet)},
		     {AMREX_D_DECL(LinOpBCType::Periodic,
				   LinOpBCType::Dirichlet,
				   LinOpBCType::Dirichlet)});

  for (int ilev = 0; ilev < nlevels; ++ilev)
    {
      amrex::Box domain(geom[ilev].Domain());
      
      for (MFIter mfi(bcdata[ilev], true); mfi.isValid(); ++mfi)
	{
	  const Box& box = mfi.tilebox();

	  amrex::BaseFab<amrex::Real> &bcdata_box = bcdata[ilev][mfi];

	  for (int i = box.loVect()[0] - bcdata[ilev].nGrow(); i<=box.hiVect()[0] + bcdata[ilev].nGrow(); i++)
	    for (int j = box.loVect()[1] - bcdata[ilev].nGrow(); j<=box.hiVect()[1] + bcdata[ilev].nGrow(); j++)
	      { 
		if (j > domain.hiVect()[1]) // Top boundary
		  {
		    bcdata_box(amrex::IntVect(i,j),0) = 0.1;
		    bcdata_box(amrex::IntVect(i,j),1) = 0.0;
		  }
		else if (i > domain.hiVect()[0]) // Right boundary
		  {
		    bcdata_box(amrex::IntVect(i,j),0) = 0.0;
		    bcdata_box(amrex::IntVect(i,j),1) = 0.0;
		  }
		else 
		  {
		    bcdata_box(amrex::IntVect(i,j),0) = 0.0;
		    bcdata_box(amrex::IntVect(i,j),1) = 0.0;
		  }
	      }
	}
      solution[ilev].setVal(0.0);
      mlabec.setLevelBC(ilev,&bcdata[ilev]);
    }

  // set coefficients

  //mlabec.setScalars(ascalar, bscalar);
  for (int ilev = 0; ilev < nlevels; ++ilev)
    {
      //mlabec.setACoeffs(ilev, acoef[ilev]);
            
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

      //mlabec.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(face_bcoef));
    }

  // configure solver

  // MLCGSolver mlcg(mlabec);
  // mlcg.setVerbose(verbose);
  // mlcg.solve(solution[0],rhs[0],tol_rel,tol_abs);

  MLMG mlmg(mlabec);
  mlmg.setMaxIter(max_iter);
  mlmg.setMaxFmgIter(max_fmg_iter);
  mlmg.setVerbose(verbose);
  mlmg.setCGVerbose(cg_verbose);
  mlmg.setBottomSolver(MLMG::BottomSolver::bicgstab);

  mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);

  //
  // WRITE PLOT FILE
  //

   const int ncomp = (acoef.empty()) ? 4 : 6;
   Vector<std::string> varname = {"solution01", "solution02", "rhs", "rhs2"};
   if (!acoef.empty()) {
     varname.emplace_back("acoef");
     varname.emplace_back("bcoef");
   }

   nlevels = max_level+1;

   Vector<MultiFab> plotmf(nlevels);
   for (int ilev = 0; ilev < nlevels; ++ilev)
     {
       plotmf[ilev].define(grids[ilev], dmap[ilev], ncomp, 0);
       MultiFab::Copy(plotmf[ilev], solution      [ilev], 0, 0, 1, 0);
       MultiFab::Copy(plotmf[ilev], solution      [ilev], 1, 1, 1, 0);
       //MultiFab::Copy(plotmf[ilev], bcdata        [ilev], 0, 1, 1, 0);
       MultiFab::Copy(plotmf[ilev], rhs           [ilev], 0, 2, 1, 0);
       MultiFab::Copy(plotmf[ilev], rhs           [ilev], 1, 3, 1, 0);
       //       MultiFab::Copy(plotmf[ilev], solution      [ilev], 0, 3, 1, 0);
       if (!acoef.empty()) {
   	MultiFab::Copy(plotmf[ilev], acoef[ilev], 0, 4, 1, 0);
   	MultiFab::Copy(plotmf[ilev], bcoef[ilev], 0, 5, 1, 0);
       }
     }

   WriteMultiLevelPlotfile("output", nlevels, amrex::GetVecOfConstPtrs(plotmf),
   			  varname, geom, 0.0, Vector<int>(nlevels, 0),
   			  Vector<IntVect>(nlevels, IntVect{ref_ratio}));

  amrex::Finalize();
}

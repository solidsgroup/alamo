
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

//#include "MyTest/MLStiffnessMatrix.H"
#include "Operator/Elastic/Isotropic/Isotropic.H"
#include "Operator/FEM/FEM.H"
#include "Model/Solid/Elastic/Elastic.H"

using namespace amrex;

int main (int argc, char* argv[])
{
  amrex::Initialize(argc, argv);


  std::array<amrex::Real,2> body_force = {0.0, 0.0};

  std::array<amrex::Real,2> disp_bc_top    = {0.1, 0.0};
  std::array<amrex::Real,2> disp_bc_left   = {0.0, 0.0};
  std::array<amrex::Real,2> disp_bc_right  = {0.0, 0.0};
  std::array<amrex::Real,2> disp_bc_bottom = {0.0, 0.0};

  LinOpBCType bc_x = LinOpBCType::Periodic; //LinOpBCType::Periodic; LinOpBCType::Neumann;
  LinOpBCType bc_y = LinOpBCType::Dirichlet;

  bool use_fsmooth = false; 

  int max_level = 2;//0;
  int ref_ratio = 2;//2
  int n_cell = 32;//128;
  int max_grid_size = 64;//64;
    
  int verbose = 2;
  int cg_verbose = 4;
  int max_iter = 1000;//100;
  int max_fmg_iter = 0;
  int linop_maxorder = 2;
  bool agglomeration = true;
  bool consolidation = false;

  const Real tol_rel = 1.0e-5;
  const Real tol_abs = 1.0e-5;


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
  pp.query("verbose", verbose);
  pp.query("cg_verbose", cg_verbose);
  pp.query("max_iter", max_iter);
  pp.query("max_fmg_iter", max_fmg_iter);
  pp.query("linop_maxorder", linop_maxorder);
  pp.query("agglomeration", agglomeration);
  pp.query("consolidation", consolidation);

  // Operator::FEM::Element<Operator::FEM::Q4> element(1.0,1.0);
  // std::ofstream out("file.dat");
  // for (amrex::Real x = 0; x<=1.0; x+=0.01)
  //   for (amrex::Real y = 0; y<=1.0; y+=0.01)
  //     {
  // 	out << x << " ";
  // 	out << y << " ";
  // 	out << element.Phi<4>({x,y}) << " ";
  // 	out << element.DPhi<4>({x,y})[0] << " ";
  // 	out << element.DPhi<4>({x,y})[1] << " ";
  // 	// out << element.Phi<2>({x,y}) << " ";
  // 	// out << element.Phi<3>({x,y}) << " ";
  // 	// out << element.Phi<4>({x,y}) << " ";
  // 	out << std::endl;
  //     }
  // exit(0);

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
  std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL((bc_x == LinOpBCType::Periodic ? 1 : 0),
							  (bc_y == LinOpBCType::Periodic ? 1 : 0),
							  1)};
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
      const Real* dx = geom[ilev].CellSize();
      rhs[ilev].setVal(body_force[0]*dx[0]*dx[1],0,1);
      rhs[ilev].setVal(body_force[1]*dx[0]*dx[1],1,1);
      solution[ilev].setVal(0.0);
    }


  //
  // SOLVE
  //
  
  LPInfo info;
  info.setAgglomeration(agglomeration);
  info.setConsolidation(consolidation);
  //const Real tol_rel = 1.e-10;
  nlevels = geom.size();
  if (!use_fsmooth) info.setMaxCoarseningLevel(0); //  <<< put in to NOT require FSmooth
  Model::Solid::Elastic model;
  Operator::FEM::FEM mlabec(model);
  //Operator::Elastic::Isotropic mlabec;
  mlabec.define(geom, grids, dmap, info);
  mlabec.setMaxOrder(linop_maxorder);
  

  // set boundary conditions

  mlabec.setDomainBC({AMREX_D_DECL(bc_x,
				   bc_y,
				   LinOpBCType::Periodic)},
		     {AMREX_D_DECL(bc_x,
				   bc_y,
				   LinOpBCType::Periodic)});

  for (int ilev = 0; ilev < nlevels; ++ilev)
    {
      amrex::Box domain(geom[ilev].Domain());
      
      solution[ilev].setVal(0.0);

      for (MFIter mfi(bcdata[ilev], true); mfi.isValid(); ++mfi)
	{
	  const Box& box = mfi.tilebox();

	  amrex::BaseFab<amrex::Real> &bcdata_box = bcdata[ilev][mfi];

	  for (int i = box.loVect()[0] - bcdata[ilev].nGrow(); i<=box.hiVect()[0] + bcdata[ilev].nGrow(); i++)
	    for (int j = box.loVect()[1] - bcdata[ilev].nGrow(); j<=box.hiVect()[1] + bcdata[ilev].nGrow(); j++)
	      { 
		if (j > domain.hiVect()[1]) // Top boundary
		  {
		    bcdata_box(amrex::IntVect(i,j),0) = disp_bc_top[0];
		    bcdata_box(amrex::IntVect(i,j),1) = disp_bc_top[1];
		  }
		else if (j < domain.loVect()[1]) // Bottom boundary
		  {
		    bcdata_box(amrex::IntVect(i,j),0) = disp_bc_bottom[0];
		    bcdata_box(amrex::IntVect(i,j),1) = disp_bc_bottom[1];
		  }
		else if (i > domain.hiVect()[0]) // Right boundary
		  {
		    bcdata_box(amrex::IntVect(i,j),0) = disp_bc_right[0];
		    bcdata_box(amrex::IntVect(i,j),1) = disp_bc_right[1];
		  }
		else if (i < domain.loVect()[0]) // Left boundary 
		  {
		    bcdata_box(amrex::IntVect(i,j),0) = disp_bc_left[0];
		    bcdata_box(amrex::IntVect(i,j),1) = disp_bc_left[1];
		  }
	      }

	}
      //solution[ilev].setVal(0.0);
      mlabec.setLevelBC(ilev,&bcdata[ilev]);
    }

  for (int ilev = 0; ilev < nlevels; ++ilev)
    {
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
    }

  // configure solver
  MLMG mlmg(mlabec);
  mlmg.setMaxIter(max_iter);
  mlmg.setMaxFmgIter(max_fmg_iter);
  mlmg.setVerbose(verbose);
  mlmg.setCGVerbose(cg_verbose);
  //mlmg.setBottomSolver(MLMG::BottomSolver::bicgstab);
  mlmg.setBottomSolver(MLMG::BottomSolver::cg);
  if (!use_fsmooth) mlmg.setFinalSmooth(0); // <<< put in to NOT require FSmooth
  if (!use_fsmooth) mlmg.setBottomSmooth(0);  // <<< put in to NOT require FSmooth
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

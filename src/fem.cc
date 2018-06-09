 
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

#include "Operator/Elastic/Cubic/Cubic.H"
#include "Operator/Elastic/Isotropic/Isotropic.H"
#include "Model/Solid/Elastic/Elastic.H"
#include "Set/Set.H"

using namespace amrex;

int main (int argc, char* argv[])
{
	amrex::Initialize(argc, argv);

	Set::Vector body_force 		= {AMREX_D_DECL(0.0, 0.0, 0.0)}; 

	Set::Vector disp_bc_top    	= {AMREX_D_DECL(0.0, 0.0, 0.1)};
	Set::Vector disp_bc_left   	= {AMREX_D_DECL(0.0, 0.0, 0.0)};
	Set::Vector disp_bc_right  	= {AMREX_D_DECL(0.0, 0.0, 0.0)};
	Set::Vector disp_bc_bottom 	= {AMREX_D_DECL(0.0, 0.0, 0.0)};
	Set::Vector disp_bc_front  	= {AMREX_D_DECL(0.0, 0.0, 0.0)};
	Set::Vector disp_bc_back 	= {AMREX_D_DECL(0.0, 0.0, 0.0)};	// configure solver

	// MLCGSolver mlcg(mlabec);
	// mlcg.setVerbose(verbose);
	// mlcg.solve(solution[0],rhs[0],tol_rel,tol_abs);


	LinOpBCType bc_x = LinOpBCType::Dirichlet; //LinOpBCType::Periodic; LinOpBCType::Neumann;
	LinOpBCType bc_y = LinOpBCType::Dirichlet;
	LinOpBCType bc_z = LinOpBCType::Dirichlet;

	//bool use_fsmooth = true; 

	int max_level = 1;//0;
	int ref_ratio = 2;//2
	int n_cell = 32;//128;
	int max_grid_size = 64;//64;
    
	bool composite_solve = true;

	int verbose 		= 2;
	int cg_verbose 		= 0;
	int max_iter		= 1000;//100;
	int max_fmg_iter 	= 0;
	int linop_maxorder 	= 2;
	bool agglomeration 	= true;
	bool consolidation 	= false;

	const Real tol_rel 	= 1.0e-5;
	const Real tol_abs 	= 1.0e-5;


	amrex::Vector<amrex::Geometry> 			geom;
	amrex::Vector<amrex::BoxArray> 			grids;
	amrex::Vector<amrex::DistributionMapping> 	dmap;
	amrex::Vector<amrex::MultiFab> 			solution;
	amrex::Vector<amrex::MultiFab> 			bcdata;	
	amrex::Vector<amrex::MultiFab> 			rhs;
	amrex::Vector<amrex::MultiFab>			stress;
	amrex::Vector<amrex::MultiFab>			energy;

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
	stress.resize(nlevels);
	energy.resize(nlevels);

	// define simulation domain
	RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
	// set periodicity
	std::array<int,AMREX_SPACEDIM> is_periodic = {AMREX_D_DECL(bc_x == LinOpBCType::Periodic ? 1 : 0,
															   bc_y == LinOpBCType::Periodic ? 1 : 0,
															   bc_z == LinOpBCType::Periodic ? 1 : 0)};
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
	int number_of_stress_components = AMREX_SPACEDIM > 1 ? (AMREX_SPACEDIM > 2 ? 6 : 3 ): 1;
	int number_of_ghost_cells = 2;
	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		dmap[ilev].define(grids[ilev]);
		solution	[ilev].define(grids[ilev], dmap[ilev], number_of_components, number_of_ghost_cells); 
		bcdata		[ilev].define(grids[ilev], dmap[ilev], number_of_components, number_of_ghost_cells);
		rhs		    [ilev].define(grids[ilev], dmap[ilev], number_of_components, number_of_ghost_cells);
		stress		[ilev].define(grids[ilev], dmap[ilev], number_of_stress_components, number_of_ghost_cells);
		energy		[ilev].define(grids[ilev], dmap[ilev], 1, number_of_ghost_cells);
	}

	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		const Real* dx = geom[ilev].CellSize();
		Set::Scalar volume = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);

		rhs[ilev].setVal(body_force[0]*volume,0,1);
#if AMREX_SPACEDIM > 1
		rhs[ilev].setVal(body_force[1]*volume,1,1);
#if AMREX_SPACEDIM > 2
		rhs[ilev].setVal(body_force[2]*volume,2,1);
#endif
#endif
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
	//if (!use_fsmooth) info.setMaxCoarseningLevel(0); //  <<< put in to NOT require FSmooth
	//Model::Solid::Elastic model;
	//Operator::FEM::FEM mlabec(model);
	Operator::Elastic::Isotropic mlabec;
	mlabec.define(geom, grids, dmap, info);
	mlabec.setMaxOrder(linop_maxorder);
  

	// set boundary conditions

	mlabec.setDomainBC({AMREX_D_DECL(bc_x, bc_y, bc_z)},
					   {AMREX_D_DECL(bc_x, bc_y, bc_z)});

	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		amrex::Box domain(geom[ilev].Domain());
      
		for (MFIter mfi(bcdata[ilev], true); mfi.isValid(); ++mfi)
		{
			const Box& box = mfi.tilebox();

			amrex::BaseFab<amrex::Real> &bcdata_box = bcdata[ilev][mfi];

			for (int i = box.loVect()[0] - bcdata[ilev].nGrow(); i<=box.hiVect()[0] + bcdata[ilev].nGrow(); i++)
				for (int j = box.loVect()[1] - bcdata[ilev].nGrow(); j<=box.hiVect()[1] + bcdata[ilev].nGrow(); j++)
#if AMREX_SPACEDIM>2
					for (int k = box.loVect()[2] - bcdata[ilev].nGrow(); k<=box.hiVect()[2] + bcdata[ilev].nGrow(); k++)
#endif
					{ 
						amrex::IntVect x(AMREX_D_DECL(i,j,k));
						if (j > domain.hiVect()[1]) // Top boundary
						{
							AMREX_D_TERM(bcdata_box(x,0) = disp_bc_top[0];,
										 bcdata_box(x,1) = disp_bc_top[1];,
										 bcdata_box(x,2) = disp_bc_top[2];)
						}
						else if (j < domain.loVect()[1]) // Bottom boundary
						{
							AMREX_D_TERM(bcdata_box(x,0) = disp_bc_bottom[0];,
										 bcdata_box(x,1) = disp_bc_bottom[1];,
										 bcdata_box(x,2) = disp_bc_bottom[2];)
						}
						else if (i > domain.hiVect()[0]) // Right boundary
						{
							AMREX_D_TERM(bcdata_box(x,0) = disp_bc_right[0];,
										 bcdata_box(x,1) = disp_bc_right[1];,
										 bcdata_box(x,2) = disp_bc_right[2];)
						}
						else if (i < domain.loVect()[0]) // Left boundary 
						{
							AMREX_D_TERM(bcdata_box(x,0) = disp_bc_left[0];,
										 bcdata_box(x,1) = disp_bc_left[1];,
										 bcdata_box(x,2) = disp_bc_left[2];)
						}
#if AMREX_SPACEDIM>2
						else if (k > domain.hiVect()[2]) // Front boundary
						{
							AMREX_D_TERM(bcdata_box(x,0) = disp_bc_front[0];,
										 bcdata_box(x,1) = disp_bc_front[1];,
										 bcdata_box(x,2) = disp_bc_front[2];)
						}
						else if (k < domain.loVect()[2]) // Back boundary 
						{
							AMREX_D_TERM(bcdata_box(x,0) = disp_bc_back[0];,
										 bcdata_box(x,1) = disp_bc_back[1];,
										 bcdata_box(x,2) = disp_bc_back[2];)
						}
#endif
					}

		}
		mlabec.setLevelBC(ilev,&bcdata[ilev]);
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
	//mlmg.setBottomSolver(MLMG::BottomSolver::bicgstab);
	mlmg.setBottomSolver(MLMG::BottomSolver::cg);
	//mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
	// if (!use_fsmooth) mlmg.setFinalSmooth(0); // <<< put in to NOT require FSmooth
	// if (!use_fsmooth) mlmg.setBottomSmooth(0);  // <<< put in to NOT require FSmooth
	mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);


	// Computing stress and energy
	for (int lev = 0; lev < nlevels; lev++)
	{
		for ( amrex::MFIter mfi(solution[lev],true); mfi.isValid(); ++mfi )
		{
			FArrayBox &ufab  = (solution[lev])[mfi];
			FArrayBox &sigmafab  = (stress[lev])[mfi];
			FArrayBox &energyfab  = (energy[lev])[mfi];
			
			mlabec.Energy(energyfab,ufab,lev,mfi);
			mlabec.Stress(sigmafab,ufab,lev,mfi);
		}
	}
		


	//
	// WRITE PLOT FILE
	//

	const int ncomp = AMREX_SPACEDIM > 2 ? 13 : 8;
#if AMREX_SPACEDIM==2
	Vector<std::string> varname = {"solution01", "solution02", "rhs01", "rhs02", "stress11", "stress22", "stress12", "energy"};
#elif AMREX_SPACEDIM>2
	Vector<std::string> varname = {"solution01", "solution02", "solution03", "rhs01", "rhs02", "rhs03",
					"stress11", "stress22", "stress33", "stress23", "stress13", "stress12", "energy"};
#endif


	nlevels = max_level+1;

	Vector<MultiFab> plotmf(nlevels);
	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		plotmf[ilev].define(grids[ilev], dmap[ilev], ncomp, 0);
#if AMREX_SPACEDIM == 2
		MultiFab::Copy(plotmf[ilev], solution      [ilev], 0, 0, 1, 0);
		MultiFab::Copy(plotmf[ilev], solution      [ilev], 1, 1, 1, 0);
		MultiFab::Copy(plotmf[ilev], rhs           [ilev], 0, 2, 1, 0);
		MultiFab::Copy(plotmf[ilev], rhs           [ilev], 1, 3, 1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 0, 4, 1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 1, 5, 1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 2, 6, 1, 0);
		MultiFab::Copy(plotmf[ilev], energy	   [ilev], 0, 7, 1, 0);
#elif AMREX_SPACEDIM == 3
		MultiFab::Copy(plotmf[ilev], solution      [ilev], 0, 0, 1, 0);
		MultiFab::Copy(plotmf[ilev], solution      [ilev], 1, 1, 1, 0);
		MultiFab::Copy(plotmf[ilev], solution      [ilev], 2, 2, 1, 0);
		MultiFab::Copy(plotmf[ilev], rhs           [ilev], 0, 3, 1, 0);
		MultiFab::Copy(plotmf[ilev], rhs           [ilev], 1, 4, 1, 0);
		MultiFab::Copy(plotmf[ilev], rhs           [ilev], 2, 5, 1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 0, 6, 1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 1, 7, 1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 2, 8, 1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 3, 9, 1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 4, 10, 1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 5, 11, 1, 0);
		MultiFab::Copy(plotmf[ilev], energy	   [ilev], 0, 12, 1, 0);
#endif 
	}

	WriteMultiLevelPlotfile("output", nlevels, amrex::GetVecOfConstPtrs(plotmf),
		varname, geom, 0.0, Vector<int>(nlevels, 0),
		Vector<IntVect>(nlevels, IntVect{ref_ratio}));

	amrex::Finalize();
}

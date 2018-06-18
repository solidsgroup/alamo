 
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
#include "IO/WriteMetaData.H"
#include "IO/FileNameParse.H"

using namespace amrex;

int main (int argc, char* argv[])
{
	amrex::Initialize(argc, argv);

	amrex::Vector<amrex::Real> body_force = {AMREX_D_DECL(0.0, 0.0, 0.0)}; 
	amrex::Vector<amrex::Real> disp_bc_top = {AMREX_D_DECL(0.0, 0.0, 0.0)}; 
	amrex::Vector<amrex::Real> disp_bc_left = {AMREX_D_DECL(0.0, 0.0, 0.0)}; 
	amrex::Vector<amrex::Real> disp_bc_right = {AMREX_D_DECL(0.0, 0.0, 0.0)}; 
	amrex::Vector<amrex::Real> disp_bc_bottom = {AMREX_D_DECL(0.0, 0.0, 0.0)}; 
	LinOpBCType AMREX_D_DECL(bc_x_lo,bc_y_lo,bc_z_lo);
	LinOpBCType AMREX_D_DECL(bc_x_hi,bc_y_hi,bc_z_hi);


#if AMREX_SPACEDIM > 2
	amrex::Vector<amrex::Real> disp_bc_front;
	amrex::Vector<amrex::Real> disp_bc_back;
#endif
	std::string AMREX_D_DECL(bc_x_lo_str, bc_y_lo_str, bc_z_lo_str);
	std::string AMREX_D_DECL(bc_x_hi_str, bc_y_hi_str, bc_z_hi_str);

	{
		ParmParse pp("bc");
		pp.queryarr("body_force",     body_force     );
		pp.queryarr("disp_bc_top",    disp_bc_top     );
		pp.queryarr("disp_bc_left",   disp_bc_left    );
		pp.queryarr("disp_bc_right",  disp_bc_right   );
		pp.queryarr("disp_bc_bottom", disp_bc_bottom  );
		
		pp.query("bc_x_lo",bc_x_lo_str);
		pp.query("bc_y_lo",bc_y_lo_str);

		if (bc_x_lo_str == "EXT_DIR" || bc_x_lo_str == "DIRICHLET")    bc_x_lo = amrex::LinOpBCType::Dirichlet;
		else if (bc_x_lo_str == "NEUMANN") bc_x_lo = amrex::LinOpBCType::Neumann;
		else                               bc_x_lo = amrex::LinOpBCType::Periodic;

		if (bc_y_lo_str == "EXT_DIR" || bc_y_lo_str == "DIRICHLET")    bc_y_lo = amrex::LinOpBCType::Dirichlet;
		else if (bc_y_lo_str == "NEUMANN") bc_y_lo = amrex::LinOpBCType::Neumann;
		else                               bc_y_lo = amrex::LinOpBCType::Periodic;

		pp.query("bc_x_hi",bc_x_hi_str);
		pp.query("bc_y_hi",bc_y_hi_str);

		if (bc_x_hi_str == "EXT_DIR" || bc_x_hi_str == "DIRICHLET")    bc_x_hi = amrex::LinOpBCType::Dirichlet;
		else if (bc_x_hi_str == "NEUMANN") bc_x_hi = amrex::LinOpBCType::Neumann;
		else                        	   bc_x_hi = amrex::LinOpBCType::Periodic;

		if (bc_y_hi_str == "EXT_DIR" || bc_y_hi_str == "DIRICHLET")    bc_y_hi = amrex::LinOpBCType::Dirichlet;
		else if (bc_y_hi_str == "NEUMANN") bc_y_hi = amrex::LinOpBCType::Neumann;
		else                               bc_y_hi = amrex::LinOpBCType::Periodic;

		if(bc_x_lo == amrex::LinOpBCType::Periodic ||
			bc_x_hi == amrex::LinOpBCType::Periodic)
		{
			std::cout << "Warning: BC for x should be periodic for both ends. Resetting to periodic." << std::endl;
			bc_x_lo_str = "PERIODIC";
			bc_x_hi_str = "PERIODIC";
			bc_x_hi = amrex::LinOpBCType::Periodic;
			bc_x_lo = amrex::LinOpBCType::Periodic;
		}

		if(bc_y_lo == amrex::LinOpBCType::Periodic ||
			bc_y_hi == amrex::LinOpBCType::Periodic)
		{
			std::cout << "Warning: BC for y should be periodic for both ends. Resetting to periodic." << std::endl;
			bc_y_lo_str = "PERIODIC";
			bc_y_hi_str = "PERIODIC";
			bc_y_hi = amrex::LinOpBCType::Periodic;
			bc_y_lo = amrex::LinOpBCType::Periodic;
		}


#if AMREX_SPACEDIM > 2
		pp.queryarr("disp_bc_front",  disp_bc_front   );
		pp.queryarr("disp_bc_back",   disp_bc_back    );
		
		pp.query("bc_z_lo",bc_z_lo_str);
		if (bc_z_lo_str == "EXT_DIR" || bc_z_lo_str == "DIRICHLET")    bc_z_lo = amrex::LinOpBCType::Dirichlet;
		else if (bc_z_lo_str == "NEUMANN") bc_z_lo = amrex::LinOpBCType::Neumann;
		else				   bc_z_lo = amrex::LinOpBCType::Periodic;

		pp.query("bc_z_hi",bc_z_lo_str);
		if (bc_z_hi_str == "EXT_DIR" || bc_z_hi_str == "DIRICHLET")    bc_z_hi = amrex::LinOpBCType::Dirichlet;
		else if (bc_z_hi_str == "NEUMANN") bc_z_hi = amrex::LinOpBCType::Neumann;
		else				   bc_z_hi = amrex::LinOpBCType::Periodic;

		if(bc_z_lo == amrex::LinOpBCType::Periodic ||
			bc_z_hi == amrex::LinOpBCType::Periodic)
		{
			std::cout << "Warning: BC for z should be periodic for both ends. Resetting to periodic." << std::endl;
			bc_z_lo_str = "PERIODIC";
			bc_z_hi_str = "PERIODIC";
			bc_z_hi = amrex::LinOpBCType::Periodic;
			bc_z_lo = amrex::LinOpBCType::Periodic;
		}
#endif
	}
	
	int max_level = 1;//0;
	int ref_ratio = 2;//2
	int n_cell = 16;//128;
	int max_grid_size = 64;//64;
	bool composite_solve = true;
	int verbose 		= 2;
	int cg_verbose 		= 0;
	int max_iter		= 1000;//100;
	int max_fmg_iter 	= 0;
	int linop_maxorder 	= 2;
	bool agglomeration 	= true;
	bool consolidation 	= false;
	Real tol_rel 	= 1.0e-5;
	Real tol_abs 	= 1.0e-5;
	{
		ParmParse pp("solver");
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
		pp.query("tol_rel", tol_rel);
		pp.query("tol_abs", tol_abs);
	}

	std::string plot_file = "output";
	{
		amrex::ParmParse pp;
		pp.query("plot_file",plot_file);
	}


	amrex::Vector<amrex::Geometry> 			geom;
	amrex::Vector<amrex::BoxArray> 			grids;
	amrex::Vector<amrex::DistributionMapping> dmap;
	amrex::Vector<amrex::MultiFab> 			u;
	amrex::Vector<amrex::MultiFab> 			bcdata;	
	amrex::Vector<amrex::MultiFab> 			rhs;
	amrex::Vector<amrex::MultiFab>			stress;
	amrex::Vector<amrex::MultiFab>			energy;

	//
	// CONSTRUCTOR
	//
	int nlevels = max_level+1;
	geom.resize(nlevels);
	grids.resize(nlevels);
	dmap.resize(nlevels);

	u.resize(nlevels);
	bcdata.resize(nlevels);
	rhs.resize(nlevels);
	stress.resize(nlevels);
	energy.resize(nlevels);

	// define simulation domain
	RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
	// set periodicity
	std::array<int,AMREX_SPACEDIM> is_periodic = {AMREX_D_DECL(	bc_x_lo == LinOpBCType::Periodic ? 1 : 0,
									bc_y_lo == LinOpBCType::Periodic ? 1 : 0,
									bc_z_lo == LinOpBCType::Periodic ? 1 : 0)};
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
		u	[ilev].define(grids[ilev], dmap[ilev], number_of_components, number_of_ghost_cells); 
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
		u[ilev].setVal(0.0);
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

	mlabec.setDomainBC({AMREX_D_DECL(bc_x_lo, bc_y_lo, bc_z_lo)},
					   {AMREX_D_DECL(bc_x_hi, bc_y_hi, bc_z_hi)});

	BC::BC *mybc;
	mybc = new BC::BC(geom,
			{AMREX_D_DECL(bc_x_hi_str,bc_y_hi_str,bc_z_hi_str)},
			{AMREX_D_DECL(bc_x_lo_str,bc_y_lo_str,bc_z_lo_str)},
			disp_bc_left,disp_bc_right,
			disp_bc_bottom,disp_bc_top,
			disp_bc_back, disp_bc_front);

	mybc->FillBoundary(u,0,0,0.0);
				
//	for (int ilev = 0; ilev < nlevels; ++ilev)
//	{
//		amrex::Box domain(geom[ilev].Domain());
//      
//		for (MFIter mfi(bcdata[ilev], true); mfi.isValid(); ++mfi)
//		{
//			const Box& box = mfi.tilebox();
//
//			amrex::BaseFab<amrex::Real> &bcdata_box = bcdata[ilev][mfi];
//
//			for (int i = box.loVect()[0] - bcdata[ilev].nGrow(); i<=box.hiVect()[0] + bcdata[ilev].nGrow(); i++)
//				for (int j = box.loVect()[1] - bcdata[ilev].nGrow(); j<=box.hiVect()[1] + bcdata[ilev].nGrow(); j++)
//#if AMREX_SPACEDIM>2
//					for (int k = box.loVect()[2] - bcdata[ilev].nGrow(); k<=box.hiVect()[2] + bcdata[ilev].nGrow(); k++)
//#endif
//					{ 
//						amrex::IntVect x(AMREX_D_DECL(i,j,k));
//						if (j > domain.hiVect()[1]) // Top boundary
//						{
//							AMREX_D_TERM(bcdata_box(x,0) = disp_bc_top[0];,
//										 bcdata_box(x,1) = disp_bc_top[1];,
//										 bcdata_box(x,2) = disp_bc_top[2];)
//						}
//						else if (j < domain.loVect()[1]) // Bottom boundary
//						{
//							AMREX_D_TERM(bcdata_box(x,0) = disp_bc_bottom[0];,
//										 bcdata_box(x,1) = disp_bc_bottom[1];,
//										 bcdata_box(x,2) = disp_bc_bottom[2];)
//						}
//						else if (i > domain.hiVect()[0]) // Right boundary
//						{
//							AMREX_D_TERM(bcdata_box(x,0) = disp_bc_right[0];,
//										 bcdata_box(x,1) = disp_bc_right[1];,
//										 bcdata_box(x,2) = disp_bc_right[2];)
//						}
//						else if (i < domain.loVect()[0]) // Left boundary 
//						{
//							AMREX_D_TERM(bcdata_box(x,0) = disp_bc_left[0];,
//										 bcdata_box(x,1) = disp_bc_left[1];,
//										 bcdata_box(x,2) = disp_bc_left[2];)
//						}
//#if AMREX_SPACEDIM>2
//						else if (k > domain.hiVect()[2]) // Front boundary
//						{
//							AMREX_D_TERM(bcdata_box(x,0) = disp_bc_front[0];,
//										 bcdata_box(x,1) = disp_bc_front[1];,
//										 bcdata_box(x,2) = disp_bc_front[2];)
//						}
//						else if (k < domain.loVect()[2]) // Back boundary 
//						{
//							AMREX_D_TERM(bcdata_box(x,0) = disp_bc_back[0];,
//										 bcdata_box(x,1) = disp_bc_back[1];,
//										 bcdata_box(x,2) = disp_bc_back[2];)
//						}
//#endif
//					}
//
//		}
//		mlabec.setLevelBC(ilev,&bcdata[ilev]);
//	}
  

	// configure solver

	// MLCGSolver mlcg(mlabec);
	// mlcg.setVerbose(verbose);
	// mlcg.solve(u[0],rhs[0],tol_rel,tol_abs);

	MLMG mlmg(mlabec);
	mlmg.setMaxIter(max_iter);
	mlmg.setMaxFmgIter(max_fmg_iter);
	mlmg.setVerbose(verbose);
	mlmg.setCGVerbose(cg_verbose);
	mlmg.setBottomSolver(MLMG::BottomSolver::bicgstab);
	//mlmg.setBottomSolver(MLMG::BottomSolver::cg);
	//mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
	// if (!use_fsmooth) mlmg.setFinalSmooth(0); // <<< put in to NOT require FSmooth
	// if (!use_fsmooth) mlmg.setBottomSmooth(0);  // <<< put in to NOT require FSmooth
	mlmg.solve(GetVecOfPtrs(u), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);

	

	// Computing stress and energy
	for (int lev = 0; lev < nlevels; lev++)
	{
		//u[lev].FillBoundary(0,AMREX_SPACEDIM, geom[lev].periodicity(),0);

		for ( amrex::MFIter mfi(u[lev],true); mfi.isValid(); ++mfi )
		{
			FArrayBox &ufab  = (u[lev])[mfi];
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
	Vector<std::string> varname = {"u01", "u02", "rhs01", "rhs02", "stress11", "stress22", "stress12", "energy"};
#elif AMREX_SPACEDIM>2
	Vector<std::string> varname = {"u01", "u02", "u03", "rhs01", "rhs02", "rhs03",
					"stress11", "stress22", "stress33", "stress23", "stress13", "stress12", "energy"};
#endif


	nlevels = max_level+1;

	Vector<MultiFab> plotmf(nlevels);
	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		plotmf[ilev].define(grids[ilev], dmap[ilev], ncomp, 0);
#if AMREX_SPACEDIM == 2
		MultiFab::Copy(plotmf[ilev], u      [ilev], 0, 0, 1, 0);
		MultiFab::Copy(plotmf[ilev], u      [ilev], 1, 1, 1, 0);
		MultiFab::Copy(plotmf[ilev], rhs           [ilev], 0, 2, 1, 0);
		MultiFab::Copy(plotmf[ilev], rhs           [ilev], 1, 3, 1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 0, 4, 1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 1, 5, 1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 2, 6, 1, 0);
		MultiFab::Copy(plotmf[ilev], energy	   [ilev], 0, 7, 1, 0);
#elif AMREX_SPACEDIM == 3
		MultiFab::Copy(plotmf[ilev], u      [ilev], 0, 0,  1, 0);
		MultiFab::Copy(plotmf[ilev], u      [ilev], 1, 1,  1, 0);
		MultiFab::Copy(plotmf[ilev], u      [ilev], 2, 2,  1, 0);
		MultiFab::Copy(plotmf[ilev], rhs           [ilev], 0, 3,  1, 0);
		MultiFab::Copy(plotmf[ilev], rhs           [ilev], 1, 4,  1, 0);
		MultiFab::Copy(plotmf[ilev], rhs           [ilev], 2, 5,  1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 0, 6,  1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 1, 7,  1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 2, 8,  1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 3, 9,  1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 4, 10, 1, 0);
		MultiFab::Copy(plotmf[ilev], stress        [ilev], 5, 11, 1, 0);
		MultiFab::Copy(plotmf[ilev], energy	       [ilev], 0, 12, 1, 0);
#endif 
	}

	IO::FileNameParse(plot_file);

	WriteMultiLevelPlotfile(plot_file, nlevels, amrex::GetVecOfConstPtrs(plotmf),
		varname, geom, 0.0, Vector<int>(nlevels, 0),
		Vector<IntVect>(nlevels, IntVect{ref_ratio}));
	
	IO::WriteMetaData(plot_file);

	amrex::Finalize();
}

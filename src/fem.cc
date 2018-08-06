#include <streambuf>
 
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

#include "Util/Util.H"
#include "Operator/Elastic/Cubic/Cubic.H"
#include "Operator/Elastic/Isotropic/Isotropic.H"
#include "Model/Solid/Elastic/Elastic.H"
#include "Set/Set.H"
#include "IO/WriteMetaData.H"
#include "IO/FileNameParse.H"
#include "IC/Eigenstrain/Sphere.H"
#include "BC/Constant.H"

using namespace amrex;

int main (int argc, char* argv[])
{
	Util::Initialize(argc, argv);

	//
	//
	// READ INPUT FILE AND INSTANTIATE PARAMETERS
	//
	//


	// Input file
	amrex::ParmParse pp;
	std::string plot_file = "output"; pp.query("plot_file",plot_file);
	

	// Boundary Conditions
	ParmParse pp_bc("bc");
	amrex::Vector<amrex::Real> body_force;     pp_bc.queryarr("body_force",     body_force     );
	amrex::Vector<amrex::Real> disp_bc_top;    pp_bc.queryarr("disp_bc_top",    disp_bc_top     );
	amrex::Vector<amrex::Real> disp_bc_left;   pp_bc.queryarr("disp_bc_left",   disp_bc_left    );
	amrex::Vector<amrex::Real> disp_bc_right;  pp_bc.queryarr("disp_bc_right",  disp_bc_right   );
	amrex::Vector<amrex::Real> disp_bc_bottom; pp_bc.queryarr("disp_bc_bottom", disp_bc_bottom  );
#if AMREX_SPACEDIM > 2
	amrex::Vector<amrex::Real> disp_bc_front;  pp_bc.queryarr("disp_bc_front",  disp_bc_front   );
	amrex::Vector<amrex::Real> disp_bc_back;   pp_bc.queryarr("disp_bc_back",   disp_bc_back    );
#endif
	std::string bc_x_lo_str; pp_bc.query("bc_x_lo",bc_x_lo_str);
	std::string bc_x_hi_str; pp_bc.query("bc_x_hi",bc_x_hi_str);
	std::string bc_y_lo_str; pp_bc.query("bc_y_lo",bc_y_lo_str);
	std::string bc_y_hi_str; pp_bc.query("bc_y_hi",bc_y_hi_str);
#if AMREX_SPACEDIM > 2
	std::string bc_z_lo_str; pp_bc.query("bc_z_lo",bc_z_lo_str);
	std::string bc_z_hi_str; pp_bc.query("bc_z_hi",bc_z_hi_str);
#endif

	
	// Solver Parameters
	ParmParse pp_solver("solver");
	std::string bottom_solver = "cg";     pp_solver.query("bottom_solver",bottom_solver);      
	int max_level             = 1;		  pp_solver.query("max_level", max_level);             
	int ref_ratio             = 2;		  pp_solver.query("ref_ratio", ref_ratio);             
	int n_cell                = 16;		  pp_solver.query("n_cell", n_cell);                   
	int max_grid_size         = 64;		  pp_solver.query("max_grid_size", max_grid_size);     
	bool composite_solve      = true;	  pp_solver.query("composite_solve", composite_solve); 
	int verbose               = 2;		  pp_solver.query("verbose", verbose);                 
	int cg_verbose            = 0;		  pp_solver.query("cg_verbose", cg_verbose);           
	int max_iter              = 100;		  pp_solver.query("max_iter", max_iter);               
	int max_fmg_iter 	        = 0;		  pp_solver.query("max_fmg_iter", max_fmg_iter);       
	int linop_maxorder 	     = 2;		  pp_solver.query("linop_maxorder", linop_maxorder);   
	bool agglomeration 	     = true;	  pp_solver.query("agglomeration", agglomeration);     
	bool consolidation 	     = false;	  pp_solver.query("consolidation", consolidation);     
	Real tol_rel 	           = 1.0e-5;	  pp_solver.query("tol_rel", tol_rel);                 
	Real tol_abs 	           = 1.0e-5;	  pp_solver.query("tol_abs", tol_abs);                 
	bool use_fsmooth          = false;	  pp_solver.query("use_fsmooth", use_fsmooth);         



	amrex::Vector<amrex::Geometry> 			geom;
	amrex::Vector<amrex::BoxArray> 			grids;
	amrex::Vector<amrex::DistributionMapping> dmap;

	amrex::Vector<amrex::MultiFab>  u;
	amrex::Vector<amrex::MultiFab>  eps0;
	amrex::Vector<amrex::MultiFab>  bcdata;	
	amrex::Vector<amrex::MultiFab>  rhs;
	amrex::Vector<amrex::MultiFab>  stress;
	amrex::Vector<amrex::MultiFab>  energy;

	//
	// CONSTRUCTOR
	//
	int nlevels = max_level+1;
	geom.resize(nlevels);
	grids.resize(nlevels);
	dmap.resize(nlevels);

	u.resize(nlevels);
	eps0.resize(nlevels);
	bcdata.resize(nlevels);
	rhs.resize(nlevels);
	stress.resize(nlevels);
	energy.resize(nlevels);


	BC::Constant *mybc;
	mybc = new BC::Constant({AMREX_D_DECL(bc_x_hi_str,bc_y_hi_str,bc_z_hi_str)},
				{AMREX_D_DECL(bc_x_lo_str,bc_y_lo_str,bc_z_lo_str)}
				,disp_bc_left
				,disp_bc_right
				,disp_bc_bottom
				,disp_bc_top
#if AMREX_SPACEDIM>2
				,disp_bc_back
				,disp_bc_front
#endif
				);
	mybc->define(geom[0]); /// \todo get rid of this line, should work without it

	// define simulation domain
	RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
	std::array<int,AMREX_SPACEDIM> is_periodic = mybc->IsPeriodic();
	std::cout << "periodicity = " << is_periodic[0] << " " << is_periodic[1] << std::endl;
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
	int number_of_stress_components = AMREX_SPACEDIM*AMREX_SPACEDIM;
	int number_of_ghost_cells = 2;
	for (int ilev = 0; ilev < nlevels; ++ilev)
		{
			dmap   [ilev].define(grids[ilev]);
			u      [ilev].define(grids[ilev], dmap[ilev], number_of_components, number_of_ghost_cells); 
			eps0   [ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM*AMREX_SPACEDIM, number_of_ghost_cells); 
			bcdata [ilev].define(grids[ilev], dmap[ilev], number_of_components, number_of_ghost_cells);
			rhs    [ilev].define(grids[ilev], dmap[ilev], number_of_components, number_of_ghost_cells);
			stress [ilev].define(grids[ilev], dmap[ilev], number_of_stress_components, number_of_ghost_cells);
			energy [ilev].define(grids[ilev], dmap[ilev], 1, number_of_ghost_cells);
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
			
			// IC::Eigenstrain::Sphere eigenstrain_ic(geom);
			// eigenstrain_ic.Initialize(ilev,eps0);
		}

	//
	// Linear Operator
	//

	LPInfo info;

	info.setAgglomeration(agglomeration);
	info.setConsolidation(consolidation);
	nlevels = geom.size();

	Operator::Elastic::Isotropic mlabec;
	mlabec.define(geom, grids, dmap, *mybc, info);
	mlabec.setMaxOrder(linop_maxorder);

	//
	// THIS STUFF IS THE OLD WAY OF SETTING BOUNDARY CONDITIONS
	//
	// mlabec.SetEigenstrain(eps0,*mybc);
	// mlabec.AddEigenstrainToRHS(rhs);
	//{AMREX_D_DECL(amrex::LinOpBCType::Dirichlet,amrex::LinOpBCType::Dirichlet,amrex::LinOpBCType::Dirichlet)}
	// mlabec.setDomainBC(mybc->GetBCTypes<amrex::LinOpBCType>()[0],
	//   		   mybc->GetBCTypes<amrex::LinOpBCType>()[1]);


	// this must be replaced...
	mlabec.setDomainBC({AMREX_D_DECL(amrex::LinOpBCType::Dirichlet,amrex::LinOpBCType::Dirichlet,amrex::LinOpBCType::Dirichlet)},
	  		   {AMREX_D_DECL(amrex::LinOpBCType::Dirichlet,amrex::LinOpBCType::Dirichlet,amrex::LinOpBCType::Dirichlet)});
	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		mybc->define(geom[ilev]);
		mybc->FillBoundary(u[ilev],0,0,0.0);
		mlabec.setLevelBC(ilev,&u[ilev]);
	}



	//
	// Solver
	//

	MLMG mlmg(mlabec);
	mlmg.setMaxIter(max_iter);
	mlmg.setMaxFmgIter(max_fmg_iter);
	mlmg.setVerbose(verbose);
	mlmg.setCGVerbose(cg_verbose);
	mlmg.setFinalFillBC(true);	
	if (bottom_solver == "cg")
		mlmg.setBottomSolver(MLMG::BottomSolver::cg);
	else if (bottom_solver == "bicgstab")
		mlmg.setBottomSolver(MLMG::BottomSolver::bicgstab);
	if (!use_fsmooth)// <<< put in to NOT require FSmooth
		{
			mlmg.setFinalSmooth(0); 
			mlmg.setBottomSmooth(0); 
		}
	mlmg.solve(GetVecOfPtrs(u), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);

	//
	// Compute post-solve values
	//

	for (int lev = 0; lev < nlevels; lev++)
		{
			for ( amrex::MFIter mfi(u[lev],true); mfi.isValid(); ++mfi )
				{
					FArrayBox &ufab  = (u[lev])[mfi];
					FArrayBox &sigmafab  = (stress[lev])[mfi];
					FArrayBox &energyfab  = (energy[lev])[mfi];
			
					mlabec.Energy(energyfab,ufab,lev,mfi);
					mlabec.Stress(sigmafab,ufab,lev,mfi);
				}
		}
		
	// // RECOMPUTE RHS
	// for (int lev = 0; lev <= max_level; lev++)
	// 	mlabec.temp_Fapply(lev, 0, rhs[lev], u[lev]);


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
			MultiFab::Copy(plotmf[ilev], rhs    [ilev], 0, 2, 1, 0);
			MultiFab::Copy(plotmf[ilev], rhs    [ilev], 1, 3, 1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 0, 4, 1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 3, 5, 1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 1, 6, 1, 0);
			MultiFab::Copy(plotmf[ilev], energy [ilev], 0, 7, 1, 0);
#elif AMREX_SPACEDIM == 3
			MultiFab::Copy(plotmf[ilev], u      [ilev], 0, 0,  1, 0);
			MultiFab::Copy(plotmf[ilev], u      [ilev], 1, 1,  1, 0);
			MultiFab::Copy(plotmf[ilev], u      [ilev], 2, 2,  1, 0);
			MultiFab::Copy(plotmf[ilev], rhs    [ilev], 0, 3,  1, 0);
			MultiFab::Copy(plotmf[ilev], rhs    [ilev], 1, 4,  1, 0);
			MultiFab::Copy(plotmf[ilev], rhs    [ilev], 2, 5,  1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 0, 6,  1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 4, 7,  1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 8, 8,  1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 5, 9,  1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 2, 10, 1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 1, 11, 1, 0);
			MultiFab::Copy(plotmf[ilev], energy	[ilev], 0, 12, 1, 0);
#endif 
		}

	IO::FileNameParse(plot_file);

	WriteMultiLevelPlotfile(plot_file, nlevels, amrex::GetVecOfConstPtrs(plotmf),
									varname, geom, 0.0, Vector<int>(nlevels, 0),
									Vector<IntVect>(nlevels, IntVect{ref_ratio}));
	
	IO::WriteMetaData(plot_file);

	Util::Finalize();
}

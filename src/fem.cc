#include <streambuf>
#include <map>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MLNodeLaplacian.H>

#include "Util/Util.H"
#include "Operator/Diagonal.H"
#include "Operator/Elastic.H"
#include "Model/Solid/Elastic/Isotropic.H"
#include "Model/Solid/Elastic/Cubic.H"
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

	// Set::Matrix R;
	// R = Eigen::AngleAxisd(30, Set::Vector::UnitZ())*
	// 	Eigen::AngleAxisd(20, Set::Vector::UnitY())*
	// 	Eigen::AngleAxisd(10, Set::Vector::UnitZ());

	//using model_type = Model::Solid::Elastic::Cubic; model_type model(10.73, 6.09, 2.830); 
	using model_type = Model::Solid::Elastic::Isotropic; model_type model(2.6,6.0); 
	
	if (amrex::ParallelDescriptor::IOProcessor())
		Util::Message(INFO,model);

	//
	// READ INPUT FILE
	//

	// Input file
	amrex::ParmParse pp;
	std::string plot_file = "output"; pp.query("plot_file",plot_file);

	// Read in boundary conditions
	ParmParse pp_bc("bc");
	amrex::Vector<amrex::Real> body_force;     pp_bc.queryarr("body_force",     body_force      ); //  
	amrex::Vector<amrex::Real> disp_bc_top;    pp_bc.queryarr("disp_bc_top",    disp_bc_top     ); // Note: we are currently
	amrex::Vector<amrex::Real> disp_bc_left;   pp_bc.queryarr("disp_bc_left",   disp_bc_left    ); // using hard-coded values and
	amrex::Vector<amrex::Real> disp_bc_right;  pp_bc.queryarr("disp_bc_right",  disp_bc_right   ); // these parameters
	amrex::Vector<amrex::Real> disp_bc_bottom; pp_bc.queryarr("disp_bc_bottom", disp_bc_bottom  ); // are not currently
#if AMREX_SPACEDIM > 2										       // being used.
	amrex::Vector<amrex::Real> disp_bc_front;  pp_bc.queryarr("disp_bc_front",  disp_bc_front   ); // 
	amrex::Vector<amrex::Real> disp_bc_back;   pp_bc.queryarr("disp_bc_back",   disp_bc_back    ); // 
#endif
	amrex::Vector<std::string> AMREX_D_DECL(bc_x_lo_str,bc_y_lo_str,bc_z_lo_str);
	amrex::Vector<std::string> AMREX_D_DECL(bc_x_hi_str,bc_y_hi_str,bc_z_hi_str);
	AMREX_D_TERM(pp_bc.queryarr("bc_x_lo",bc_x_lo_str);,
		     pp_bc.queryarr("bc_y_lo",bc_y_lo_str);,
		     pp_bc.queryarr("bc_z_lo",bc_z_lo_str););
	AMREX_D_TERM(pp_bc.queryarr("bc_x_hi",bc_x_hi_str);,
		     pp_bc.queryarr("bc_y_hi",bc_y_hi_str);,
		     pp_bc.queryarr("bc_z_hi",bc_z_hi_str););
	std::map<std::string,Operator::Elastic<model_type>::BC> bc;
	bc["displacement"] = Operator::Elastic<model_type>::BC::Displacement;
	bc["disp"] = Operator::Elastic<model_type>::BC::Displacement;
	bc["traction"] = Operator::Elastic<model_type>::BC::Traction;
	bc["trac"] = Operator::Elastic<model_type>::BC::Traction;
	bc["periodic"] = Operator::Elastic<model_type>::BC::Periodic;
	std::array<Operator::Elastic<model_type>::BC,AMREX_SPACEDIM> bc_x_lo = {AMREX_D_DECL(bc[bc_x_lo_str[0]], bc[bc_x_lo_str[1]], bc[bc_x_lo_str[2]])};
	std::array<Operator::Elastic<model_type>::BC,AMREX_SPACEDIM> bc_x_hi = {AMREX_D_DECL(bc[bc_x_hi_str[0]], bc[bc_x_hi_str[1]], bc[bc_x_hi_str[2]])};
#if AMREX_SPACEDIM > 1
	std::array<Operator::Elastic<model_type>::BC,AMREX_SPACEDIM> bc_y_lo = {AMREX_D_DECL(bc[bc_y_lo_str[0]], bc[bc_y_lo_str[1]], bc[bc_y_lo_str[2]])};
	std::array<Operator::Elastic<model_type>::BC,AMREX_SPACEDIM> bc_y_hi = {AMREX_D_DECL(bc[bc_y_hi_str[0]], bc[bc_y_hi_str[1]], bc[bc_y_hi_str[2]])};
#endif
#if AMREX_SPACEDIM > 2
	std::array<Operator::Elastic<model_type>::BC,AMREX_SPACEDIM> bc_z_lo = {AMREX_D_DECL(bc[bc_z_lo_str[0]], bc[bc_z_lo_str[1]], bc[bc_z_lo_str[2]])};
	std::array<Operator::Elastic<model_type>::BC,AMREX_SPACEDIM> bc_z_hi = {AMREX_D_DECL(bc[bc_z_hi_str[0]], bc[bc_z_hi_str[1]], bc[bc_z_hi_str[2]])};
#endif

	// Read in solver parameters
	ParmParse pp_solver("solver");
	std::string bottom_solver = "cg";     pp_solver.query("bottom_solver",bottom_solver);      
	int max_level             = 1;		  pp_solver.query("max_level", max_level);             
	int ref_ratio             = 2;		  pp_solver.query("ref_ratio", ref_ratio);             
	int n_cell                = 16;		  pp_solver.query("n_cell", n_cell);                   
	int max_grid_size         = 64;		  pp_solver.query("max_grid_size", max_grid_size);     
	bool composite_solve      = true;	  pp_solver.query("composite_solve", composite_solve); 
	int verbose               = 2;		  pp_solver.query("verbose", verbose);                 
	int cg_verbose            = 0;		  pp_solver.query("cg_verbose", cg_verbose);           
	int max_iter              = 100;	  pp_solver.query("max_iter", max_iter);               
	int max_fmg_iter 	  = 0;		  pp_solver.query("max_fmg_iter", max_fmg_iter);       
	int max_mg_level          = 4;            pp_solver.query("max_mg_level", max_mg_level);
	int linop_maxorder 	  = 2;		  pp_solver.query("linop_maxorder", linop_maxorder);   
	bool agglomeration 	  = true;	  pp_solver.query("agglomeration", agglomeration);     
	bool consolidation 	  = false;	  pp_solver.query("consolidation", consolidation);     
	Real tol_rel 	          = 1.0e-5;	  pp_solver.query("tol_rel", tol_rel);                 
	Real tol_abs 	          = 1.0e-5;	  pp_solver.query("tol_abs", tol_abs);                 
	bool use_fsmooth          = false;	  pp_solver.query("use_fsmooth", use_fsmooth);         


	//
	// Define problem domain and field variables
	//
	// Note: both cell centered (cgrids, etc) and node centered (ngrids, etc)
	// are defined. Only node centered is used currently.
	//

	amrex::Vector<amrex::Geometry> 			geom;
	amrex::Vector<amrex::BoxArray> 			/*cgrids,*/ ngrids;
	amrex::Vector<amrex::DistributionMapping>       /*cdmap,*/ ndmap;

	amrex::Vector<amrex::MultiFab>  u;
	amrex::Vector<amrex::MultiFab>  res;
	amrex::Vector<amrex::MultiFab>  eps0;
	amrex::Vector<amrex::MultiFab>  bcdata;	
	amrex::Vector<amrex::MultiFab>  rhs;
	amrex::Vector<amrex::MultiFab>  stress;
	amrex::Vector<amrex::MultiFab>  energy;
	amrex::Vector<amrex::MultiFab>  verify;

	amrex::Vector<amrex::FabArray<amrex::BaseFab<model_type> > > modelfab;

	//
	// Define domain
	//
	int nlevels = max_level+1;
	geom.resize(nlevels);
	//cgrids.resize(nlevels);
	ngrids.resize(nlevels);
	//cdmap.resize(nlevels);
	ndmap.resize(nlevels);

	u.resize(nlevels);
	res.resize(nlevels);
	eps0.resize(nlevels);
	bcdata.resize(nlevels);
	rhs.resize(nlevels);
	stress.resize(nlevels);
	energy.resize(nlevels);
	verify.resize(nlevels);
	modelfab.resize(nlevels);

	RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
	Geometry::Setup(&rb, 0);

	Box NDomain(IntVect{AMREX_D_DECL(0,0,0)},
		    IntVect{AMREX_D_DECL(n_cell,n_cell,n_cell)},
	 	    IntVect::TheNodeVector());
	Box CDomain = amrex::convert(NDomain, IntVect::TheCellVector());

	Box domain = CDomain;
	for (int ilev = 0; ilev < nlevels; ++ilev)
		{
			geom[ilev].define(domain);
			domain.refine(ref_ratio);
		}
	Box cdomain = CDomain, ndomain = NDomain;
	for (int ilev = 0; ilev < nlevels; ++ilev)
		{
			// cgrids[ilev].define(cdomain);
			// cgrids[ilev].maxSize(max_grid_size);

			ngrids[ilev].define(ndomain);
			ngrids[ilev].maxSize(max_grid_size);

			//cdomain.grow(IntVect(AMREX_D_DECL(-n_cell/4,0,0))); 
			//cdomain.grow(IntVect(0,-n_cell/4)); 
			//cdomain.grow(-n_cell/4); 
			//cdomain.refine(ref_ratio); 

			ndomain.grow(IntVect(-n_cell/4,0));
			ndomain.refine(ref_ratio);
			//amrex::convert(cdomain,IntVect::TheNodeVector());
		}

	//
	// Initialize fields
	//

	int number_of_components = AMREX_SPACEDIM;
	int number_of_stress_components = AMREX_SPACEDIM*AMREX_SPACEDIM;
	int number_of_ghost_cells = 1; // \todo Reduce number of ghost cells to 0 |||| or not? 
	for (int ilev = 0; ilev < nlevels; ++ilev)
		{
			//cdmap   [ilev].define(cgrids[ilev]);
			ndmap   [ilev].define(ngrids[ilev]);
			u       [ilev].define(ngrids[ilev], ndmap[ilev], number_of_components, number_of_ghost_cells); 
			res     [ilev].define(ngrids[ilev], ndmap[ilev], number_of_components, number_of_ghost_cells); 
			eps0    [ilev].define(ngrids[ilev], ndmap[ilev], AMREX_SPACEDIM*AMREX_SPACEDIM, number_of_ghost_cells); 
			bcdata  [ilev].define(ngrids[ilev], ndmap[ilev], number_of_components, number_of_ghost_cells);
			rhs     [ilev].define(ngrids[ilev], ndmap[ilev], number_of_components, number_of_ghost_cells);
			stress  [ilev].define(ngrids[ilev], ndmap[ilev], number_of_stress_components, number_of_ghost_cells);
			energy  [ilev].define(ngrids[ilev], ndmap[ilev], 1, number_of_ghost_cells);
			verify	[ilev].define(ngrids[ilev], ndmap[ilev], number_of_components, number_of_ghost_cells);
			modelfab[ilev].define(ngrids[ilev], ndmap[ilev], 1, number_of_ghost_cells);
		}

	
	//res
	// Initialize fields, rhs, coefficients, etc.
	// Note: everything is currently hard-coded for testing purposes.
	//

	//Model::Solid::Elastic::Isotropic::Isotropic modeltype(2.6,6.0);
	//Model::Solid::Elastic::Cubic::Cubic modeltype(10.73, 6.09, 2.830);
	
	for (int ilev = 0; ilev < nlevels; ++ilev)
		{
			// const Real* dx = geom[ilev].CellSize();
			// Set::Scalar volume = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);

			
 			// AMREX_D_TERM(rhs[ilev].setVal(body_force[0]*volume,0,1);,
			// 	     rhs[ilev].setVal(body_force[1]*volume,1,1);,
			// 	     rhs[ilev].setVal(body_force[2]*volume,2,1);)

			rhs[ilev].setVal(0.0000);
			rhs[ilev].setVal(1.0,0);
			// rhs[ilev].setVal(100.0,1);
			u[ilev].setVal(0.0);
			stress[ilev].setVal(0.0);
			verify[ilev].setVal(0.0);
			modelfab[ilev].setVal(model);

			for (amrex::MFIter mfi(rhs[ilev],true); mfi.isValid(); ++mfi)
			{
			 	const amrex::Box& box = mfi.tilebox();

			 	amrex::BaseFab<amrex::Real> &rhsfab = (rhs[ilev])[mfi];

			 	AMREX_D_TERM(for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++),
			 		     for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++),
			 		     for (int k = box.loVect()[2]; k<=box.hiVect()[2]; k++))
			 	{
			 		if (AMREX_D_TERM(i == geom[ilev].Domain().loVect()[0],
							 || j == geom[ilev].Domain().loVect()[1],
							 || k == geom[ilev].Domain().loVect()[2] ) ||    
					    AMREX_D_TERM(i == geom[ilev].Domain().hiVect()[0]+1,
							 || j == geom[ilev].Domain().hiVect()[1]+1,  
							 || k == geom[ilev].Domain().hiVect()[2]+1))
			 		{
			 			AMREX_D_TERM(rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = 0.0;,
							     rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),1) = 0.0;,
							     rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),2) = 0.0;);
			 		}						
					if (i == geom[ilev].Domain().loVect()[0])
					    	rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = 0.1;
					// if (i == geom[ilev].Domain().hiVect()[0]+1)
					//    	rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = 0.1;
					// if (j == geom[ilev].Domain().hiVect()[0])
					//    	rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = 0.1;
			 	}
			}
			for (MFIter mfi(u[ilev], true); mfi.isValid(); ++mfi)
			{
				amrex::BaseFab<model_type> &C = modelfab[ilev][mfi];
				const amrex::FArrayBox &ufab    = u[ilev][mfi];
				Util::Message(INFO,"u lovect = ",amrex::IntVect(ufab.loVect()), ", C lovect = ",amrex::IntVect(C.loVect()));
				Util::Message(INFO,"u hivect = ",amrex::IntVect(ufab.hiVect()), ", C hivect = ",amrex::IntVect(C.hiVect()));

			}
		}

	//
	// Linear Operator
	//

	LPInfo info;
	info.setAgglomeration(agglomeration);
	info.setConsolidation(consolidation);
	//info.setMaxCoarseningLevel(0); // Multigrid does not work yet
	//info.setMaxCoarseningLevel(1); // Multigrid does not work yet
	info.setMaxCoarseningLevel(max_mg_level);
	nlevels = geom.size();

	// Operator::Diagonal mlabec;
	// mlabec.define(geom, cgrids, cdmap, info);
	// mlabec.setMaxOrder(linop_maxorder);

	Operator::Elastic<model_type> mlabec;
	mlabec.define(geom, ngrids, ndmap, info);
	//mlabec.setMaxOrder(linop_maxorder);
	for (int ilev = 0; ilev < nlevels; ++ilev) mlabec.SetModel(ilev,modelfab[ilev]);
	mlabec.SetBC({{AMREX_D_DECL(bc_x_lo,bc_y_lo,bc_z_lo)}},
		     {{AMREX_D_DECL(bc_x_hi,bc_y_hi,bc_z_hi)}});
	


	//mlabec.FsmoothTest(0, 0, u[0], rhs[0], res[0]);
	//mlabec.Diagonal(0, 0, res[0]);


	//
	// Solver
	//

	MLMG mlmg(mlabec);
	mlmg.setMaxIter(max_iter);
	mlmg.setMaxFmgIter(max_fmg_iter);
	mlmg.setVerbose(verbose);
	mlmg.setCGVerbose(cg_verbose);
	mlmg.setBottomMaxIter(200);
	mlmg.setFinalFillBC(true);	
	if (bottom_solver == "cg")
		mlmg.setBottomSolver(MLMG::BottomSolver::cg);
	else if (bottom_solver == "bicgstab")
		mlmg.setBottomSolver(MLMG::BottomSolver::bicgstab);
	else if (bottom_solver == "smoother")
		mlmg.setBottomSolver(MLMG::BottomSolver::smoother);
	//mlmg.setFinalSmooth(100); 
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
			mlabec.Stress(lev,stress[lev],u[lev]);
			mlabec.Energy(lev,energy[lev],u[lev]);
		}
		
	//
	// WRITE PLOT FILE
	//

#if AMREX_SPACEDIM==2
	const int ncomp = 10;
	Vector<std::string> varname = {"u01", "u02", "rhs01", "rhs02", "res01", "res02", "stress11", "stress22", "stress12", "energy"};
#elif AMREX_SPACEDIM==3
	const int ncomp = 16;
	Vector<std::string> varname = {"u01", "u02", "u03", "rhs01", "rhs02", "rhs03", "res01", "res02", "res03",
				       "stress11", "stress22", "stress33", "stress23", "stress13", "stress12", "energy"};
#endif

	nlevels = max_level+1;

	Vector<MultiFab> plotmf(nlevels);
	for (int ilev = 0; ilev < nlevels; ++ilev)
		{
			plotmf[ilev].define(ngrids[ilev], ndmap[ilev], ncomp, 0);
#if AMREX_SPACEDIM == 2
			MultiFab::Copy(plotmf[ilev], u      [ilev], 0, 0, 1, 0);
			MultiFab::Copy(plotmf[ilev], u      [ilev], 1, 1, 1, 0);
			MultiFab::Copy(plotmf[ilev], rhs    [ilev], 0, 2, 1, 0);
			MultiFab::Copy(plotmf[ilev], rhs    [ilev], 1, 3, 1, 0);
			MultiFab::Copy(plotmf[ilev], res    [ilev], 0, 4, 1, 0);
			MultiFab::Copy(plotmf[ilev], res    [ilev], 1, 5, 1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 0, 6, 1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 3, 7, 1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 1, 8, 1, 0);
			MultiFab::Copy(plotmf[ilev], energy [ilev], 0, 9, 1, 0);
#elif AMREX_SPACEDIM == 3
			MultiFab::Copy(plotmf[ilev], u      [ilev], 0, 0,  1, 0);
			MultiFab::Copy(plotmf[ilev], u      [ilev], 1, 1,  1, 0);
			MultiFab::Copy(plotmf[ilev], u      [ilev], 2, 2,  1, 0);
			MultiFab::Copy(plotmf[ilev], rhs    [ilev], 0, 3,  1, 0);
			MultiFab::Copy(plotmf[ilev], rhs    [ilev], 1, 4,  1, 0);
			MultiFab::Copy(plotmf[ilev], rhs    [ilev], 2, 5,  1, 0);
			MultiFab::Copy(plotmf[ilev], res    [ilev], 0, 6,  1, 0);
			MultiFab::Copy(plotmf[ilev], res    [ilev], 1, 7,  1, 0);
			MultiFab::Copy(plotmf[ilev], res    [ilev], 2, 8,  1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 0, 9,  1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 4, 10,  1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 8, 11,  1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 5, 12,  1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 2, 13, 1, 0);
			MultiFab::Copy(plotmf[ilev], stress [ilev], 1, 14, 1, 0);
			MultiFab::Copy(plotmf[ilev], energy	[ilev], 0, 15, 1, 0);
#endif 
		}

	IO::FileNameParse(plot_file);

	Util::Message(INFO,"varname size = ", varname.size());
	Util::Message(INFO,"mf->nComp() = ", plotmf[0].nComp());
	
	IO::WriteMetaData(plot_file);

	WriteMultiLevelPlotfile(plot_file, nlevels, amrex::GetVecOfConstPtrs(plotmf),
	 			varname, geom, 0.0, Vector<int>(nlevels, 0),
	 			Vector<IntVect>(nlevels, IntVect{ref_ratio}));
	

	Util::Finalize();
}

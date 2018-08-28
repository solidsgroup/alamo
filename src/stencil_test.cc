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
#include "BC/Elastic/Elastic.H"

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


	//BC::Constant *mybc;
	BC::Elastic *mybc;
	//mybc = new BC::Constant({AMREX_D_DECL(bc_x_hi_str,bc_y_hi_str,bc_z_hi_str)},
	//			{AMREX_D_DECL(bc_x_lo_str,bc_y_lo_str,bc_z_lo_str)},
	//			AMREX_D_DECL(disp_bc_left,disp_bc_bottom,disp_bc_back),
	//			AMREX_D_DECL(disp_bc_right,disp_bc_top,disp_bc_front));
	mybc = new BC::Elastic({AMREX_D_DECL(bc_x_hi_str,bc_y_hi_str,bc_z_hi_str)},
				{AMREX_D_DECL(bc_x_lo_str,bc_y_lo_str,bc_z_lo_str)},
				AMREX_D_DECL(disp_bc_left,disp_bc_bottom,disp_bc_back),
				AMREX_D_DECL(disp_bc_right,disp_bc_top,disp_bc_front));
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

	mybc->SetElasticOperator(&mlabec);

	amrex::Vector<Set::Vector> stencil;
	//amrex::Vector<Set::Vector> traction;
	///amrex::Vector<int> points;
	amrex::IntVect m(AMREX_D_DECL(0,6,5));

#if AMREX_SPACEDIM == 1
	static amrex::IntVect dx(1);
#elif AMREX_SPACEDIM == 2
	static amrex::IntVect dx(1,0), dy(0,1);
#elif AMREX_SPACEDIM == 3
	static amrex::IntVect dx(1,0,0), dy(0,1,0), dz(0,0,1);
#endif 

	/* Testing only one point unknown in the stencil - left boundary*/
	stencil.push_back(Set::Vector(AMREX_D_DECL(0.0008545907028,0.003169449969,-0.0001377147974)));
	AMREX_D_TERM(	stencil.push_back(Set::Vector(AMREX_D_DECL(0,0,0)));
			stencil.push_back(Set::Vector(AMREX_D_DECL(0.0005749858986,0.005498521306,-0.0001355297027)));
			,
			stencil.push_back(Set::Vector(AMREX_D_DECL(0.0006662620317,0.001609323344,-0.0002173323415)));
			stencil.push_back(Set::Vector(AMREX_D_DECL(0.001608359934,0.006584905797,3.210572128e-05)));
			,
			stencil.push_back(Set::Vector(AMREX_D_DECL(0.0008594510531,0.003207860629,-1.956894045e-05)));
			stencil.push_back(Set::Vector(AMREX_D_DECL(0.0008384021923,0.002920127078,-0.0003646250717)));
			);
	amrex::Vector<Set::Vector> trac;
	trac.push_back(Set::Vector(AMREX_D_DECL(disp_bc_left[0],disp_bc_left[1],disp_bc_left[2])));

	amrex::Vector<int> pts;
	pts.push_back(1);

	amrex::MFIter mfi (u[0], true);

	//std::cout << "calling StencilFill routine. at (" << i << "," << j << "," << k << "). point size = " << points.size() << std::endl;
	for (int p = 0; p < 7; p++)
		std::cout << "Stencil before at p = " << p+1 << ", " << stencil[p](0) << "," << stencil[p](1) << "," << stencil[p](2) << std::endl;
	std::cout << std::endl;
	
	mybc->StencilFill(stencil, trac, pts, m, 0, 0, mfi,true);

	for (int p = 0; p < 7; p++)
		std::cout << "Stencil after at p = " << p+1 << ", " << stencil[p](0) << "," << stencil[p](1) << "," << stencil[p](2) << std::endl;
	std::cout << std::endl;
}
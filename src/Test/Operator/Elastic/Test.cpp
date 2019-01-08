#include <streambuf>
#include <map>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MLNodeLaplacian.H>


#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLCGSolver.H>

#include "Test.H"
#include "Set/Set.H"
#include "IC/Trig.H"
#include "Operator/Elastic.H"
#include "Model/Solid/LinearElastic/Laplacian.H"

namespace Test
{
namespace Operator
{
namespace Elastic
{

void Test::Define(int _ncells,
		  int _nlevels)
{

	ncells = _ncells;
 	nlevels = _nlevels;
	int max_grid_size = 2000;//16;
	std::string orientation = "h";
 	geom.resize(nlevels);
 	cgrids.resize(nlevels);
 	ngrids.resize(nlevels);
 	dmap.resize(nlevels);

 	solution_exact.resize(nlevels);
 	solution_numeric.resize(nlevels);
 	solution_error.resize(nlevels);
 	rhs_prescribed.resize(nlevels);
 	rhs_numeric.resize(nlevels);
 	rhs_exact.resize(nlevels);
 	res_numeric.resize(nlevels);
 	res_exact.resize(nlevels);

	amrex::RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
	amrex::Geometry::Setup(&rb, 0);

	amrex::Box NDomain(amrex::IntVect{AMREX_D_DECL(0,0,0)},
			   amrex::IntVect{AMREX_D_DECL(ncells,ncells,ncells)},
			   amrex::IntVect::TheNodeVector());
	amrex::Box CDomain = amrex::convert(NDomain, amrex::IntVect::TheCellVector());

	amrex::Box domain = CDomain;
 	for (int ilev = 0; ilev < nlevels; ++ilev)
 		{
 			geom[ilev].define(domain);
 			domain.refine(ref_ratio);
 		}
	amrex::Box cdomain = CDomain;
 	for (int ilev = 0; ilev < nlevels; ++ilev)
 		{
 			cgrids[ilev].define(cdomain);
 			cgrids[ilev].maxSize(max_grid_size);
 			// if (orientation == "h")
			//cdomain.grow(amrex::IntVect(AMREX_D_DECL(0,-ncells/4,0))); 
 			// else if (orientation == "v")
 			// 	cdomain.grow(amrex::IntVect(AMREX_D_DECL(-ncells/4,0,0))); 
			cdomain.grow(amrex::IntVect(-ncells/4)); 
 			cdomain.refine(ref_ratio); 

 			ngrids[ilev] = cgrids[ilev];
 			ngrids[ilev].convert(amrex::IntVect::TheNodeVector());
 		}

 	int number_of_components = AMREX_SPACEDIM;
 	for (int ilev = 0; ilev < nlevels; ++ilev)
 		{
 			dmap   [ilev].define(cgrids[ilev]);
 			solution_numeric[ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 1); 
 			solution_exact  [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 1);
 			solution_error  [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 1);
 			rhs_prescribed  [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 1);
 			rhs_numeric     [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 1);
 			rhs_exact       [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 1);
 			res_numeric     [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 1); 
 			res_exact       [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 1); 
 		}

}

//
//
// Let Omega = [0, 1]^2
//
// Let u_i = u_i^{mn} sin(pi x m) sin(pi y n) and
//     b_i = b_i^{mn} sin(pi x m) sin(pi y n) and
// Then u_{i,jj} = - u_i^{mn} pi^2 (m^2 + n^2) sin(pi x m) sin(pi y n) 
//
// Governing equation
//   C_{ijkl} u_{k,jl} + b_i = 0
//
// Let C_{ijkl} = alpha delta_{ik} delta_{jl}
// then
//   C_{ijkl} u_{k,jl} = alpha delta_{ik} delta_{jl} u_{k,jl}
//                     = alpha u_{i,jj}
// so
//   - alpha u_i^{mn} pi^2 (m^2 + n^2) sin(pi x m) sin(pi y n) + b_i^{mn} sin(pi x m) sin(pi y n) = 0
// or
//     alpha u_i^{mn} pi^2 (m^2 + n^2) sin(pi x m) sin(pi y n) = b_i^{mn} sin(pi x m) sin(pi y n)
// using orthogonality
//     alpha u_i^{mn} pi^2 (m^2 + n^2) = b_i^{mn}
// or
//     u_i^{mn}  = b_i^{mn} / (alpha * pi^2 * (m^2 + n^2))
// 
//
int
Test::TrigTest(bool verbose, int component, int n, std::string plotfile)
{
	Set::Scalar tolerance = 0.001;

	int failed = 0;

	Set::Scalar alpha = 1.0;
	Model::Solid::LinearElastic::Laplacian model(alpha);
	amrex::Vector<amrex::FabArray<amrex::BaseFab<Model::Solid::LinearElastic::Laplacian> > >
		modelfab(nlevels); 

 	for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev].define(ngrids[ilev], dmap[ilev], 1, 1);
 	for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev].setVal(model);

	IC::Trig icrhs(geom,1.0,n,n);
	icrhs.SetComp(component);
	IC::Trig icexact(geom,-(1./2./Set::Constant::Pi/Set::Constant::Pi),n,n);
	icexact.SetComp(component);

	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		solution_exact  [ilev].setVal(0.0);
		solution_numeric[ilev].setVal(0.0);
		solution_error  [ilev].setVal(0.0);
		rhs_prescribed  [ilev].setVal(0.0);
		rhs_exact       [ilev].setVal(0.0);
		rhs_numeric     [ilev].setVal(0.0);
		res_exact       [ilev].setVal(0.0);
		res_numeric     [ilev].setVal(0.0);

		icrhs.Initialize(ilev,rhs_prescribed);
		icexact.Initialize(ilev,solution_exact);
		///icexact.Initialize(ilev,solution_numeric);
	}

	amrex::LPInfo info;
 	info.setAgglomeration(1);
 	info.setConsolidation(1);
 	info.setMaxCoarseningLevel(0);
 	nlevels = geom.size();

	::Operator::Elastic<Model::Solid::LinearElastic::Laplacian> elastic;

 	elastic.define(geom, cgrids, dmap, info);
 	//elastic.setMaxOrder(linop_maxorder);
	
 	using bctype = ::Operator::Elastic<Model::Solid::LinearElastic::Laplacian>::BC;

 	std::array<bctype,AMREX_SPACEDIM> bc_x_lo
 		= {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)};
 	std::array<bctype,AMREX_SPACEDIM> bc_x_hi 
 		= {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)};
#if AMREX_SPACEDIM > 1
 	std::array<bctype,AMREX_SPACEDIM> bc_y_lo
 		= {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)};
 	std::array<bctype,AMREX_SPACEDIM> bc_y_hi
 		= {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)};
#endif
#if AMREX_SPACEDIM > 2
 	std::array<bctype,AMREX_SPACEDIM> bc_z_lo
 		= {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)};
 	std::array<bctype,AMREX_SPACEDIM> bc_z_hi
 		= {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)};
#endif

 	for (int ilev = 0; ilev < nlevels; ++ilev) elastic.SetModel(ilev,modelfab[ilev]);
 	elastic.SetBC({{AMREX_D_DECL(bc_x_lo,bc_y_lo,bc_z_lo)}},
 		      {{AMREX_D_DECL(bc_x_hi,bc_y_hi,bc_z_hi)}});


	amrex::MLMG mlmg(elastic);
	mlmg.setMaxIter(20);
	mlmg.setMaxFmgIter(20);
 	if (verbose)
 	{
 		mlmg.setVerbose(2);
 		mlmg.setCGVerbose(2);
 	}
 	else
 	{
 		mlmg.setVerbose(0);
 		mlmg.setCGVerbose(0);
	}
 	mlmg.setBottomMaxIter(20);
 	mlmg.setFinalFillBC(false);	
 	mlmg.setBottomSolver(MLMG::BottomSolver::bicgstab);

	// Solution	
	Set::Scalar tol_rel = 1E-8;
	Set::Scalar tol_abs = 0.0;
 	mlmg.solve(GetVecOfPtrs(solution_numeric), GetVecOfConstPtrs(rhs_prescribed), tol_rel, tol_abs);

	// Compute solution error
	for (int i = 0; i < nlevels; i++)
	{
		amrex::MultiFab::Copy(solution_error[i],solution_numeric[i],component,component,1,0);
		amrex::MultiFab::Subtract(solution_error[i],solution_exact[i],component,component,1,0);
	}
	

	// Compute numerical right hand side
	for (int ilev = 0; ilev < nlevels; ilev++) elastic.FApply(ilev,0,rhs_numeric[ilev],solution_numeric[ilev]);

	// Compute exact right hand side
	for (int ilev = 0; ilev < nlevels; ilev++) elastic.FApply(ilev,0,rhs_exact[ilev]  ,solution_exact[ilev]);

	
	// Compute the numeric residual
	for (int i = 0; i < nlevels; i++)
	{
		amrex::MultiFab::Copy(res_numeric[i],rhs_numeric[i],component,component,1,1);
		amrex::MultiFab::Subtract(res_numeric[i],rhs_prescribed[i],component,component,1,1);
	}
	for (int ilev = nlevels-1; ilev > 0; ilev--)
	 	elastic.Reflux(0,
	 		       res_numeric[ilev-1], solution_numeric[ilev-1], rhs_prescribed[ilev-1],
	 		       res_numeric[ilev],   solution_numeric[ilev],   rhs_prescribed[ilev]);

	// Compute the exact residual
	for (int i = 0; i < nlevels; i++)
	{
		amrex::MultiFab::Copy(res_exact[i],rhs_exact[i],component,component,1,1);
		amrex::MultiFab::Subtract(res_exact[i],rhs_prescribed[i],component,component,1,1);
	}
	for (int ilev = nlevels-1; ilev > 0; ilev--)
	 	elastic.Reflux(0,
	 		       res_exact[ilev-1], solution_exact[ilev-1], rhs_prescribed[ilev-1],
	 		       res_exact[ilev],   solution_exact[ilev],   rhs_prescribed[ilev]);


	if (plotfile != "")
	{
		Util::Message(INFO,"Printing plot file to ",plotfile);
		const int output_comp = 16;
		Vector<std::string> varname = {"solution_exact1", "solution_exact2",
					       "solution_numeric1", "solution_numeric2",
					       "solution_error1", "solution_error2",
					       "rhs_prescribed1", "rhs_prescribed2",
					       "rhs_exact1","rhs_exact2",
					       "rhs_numeric1","rhs_numeric2",
					       "res_exact1","res_exact2",
					       "res_numeric1","res_numeric2"};
		Vector<MultiFab> plotmf(nlevels);
		for (int ilev = 0; ilev < nlevels; ++ilev)
		{
			plotmf			[ilev].define(ngrids[ilev], dmap[ilev], output_comp, 0);
			MultiFab::Copy(plotmf	[ilev], solution_exact [ilev], 0, 0,  2, 0); // 
			MultiFab::Copy(plotmf	[ilev], solution_numeric[ilev],0, 2,  2, 0); // 
			MultiFab::Copy(plotmf	[ilev], solution_error [ilev], 0, 4,  2, 0); // 
			MultiFab::Copy(plotmf	[ilev], rhs_prescribed [ilev], 0, 6,  2, 0); // 
			MultiFab::Copy(plotmf	[ilev], rhs_exact      [ilev], 0, 8,  2, 0); // 
			MultiFab::Copy(plotmf	[ilev], rhs_numeric    [ilev], 0, 10, 2, 0); // 
			MultiFab::Copy(plotmf	[ilev], res_exact      [ilev], 0, 12, 2, 0); // 
			MultiFab::Copy(plotmf	[ilev], res_numeric    [ilev], 0, 14, 2, 0); // 

		}


		amrex::WriteMultiLevelPlotfile(plotfile, nlevels, amrex::GetVecOfConstPtrs(plotmf),
					       varname, geom, 0.0, Vector<int>(nlevels, 0),
					       Vector<IntVect>(nlevels, IntVect{ref_ratio}));
	}

	// Find maximum solution error
	std::vector<Set::Scalar> error_norm(nlevels);
	for (int i = 0; i < nlevels; i++) error_norm[i] = solution_error[0].norm0(component,0,false) / solution_exact[0].norm0(component,0,false);
	Set::Scalar maxnorm = fabs(*std::max_element(error_norm.begin(),error_norm.end()));

	if (verbose) Util::Message(INFO,"relative error = ", 100*maxnorm, " %");

	if (maxnorm > tolerance) failed += 1;

	return failed;
}

}
}
}
	     

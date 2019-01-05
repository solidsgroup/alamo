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
	int max_grid_size = 100000;
	std::string orientation = "h";
 	geom.resize(nlevels);
 	cgrids.resize(nlevels);
 	ngrids.resize(nlevels);
 	dmap.resize(nlevels);

 	u.resize(nlevels);
 	res.resize(nlevels);
 	rhs.resize(nlevels);
 	exact.resize(nlevels);

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
 			u       [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 1); 
 			res     [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 1); 
 			rhs     [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 1);
 			exact   [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 1);
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
		u[ilev].setVal(0.0);
		res[ilev].setVal(0.0);

		rhs[ilev].setVal(0.0);
		icrhs.Initialize(ilev,rhs);
		exact[ilev].setVal(0.0);
		icexact.Initialize(ilev,exact);
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
 	mlmg.setBottomMaxIter(200);
 	mlmg.setFinalFillBC(false);	
 	mlmg.setBottomSolver(MLMG::BottomSolver::bicgstab);

	std::vector<Set::Scalar> initial_error_norm(nlevels);

	for (int i = 0; i < nlevels; i++)
	{
		amrex::MultiFab::Copy(res[i],exact[i],component,component,1,0);
		amrex::MultiFab::Subtract(res[i],u[i],component,component,1,0);
		initial_error_norm[i] = res[0].norm0(component,0,false);
	}
	

	Set::Scalar tol_rel = 1E-8;
	Set::Scalar tol_abs = 0.0;

 	mlmg.solve(GetVecOfPtrs(u), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);

	if (plotfile != "")
	{
		Util::Message(INFO,"Printing plot file to ",plotfile);
		const int output_comp = 8;
		Vector<std::string> varname = {"u1", "u2", "rhs1", "rhs2", "res1", "res2", "exact1", "exact2"};
		Vector<MultiFab> plotmf(nlevels);
		for (int ilev = 0; ilev < nlevels; ++ilev)
		{
			plotmf[ilev].define(ngrids[ilev], dmap[ilev], output_comp, 0);
			MultiFab::Copy(plotmf[ilev], u      [ilev], 0, 0, 1, 0);
			MultiFab::Copy(plotmf[ilev], u      [ilev], 1, 1, 1, 0);
			MultiFab::Copy(plotmf[ilev], rhs    [ilev], 0, 2, 1, 0);
			MultiFab::Copy(plotmf[ilev], rhs    [ilev], 1, 3, 1, 0);
			MultiFab::Copy(plotmf[ilev], res    [ilev], 0, 4, 1, 0);
			MultiFab::Copy(plotmf[ilev], res    [ilev], 1, 5, 1, 0);
			MultiFab::Copy(plotmf[ilev], exact  [ilev], 0, 6, 1, 0);
			MultiFab::Copy(plotmf[ilev], exact  [ilev], 1, 7, 1, 0);
		}


		amrex::WriteMultiLevelPlotfile(plotfile, nlevels, amrex::GetVecOfConstPtrs(plotmf),
					       varname, geom, 0.0, Vector<int>(nlevels, 0),
					       Vector<IntVect>(nlevels, IntVect{ref_ratio}));
	}

	std::vector<Set::Scalar> error_norm(nlevels);

	for (int i = 0; i < nlevels; i++)
	{
		amrex::MultiFab::Copy(res[i],exact[i],component,component,1,0);
		amrex::MultiFab::Subtract(res[i],u[i],component,component,1,0);
		error_norm[i] = res[0].norm0(component,0,false) / initial_error_norm[i];
	}

	Set::Scalar maxnorm = *std::max_element(error_norm.begin(),error_norm.end());

	if (maxnorm > tolerance)
	{
		if (verbose) Util::Message(INFO,"relative error = ", maxnorm);
		return 1;
	}
	else return 0;

}

}
}
}
	     

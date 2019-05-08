#include <AMReX_MLMG.H>

#include "Test/Operator/Elastic.H"
#include "Model/Solid/LinearElastic/Isotropic.H"
#include "IC/Affine.H"
#include "IC/Trig.H"
#include "Operator/Elastic.H"

namespace Test
{
namespace Operator
{
int Elastic::UniaxialTest(int verbose,
		 int component,
		 int n,
		 std::string plotfile)
{
        Set::Scalar tolerance = 0.01;

	int failed = 0;

	//using model_type = Model::Solid::LinearElastic::Laplacian;
	//Set::Scalar alpha = 1.0;
	//model_type model(alpha);

	using model_type = Model::Solid::LinearElastic::Isotropic;
	Set::Scalar lame=2.6, shear=6.0;
	model_type model(lame,shear);
	//model.Randomize();

	amrex::Vector<amrex::FabArray<amrex::BaseFab<model_type> > > modelfab(nlevels); 

 	for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev].define(ngrids[ilev], dmap[ilev], 1, 2);
 	for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev].setVal(model);

	Set::Vector myvec(AMREX_D_DECL(0.0, 0.0, 0.0));
	Set::Vector myvec1(AMREX_D_DECL(0.0, 0.1, 0.0));
	
	//rhs_prescribed - BODY FORCE
	std::complex<int> i(0,1);
	//IC::Trig icrhs(geom,1.0,AMREX_D_DECL(i,i,i));
	IC::Affine icrhs(geom, myvec, 1.0, Set::Vector::Zero() , false, 1.0);
	icrhs.SetComp(component);

	//Trig(Geometry , _alpha , AMREX_D_DECL(phi1,phi2,phi3) ) --- (Vector::geom, Scalar::Alpha , Scalar::AMREX)
	//Affine(IC(_geom), n(a_n), alpha(a_alpha), b(a_b), halfspace(a_halfspace), m(a_m)
	//Vector::geom , Vector::a_n, Scalar::a_alpha, Vector::a_b = Set::Vector::Zero(), bool::a_halfspace=false, Scalar::a_m = 1.0 

	//Set::Scalar dim = (Set::Scalar)(AMREX_SPACEDIM); //dim=2 -> 2D ; dim=3 ->3D
	//IC::Trig icexact(geom,-(1./dim/Set::Constant::Pi/Set::Constant::Pi),AMREX_D_DECL(n*i,n*i,n*i));
	IC::Affine icexact(geom, myvec1, 1.0, Set::Vector::Zero(), false, 1.0);
	icexact.SetComp(1); //exact solution
	
	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		icrhs.Initialize(ilev,rhs_prescribed);
		icexact.Initialize(ilev,solution_exact);
	}

	amrex::LPInfo info;
 	info.setAgglomeration(1);
 	info.setConsolidation(1);
 	//info.setMaxCoarseningLevel(2);
 	nlevels = geom.size();

	::Operator::Elastic<model_type> elastic;
	elastic.SetHomogeneous(false);
 	elastic.define(geom, cgrids, dmap, info);
	for (int ilev = 0; ilev < nlevels; ++ilev) elastic.SetModel(ilev,modelfab[ilev]);
	BC::Operator::Elastic<model_type> bc;

	bc.Set(bc.Face::XHI, bc.Direction::X, bc.Type::Displacement, 0.01, rhs_prescribed, geom);
	bc.Set(bc.Face::YLO, bc.Direction::X, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
	bc.Set(bc.Face::YLO, bc.Direction::Y, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
	bc.Set(bc.Face::YHI, bc.Direction::X, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
	bc.Set(bc.Face::YHI, bc.Direction::Y, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
	
	elastic.SetBC(&bc);



	amrex::MLMG mlmg(elastic);
	// mlmg.setMaxIter(100);
	// mlmg.setMaxFmgIter(20);
 	if (verbose)
 	{
 		mlmg.setVerbose(verbose);
		if (verbose > 4) mlmg.setCGVerbose(verbose);
 	}
 	else
 	{
 		mlmg.setVerbose(0);
 		mlmg.setCGVerbose(0);
	}
 	mlmg.setBottomMaxIter(50);
 	mlmg.setFinalFillBC(false);	
 	mlmg.setBottomSolver(MLMG::BottomSolver::bicgstab);

	// Solution
	
	Set::Scalar tol_rel = 1E-8;
	Set::Scalar tol_abs = 0;
 	mlmg.solve(GetVecOfPtrs(solution_numeric), GetVecOfConstPtrs(rhs_prescribed), tol_rel,tol_abs);

	// Compute solution error
	for (int i = 0; i < nlevels; i++)
	{
	        //amrex::MultiFab::Copy(solution_error[i],solution_numeric[i],component,component,1,1);
		amrex::MultiFab::Copy(solution_error[i],solution_numeric[i],0,0,AMREX_SPACEDIM,1);
		amrex::MultiFab::Subtract(solution_error[i],solution_exact[i],0,0,AMREX_SPACEDIM,1);
	}
	

	// Compute numerical right hand side
	mlmg.apply(GetVecOfPtrs(rhs_numeric),GetVecOfPtrs(solution_numeric));

	// Compute exact right hand side
	mlmg.apply(GetVecOfPtrs(rhs_exact),GetVecOfPtrs(solution_exact));

	
	// Compute the numeric residual
	for (int i = 0; i < nlevels; i++)
	{
		amrex::MultiFab::Copy(res_numeric[i],rhs_numeric[i],0,0,AMREX_SPACEDIM,0);
		amrex::MultiFab::Subtract(res_numeric[i],rhs_prescribed[i],0,0,AMREX_SPACEDIM,0);
	}
	for (int ilev = nlevels-1; ilev > 0; ilev--)
	 	elastic.Reflux(0,
	 		       res_numeric[ilev-1], solution_numeric[ilev-1], rhs_prescribed[ilev-1],
	 		       res_numeric[ilev],   solution_numeric[ilev],   rhs_prescribed[ilev]);

	// Compute the exact residual
	for (int i = 0; i < nlevels; i++)
	{
		amrex::MultiFab::Copy(res_exact[i],rhs_exact[i],0,0,AMREX_SPACEDIM,0);
		amrex::MultiFab::Subtract(res_exact[i],rhs_prescribed[i],0,0,AMREX_SPACEDIM,0);
	}
	for (int ilev = nlevels-1; ilev > 0; ilev--)
	 	elastic.Reflux(0,
	 		       res_exact[ilev-1], solution_exact[ilev-1], rhs_prescribed[ilev-1],
	 		       res_exact[ilev],   solution_exact[ilev],   rhs_prescribed[ilev]);

	// Compute the "ghost force" that introduces the error
	mlmg.apply(GetVecOfPtrs(ghost_force),GetVecOfPtrs(solution_error));

	if (plotfile != "")
	{
		Util::Message(INFO,"Printing plot file to ",plotfile);
		WritePlotFile(plotfile);

	}

	// Find maximum solution error
	Set::Scalar norm = 0.0, total = 0.0;
	for (int i = 0; i < nlevels; i++)
	  {
	    for (int j = 0; j < AMREX_SPACEDIM; j++)
	      {
		total += solution_exact[i].norm0(j);
		norm += solution_error[i].norm0(j) ;
	      }
	  }
	Set::Scalar maxnorm = fabs(norm/total);

	if (verbose) Util::Message(INFO,"relative error = ", 100*maxnorm, " %");
	          
	if (maxnorm > tolerance) failed += 1;

	return failed;
	
}
}
}

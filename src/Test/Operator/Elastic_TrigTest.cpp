#include "Elastic.H"
#include "Set/Set.H"
#include "IC/Trig.H"
#include "IC/Affine.H"
#include "IC/Random.H"
#include "Operator/Elastic.H"
#include "Model/Solid/LinearElastic/Laplacian.H"
#include "BC/Operator/Elastic.H"

namespace Test
{
namespace Operator
{
// This is a test that solves the equation
//     Lap(u) = 0
// with u = 0 on the boundary.
// Only one component is solved for, specified by the "component" argument.
int
Elastic::TrigTest(int verbose, int component, int n, std::string plotfile)
{
	Set::Scalar tolerance = 0.02;

	int failed = 0;

	// Define the "model" fab to be a Laplacian, so that this
	// elastic operator acts as a Laplacian on the "component-th" component of the fab.
	Set::Scalar alpha = 1.0;
	//using model_type = Model::Solid::LinearElastic::Isotropic; model_type model(2.6,6.0); 
	using model_type = Model::Solid::LinearElastic::Laplacian; model_type model(alpha);
	amrex::Vector<amrex::FabArray<amrex::BaseFab<model_type> > > modelfab(nlevels); 
 	for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev].define(ngrids[ilev], dmap[ilev], 1, 2);
 	for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev].setVal(model);

	// Initialize: set the rhs_prescribed to sin(n pi x1 / L) * sin(n pi x2 / L), so that
	// the exact solution is sin(n pi x1 / L) * sin(n pi x2 / L) / pi / 2.
	// Set everything else to zero.
	std::complex<int> i(0,1);
	IC::Trig icrhs(geom,1.0,AMREX_D_DECL(n*i,n*i,n*i),dim);
	icrhs.SetComp(component);
	//Set::Scalar dim = (Set::Scalar)(AMREX_SPACEDIM);
	IC::Trig icexact(geom,-(1./dim/Set::Constant::Pi/Set::Constant::Pi/n/n),AMREX_D_DECL(n*i,n*i,n*i),dim);
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
		ghost_force     [ilev].setVal(0.0);

		icrhs.Initialize(ilev,rhs_prescribed);
		icexact.Initialize(ilev,solution_exact);


		// Zero out the boundaries of the presecribed RHS
		for (amrex::MFIter mfi(rhs_prescribed[ilev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.tilebox();

			amrex::BaseFab<amrex::Real> &rhs = (rhs_prescribed[ilev])[mfi];

			AMREX_D_TERM(for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++),
				     for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++),
				     for (int k = box.loVect()[2]; k<=box.hiVect()[2]; k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				for (int p = 0; p<AMREX_SPACEDIM; p++)
				{
					AMREX_D_TERM( if (i == geom[ilev].Domain().loVect()[0]) rhs(m,p) = 0.0;,
						      if (j == geom[ilev].Domain().loVect()[1]) rhs(m,p) = 0.0;,
						      if (k == geom[ilev].Domain().loVect()[2]) rhs(m,p) = 0.0; );
					AMREX_D_TERM( if (i == geom[ilev].Domain().hiVect()[0]+1) rhs(m,p) = 0.0;,
						      if (j == geom[ilev].Domain().hiVect()[1]+1) rhs(m,p) = 0.0;,
						      if (k == geom[ilev].Domain().hiVect()[2]+1) rhs(m,p) = 0.0; );
				}
			}
		}

		// Uncomment this to initialize the numeric solution with the
		// exact solution

		// Util::Warning(INFO,"icexact.Initialize(ilev,solution_numeric);");
		// icexact.Initialize(ilev,solution_numeric); 
	}

	// Create and configure the Elastic operator
	// See:
	//   - src/Operator/Elastic.cpp   <<  the actual operator
	//   - src/Operator/Operator.cpp  <<  the abstract operator that Elastic inherits from. (where reflux is)
	amrex::LPInfo info;
 	info.setAgglomeration(1);
 	info.setConsolidation(1);
 	//info.setMaxCoarseningLevel(1);
 	nlevels = geom.size();
	::Operator::Elastic<model_type> elastic;
	elastic.SetHomogeneous(false);
 	elastic.define(geom, cgrids, dmap, info);
 	using bctype = ::Operator::Elastic<model_type>::BC;
 	for (int ilev = 0; ilev < nlevels; ++ilev) elastic.SetModel(ilev,modelfab[ilev]);

	Set::Vector disp;
	BC::Operator::Elastic<model_type> bc;
	bc.SetBC(bc.min,bc.min,bc.min,
		 BC::Operator::Elastic<model_type>::Displacement,
		 disp,
		 rhs_prescribed,
		 geom);
	//
	// Set up different versions of the problem based on the problem dimension (dim).
	// Note that this is not the actual dimension, just the dimension of the problem.
	//

	if (dim == 1)
	{
		elastic.SetBC(
			      {{AMREX_D_DECL({AMREX_D_DECL(bctype::Displacement,bctype::Traction,    bctype::Traction)},
					     {AMREX_D_DECL(bctype::Traction,    bctype::Displacement,bctype::Traction)},
					     {AMREX_D_DECL(bctype::Traction,    bctype::Traction,    bctype::Displacement)})}},
			      {{AMREX_D_DECL({AMREX_D_DECL(bctype::Displacement,bctype::Traction,    bctype::Traction)},
					     {AMREX_D_DECL(bctype::Traction,    bctype::Displacement,bctype::Traction)},
					     {AMREX_D_DECL(bctype::Traction,    bctype::Traction,    bctype::Displacement)})}});

	}
	if (dim == 2)
	{
		elastic.SetBC(
			      {{AMREX_D_DECL({AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Traction)},
					     {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Traction)},
					     {AMREX_D_DECL(bctype::Traction,    bctype::Traction,    bctype::Displacement)})}},
			      {{AMREX_D_DECL({AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Traction)},
					     {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Traction)},
					     {AMREX_D_DECL(bctype::Traction,    bctype::Traction,    bctype::Displacement)})}});
	}
	if (dim == 3)
	{
		elastic.SetBC(
			      {{AMREX_D_DECL({AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)},
					     {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)},
					     {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)})}},
			      {{AMREX_D_DECL({AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)},
					     {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)},
					     {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)})}});
	}

	// Create MLMG solver and solve
	amrex::MLMG mlmg(elastic);
	//mlmg.setFixedIter(1);
	//mlmg.setMaxIter(1000);
	//mlmg.setMaxFmgIter(50);
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
	// mlmg.setPreSmooth(4);
	// mlmg.setPostSmooth(4);

	Set::Scalar tol_rel = 1E-8;
	Set::Scalar tol_abs = 0;


	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		solution_error  [ilev].setVal(0.0);
		res_numeric     [ilev].setVal(0.0);
	}


 	mlmg.solve(GetVecOfPtrs(solution_numeric), GetVecOfConstPtrs(rhs_prescribed), tol_rel,tol_abs);

	// Call out specific value at C/F point
	// if (nlevels > 1)
	// 	Util::Message(INFO,
	// 		      "crse=",solution_numeric[0][amrex::MFIter(solution_numeric[0])](amrex::IntVect(AMREX_D_DECL(4,8,8))), " ",
	// 		      "fine=",solution_numeric[1][amrex::MFIter(solution_numeric[1])](amrex::IntVect(AMREX_D_DECL(8,16,16))));

	// Compute solution error
	for (int i = 0; i < nlevels; i++)
	{
	  	amrex::MultiFab::Copy(solution_error[i],solution_numeric[i],component,component,1,2);
	  	amrex::MultiFab::Subtract(solution_error[i],solution_exact[i],component,component,1,2);
	}

	//Compute numerical right hand side
	mlmg.apply(GetVecOfPtrs(rhs_numeric),GetVecOfPtrs(solution_numeric));

	// Check to make sure that point didn't change
	// if (nlevels > 1)
	// 	Util::Message(INFO,
	// 		      "crse=",solution_numeric[0][amrex::MFIter(solution_numeric[0])](amrex::IntVect(AMREX_D_DECL(4,8,8))), " ",
	// 		      "fine=",solution_numeric[1][amrex::MFIter(solution_numeric[1])](amrex::IntVect(AMREX_D_DECL(8,16,16))));

	// Compute exact right hand side
	mlmg.apply(GetVecOfPtrs(rhs_exact),GetVecOfPtrs(solution_exact));
	// elastic.Fapply(0,0,rhs_exact[0],solution_exact[0]);
	// elastic.Fapply(1,0,rhs_exact[1],solution_exact[1]);
	// elastic.reflux(0, rhs_exact[0], rhs_exact[0], rhs_exact[0], rhs_exact[1], rhs_exact[1], rhs_exact[1]);

	// amrex::MLCGSolver mlcg(&mlmg,elastic);
	// elastic.prepareForSolve();
	// mlcg.solve(solution_numeric[0],rhs_prescribed[0],tol_rel,tol_abs);


	// Compute numerical residual
	mlmg.compResidual(GetVecOfPtrs(res_numeric),GetVecOfPtrs(solution_numeric),GetVecOfConstPtrs(rhs_prescribed));

	// Compute exact residual
	mlmg.compResidual(GetVecOfPtrs(res_exact),GetVecOfPtrs(solution_exact),GetVecOfConstPtrs(rhs_prescribed));

	// Compute the "ghost force" that introduces the error
	// mlmg.apply(GetVecOfPtrs(ghost_force),GetVecOfPtrs(solution_error));
	
	// for (int i = 0; i < nlevels; i++)
	// {
	//  	amrex::MultiFab::Copy(ghost_force[i],rhs_prescribed[i],component,component,1,2);
	//  	amrex::MultiFab::Subtract(ghost_force[i],rhs_numeric[i],component,component,1,2);
	// }

	// // Normalize so that they are all per unit area
	// for (int i = 0; i < nlevels; i++)
	// {
	// 	const Real* DX = geom[i].CellSize();
	// 	ghost_force[i].mult(1.0/DX[0]/DX[1]);
	// }	


	


	// for (int alev=0; alev < nlevels; alev++)
	// {
	//  	varname[14] = "res";    MultiFab::Copy(res_numeric[alev], mlmg.res[alev][0], 0, 0, AMREX_SPACEDIM, 2);
	//  	varname[4] = "cor";     MultiFab::Copy(solution_error[alev], *mlmg.cor[alev][0], 0, 0, AMREX_SPACEDIM, 2);    
	//  	varname[6] = "rhs";     MultiFab::Copy(rhs_prescribed[alev], mlmg.rhs[alev], 0, 0, AMREX_SPACEDIM, 2);
	//  	varname[16] = "rescor"; MultiFab::Copy(ghost_force[alev], mlmg.rescor[alev][0], 0, 0, AMREX_SPACEDIM, 2);
	// }
	


	// Output plot file
	if (plotfile != "")
	{
		Util::Message(INFO,"Printing plot file to ",plotfile);
		WritePlotFile(plotfile);
		// std::vector<int> nghosts(nlevels); 
		// for (int i = 0; i < nlevels; i++) nghosts[i] = 2; nghosts[0] = 0;
		// WritePlotFile(plotfile,{0,2,2});
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

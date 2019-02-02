#include <streambuf>
#include <map>
#include <complex>

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

#include "Elastic.H"
#include "Set/Set.H"
#include "IC/Trig.H"
#include "IC/Affine.H"
#include "IC/Random.H"
#include "Operator/Elastic.H"
#include "Model/Solid/LinearElastic/Laplacian.H"

namespace Test
{
namespace Operator
{
void Elastic::Define(int _ncells,
		     int _nlevels,
		     int _dim,
		     Grid _config)
{
	dim = _dim;
	ncells = _ncells;
 	nlevels = _nlevels;
	int max_grid_size = 100000;//ncells/2;
	//std::string orientation = "h";
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
	ghost_force.resize(nlevels);

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
			if (_config == Grid::XYZ)
				cdomain.grow(amrex::IntVect(-ncells/4)); 
			else if (_config == Grid::X)
				cdomain.grow(amrex::IntVect(AMREX_D_DECL(-ncells/4,0,0)));
			else if (_config == Grid::Y)
				cdomain.grow(amrex::IntVect(AMREX_D_DECL(0,-ncells/4,0)));
			else if (_config == Grid::Z)
				cdomain.grow(amrex::IntVect(AMREX_D_DECL(0,0,-ncells/4)));
			else if (_config == Grid::YZ)
				cdomain.grow(amrex::IntVect(AMREX_D_DECL(0,-ncells/4,-ncells/4)));
			else if (_config == Grid::ZX)
				cdomain.grow(amrex::IntVect(AMREX_D_DECL(-ncells/4,0,-ncells/4)));
			else if (_config == Grid::XY)
				cdomain.grow(amrex::IntVect(AMREX_D_DECL(-ncells/4,-ncells/4,0)));
			
			//cdomain.growHi(1,-ncells/2); 

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
 			ghost_force     [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 1); 
 		}

}

int
Elastic::SPD(bool verbose,std::string plotfile)
{
	int failed = 0;

	using model_type = Model::Solid::LinearElastic::Isotropic; model_type model(2.6,6.0); 
	//using model_type = Model::Solid::LinearElastic::Laplacian; model_type model(1.0); 

	amrex::Vector<amrex::FabArray<amrex::BaseFab<model_type> > > modelfab(nlevels); 
 	for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev].define(ngrids[ilev], dmap[ilev], 1, 1);
 	for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev].setVal(model);

	// Clear all multifabs
	for (int ilev = 0; ilev < nlevels; ilev++) solution_exact  [ilev].setVal(0.0);
	for (int ilev = 0; ilev < nlevels; ilev++) solution_numeric[ilev].setVal(0.0);
	for (int ilev = 0; ilev < nlevels; ilev++) solution_error  [ilev].setVal(0.0);
	for (int ilev = 0; ilev < nlevels; ilev++) rhs_prescribed  [ilev].setVal(0.0);
	for (int ilev = 0; ilev < nlevels; ilev++) rhs_exact       [ilev].setVal(0.0);
	for (int ilev = 0; ilev < nlevels; ilev++) rhs_numeric     [ilev].setVal(0.0);
	for (int ilev = 0; ilev < nlevels; ilev++) res_exact       [ilev].setVal(0.0);
	for (int ilev = 0; ilev < nlevels; ilev++) res_numeric     [ilev].setVal(0.0);
	for (int ilev = 0; ilev < nlevels; ilev++) ghost_force     [ilev].setVal(0.0);


	amrex::Vector<amrex::MultiFab> &u   = solution_exact;
	amrex::Vector<amrex::MultiFab> &v   = solution_numeric;
	amrex::Vector<amrex::MultiFab> &Au  = solution_error;
	amrex::Vector<amrex::MultiFab> &Av  = rhs_prescribed;
	amrex::Vector<amrex::MultiFab> &vAu = rhs_exact;
	amrex::Vector<amrex::MultiFab> &uAv = rhs_numeric;
	amrex::Vector<amrex::MultiFab> &diff = res_exact;

	//IC::Random ic(geom);
	std::complex<int> I(0,1);

	IC::Trig   ic(geom);//,1.0,i,i,i);
	for (int d=0; d<AMREX_SPACEDIM; d++)
	{
		ic.SetComp(d);
		AMREX_D_TERM(for (int i = 1; i < 2; i++),
			     for (int j = 1; j < 2; j++),
			     for (int k = 0; k < 1; k++))
		{
			ic.Define(Util::Random(),AMREX_D_DECL(i*I,j*I,k*I));
			for (int ilev = 0; ilev < nlevels; ++ilev)  ic.Add(ilev, u);
			ic.Define(Util::Random(),AMREX_D_DECL(i*I,j*I,k*I));
			for (int ilev = 0; ilev < nlevels; ++ilev)  ic.Add(ilev, v);
		}
	}
	for (int ilev = 0; ilev < nlevels; ++ilev)  Au[ilev].setVal(0.0);
	for (int ilev = 0; ilev < nlevels; ++ilev)  Av[ilev].setVal(0.0);
		
 	LPInfo info;
 	info.setAgglomeration(1);
 	info.setConsolidation(1);
 	info.setMaxCoarseningLevel(0);
 	nlevels = geom.size();

	::Operator::Elastic<model_type> elastic;
 	elastic.define(geom, cgrids, dmap, info);

	using bctype = ::Operator::Elastic<model_type>::BC;
	elastic.SetBC({{AMREX_D_DECL({AMREX_D_DECL(bctype::Displacement,bctype::Traction,bctype::Displacement)},
				     {AMREX_D_DECL(bctype::Traction,bctype::Displacement,bctype::Displacement)},
				     {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)})}},
		      {{AMREX_D_DECL({AMREX_D_DECL(bctype::Displacement,bctype::Traction,bctype::Displacement)},
				     {AMREX_D_DECL(bctype::Traction,bctype::Displacement,bctype::Displacement)},
				     {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)})}});


 	for (int ilev = 0; ilev < nlevels; ++ilev) elastic.SetModel(ilev,modelfab[ilev]);

	amrex::MLMG mlmg(elastic);

	
	mlmg.apply(GetVecOfPtrs(Au),GetVecOfPtrs(u));
	mlmg.apply(GetVecOfPtrs(Av),GetVecOfPtrs(v));

	for (int ilev=0; ilev < nlevels; ilev++)
	{
		Util::Message(INFO,amrex::MultiFab::Dot(u[ilev],0,Av[ilev],0,AMREX_SPACEDIM,0));
		Util::Message(INFO,amrex::MultiFab::Dot(v[ilev],0,Au[ilev],0,AMREX_SPACEDIM,0));

		amrex::MultiFab::Copy(uAv[ilev],u[ilev],0,0,AMREX_SPACEDIM,0); // uAv = u
		amrex::MultiFab::Multiply(uAv[ilev],Av[ilev],0,0,AMREX_SPACEDIM,0); // uAv *= Av

		amrex::MultiFab::Copy(vAu[ilev],v[ilev],0,0,AMREX_SPACEDIM,0); // vAu = v
		amrex::MultiFab::Multiply(vAu[ilev],Au[ilev],0,0,AMREX_SPACEDIM,0); // vAu *= Au

		amrex::MultiFab::Copy(diff[ilev],uAv[ilev],0,0,AMREX_SPACEDIM,0);
		amrex::MultiFab::Subtract(diff[ilev],vAu[ilev],0,0,AMREX_SPACEDIM,0);
	}

	if (plotfile != "")
	{
		if (verbose) Util::Message(INFO,"Printing plot file to ",plotfile);
		WritePlotFile(plotfile);
	}



	
	return failed;
}


// This is a test that solves the equation
//     Lap(u) = 0
// with u = 0 on the boundary.
// Only one component is solved for, specified by the "component" argument.
int
Elastic::TrigTest(bool verbose, int component, int n, std::string plotfile)
{
	Set::Scalar tolerance = 0.01;

	int failed = 0;

	// Define the "model" fab to be a Laplacian, so that this
	// elastic operator acts as a Laplacian on the "component-th" component of the fab.
	Set::Scalar alpha = 1.0;
	using model_type = Model::Solid::LinearElastic::Laplacian;
	model_type model(alpha);
	amrex::Vector<amrex::FabArray<amrex::BaseFab<model_type> > > modelfab(nlevels); 
 	for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev].define(ngrids[ilev], dmap[ilev], 1, 1);
 	for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev].setVal(model);

	// Initialize: set the rhs_prescribed to sin(n pi x1 / L) * sin(n pi x2 / L), so that
	// the exact solution is sin(n pi x1 / L) * sin(n pi x2 / L) / pi / 2.
	// Set everything else to zero.
	std::complex<int> i(0,1);
	IC::Trig icrhs(geom,1.0,AMREX_D_DECL(n*i,n*i,n*i),dim);
	icrhs.SetComp(component);
	//Set::Scalar dim = (Set::Scalar)(AMREX_SPACEDIM);
	IC::Trig icexact(geom,-(1./dim/Set::Constant::Pi/Set::Constant::Pi),AMREX_D_DECL(n*i,n*i,n*i),dim);
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

		// icexact.Initialize(ilev,solution_numeric); 
	}







	// Create and configure the Elastic operator
	// See:
	//   - src/Operator/Elastic.cpp   <<  the actual operator
	//   - src/Operator/Operator.cpp  <<  the abstract operator that Elastic inherits from. (where reflux is)
	amrex::LPInfo info;
 	info.setAgglomeration(1);
 	info.setConsolidation(1);
 	info.setMaxCoarseningLevel(0);
 	nlevels = geom.size();
	::Operator::Elastic<model_type> elastic;
 	elastic.define(geom, cgrids, dmap, info);
 	using bctype = ::Operator::Elastic<model_type>::BC;
 	for (int ilev = 0; ilev < nlevels; ++ilev) elastic.SetModel(ilev,modelfab[ilev]);
	if (dim == 1)
	{
		elastic.SetBC(
			      {{AMREX_D_DECL({AMREX_D_DECL(bctype::Displacement,bctype::Traction,bctype::Traction)},
					     {AMREX_D_DECL(bctype::Traction,bctype::Displacement,bctype::Traction)},
					     {AMREX_D_DECL(bctype::Traction,bctype::Traction,bctype::Displacement)})}},
			      {{AMREX_D_DECL({AMREX_D_DECL(bctype::Displacement,bctype::Traction,bctype::Traction)},
					     {AMREX_D_DECL(bctype::Traction,bctype::Displacement,bctype::Traction)},
					     {AMREX_D_DECL(bctype::Traction,bctype::Traction,bctype::Displacement)})}});

	}
	if (dim == 2)
	{
		elastic.SetBC(
			      {{AMREX_D_DECL({AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)},
					     {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)},
					     {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)})}},
			      {{AMREX_D_DECL({AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)},
					     {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)},
					     {AMREX_D_DECL(bctype::Displacement,bctype::Displacement,bctype::Displacement)})}});
	}





	// Call out specific value at a C/F point
	if (nlevels > 1)
		Util::Message(INFO,
			      "crse=",solution_numeric[0][amrex::MFIter(solution_numeric[0])](amrex::IntVect(AMREX_D_DECL(4,8,8))), " ",
			      "fine=",solution_numeric[1][amrex::MFIter(solution_numeric[1])](amrex::IntVect(AMREX_D_DECL(8,16,16))));


	// Create MLMG solver and solve
	amrex::MLMG mlmg(elastic);
	mlmg.setFixedIter(20);
	mlmg.setMaxIter(1000);
	mlmg.setMaxFmgIter(20);
 	if (verbose)
 	{
 		mlmg.setVerbose(4);
 		mlmg.setCGVerbose(0);
 	}
 	else
 	{
 		mlmg.setVerbose(0);
 		mlmg.setCGVerbose(0);
	}
 	mlmg.setBottomMaxIter(50);
 	mlmg.setFinalFillBC(false);	
 	mlmg.setBottomSolver(MLMG::BottomSolver::bicgstab);
	Set::Scalar tol_rel = 1E-8;
	Set::Scalar tol_abs = 0;
 	mlmg.solve(GetVecOfPtrs(solution_numeric), GetVecOfConstPtrs(rhs_prescribed), tol_rel,tol_abs);

	// Call out specific value at C/F point
	if (nlevels > 1)
		Util::Message(INFO,
			      "crse=",solution_numeric[0][amrex::MFIter(solution_numeric[0])](amrex::IntVect(AMREX_D_DECL(4,8,8))), " ",
			      "fine=",solution_numeric[1][amrex::MFIter(solution_numeric[1])](amrex::IntVect(AMREX_D_DECL(8,16,16))));

	// Compute solution error
	for (int i = 0; i < nlevels; i++)
	{
	 	amrex::MultiFab::Copy(solution_error[i],solution_numeric[i],component,component,1,1);
	 	amrex::MultiFab::Subtract(solution_error[i],solution_exact[i],component,component,1,1);
	}

	//Compute numerical right hand side
	mlmg.apply(GetVecOfPtrs(rhs_numeric),GetVecOfPtrs(solution_numeric));

	// Check to make sure that point didn't change
	if (nlevels > 1)
		Util::Message(INFO,
			      "crse=",solution_numeric[0][amrex::MFIter(solution_numeric[0])](amrex::IntVect(AMREX_D_DECL(4,8,8))), " ",
			      "fine=",solution_numeric[1][amrex::MFIter(solution_numeric[1])](amrex::IntVect(AMREX_D_DECL(8,16,16))));

	// Compute exact right hand side
	mlmg.apply(GetVecOfPtrs(rhs_exact),GetVecOfPtrs(solution_exact));

	// Compute numerical residual
	mlmg.compResidual(GetVecOfPtrs(res_numeric),GetVecOfPtrs(solution_numeric),GetVecOfConstPtrs(rhs_prescribed));

	// Compute exact residual
	mlmg.compResidual(GetVecOfPtrs(res_exact),GetVecOfPtrs(solution_exact),GetVecOfConstPtrs(rhs_prescribed));

	// Compute the "ghost force" that introduces the error
	//mlmg.apply(GetVecOfPtrs(ghost_force),GetVecOfPtrs(solution_error));
	
	for (int i = 0; i < nlevels; i++)
	{
	 	amrex::MultiFab::Copy(ghost_force[i],rhs_prescribed[i],component,component,1,1);
	 	amrex::MultiFab::Subtract(ghost_force[i],rhs_numeric[i],component,component,1,1);
	}

	// // Normalize so that they are all per unit area
	// for (int i = 0; i < nlevels; i++)
	// {
	// 	const Real* DX = geom[i].CellSize();
	// 	ghost_force[i].mult(1.0/DX[0]/DX[1]);
	// }	

	// Output plot file
	if (plotfile != "")
	{
		Util::Message(INFO,"Printing plot file to ",plotfile);
		WritePlotFile(plotfile);
	}

	// Find maximum solution error
	std::vector<Set::Scalar> error_norm(nlevels);
	for (int i = 0; i < nlevels; i++) error_norm[i] = solution_error[0].norm0(component,0,false) / solution_exact[0].norm0(component,0,false);
	Set::Scalar maxnorm = fabs(*std::max_element(error_norm.begin(),error_norm.end()));

	if (verbose) Util::Message(INFO,"relative error = ", 100*maxnorm, " %");

	if (maxnorm > tolerance) failed += 1;

	return failed;
}

int Elastic::RefluxTest(int verbose)
{
	int failed = 0;

	using model_type = Model::Solid::LinearElastic::Isotropic; model_type model(2.6,6.0); 

	amrex::Vector<amrex::FabArray<amrex::BaseFab<Model::Solid::LinearElastic::Isotropic> > >
		modelfab(nlevels); 
 	for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev].define(ngrids[ilev], dmap[ilev], 1, 1);
 	for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev].setVal(model);


 	LPInfo info;
 	info.setAgglomeration(1);
 	info.setConsolidation(1);
 	info.setMaxCoarseningLevel(0);
 	nlevels = geom.size();


	::Operator::Elastic<Model::Solid::LinearElastic::Isotropic> elastic;
 	elastic.define(geom, cgrids, dmap, info);
 	for (int ilev = 0; ilev < nlevels; ++ilev) elastic.SetModel(ilev,modelfab[ilev]);


	std::vector<int> comps = {1};
	std::vector<Set::Scalar> alphas = {1.0};
	std::vector<Set::Scalar> ms = {{1.0,2.0}};
	std::vector<Set::Vector> ns = {{Set::Vector(AMREX_D_DECL(1,0,0)),
					Set::Vector(AMREX_D_DECL(-1,0,0)),
					Set::Vector(AMREX_D_DECL(0,1,0)),
					Set::Vector(AMREX_D_DECL(0,-1,0)),
					Set::Vector(AMREX_D_DECL(1,1,0))}};
	std::vector<Set::Vector> bs = {{Set::Vector(AMREX_D_DECL(0,0.25,0))}};

	for (std::vector<int>::iterator comp = comps.begin(); comp != comps.end(); comp++)
	for (std::vector<Set::Scalar>::iterator m = ms.begin(); m != ms.end(); m++)
	for (std::vector<Set::Vector>::iterator n = ns.begin(); n != ns.end(); n++)
	for (std::vector<Set::Vector>::iterator b = bs.begin(); b != bs.end(); b++)
	{
		::Operator::Elastic<model_type> mlabec;
	
		mlabec.define(geom, cgrids, dmap, info);
		mlabec.setMaxOrder(2);
		for (int ilev = 0; ilev < nlevels; ++ilev) mlabec.SetModel(ilev,modelfab[ilev]);
		mlabec.SetBC({{{AMREX_D_DECL(::Operator::Elastic<model_type>::BC::Displacement,
					    ::Operator::Elastic<model_type>::BC::Displacement,
					     ::Operator::Elastic<model_type>::BC::Displacement)}}},
			{{{AMREX_D_DECL(::Operator::Elastic<model_type>::BC::Displacement,
				       ::Operator::Elastic<model_type>::BC::Displacement,
					::Operator::Elastic<model_type>::BC::Displacement)}}});


		solution_exact[0].setVal(0.0);
		solution_exact[1].setVal(0.0);
		rhs_prescribed[0].setVal(0.0);
		rhs_prescribed[1].setVal(0.0);

		IC::Affine ic(geom,*n,1.0,*b,true,*m);

		ic.SetComp(*comp);
		ic.Initialize(0,solution_exact);
		ic.Initialize(1,solution_exact);

		mlabec.FApply(0,0,rhs_prescribed[0],solution_exact[0]);
		mlabec.FApply(1,0,rhs_prescribed[1],solution_exact[1]);

		res_numeric[0].setVal(0.0);
		res_numeric[1].setVal(0.0);

		mlabec.BuildMasks();
		mlabec.Reflux(0,
		 	      res_numeric[0], solution_exact[0], rhs_prescribed[0],
		 	      res_numeric[1], solution_exact[1], rhs_prescribed[1]);

		
		Set::Scalar residual = res_numeric[0].norm0();

		if (rhs_prescribed[0].norm0() > 1E-15) residual /= rhs_prescribed[0].norm0();

		std::stringstream ss;
		ss << "n=["<< std::setw(5) << (n->transpose()) << "], "
		   << "b=["<< std::setw(5) << (b->transpose()) << "], "
		   << "m=" << *m<<", "
		   << "comp="<<(*comp);

		bool pass = fabs(residual) < 1E-12;

		if (verbose > 0) Util::Test::SubMessage(ss.str(), !pass);

		if (!pass) failed++;
	}
	
	return failed;

}


int SpatiallyVaryingCTest(int /*verbose*/)
{
	/*
	using model_type = Model::Solid::LinearElastic::Degradable::Isotropic;
	model_type model(400,300);

	amrex::Vector<amrex::Geometry>		geom;
	amrex::Vector<amrex::BoxArray> 		ngrids,cgrids;
	amrex::Vector<amrex::DistributionMapping>	dmap;

	amrex::Vector<amrex::MultiFab> 	u;
	amrex::Vector<amrex::MultiFab> 	rhs;
	amrex::Vector<amrex::MultiFab> 	res;
	amrex::Vector<amrex::MultiFab> 	stress;
	amrex::Vector<amrex::MultiFab> 	strain;
	amrex::Vector<amrex::MultiFab> 	energy;

	amrex::Vector<amrex::FabArray<amrex::BaseFab<model_type> > > modelfab;
	std::map<std::string,::Operator::Elastic<model_type>::BC > 			bc_map;
	amrex::Vector<std::string> AMREX_D_DECL(bc_x_lo_str,bc_y_lo_str,bc_z_lo_str);
	amrex::Vector<std::string> AMREX_D_DECL(bc_x_hi_str,bc_y_hi_str,bc_z_hi_str);
	std::array<::Operator::Elastic<model_type>::BC,AMREX_SPACEDIM> AMREX_D_DECL(bc_x_lo,bc_y_lo,bc_z_lo);
	std::array<::Operator::Elastic<model_type>::BC,AMREX_SPACEDIM> AMREX_D_DECL(bc_x_hi,bc_y_hi,bc_z_hi);

	bc_map["displacement"] = 	::Operator::Elastic<model_type>::BC::Displacement;
	bc_map["disp"] = 		::Operator::Elastic<model_type>::BC::Displacement;
	bc_map["traction"] = 		::Operator::Elastic<model_type>::BC::Traction;
	bc_map["trac"] = 		::Operator::Elastic<model_type>::BC::Traction;
	bc_map["periodic"] = 		::Operator::Elastic<model_type>::BC::Periodic;

	AMREX_D_TERM(	bc_x_lo_str = {AMREX_D_DECL("disp", "disp", "disp")};
					bc_x_hi_str = {AMREX_D_DECL("disp", "disp", "disp")};
					,
					bc_y_lo_str = {AMREX_D_DECL("trac", "trac", "trac")};
					bc_y_hi_str = {AMREX_D_DECL("trac", "trac", "trac")};
					,
					bc_z_lo_str = {AMREX_D_DECL("trac", "trac", "trac")};
					bc_z_hi_str = {AMREX_D_DECL("trac", "trac", "trac")};
				);

	AMREX_D_TERM(	bc_x_lo = {AMREX_D_DECL(bc_map[bc_x_lo_str[0]], bc_map[bc_x_lo_str[1]], bc_map[bc_x_lo_str[2]])};
					bc_x_hi = {AMREX_D_DECL(bc_map[bc_x_hi_str[0]], bc_map[bc_x_hi_str[1]], bc_map[bc_x_hi_str[2]])};
					,
					bc_y_lo = {AMREX_D_DECL(bc_map[bc_y_lo_str[0]], bc_map[bc_y_lo_str[1]], bc_map[bc_y_lo_str[2]])};
					bc_y_hi = {AMREX_D_DECL(bc_map[bc_y_hi_str[0]], bc_map[bc_y_hi_str[1]], bc_map[bc_y_hi_str[2]])};
					,
					bc_z_lo = {AMREX_D_DECL(bc_map[bc_z_lo_str[0]], bc_map[bc_z_lo_str[1]], bc_map[bc_z_lo_str[2]])};
					bc_z_hi = {AMREX_D_DECL(bc_map[bc_z_hi_str[0]], bc_map[bc_z_hi_str[1]], bc_map[bc_z_hi_str[2]])};);

	Set::Vector AMREX_D_DECL(bc_left,bc_bottom,bc_back);
	Set::Vector AMREX_D_DECL(bc_right,bc_top,bc_front);

	AMREX_D_TERM(	bc_left = {AMREX_D_DECL(0.,0.,0.)};
					bc_right = {AMREX_D_DECL(1.,0.,0.)};,
					bc_bottom = {AMREX_D_DECL(0.,0.,0.)};
					bc_top = {AMREX_D_DECL(0.,0.,0.)};,
					bc_back = {AMREX_D_DECL(0.,0.,0.)};
					bc_front = {AMREX_D_DECL(0.,0.,0.)};);

	
	int n_cell = 16;
 	int nlevels = 1;
	int ref_ratio = 2;
	int max_grid_size = 100000;
	Set::Vector body_force(AMREX_D_DECL(0.,-0.0001,0.));

 	
 	geom.resize(nlevels);	

 	ngrids.resize(nlevels);
 	cgrids.resize(nlevels);
 	dmap.resize(nlevels);

 	u.resize(nlevels);
 	res.resize(nlevels);
 	rhs.resize(nlevels);
 	stress.resize(nlevels);
 	strain.resize(nlevels);
 	energy.resize(nlevels);
	modelfab.resize(nlevels);

	amrex::RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
	amrex::Geometry::Setup(&rb, 0);

	amrex::Box NDomain(amrex::IntVect{AMREX_D_DECL(0,0,0)},
			   amrex::IntVect{AMREX_D_DECL(n_cell,n_cell,n_cell)},
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
 		cdomain.refine(ref_ratio); 

		ngrids[ilev] = cgrids[ilev];
		ngrids[ilev].convert(amrex::IntVect::TheNodeVector());
	}

	int number_of_components = AMREX_SPACEDIM;
	int number_of_stress_components = AMREX_SPACEDIM*AMREX_SPACEDIM;
 	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
 		dmap   [ilev].define(cgrids[ilev]);
 		u       [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 1); 
 		res     [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 1); 
 		rhs     [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 1);
 		stress     [ilev].define(ngrids[ilev], dmap[ilev], number_of_stress_components, 1);
 		strain     [ilev].define(ngrids[ilev], dmap[ilev], number_of_stress_components, 1);
 		energy     [ilev].define(ngrids[ilev], dmap[ilev], 1, 1);
		modelfab[ilev].define(ngrids[ilev], dmap[ilev], 1, 1);
 	}

 	for (int ilev = 0; ilev < nlevels; ++ilev)
 	{
 		u[ilev].setVal(0.0);
 		res[ilev].setVal(0.0);
 		rhs[ilev].setVal(0.0);
 		stress[ilev].setVal(0.0);
 		strain[ilev].setVal(0.0);
 		energy[ilev].setVal(0.0);
		modelfab[ilev].setVal(model);
 	}

 	LPInfo info;
 	info.setAgglomeration(1);
 	info.setConsolidation(1);
 	info.setMaxCoarseningLevel(0);
 	nlevels = geom.size();

 	std::vector<Set::Vector> normal = {{Set::Vector(AMREX_D_DECL(0.,1.,0.))}};
 	std::vector<Set::Vector> point = {{Set::Vector(AMREX_D_DECL(0.95,0.95,0.95))}};
 	std::vector<Set::Scalar> value = {0.1};
 	std::vector<Set::Scalar> exponent = {0.};

 	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		for (amrex::MFIter mfi(modelfab[ilev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.tilebox();
	 		amrex::BaseFab<model_type> &modelbox = (modelfab[ilev])[mfi];
	 		AMREX_D_TERM(for (int i = box.loVect()[0]-1; i<=box.hiVect()[0]+1; i++),
		 		     for (int j = box.loVect()[1]-1; j<=box.hiVect()[1]+1; j++),
		 		     for (int k = box.loVect()[2]-1; k<=box.hiVect()[2]+1; k++))
	 		{
	 			amrex::IntVect m(AMREX_D_DECL(i,j,k));
	 			std::vector<Set::Vector> corners = {	{AMREX_D_DECL(geom[ilev].ProbLo()[0],geom[ilev].ProbLo()[1],geom[ilev].ProbLo()[2])},
	 													{AMREX_D_DECL(geom[ilev].ProbHi()[0],geom[ilev].ProbLo()[1],geom[ilev].ProbLo()[2])}
	 													#if AMREX_SPACEDIM > 1
	 													,{AMREX_D_DECL(geom[ilev].ProbLo()[0],geom[ilev].ProbHi()[1],geom[ilev].ProbLo()[2])}
	 													,{AMREX_D_DECL(geom[ilev].ProbHi()[0],geom[ilev].ProbHi()[1],geom[ilev].ProbLo()[2])}
	 													#if AMREX_SPACEDIM > 2
	 													,{AMREX_D_DECL(geom[ilev].ProbLo()[0],geom[ilev].ProbLo()[1],geom[ilev].ProbHi()[2])}
	 													,{AMREX_D_DECL(geom[ilev].ProbHi()[0],geom[ilev].ProbLo()[1],geom[ilev].ProbHi()[2])}
	 													,{AMREX_D_DECL(geom[ilev].ProbLo()[0],geom[ilev].ProbHi()[1],geom[ilev].ProbHi()[2])}
	 													,{AMREX_D_DECL(geom[ilev].ProbHi()[0],geom[ilev].ProbHi()[1],geom[ilev].ProbHi()[2])}
	 													#endif
	 													#endif
	 													};
	 			Set::Scalar max_dist = 0.0;
	 			AMREX_D_TERM(	amrex::Real x1 = geom[ilev].ProbLo()[0] + ((amrex::Real)(i)) * geom[ilev].CellSize()[0];,
					     		amrex::Real x2 = geom[ilev].ProbLo()[1] + ((amrex::Real)(j)) * geom[ilev].CellSize()[1];,
					     		amrex::Real x3 = geom[ilev].ProbLo()[2] + ((amrex::Real)(k)) * geom[ilev].CellSize()[2];);
	 			
	 			Set::Vector x(AMREX_D_DECL(x1,x2,x3));

	 			for (std::vector<Set::Vector>::iterator corner = corners.begin(); corner != corners.end(); corner++)
	 				max_dist = std::max (max_dist, (*corner-point[0]).dot(normal[0]));
	 			
	 			Set::Scalar dist = (x-point[0]).dot(normal[0]);
	 			if(dist >= 0.0)
	 				modelbox(m).DegradeModulus({value[0]*std::pow((dist/max_dist),exponent[0])});

	 		}
		}
	}

	Elastic<model_type> mlabec;
	
	mlabec.define(geom, cgrids, dmap, info);
	mlabec.setMaxOrder(2);
	mlabec.SetBC({{AMREX_D_DECL(bc_x_lo,bc_y_lo,bc_z_lo)}},
		     				{{AMREX_D_DECL(bc_x_hi,bc_y_hi,bc_z_hi)}});

	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		const Real* DX = geom[ilev].CellSize();
		Set::Scalar volume = AMREX_D_TERM(DX[0],*DX[1],*DX[2]);

		mlabec.SetModel(ilev,modelfab[ilev]);
		AMREX_D_TERM(rhs[ilev].setVal(body_force[0]*volume,0,1);,
					rhs[ilev].setVal(body_force[1]*volume,1,1);,
					rhs[ilev].setVal(body_force[2]*volume,2,1););
		for (amrex::MFIter mfi(rhs[ilev],true); mfi.isValid(); ++mfi)
		{
		 	const amrex::Box& box = mfi.tilebox();

		 	amrex::BaseFab<amrex::Real> &rhsfab = (rhs[ilev])[mfi];

		 	AMREX_D_TERM(for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++),
		 		     for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++),
		 		     for (int k = box.loVect()[2]; k<=box.hiVect()[2]; k++))
		 	{
		 		bool AMREX_D_DECL(xmin = false, ymin = false, zmin = false);
				bool AMREX_D_DECL(xmax = false, ymax = false, zmax = false);

		 		AMREX_D_TERM(	xmin = (i == geom[ilev].Domain().loVect()[0]);
		 						xmax = (i == geom[ilev].Domain().hiVect()[0]+1);
		 						,
		 						ymin = (j == geom[ilev].Domain().loVect()[1]);
		 						ymax = (j == geom[ilev].Domain().hiVect()[1]+1);
		 						,
		 						zmin = (k == geom[ilev].Domain().loVect()[2]);
		 						zmax = (k == geom[ilev].Domain().hiVect()[2]+1););
		 	
		 
		 		if (false || AMREX_D_TERM(xmin || xmax, || ymin || ymax, || zmin || zmax))
		 		{
		 			AMREX_D_TERM(	rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = 0.0;,
 									rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),1) = 0.0;,
 									rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),2) = 0.0;);
		 			for(int l = 0; l<AMREX_SPACEDIM; l++)
					{
						AMREX_D_TERM(
							if(xmin && bc_x_lo[l]==Elastic<model_type>::BC::Displacement)
								rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = bc_left[l];
							if(xmax && bc_x_hi[l]==Elastic<model_type>::BC::Displacement)
								rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = bc_right[l];
							,
							if(ymin && bc_y_lo[l]==Elastic<model_type>::BC::Displacement)
								rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = bc_bottom[l];
							if(ymax && bc_y_hi[l]==Elastic<model_type>::BC::Displacement)
								rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = bc_top[l];
							,
							if(zmin && bc_z_lo[l]==Elastic<model_type>::BC::Displacement)
								rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = bc_back[l];
							if(zmax && bc_z_hi[l]==Elastic<model_type>::BC::Displacement)
								rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = bc_front[l];
							);
					}
		 		}
		 	}
		}
	}

	amrex::MLMG solver(mlabec);
	solver.setMaxIter(10000);
	solver.setMaxFmgIter(10000);
	solver.setVerbose(4);
	solver.setCGVerbose(4);
	solver.setBottomMaxIter(10000);
	solver.setBottomSolver(MLMG::BottomSolver::bicgstab);
	
	solver.solve(GetVecOfPtrs(u),
		GetVecOfConstPtrs(rhs),
		1E-8,
		1E-8);

	for (int lev = 0; lev < nlevels; lev++)
	{
		mlabec.Strain(lev,strain[lev],u[lev]);
		mlabec.Stress(lev,stress[lev],u[lev]);
		mlabec.Energy(lev,energy[lev],u[lev]);
	}
	*/
	Util::Abort(INFO,"Test is not fully implemente yet");
	return 1;
}


}
}
	     

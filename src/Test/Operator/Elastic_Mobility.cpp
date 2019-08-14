#include "Test/Operator/Elastic.H"
#include "Model/Solid/LinearElastic/Isotropic.H"
#include "Model/Solid/LinearElastic/Cubic.H"
#include "Model/Solid/LinearElastic/MultiWell.H"
#include "IC/Affine.H"
#include "IC/Trig.H"
#include "Operator/Elastic2.H"
#include "Solver/Nonlocal/Linear.H"

namespace Test
{
namespace Operator
{
int Elastic::Mobility(int verbose, int component, std::string plotfile)
{/*
	Generate();

        Set::Scalar tolerance = 0.01;

	int failed = 0;

	using model_type = Model::Solid::LinearElastic::Multiwell;
	Set::Scalar lame=2.6, shear=6.0;
	model_type model(lame,shear);//(lame, shear);
    model.Randomize();
	//Use this instead to run for Cubic elastic case.
	//using model_type = Model::Solid::LinearElastic::Cubic; model_type model; model.Randomize(); 

	amrex::Vector<amrex::FabArray<amrex::BaseFab<model_type> > > modelfab(nlevels); 

 	for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev].define(ngrids[ilev], dmap[ilev], 1, 2);
 	for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev].setVal(model);

	Set::Vector vec = Set::Vector::Zero();
	vec[component]=0.1;

	IC::Affine icrhs(geom, {0,1}, 0.2, {0.0,0.5} , true, 1.0);
	icrhs.SetComp(0);

	IC::Affine icexact(geom, vec, 1.0, Set::Vector::Zero(), false, 1.0);
	icexact.SetComp(component);
	
	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		icrhs.Initialize(ilev,rhs_prescribed);
		icexact.Initialize(ilev,solution_exact);
	}

	amrex::LPInfo info;
 	info.setAgglomeration(m_agglomeration);
 	info.setConsolidation(m_consolidation);
 	if (m_maxCoarseningLevel > -1) info.setMaxCoarseningLevel(m_maxCoarseningLevel);
 	nlevels = geom.size();

	::Operator::Elastic2<model_type> elastic;
	elastic.SetHomogeneous(false);
 	elastic.define(geom, cgrids, dmap, info);
	for (int ilev = 0; ilev < nlevels; ++ilev) elastic.SetModel(ilev,modelfab[ilev]);
	
	BC::Operator::Elastic<model_type> bc;
    bc.Set(bc.Face::XLO, bc.Direction::X, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
    bc.Set(bc.Face::XLO, bc.Direction::Y, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
    bc.Set(bc.Face::XHI, bc.Direction::X, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
    bc.Set(bc.Face::XHI, bc.Direction::Y, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
    bc.Set(bc.Face::YLO, bc.Direction::X, bc.Type::Displacement, 0.0, rhs_prescribed, geom);
    bc.Set(bc.Face::YLO, bc.Direction::Y, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
    bc.Set(bc.Face::YHI, bc.Direction::X, bc.Type::Displacement, 0.0, rhs_prescribed, geom);
    bc.Set(bc.Face::YHI, bc.Direction::Y, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
	elastic.SetBC(&bc);

	Solver::Nonlocal::Linear mlmg(elastic);
	if (m_fixedIter > -1)     mlmg.setFixedIter(m_fixedIter);
	if (m_maxIter > -1 )      mlmg.setMaxIter(m_maxIter);
	if (m_maxFmgIter > -1)    mlmg.setMaxFmgIter(m_maxFmgIter);
 	mlmg.setVerbose(verbose);
 	if (m_bottomMaxIter > -1) mlmg.setBottomMaxIter(m_bottomMaxIter);
	
	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		icrhs.Initialize(ilev,rhs_exact);
		icexact.Initialize(ilev,solution_exact);
	}

	mlmg.apply(GetVecOfPtrs(rhs_prescribed),GetVecOfPtrs(rhs_exact));
	bc.Set(bc.Face::XLO, bc.Direction::X, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
    bc.Set(bc.Face::XLO, bc.Direction::Y, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
    bc.Set(bc.Face::XHI, bc.Direction::X, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
    bc.Set(bc.Face::XHI, bc.Direction::Y, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
    bc.Set(bc.Face::YLO, bc.Direction::X, bc.Type::Displacement, 0.0, rhs_prescribed, geom);
    bc.Set(bc.Face::YLO, bc.Direction::Y, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
    bc.Set(bc.Face::YHI, bc.Direction::X, bc.Type::Displacement, 0.1, rhs_prescribed, geom);
    bc.Set(bc.Face::YHI, bc.Direction::Y, bc.Type::Neumann, 0.0, rhs_prescribed, geom);

	
    if (component!=0) component-=component;
         
 	mlmg.solve(GetVecOfPtrs(solution_numeric),
		   GetVecOfConstPtrs(rhs_prescribed),
		   m_tol_rel,m_tol_abs);

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

	return failed;*/

}
}
}

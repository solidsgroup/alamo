#include "Test/Operator/Elastic.H"
#include "Model/Solid/LinearElastic/Isotropic.H"
#include "IC/Affine.H"
#include "IC/Trig.H"
#include "Operator/Elastic.H"
#include "Solver/Linear.H"

namespace Test
{
namespace Operator
{
int Elastic::UniaxialTest(int verbose, int component, std::string plotfile)
{
	Generate();

        Set::Scalar tolerance = 0.01;

	int failed = 0;

	using model_type = Model::Solid::LinearElastic::Isotropic;
	Set::Scalar lame=2.6, shear=6.0;
	model_type model(lame, shear);

	amrex::Vector<amrex::FabArray<amrex::BaseFab<model_type> > > modelfab(nlevels); 

 	for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev].define(ngrids[ilev], dmap[ilev], 1, 2);
 	for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev].setVal(model);

	Set::Vector vec = Set::Vector::Zero();
	vec[component]=0.1;

	IC::Affine icrhs(geom, Set::Vector::Zero(), 1.0, Set::Vector::Zero() , false, 1.0);
	icrhs.SetComp(component);

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

	::Operator::Elastic<model_type> elastic;
	elastic.SetHomogeneous(false);
 	elastic.define(geom, cgrids, dmap, info);
	for (int ilev = 0; ilev < nlevels; ++ilev) elastic.SetModel(ilev,modelfab[ilev]);
	BC::Operator::Elastic<model_type> bc;

	
	if (component==0)
	{
	  AMREX_D_TERM(,
	           bc.Set(bc.Face::XHI, bc.Direction::X, bc.Type::Displacement, 0.1, rhs_prescribed, geom);
		       bc.Set(bc.Face::XHI, bc.Direction::Y, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
		       bc.Set(bc.Face::YLO, bc.Direction::X, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
		       bc.Set(bc.Face::YHI, bc.Direction::X, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
		       ,
		       bc.Set(bc.Face::XHI, bc.Direction::Z, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
		       bc.Set(bc.Face::ZLO, bc.Direction::X, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
		       bc.Set(bc.Face::ZHI, bc.Direction::X, bc.Type::Neumann, 0.0, rhs_prescribed, geom););
	}
	
#if AMREX_SPACEDIM > 1	
	else if (component==1)
	{
	  AMREX_D_TERM(,
		       bc.Set(bc.Face::XLO, bc.Direction::Y, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
		       bc.Set(bc.Face::XHI, bc.Direction::Y, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
		       bc.Set(bc.Face::YHI, bc.Direction::X, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
		       bc.Set(bc.Face::YHI, bc.Direction::Y, bc.Type::Displacement, 0.1, rhs_prescribed, geom);
		       ,
		       bc.Set(bc.Face::YHI, bc.Direction::Z, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
		       bc.Set(bc.Face::ZLO, bc.Direction::Y, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
		       bc.Set(bc.Face::ZHI, bc.Direction::Y, bc.Type::Neumann, 0.0, rhs_prescribed, geom););
	}
#endif

#if AMREX_SPACEDIM > 2
	if (component==2) 
	{
	  bc.Set(bc.Face::XLO, bc.Direction::Z, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
	  bc.Set(bc.Face::XHI, bc.Direction::Z, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
	  bc.Set(bc.Face::YLO, bc.Direction::Z, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
	  bc.Set(bc.Face::YHI, bc.Direction::Z, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
	  bc.Set(bc.Face::ZHI, bc.Direction::X, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
	  bc.Set(bc.Face::ZHI, bc.Direction::Y, bc.Type::Neumann, 0.0, rhs_prescribed, geom);
	  bc.Set(bc.Face::ZHI, bc.Direction::Z, bc.Type::Displacement, 0.1, rhs_prescribed, geom);
	}
#endif	
	
	elastic.SetBC(&bc);

	Solver::Linear mlmg(elastic);
	if (m_fixedIter > -1)     mlmg.setFixedIter(m_fixedIter);
	if (m_maxIter > -1 )      mlmg.setMaxIter(m_maxIter);
	if (m_maxFmgIter > -1)    mlmg.setMaxFmgIter(m_maxFmgIter);
 	mlmg.setVerbose(verbose);
 	if (m_bottomMaxIter > -1) mlmg.setBottomMaxIter(m_bottomMaxIter);
	
  if (component!=0)
	{
		component-=component;
	}
         
 	mlmg.solve(GetVecOfPtrs(solution_numeric),
		   GetVecOfConstPtrs(rhs_prescribed),
		   m_tol_rel,m_tol_abs);

	// Compute solution error
	for (int i = 0; i < nlevels; i++)
	{
	    amrex::MultiFab::Copy(solution_error[i],solution_numeric[i],component,component,AMREX_SPACEDIM,1);
		amrex::MultiFab::Subtract(solution_error[i],solution_exact[i],component,component,AMREX_SPACEDIM,1);
	}
	

	// Compute numerical right hand side
	mlmg.apply(GetVecOfPtrs(rhs_numeric),GetVecOfPtrs(solution_numeric));

	// Compute exact right hand side
	mlmg.apply(GetVecOfPtrs(rhs_exact),GetVecOfPtrs(solution_exact));

	
	// Compute the numeric residual
	for (int i = 0; i < nlevels; i++)
	{
		amrex::MultiFab::Copy(res_numeric[i],rhs_numeric[i],component,component,AMREX_SPACEDIM,0);
		amrex::MultiFab::Subtract(res_numeric[i],rhs_prescribed[i],component,component,AMREX_SPACEDIM,0);
	}
	for (int ilev = nlevels-1; ilev > 0; ilev--)
	 	elastic.Reflux(0,
	 		       res_numeric[ilev-1], solution_numeric[ilev-1], rhs_prescribed[ilev-1],
	 		       res_numeric[ilev],   solution_numeric[ilev],   rhs_prescribed[ilev]);

	// Compute the exact residual
	for (int i = 0; i < nlevels; i++)
	{
		amrex::MultiFab::Copy(res_exact[i],rhs_exact[i],component,component,AMREX_SPACEDIM,0);
		amrex::MultiFab::Subtract(res_exact[i],rhs_prescribed[i],component,component,AMREX_SPACEDIM,0);
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

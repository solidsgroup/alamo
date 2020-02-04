#include <AMReX_MLMG.H>

#include "Operator/Elastic.H"
#include "IC/Affine.H"
#include "Model/Solid/Linear/Isotropic.H"
#include "Elastic.H"

namespace Test
{
namespace Operator
{
int Elastic::RefluxTest(int verbose)
{
	Generate();
	int failed = 0;

	int nghost = 2;

	using model_type = Model::Solid::Linear::Isotropic; model_type model(2.6,6.0); 

	Set::Field<Model::Solid::Linear::Isotropic> modelfab(nlevels,ngrids,dmap,1,nghost); 
 	for (int ilev = 0; ilev < nlevels; ++ilev) modelfab[ilev]->setVal(model);


 	LPInfo info;
 	info.setAgglomeration(1);
 	info.setConsolidation(1);
 	info.setMaxCoarseningLevel(0);
 	nlevels = geom.size();



	::Operator::Elastic<Model::Solid::Linear::Isotropic> elastic;
 	elastic.define(geom, cgrids, dmap, info);

 	elastic.SetModel(modelfab);

	BC::Operator::Elastic<model_type> bc;
	elastic.SetBC(&bc);


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
		mlabec.SetModel(modelfab);

		BC::Operator::Elastic<model_type> bc;
		mlabec.SetBC(&bc);


		solution_exact[0]->setVal(0.0);
		solution_exact[1]->setVal(0.0);
		rhs_prescribed[0]->setVal(0.0);
		rhs_prescribed[1]->setVal(0.0);

		IC::Affine ic(geom,*n,1.0,*b,true,*m);

		ic.SetComp(*comp);
		ic.Initialize(0,solution_exact);
		ic.Initialize(1,solution_exact);

		mlabec.Fapply(0,0,*rhs_prescribed[0],*solution_exact[0]);
		mlabec.Fapply(1,0,*rhs_prescribed[1],*solution_exact[1]);

		res_numeric[0]->setVal(0.0);
		res_numeric[1]->setVal(0.0);

		mlabec.buildMasks();
		mlabec.reflux(0,
		 	      *res_numeric[0], *solution_exact[0], *rhs_prescribed[0],
		 	      *res_numeric[1], *solution_exact[1], *rhs_prescribed[1]);

		
		Set::Scalar residual = res_numeric[0]->norm0();

		if (rhs_prescribed[0]->norm0() > 1E-15) residual /= rhs_prescribed[0]->norm0();

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
}
}

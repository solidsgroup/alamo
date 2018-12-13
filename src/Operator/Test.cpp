#include <AMReX.H>

#include "Operator/Elastic.H"
#include "Model/Solid/LinearElastic/Isotropic.H"
#include "Test.H"
#include "IC/Affine.H"
#include "IC/Trig.H"
#include "IC/Random.H"

namespace Operator
{

template<>
int Test<Elastic<Model::Solid::LinearElastic::Isotropic> >::RefluxTest(int verbose)
{
	int failed = 0;

	using model_type = Model::Solid::LinearElastic::Isotropic; model_type model(2.6,6.0); 

 	amrex::Vector<amrex::Geometry> 			geom;
 	amrex::Vector<amrex::BoxArray> 			cgrids, ngrids;
 	amrex::Vector<amrex::DistributionMapping>       dmap;

 	amrex::Vector<amrex::MultiFab>  u;
 	amrex::Vector<amrex::MultiFab>  res;
 	amrex::Vector<amrex::MultiFab>  rhs;

	amrex::Vector<amrex::FabArray<amrex::BaseFab<model_type> > > modelfab;

	int n_cell = 16;
 	int nlevels = 2;
	int ref_ratio = 2;
	int max_grid_size = 100000;
	std::string orientation = "h";
 	geom.resize(nlevels);
 	cgrids.resize(nlevels);
 	ngrids.resize(nlevels);
 	dmap.resize(nlevels);

 	u.resize(nlevels);
 	res.resize(nlevels);
 	rhs.resize(nlevels);
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
 			if (orientation == "h")
 				cdomain.grow(amrex::IntVect(AMREX_D_DECL(0,-n_cell/4,0))); 
 			else if (orientation == "v")
 				cdomain.grow(amrex::IntVect(AMREX_D_DECL(-n_cell/4,0,0))); 

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
			modelfab[ilev].define(ngrids[ilev], dmap[ilev], 1, 1);
 		}

 	for (int ilev = 0; ilev < nlevels; ++ilev)
 		{
 			u[ilev].setVal(0.0);
 			res[ilev].setVal(0.0);
 			rhs[ilev].setVal(0.0);
			modelfab[ilev].setVal(model);
 		}

 	LPInfo info;
 	info.setAgglomeration(1);
 	info.setConsolidation(1);
 	info.setMaxCoarseningLevel(0);
 	nlevels = geom.size();



 	// int test_on = 5;
 	// std::string testtype = "random"; 
 	// amrex::Vector<amrex::Real> tmpn; 
 	// amrex::Vector<amrex::Real> tmpb; 
 	// amrex::Real tmpalpha; 
 	// amrex::Real tmpm = 1.0; 
 	// int comp; 

	// Set::Scalar alpha = tmpalpha, m = tmpm;
	// //Set::Vector n(1.0,1.0); Set::Vector b(0.5,0.0);
	// Set::Vector n(tmpn[0],tmpn[1]);
	// Set::Vector b(tmpb[0],tmpb[1]);

	// Util::Message(INFO,"alpha=",alpha);
	// Util::Message(INFO,"n=",n.transpose());
	// Util::Message(INFO,"b=",b.transpose());
	// Util::Message(INFO,"comp=",comp);

	std::vector<int> comps = {{0,1}};
	std::vector<Set::Scalar> alphas = {{1.0}};
	std::vector<Set::Scalar> ms = {{1.0,2.0,3.0}};
	std::vector<Set::Vector> ns = {{Set::Vector(AMREX_D_DECL(1,0,0)),
					Set::Vector(AMREX_D_DECL(-1,0,0)),
					Set::Vector(AMREX_D_DECL(0,1,0)),
					Set::Vector(AMREX_D_DECL(0,-1,0)),
					Set::Vector(AMREX_D_DECL(1,1,0))}};
	std::vector<Set::Vector> bs = {{Set::Vector(AMREX_D_DECL(0,0.25,0)),
					Set::Vector(AMREX_D_DECL(0,0.75,0)),
					Set::Vector(AMREX_D_DECL(0.5,0.5,0))}};

	for (std::vector<int>::iterator comp = comps.begin(); comp != comps.end(); comp++)
	for (std::vector<Set::Scalar>::iterator m = ms.begin(); m != ms.end(); m++)
	for (std::vector<Set::Vector>::iterator n = ns.begin(); n != ns.end(); n++)
	for (std::vector<Set::Vector>::iterator b = bs.begin(); b != bs.end(); b++)
	{
		Elastic<model_type> mlabec;
	
		mlabec.define(geom, cgrids, dmap, info);
		mlabec.setMaxOrder(2);
		for (int ilev = 0; ilev < nlevels; ++ilev) mlabec.SetModel(ilev,modelfab[ilev]);
		mlabec.SetBC({{AMREX_D_DECL(Elastic<model_type>::BC::Displacement,
					    Elastic<model_type>::BC::Displacement,
					    Elastic<model_type>::BC::Displacement)}},
			{{AMREX_D_DECL(Elastic<model_type>::BC::Displacement,
				       Elastic<model_type>::BC::Displacement,
				       Elastic<model_type>::BC::Displacement)}});


		u[0].setVal(0.0);
		u[1].setVal(0.0);
		rhs[0].setVal(0.0);
		rhs[1].setVal(0.0);

		IC::Affine ic(geom,*n,1.0,*b,true,*m);

		ic.SetComp(*comp);
		ic.Initialize(0,u);
		ic.Initialize(1,u);

		mlabec.FApply(0,0,rhs[0],u[0]);
		mlabec.FApply(1,0,rhs[1],u[1]);

		res[0].setVal(0.0);
		res[1].setVal(0.0);

		mlabec.BuildMasks();
		mlabec.Reflux(0,
		 	      res[0], u[0], rhs[0],
		 	      res[1], u[1], rhs[1]);

		
		Set::Scalar residual = res[0].norm0();

		if (rhs[0].norm0() > 1E-15) residual /= rhs[0].norm0();

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

	// if (testtype=="affine")
	// {
	// 	Util::Message(INFO,"affine, comp = ", comp);
	// 	IC::Affine ic(geom,n,alpha,b,true,m);
	// 	if (comp == 0) {ic.SetComp(0); ic.Initialize(0,u); ic.Initialize(1,u);}
	// 	if (comp == 1) {ic.SetComp(1); ic.Initialize(0,u); ic.Initialize(1,u);}
	// }
	// else if (testtype=="trig")
	// {
	// 	Util::Message(INFO,"trig");
	// 	IC::Trig ic(geom);
	// 	if (comp == 0) {ic.SetComp(0); ic.Initialize(0,u); ic.Initialize(1,u);}
	// 	if (comp == 1) {ic.SetComp(1); ic.Initialize(0,u); ic.Initialize(1,u);}
	// }
	// else if (testtype=="random")
	// {
	// 	Util::Message(INFO,"random");
	// 	IC::Random ic(geom);
	// 	if (comp == 0) {ic.SetComp(0); ic.Initialize(0,u); ic.Initialize(1,u);}
	// 	if (comp == 1) {ic.SetComp(1); ic.Initialize(0,u); ic.Initialize(1,u);}
	// }
	// else
	// 	Util::Abort(INFO,"invalid test type");



	// mlabec.FApply(0,0,rhs[0],u[0]);
	// mlabec.FApply(1,0,rhs[1],u[1]);

	// res[0].setVal(0.0);
	// res[1].setVal(0.0);

	// mlabec.BuildMasks();
	// mlabec.Reflux(0,
	// 	      res[0], u[0], rhs[0],
	// 	      res[1], u[1], rhs[1]);
}


}

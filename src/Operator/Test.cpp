#include <AMReX.H>

#include "Test.H"
#include "Operator/Elastic.H"
#include "Model/Solid/LinearElastic/Isotropic.H"


namespace Operator
{

int Test::RefluxTest(int verbose)
{
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
	amrex::Box cdomain = CDomain, ndomain = NDomain;
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
 			u       [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 0); 
 			res     [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 0); 
 			rhs     [ilev].define(ngrids[ilev], dmap[ilev], number_of_components, 0);
 		}

 	for (int ilev = 0; ilev < nlevels; ++ilev)
 		{
 			u[ilev].setVal(0.0);
 			res[ilev].setVal(0.0);
 			rhs[ilev].setVal(0.0);
 		}

 	LPInfo info;
 	info.setAgglomeration(1);
 	info.setConsolidation(1);
 	info.setMaxCoarseningLevel(0);
 	nlevels = geom.size();

 	Elastic<model_type> mlabec;

 	mlabec.define(geom, cgrids, dmap, info);
 	mlabec.setMaxOrder(2);
 	for (int ilev = 0; ilev < nlevels; ++ilev) mlabec.SetModel(ilev,modelfab[ilev]);
	// mlabec.SetBC({{AMREX_D_DECL(Elastic<model_type>::BC::Dirichlet,
	// 			    Elastic<model_type>::BC::Dirichlet,
	// 			    Elastic<model_type>::BC::Dirichlet)}},
	// 	{{AMREX_D_DECL(Elastic<model_type>::BC::Dirichlet,
	// 		       Elastic<model_type>::BC::Dirichlet,
	// 		       Elastic<model_type>::BC::Dirichlet)}});
	
	

// 	int test_on = 5; pptest.query("on",test_on);
// 	std::string testtype = "random"; pptest.query("type",testtype);
// 	amrex::Vector<amrex::Real> tmpn; pptest.queryarr("n",tmpn);
// 	amrex::Vector<amrex::Real> tmpb; pptest.queryarr("b",tmpb);
// 	amrex::Real tmpalpha; pptest.query("alpha",tmpalpha);
// 	amrex::Real tmpm = 1.0; pptest.query("m",tmpm);
// 	int comp; pptest.query("comp",comp);

// 	if (test_on == 0)
// 	{
// 		mlmg.solve(GetVecOfPtrs(u), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);
// 	}
// 	else
// 	{
// 		Set::Scalar alpha = tmpalpha, m = tmpm;
// 		//Set::Vector n(1.0,1.0); Set::Vector b(0.5,0.0);
// 		Set::Vector n(tmpn[0],tmpn[1]);
// 		Set::Vector b(tmpb[0],tmpb[1]);

// 		Util::Message(INFO,"alpha=",alpha);
// 		Util::Message(INFO,"n=",n.transpose());
// 		Util::Message(INFO,"b=",b.transpose());
// 		Util::Message(INFO,"comp=",comp);

		
// 		u[0].setVal(0.0);
// 		u[1].setVal(0.0);
// 		rhs[0].setVal(0.0);
// 		rhs[1].setVal(0.0);

// 		if (testtype=="affine")
// 		{
// 			Util::Message(INFO,"affine, comp = ", comp);
// 			IC::Affine ic(geom,n,alpha,b,true,m);
// 			if (comp == 0) {ic.SetComp(0); ic.Initialize(0,u); ic.Initialize(1,u);}
// 			if (comp == 1) {ic.SetComp(1); ic.Initialize(0,u); ic.Initialize(1,u);}
// 		}
// 		else if (testtype=="trig")
// 		{
// 			Util::Message(INFO,"trig");
// 			IC::Trig ic(geom);
// 			if (comp == 0) {ic.SetComp(0); ic.Initialize(0,u); ic.Initialize(1,u);}
// 			if (comp == 1) {ic.SetComp(1); ic.Initialize(0,u); ic.Initialize(1,u);}
// 		}
// 		else if (testtype=="random")
// 		{
// 			Util::Message(INFO,"random");
// 			IC::Random ic(geom);
// 			if (comp == 0) {ic.SetComp(0); ic.Initialize(0,u); ic.Initialize(1,u);}
// 			if (comp == 1) {ic.SetComp(1); ic.Initialize(0,u); ic.Initialize(1,u);}
// 		}
// 		else
// 			Util::Abort(INFO,"invalid test type");



// 		mlabec.FApply(0,0,rhs[0],u[0]);
// 		mlabec.FApply(1,0,rhs[1],u[1]);

// 		res[0].setVal(0.0);
// 		res[1].setVal(0.0);

// 		mlabec.BuildMasks();
// 		mlabec.Reflux(0,
// 			      res[0], u[0], rhs[0],
// 			      res[1], u[1], rhs[1]);
// 	}
	


// 	// 	mlabec.FSmooth(0,0,u[0],rhs[0]);
// 	// 	mlabec.FSmooth(1,0,u[1],rhs[1]);
// 	// }




// 	//

// 	//
// 	// Compute post-solve values
// 	//

// 	for (int lev = 0; lev < nlevels; lev++)
// 		{
// 			mlabec.Stress(lev,stress[lev],u[lev]);
// 			mlabec.Energy(lev,energy[lev],u[lev]);
// 		}
		
// 	//
// 	// WRITE PLOT FILE
// 	//

// #if AMREX_SPACEDIM==2
// 	const int ncomp = 10;
// 	Vector<std::string> varname = {"u01", "u02", "rhs01", "rhs02", "res01", "res02", "stress11", "stress22", "stress12", "energy"};
// #elif AMREX_SPACEDIM==3
// 	const int ncomp = 16;
// 	Vector<std::string> varname = {"u01", "u02", "u03", "rhs01", "rhs02", "rhs03", "res01", "res02", "res03",
// 				       "stress11", "stress22", "stress33", "stress23", "stress13", "stress12", "energy"};
// #endif

// 	nlevels = max_level+1;

// 	Vector<MultiFab> plotmf(nlevels);
// 	for (int ilev = 0; ilev < nlevels; ++ilev)
// 		{
// 			plotmf[ilev].define(ngrids[ilev], dmap[ilev], ncomp, 0);
// #if AMREX_SPACEDIM == 2
// 			MultiFab::Copy(plotmf[ilev], u      [ilev], 0, 0, 1, 0);
// 			MultiFab::Copy(plotmf[ilev], u      [ilev], 1, 1, 1, 0);
// 			MultiFab::Copy(plotmf[ilev], rhs    [ilev], 0, 2, 1, 0);
// 			MultiFab::Copy(plotmf[ilev], rhs    [ilev], 1, 3, 1, 0);
// 			MultiFab::Copy(plotmf[ilev], res    [ilev], 0, 4, 1, 0);
// 			MultiFab::Copy(plotmf[ilev], res    [ilev], 1, 5, 1, 0);
// 			MultiFab::Copy(plotmf[ilev], stress [ilev], 0, 6, 1, 0);
// 			MultiFab::Copy(plotmf[ilev], stress [ilev], 3, 7, 1, 0);
// 			MultiFab::Copy(plotmf[ilev], stress [ilev], 1, 8, 1, 0);
// 			MultiFab::Copy(plotmf[ilev], energy [ilev], 0, 9, 1, 0);
// #elif AMREX_SPACEDIM == 3
// 			MultiFab::Copy(plotmf[ilev], u      [ilev], 0, 0,  1, 0);
// 			MultiFab::Copy(plotmf[ilev], u      [ilev], 1, 1,  1, 0);
// 			MultiFab::Copy(plotmf[ilev], u      [ilev], 2, 2,  1, 0);
// 			MultiFab::Copy(plotmf[ilev], rhs    [ilev], 0, 3,  1, 0);
// 			MultiFab::Copy(plotmf[ilev], rhs    [ilev], 1, 4,  1, 0);
// 			MultiFab::Copy(plotmf[ilev], rhs    [ilev], 2, 5,  1, 0);
// 			MultiFab::Copy(plotmf[ilev], res    [ilev], 0, 6,  1, 0);
// 			MultiFab::Copy(plotmf[ilev], res    [ilev], 1, 7,  1, 0);
// 			MultiFab::Copy(plotmf[ilev], res    [ilev], 2, 8,  1, 0);
// 			MultiFab::Copy(plotmf[ilev], stress [ilev], 0, 9,  1, 0);
// 			MultiFab::Copy(plotmf[ilev], stress [ilev], 4, 10,  1, 0);
// 			MultiFab::Copy(plotmf[ilev], stress [ilev], 8, 11,  1, 0);
// 			MultiFab::Copy(plotmf[ilev], stress [ilev], 5, 12,  1, 0);
// 			MultiFab::Copy(plotmf[ilev], stress [ilev], 2, 13, 1, 0);
// 			MultiFab::Copy(plotmf[ilev], stress [ilev], 1, 14, 1, 0);
// 			MultiFab::Copy(plotmf[ilev], energy	[ilev], 0, 15, 1, 0);
// #endif 
// 		}

// 	IO::FileNameParse(plot_file);

// 	Util::Message(INFO,"varname size = ", varname.size());
// 	Util::Message(INFO,"mf->nComp() = ", plotmf[0].nComp());
	
// 	IO::WriteMetaData(plot_file);

// 	pptest.query("nlevels",nlevels);
// 	Util::Warning(INFO,"TODO: change nlevels back!");// nlevels=1;
// 	WriteMultiLevelPlotfile(plot_file, nlevels, amrex::GetVecOfConstPtrs(plotmf),
// 	 			varname, geom, 0.0, Vector<int>(nlevels, 0),
// 	 			Vector<IntVect>(nlevels, IntVect{ref_ratio}));
	


	
}


}

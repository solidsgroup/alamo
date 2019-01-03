#include <AMReX.H>
#include <AMReX_MLMG.H>
#include <AMReX_PlotFileUtil.H>
#include "Operator/Elastic.H"
#include "Model/Solid/LinearElastic/Isotropic.H"
#include "Model/Solid/LinearElastic/Degradable/Isotropic.H"
#include "Test.H"
#include "IC/Affine.H"
#include "IC/Trig.H"
#include "IC/Random.H"
#include "IO/WriteMetaData.H"
#include "IO/FileNameParse.H"


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
// Let b(x,y) = b0 * sin(pi x) sin(pi y)
// so the solution is
//     u_i^{mn}  = b0 / (alpha * pi^2 * 2) if m=n=1,    0 else
// 
//




namespace Operator
{

/// \todo Add verbosity messages 
template<>
int Test<Elastic<Model::Solid::LinearElastic::Degradable::Isotropic> >::SpatiallyVaryingCTest(int /*verbose*/)
{
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
	std::map<std::string,Elastic<model_type>::BC > 			bc_map;
	amrex::Vector<std::string> AMREX_D_DECL(bc_x_lo_str,bc_y_lo_str,bc_z_lo_str);
	amrex::Vector<std::string> AMREX_D_DECL(bc_x_hi_str,bc_y_hi_str,bc_z_hi_str);
	std::array<Elastic<model_type>::BC,AMREX_SPACEDIM> AMREX_D_DECL(bc_x_lo,bc_y_lo,bc_z_lo);
	std::array<Elastic<model_type>::BC,AMREX_SPACEDIM> AMREX_D_DECL(bc_x_hi,bc_y_hi,bc_z_hi);

	bc_map["displacement"] = 	Elastic<model_type>::BC::Displacement;
	bc_map["disp"] = 			Elastic<model_type>::BC::Displacement;
	bc_map["traction"] = 		Elastic<model_type>::BC::Traction;
	bc_map["trac"] = 			Elastic<model_type>::BC::Traction;
	bc_map["periodic"] = 		Elastic<model_type>::BC::Periodic;

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
					bc_front = {AMREX_D_DECL(0.,0.,0.)};
				);

	
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
	return 1;
}

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

	std::vector<int> comps = {{1}};
	std::vector<Set::Scalar> alphas = {{1.0}};
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

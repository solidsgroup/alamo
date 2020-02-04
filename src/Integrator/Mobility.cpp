#include "Mobility.H"
#include "BC/Constant.H"
#include "Operator/Elastic.H"
#include "Solver/Nonlocal/Linear.H"
//#include "Model/Solid/LinearElastic/Isotropic.H"
#include "Model/Solid/LinearElastic/MultiWell.H"
#include "BC/Operator/Elastic.H"
#include "IC/Sphere.H"
#include "IC/Affine.H"
#include "IC/Constant.H"
#include "IC/TabulatedInterface.H"
#include "Numeric/Stencil.H"

namespace Integrator
{
Mobility::Mobility() :
	Integrator()
{
	amrex::ParmParse pp("heat");
	pp.query("alpha", alpha);
	pp.query("refinement_threshold", refinement_threshold);
	//pp.query("ic_type", ic_type);

	{
		amrex::ParmParse pp("amr");
		pp.query("ref_criterion", amr.ref_criterion);
	}

	{
		amrex::ParmParse pp("ic");
		std::string type;
		pp.query("type", type);

		// // Determine initial condition
		if (type=="tabulated")
		{
			std::vector<Set::Scalar> xs, ys;
			pp.queryarr("xs",xs);
			pp.queryarr("ys",ys);
			Set::Scalar gammagb1 = -0.1, gammagb2 = 0.1;
			pp.query("gammagb1",gammagb1);
			pp.query("gammagb2",gammagb2);
			ic = new IC::TabulatedInterface(geom,xs,ys,gammagb1,gammagb2);
		}
		else
			ic = new IC::Constant(geom);
	}

	{
		amrex::ParmParse pp("bc");
		amrex::Vector<std::string> bc_hi_str(AMREX_SPACEDIM);
		amrex::Vector<std::string> bc_lo_str(AMREX_SPACEDIM);
		pp.queryarr("lo",bc_lo_str,0,BL_SPACEDIM);
		pp.queryarr("hi",bc_hi_str,0,BL_SPACEDIM);
		amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
		if (pp.countval("lo_1")) pp.getarr("lo_1",bc_lo_1);
		if (pp.countval("hi_1")) pp.getarr("hi_1",bc_hi_1);
		amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
		if (pp.countval("lo_2")) pp.getarr("lo_2",bc_lo_2);
		if (pp.countval("hi_2")) pp.getarr("hi_2",bc_hi_2);
		amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;
		if (pp.countval("lo_3")) pp.getarr("lo_3",bc_lo_3);
		if (pp.countval("hi_3")) pp.getarr("hi_3",bc_hi_3);

		mybc = new BC::Constant(1,bc_hi_str, bc_lo_str,
					AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3),
					AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3));
	}
	{
		amrex::ParmParse pp("bc.disp");
		pp.query("yhi_x",bc.disp.yhi_x);
	}
	{
		amrex::ParmParse pp("solver");
		pp.query("interval",solver.interval);
		pp.query("bottom_max_iter",solver.bottom_max_iter);
		pp.query("bottomsolver",solver.bottomsolver);
		pp.query("fixed_iter",solver.fixed_iter);
		pp.query("max_coarsening_level",solver.max_coarsening_level);
		pp.query("verbose",solver.verbose);
	}

	{
		amrex::ParmParse pp("physics");
		pp.query("gamma",physics.gamma);
		pp.query("elastic_mult",physics.elastic_mult);
		pp.query("kappa",physics.kappa);
		pp.query("L",physics.L);
		pp.query("gammagb0",physics.gammagb0);
	}

	RegisterNewFab(gammagb_mf,    mybc, 1, 3, "gammagb",true);
	RegisterNewFab(gammagbold_mf, mybc, 1, 3, "gammagbold",false);

	RegisterNodalFab(disp, AMREX_SPACEDIM, number_of_ghost_cells, "Disp",true);
	RegisterNodalFab(rhs,  AMREX_SPACEDIM, number_of_ghost_cells, "RHS",true);
	RegisterNodalFab(res,  AMREX_SPACEDIM, number_of_ghost_cells, "res",true);
	RegisterNodalFab(energy_mf, 1, number_of_ghost_cells, "energy",true);
	RegisterNodalFab(sigma,  AMREX_SPACEDIM*AMREX_SPACEDIM, number_of_ghost_cells, "sigma",true);
}

Mobility::~Mobility()
{
}

void
Mobility::Initialize (int lev)
{
	//ic->Initialize(lev,TempFab);
	disp[lev]->setVal(0.0);
	rhs[lev]->setVal(0.0);
	
	ic->Initialize(lev,gammagb_mf);
}

void 
Mobility::TimeStepBegin(amrex::Real /*time*/, int iter)
{
	Util::Message(INFO);
	if (!solver.interval || iter%solver.interval) return;
	Util::Message(INFO);

	for (int lev = 0; lev < disp.size(); ++lev)
	{
		//disp[lev]->setVal(0.0);
		rhs[lev]->setVal(0.0);
	}

	Set::Scalar lame = 2.6, shear = 6.0;
	using model_type = Model::Solid::LinearElastic::Multiwell;
	Operator::Elastic<model_type> elastic;
	elastic.SetUniform(false);

	amrex::LPInfo info;
	if (solver.max_coarsening_level > -1) info.setMaxCoarseningLevel(solver.max_coarsening_level);
	elastic.define(geom,grids,dmap,info);

	model_type mymodel(lame,shear,Set::Matrix::Zero());
	
	amrex::Vector<amrex::FabArray<amrex::BaseFab<model_type> > > model_mf;
	model_mf.resize(disp.size());


	for (int lev = 0; lev < disp.size(); ++lev)
	{
		amrex::Box domain(geom[lev].Domain());
		domain.convert(amrex::IntVect::TheNodeVector());
		model_mf[lev].define(disp[lev]->boxArray(), disp[lev]->DistributionMap(), 1, 2);
		model_mf[lev].setVal(mymodel);
		gammagb_mf[lev]->FillBoundary();
		
		const amrex::Real* DX  = geom[lev].CellSize();

		for (MFIter mfi(model_mf[lev],amrex::TilingIfNotGPU());mfi.isValid();++mfi)
		{
			//amrex::Box bx = mfi.growntilebox(2);
			amrex::Box bx = mfi.tilebox();
			bx.grow(2);
			//bx = bx & domain;

			amrex::Array4<model_type> const & model = model_mf[lev].array(mfi);
			amrex::Array4<const Set::Scalar> const & gammagb = gammagb_mf[lev]->array(mfi);
			amrex::Array4<Set::Scalar> const & RHS = rhs[lev]->array(mfi);

			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) 
			{
				Set::Matrix Fgb = Set::Matrix::Zero();
				Fgb(0,1) = 0.25*(gammagb(i,j,k) + gammagb(i-1,j,k) + gammagb(i,j-1,k)+ gammagb(i-1,j-1,k));
				RHS(i,j,k,0)=gammagb(i,j,k);
				model(i,j,k) = model_type(lame,shear,Fgb);
			});
			bx = mfi.tilebox();
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) 
			{
				model(i,j,k).gradFgb[0](1,0) = ((model(i+1,j,k).Fgb - model(i-1,j,k).Fgb)/2./DX[0])(0,1);
				model(i,j,k).gradFgb[0](1,1) = ((model(i,j+1,k).Fgb - model(i,j-1,k).Fgb)/2./DX[1])(0,1);
			});
		}
	}


	elastic.SetModel(model_mf);

	BC::Operator::Elastic<model_type> bc;
	for (int lev = 0; lev < rhs.size(); lev++) rhs[lev]->setVal(0.0);
	for (int lev = 0; lev < rhs.size(); lev++) disp[lev]->setVal(0.0);
    bc.Set(bc.Face::XLO, bc.Direction::X, bc.Type::Neumann, 0.0,      rhs, geom);
    bc.Set(bc.Face::XLO, bc.Direction::Y, bc.Type::Neumann, 0.0,      rhs, geom);
    bc.Set(bc.Face::XHI, bc.Direction::X, bc.Type::Neumann, 0.0,      rhs, geom);
    bc.Set(bc.Face::XHI, bc.Direction::Y, bc.Type::Neumann, 0.0,      rhs, geom);
    bc.Set(bc.Face::YLO, bc.Direction::X, bc.Type::Displacement, 0.0, rhs, geom);
    bc.Set(bc.Face::YLO, bc.Direction::Y, bc.Type::Displacement, 0.0, rhs, geom);
    bc.Set(bc.Face::YHI, bc.Direction::X, bc.Type::Displacement, this->bc.disp.yhi_x, rhs, geom);
    bc.Set(bc.Face::YHI, bc.Direction::Y, bc.Type::Displacement, 0.0, rhs, geom);

	elastic.SetBC(&bc);
	Solver::Nonlocal::Linear linearsolver(elastic);
	if (solver.verbose)                    linearsolver.setVerbose(solver.verbose);
	if (solver.bottom_max_iter > -1)       linearsolver.setBottomMaxIter(solver.bottom_max_iter);
	if (solver.bottomsolver == "smoother") linearsolver.setBottomSolver(amrex::MLMG::BottomSolver::smoother);
	if (solver.fixed_iter > -1)            linearsolver.setFixedIter(solver.fixed_iter);
	
	



	Set::Scalar tol_rel = 1E-8, tol_abs = 1E-8;
	linearsolver.solveaffine(disp,rhs,tol_rel,tol_abs,true);
	
	//elastic.SetHomogeneous(false);
	linearsolver.compResidual(GetVecOfPtrs(res),GetVecOfPtrs(disp),GetVecOfConstPtrs(rhs));

	//for (int lev = 0; lev < disp.size(); ++lev)
	//{
	//	amrex::Box domain(geom[lev].Domain());
	//	domain.convert(amrex::IntVect::TheNodeVector());
	//	for (MFIter mfi(model_mf[lev],amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	//	{
	//		amrex::Box bx = mfi.tilebox();
	//		bx.grow(model_mf[lev].nGrow());
	//		bx = bx & domain;
	//		//amrex::Array4<model_type> const & model = model_mf[lev].array(mfi);
	//		//amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) { model(i,j,k).SetHomogeneous(false);} );
	//	}
	//}
	elastic.SetModel(model_mf);
	
	for (int lev = 0; lev < sigma.size(); lev++)
	{
		elastic.Stress(lev,*sigma[lev],*disp[lev]);
		elastic.Energy(lev,*energy_mf[lev],*disp[lev]);
	}
	
}

void
Mobility::Advance (int lev, amrex::Real /*time*/, amrex::Real dt)
{
	std::swap(gammagbold_mf[lev], gammagb_mf[lev]);
	const amrex::Real* DX = geom[lev].CellSize();

	for (amrex::MFIter mfi(*gammagb_mf[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		amrex::Box bx = mfi.tilebox();
		bx.grow(2);

		amrex::Array4<Set::Scalar> const & gammagb = gammagb_mf[lev]->array(mfi);
		amrex::Array4<const Set::Scalar> const & gammagbold = gammagbold_mf[lev]->array(mfi);
		amrex::Array4<const Set::Scalar> const & Sigma = sigma[lev]->array(mfi);
		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) 
		{
			Set::Scalar df = 0.0;
			
			
			df += - physics.gamma *(Set::Constant::Pi/physics.gammagb0) * sin(Set::Constant::Pi * gammagbold(i,j,k) / physics.gammagb0);
			Set::Scalar sig12 = 0.25*(Sigma(i,j,k,1) + Sigma(i+1,j,k,1) + Sigma(i,j+1,k,1)+ Sigma(i+1,j+1,k,1));
			df += - physics.elastic_mult*sig12;
			df += - physics.kappa*Numeric::Laplacian(gammagbold,i,j,k,0,DX);
			gammagb(i,j,k) = gammagbold(i,j,k) - dt * physics.L * df;
		});
	}
}

void
Mobility::TagCellsForRefinement (int lev, amrex::TagBoxArray& a_tags, amrex::Real /*time*/, int /*ngrow*/)
{
	const amrex::Real* DX      = geom[lev].CellSize();
	const Set::Scalar dxnorm = sqrt(DX[0]*DX[0] +  DX[1]*DX[1]);

	for (amrex::MFIter mfi(*gammagb_mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		amrex::Box bx = mfi.tilebox();
		//bx.grow(-1);
		//amrex::Array4<const amrex::Real> const& RHS    = (*rhs[lev]).array(mfi);
		//amrex::Array4<const amrex::Real> const& u    = (*disp[lev]).array(mfi);
		amrex::Array4<char> const& tags    = a_tags.array(mfi);
		amrex::Array4<Set::Scalar> const & gammagb = gammagb_mf[lev]->array(mfi);

		for (int n = 0; n < 2 ; n++)
		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k)
		{
			Set::Vector grad = Numeric::Gradient(gammagb,i,j,k,0,DX);
			//tags(i,j,k) = amrex::TagBox::SET;
			if (grad.lpNorm<2>()*dxnorm > amr.ref_criterion) tags(i,j,k) = amrex::TagBox::SET;
		});
	}
	
}

}

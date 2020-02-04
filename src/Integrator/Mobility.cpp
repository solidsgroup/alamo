#include "Mobility.H"
#include "BC/Constant.H"
#include "Operator/Elastic.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"
#include "Model/Solid/Affine/Cubic.H"
#include "Model/Solid/LinearElastic/Isotropic.H"
#include "BC/Operator/Elastic.H"
//#include "BC/Step.H"
#include "IC/Sphere.H"
#include "IC/Affine.H"
#include "IC/Constant.H"
#include "IC/TabulatedInterface.H"
#include "IC/DoubleNotch.H"
#include "Numeric/Stencil.H"

#include "Util/Util.H"
//#include "BC/Step.H"

#include "IO/ParmParse.H"

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
		amrex::ParmParse pp("physics");
		pp.query("gamma",physics.gamma);
		pp.query("elastic_mult",physics.elastic_mult);
		pp.query("kappa",physics.kappa);
		pp.query("L",physics.L);
		pp.query("gammagb1",physics.gammagb1);
		pp.query("gammagb2",physics.gammagb2);
	}

	{
		amrex::ParmParse pp("ic");
		std::string type;
		pp.query("type", type);

		// // Determine initial condition
		if (type=="tabulated")
		{
			pp.queryarr("xs",xs);
			pp.queryarr("ys",ys);
			ic = new IC::TabulatedInterface(geom,xs,ys,physics.gammagb1,physics.gammagb2,IC::TabulatedInterface::Type::Values);
		}
		else if (type == "circle")
		{
			amrex::Vector<Set::Scalar> center;
			pp.queryarr("center",center);

			ic = new IC::Sphere(geom, 0.1, center, IC::Sphere::Type::XYZ ,physics.gammagb1, physics.gammagb2);

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
//		mybc = new BC::Step(bc_hi_str, bc_lo_str,
//					AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3),
//					AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3),ys[0],ys[ys.size()-1],gammagb1,gammagb2);

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
		IO::ParmParse pp("elastic");
		pp.queryclass("model",elastic.model);
		pp.queryclass("bc",elastic.bc);
	}

	RegisterNewFab(gammagb_mf,    mybc, 1, 3, "gammagb",true);
	RegisterNewFab(gammagbold_mf, mybc, 1, 3, "gammagbold",false);
	RegisterNewFab(L_mf, mybc, 1, number_of_ghost_cells, "L",true);

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

	IO::ParmParse pp("io");
	IC::DoubleNotch icdn(geom);
	pp.queryclass("doublenotch",icdn);
	icdn.Initialize(lev,L_mf);
}

void 
Mobility::TimeStepBegin(amrex::Real /*time*/, int iter)
{
	if (!solver.interval || iter%solver.interval) return;

	//elastic.op.SetUniform(false);

	amrex::LPInfo info;
	//elastic.op = Operator::Elastic<model_type>(geom,grids,dmap,info);
	Operator::Elastic<model_type> op;
	op.define(geom,grids,dmap,info);
	op.SetUniform(false);

	Set::Field<model_type> model_mf;
	model_mf.resize(disp.size());

	for (int lev = 0; lev < disp.size(); ++lev)
	{
		rhs[lev]->setVal(0.0);
		disp[lev]->setVal(0.0);
		amrex::Box domain(geom[lev].Domain());

		domain.convert(amrex::IntVect::TheNodeVector());
		model_mf.Define(lev,disp[lev]->boxArray(), disp[lev]->DistributionMap(), 1, 2);
		model_mf[lev]->setVal(elastic.model);
		gammagb_mf[lev]->FillBoundary();
		
		const amrex::Real* DX  = geom[lev].CellSize();

		for (MFIter mfi(*model_mf[lev],amrex::TilingIfNotGPU());mfi.isValid();++mfi)
		{
			amrex::Box bx = mfi.tilebox();
			bx.grow(2);

			amrex::Array4<model_type> const & model = model_mf[lev]->array(mfi);
			amrex::Array4<const Set::Scalar> const & gammagb = gammagb_mf[lev]->array(mfi);

			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) 
			{
				Set::Matrix Fgb = Set::Matrix::Zero();
				Fgb(0,1) = 0.25*(gammagb(i,j,k) + gammagb(i-1,j,k) + gammagb(i,j-1,k)+ gammagb(i-1,j-1,k));
				model(i,j,k).F0 = Fgb;
			});

			//Util::RealFillBoundary(model_mf[lev],geom[lev]); // this is causing the problem..
		}
	}
	op.SetModel(model_mf);


	elastic.bc.Init(rhs,geom);
	op.SetBC(&elastic.bc);
	
	elastic.newton = new Solver::Nonlocal::Newton<model_type>(op);
	IO::ParmParse pp("elastic");
	pp.queryclass("newton",*elastic.newton);
	
	Set::Scalar tol_rel = 1E-8, tol_abs = 1E-8;
	elastic.newton->solve(disp,rhs,model_mf,tol_rel,tol_abs);
	elastic.newton->compResidual(res,disp,rhs,model_mf);
	
	for (int lev = 0; lev < sigma.size(); lev++)
	{
		op.Stress(lev,*sigma[lev],*disp[lev]);
		op.Energy(lev,*energy_mf[lev],*disp[lev]);
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
		amrex::Array4<const Set::Scalar> const & L = L_mf[lev]->array(mfi);
		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) 
		{
			Set::Scalar x1 = geom[lev].ProbLo()[0] + ((Set::Scalar)(i)) * geom[lev].CellSize()[0];

			//if (x1 >= xs[1] && x1 <= xs[xs.size()-2]) 
			{
				Set::Scalar df = 0.0;

				Set::Scalar gammagb0 = 0.5*fabs(physics.gammagb2 - physics.gammagb1);

				df += - physics.gamma *(Set::Constant::Pi/gammagb0) * sin(Set::Constant::Pi * gammagbold(i,j,k) / gammagb0);
				Set::Scalar sig12 = 0.25*(Sigma(i,j,k,1) + Sigma(i+1,j,k,1) + Sigma(i,j+1,k,1)+ Sigma(i+1,j+1,k,1));
				df += - physics.elastic_mult*sig12;
				df += - physics.kappa*Numeric::Laplacian(gammagbold,i,j,k,0,DX);
				gammagb(i,j,k) = gammagbold(i,j,k) - dt * L(i,j,k) * df;
				//gammagb(i,j,k) = 0.0;
			}
			//else
			//{
			//	gammagb(i,j,k) = gammagbold(i,j,k);
			//}
		});
	}
}

void
Mobility::TagCellsForRefinement (int lev, amrex::TagBoxArray& a_tags, amrex::Real /*time*/, int /*ngrow*/)
{
	const amrex::Real* DX      = geom[lev].CellSize();
	const Set::Scalar dxnorm = std::sqrt(DX[0]*DX[0] +  DX[1]*DX[1]);

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

#include "HeatConduction.H"
#include "BC/Constant.H"
#include "Operator/Elastic.H"
#include "Solver/Nonlocal/Linear.H"
#include "Model/Solid/LinearElastic/Isotropic.H"
#include "Model/Solid/LinearElastic/MultiWell.H"
#include "BC/Operator/Elastic.H"
#include "IC/Sphere.H"
#include "IC/Affine.H"
//#include "IC/PS.H"
#include "Numeric/Stencil.H"

namespace Integrator
{
/// \fn    HeatConduction::Integrator::Integrator
///
/// Read in the following simulation parameters
///
///     heat.alpha                (default 1.0)
///     heat.refinement_threshold (default 0.01)
///     ic.type
///

/// Initialize initial condition pointer #ic, and register
/// the #Temp, #Temp_old Multifab arrays.

HeatConduction::HeatConduction() :
	Integrator()
{
	amrex::ParmParse pp("heat");
	pp.query("alpha", alpha);
	pp.query("refinement_threshold", refinement_threshold);
	pp.query("ic_type", ic_type);

	// // Determine initial condition
	if (ic_type == "cylinder")
		ic = new IC::Cylinder(geom);
	else if (ic_type == "sphere")
		ic = new IC::Sphere(geom);
	else if (ic_type == "constant")
		ic = new IC::Constant(geom);
//	else if (ic_type == "packed_spheres")
//		ic = new IC::PS(geom,20,0.0,1.0);
	else
		ic = new IC::Constant(geom);
    
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

		mybc = new BC::Constant(bc_hi_str, bc_lo_str,
					AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3),
					AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3));
	}

	RegisterNewFab(gammagb_mf,    mybc, 1, 2, "gammagb");
	RegisterNewFab(gammagbold_mf, mybc, 1, 2, "gammagbold");

	RegisterNodalFab(disp, AMREX_SPACEDIM, number_of_ghost_cells, "Disp");
	RegisterNodalFab(rhs,  AMREX_SPACEDIM, number_of_ghost_cells, "RHS");
	RegisterNodalFab(res,  AMREX_SPACEDIM, number_of_ghost_cells, "res");
	RegisterNodalFab(energy_mf, 1, number_of_ghost_cells, "energy");
	RegisterNodalFab(sigma,  AMREX_SPACEDIM*AMREX_SPACEDIM, number_of_ghost_cells, "sigma");
}

HeatConduction::~HeatConduction()
{
}

/// \fn HeatConduction::Integrator::Initialize
///
/// Use the #ic object to initialize #Temp
void
HeatConduction::Initialize (int lev)
{
	//ic->Initialize(lev,TempFab);
	disp[lev]->setVal(0.0);
	rhs[lev]->setVal(0.0);

	for (amrex::MFIter mfi(*gammagb_mf[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		amrex::Box bx = mfi.growntilebox(2);

		amrex::Array4<Set::Scalar> const & gammagb = gammagb_mf[lev]->array(mfi);
		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) 
		{
			Set::Scalar AMREX_D_DECL(x,y,z);
			AMREX_D_TERM(x = geom[lev].ProbLo()[0] + (amrex::Real)(i) * geom[lev].CellSize()[0];,
						 y = geom[lev].ProbLo()[1] + (amrex::Real)(j) * geom[lev].CellSize()[1];,
						 z = geom[lev].ProbLo()[2] + (amrex::Real)(k) * geom[lev].CellSize()[2];);

			// Simple bicrystal
			if (y >= 0.5)  gammagb(i,j,k) = 0.1;
			else gammagb(i,j,k)=-0.1;

			//// Single Disconnection
			//if (x < 0.5 && y > 0.6)  gammagb(i,j,k) = 0.1;
			//else if (x >= 0.5 && y > 0.4) gammagb(i,j,k) = 0.1;
			//else gammagb(i,j,k)=-0.1;

			//// Sinusoidal perturbation
			//if (y > 0.7 + 0.01*sin(8. * Set::Constant::Pi * x))  gammagb(i,j,k) = 0.1;
			//else gammagb(i,j,k)=-0.1;

			// // Ramp
			//if (y > (0.55 - 0.1 * x)) gammagb(i,j,k) = 0.1;
			//else gammagb(i,j,k) = -0.1;

			// // Single Disconnection
			//if (x < 0.25 && y > 0.525)  gammagb(i,j,k) = 0.1;
			//else if (0.25 <= x && x < 0.75 && y > 0.475) gammagb(i,j,k) = 0.1;
			//else if (0.75 <= x && y > 0.525)  gammagb(i,j,k) = 0.1;
			//else gammagb(i,j,k)=0;

		});
	}
//	//ic->Initialize(lev,TempFab);
//	disp[lev]->setVal(0.0);
//	rhs[lev]->setVal(0.0);
//
//	for (amrex::MFIter mfi(*gammagb_mf[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
//	{
//		amrex::Box bx = mfi.growntilebox(2);
//
//		amrex::Array4<Set::Scalar> const & gammagb = gammagb_mf[lev]->array(mfi);
//		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) 
//		{
//			Set::Scalar AMREX_D_DECL(x,y,z);
//			AMREX_D_TERM(x = geom[lev].ProbLo()[0] + (amrex::Real)(i) * geom[lev].CellSize()[0];,
//						 y = geom[lev].ProbLo()[1] + (amrex::Real)(j) * geom[lev].CellSize()[1];,
//						 z = geom[lev].ProbLo()[2] + (amrex::Real)(k) * geom[lev].CellSize()[2];);
//
//			Set::Matrix Fgb = Set::Matrix::Zero();
//			if (x < 0.5 && y > 0.525)  gammagb(i,j,k) = 0.1;
//			else if (x >= 0.5 && y > 0.475) gammagb(i,j,k) = 0.1;
//			else gammagb(i,j,k)=-0.1;
//		});
//	}

	

}

void 
HeatConduction::TimeStepBegin(amrex::Real /*time*/, int iter)
{
	//for (int lev = 0; lev < disp.size(); ++lev)
	//{
	//	disp[lev]->setVal(0.0);
	//	rhs[lev]->setVal(0.0);
	//	energy_mf[lev]->setVal(0.0);
	//	sigma[lev]->setVal(0.0);
	//}
	//return;

	//for (int lev = 0; lev < disp.size(); ++lev)	{disp[lev]->setVal(0.0);rhs[lev]->setVal(0.0);energy_mf[lev]->setVal(0.0);sigma[lev]->setVal(0.0);} Util::Warning(INFO,"Note, settting everything to zero!"); return; 
	if (iter == 0) return;
	if (iter%plot_int) return;

Util::Message(INFO);

	for (int lev = 0; lev < disp.size(); ++lev)
	{
		disp[lev]->setVal(0.0);
		rhs[lev]->setVal(0.0);
	}

Util::Message(INFO);
	Set::Scalar lame = 2.6, shear = 6.0;
	using model_type = Model::Solid::LinearElastic::Multiwell;
	Util::Message(INFO);
	Operator::Elastic<model_type> elastic;
	elastic.SetHomogeneous(false);

	amrex::LPInfo info;
	info.setMaxCoarseningLevel(0);
	elastic.define(geom,grids,dmap,info);
	

Util::Message(INFO);
	model_type mymodel(lame,shear,Set::Matrix::Zero());
	
Util::Message(INFO);
	amrex::Vector<amrex::FabArray<amrex::BaseFab<model_type> > > model_mf;
	model_mf.resize(disp.size());

Util::Message(INFO);
	for (int lev = 0; lev < disp.size(); ++lev)
	{
		amrex::Box domain(geom[lev].Domain());
		domain.convert(amrex::IntVect::TheNodeVector());
		model_mf[lev].define(disp[lev]->boxArray(), disp[lev]->DistributionMap(), 1, 2);
		model_mf[lev].setVal(mymodel);
		
		const amrex::Real* DX  = geom[lev].CellSize();

		for (MFIter mfi(model_mf[lev],amrex::TilingIfNotGPU());mfi.isValid();++mfi)
		{
			//amrex::Box bx = mfi.growntilebox(2);
			amrex::Box bx = mfi.tilebox();
			bx.grow(1);
			//bx = bx & domain;

			amrex::Array4<model_type> const & model = model_mf[lev].array(mfi);
			amrex::Array4<const Set::Scalar> const & gammagb = gammagb_mf[lev]->array(mfi);
			amrex::Array4<Set::Scalar> const & RHS = rhs[lev]->array(mfi);

			amrex::Dim3 lo= amrex::lbound(bx), hi = amrex::ubound(bx);
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) 
			{
				Set::Matrix Fgb = Set::Matrix::Zero();
				//if      (i == lo.x && j == lo.y) Fgb(0,1) = gammagb(i,j,k);
				//else if (i == lo.x && j == hi.y) Fgb(0,1) = gammagb(i,j-1,k);
				//else if (i == hi.x && j == lo.y) Fgb(0,1) = gammagb(i-1,j,k);
				//else if (i == hi.x && j == hi.y) Fgb(0,1) = gammagb(i-1,j-1,k);
				//else if (i == lo.x) Fgb(0,1) = 0.5*(gammagb(i,j,k) + gammagb(i,j-1,k));
				//else if (i == hi.x) Fgb(0,1) = 0.5*(gammagb(i-1,j,k) + gammagb(i-1,j-1,k));
				//else if (i == lo.y) Fgb(0,1) = 0.5*(gammagb(i,j,k) + gammagb(i-1,j,k));
				//else if (i == hi.y) Fgb(0,1) = 0.5*(gammagb(i,j-1,k) + gammagb(i-1,j-1,k));
				//else 
				
				Fgb(0,1) = 0.25*(gammagb(i,j,k) + gammagb(i-1,j,k) + gammagb(i,j-1,k)+ gammagb(i-1,j-1,k));

				//Util::Message(INFO,Fgb(0,1));
				//RHS(i,j,k,0)=Fgb(0,1);
				RHS(i,j,k,0)=gammagb(i,j,k);
				
				model(i,j,k) = model_type(lame,shear,Fgb);

			});
			bx = mfi.tilebox();
			Util::Message(INFO);
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) 
			{
				model(i,j,k).gradFgb[0](1,0) = ((model(i+1,j,k).Fgb - model(i-1,j,k).Fgb)/2./DX[0])(0,1);
				model(i,j,k).gradFgb[0](1,1) = ((model(i,j+1,k).Fgb - model(i,j-1,k).Fgb)/2./DX[1])(0,1);
			});
			Util::Message(INFO);
		}
	}

	elastic.SetModel(model_mf);

	BC::Operator::Elastic<model_type> bc;
	for (int lev = 0; lev < rhs.size(); lev++) rhs[lev]->setVal(0.0);
    bc.Set(bc.Face::XLO, bc.Direction::X, bc.Type::Neumann, 0.0,      rhs, geom);
    bc.Set(bc.Face::XLO, bc.Direction::Y, bc.Type::Neumann, 0.0,      rhs, geom);
    bc.Set(bc.Face::XHI, bc.Direction::X, bc.Type::Neumann, 0.0,      rhs, geom);
    bc.Set(bc.Face::XHI, bc.Direction::Y, bc.Type::Neumann, 0.0,      rhs, geom);
    bc.Set(bc.Face::YLO, bc.Direction::X, bc.Type::Displacement, 0.0, rhs, geom);
    bc.Set(bc.Face::YLO, bc.Direction::Y, bc.Type::Displacement, 0.0, rhs, geom);
    bc.Set(bc.Face::YHI, bc.Direction::X, bc.Type::Displacement, 0.0, rhs, geom);
    bc.Set(bc.Face::YHI, bc.Direction::Y, bc.Type::Displacement, 0.0, rhs, geom);

	elastic.SetBC(&bc);
	Solver::Nonlocal::Linear solver(elastic);
	solver.setVerbose(4);
	solver.setBottomMaxIter(20);
	//solver.setBottomSolver(amrex::MLMG::BottomSolver::smoother);
	solver.setFixedIter(10000);

	solver.apply(GetVecOfPtrs(rhs),GetVecOfPtrs(disp));

	bc.Set(bc.Face::XLO, bc.Direction::X, bc.Type::Neumann, 0.0,      rhs, geom);
    bc.Set(bc.Face::XLO, bc.Direction::Y, bc.Type::Neumann, 0.0,      rhs, geom);
    bc.Set(bc.Face::XHI, bc.Direction::X, bc.Type::Neumann, 0.0,      rhs, geom);
    bc.Set(bc.Face::XHI, bc.Direction::Y, bc.Type::Neumann, 0.0,      rhs, geom);
    bc.Set(bc.Face::YLO, bc.Direction::X, bc.Type::Displacement, 0.0, rhs, geom);
    bc.Set(bc.Face::YLO, bc.Direction::Y, bc.Type::Displacement, 0.0, rhs, geom);
    bc.Set(bc.Face::YHI, bc.Direction::X, bc.Type::Displacement, 0.0, rhs, geom);
    bc.Set(bc.Face::YHI, bc.Direction::Y, bc.Type::Displacement, 0.0, rhs, geom);


	for (int lev = 0; lev < disp.size(); ++lev)
	{
		amrex::Box domain(geom[lev].Domain());
		domain.convert(amrex::IntVect::TheNodeVector());

		for (MFIter mfi(model_mf[lev],amrex::TilingIfNotGPU());mfi.isValid();++mfi)
		{
			amrex::Box bx = mfi.tilebox();
			bx.grow(model_mf[lev].nGrow());
			bx = bx & domain;

			amrex::Array4<model_type> const & model = model_mf[lev].array(mfi);
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) { model(i,j,k).SetHomogeneous(true);} );
		}
	}
	elastic.SetModel(model_mf);
	Set::Scalar tol_rel = 1E-8, tol_abs = 1E-8;
	solver.solve(GetVecOfPtrs(disp),GetVecOfConstPtrs(rhs),tol_rel,tol_abs);
	
	solver.compResidual(GetVecOfPtrs(res),GetVecOfPtrs(disp),GetVecOfConstPtrs(rhs));

	for (int lev = 0; lev < disp.size(); ++lev)
	{
		amrex::Box domain(geom[lev].Domain());
		domain.convert(amrex::IntVect::TheNodeVector());
		for (MFIter mfi(model_mf[lev],amrex::TilingIfNotGPU());mfi.isValid();++mfi)
		{
			amrex::Box bx = mfi.tilebox();
			bx.grow(model_mf[lev].nGrow());
			bx = bx & domain;
			amrex::Array4<model_type> const & model = model_mf[lev].array(mfi);
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) { model(i,j,k).SetHomogeneous(false);} );
		}
	}
	elastic.SetModel(model_mf);
	
	for (int lev = 0; lev < sigma.size(); lev++)
	{
		elastic.Stress(lev,*sigma[lev],*disp[lev]);
		elastic.Energy(lev,*energy_mf[lev],*disp[lev]);
	}
	
}


/// \fn    Integrator::HeatConduction::Advance
///
/// Integrate the heat diffusion equation
/// \f[\nabla^2T = \alpha \frac{\partial T}{\partial t}\f]
/// using an explicit forward Euler method.
/// \f$\alpha\f$ is stored in #alpha
void
HeatConduction::Advance (int lev, amrex::Real /*time*/, amrex::Real dt)
{
	std::swap(gammagbold_mf[lev], gammagb_mf[lev]);
	Set::Scalar L = 0.1;
	const amrex::Real* DX = geom[lev].CellSize();

	for (amrex::MFIter mfi(*gammagb_mf[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		amrex::Box bx = mfi.tilebox();
		bx.grow(1);

		amrex::Array4<Set::Scalar> const & gammagb = gammagb_mf[lev]->array(mfi);
		amrex::Array4<const Set::Scalar> const & gammagbold = gammagbold_mf[lev]->array(mfi);
		amrex::Array4<const Set::Scalar> const & Sigma = sigma[lev]->array(mfi);
		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) 
		{
			Set::Scalar df = 0.0;
			
			df += - 10.*(Set::Constant::Pi/0.1) * sin(Set::Constant::Pi * gammagbold(i,j,k) / 0.1);
			//df += 2.0 * (Set::Constant::Pi/0.1) * cos(2*Set::Constant::Pi * gammagbold(i,j,k) / 0.1);

			//Set::Scalar sig12 = 0.25*(Sigma(i,j,k,1) + Sigma(i+1,j,k,1) + Sigma(i,j+1,k,1)+ Sigma(i+1,j+1,k,1));
			//df += - 200.0*sig12;
			//df += - 300.0*Sigma(i,j,k,1);
			//df += - 100.0*Sigma(i,j,k,1);

			df += - 1.0*Numeric::Laplacian(gammagbold,i,j,k,0,DX);
			//df += -1*Numeric::Laplacian(gammagbold,i,j,k,0,DX);
			
			gammagb(i,j,k) = gammagbold(i,j,k) - dt * L * df;
			
			//Set::Matrix Fgb = Set::Matrix::Zero();
			//if (x < 0.5 && y > 0.525)  gammagb(i,j,k) = 0.1;
			//if (x >= 0.5 && y > 0.475) gammagb(i,j,k) = 0.1;
		});
	}
}

void
HeatConduction::TagCellsForRefinement (int lev, amrex::TagBoxArray& a_tags, amrex::Real time, int /*ngrow*/)
{
	const amrex::Real* DX      = geom[lev].CellSize();
	const Set::Scalar dxnorm = sqrt(DX[0]*DX[0] +  DX[1]*DX[1]);

	for (amrex::MFIter mfi(*gammagb_mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		amrex::Box bx = mfi.tilebox();
		bx.grow(-1);
		//amrex::Array4<const amrex::Real> const& RHS    = (*rhs[lev]).array(mfi);
		//amrex::Array4<const amrex::Real> const& u    = (*disp[lev]).array(mfi);
		amrex::Array4<char> const& tags    = a_tags.array(mfi);
		amrex::Array4<Set::Scalar> const & gammagb = gammagb_mf[lev]->array(mfi);

		for (int n = 0; n < 2 ; n++)
		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k)
		{
			Set::Vector grad = Numeric::Gradient(gammagb,i,j,k,0,DX);
			//tags(i,j,k) = amrex::TagBox::SET;
			if (grad.lpNorm<2>()*dxnorm > 0.01) tags(i,j,k) = amrex::TagBox::SET;
		});
	}
	
}

}

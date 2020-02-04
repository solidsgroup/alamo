
#include <eigen3/Eigen/Eigenvalues>

#include <omp.h>
#include <cmath>

#include <AMReX_SPACE.H>

#include "PhaseFieldMicrostructure.H"
#include "BC/Constant.H"
#include "Set/Set.H"
#include "Util/Util.H"
#include "IC/Random.H"
#include "IC/Trig.H"
#include "IC/Sphere.H"
//#include "Model/Solid/LinearElastic/Isotropic.H"
#include "Model/Solid/LinearElastic/MultiWell.H"
//#include "Model/Solid/LinearElastic/Laplacian.H"
#include "Model/Interface/GB/SH.H"
#include "Numeric/Stencil.H"
#include "Solver/Nonlocal/Linear.H"
#include "IC/Trig.H"
namespace Integrator
{
PhaseFieldMicrostructure::PhaseFieldMicrostructure() : Integrator()
{
	BL_PROFILE("PhaseFieldMicrostructure::PhaseFieldMicrostructure()");
	//
	// READ INPUT PARAMETERS
	//
	{
		amrex::ParmParse pp("pf"); // Phase-field model parameters
		pp.query("number_of_grains", number_of_grains);
		pp.query("M", pf.M);
		if (pp.contains("mu"))
			pp.query("mu", pf.mu);
		pp.query("gamma", pf.gamma);
		pp.query("sigma0", pf.sigma0);
		pp.query("l_gb", pf.l_gb);
		pp.query("elastic_mult",pf.elastic_mult);
	}
	{
		amrex::ParmParse pp("amr");
		pp.query("max_level", max_level);
		pp.query("ref_threshold", ref_threshold);
	}
	{
		amrex::ParmParse pp("lagrange");
		pp.query("on", lagrange.on);
		if (lagrange.on)
		{
			pp.query("lambda", lagrange.lambda);
			pp.query("vol0", lagrange.vol0);
			pp.query("tstart", lagrange.tstart);
			SetThermoInt(1);
		}
	}
	{
		amrex::ParmParse pp("anisotropy"); // Phase-field model parameters
		pp.query("on", anisotropy.on);
		pp.query("theta0", anisotropy.theta0);
		pp.query("phi0", anisotropy.phi0);
		pp.query("filename", filename);
		pp.query("gb_type", gb_type);
		anisotropy.theta0 *= 0.01745329251; // convert degrees into radians
		anisotropy.phi0 *= 0.01745329251;   // convert degrees into radians
		pp.query("sigma0", anisotropy.sigma0);
		pp.query("sigma1", anisotropy.sigma1);
		pp.query("beta", anisotropy.beta);
		pp.query("tstart", anisotropy.tstart);
		anisotropy.timestep = timestep;
		pp.query("timestep", anisotropy.timestep);
		anisotropy.plot_int = plot_int;
		pp.query("plot_int", anisotropy.plot_int);
		anisotropy.plot_dt = plot_dt;
		pp.query("plot_int", anisotropy.plot_dt);
		pp.query("thermo_int", anisotropy.thermo_int);
		pp.query("thermo_plot_int", anisotropy.thermo_plot_int);

		std::map<std::string, RegularizationType> regularization_type;
		regularization_type["wilmore"] = RegularizationType::Wilmore;
		regularization_type["k12"] = RegularizationType::K12;
		std::string regularization_type_input = "k12";
		pp.query("regularization", regularization_type_input);
		regularization = regularization_type[regularization_type_input];

		if (gb_type == "abssin")
			boundary = new Model::Interface::GB::AbsSin(anisotropy.theta0,
														anisotropy.sigma0,
														anisotropy.sigma1);
		else if (gb_type == "sin")
			boundary = new Model::Interface::GB::Sin(anisotropy.theta0,
													 anisotropy.sigma0,
													 anisotropy.sigma1);
		else if (gb_type == "read")
			boundary = new Model::Interface::GB::Read(filename);

		else if (gb_type == "sh")
		{
			//Need to make this check for other gb_types as well.
			if (AMREX_SPACEDIM < 2)
				Util::Abort(INFO, "SH model is only for 3D");
			boundary = new Model::Interface::GB::SH(anisotropy.theta0,
													anisotropy.phi0,
													anisotropy.sigma0,
													anisotropy.sigma1);
		}
		else
			boundary = new Model::Interface::GB::Sin(anisotropy.theta0, anisotropy.sigma0, anisotropy.sigma1);
	}

	{
		amrex::ParmParse pp("bc");
		amrex::Vector<std::string> bc_hi_str(AMREX_SPACEDIM);
		amrex::Vector<std::string> bc_lo_str(AMREX_SPACEDIM);
		pp.queryarr("lo", bc_lo_str, 0, BL_SPACEDIM);
		pp.queryarr("hi", bc_hi_str, 0, BL_SPACEDIM);
		amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
		if (pp.countval("lo_1"))
			pp.getarr("lo_1", bc_lo_1);
		if (pp.countval("hi_1"))
			pp.getarr("hi_1", bc_hi_1);
		amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
		if (pp.countval("lo_2"))
			pp.getarr("lo_2", bc_lo_2);
		if (pp.countval("hi_2"))
			pp.getarr("hi_2", bc_hi_2);
		amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;
		if (pp.countval("lo_3"))
			pp.getarr("lo_3", bc_lo_3);
		if (pp.countval("hi_3"))
			pp.getarr("hi_3", bc_hi_3);

		mybc = new BC::Constant(bc_hi_str, bc_lo_str,
								AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3),
								AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3));
	}

	{
		IO::ParmParse pp("ic"); // Phase-field model parameters
		pp.query("type", ic_type);
		if (ic_type == "perturbed_interface") 
		{
			ic = new IC::PerturbedInterface(geom);
			pp.queryclass(static_cast<IC::PerturbedInterface*>(ic));
		}
		else if (ic_type == "tabulated_interface")
			ic = new IC::TabulatedInterface(geom);
		else if (ic_type == "voronoi")
		{
			int total_grains = number_of_grains;
			pp.query("voronoi.number_of_grains", total_grains);
			ic = new IC::Voronoi(geom, total_grains);
		}
		else if (ic_type == "sphere")
			ic = new IC::Sphere(geom);
		else
			Util::Abort(INFO, "No valid initial condition specified");
	}

	eta_new_mf.resize(maxLevel() + 1);
	RegisterNewFab(eta_new_mf, mybc, number_of_grains, number_of_ghost_cells, "Eta",true);
	RegisterNewFab(eta_old_mf, mybc, number_of_grains, number_of_ghost_cells, "Eta old",false);

	volume = 1.0;
	RegisterIntegratedVariable(&volume, "volume");
	RegisterIntegratedVariable(&area, "area");
	RegisterIntegratedVariable(&gbenergy, "gbenergy");
	RegisterIntegratedVariable(&realgbenergy, "realgbenergy");
	RegisterIntegratedVariable(&regenergy, "regenergy");

	// Elasticity
	{
		IO::ParmParse pp("elastic");
		pp.query("on", elastic.on);
		if (elastic.on)
		{
			RegisterNodalFab(disp_mf, AMREX_SPACEDIM, 2, "disp",true);
			RegisterNodalFab(rhs_mf, AMREX_SPACEDIM, 2, "rhs",true);
			RegisterNodalFab(stress_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "stress",true);
			RegisterNodalFab(energies_mf, number_of_grains, 2, "energies",false);

			pp.query("interval", elastic.interval);
			pp.query("max_coarsening_level", elastic.max_coarsening_level);
			pp.query("tol_rel", elastic.tol_rel);
			pp.query("tol_abs", elastic.tol_abs);
			pp.query("tstart", elastic.tstart);

			pp.queryclass("op",elastic.op);

			pp.queryclass("bc",elastic.bc);


			elastic.model.resize(number_of_grains);
			for (int i = 0; i < number_of_grains; i++)
			{
				model_type mymodel;
				pp.queryclass("model",elastic.model[i]);
			}
		}
	}
}

#define ETA(i, j, k, n) eta_old(amrex::IntVect(AMREX_D_DECL(i, j, k)), n)

void PhaseFieldMicrostructure::Advance(int lev, amrex::Real time, amrex::Real dt)
{
	BL_PROFILE("PhaseFieldMicrostructure::Advance");
	/// TODO Make this optional
	//if (lev != max_level) return;
	std::swap(eta_old_mf[lev], eta_new_mf[lev]);
	const amrex::Real *DX = geom[lev].CellSize();

	Model::Interface::GB::SH gbmodel(0.0, 0.0, anisotropy.sigma0, anisotropy.sigma1);

	for (amrex::MFIter mfi(*eta_new_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		const amrex::Box &bx = mfi.tilebox();
		amrex::Array4<const amrex::Real> const &eta = (*eta_old_mf[lev]).array(mfi);
		amrex::Array4<amrex::Real> const &etanew = (*eta_new_mf[lev]).array(mfi);
		amrex::Array4<amrex::Real> const *energies;
		if (elastic.on)
		{
			amrex::Array4<amrex::Real> const &tmp_energies = (*energies_mf[lev]).array(mfi);
			energies = &tmp_energies;
		}
		
		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
			for (int m = 0; m < number_of_grains; m++)
			{
				Set::Scalar driving_force = 0.0;

				Set::Scalar kappa = NAN, mu = NAN;

				//
				// BOUNDARY TERM and SECOND ORDER REGULARIZATION
				//

				Set::Vector Deta = Numeric::Gradient(eta, i, j, k, m, DX);
				Set::Scalar normgrad = Deta.lpNorm<2>();
				if (normgrad < 1E-4)
					continue; // This ought to speed things up.

				Set::Matrix DDeta = Numeric::Hessian(eta, i, j, k, m, DX);
				Set::Scalar laplacian = DDeta.trace();

				if (!anisotropy.on || time < anisotropy.tstart)
				{
					kappa = pf.l_gb * 0.75 * pf.sigma0;
					mu = 0.75 * (1.0 / 0.23) * pf.sigma0 / pf.l_gb;
					driving_force += -kappa * laplacian;
				}
				else
				{
					Set::Vector normal = Deta / normgrad;
					Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Full> DDDDEta = Numeric::DoubleHessian<AMREX_SPACEDIM>(eta, i, j, k, m, DX);

#if AMREX_SPACEDIM == 1
					Util::Abort(INFO, "Anisotropy is enabled but works in 2D/3D ONLY");
#elif AMREX_SPACEDIM == 2
						Set::Vector tangent(normal[1],-normal[0]);
						Set::Scalar Theta = atan2(Deta(1),Deta(0));
						Set::Scalar kappa = pf.l_gb*0.75*boundary->W(Theta);
						Set::Scalar Dkappa = pf.l_gb*0.75*boundary->DW(Theta);
						Set::Scalar DDkappa = pf.l_gb*0.75*boundary->DDW(Theta);
						mu = 0.75 * (1.0/0.23) * boundary->W(Theta) / pf.l_gb;
						Set::Scalar sinTheta = sin(Theta);
						Set::Scalar cosTheta = cos(Theta);
			
						Set::Scalar Curvature_term =
							DDDDEta(0,0,0,0)*(    sinTheta*sinTheta*sinTheta*sinTheta) +
							DDDDEta(0,0,0,1)*(4.0*sinTheta*sinTheta*sinTheta*cosTheta) +
							DDDDEta(0,0,1,1)*(6.0*sinTheta*sinTheta*cosTheta*cosTheta) +
							DDDDEta(0,1,1,1)*(4.0*sinTheta*cosTheta*cosTheta*cosTheta) +
							DDDDEta(1,1,1,1)*(    cosTheta*cosTheta*cosTheta*cosTheta);

						Set::Scalar Boundary_term =
							kappa*laplacian +
							Dkappa*(cos(2.0*Theta)*DDeta(0,1) + 0.5*sin(2.0*Theta)*(DDeta(1,1) - DDeta(0,0)))
							+ 0.5*DDkappa*(sinTheta*sinTheta*DDeta(0,0) - 2.*sinTheta*cosTheta*DDeta(0,1) + cosTheta*cosTheta*DDeta(1,1));
						if (std::isnan(Boundary_term)) Util::Abort(INFO,"nan at m=",i,",",j,",",k);
			
						driving_force += - (Boundary_term) + anisotropy.beta*(Curvature_term);
						if (std::isnan(driving_force)) Util::Abort(INFO,"nan at m=",i,",",j,",",k);

#elif AMREX_SPACEDIM == 3
						// GRAHM-SCHMIDT PROCESS 
						const Set::Vector e1(1,0,0), e2(0,1,0), e3(0,0,1);
						Set::Vector _t2, _t3;
						if      (fabs(normal(0)) > fabs(normal(1)) && fabs(normal(0)) > fabs(normal(2)))
						{
	 						_t2 = e2 - normal.dot(e2)*normal; _t2 /= _t2.lpNorm<2>();
							_t3 = e3 - normal.dot(e3)*normal - _t2.dot(e3)*_t2; _t3 /= _t3.lpNorm<2>();
						}
						else if (fabs(normal(1)) > fabs(normal(0)) && fabs(normal(1)) > fabs(normal(2)))
						{
	 						_t2 = e1 - normal.dot(e1)*normal; _t2 /= _t2.lpNorm<2>();
							_t3 = e3 - normal.dot(e3)*normal - _t2.dot(e3)*_t2; _t3 /= _t3.lpNorm<2>();
						}
						else
						{
	 						_t2 = e1 - normal.dot(e1)*normal; _t2 /= _t2.lpNorm<2>();
							_t3 = e2 - normal.dot(e2)*normal - _t2.dot(e2)*_t2; _t3 /= _t3.lpNorm<2>();
						}
												
						// Compute Hessian projected into tangent space (spanned by _t1,_t2)
						Eigen::Matrix2d DDeta2D;
						DDeta2D <<
							_t2.dot(DDeta*_t2) , _t2.dot(DDeta*_t3),
							_t3.dot(DDeta*_t2) , _t3.dot(DDeta*_t3);
						Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(2);
						eigensolver.computeDirect(DDeta2D);
						Eigen::Matrix2d eigenvecs = eigensolver.eigenvectors();

						// Compute tangent vectors embedded in R^3
						Set::Vector t2 = _t2*eigenvecs(0,0) + _t3*eigenvecs(0,1),
									t3 = _t2*eigenvecs(1,0) + _t3*eigenvecs(1,1);

						// Compute components of second Hessian in t2,t3 directions
						Set::Scalar DH2 = 0.0, DH3 = 0.0;
						Set::Scalar DH23 = 0.0;
						for (int p = 0; p < 3; p++)
							for (int q = 0; q < 3; q++)
								for (int r = 0; r < 3; r++)
									for (int s = 0; s < 3; s++)
									{
										DH2 += DDDDEta(p,q,r,s)*t2(p)*t2(q)*t2(r)*t2(s);
										DH3 += DDDDEta(p,q,r,s)*t3(p)*t3(q)*t3(r)*t3(s);
										DH23 += DDDDEta(p,q,r,s)*t2(p)*t2(q)*t3(r)*t3(s);
									}

						Set::Scalar gbe = gbmodel.W(normal);
						//Set::Scalar kappa = l_gb*0.75*gbe;
						kappa = pf.l_gb*0.75*gbe;
						mu = 0.75 * (1.0/0.23) * gbe / pf.l_gb;
						Set::Scalar DDK2 = gbmodel.DDW(normal,_t2) * pf.l_gb * 0.75;
						Set::Scalar DDK3 = gbmodel.DDW(normal,_t3) * pf.l_gb * 0.75;

						// GB energy anisotropy term
						Set::Scalar gbenergy_df = - kappa*laplacian - DDK2*DDeta2D(0,0) - DDK3*DDeta2D(1,1);
						driving_force += gbenergy_df;
								  
						// Second order curvature term
						Set::Scalar reg_df = NAN;
						switch(regularization)
						{
							case Wilmore:
								reg_df = anisotropy.beta*(DH2 + DH3 + 2.0*DH23);
								break;
							case K12:
								reg_df = anisotropy.beta*(DH2+DH3);
								break;
						}
						driving_force += reg_df;

						if (std::isnan(driving_force) || std::isinf(driving_force))
						{
							for (int p = 0; p < 3; p++)
							for (int q = 0; q < 3; q++)
								for (int r = 0; r < 3; r++)
									for (int s = 0; s < 3; s++)
									{
										Util::Message(INFO,p,q,r,s," ",DDDDEta(p,q,r,s));
									}
							Util::Abort(INFO,"nan/inf detected at amrlev = ", lev," i=",i," j=",j," k=",k);
						}
#endif
				}

				//
				// CHEMICAL POTENTIAL
				//

				Set::Scalar sum_of_squares = 0.;
				for (int n = 0; n < number_of_grains; n++)
				{
					if (m == n)
						continue;
					sum_of_squares += eta(i, j, k, n) * eta(i, j, k, n);
				}
				driving_force += mu * (eta(i, j, k, m) * eta(i, j, k, m) - 1.0 + 2.0 * pf.gamma * sum_of_squares) * eta(i, j, k, m);

				//
				// SYNTHETIC DRIVING FORCE
				//
				if (lagrange.on && m == 0 && time > lagrange.tstart)
				{
					driving_force += lagrange.lambda * (volume - lagrange.vol0);
				}

				//
				// ELASTIC DRIVING FORCE
				//

				if (elastic.on && time > elastic.tstart)
				{
					driving_force += (*energies)(i,j,k,m);
				}

				//
				// EVOLVE ETA
				//
				etanew(i, j, k, m) = eta(i, j, k, m) - pf.M * dt * driving_force;
				if (std::isnan(driving_force))
					Util::Abort(INFO, i, " ", j, " ", k, " ", m);
			}
		});
	}
}

void PhaseFieldMicrostructure::Initialize(int lev)
{
	BL_PROFILE("PhaseFieldMicrostructure::Initialize");
	eta_new_mf[lev]->setVal(0.0);
	eta_old_mf[lev]->setVal(0.0);

	ic->Initialize(lev, eta_new_mf);
	ic->Initialize(lev, eta_old_mf);

	if (elastic.on)
	{
		disp_mf[lev].get()->setVal(0.0);
		rhs_mf[lev].get()->setVal(0.0);
		//res_mf[lev].get()->setVal(0.0);
		stress_mf[lev].get()->setVal(0.0);
	}
}

void PhaseFieldMicrostructure::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/)
{
	BL_PROFILE("PhaseFieldMicrostructure::TagCellsForRefinement");
	const amrex::Real *DX = geom[lev].CellSize();
	const Set::Vector dx(DX);
	const Set::Scalar dxnorm = dx.lpNorm<2>();

	for (amrex::MFIter mfi(*eta_new_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		const amrex::Box &bx = mfi.tilebox();
		amrex::Array4<const amrex::Real> const &etanew = (*eta_new_mf[lev]).array(mfi);
		amrex::Array4<char> const &tags = a_tags.array(mfi);

		for (int n = 0; n < number_of_grains; n++)
			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
				Set::Vector grad = Numeric::Gradient(etanew, i, j, k, n, DX);

				if (dxnorm * grad.lpNorm<2>() > ref_threshold)
					tags(i, j, k) = amrex::TagBox::SET;
			});
	}
}

void PhaseFieldMicrostructure::TimeStepComplete(amrex::Real /*time*/, int /*iter*/)
{
	// TODO: remove this function, it is no longer needed.
}

void PhaseFieldMicrostructure::TimeStepBegin(amrex::Real time, int iter)
{
	if (!elastic.on) return;
	if (time < elastic.tstart)   return;
	if (iter % elastic.interval) return;


	if (finest_level != rhs_mf.size() - 1)
	{
		Util::Abort(INFO, "amr.max_level is larger than necessary. Set to ", finest_level, " or less");
	}
	for (int lev = 0; lev < rhs_mf.size(); lev++)
		rhs_mf[lev]->setVal(0.0);

	//Operator::Elastic<model_type> elasticop;
	elastic.op.SetUniform(false);
	amrex::LPInfo info;
	//info.setMaxCoarseningLevel(0);
	elastic.op.define(geom, grids, dmap, info);

	// Set linear elastic model
	amrex::Vector<amrex::FabArray<amrex::BaseFab<model_type>>> model_mf;
	model_mf.resize(disp_mf.size());
	for (int lev = 0; lev < rhs_mf.size(); ++lev)
	{
		amrex::Box domain(geom[lev].Domain());
		domain.convert(amrex::IntVect::TheNodeVector());
		model_mf[lev].define(disp_mf[lev]->boxArray(), disp_mf[lev]->DistributionMap(), 1, 2);
		//model_mf[lev].setVal(mymodel);

		eta_new_mf[lev]->FillBoundary();

		Set::Vector DX(geom[lev].CellSize());

		for (MFIter mfi(model_mf[lev], false); mfi.isValid(); ++mfi)
		{
			amrex::Box bx = mfi.growntilebox(2);

			amrex::Array4<model_type> const &model = model_mf[lev].array(mfi);
			amrex::Array4<const Set::Scalar> const &eta = eta_new_mf[lev]->array(mfi);

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
				model(i, j, k) = elastic.model[0] * 0.0;
				for (int n = 0; n < number_of_grains; n++)
					model(i, j, k) += elastic.model[n] * eta(i, j, k, n);
			});
		}

		Util::RealFillBoundary(model_mf[lev],elastic.op.Geom(lev));
	}
	elastic.op.SetModel(model_mf);

	elastic.bc.Init(rhs_mf,geom);
	elastic.op.SetBC(&elastic.bc);

	Solver::Nonlocal::Linear linearsolver(elastic.op);
	IO::ParmParse pp("elastic");
	pp.queryclass("solver",linearsolver);
	linearsolver.solve(disp_mf, rhs_mf);

	for (int lev = 0; lev < disp_mf.size(); lev++)
	{
		elastic.op.Stress(lev, *stress_mf[lev], *disp_mf[lev]);
		elastic.op.Energy(lev,*energies_mf[lev],*disp_mf[lev],elastic.model);
	}
}

void PhaseFieldMicrostructure::Integrate(int amrlev, Set::Scalar time, int /*step*/,
										 const amrex::MFIter &mfi, const amrex::Box &box)
{
	Model::Interface::GB::SH gbmodel(0.0, 0.0, anisotropy.sigma0, anisotropy.sigma1);

	BL_PROFILE("PhaseFieldMicrostructure::Integrate");
	const amrex::Real *DX = geom[amrlev].CellSize();
	amrex::Array4<amrex::Real> const &eta = (*eta_new_mf[amrlev]).array(mfi);
	amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
		Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);

		volume += eta(i, j, k, 0) * dv;

		Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
		Set::Scalar normgrad = grad.lpNorm<2>();

		if (normgrad > 1E-8)
		{
			Set::Vector normal = grad / normgrad;

			Set::Scalar da = normgrad * dv;
			area += da;

			if (!anisotropy.on || time < anisotropy.tstart)
			{
				gbenergy += pf.sigma0 * da;

				Set::Scalar k = 0.75 * pf.sigma0 * pf.l_gb;
				realgbenergy += 0.5 * k * normgrad * normgrad * dv;
				regenergy = 0.0;
			}
			else
			{
#if AMREX_SPACEDIM == 2
				Set::Scalar theta = atan2(grad(1), grad(0));
				Set::Scalar sigma = boundary->W(theta);
				gbenergy += sigma * da;

				Set::Scalar k = 0.75 * sigma * pf.l_gb;
				realgbenergy += 0.5 * k * normgrad * normgrad * dv;

				Set::Matrix DDeta = Numeric::Hessian(eta, i, j, k, 0, DX);
				Set::Vector tangent(normal[1], -normal[0]);
				Set::Scalar k2 = (DDeta * tangent).dot(tangent);
				regenergy += 0.5 * anisotropy.beta * k2 * k2;
#elif AMREX_SPACEDIM == 3
				gbenergy += gbmodel.W(normal) * da;
#endif
			}
		}
	});
}

} // namespace Integrator


#include <eigen3/Eigen/Eigenvalues>

#include <omp.h>
#include <cmath>
#include <chrono>
#include <random>

#include <AMReX_SPACE.H>
#include <AMReX_ParallelReduce.H>

#include "PhaseFieldMicrostructure.H"
#include "BC/Constant.H"
#include "BC/Step.H"
#include "Set/Set.H"
#include "Util/Util.H"
#include "IC/Random.H"
#include "IC/Trig.H"
#include "IC/Sphere.H"
#include "Model/Interface/GB/SH.H"
#include "Numeric/Stencil.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"
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
		pp.query("elastic_threshold",pf.elastic_threshold);
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
		IO::ParmParse pp("anisotropy"); // Phase-field model parameters
		pp.query("on", anisotropy.on);

		pp.query("beta", anisotropy.beta);
		pp.query("tstart", anisotropy.tstart);
		anisotropy.timestep = timestep;
		pp.query("timestep", anisotropy.timestep);
		anisotropy.plot_int = plot_int;
		pp.query("plot_int", anisotropy.plot_int);
		anisotropy.plot_dt = plot_dt;
		pp.query("plot_dt", anisotropy.plot_dt);
		pp.query("thermo_int", anisotropy.thermo_int);
		pp.query("thermo_plot_int", anisotropy.thermo_plot_int);
		pp.query("elastic_int",anisotropy.elastic_int);

		std::map<std::string, RegularizationType> regularization_type;
		regularization_type["wilmore"] = RegularizationType::Wilmore;
		regularization_type["k12"] = RegularizationType::K12;
		std::string regularization_type_input = "k12";
		pp.query("regularization", regularization_type_input);
		regularization = regularization_type[regularization_type_input];

		pp.query("gb_type", gb_type);
		if (gb_type == "abssin")
		{
			boundary = new Model::Interface::GB::AbsSin();
			pp.queryclass(*static_cast<Model::Interface::GB::AbsSin *>(boundary));
		}
		else if (gb_type == "sin")
		{
			boundary = new Model::Interface::GB::Sin();
			pp.queryclass(*static_cast<Model::Interface::GB::Sin *>(boundary));
			
		}
		else if (gb_type == "read")
		{
			boundary = new Model::Interface::GB::Read();
			pp.queryclass(*static_cast<Model::Interface::GB::Read *>(boundary));
		}

		else if (gb_type == "sh")
		{
			Util::Assert(INFO, TEST(AMREX_SPACEDIM == 3));
			boundary = new Model::Interface::GB::SH();
			pp.queryclass(*static_cast<Model::Interface::GB::SH *>(boundary));
		}
		else if (anisotropy.on)
		{
			Util::Abort(INFO,"A GB model must be specified");
		}
	}

	{
		IO::ParmParse pp("fluctuation");
		pp.query("on",fluctuation.on);
		pp.query("amp",fluctuation.amp);
		pp.query("sd",fluctuation.sd);
		pp.query("tstart", fluctuation.tstart);
		fluctuation.norm_dist = std::normal_distribution<double>(0.0,fluctuation.sd);
		
	}

	{
		IO::ParmParse pp("disconnection");
		pp.query("on",disconnection.on);
		pp.query("tstart", disconnection.tstart);
		pp.query("range",disconnection.range);
		pp.query("nucleation_energy",disconnection.nucleation_energy);
		pp.query("temp",disconnection.temp);
		pp.query("box_size",disconnection.box_size);
		pp.query("interval",disconnection.interval);

		disconnection.unif_dist = std::uniform_real_distribution<double>(0.0,1.0);
		disconnection.int_dist = std::uniform_int_distribution<int>(0,1);
		disconnection.p = exp(-disconnection.nucleation_energy/(disconnection.K_b*disconnection.temp));
	}

	{
		IO::ParmParse pp("bc");
		std::string bc_type = "constant";
		pp.query("eta.type",bc_type);
		if (bc_type == "constant")
		{
			mybc = new BC::Constant(number_of_grains);
			pp.queryclass("eta",*static_cast<BC::Constant *>(mybc));
		}
		else if (bc_type == "step")
		{
			mybc = new BC::Step();
			pp.queryclass("eta",*static_cast<BC::Step *>(mybc));
		}
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
		{
			Util::Assert(INFO,TEST(number_of_grains == 2));
			ic = new IC::TabulatedInterface(geom);
			pp.queryclass(static_cast<IC::TabulatedInterface*>(ic));
		}
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
	
 	RegisterNewFab(fluct_mf, new BC::Nothing(), 1, number_of_ghost_cells, "fluct",true);
	RegisterNewFab(disc_mf, new BC::Nothing(), 1, number_of_ghost_cells, "disc",true);  // see box
	
	volume = 1.0;
	RegisterIntegratedVariable(&volume, "volume");
	RegisterIntegratedVariable(&area, "area");
	RegisterIntegratedVariable(&gbenergy, "gbenergy");
	RegisterIntegratedVariable(&realgbenergy, "realgbenergy");
	RegisterIntegratedVariable(&regenergy, "regenergy");
	RegisterIntegratedVariable(&elastic.strainenergy, "strainenergy");
	RegisterIntegratedVariable(&elastic.force, "force");
	RegisterIntegratedVariable(&elastic.disp, "disp");

	RegisterGeneralFab(model_mf,1,2);

	// Elasticity
	{
		IO::ParmParse pp("elastic");
		pp.query("on", elastic.on);
		if (elastic.on)
		{
			RegisterNodalFab(disp_mf, AMREX_SPACEDIM, 2, "disp",true);
			RegisterNodalFab(rhs_mf, AMREX_SPACEDIM, 2, "rhs",true);
			RegisterNodalFab(stress_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "stress",true);
			RegisterNodalFab(energy_mf, 1, 2, "energy",true);

			pp.query("interval", elastic.interval);
			pp.query("max_coarsening_level", elastic.max_coarsening_level);
			pp.query("tol_rel", elastic.tol_rel);
			pp.query("tol_abs", elastic.tol_abs);
			pp.query("tstart", elastic.tstart);

			pp.queryclass("bc",elastic.bc);


			elastic.model.resize(number_of_grains);
			for (int i = 0; i < number_of_grains; i++)
			{
				model_type mymodel;
				pp.queryclass("model",elastic.model[i]);
			}
			pp.queryclass("model1",elastic.model[0]);
			pp.queryclass("model2",elastic.model[1]);
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
		amrex::Array4<amrex::Real> const &fluct = (*fluct_mf[lev]).array(mfi);

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
				// EVOLVE ETA
				//
				etanew(i, j, k, m) = eta(i, j, k, m) - pf.M * dt * driving_force;

				//
				// FLUCTUATION TERM
				//

				if (fluctuation.on && time > fluctuation.tstart)
				{
					fluct(i,j,k,0) = fluctuation.amp * fluctuation.norm_dist(fluctuation.rand_num_gen) / DX[0];
					etanew(i,j,k,m) += fluct(i,j,k,0) * dt;
				}

				if (std::isnan(driving_force))
					Util::Abort(INFO, i, " ", j, " ", k, " ", m);
			}
		});

		//
		// ELASTIC DRIVING FORCE
		//
		if (elastic.on && time > elastic.tstart)
		{
			const amrex::Box &bx = mfi.tilebox();
			amrex::Array4<const amrex::Real> const &eta = (*eta_old_mf[lev]).array(mfi);
			amrex::Array4<amrex::Real> const &etanew = (*eta_new_mf[lev]).array(mfi);
			amrex::Array4<const amrex::Real> const &sigma = (*stress_mf[lev]).array(mfi);

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
			{
				for (int m = 0; m < number_of_grains; m++)
				{
					Set::Scalar driving_force = 0.0;
					Set::Scalar etasum = 0.0;
					Set::Matrix F0avg = Set::Matrix::Zero();
					
					for (int n = 0; n < number_of_grains; n++)
					{
						etasum += eta(i,j,k,n);
						F0avg  += eta(i,j,k,n) * elastic.model[n].F0;
					} 
					
					Set::Matrix dF0deta = elastic.model[m].F0;//(etasum * elastic.model[m].F0 - F0avg) / (etasum * etasum);


					Set::Matrix sig;
					#if AMREX_SPACEDIM == 2
					sig(0,0) = Numeric::Interpolate::CellToNodeAverage(sigma,i,j,k,0);
					sig(0,1) = Numeric::Interpolate::CellToNodeAverage(sigma,i,j,k,1);
					sig(1,0) = Numeric::Interpolate::CellToNodeAverage(sigma,i,j,k,2);
					sig(1,1) = Numeric::Interpolate::CellToNodeAverage(sigma,i,j,k,3);
					#elif AMREX_SPACEDIM == 3
					sig(0,0) = Numeric::Interpolate::CellToNodeAverage(sigma,i,j,k,0);
					sig(0,1) = Numeric::Interpolate::CellToNodeAverage(sigma,i,j,k,1);
					sig(0,2) = Numeric::Interpolate::CellToNodeAverage(sigma,i,j,k,2);
					sig(1,0) = Numeric::Interpolate::CellToNodeAverage(sigma,i,j,k,3);
					sig(1,1) = Numeric::Interpolate::CellToNodeAverage(sigma,i,j,k,4);
					sig(1,2) = Numeric::Interpolate::CellToNodeAverage(sigma,i,j,k,5);
					sig(2,0) = Numeric::Interpolate::CellToNodeAverage(sigma,i,j,k,6);
					sig(2,1) = Numeric::Interpolate::CellToNodeAverage(sigma,i,j,k,7);
					sig(2,2) = Numeric::Interpolate::CellToNodeAverage(sigma,i,j,k,8);
					#endif

					Set::Scalar tmpdf = (dF0deta.transpose() * sig).trace();

					if (tmpdf > pf.elastic_threshold)
					{
						driving_force -= pf.elastic_mult * (tmpdf-pf.elastic_threshold);
					}
					else if (tmpdf < -pf.elastic_threshold)
					{
						driving_force -= pf.elastic_mult * (tmpdf+pf.elastic_threshold);
					}

					etanew(i, j, k, m) -= pf.M * dt * driving_force;
				}
			});

		}
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
	//
	// Manual Disconnection Nucleation
	//

	if (disconnection.on && time > disconnection.tstart && !(iter % disconnection.interval))
	{
		disconnection.nucleation_sites.clear();

		// iterate over all AMR levels
		//for (int lev = 0; lev <= max_level; lev ++)  
		{
			int lev = max_level;
			const amrex::Real *DX = geom[lev].CellSize();
		
			// iterate over all Patches
			for (amrex::MFIter mfi(*eta_new_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
			{

				const amrex::Box &bx = mfi.tilebox();
				amrex::Array4<amrex::Real> const &eta = (*eta_old_mf[lev]).array(mfi);

				// iterate over the GRID (index i,j,k)
				amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

					if (eta(i,j,k,0) >= (0.5 - disconnection.range) && eta(i,j,k,0) <= (0.5 + disconnection.range))
					{
						amrex::Real q = 0.0; // fixed nucleation site
						//amrex::Real q = disconnection.unif_dist(disconnection.rand_num_gen);
						
						if (q < disconnection.p)
						{
							Set::Vector x;

							AMREX_D_TERM(
								x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * DX[0];,
								x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * DX[1];,
								x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k)) * DX[2];
							);

							if (x(0) == 4.0 && x(1) == 0.0){ // fixed nucleation site
								disconnection.nucleation_sites.push_back(x);
							};

							//disconnection.nucleation_sites.push_back(x);
						}
					}
				});
			}
		}
		for (int lev = 0; lev <= max_level; lev++)
		{
			const amrex::Real *DX = geom[lev].CellSize();
			for (amrex::MFIter mfi(*eta_new_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
			{
				const amrex::Box &bx = mfi.tilebox();
				amrex::Array4<amrex::Real> const &etanew = (*eta_new_mf[lev]).array(mfi);
				amrex::Array4<amrex::Real> const &disc = (*disc_mf[lev]).array(mfi);

				//iterate again
				amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
					Set::Vector x;

					AMREX_D_TERM(
						x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * DX[0];,
						x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * DX[1];,
						x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k)) * DX[2];
					);

					for (unsigned int m = 0; m < disconnection.nucleation_sites.size(); m++){	
						amrex::Real r_squared = 0;
						for (int n = 0; n < AMREX_SPACEDIM; n++){
							amrex::Real dist = disconnection.nucleation_sites[m](n) - x(n);
							r_squared += dist * dist;

							if (sqrt(r_squared) > disconnection.box_size / 2) break;
							if (n == AMREX_SPACEDIM - 1){
								//int phase = disconnection.int_dist(disconnection.rand_num_gen);
								int phase = 1; // for fixed nucleation site, 0 for up, 1 for down
								amrex::Real bump = exp(1 - 1 / (1 - 2/disconnection.box_size * r_squared));
								
								disc(i,j,k,0) = bump;
								etanew(i,j,k,phase) = bump * (1-etanew(i,j,k,phase)) + etanew(i,j,k,phase);
								etanew(i,j,k,1-phase) = 1 - etanew(i,j,k,phase);
							}
						}
					}
				});
			}
		}

		// Note: this code will need to be implemented before
		// we can run in parallel.
		//std::vector<Set::Vector> buffer;
		//amrex::ParallelAllGather::AllGather(disconnection.nucleation_sites.data(),
		//						 disconnection.nucleation_sites.size(),
		//						 buffer.data(),
		//						 amrex::ParallelContext::CommunicatorAll());
	}

	// 
	// Elastic solve
	//
	if (anisotropy.on && time >= anisotropy.tstart)
	{
		SetTimestep(anisotropy.timestep);
		if (anisotropy.elastic_int > 0) 
			if (iter % anisotropy.elastic_int) return;
	}
	
	if (!elastic.on) return;
	if (time < elastic.tstart)   return;
	if (iter % elastic.interval) return;

	if (finest_level != rhs_mf.size() - 1)
	{
		Util::Abort(INFO, "amr.max_level is larger than necessary. Set to ", finest_level, " or less");
	}
	for (int lev = 0; lev < rhs_mf.size(); lev++)
		rhs_mf[lev]->setVal(0.0);

	Operator::Elastic<model_type::sym> elasticop;
	elasticop.SetUniform(false);
	amrex::LPInfo info;
	//info.setMaxCoarseningLevel(0);
	elasticop.define(geom, grids, dmap, info);

	// Set linear elastic model
	for (int lev = 0; lev < rhs_mf.size(); ++lev)
	{
		amrex::Box domain(geom[lev].Domain());
		domain.convert(amrex::IntVect::TheNodeVector());

		eta_new_mf[lev]->FillBoundary();

		Set::Vector DX(geom[lev].CellSize());

		for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
		{
			amrex::Box bx = mfi.grownnodaltilebox(-1,2);

			amrex::Array4<model_type> const &model = model_mf[lev]->array(mfi);
			amrex::Array4<const Set::Scalar> const &eta = eta_new_mf[lev]->array(mfi);

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
				std::vector<Set::Scalar> etas(number_of_grains);
				for (int n = 0; n < number_of_grains; n++) etas[n] = 0.25*(eta(i,j,k,n) + eta(i,j-1,k,n) + eta(i-1,j,k,n) + eta(i-1,j-1,k,n));
				model(i, j, k) = model_type::Combine(elastic.model,etas);
			});
		}

		Util::RealFillBoundary(*model_mf[lev],elasticop.Geom(lev));
	}

	elastic.bc.SetTime(time);
	elastic.bc.Init(rhs_mf,geom);
	elasticop.SetBC(&elastic.bc);

	Solver::Nonlocal::Newton<model_type> linearsolver(elasticop);
	IO::ParmParse pp("elastic");
	pp.queryclass("solver",linearsolver);
	linearsolver.solve(disp_mf, rhs_mf, model_mf, 1E-8, 1E-8);

	linearsolver.W(energy_mf,disp_mf,model_mf);
	linearsolver.DW(stress_mf,disp_mf,model_mf);
}

void PhaseFieldMicrostructure::Integrate(int amrlev, Set::Scalar time, int /*step*/,
										 const amrex::MFIter &mfi, const amrex::Box &box)
{
	Model::Interface::GB::SH gbmodel(0.0, 0.0, anisotropy.sigma0, anisotropy.sigma1);
	const amrex::Real *DX = geom[amrlev].CellSize();
	Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);

	BL_PROFILE("PhaseFieldMicrostructure::Integrate");
	amrex::Array4<amrex::Real> const &eta = (*eta_new_mf[amrlev]).array(mfi);
	amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

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
	if (elastic.on)
	{
		amrex::Array4<amrex::Real> const &w        = (*energy_mf[amrlev]).array(mfi);
		amrex::Array4<amrex::Real> const &stress   = (*stress_mf[amrlev]).array(mfi);
		amrex::Array4<amrex::Real> const &u        = (*disp_mf[amrlev])  .array(mfi);
		amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
		{
			if (j == geom[amrlev].Domain().hiVect()[1])
			{
				elastic.force += 0.5*(stress(i,j+1,k,1) + stress(i+1,j+1,k,1)) * DX[0];
				elastic.disp  += 0.5*(u(i,j+1,k,0)      + u(i+1,j+1,k,0)     ) * DX[0];
			}
			elastic.strainenergy += 0.25 * (w(i,j,k) + w(i+1,j,k) + w(i,j+1,k) + w(i+1,j+1,k)) * volume;
		});
	}
}

} // namespace Integrator

#include "PhaseFieldMicrostructure.H"
#include "BC/Constant.H"
#include "Set/Set.H"
#include "Util/Util.H"
namespace Integrator
{
PhaseFieldMicrostructure::PhaseFieldMicrostructure() : Integrator()
{
	//
	// READ INPUT PARAMETERS
	//
	{
		amrex::ParmParse pp("pf"); // Phase-field model parameters
		pp.query("number_of_grains", number_of_grains);
		pp.query("M", M);
		if (pp.contains("mu")) pp.query("mu", mu);
		pp.query("gamma", gamma);
		pp.query("sigma0", sigma0);
		pp.query("l_gb", l_gb);
		pp.query("sdf_on", sdf_on);
		if (sdf_on)
		{
			sdf.resize(number_of_grains);
			pp.queryarr("sdf",sdf);
		}
	}
	{
		amrex::ParmParse pp("amr");
		pp.query("max_level",max_level);
	}
	{
		amrex::Real theta0,sigma0,sigma1;

		amrex::ParmParse pp("anisotropy"); // Phase-field model parameters
		pp.query("on", anisotropy);
		pp.query("theta0", theta0);
		pp.query("filename", filename);
		pp.query("gb_type", gb_type);
		theta0 *= 0.01745329251; // convert degrees into radians
		pp.query("sigma0", sigma0);
		pp.query("sigma1", sigma1);
		pp.query("beta", beta);
		pp.query("tstart", anisotropy_tstart);


		if(gb_type=="abssin")
			boundary = new Model::Interface::GrainBoundary::AbsSin(theta0,sigma0,sigma1);
		else if(gb_type=="sin")
			boundary = new Model::Interface::GrainBoundary::Sin(theta0,sigma0,sigma1);
		else if(gb_type=="read")
			boundary = new Model::Interface::GrainBoundary::Read(filename);
		else
			boundary = new Model::Interface::GrainBoundary::Sin(theta0,sigma0,sigma1);

    
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

		mybc = new BC::Constant(bc_hi_str, bc_lo_str,
					AMREX_D_DECL(bc_lo_1,bc_lo_2,bc_lo_3),
					AMREX_D_DECL(bc_hi_1,bc_hi_2,bc_hi_3));
	}


	{
		amrex::ParmParse pp("ic"); // Phase-field model parameters
		pp.query("type", ic_type);
		if (ic_type == "perturbed_interface")
			ic = new IC::PerturbedInterface(geom);
		else if (ic_type == "tabulated_interface")
			ic = new IC::TabulatedInterface(geom);
		else if (ic_type == "voronoi")
			ic = new IC::Voronoi(geom,number_of_grains);
		else if (ic_type == "circle")
			ic = new IC::Circle(geom);
		else
			Util::Abort(INFO, "No valid initial condition specified");
	}
	/*
	 */
  
	eta_new_mf.resize(maxLevel()+1);

	// amrex::ParallelDescriptor::Barrier();
	// amrex::Abort("No error, just exiting here");

	RegisterNewFab(eta_new_mf, mybc, number_of_grains, number_of_ghost_cells, "Eta");
	//eta_old_mf.resize(maxLevel()+1);
	RegisterNewFab(eta_old_mf, mybc, number_of_grains, number_of_ghost_cells, "Eta old");
	RegisterNewFab(etas_mf, 1, "Etas");

	volume = 10;
	RegisterIntegratedVariable(&volume, "volume");
  
	// Elasticity
	{
		amrex::ParmParse pp("elastic");
		pp.query("on",elastic.on);

		if (elastic.on)
		{
			pp.query("interval",elastic.interval);
			//pp.query("type",elastic_type);
			pp.query("max_iter",elastic.max_iter);
			pp.query("max_fmg_iter",elastic.max_fmg_iter);
			pp.query("verbose",elastic.verbose);
			pp.query("cgverbose",elastic.cgverbose);
			pp.query("tol_rel",elastic.tol_rel);
			pp.query("tol_abs",elastic.tol_abs);
			pp.query("tstart",elastic.tstart);


			amrex::ParmParse pp_bc("elastic.bc");

			// Read in boundary types as strings, then convert to Operator::Elastic::BC types and store for future use.
			amrex::Vector<std::string> AMREX_D_DECL(bctype_xlo_str,bctype_ylo_str,bctype_zlo_str);
			amrex::Vector<std::string> AMREX_D_DECL(bctype_xhi_str,bctype_yhi_str,bctype_zhi_str);
			AMREX_D_TERM(pp_bc.queryarr("type_xlo",bctype_xlo_str);, pp_bc.queryarr("type_ylo",bctype_ylo_str);, pp_bc.queryarr("type_zlo",bctype_zlo_str););
			AMREX_D_TERM(pp_bc.queryarr("type_xhi",bctype_xhi_str);, pp_bc.queryarr("type_yhi",bctype_yhi_str);, pp_bc.queryarr("type_zhi",bctype_zhi_str););
			if ( AMREX_D_TERM(bctype_xlo_str.size() < AMREX_SPACEDIM, || bctype_ylo_str.size() < AMREX_SPACEDIM, || bctype_zlo_str.size() < AMREX_SPACEDIM)
			     ||
			     AMREX_D_TERM(bctype_xhi_str.size() < AMREX_SPACEDIM, || bctype_yhi_str.size() < AMREX_SPACEDIM, || bctype_zhi_str.size() < AMREX_SPACEDIM))
				Util::Abort(INFO, "incorrect number of terms specified in bctype");
			std::map<std::string,Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC> bc;
			bc["displacement"] = Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC::Displacement;
			bc["disp"]         = Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC::Displacement;
			bc["traction"]     = Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC::Traction;
			bc["trac"]         = Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC::Traction;
			bc["periodic"]     = Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC::Periodic;
			AMREX_D_TERM(elastic.bctype_xlo = {AMREX_D_DECL(bc[bctype_xlo_str[0]], bc[bctype_xlo_str[1]], bc[bctype_xlo_str[2]])};,
				     elastic.bctype_ylo = {AMREX_D_DECL(bc[bctype_ylo_str[0]], bc[bctype_ylo_str[1]], bc[bctype_ylo_str[2]])};,
				     elastic.bctype_zlo = {AMREX_D_DECL(bc[bctype_zlo_str[0]], bc[bctype_zlo_str[1]], bc[bctype_zlo_str[2]])};);
			AMREX_D_TERM(elastic.bctype_xhi = {AMREX_D_DECL(bc[bctype_xhi_str[0]], bc[bctype_xhi_str[1]], bc[bctype_xhi_str[2]])};,
				     elastic.bctype_yhi = {AMREX_D_DECL(bc[bctype_yhi_str[0]], bc[bctype_yhi_str[1]], bc[bctype_yhi_str[2]])};,
				     elastic.bctype_zhi = {AMREX_D_DECL(bc[bctype_zhi_str[0]], bc[bctype_zhi_str[1]], bc[bctype_zhi_str[2]])};);


			AMREX_D_TERM(pp_bc.queryarr("xlo",elastic.bc_xlo);, pp_bc.queryarr("ylo",elastic.bc_ylo);, pp_bc.queryarr("zlo",elastic.bc_zlo););
			AMREX_D_TERM(pp_bc.queryarr("xhi",elastic.bc_xhi);, pp_bc.queryarr("yhi",elastic.bc_yhi);, pp_bc.queryarr("zhi",elastic.bc_zhi););


			RegisterNodalFab(displacement, AMREX_SPACEDIM,"disp");
			RegisterNodalFab(body_force, AMREX_SPACEDIM,"rhs");
			//RegisterNodalFab(strain, AMREX_SPACEDIM,"strain");
			RegisterNodalFab(stress, AMREX_SPACEDIM*AMREX_SPACEDIM,"stress");

		}

		//
		// Initialize elastic models
		//
		// for (int n = 0; n <  number_of_grains; n++) 
		// 	models.push_back(new OperatorCell::Elastic::PolyCrystal::Cubic(107.3, 60.9, 28.30)); // randomized angles

	}
}


#define ETA(i,j,k,n) eta_old(amrex::IntVect(AMREX_D_DECL(i,j,k)),n)

void
PhaseFieldMicrostructure::Advance (int lev, amrex::Real time, amrex::Real dt)
{

	/// TODO Make this optional
	if (lev != max_level) return;
	std::swap(eta_old_mf[lev], eta_new_mf[lev]);
	const amrex::Real* DX = geom[lev].CellSize();

	static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
					   dy(AMREX_D_DECL(0,1,0)),
					   dz(AMREX_D_DECL(0,0,1)));

	for ( amrex::MFIter mfi(*eta_new_mf[lev],true); mfi.isValid(); ++mfi )
	{
		const amrex::Box& bx = mfi.tilebox();
		FArrayBox &eta_new_box     = (*eta_new_mf[lev])[mfi];
		FArrayBox &eta_old     = (*eta_old_mf[lev])[mfi];

		AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
			     for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
			     for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
		{
			amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
			
			for (int i = 0; i < number_of_grains; i++)
			{
		  
				//
				// CHEMICAL POTENTIAL
				//

				amrex::Real sum_of_squares = 0.;
				for (int j = 0; j < number_of_grains; j++)
				{
					if (i==j) continue;
					sum_of_squares += eta_old(m,j)*eta_old(m,j);
				}

				//
				// BOUNDARY TERM
				// and
				// CURVATURE PENALTY
				//

				// amrex::Real grad1 =  grad1_normal;
				// amrex::Real grad2 =  grad2_normal;
				// amrex::Real grad12 = (eta_old(m+dx+dy,i) - eta_old(m-dx+dy,i) - eta_old(m+dx-dy,i) + eta_old(m-dx-dy))/(4.*DX[0]*DX[1]);
				amrex::Real grad11 = (eta_old(m+dx,i) - 2.*eta_old(m,i) + eta_old(m-dx,i))/DX[0]/DX[0]; // 3 point
				amrex::Real grad22 = (eta_old(m+dy,i) - 2.*eta_old(m,i) + eta_old(m-dy,i))/DX[1]/DX[1]; // 3 point
		      
				amrex::Real laplacian = grad11 + grad22;

				amrex::Real kappa = l_gb*0.75*sigma0;
				mu = 0.75 * (1.0/0.23) * sigma0 / l_gb;

				if (anisotropy && time > anisotropy_tstart)
				{
#if AMREX_SPACEDIM == 2
					amrex::Real grad1 = (eta_old(m+dx,i) - eta_old(m-dx,i))/(2*DX[0]);
					amrex::Real grad2 = (eta_old(m+dy,i) - eta_old(m-dy,i))/(2*DX[1]);

					amrex::Real grad12 = (eta_old(m+dx+dy,i) - eta_old(m-dx+dy,i) - eta_old(m+dx-dy,i) + eta_old(m-dx-dy))/(4.*DX[0]*DX[1]);

					amrex::Real Theta = atan2(grad2,grad1);
					amrex::Real Kappa = l_gb*0.75*boundary->W(Theta);
					amrex::Real DKappa = l_gb*0.75*boundary->DW(Theta);
					amrex::Real DDKappa = l_gb*0.75*boundary->DDW(Theta);
					amrex::Real Mu = 0.75 * (1.0/0.23) * boundary->W(Theta) / l_gb;
					amrex::Real Sine_theta = sin(Theta);
					amrex::Real Cos_theta = cos(Theta);
		
					// amrex::Real norm_grad = grad1*grad1+grad2*grad2; //(UNUSED)
		
					amrex::Real grad1111=(1/(DX[0]*DX[0]*DX[0]*DX[0])) * ((ETA(m1+2,m2,m3,i)) - 4*(ETA(m1+1,m2,m3,i))+6*(ETA(m1,m2,m3,i)) - 4*(ETA(m1-1,m2,m3,i)) + (ETA(m1-2,m2,m3,i)));
			
					amrex::Real grad2222=(1/(DX[1]*DX[1]*DX[1]*DX[1])) * ((ETA(m1,m2+2,m3,i)) - 4*(ETA(m1,m2+1,m3,i))+6*(ETA(m1,m2,m3,i)) - 4*(ETA(m1,m2-1,m3,i)) + (ETA(m1,m2-2,m3,i)));
			
					amrex::Real grad1112= (1/(24*DX[0]*DX[0]*DX[0]*DX[1]))*
						((-ETA(m1+2,m2+2,m3,i)+8*ETA(m1+2,m2+1,m3,i)
						  -8*ETA(m1+2,m2-1,m3,i)+ETA(m1+2,m2-2,m3,i))
						 -2*(-ETA(m1+1,m2+2,m3,i)+8*ETA(m1+1,m2+1,m3,i)
						     -8*ETA(m1+1,m2-1,m3,i)+ETA(m1+1,m2-2,m3,i))
						 +2*(-ETA(m1-1,m2+2,m3,i)+8*ETA(m1-1,m2+1,m3,i)
						     -8*ETA(m1-1,m2-1,m3,i)+ETA(m1-1,m2-2,m3,i))
						 -(-ETA(m1-2,m2+2,m3,i)+8*ETA(m1-2,m2+1,m3,i)
						   -8*ETA(m1-2,m2-1,m3,i)+ETA(m1-2,m2-2,m3,i)));
			
					amrex::Real grad1222= (1/(24*DX[0]*DX[1]*DX[1]*DX[1]))*
						((-ETA(m1+2,m2+2,m3,i)+8*ETA(m1+1,m2+2,m3,i)
						  -8*ETA(m1-1,m2+2,m3,i)+ETA(m1-2,m2+2,m3,i))
						 -2*(-ETA(m1+2,m2+1,m3,i)+8*ETA(m1+1,m2+1,m3,i)
						     -8*ETA(m1-1,m2+1,m3,i)+ETA(m1-2,m2+1,m3,i))
						 +2*(-ETA(m1+2,m2-1,m3,i)+8*ETA(m1+1,m2-1,m3,i)
						     -8*ETA(m1-1,m2-1,m3,i)+ETA(m1-2,m2-1,m3,i))
						 -(-ETA(m1+2,m2-2,m3,i)+8*ETA(m1+1,m2-2,m3,i)
						   -8*ETA(m1-1,m2-2,m3,i)+ETA(m1-2,m2-2,m3,i)));
			
					amrex::Real grad1122= (1/(144*DX[0]*DX[0]*DX[1]*DX[1]))*
						(-(-ETA(m1+2,m2+2,m3,i)+16*ETA(m1+1,m2+2,m3,i)-30*ETA(m1,m2+2,m3,i)
						   +16*ETA(m1-1,m2+2,m3,i)-ETA(m1-2,m2+2,m3,i))
						 +16*(-ETA(m1+2,m2+1,m3,i)+16*ETA(m1+1,m2+1,m3,i)-30*ETA(m1,m2+1,m3,i)
						      +16*ETA(m1-1,m2+1,m3,i)-ETA(m1-2,m2+1,m3,i))
						 -30*(-ETA(m1+2,m2,m3,i)+16*ETA(m1+1,m2,m3,i)-30*ETA(m1,m2,m3,i)
						      +16*ETA(m1-1,m2,m3,i)-ETA(m1-2,m2,m3,i))
						 +16*(-ETA(m1+2,m2-1,m3,i)+16*ETA(m1+1,m2-1,m3,i)-30*ETA(m1,m2-1,m3,i)
						      +16*ETA(m1-1,m2-1,m3,i)-ETA(m1-2,m2-1,m3,i))
						 -(-ETA(m1+2,m2-2,m3,i)+16*ETA(m1+1,m2-2,m3,i)-30*ETA(m1,m2-2,m3,i)
						   +16*ETA(m1-1,m2-2,m3,i)-ETA(m1-2,m2-2,m3,i)));
			
					amrex::Real Curvature_term =
						grad1111*(Sine_theta*Sine_theta*Sine_theta*Sine_theta)
						+grad1112*(-4*Sine_theta*Sine_theta*Sine_theta*Cos_theta)
						+grad1122*(6*Sine_theta*Sine_theta*Cos_theta*Cos_theta)
						+grad1222*(-4*Sine_theta*Cos_theta*Cos_theta*Cos_theta)
						+grad2222*(Cos_theta*Cos_theta*Cos_theta*Cos_theta);

					amrex::Real W =
						Mu*(ETA(m1,m2,m3,i)*ETA(m1,m2,m3,i)
						    - 1.0 + 2.0*gamma*sum_of_squares)*ETA(m1,m2,m3,i);

					amrex::Real Boundary_term =
						Kappa*laplacian +
						DKappa*(cos(2.0*Theta)*grad12 + 0.5*sin(2.0*Theta)*(grad22-grad11))
						+ 0.5*DDKappa*(Sine_theta*Sine_theta*grad11 - 2.*Sine_theta*Cos_theta*grad12 + Cos_theta*Cos_theta*grad22);
			
			
					eta_new_box(amrex::IntVect(AMREX_D_DECL(m1,m2,m3)),i) =
						ETA(m1,m2,m3,i) -
						M*dt*(W - (Boundary_term) + beta*(Curvature_term));
#else
					Util::Abort(INFO, "Anisotropy is enabled but works in 2D ONLY");
#endif
				}
				else // Isotropic response if less than anisotropy_tstart
				{
					eta_new_box(amrex::IntVect(AMREX_D_DECL(m1,m2,m3)),i) =
						ETA(m1,m2,m3,i) -
						M*dt*(mu*(ETA(m1,m2,m3,i)*ETA(m1,m2,m3,i)
							  - 1.0 +
							  2.0*gamma*sum_of_squares)*ETA(m1,m2,m3,i)
						      - kappa*laplacian);
				}

				//
				// SYNTHETIC DRIVING FORCE
				//
				if (sdf_on)
				{
					eta_new_box(amrex::IntVect(AMREX_D_DECL(m1,m2,m3)),i) -=  M*dt*(sdf[i]);
				}

				//
				// ELASTIC DRIVING FORCE
				//
				Util::Warning(INFO,"elastic driving force disabled!");
				// if (elastic.on && time > elastic.tstart)
				// {
				// 	FArrayBox &energiesfab     = (*energies[lev])[mfi];

				// 	eta_new_box(amrex::IntVect(AMREX_D_DECL(m1,m2,m3)),i)
				// 		-= M*dt*( energiesfab(amrex::IntVect(AMREX_D_DECL(m1,m2,m3)),i));
				// }

			}
	      
		}
	}

}

void
PhaseFieldMicrostructure::Initialize (int lev)
{
	ic->Initialize(lev,eta_new_mf);
	ic->Initialize(lev,eta_old_mf);
  
  
	if (elastic.on)
	{
		// //displacement[lev].get()->setVal(0.0); ///TODO
		// strain[lev].get()->setVal(0.0); 
		// stress[lev].get()->setVal(0.0); 
		// stress_vm[lev].get()->setVal(0.0);
		// //body_force[lev].get()->setVal(0.0);
		// energy[lev].get()->setVal(0.0); 
		// energies[lev].get()->setVal(0.0); 
	}
}


void
PhaseFieldMicrostructure::TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{
	const amrex::Real* DX      = geom[lev].CellSize();

	amrex::Vector<int>  itags;

	for (amrex::MFIter mfi(*eta_new_mf[lev],true); mfi.isValid(); ++mfi)
	{
		const amrex::Box&  bx  = mfi.tilebox();
		amrex::TagBox&     tag  = tags[mfi];
		amrex::BaseFab<amrex::Real> &eta_new_box = (*eta_new_mf[lev])[mfi];

#if BL_SPACEDIM==2
		for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
			for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
			{
				for (int n = 0; n < number_of_grains; n++)
				{
					amrex::Real gradx = (eta_new_box(amrex::IntVect(i+1,j),n) - eta_new_box(amrex::IntVect(i-1,j),n))/(2.*DX[0]);
					amrex::Real grady = (eta_new_box(amrex::IntVect(i,j+1),n) - eta_new_box(amrex::IntVect(i,j-1),n))/(2.*DX[1]);
					if (DX[0]*sqrt(gradx*gradx + grady*grady)>0.1) tag(amrex::IntVect(i,j)) = amrex::TagBox::SET;
				}
			}

#elif BL_SPACEDIM==3

		for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
			for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
				for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
				{
					for (int n = 0; n < number_of_grains; n++)
					{
						amrex::Real gradx = (eta_new_box(amrex::IntVect(i+1,j,k),n) - eta_new_box(amrex::IntVect(i-1,j,k),n))/(2.*DX[0]);
						amrex::Real grady = (eta_new_box(amrex::IntVect(i,j+1,k),n) - eta_new_box(amrex::IntVect(i,j-1,k),n))/(2.*DX[1]);
						amrex::Real gradz = (eta_new_box(amrex::IntVect(i,j,k+1),n) - eta_new_box(amrex::IntVect(i,j,k-1),n))/(2.*DX[2]);
						if (DX[0]*sqrt(gradx*gradx + grady*grady + gradz*gradz)>0.1) tag(amrex::IntVect(i,j,k)) = amrex::TagBox::SET;
					}
				}
#endif

	}
}


void PhaseFieldMicrostructure::TimeStepComplete(amrex::Real /*time*/, int iter)
{
	if (!(iter % plot_int))
	{
		Util::Warning(INFO,"etas calculation diabled");
		// for (int ilev = 0; ilev < displacement.size(); ilev++)
		// {
		// 	for ( amrex::MFIter mfi(*strain[ilev],true); mfi.isValid(); ++mfi )
		// 	{
		// 		const amrex::Box& box = mfi.tilebox();
		// 		amrex::FArrayBox &etas  = (*etas_mf[ilev])[mfi];
		// 		amrex::FArrayBox &eta_new  = (*eta_new_mf[ilev])[mfi];
		// 		AMREX_D_TERM(for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++),
		// 			     for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++),
		// 			     for (int k = box.loVect()[2]; k<=box.hiVect()[2]; k++))
		// 		{
		// 			amrex::IntVect m(AMREX_D_DECL(i,j,k));
					
		// 			etas(m) = 0.0;
		// 			for (int n = 0; n < number_of_grains; n++)
		// 				etas(m) += ((Set::Scalar)n)*eta_new(m,n);
		// 		}
		// 	}
		// }

	}
}


void PhaseFieldMicrostructure::TimeStepBegin(amrex::Real time, int iter)
{
	if (anisotropy && time > anisotropy_tstart)
	{
		SetTimestep(anisotropy_timestep);
	}

	if (!elastic.on) return;
	if (iter%elastic.interval) return;
	if (time < elastic.tstart) return;

	LPInfo info;
	info.setAgglomeration(true);
	info.setConsolidation(true);
	int max_mg_level = 0;
	info.setMaxCoarseningLevel(max_mg_level);



	int nlevels = maxLevel() + 1;
	amrex::Vector<amrex::BoxArray> 	ngrids;
	amrex::Vector<amrex::FabArray<amrex::BaseFab<Model::Solid::LinearElastic::Cubic> > > model;

	ngrids.resize(nlevels);
	// displacement.resize(nlevels);
	// body_force.resize(nlevels);
	model.resize(nlevels);

	Model::Solid::LinearElastic::Cubic grain1; grain1.Randomize(); //(10.73, 6.09, 2.830); //testmodel.Randomize();
	Model::Solid::LinearElastic::Cubic grain2 = grain1 * 2.0;
	// Util::Message(INFO,grain1);
	// Util::Message(INFO,grain2);
	

	// amrex::Vector<std::unique_ptr<amrex::MultiFab> > residual;
	// residual.resize(nlevels);

	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		//ngeom[ilev].define(amrex::convert(geom[ilev],amrex::IntVect::TheNodeVector()));
		
		ngrids[ilev] = grids[ilev];
		ngrids[ilev].convert(amrex::IntVect::TheNodeVector());

		// displacement[ilev].reset(new amrex::MultiFab(ngrids[ilev],dmap[ilev],AMREX_SPACEDIM,1));
		// body_force[ilev]  .reset(new amrex::MultiFab(ngrids[ilev],dmap[ilev],AMREX_SPACEDIM,1));
		// residual[ilev]    .reset(new amrex::MultiFab(ngrids[ilev],dmap[ilev],AMREX_SPACEDIM,1));

		model[ilev].define(ngrids[ilev],dmap[ilev],1,1);

		Util::Message(INFO,"displacement size = " , displacement.size());
		Util::Message(INFO,"body_force size = " , body_force.size());
		displacement[ilev]->setVal(0.0);
		body_force[ilev]->setVal(0.000000001);

		//model[ilev].setVal(grain1); //TODO


		for (amrex::MFIter mfi(*body_force[ilev],true); mfi.isValid(); ++mfi)
		{
		 	const amrex::Box& box = mfi.tilebox();

			amrex::BaseFab<amrex::Real> &etafab = (*(eta_old_mf[ilev]))[mfi];
		 	amrex::BaseFab<amrex::Real> &rhsfab = (*(body_force[ilev]))[mfi];
			amrex::BaseFab<Model::Solid::LinearElastic::Cubic> &modelfab = (model[ilev])[mfi];

			AMREX_D_TERM(for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++),
				     for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++),
				     for (int k = box.loVect()[2]; k<=box.hiVect()[2]; k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));

				// Set boundary conditions in RHS
				for (int p = 0; p<AMREX_SPACEDIM; p++)
				{
					AMREX_D_TERM( if (i == geom[ilev].Domain().loVect()[0]) rhsfab(m,p)   = elastic.bc_xlo[p];,
						      if (j == geom[ilev].Domain().loVect()[1]) rhsfab(m,p)   = elastic.bc_ylo[p];,
						      if (k == geom[ilev].Domain().loVect()[2]) rhsfab(m,p)   = elastic.bc_zlo[p]; );
					AMREX_D_TERM( if (i == geom[ilev].Domain().hiVect()[0]+1) rhsfab(m,p) = elastic.bc_xhi[p];,
						      if (j == geom[ilev].Domain().hiVect()[1]+1) rhsfab(m,p) = elastic.bc_yhi[p];,
						      if (k == geom[ilev].Domain().hiVect()[2]+1) rhsfab(m,p) = elastic.bc_zhi[p]; );
				}

				
				modelfab(m) =
					grain1 * 0.25*(etafab(m-dx[0]-dx[1], 0) + etafab(m-dx[0], 0) + etafab(m-dx[1], 0) + etafab(m, 0)) +
					grain2 * 0.25*(etafab(m-dx[0]-dx[1], 1) + etafab(m-dx[0], 1) + etafab(m-dx[1], 1) + etafab(m, 1));



				
			}
		}
	}
	



	elastic.op = new Operator::Elastic<Model::Solid::LinearElastic::Cubic>();
	// amrex::Vector<amrex::Geometry> elgeom(nlevels);
	// for (int i = 0; i < nlevels; i++)
	// {

	//        elgeom = geom;
	elastic.op->define(geom,grids,dmap,info);

	elastic.op->SetBC({{AMREX_D_DECL(elastic.bctype_xlo,elastic.bctype_ylo,elastic.bctype_zlo)}},
			  {{AMREX_D_DECL(elastic.bctype_xhi,elastic.bctype_yhi,elastic.bctype_zhi)}});

	for (int ilev = 0; ilev < nlevels; ++ilev) elastic.op->SetModel(ilev,model[ilev]);

	amrex::MLMG solver(*elastic.op);
	solver.setMaxIter(elastic.max_iter);
	solver.setMaxFmgIter(elastic.max_fmg_iter);
	solver.setVerbose(elastic.verbose);
	solver.setCGVerbose(elastic.cgverbose);

      	solver.solve(GetVecOfPtrs(displacement),
	  	     GetVecOfConstPtrs(body_force),
	  	     elastic.tol_rel,
	  	     elastic.tol_abs);

 	for (int lev = 0; lev < displacement.size(); lev++)
 	{
 		const amrex::Real* DX = geom[lev].CellSize();

 		elastic.op->Stress(lev,*stress[lev],*displacement[lev]);
 		//elastic_operator->Energy(lev,*energy[lev],*displacement[lev]);

// 		for ( amrex::MFIter mfi(*strain[lev],true); mfi.isValid(); ++mfi )
// 		{
// 			const Box& bx = mfi.tilebox();

// 			FArrayBox &ufab  = (*displacement[lev])[mfi];
// 			FArrayBox &epsfab  = (*strain[lev])[mfi];
// 			FArrayBox &sigmafab  = (*stress[lev])[mfi];
// 			//FArrayBox &energyfab  = (*energy[lev])[mfi];
// 			FArrayBox &energiesfab  = (*energies[lev])[mfi];
// 			FArrayBox &sigmavmfab  = (*stress_vm[lev])[mfi];


// 			AMREX_D_TERM(for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++),
// 				     for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++),
// 				     for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++))
// 			 	{
// 			 		amrex::IntVect m(AMREX_D_DECL(i,j,k));
// #if AMREX_SPACEDIM == 2
// 			 		epsfab(m,0) = (ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0]);
// 			 		epsfab(m,1) = 0.5*(ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]) +
// 			 		 	0.5*(ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1]);
// 					epsfab(m,2) = epsfab(m,1);
// 			 		epsfab(m,3) = (ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1]);
// #elif AMREX_SPACEDIM == 3
// 			 		epsfab(m,0) = (ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0]);
// 			 		epsfab(m,1) = 0.5*(ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]) + 0.5*(ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1]);
// 			 		epsfab(m,2) = 0.5*(ufab(m+dx,2) - ufab(m-dx,2))/(2.0*DX[0]) + 0.5*(ufab(m+dz,0) - ufab(m-dz,0))/(2.0*DX[2]);
// 					epsfab(m,3) = epsfab(m,1);
// 					epsfab(m,4) = (ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1]);
// 					epsfab(m,5) = 0.5*(ufab(m+dy,2) - ufab(m-dy,2))/(2.0*DX[1]) + 0.5*(ufab(m+dz,1) - ufab(m-dz,1))/(2.0*DX[2]);
// 					epsfab(m,6) = epsfab(m,2);
// 					epsfab(m,7) = epsfab(m,5);
// 			 		epsfab(m,8) = (ufab(m+dz,2) - ufab(m-dz,2))/(2.0*DX[2]);
// #endif

// 					// elastic_operator->Stress((*stress[lev])[mfi],(*displacement[lev])[mfi],lev,mfi);
// 					// elastic_operator->Energy((*stress[lev])[mfi],(*displacement[lev])[mfi],lev,mfi);


// 			 		sigmavmfab(m,0) = 0.0;

// #if AMREX_SPACEDIM == 2
// 			 		sigmavmfab(m) =
// 			 		 	sqrt(0.5*(sigmafab(m,0) - sigmafab(m,1)*(sigmafab(m,0) - sigmafab(m,1))
// 			 		 		  + sigmafab(m,0)*sigmafab(m,0)
// 			 		 		  + sigmafab(m,1)*sigmafab(m,1)
// 			 		 		  + 6.0*sigmafab(m,2)*sigmafab(m,2)));
// #elif AMREX_SPACEDIM == 3
// 			 		sigmavmfab(m) =
// 			 		 	sqrt(0.5*((sigmafab(m,0) - sigmafab(m,4))*(sigmafab(m,0) - sigmafab(m,4)) +
// 							  (sigmafab(m,4) - sigmafab(m,8))*(sigmafab(m,4) - sigmafab(m,8)) +
// 							  (sigmafab(m,8) - sigmafab(m,0))*(sigmafab(m,8) - sigmafab(m,0)))+
// 						     + 3.0 * (sigmafab(m,1)*sigmafab(m,1) +
// 							      sigmafab(m,2)*sigmafab(m,2) +
// 							      sigmafab(m,5)*sigmafab(m,5)));
// #endif

// 			 	}

			
// 			//elastic_operator->Energies(energiesfab,ufab,lev,mfi);
// 		}
 	}
}
}


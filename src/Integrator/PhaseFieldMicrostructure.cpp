#include "PhaseFieldMicrostructure.H"
#include "BC/Constant.H"
#include "Set/Set.H"
#include "Util/Util.H"
#include "IC/Affine.H"
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
		anisotropy_timestep = timestep;
		pp.query("timestep",anisotropy_timestep);


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
		pp.query("on",elastic_on);
		pp.query("int",elastic_int);
		pp.query("type",elastic_type);
		pp.query("max_iter",elastic_max_iter);
		pp.query("max_fmg_iter",elastic_max_fmg_iter);
		pp.query("verbose",elastic_verbose);
		pp.query("cgverbose",elastic_cgverbose);
		pp.query("tol_rel",elastic_tol_rel);
		pp.query("tol_abs",elastic_tol_abs);
		pp.query("tstart",elastic_tstart);

		pp.queryarr("load_t",elastic_load_t);
		pp.queryarr("load_disp",elastic_load_disp);
		if (elastic_load_t.size() != elastic_load_disp.size())
			Util::Abort(INFO, "load_t and load_disp must have the same number of entries");

		if (elastic_on)
		{
			//RegisterNewFab(displacement, mybc,AMREX_SPACEDIM,1,"u");
			//RegisterNewFab(body_force,mybc,AMREX_SPACEDIM,0,"b");
			RegisterNewFab(strain,    mybc,AMREX_SPACEDIM*AMREX_SPACEDIM,0,"eps");
			RegisterNewFab(stress,    mybc,AMREX_SPACEDIM*AMREX_SPACEDIM,0,"sig");
			RegisterNewFab(stress_vm, mybc,1,0,"sig_VM");
			RegisterNewFab(energy,    mybc,1,0,"W");
			RegisterNewFab(energies,  mybc,number_of_grains,0,"W");


			amrex::ParmParse pp("elastic.bc");
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

			elastic_bc = new BC::Constant(bc_hi_str, bc_lo_str,
						      AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3),
						      AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3));
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
				if (elastic_on && time > elastic_tstart)
				{
					FArrayBox &energiesfab     = (*energies[lev])[mfi];

					eta_new_box(amrex::IntVect(AMREX_D_DECL(m1,m2,m3)),i)
						-= M*dt*( energiesfab(amrex::IntVect(AMREX_D_DECL(m1,m2,m3)),i));
				}

			}
	      
		}
	}

}

void
PhaseFieldMicrostructure::Initialize (int lev)
{
	ic->Initialize(lev,eta_new_mf);
	ic->Initialize(lev,eta_old_mf);
  
  
	if (elastic_on)
	{
		//displacement[lev].get()->setVal(0.0);
		strain[lev].get()->setVal(0.0); 
		stress[lev].get()->setVal(0.0); 
		stress_vm[lev].get()->setVal(0.0);
		//body_force[lev].get()->setVal(0.0);
		energy[lev].get()->setVal(0.0); 
		energies[lev].get()->setVal(0.0); 
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
		for (int ilev = 0; ilev < displacement.size(); ilev++)
		{
			for ( amrex::MFIter mfi(*strain[ilev],true); mfi.isValid(); ++mfi )
			{
				const amrex::Box& box = mfi.tilebox();
				amrex::FArrayBox &etas  = (*etas_mf[ilev])[mfi];
				amrex::FArrayBox &eta_new  = (*eta_new_mf[ilev])[mfi];
				AMREX_D_TERM(for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++),
					     for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++),
					     for (int k = box.loVect()[2]; k<=box.hiVect()[2]; k++))
				{
					amrex::IntVect m(AMREX_D_DECL(i,j,k));
					
					etas(m) = 0.0;
					for (int n = 0; n < number_of_grains; n++)
						etas(m) += ((Set::Scalar)n)*eta_new(m,n);
				}
			}
		}

	}

	// Set::Scalar volume = 0.0;
	// for (int ilev = 0; ilev < max_level; ilev++)
	// {
	// 	const amrex::Real* DX = geom[ilev].CellSize();

	// 	const BoxArray& cfba = amrex::coarsen(grids[ilev+1], refRatio(ilev));

	// 	for ( amrex::MFIter mfi(grids[ilev],dmap[ilev],true); mfi.isValid(); ++mfi )
	// 	{
	// 		const amrex::Box& box = mfi.tilebox();
	// 		const::BoxArray & comp = amrex::complementIn(box,cfba);
	// 		amrex::FArrayBox &eta_new  = (*eta_new_mf[ilev])[mfi];

	// 		for (int i = 0; i < comp.size(); i++)
	// 		{
	// 			amrex::Box mybox = comp[i];
	// 			AMREX_D_TERM(for (int m1 = mybox.loVect()[0]; m1<=mybox.hiVect()[0]; m1++),
	// 				     for (int m2 = mybox.loVect()[1]; m2<=mybox.hiVect()[1]; m2++),
	// 				     for (int m3 = mybox.loVect()[2]; m3<=mybox.hiVect()[2]; m3++))
	// 			{
	// 				amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
	// 				volume += AMREX_D_TERM(DX[0],*DX[1],*DX[2]);
	// 			}
	// 		}
	// 	}
	// }
	// {
	// 	const amrex::Real* DX = geom[max_level].CellSize();

	// 	for ( amrex::MFIter mfi(grids[max_level],dmap[max_level],true); mfi.isValid(); ++mfi )
	// 	{
	// 		const amrex::Box& box = mfi.tilebox();
	// 		//amrex::FArrayBox &etas  = (*etas_mf[ilev])[mfi];
	// 		amrex::FArrayBox &eta_new  = (*eta_new_mf[max_level])[mfi];
	// 		AMREX_D_TERM(for (int m1 = box.loVect()[0]; m1<=box.hiVect()[0]; m1++),
	// 			     for (int m2 = box.loVect()[1]; m2<=box.hiVect()[1]; m2++),
	// 			     for (int m3 = box.loVect()[2]; m3<=box.hiVect()[2]; m3++))
	// 		{
	// 			amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
					
	// 			volume += AMREX_D_TERM(DX[0],*DX[1],*DX[2]);
				
	// 		}
	// 	}
	// }

	// std::cout << "volbefore = " << volume << std::endl;
	// amrex::ParallelDescriptor::ReduceRealSum(volume);

	// std::cout << "volafter = " << volume << std::endl;
	// Util::Message(INFO,"volume = ", volume);
	// //(*etas_mf[max_level]).setVal(max_level);
	
}


void PhaseFieldMicrostructure::TimeStepBegin(amrex::Real time, int iter)
{
	if (anisotropy && time > anisotropy_tstart)
	{
		SetTimestep(anisotropy_timestep);
	}

	if (!elastic_on) return;
	if (iter%elastic_int) return;
	if (time < elastic_tstart) return;

	LPInfo info;
	info.setAgglomeration(true);
	info.setConsolidation(true);
	int max_mg_level = 0;
	info.setMaxCoarseningLevel(max_mg_level);



	int nlevels = maxLevel() + 1;

	//amrex::Vector<amrex::Geometry> 	ngeom;
	amrex::Vector<amrex::BoxArray> 	ngrids;
	//ngeom.resize(nlevels);
	ngrids.resize(nlevels);
	displacement.resize(nlevels);
	body_force.resize(nlevels);
	amrex::Vector<amrex::FabArray<amrex::BaseFab<Model::Solid::LinearElastic::Cubic> > > modelfab;
	modelfab.resize(nlevels);
	Model::Solid::LinearElastic::Cubic testmodel;
	testmodel.Randomize();

	amrex::Vector<std::unique_ptr<amrex::MultiFab> > residual;
	residual.resize(nlevels);

	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		//ngeom[ilev].define(amrex::convert(geom[ilev],amrex::IntVect::TheNodeVector()));
		
		ngrids[ilev] = grids[ilev];
		ngrids[ilev].convert(amrex::IntVect::TheNodeVector());

		displacement[ilev].reset(new amrex::MultiFab(ngrids[ilev],dmap[ilev],AMREX_SPACEDIM,1));
		body_force[ilev]  .reset(new amrex::MultiFab(ngrids[ilev],dmap[ilev],AMREX_SPACEDIM,1));
		residual[ilev].reset(new amrex::MultiFab(ngrids[ilev],dmap[ilev],AMREX_SPACEDIM,1));

		modelfab[ilev].define(ngrids[ilev],dmap[ilev],1,1);

		displacement[ilev]->setVal(0.0);
		body_force[ilev]->setVal(0.0,0);
		body_force[ilev]->setVal(0.00001,1);
		modelfab[ilev].setVal(testmodel);


		for (amrex::MFIter mfi(*body_force[ilev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.tilebox();

			amrex::BaseFab<amrex::Real> &rhsfab = (*(body_force[ilev]))[mfi];

			AMREX_D_TERM(for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++),
				     for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++),
				     for (int k = box.loVect()[2]; k<=box.hiVect()[2]; k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));

				for (int p = 0; p<AMREX_SPACEDIM; p++)
				{
					if (j == geom[ilev].Domain().hiVect()[1]) rhsfab(m,p) = 0.1;
				}
			}
		}
	}
	
	Set::Vector n(1.0,0);
	Set::Vector b(1.0,0);
	Set::Scalar alpha = 1.0;
	Set::Scalar m = 2.0;
	IC::Affine ic(geom,n,alpha,b,true,m);
	ic.SetComp(0);
	for (int ilev = 0; ilev < nlevels; ilev++) ic.Initialize(ilev,displacement);

		
	//for (int ilev = 0; ilev < nlevels; ilev++) res[ilev].minus(rhs[ilev], 0, 2, 0);




	elastic_operator = new Operator::Elastic<Model::Solid::LinearElastic::Cubic>();
	elastic_operator->define(geom,grids,dmap,info);



	std::array<Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC,AMREX_SPACEDIM>
		bc_x_lo = {AMREX_D_DECL(Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC::Displacement,
					Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC::Displacement,
					Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC::Displacement)};
	std::array<Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC,AMREX_SPACEDIM>
		bc_x_hi = {AMREX_D_DECL(Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC::Displacement,
					Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC::Displacement,
					Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC::Displacement)};
	std::array<Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC,AMREX_SPACEDIM>
		bc_y_lo = {AMREX_D_DECL(Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC::Displacement,
					Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC::Displacement,
					Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC::Displacement)};
	std::array<Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC,AMREX_SPACEDIM>
		bc_y_hi = {AMREX_D_DECL(Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC::Displacement,
					Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC::Displacement,
					Operator::Elastic<Model::Solid::LinearElastic::Cubic>::BC::Displacement)};
	
	elastic_operator->SetBC({{AMREX_D_DECL(bc_x_lo,bc_y_lo,bc_z_lo)}},
				{{AMREX_D_DECL(bc_x_hi,bc_y_hi,bc_z_hi)}});

	for (int ilev = 0; ilev < nlevels; ++ilev) elastic_operator->SetModel(ilev,modelfab[ilev]);

	amrex::MLMG solver(*elastic_operator);
	solver.setMaxIter(elastic_max_iter);
	solver.setMaxFmgIter(elastic_max_fmg_iter);
	solver.setVerbose(elastic_verbose);
	solver.setCGVerbose(elastic_cgverbose);

	
	
	for (int ilev = 0; ilev < nlevels; ilev++)
	{
		Util::Message(INFO,"level = ", ilev, " norm = ", body_force[ilev]->norm0());
		elastic_operator->FApply(ilev,0,*(body_force[ilev]),*(displacement[ilev]));
		Util::Message(INFO,"level = ", ilev, " norm = ", body_force[ilev]->norm0());
		residual[ilev]->setVal(0.0);

	}
	elastic_operator->BuildMasks();

	for (int ilev = nlevels-1; ilev > 0; ilev--)
	{
		elastic_operator->Reflux(0, *residual[ilev-1], *displacement[ilev-1], *body_force[ilev-1], *residual[ilev], *displacement[ilev], *body_force[ilev]);
	}

	// solver.solve(GetVecOfPtrs(displacement),
	//  	     GetVecOfConstPtrs(body_force),
	//  	     elastic_tol_rel,
	//  	     elastic_tol_abs);



	amrex::Vector<std::string> complete_name_array;
	complete_name_array.push_back("1");
	complete_name_array.push_back("2");
	amrex::WriteMultiLevelPlotfile("residual", nlevels, amrex::GetVecOfConstPtrs(residual), complete_name_array, Geom(), t_new[0], istep, refRatio());
	amrex::WriteMultiLevelPlotfile("bodyforce", nlevels, amrex::GetVecOfConstPtrs(body_force), complete_name_array, Geom(), t_new[0], istep, refRatio());



	for (int lev = 0; lev < displacement.size(); lev++)
	{
		const amrex::Real* DX = geom[lev].CellSize();
		const amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
						  dy(AMREX_D_DECL(0,1,0)),
						  dz(AMREX_D_DECL(0,0,1)));

		// elastic_operator->Stress(lev,*stress[lev],*displacement[lev]);
		// elastic_operator->Energy(lev,*energy[lev],*displacement[lev]);

		for ( amrex::MFIter mfi(*strain[lev],true); mfi.isValid(); ++mfi )
		{
			const Box& bx = mfi.tilebox();

			FArrayBox &ufab  = (*displacement[lev])[mfi];
			FArrayBox &epsfab  = (*strain[lev])[mfi];
			FArrayBox &sigmafab  = (*stress[lev])[mfi];
			//FArrayBox &energyfab  = (*energy[lev])[mfi];
			FArrayBox &energiesfab  = (*energies[lev])[mfi];
			FArrayBox &sigmavmfab  = (*stress_vm[lev])[mfi];


			AMREX_D_TERM(for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++),
				     for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++),
				     for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++))
			 	{
			 		amrex::IntVect m(AMREX_D_DECL(i,j,k));
#if AMREX_SPACEDIM == 2
			 		epsfab(m,0) = (ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0]);
			 		epsfab(m,1) = 0.5*(ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]) +
			 		 	0.5*(ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1]);
					epsfab(m,2) = epsfab(m,1);
			 		epsfab(m,3) = (ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1]);
#elif AMREX_SPACEDIM == 3
			 		epsfab(m,0) = (ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0]);
			 		epsfab(m,1) = 0.5*(ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]) + 0.5*(ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1]);
			 		epsfab(m,2) = 0.5*(ufab(m+dx,2) - ufab(m-dx,2))/(2.0*DX[0]) + 0.5*(ufab(m+dz,0) - ufab(m-dz,0))/(2.0*DX[2]);
					epsfab(m,3) = epsfab(m,1);
					epsfab(m,4) = (ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1]);
					epsfab(m,5) = 0.5*(ufab(m+dy,2) - ufab(m-dy,2))/(2.0*DX[1]) + 0.5*(ufab(m+dz,1) - ufab(m-dz,1))/(2.0*DX[2]);
					epsfab(m,6) = epsfab(m,2);
					epsfab(m,7) = epsfab(m,5);
			 		epsfab(m,8) = (ufab(m+dz,2) - ufab(m-dz,2))/(2.0*DX[2]);
#endif

					// elastic_operator->Stress((*stress[lev])[mfi],(*displacement[lev])[mfi],lev,mfi);
					// elastic_operator->Energy((*stress[lev])[mfi],(*displacement[lev])[mfi],lev,mfi);


			 		sigmavmfab(m,0) = 0.0;

#if AMREX_SPACEDIM == 2
			 		sigmavmfab(m) =
			 		 	sqrt(0.5*(sigmafab(m,0) - sigmafab(m,1)*(sigmafab(m,0) - sigmafab(m,1))
			 		 		  + sigmafab(m,0)*sigmafab(m,0)
			 		 		  + sigmafab(m,1)*sigmafab(m,1)
			 		 		  + 6.0*sigmafab(m,2)*sigmafab(m,2)));
#elif AMREX_SPACEDIM == 3
			 		sigmavmfab(m) =
			 		 	sqrt(0.5*((sigmafab(m,0) - sigmafab(m,4))*(sigmafab(m,0) - sigmafab(m,4)) +
							  (sigmafab(m,4) - sigmafab(m,8))*(sigmafab(m,4) - sigmafab(m,8)) +
							  (sigmafab(m,8) - sigmafab(m,0))*(sigmafab(m,8) - sigmafab(m,0)))+
						     + 3.0 * (sigmafab(m,1)*sigmafab(m,1) +
							      sigmafab(m,2)*sigmafab(m,2) +
							      sigmafab(m,5)*sigmafab(m,5)));
#endif

			 	}

			
			//elastic_operator->Energies(energiesfab,ufab,lev,mfi);
		}
	}
}


void
PhaseFieldMicrostructure::Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/,
				    const amrex::MFIter &mfi, const amrex::Box &box)
{
	const amrex::Real* DX = geom[amrlev].CellSize();

	amrex::FArrayBox &eta_new  = (*eta_new_mf[amrlev])[mfi];

	
	AMREX_D_TERM(for (int m1 = box.loVect()[0]; m1<=box.hiVect()[0]; m1++),
		     for (int m2 = box.loVect()[1]; m2<=box.hiVect()[1]; m2++),
		     for (int m3 = box.loVect()[2]; m3<=box.hiVect()[2]; m3++))
	{
		amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
		volume += eta_new(m,0)*AMREX_D_TERM(DX[0],*DX[1],*DX[2]);
	}

}

}

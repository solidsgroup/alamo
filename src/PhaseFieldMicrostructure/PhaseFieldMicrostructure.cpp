#include "PhaseFieldMicrostructure/PhaseFieldMicrostructure.H"

#if BL_SPACEDIM == 2

PhaseFieldMicrostructure::PhaseFieldMicrostructure() :
  GeneralAMRIntegrator(), 
  mybc(geom),
  linop(geom,grids,dmap,info)
{

  //
  // READ INPUT PARAMETERS
  //

  {
    amrex::ParmParse pp("pf"); // Phase-field model parameters
    pp.query("number_of_grains", number_of_grains);
    pp.query("M", M);
    pp.query("mu", mu);
    pp.query("gamma", gamma);
    pp.query("sigma0", sigma0);
    pp.query("l_gb", l_gb);
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
    pp.query("damp", damp);
    pp.query("tstart", anisotropy_tstart);

    // if (amrex::Verbose()) std::cout << "should only print if run with -v flag" << std::cout;
  
    if(gb_type=="abssin") 
 	boundary = new PFBoundaryAbsSin(theta0,sigma0,sigma1);
    else if(gb_type=="sin")
	boundary = new PFBoundarySin(theta0,sigma0,sigma1);
    else if(gb_type=="read")
	boundary = new PFBoundaryRead(filename);
    else
	boundary = new PFBoundarySin(theta0,sigma0,sigma1);
 

    
    // if(ParallelDescriptor::IOProcessor())
    //   if (!boundary->Test()) amrex::Error("Boundary functor does not pass derivative test");
  }

  {
    amrex::ParmParse pp("ic"); // Phase-field model parameters
    pp.query("type", ic_type);
    if (ic_type == "perturbed_interface")
      ic = new PerturbedInterface(geom);
  }

  RegisterNewFab(eta_new, mybc, number_of_grains, number_of_ghost_cells, "Eta");
  RegisterNewFab(eta_old, mybc, number_of_grains, number_of_ghost_cells, "Eta old");
  /// \todo Replace `mybc` with new BC object
  RegisterNewFab(displacement, mybc, AMREX_SPACEDIM, 1, "u");
  RegisterNewFab(body_force, mybc, AMREX_SPACEDIM, 1, "b");
}


#define ETA(i,j,k,n) eta_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n)

void
PhaseFieldMicrostructure::Advance (int lev, Real /*time*/, Real dt)
{
  std::swap(eta_old[lev], eta_new[lev]);
  const Real* dx = geom[lev].CellSize();

  for ( amrex::MFIter mfi(*eta_new[lev],true); mfi.isValid(); ++mfi )
    {
      const amrex::Box& bx = mfi.tilebox();

      amrex::BaseFab<Real> &eta_new_box = (*eta_new[lev])[mfi];
      amrex::BaseFab<Real> &eta_old_box = (*eta_old[lev])[mfi];

      for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
	for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
#if AMREX_SPACEDIM > 2
	  for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
#endif
	  {
	    for (int m = 0; m < number_of_grains; m++)
	      {

		amrex::Real sum_of_squares = 0.;
		for (int n = 0; n < number_of_grains; n++)
		  {
		    if (n==m) continue;
		    sum_of_squares += ETA(i,j,k,n)*ETA(i,j,k,n);
		  }
	       
	       
		amrex::Real grad1_normal = (ETA(i+1,j,k,m) - ETA(i-1,j,k,m))/(2*dx[0]);
	       
		amrex::Real grad2_normal = (ETA(i,j+1,k,m) - ETA(i,j-1,k,m))/(2*dx[1]);

		amrex::Real grad1 =  grad1_normal;

		amrex::Real grad2 =  grad2_normal;
	       
		amrex::Real grad12 = (ETA(i+1,j+1,k,m) - ETA(i-1,j+1,k,m)
					     - ETA(i+1,j-1,k,m) + ETA(i-1,j-1,k,m))/(4.*dx[0]*dx[1]);
	       
		amrex::Real grad11 =  (ETA(i+1,j,k,m) - 2.*ETA(i,j,k,m)
					      + ETA(i-1,j,k,m))/dx[0]/dx[0]; // 3 point
	       
		amrex::Real grad22 = (ETA(i,j+1,k,m) - 2.*ETA(i,j,k,m)
					     + ETA(i,j-1,k,m))/dx[1]/dx[1]; // 3 point
	       
		//Second curvature variational terms
		amrex::Real grad1111=(1/((dx[0]*dx[0]*dx[0]*dx[0]))) * ((ETA(i+2,j,k,m)) - 4*(ETA(i+1,j,k,m))+6*(ETA(i,j,k,m)) - 4*(ETA(i-1,j,k,m)) + (ETA(i-2,j,k,m)));

		amrex::Real grad2222=(1/((dx[1]*dx[1]*dx[1]*dx[1]))) * ((ETA(i,j+2,k,m)) - 4*(ETA(i,j+1,k,m))+6*(ETA(i,j,k,m)) - 4*(ETA(i,j-1,k,m)) + (ETA(i,j-2,k,m)));
	       
		amrex::Real grad1112= (1/(24*dx[0]*dx[0]*dx[0]*dx[1]))*
		  ((-ETA(i+2,j+2,k,m)+8*ETA(i+2,j+1,k,m)
		    -8*ETA(i+2,j-1,k,m)+ETA(i+2,j-2,k,m))
		   -2*(-ETA(i+1,j+2,k,m)+8*ETA(i+1,j+1,k,m)
		       -8*ETA(i+1,j-1,k,m)+ETA(i+1,j-2,k,m))
		   +2*(-ETA(i-1,j+2,k,m)+8*ETA(i-1,j+1,k,m)
		       -8*ETA(i-1,j-1,k,m)+ETA(i-1,j-2,k,m))
		   -(-ETA(i-2,j+2,k,m)+8*ETA(i-2,j+1,k,m)
		     -8*ETA(i-2,j-1,k,m)+ETA(i-2,j-2,k,m)));

		amrex::Real grad1222= (1/(24*dx[0]*dx[1]*dx[1]*dx[1]))*
		  ((-ETA(i+2,j+2,k,m)+8*ETA(i+1,j+2,k,m)
		    -8*ETA(i-1,j+2,k,m)+ETA(i-2,j+2,k,m))
		   -2*(-ETA(i+2,j+1,k,m)+8*ETA(i+1,j+1,k,m)
		       -8*ETA(i-1,j+1,k,m)+ETA(i-2,j+1,k,m))
		   +2*(-ETA(i+2,j-1,k,m)+8*ETA(i+1,j-1,k,m)
		       -8*ETA(i-1,j-1,k,m)+ETA(i-2,j-1,k,m))
		   -(-ETA(i+2,j-2,k,m)+8*ETA(i+1,j-2,k,m)
		     -8*ETA(i-1,j-2,k,m)+ETA(i-2,j-2,k,m)));
		
		amrex::Real grad1122= (1/(144*dx[0]*dx[0]*dx[1]*dx[1]))*
		  (-(-ETA(i+2,j+2,k,m)+16*ETA(i+1,j+2,k,m)-30*ETA(i,j+2,k,m)
		     +16*ETA(i-1,j+2,k,m)-ETA(i-2,j+2,k,m))
		   +16*(-ETA(i+2,j+1,k,m)+16*ETA(i+1,j+1,k,m)-30*ETA(i,j+1,k,m)
			+16*ETA(i-1,j+1,k,m)-ETA(i-2,j+1,k,m))
		   -30*(-ETA(i+2,j,k,m)+16*ETA(i+1,j,k,m)-30*ETA(i,j,k,m)
			+16*ETA(i-1,j,k,m)-ETA(i-2,j,k,m))
		   +16*(-ETA(i+2,j-1,k,m)+16*ETA(i+1,j-1,k,m)-30*ETA(i,j-1,k,m)
			+16*ETA(i-1,j-1,k,m)-ETA(i-2,j-1,k,m))
		   -(-ETA(i+2,j-2,k,m)+16*ETA(i+1,j-2,k,m)-30*ETA(i,j-2,k,m)
		     +16*ETA(i-1,j-2,k,m)-ETA(i-2,j-2,k,m)));

		amrex::Real laplacian = grad11 + grad22;
	       
		amrex::Real kappa = l_gb*0.75*sigma0;
		mu = 0.75 * (1.0/0.23) * sigma0 / l_gb;

		if (anisotropy)
		  {
		    amrex::Real Theta = atan2(grad2,grad1);
		    amrex::Real Kappa = l_gb*0.75*boundary->W(Theta);
		    amrex::Real DKappa = l_gb*0.75*boundary->DW(Theta);
		    amrex::Real DDKappa = l_gb*0.75*boundary->DDW(Theta);
		    amrex::Real Mu = 0.75 * (1.0/0.23) * boundary->W(Theta) / l_gb;


		    eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),m) =
		      ETA(i,j,k,m) -
		      M*dt*(Mu*(ETA(i,j,k,m)*ETA(i,j,k,m)
				- 1.0 +
				2.0*gamma*sum_of_squares)*ETA(i,j,k,m)
			    - (Kappa*laplacian
			       + DKappa*(cos(2.0*Theta)*grad12 + 0.5*sin(2.0*Theta)*(grad22-grad11))
			       + damp*0.5*DDKappa*(sin(Theta)*sin(Theta)*grad11 - 2.*sin(Theta)*cos(Theta)*grad12 + cos(Theta)*cos(Theta)*grad22))+
			    beta*(grad1111*(sin(Theta)*sin(Theta)*sin(Theta)*sin(Theta))
				  +grad1112*(-6*sin(Theta)*sin(Theta)*sin(Theta)*cos(Theta))
				  +grad1122*(10*sin(Theta)*sin(Theta)*cos(Theta)*cos(Theta))
				  +grad1222*(-6*sin(Theta)*cos(Theta)*cos(Theta)*cos(Theta))
				  +grad2222*(cos(Theta)*cos(Theta)*cos(Theta)*cos(Theta))) 
			    );
		  }
		// Isotropic response if less than anisotropy_tstart
		else 
		  {
		    eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),m) =
		      ETA(i,j,k,m) -
		      M*dt*(mu*(ETA(i,j,k,m)*ETA(i,j,k,m)
				- 1.0 +
				2.0*gamma*sum_of_squares)*ETA(i,j,k,m)
			    - kappa*laplacian);
		  }

	      }

	  }
    }
}



void
PhaseFieldMicrostructure::Initialize (int lev)
{
  ic->Initialize(lev,eta_new);
  ic->Initialize(lev,eta_old);
}


void
PhaseFieldMicrostructure::TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{
  const amrex::Real* dx      = geom[lev].CellSize();

  amrex::Array<int>  itags;
	
  for (amrex::MFIter mfi(*eta_new[lev],true); mfi.isValid(); ++mfi)
    {
      const amrex::Box&  bx  = mfi.tilebox();

      amrex::TagBox&     tag  = tags[mfi];
	    
      amrex::BaseFab<Real> &eta_new_box = (*eta_new[lev])[mfi];

#if BL_SPACEDIM==2
      for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
	for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
	  {
	    for (int n = 0; n < number_of_grains; n++)
	      {
		amrex::Real gradx = (eta_new_box(amrex::IntVect(i+1,j),n) - eta_new_box(amrex::IntVect(i-1,j),n))/(2.*dx[0]);
		amrex::Real grady = (eta_new_box(amrex::IntVect(i,j+1),n) - eta_new_box(amrex::IntVect(i,j-1),n))/(2.*dx[1]);
		if (dx[0]*sqrt(gradx*gradx + grady*grady)>0.1) tag(amrex::IntVect(i,j)) = amrex::TagBox::SET;
	      }
	  }
#elif BL_SPACEDIM==3
      for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
	for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
	  for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
	    {
	      for (int n = 0; n < number_of_grains; n++)
		{
		  amrex::Real gradx = (eta_new_box(amrex::IntVect(i+1,j,k),n) - eta_new_box(amrex::IntVect(i-1,j,k),n))/(2.*dx[0]);
		  amrex::Real grady = (eta_new_box(amrex::IntVect(i,j+1,k),n) - eta_new_box(amrex::IntVect(i,j-1,k),n))/(2.*dx[1]);
		  amrex::Real gradz = (eta_new_box(amrex::IntVect(i,j,k+1),n) - eta_new_box(amrex::IntVect(i,j,k-1),n))/(2.*dx[2]);
		  if (dx[0]*sqrt(gradx*gradx + grady*grady + gradz*gradz)>0.1) tag(amrex::IntVect(i,j,k)) = amrex::TagBox::SET;
		}
	    }
#endif

    }

}

void PhaseFieldMicrostructure::TimeStepComplete(amrex::Real time)
{
  std::cout << "Starting Linear Solve now" << std::endl;
  std::cout << __FILE__ << ":" << __LINE__ << std::endl;

  // Hard code BCs for now.
  LPInfo info;
  info.setAgglomeration(true);
  info.setConsolidation(true);
  //const Real tol_rel = 1.e-10;
  //const Real tol_rel = 1e-10;//1.e-10;
  //const Real tol_abs = 0.0;
  MLStiffnessMatrix linop(geom,grids,dmap,info);

  linop.setMaxOrder(2);
  
  linop.setDomainBC(
		    {AMREX_D_DECL(LinOpBCType::Periodic,
				  LinOpBCType::Dirichlet,
				  LinOpBCType::Dirichlet)},
		    {AMREX_D_DECL(LinOpBCType::Periodic,
				  LinOpBCType::Dirichlet,
				  LinOpBCType::Dirichlet)});

  /// \todo This is not really necessary, should be handled with BC object.
  for (int ilev = 0; ilev < displacement.size(); ++ilev)
    {
      amrex::Box domain(geom[ilev].Domain());
      for (amrex::MFIter mfi(*displacement[ilev], true); mfi.isValid(); ++mfi)
	{
	  const amrex::Box& box = mfi.tilebox();
	  amrex::BaseFab<amrex::Real> &disp_box = (*displacement[ilev])[mfi];

	  for (int i = box.loVect()[0] - displacement[ilev]->nGrow(); i<=box.hiVect()[0] + displacement[ilev]->nGrow(); i++)
	    for (int j = box.loVect()[1] - displacement[ilev]->nGrow(); j<=box.hiVect()[1] + displacement[ilev]->nGrow(); j++)
	      { 
		if (j > domain.hiVect()[1]) // Top boundary
		  {
		    disp_box(amrex::IntVect(i,j),0) = 0.1;
		    disp_box(amrex::IntVect(i,j),1) = 0.0;
		  }
		else if (i > domain.hiVect()[0]) // Right boundary
		  {
		    disp_box(amrex::IntVect(i,j),0) = 0.0;
		    disp_box(amrex::IntVect(i,j),1) = 0.0;
		  }
		else if (j < domain.loVect()[1]) // Bottom
		  {
		    disp_box(amrex::IntVect(i,j),0) = 0.0;
		    disp_box(amrex::IntVect(i,j),1) = 0.0;
		  }
		else if (i < domain.loVect()[0]) // Left boundary
		  {
		    disp_box(amrex::IntVect(i,j),0) = 0.0;
		    disp_box(amrex::IntVect(i,j),1) = 0.0;
		  }
	      }
	}
      linop.setLevelBC(ilev,displacement[ilev].get());

      /// \todo Replace with proper driving force initialization
      body_force[ilev]->setVal(0.0);
    }
  
  amrex::Real tol_rel = 0, tol_abs = 1E-10;
  amrex::MLMG solver(linop);
  solver.setMaxIter(20);
  solver.setMaxFmgIter(0);
  solver.setVerbose(1);
  solver.solve(GetVecOfPtrs(displacement),
   	       GetVecOfConstPtrs(body_force),
   	       tol_rel,
   	       tol_abs);
  
  std::cout << "Done with Linear Solve" << std::endl;
}


#endif

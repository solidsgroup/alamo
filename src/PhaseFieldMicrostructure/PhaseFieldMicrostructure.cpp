#include "PhaseFieldMicrostructure/PhaseFieldMicrostructure.H"

#if BL_SPACEDIM == 2

PhaseFieldMicrostructure::PhaseFieldMicrostructure() :
  GeneralAMRIntegrator(), 
  mybc(geom)
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
  //  RegisterNewFab(Curvature, mybc, number_of_grains, 0, "Curvature");
  //  RegisterNewFab(Boundary_energy, mybc, number_of_grains, 0, "Boundary Energy");
}

// lev = max_level+1;
#define ETA(i,j,k,n) eta_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n)

void	
PhaseFieldMicrostructure::Advance ( int lev , Real time, Real dt)
{
  if (lev != max_level)
    {
      //      std::cout << (lev)<< " ";
      //std::cout << (max_level) << std::endl;
      
      return;
    }
  std::swap(eta_old[lev], eta_new[lev]);
  const Real* dx = geom[lev].CellSize();
  
  for ( amrex::MFIter mfi(*eta_new[lev],true); mfi.isValid(); ++mfi )
    {
      const amrex::Box& bx = mfi.tilebox();
      
      amrex::BaseFab<Real> &eta_new_box = (*eta_new[lev])[mfi];
      amrex::BaseFab<Real> &eta_old_box = (*eta_old[lev])[mfi];
	  // amrex::BaseFab<Real> &Curvature_box = (*Curvature[lev])[mfi];
	  // amrex::BaseFab<Real> &Boundary_energy_box = (*Boundary_energy[lev])[mfi];
      
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
		      
		  amrex::Real laplacian = grad11 + grad22;
		      
		  amrex::Real kappa = l_gb*0.75*sigma0;
		  mu = 0.75 * (1.0/0.23) * sigma0 / l_gb;
		  amrex::Real Theta = atan2(grad2,grad1);
		  amrex::Real Kappa = l_gb*0.75*boundary->W(Theta);
		  amrex::Real DKappa = l_gb*0.75*boundary->DW(Theta);
		  amrex::Real DDKappa = l_gb*0.75*boundary->DDW(Theta);
		  amrex::Real Mu = 0.75 * (1.0/0.23) * boundary->W(Theta) / l_gb;
		  amrex::Real Sine_theta = sin(Theta);
		  amrex::Real Cos_theta = cos(Theta);
	       
		
		  amrex::Real norm_grad = grad1*grad1+grad2*grad2;
		
		  //	Curvature_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),m) = Curvature_term;
		  //	Boundary_energy_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),m) = W - Boundary_term + beta*Curvature_term; 	     
	        
		    
		  if (time > anisotropy_tstart)
		    {
		      amrex::Real grad1111=(1/(dx[0]*dx[0]*dx[0]*dx[0])) * ((ETA(i+2,j,k,m)) - 4*(ETA(i+1,j,k,m))+6*(ETA(i,j,k,m)) - 4*(ETA(i-1,j,k,m)) + (ETA(i-2,j,k,m)));
			
		      amrex::Real grad2222=(1/(dx[1]*dx[1]*dx[1]*dx[1])) * ((ETA(i,j+2,k,m)) - 4*(ETA(i,j+1,k,m))+6*(ETA(i,j,k,m)) - 4*(ETA(i,j-1,k,m)) + (ETA(i,j-2,k,m)));
			
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
			
		      amrex::Real Curvature_term =
			grad1111*(Sine_theta*Sine_theta*Sine_theta*Sine_theta)
			+grad1112*(-4*Sine_theta*Sine_theta*Sine_theta*Cos_theta)
			+grad1122*(6*Sine_theta*Sine_theta*Cos_theta*Cos_theta)
			+grad1222*(-4*Sine_theta*Cos_theta*Cos_theta*Cos_theta)
			+grad2222*(Cos_theta*Cos_theta*Cos_theta*Cos_theta);
		      amrex::Real W =
			Mu*(ETA(i,j,k,m)*ETA(i,j,k,m)
			    - 1.0 + 2.0*gamma*sum_of_squares)*ETA(i,j,k,m);
		      amrex::Real Boundary_term =
			Kappa*laplacian +
			DKappa*(cos(2.0*Theta)*grad12 + 0.5*sin(2.0*Theta)*(grad22-grad11))
			+ damp*0.5*DDKappa*(Sine_theta*Sine_theta*grad11 - 2.*Sine_theta*Cos_theta*grad12 + Cos_theta*Cos_theta*grad22);
			
			
		      eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),m) =
			ETA(i,j,k,m) -
			M*dt*(W
			      - (Boundary_term)+
			      beta*(Curvature_term) 
			      );
			
		    }
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
		  
	      // Isotropic response if less than anisotropy_tstart
		  
		  
		  
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
		if (dx[0]*sqrt(gradx*gradx + grady*grady)>0.005) tag(amrex::IntVect(i,j)) = amrex::TagBox::SET;
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

#endif

#include "PhaseFieldMicrostructure.H"

#if BL_SPACEDIM == 2

PhaseFieldMicrostructure::PhaseFieldMicrostructure() :
  GeneralAMRIntegrator(), 
  mybc(geom)
{

  //
  // READ INPUT PARAMETERS
  //
  { 
    amrex::ParmParse pp;   // Basic run parameters
    pp.query("newOld",newOld);
  }

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
    //std::string filename;
    
    amrex::ParmParse pp("anisotropy"); // Phase-field model parameters
    pp.query("on", anisotropy);
    pp.query("theta0", theta0);
    //pp.query("filename", filename);
    theta0 *= 0.01745329251; // convert degrees into radians
    pp.query("sigma0", sigma0);
    pp.query("sigma1", sigma1);
    pp.query("beta", beta);
    pp.query("damp", damp);
    pp.query("tstart", anisotropy_tstart);

    boundary = new PFBoundarySin(theta0,sigma0,sigma1);    

    // if(filename == "no")
    //   {
    // 	if(newOld=="new")
    // 	  {  
    // 	    std::cout << "\n**USING NEW VERSION**\n" << std::endl;
    // 	    boundary = new PFBoundaryAbsSin(theta0,sigma0,sigma1);
    // 	  }
    // 	else if(newOld=="old")
    // 	  {
    // 	    std::cout << "\n**USING OLD VERSION**\n" << std::endl;
    // 	    boundary = new PFBoundarySin(theta0,sigma0,sigma1);
    // 	  }
    // 	else
    // 	  {
    // 	    std::cout << "\n**NO VERSION SELECTED**\n" << std::endl;
    // 	  }
    //   }
    // else
    //   {
    // 	std::cout << "\n**READING FILE "+filename+"**\n" << std::endl;
    // 	boundary = new PFBoundaryRead(filename);
    //   }

    // bool a = boundary -> Test();
    // if (a==false)
    //   {
    //  	std::cout << "\n**BOUNDARY DERIVATIVE TEST DID NOT PASS**\n" << std::endl;
    //  	//amrex::Abort("Boundary derivative test did not pass");
    //   }
    // else
    //   {
    //  	std::cout <<"\n**TEST SUCCESSFULLY PASSED**\n" << std::endl;
    //   }

    
    // if(ParallelDescriptor::IOProcessor())
    //   if (!boundary->Test()) amrex::Error("Boundary functor does not pass derivative test");
  }

  {
    amrex::ParmParse pp("ic"); // Phase-field model parameters
    pp.query("type", ic_type);
  }

  RegisterNewFab(eta_new, mybc, number_of_grains, number_of_ghost_cells, "Eta");
  RegisterNewFab(eta_old, mybc, number_of_grains, number_of_ghost_cells, "Eta old");
}


void
PhaseFieldMicrostructure::Advance (int lev, Real time, Real dt)
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
	  {

	    for (int m = 0; m < number_of_grains; m++)
	      {

		amrex::Real sum_of_squares = 0.;
		for (int n = 0; n < number_of_grains; n++)
		  {
		    if (n==m) continue;
		    sum_of_squares += eta_old_box(amrex::IntVect(i,j),n)*eta_old_box(amrex::IntVect(i,j),n);
		  }
	       
	       
		amrex::Real grad1_normal = (eta_old_box(amrex::IntVect(i+1,j),m) - eta_old_box(amrex::IntVect(i-1,j),m))/(2*dx[0]);
	       
		amrex::Real grad2_normal = (eta_old_box(amrex::IntVect(i,j+1),m) - eta_old_box(amrex::IntVect(i,j-1),m))/(2*dx[1]);

		amrex::Real grad1 =  grad1_normal;

		amrex::Real grad2 =  grad2_normal;

	       
		amrex::Real grad12_normal = (eta_old_box(amrex::IntVect(i+1,j+1),m) - eta_old_box(amrex::IntVect(i-1,j+1),m) - eta_old_box(amrex::IntVect(i+1,j-1),m) + eta_old_box(amrex::IntVect(i-1,j-1),m))/(4.*dx[0]*dx[1]);
	       
	       
		amrex::Real grad11_normal =  (eta_old_box(amrex::IntVect(i+1,j),m) - 2.*eta_old_box(amrex::IntVect(i,j),m) + eta_old_box(amrex::IntVect(i-1,j),m))/dx[0]/dx[0]; // 3 point
	       
		amrex::Real grad22_normal = (eta_old_box(amrex::IntVect(i,j+1),m) - 2.*eta_old_box(amrex::IntVect(i,j),m) + eta_old_box(amrex::IntVect(i,j-1),m))/dx[1]/dx[1]; // 3 point
		//Doing the change of basis and averaging them.
		amrex::Real grad11 = grad11_normal;
		amrex::Real grad22 = grad22_normal;
		amrex::Real grad12 = grad12_normal;
	       
		//Second curvature variational terms
		amrex::Real grad1111=(1/((dx[0]*dx[0]*dx[0]*dx[0])))*
		  ((eta_old_box(amrex::IntVect(i+2,j),m))-4*(eta_old_box(amrex::IntVect(i+1,j),m))+6*(eta_old_box(amrex::IntVect(i,j),m))
		   -4*(eta_old_box(amrex::IntVect(i-1,j),m))+(eta_old_box(amrex::IntVect(i-2,j),m)));

		amrex::Real grad2222=(1/((dx[1]*dx[1]*dx[1]*dx[1])))*
		  ((eta_old_box(amrex::IntVect(i,j+2),m))-4*(eta_old_box(amrex::IntVect(i,j+1),m))+6*(eta_old_box(amrex::IntVect(i,j),m))
		   -4*(eta_old_box(amrex::IntVect(i,j-1),m))+(eta_old_box(amrex::IntVect(i,j-2),m)));
	       
		amrex::Real grad1112= (1/(24*dx[0]*dx[0]*dx[0]*dx[1]))*
		  ((-eta_old_box(amrex::IntVect(i+2,j+2),m)+8*eta_old_box(amrex::IntVect(i+2,j+1),m)
		    -8*eta_old_box(amrex::IntVect(i+2,j-1),m)+eta_old_box(amrex::IntVect(i+2,j-2),m))
		   -2*(-eta_old_box(amrex::IntVect(i+1,j+2),m)+8*eta_old_box(amrex::IntVect(i+1,j+1),m)
		       -8*eta_old_box(amrex::IntVect(i+1,j-1),m)+eta_old_box(amrex::IntVect(i+1,j-2),m))
		   +2*(-eta_old_box(amrex::IntVect(i-1,j+2),m)+8*eta_old_box(amrex::IntVect(i-1,j+1),m)
		       -8*eta_old_box(amrex::IntVect(i-1,j-1),m)+eta_old_box(amrex::IntVect(i-1,j-2),m))
		   -(-eta_old_box(amrex::IntVect(i-2,j+2),m)+8*eta_old_box(amrex::IntVect(i-2,j+1),m)
		     -8*eta_old_box(amrex::IntVect(i-2,j-1),m)+eta_old_box(amrex::IntVect(i-2,j-2),m)));

		amrex::Real grad1222= (1/(24*dx[0]*dx[1]*dx[1]*dx[1]))*
		  ((-eta_old_box(amrex::IntVect(i+2,j+2),m)+8*eta_old_box(amrex::IntVect(i+1,j+2),m)
		    -8*eta_old_box(amrex::IntVect(i-1,j+2),m)+eta_old_box(amrex::IntVect(i-2,j+2),m))
		   -2*(-eta_old_box(amrex::IntVect(i+2,j+1),m)+8*eta_old_box(amrex::IntVect(i+1,j+1),m)
		       -8*eta_old_box(amrex::IntVect(i-1,j+1),m)+eta_old_box(amrex::IntVect(i-2,j+1),m))
		   +2*(-eta_old_box(amrex::IntVect(i+2,j-1),m)+8*eta_old_box(amrex::IntVect(i+1,j-1),m)
		       -8*eta_old_box(amrex::IntVect(i-1,j-1),m)+eta_old_box(amrex::IntVect(i-2,j-1),m))
		   -(-eta_old_box(amrex::IntVect(i+2,j-2),m)+8*eta_old_box(amrex::IntVect(i+1,j-2),m)
		     -8*eta_old_box(amrex::IntVect(i-1,j-2),m)+eta_old_box(amrex::IntVect(i-2,j-2),m)));
		
		amrex::Real grad1122= (1/(144*dx[0]*dx[0]*dx[1]*dx[1]))*
		  (-(-eta_old_box(amrex::IntVect(i+2,j+2),m)+16*eta_old_box(amrex::IntVect(i+1,j+2),m)-30*eta_old_box(amrex::IntVect(i,j+2),m)
		     +16*eta_old_box(amrex::IntVect(i-1,j+2),m)-eta_old_box(amrex::IntVect(i-2,j+2),m))
		   +16*(-eta_old_box(amrex::IntVect(i+2,j+1),m)+16*eta_old_box(amrex::IntVect(i+1,j+1),m)-30*eta_old_box(amrex::IntVect(i,j+1),m)
			+16*eta_old_box(amrex::IntVect(i-1,j+1),m)-eta_old_box(amrex::IntVect(i-2,j+1),m))
		   -30*(-eta_old_box(amrex::IntVect(i+2,j),m)+16*eta_old_box(amrex::IntVect(i+1,j),m)-30*eta_old_box(amrex::IntVect(i,j),m)
			+16*eta_old_box(amrex::IntVect(i-1,j),m)-eta_old_box(amrex::IntVect(i-2,j),m))
		   +16*(-eta_old_box(amrex::IntVect(i+2,j-1),m)+16*eta_old_box(amrex::IntVect(i+1,j-1),m)-30*eta_old_box(amrex::IntVect(i,j-1),m)
			+16*eta_old_box(amrex::IntVect(i-1,j-1),m)-eta_old_box(amrex::IntVect(i-2,j-1),m))
		   -(-eta_old_box(amrex::IntVect(i+2,j-2),m)+16*eta_old_box(amrex::IntVect(i+1,j-2),m)-30*eta_old_box(amrex::IntVect(i,j-2),m)
		     +16*eta_old_box(amrex::IntVect(i-1,j-2),m)-eta_old_box(amrex::IntVect(i-2,j-2),m)));

	       
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
		    //amrex::Real beta= 0.00002;
		    amrex::Real beta2 = 1*beta;
		    amrex::Real timebeta = 30;
		    //amrex::Real Mu = mu;
	

		    if (time>anisotropy_tstart)
		      {
			eta_new_box(amrex::IntVect(i,j),m) =
			  eta_old_box(amrex::IntVect(i,j),m) -
			  M*dt*(Mu*(eta_old_box(amrex::IntVect(i,j),m)*eta_old_box(amrex::IntVect(i,j),m)
				    - 1.0 +
				    2.0*gamma*sum_of_squares)*eta_old_box(amrex::IntVect(i,j),m)
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
		    else if (time>timebeta)// beta change
		      {
			eta_new_box(amrex::IntVect(i,j),m) =
			  eta_old_box(amrex::IntVect(i,j),m) -
			  M*dt*(Mu*(eta_old_box(amrex::IntVect(i,j),m)*eta_old_box(amrex::IntVect(i,j),m)
				    - 1.0 +
				    2.0*gamma*sum_of_squares)*eta_old_box(amrex::IntVect(i,j),m)
				- (Kappa*laplacian
				   + DKappa*(cos(2.0*Theta)*grad12 + 0.5*sin(2.0*Theta)*(grad22-grad11))
				   + damp*0.5*DDKappa*(sin(Theta)*sin(Theta)*grad11 - 2.*sin(Theta)*cos(Theta)*grad12 + cos(Theta)*cos(Theta)*grad22))+
				beta2*(grad1111*(sin(Theta)*sin(Theta)*sin(Theta)*sin(Theta))
				       +grad1112*(-6*sin(Theta)*sin(Theta)*sin(Theta)*cos(Theta))
				       +grad1122*(10*sin(Theta)*sin(Theta)*cos(Theta)*cos(Theta))
				       +grad1222*(-6*sin(Theta)*cos(Theta)*cos(Theta)*cos(Theta))
				       +grad2222*(cos(Theta)*cos(Theta)*cos(Theta)*cos(Theta))) 
				);
		      }
		    else
		      {
			eta_new_box(amrex::IntVect(i,j),m) =
			  eta_old_box(amrex::IntVect(i,j),m) -
			  M*dt*(Mu*(eta_old_box(amrex::IntVect(i,j),m)*eta_old_box(amrex::IntVect(i,j),m)
				    - 1.0 +
				    2.0*gamma*sum_of_squares)*eta_old_box(amrex::IntVect(i,j),m)
				- kappa*laplacian);
		      }


		  }
	      }

	  }
    }
}



void
PhaseFieldMicrostructure::Initialize (int lev)
{
  const amrex::Real width = geom[lev].ProbHi()[0] - geom[lev].ProbHi()[1];

  for (amrex::MFIter mfi(*eta_new[lev],true); mfi.isValid(); ++mfi)
    {
      const amrex::Box& box = mfi.tilebox();

      amrex::BaseFab<Real> &eta_new_box = (*eta_new[lev])[mfi];

      for (int i = box.loVect()[0]-number_of_ghost_cells; i<=box.hiVect()[0]+number_of_ghost_cells; i++) 
       	for (int j = box.loVect()[1]-number_of_ghost_cells; j<=box.hiVect()[1]+number_of_ghost_cells; j++)
#if BL_SPACEDIM==3
	  for (int k = box.loVect()[2]-number_of_ghost_cells; k<=box.hiVect()[2]+number_of_ghost_cells; k++)
#endif
	    {
	      amrex::Real x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];
	      amrex::Real y = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];
#if BL_SPACEDIM==3
	      amrex::Real z = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];
#endif

	      if (anisotropy)
		{
		  if (ic_type == "circle")
		    {
		      //
		      // circle IC
		      //
		      if (sqrt(x*x+y*y)<=0.5)

			{
			  eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = 1.;     
			  eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),1) = 0.;     
			}
		      else
			{
			  eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = 0.;     
			  eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),1) = 1.;     
			}
		    }
		  else if (ic_type == "perturbed_bar")
		    {
		      //
		      // perturbed bar IC
		      //
		      amrex::Real pi = 3.14159265359;
		      amrex::Real bdry  =
			(+ 0.1 * sin(x*(2*pi)/width)
			 + 0.1 * cos(2.0*x*(2*pi)/width)
			 + 0.1 * sin(3.0*x*(2*pi)/width)
			 + 0.1 * cos(4.0*x*(2*pi)/width)
			 + 0.1 * sin(5.0*x*(2*pi)/width)
			 + 0.15 * cos(6.0*x*(2*pi)/width)
			 + 0.15 * sin(7.0*x*(2*pi)/width)
			 + 0.15 * cos(8.0*x*(2*pi)/width)
			 + 0.15 * sin(9.0*x*(2*pi)/width)
			 + 0.15 * cos(10.0*x*(2*pi)/width))*0.1
			;

		      if (y < bdry)

			{
			  eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = 1.;     
			  eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),1) = 0.;     
			}
		      else
			{
			  eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = 0.;     
			  eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),1) = 1.;     
			}
		    }
		  else if (ic_type == "random")
		    {
		      eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = ((amrex::Real)rand())/((amrex::Real)RAND_MAX);
		      eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),1) = 1.0 - eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),0);
		    }
		}	      
	    }
    }
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

#endif

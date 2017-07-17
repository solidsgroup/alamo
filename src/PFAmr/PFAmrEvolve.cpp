
#include <limits>

#include <AMReX_MultiFabUtil.H>
#include <PFAmr.H>

using namespace amrex;

void
PFAmr::Evolve ()
{
  Real cur_time = t_new[0];
  int last_plot_file_step = 0;

  for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "\nSTEP " << step+1 << " starts ..." << std::endl;
      }
      int lev = 0;
      int iteration = 1;
      TimeStep(lev, cur_time, iteration);
      cur_time += dt[0];

      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "STEP " << step+1 << " ends."
		  << " TIME = " << cur_time << " DT = " << dt[0]
		  << std::endl;
      }

      // sync up time
      for (int lev = 0; lev <= finest_level; ++lev) {
	t_new[lev] = cur_time;
      }

      if (plot_int > 0 && (step+1) % plot_int == 0) {
	last_plot_file_step = step+1;
	WritePlotFile();
      }

      if (cur_time >= stop_time - 1.e-6*dt[0]) break;
    }

  if (plot_int > 0 && istep[0] > last_plot_file_step) {
    WritePlotFile();
  }
}

void
PFAmr::TimeStep (int lev, Real time, int iteration)
{

  if (regrid_int > 0)  // We may need to regrid
    {
      static Array<int> last_regrid_step(max_level+1, 0);

      // regrid doesn't change the base level, so we don't regrid on max_level
      if (lev < max_level && istep[lev] > last_regrid_step[lev])
	{
          if (istep[lev] % regrid_int == 0)
	    {
	      int old_finest = finest_level; // regrid changes finest_level
	      regrid(lev, time, false); // <--segfault is here.
	      for (int k = lev; k <= finest_level; ++k) {
		last_regrid_step[k] = istep[k];
	      }
	      for (int k = old_finest+1; k <= finest_level; ++k) {
		dt[k] = dt[k-1] / MaxRefRatio(k-1);
	      }
  	    }
  	}
    }

  if (Verbose() && ParallelDescriptor::IOProcessor()) {
    std::cout << "[Level " << lev 
	      << " step " << istep[lev]+1 << "] ";
    std::cout << "ADVANCE with dt = "
	      << dt[lev]
	      << std::endl;
  }
  Advance(lev, time, dt[lev]);
  ++istep[lev];

  if (Verbose() && ParallelDescriptor::IOProcessor())
    {
      std::cout << "[Level " << lev
		<< " step " << istep[lev] << "] ";
      std::cout << "Advanced "
		<< CountCells(lev)
		<< " cells"
		<< std::endl;
    }

  if (lev < finest_level)
    {
      // for (int i = 1; i <= nsubsteps[lev+1]; ++i)
      // {
      //    TimeStep(lev+1, time+(i-1)*dt[lev+1], i);
      // }
      TimeStep(lev+1, time, 1);
      for (int n = 0; n < number_of_fabs; n++)
	{
	  amrex::average_down(*phi_new[n][lev+1], *phi_new[n][lev],
			      geom[lev+1], geom[lev],
			      0, phi_new[n][lev]->nComp(), refRatio(lev));
	}
    }
}

void
PFAmr::Advance (int lev, Real time, Real dt)
{
 std::swap(phi_old[0][lev], phi_new[0][lev]);
 const Real* dx = geom[lev].CellSize();

 amrex::Array<std::unique_ptr<amrex::MultiFab> > Sborder(number_of_fabs);
 for (int n=0; n<number_of_fabs; n++)
   Sborder[n].reset(new amrex::MultiFab(grids[lev], dmap[lev], number_of_grains+2,nghost)); 
 FillPatch(lev,t_old[lev],Sborder,0);

 for ( MFIter mfi(*phi_new[0][lev],true); mfi.isValid(); ++mfi )
   {
     const Box& bx = mfi.tilebox();

     amrex::BaseFab<Real> &old_phi = (*Sborder[0])[mfi];
     amrex::BaseFab<Real> &new_phi = (*phi_new[0][lev])[mfi];

#if BL_SPACEDIM == 2
     amrex::Array<amrex::Real> Laplacian(number_of_fabs);
     for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
       for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
	 {
	   new_phi(amrex::IntVect(i,j),number_of_grains) = 0.;
	   new_phi(amrex::IntVect(i,j),number_of_grains+1) = 0.;

	   for (int m = 0; m < number_of_grains; m++)
	     {

	       amrex::Real sum_of_squares = 0.;
	       for (int n = 0; n < number_of_grains; n++)
		 {
		   if (n==m) continue;
		   sum_of_squares += old_phi(amrex::IntVect(i,j),n)*old_phi(amrex::IntVect(i,j),n);
		 }
	       
	       // THREE POINT STENCILS
	       //amrex::Real grad1  = (old_phi(amrex::IntVect(i+1,j),m) - old_phi(amrex::IntVect(i-1,j),m))/(2*dx[0]); // 3 point
	       //amrex::Real grad2  = (old_phi(amrex::IntVect(i,j+1),m) - old_phi(amrex::IntVect(i,j-1),m))/(2*dx[1]); // 3 point
	       //amrex::Real grad11 = (old_phi(amrex::IntVect(i+1,j),m) - 2.*old_phi(amrex::IntVect(i,j),m) + old_phi(amrex::IntVect(i-1,j),m))/dx[0]/dx[0]; // 3 point
	       // amrex::Real grad22 = (old_phi(amrex::IntVect(i,j+1),m) - 2.*old_phi(amrex::IntVect(i,j),m) + old_phi(amrex::IntVect(i,j-1),m))/dx[1]/dx[1]; // 3 point
	       // amrex::Real grad12 = (old_phi(amrex::IntVect(i+1,j+1),m) - old_phi(amrex::IntVect(i-1,j+1),m) - old_phi(amrex::IntVect(i+1,j-1),m) + old_phi(amrex::IntVect(i-1,j-1),m))/(4.*dx[0]*dx[1]);

	       // // FIVE POINT STENCILS
	       // amrex::Real grad1  = (-old_phi(amrex::IntVect(i+2,j),m)+8*old_phi(amrex::IntVect(i+1,j),m) -8*old_phi(amrex::IntVect(i-1,j),m)+old_phi(amrex::IntVect(i-2,j),m))/(12*dx[0]); // 5 point
	       // amrex::Real grad2  = (-old_phi(amrex::IntVect(i,j+2),m)+8*old_phi(amrex::IntVect(i,j+1),m) -8*old_phi(amrex::IntVect(i,j-1),m)+old_phi(amrex::IntVect(i,j-2),m))/(12*dx[0]); // 5 point
	       // amrex::Real grad11 = (-old_phi(amrex::IntVect(i+2,j),m)+16*old_phi(amrex::IntVect(i+1,j),m)-30*old_phi(amrex::IntVect(i,j),m)+16*old_phi(amrex::IntVect(i-1,j),m)-old_phi(amrex::IntVect(i-2,j),m))/(12*dx[0]*dx[0]); 
	       // amrex::Real grad22 = (-old_phi(amrex::IntVect(i,j+2),m)+16*old_phi(amrex::IntVect(i,j+1),m)-30*old_phi(amrex::IntVect(i,j),m)+16*old_phi(amrex::IntVect(i,j-1),m)-old_phi(amrex::IntVect(i,j-2),m))/(12*dx[1]*dx[1]); 
	       amrex::Real grad12 =
		 (+ 8*(old_phi(amrex::IntVect(i+1,j-2),m)+old_phi(amrex::IntVect(i+2,j-1),m)+old_phi(amrex::IntVect(i-2,j+1),m)+old_phi(amrex::IntVect(i-1,j+2),m))
		  - 8*(old_phi(amrex::IntVect(i-1,j-2),m)+old_phi(amrex::IntVect(i-2,j-1),m)+old_phi(amrex::IntVect(i+1,j+2),m)+old_phi(amrex::IntVect(i+2,j+1),m))
		  - 1*(old_phi(amrex::IntVect(i+2,j-2),m)+old_phi(amrex::IntVect(i-2,j+2),m)-old_phi(amrex::IntVect(i-2,j-2),m)-old_phi(amrex::IntVect(i+2,j+2),m))
		  +64*(old_phi(amrex::IntVect(i-1,j-1),m)+old_phi(amrex::IntVect(i+1,j+1),m)-old_phi(amrex::IntVect(i+1,j-1),m)-old_phi(amrex::IntVect(i-1,j+1),m)))
		 /(144*dx[0]*dx[1]);
	       
	       // Rotated stencils (5 Points)
	       
	       amrex::Real grad1  =
		 (0.5*(-old_phi(amrex::IntVect(i+2,j),m)+8*old_phi(amrex::IntVect(i+1,j),m) -8*old_phi(amrex::IntVect(i-1,j),m)+old_phi(amrex::IntVect(i-2,j),m))
		  +0.25*(-old_phi(amrex::IntVect(i+2,j+2),m)+8*old_phi(amrex::IntVect(i+1,j+1),m) -8*old_phi(amrex::IntVect(i-1,j-1),m)+old_phi(amrex::IntVect(i-2,j-2),m)+old_phi(amrex::IntVect(i-2,j+2),m)-8*old_phi(amrex::IntVect(i-1,j+1),m) +8*old_phi(amrex::IntVect(i+1,j-1),m)-old_phi(amrex::IntVect(i+2,j-2),m)))/(12*dx[0]); // 5 point
	       
	       amrex::Real grad2  =
		 (0.5*(-old_phi(amrex::IntVect(i,j+2),m)+8*old_phi(amrex::IntVect(i,j+1),m) -8*old_phi(amrex::IntVect(i,j-1),m)+old_phi(amrex::IntVect(i,j-2),m))
		  +0.25*(-old_phi(amrex::IntVect(i+2,j+2),m)+8*old_phi(amrex::IntVect(i+1,j+1),m) -8*old_phi(amrex::IntVect(i-1,j-1),m)+old_phi(amrex::IntVect(i-2,j-2),m)-old_phi(amrex::IntVect(i-2,j+2),m)+8*old_phi(amrex::IntVect(i-1,j+1),m) -8*old_phi(amrex::IntVect(i+1,j-1),m)+old_phi(amrex::IntVect(i+2,j-2),m)))/(12*dx[0]); // 5 point
	       
	       amrex::Real grad11 =
		 (0.5*(-old_phi(amrex::IntVect(i+2,j),m)+16*old_phi(amrex::IntVect(i+1,j),m)-30*old_phi(amrex::IntVect(i,j),m)+16*old_phi(amrex::IntVect(i-1,j),m)-old_phi(amrex::IntVect(i-2,j),m))+0.25*sqrt(2)*(-old_phi(amrex::IntVect(i+2,j+2),m)+16*old_phi(amrex::IntVect(i+1,j+1),m)+16*old_phi(amrex::IntVect(i-1,j-1),m)-old_phi(amrex::IntVect(i-2,j-2),m)+old_phi(amrex::IntVect(i-2,j+2),m)-16*old_phi(amrex::IntVect(i-1,j+1),m)-16*old_phi(amrex::IntVect(i+1,j-1),m)+old_phi(amrex::IntVect(i+2,j-2),m)))/(12*dx[0]*dx[0]);
	       
	       amrex::Real grad22 = (0.5*(-old_phi(amrex::IntVect(i,j+2),m)+16*old_phi(amrex::IntVect(i,j+1),m)-30*old_phi(amrex::IntVect(i,j),m)+16*old_phi(amrex::IntVect(i,j-1),m)-old_phi(amrex::IntVect(i,j-2),m))+0.25*sqrt(2)*(-old_phi(amrex::IntVect(i+2,j+2),m)+16*old_phi(amrex::IntVect(i+1,j+1),m)+16*old_phi(amrex::IntVect(i-1,j-1),m)-old_phi(amrex::IntVect(i-2,j-2),m)-old_phi(amrex::IntVect(i-2,j+2),m)+16*old_phi(amrex::IntVect(i-1,j+1),m)-60*old_phi(amrex::IntVect(i,j),m)+16*old_phi(amrex::IntVect(i+1,j-1),m)-old_phi(amrex::IntVect(i+2,j-2),m)))/(12*dx[1]*dx[1]); 
	       
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
		   //amrex::Real Mu = mu;
		   

		   // amrex::Real DDKappa1 = l_gb*0.75*(boundary->DW(Theta_E) - boundary->DW(Theta_W))/(2*dx[0]);
		   // amrex::Real DDKappa2 = l_gb*0.75*(boundary->DW(Theta_N) - boundary->DW(Theta_S))/(2*dx[0]);


		   if (time>tstart)
		     {
		       new_phi(amrex::IntVect(i,j),m) =
			 old_phi(amrex::IntVect(i,j),m) -
			 M*dt*(Mu*(old_phi(amrex::IntVect(i,j),m)*old_phi(amrex::IntVect(i,j),m)
				   - 1.0 +
				   2.0*gamma*sum_of_squares)*old_phi(amrex::IntVect(i,j),m)
			       - (Kappa*laplacian
				  + DKappa*(cos(2.0*Theta)*grad12 + 0.5*sin(2.0*Theta)*(grad22-grad11))
				  + damp*0.5*DDKappa*(sin(Theta)*sin(Theta)*grad11 - 2.*sin(Theta)*cos(Theta)*grad12 + cos(Theta)*cos(Theta)*grad22)
				  )
			       );
		     }
		   else
		     {
		       new_phi(amrex::IntVect(i,j),m) =
			 old_phi(amrex::IntVect(i,j),m) -
			 M*dt*(mu*(old_phi(amrex::IntVect(i,j),m)*old_phi(amrex::IntVect(i,j),m)
				   - 1.0
				   + 2.0*gamma*sum_of_squares)*old_phi(amrex::IntVect(i,j),m)
			       - kappa*laplacian);
		     }


		   new_phi(amrex::IntVect(i,j),number_of_grains) = DDKappa*(sin(Theta)*sin(Theta)*grad11 - 2.*sin(Theta)*cos(Theta)*grad12 + cos(Theta)*cos(Theta)*grad22);//*(cos(2.0*Theta)*grad12 + 0.5*sin(2.0*Theta)*(grad22-grad11));
		   new_phi(amrex::IntVect(i,j),number_of_grains+1) = sin(Theta)*sin(Theta)*grad11 - 2.*sin(Theta)*cos(Theta)*grad12 + cos(Theta)*cos(Theta)*grad22;//*(sin(Theta)*sin(Theta)*grad11 - 2.*sin(Theta)*cos(Theta)*grad12 + cos(Theta)*cos(Theta)*grad22);
		 }
	       else
		 {
		   new_phi(amrex::IntVect(i,j),m) =
		     old_phi(amrex::IntVect(i,j),m) -
		     M*dt*(mu*(old_phi(amrex::IntVect(i,j),m)*old_phi(amrex::IntVect(i,j),m)
			       - 1.0 +
			       2.0*gamma*sum_of_squares)*old_phi(amrex::IntVect(i,j),m)
			   - kappa*laplacian);
		 }


	       //if (new_phi(amrex::IntVect(i,j),m)>0.5) new_phi(amrex::IntVect(i,j),number_of_grains) = (amrex::Real)m;

	       // Boundary field

	       //new_phi(amrex::IntVect(i,j),number_of_grains+1) += sqrt(grad1*grad1 + grad2*grad2);

	     }

	 }
#endif
   }
}


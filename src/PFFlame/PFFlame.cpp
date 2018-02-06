#include "PFFlame.H"

PFFlame::PFFlame () : GeneralAMRIntegrator(), TempBC(geom,"TempBC"), EtaBC(geom,"EtaBC")
{
  amrex::ParmParse pp("physics"); 
  pp.query("M",M);
  pp.query("kappa",kappa);
  pp.query("w1",w1);
  pp.query("w12",w12);
  pp.query("w0",w0);
  pp.query("rho1",rho1);
  pp.query("rho0",rho0);
  pp.query("k1",k1);
  pp.query("k0",k0);
  pp.query("cp1",cp1);
  pp.query("cp0",cp0);
  pp.query("qdotburn",qdotburn);

  RegisterNewFab(Temp,     TempBC, 1, 1, "Temp");
  RegisterNewFab(Temp_old, TempBC, 1, 1, "Temp_old");
  RegisterNewFab(Eta,      EtaBC,  1, 1, "Eta");
  RegisterNewFab(Eta_old,  EtaBC,  1, 1, "Eta_old");

}

void PFFlame::Initialize (int lev)
{
  for (amrex::MFIter mfi(*Temp[lev],true); mfi.isValid(); ++mfi)
       {
	 const amrex::Box& box = mfi.tilebox();

	 amrex::BaseFab<Real> &Eta_box		= (*Eta[lev])[mfi];
	 amrex::BaseFab<Real> &Eta_old_box	= (*Eta_old[lev])[mfi];
	 amrex::BaseFab<Real> &Temp_box		= (*Temp[lev])[mfi];
	 amrex::BaseFab<Real> &Temp_old_box	= (*Temp_old[lev])[mfi];

	 for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++) 
	   for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++)
#if BL_SPACEDIM > 2
	   for (int k = box.loVect()[2]; k<=box.hiVect()[2]; k++)
#endif
	     {
	       Eta_box     (amrex::IntVect(AMREX_D_DECL(i,j,k))) =  1.0;
	       Eta_old_box (amrex::IntVect(AMREX_D_DECL(i,j,k))) =  1.0;
	       Temp_box    (amrex::IntVect(AMREX_D_DECL(i,j,k))) =  0.0;
	       Temp_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) =  0.0;
	     }
    }
}



void PFFlame::Advance (int lev, Real time, Real dt)
{
  std::swap(Eta_old [lev], Eta [lev]);
  std::swap(Temp_old[lev], Temp[lev]);

  const Real* dx = geom[lev].CellSize();

  amrex::Real a0=w0, a1=0.0, a2= -5*w1 + 16*w12 - 11*a0, a3=14*w1 - 32*w12 + 18*a0, a4=-8*w1 + 16*w12 - 8*a0;

  for ( amrex::MFIter mfi(*Temp[lev],true); mfi.isValid(); ++mfi )
    {
      const amrex::Box& bx = mfi.tilebox();

      amrex::BaseFab<Real> &Eta_box		= (*Eta[lev])[mfi];
      amrex::BaseFab<Real> &Eta_old_box		= (*Eta_old[lev])[mfi];
      amrex::BaseFab<Real> &Temp_box		= (*Temp[lev])[mfi];
      amrex::BaseFab<Real> &Temp_old_box	= (*Temp_old[lev])[mfi];


      for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
	for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
#if BL_SPACEDIM>2
	for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
#endif
	  {
	    //
	    // Phase field evolution
	    //

	    amrex::Real eta_lap = 
	      (Eta_old_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k))) - 2.*Eta_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) + Eta_old_box(amrex::IntVect(AMREX_D_DECL(i-1,j,k))))/dx[0]/dx[0] +
	      (Eta_old_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k))) - 2.*Eta_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) + Eta_old_box(amrex::IntVect(AMREX_D_DECL(i,j-1,k))))/dx[1]/dx[1];
	     
	     
	    amrex::Real oldeta = Eta_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k)));
	    Eta_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) = oldeta -
	      M * dt * (a1 + 2*a2*oldeta + 3*a3*oldeta*oldeta + 4*a4*oldeta*oldeta*oldeta
			- kappa*eta_lap);

	    //
	    // Temperature evolution
	    //
	    
	    //if (0) // todo ignore temperature for the time being
	      //if (Eta_box(amrex::IntVect(i,j))>0.01)
	    //{
	    amrex::Real temperature_delay = 0.05;
	    if (time<temperature_delay) continue;

	    amrex::Real eta_gradx = (Eta_old_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k))) - Eta_old_box(amrex::IntVect(AMREX_D_DECL(i-1,j,k))))/(2*dx[0]);
	    amrex::Real eta_grady = (Eta_old_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k))) - Eta_old_box(amrex::IntVect(AMREX_D_DECL(i,j-1,k))))/(2*dx[1]);
	    amrex::Real T_gradx = (Temp_old_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k))) - Temp_old_box(amrex::IntVect(AMREX_D_DECL(i-1,j,k))))/(2*dx[0]);
	    amrex::Real T_grady = (Temp_old_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k))) - Temp_old_box(amrex::IntVect(AMREX_D_DECL(i,j-1,k))))/(2*dx[1]);

	    amrex::Real eta_grad_mag = sqrt(eta_gradx*eta_gradx + eta_grady*eta_grady);

	    amrex::Real T_lap = 
	      (Temp_old_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k))) - 2.*Temp_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) + Temp_old_box(amrex::IntVect(AMREX_D_DECL(i-1,j,k))))/dx[0]/dx[0] +
	      (Temp_old_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k))) - 2.*Temp_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) + Temp_old_box(amrex::IntVect(AMREX_D_DECL(i,j-1,k))))/dx[1]/dx[1];
	    
	    amrex::Real rho = (rho1-rho0)*Eta_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) + rho0;
	    amrex::Real K   = (k1-k0)*Eta_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) + k0;
	    amrex::Real cp  = (cp1-cp0)*Eta_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) + cp0;


	    Temp_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) =
	      Temp_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) + (dt/rho/cp) * ((k1-k0)*(eta_gradx*T_gradx + eta_grady*T_grady)  + K*T_lap  + (w1 - w0 - qdotburn)*eta_grad_mag);

	    if (std::isnan(Temp_box(amrex::IntVect(AMREX_D_DECL(i,j,k)))))
	      amrex::Abort("NaN encountered");
	  }
    }
}



void PFFlame::TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{

  const Real* dx      = geom[lev].CellSize();

  const amrex::MultiFab& state = *Temp[lev];
  amrex::Array<int>  itags;
	
  for (amrex::MFIter mfi(state,true); mfi.isValid(); ++mfi)
    {
      const amrex::Box&  bx  = mfi.tilebox();

      amrex::TagBox&     tag  = tags[mfi];
	    
      //amrex::BaseFab<Real> &Eta_box = (*Temp[lev])[mfi];
      amrex::BaseFab<Real> &Eta_box = (*Eta[lev])[mfi];

      for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
	for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
#if BL_SPACEDIM>2
	  for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
#endif
	    {
	      amrex::Real gradx = (Eta_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k))) - Eta_box(amrex::IntVect(AMREX_D_DECL(i-1,j,k))))/(2.*dx[0]);
	      amrex::Real grady = (Eta_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k))) - Eta_box(amrex::IntVect(AMREX_D_DECL(i,j-1,k))))/(2.*dx[1]);
	      amrex::Real gradz = (Eta_box(amrex::IntVect(AMREX_D_DECL(i,j,k+1))) - Eta_box(amrex::IntVect(AMREX_D_DECL(i,j,k-1))))/(2.*dx[2]);
	      if (dx[0]*dx[1]*dx[2]*(gradx*gradx + grady*grady + gradz*gradz)>0.001) tag(amrex::IntVect(AMREX_D_DECL(i,j,k))) = amrex::TagBox::SET;
	    }
    }

}


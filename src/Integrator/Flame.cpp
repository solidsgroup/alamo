#include "Flame.H"
#include "BC/Constant.H"

namespace Integrator
{
Flame::Flame () : Integrator()
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
  pp.query("fs_number",fs_number);
  pp.query("fs_min",fs_min);
  pp.query("fs_max",fs_max);

  {
    amrex::ParmParse pp("TempBC");
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

    TempBC = new BC::Constant(bc_hi_str, bc_lo_str,
			      AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3),
			      AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3));
  }
  {
    amrex::ParmParse pp("EtaBC");
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

    EtaBC = new BC::Constant(bc_hi_str, bc_lo_str,
			     AMREX_D_DECL(bc_lo_1,bc_lo_2,bc_lo_3),
			     AMREX_D_DECL(bc_hi_1,bc_hi_2,bc_hi_3));
  }

  VoronoiIC = new IC::Voronoi(geom);
  std::vector<Set::Scalar> fs(fs_number);
  for (int i = 0; i < fs_number; i++) fs[i] = 0.5*(1.0 + Util::Random());
  VoronoiIC->Define(fs_number,fs,IC::Voronoi::Type::Values);

  EtaIC = new IC::Wedge(geom);
  // EtaIC = new IC::Constant(geom,value);

  RegisterNewFab(Temp,     TempBC, 1, 1, "Temp", true);
  RegisterNewFab(Temp_old, TempBC, 1, 1, "Temp_old", false);
  RegisterNewFab(Eta,      EtaBC,  1, 1, "Eta", true);
  RegisterNewFab(Eta_old,  EtaBC,  1, 1, "Eta_old", false);
  RegisterNewFab(FlameSpeedFab, EtaBC,  1, 1, "FlameSpeed",true);
}

void Flame::Initialize (int lev)
{
	for (amrex::MFIter mfi(*Temp[lev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.tilebox();

			amrex::BaseFab<amrex::Real> &Eta_box		= (*Eta[lev])[mfi];
			amrex::BaseFab<amrex::Real> &Eta_old_box	= (*Eta_old[lev])[mfi];
			amrex::BaseFab<amrex::Real> &Temp_box		= (*Temp[lev])[mfi];
			amrex::BaseFab<amrex::Real> &Temp_old_box	= (*Temp_old[lev])[mfi];

			AMREX_D_TERM(for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++),
							 for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++),
							 for (int k = box.loVect()[2]; k<=box.hiVect()[2]; k++))
				{
					Eta_box     (amrex::IntVect(AMREX_D_DECL(i,j,k))) =  1.0;
					Eta_old_box (amrex::IntVect(AMREX_D_DECL(i,j,k))) =  1.0;
					Temp_box    (amrex::IntVect(AMREX_D_DECL(i,j,k))) =  0.0;
					Temp_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) =  0.0;
				}
		}
	EtaIC->Initialize(lev,Eta);
	EtaIC->Initialize(lev,Eta_old);
	
	VoronoiIC->Initialize(lev,FlameSpeedFab);
}



void Flame::Advance (int lev, amrex::Real time, amrex::Real dt)
{
  std::swap(Eta_old [lev], Eta [lev]);
  std::swap(Temp_old[lev], Temp[lev]);

  static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
												 dy(AMREX_D_DECL(0,1,0)),
												 dz(AMREX_D_DECL(0,0,1)));

  const amrex::Real* DX = geom[lev].CellSize();

  amrex::Real a0=w0, a1=0.0, a2= -5*w1 + 16*w12 - 11*a0, a3=14*w1 - 32*w12 + 18*a0, a4=-8*w1 + 16*w12 - 8*a0;

  for ( amrex::MFIter mfi(*Temp[lev],true); mfi.isValid(); ++mfi )
    {
      const amrex::Box& bx = mfi.tilebox();

      amrex::FArrayBox &Eta_box		= (*Eta[lev])[mfi];
      amrex::FArrayBox &Eta_old_box		= (*Eta_old[lev])[mfi];
      amrex::FArrayBox &Temp_box		= (*Temp[lev])[mfi];
      amrex::FArrayBox &Temp_old_box	= (*Temp_old[lev])[mfi];
      amrex::FArrayBox &FlameSpeed	= (*FlameSpeedFab[lev])[mfi];


		AMREX_D_TERM(for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++),
						 for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++),
						 for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				//
				// Phase field evolution
				//

				amrex::Real M_dev = fs_min + FlameSpeed(m)*(fs_max-fs_min)/(amrex::Real)fs_number;

				amrex::Real eta_lap = AMREX_D_TERM((Eta_old_box(m+dx) - 2.*Eta_old_box(m) + Eta_old_box(m-dx))/DX[0]/DX[0],
															  + (Eta_old_box(m+dy) - 2.*Eta_old_box(m) + Eta_old_box(m-dy))/DX[1]/DX[1],
															  + (Eta_old_box(m+dz) - 2.*Eta_old_box(m) + Eta_old_box(m-dz))/DX[2]/DX[2]);
	     
	     
				amrex::Real oldeta = Eta_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k)));
				Eta_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) = oldeta -
					(M + M_dev) * dt * (a1 + 2*a2*oldeta + 3*a3*oldeta*oldeta + 4*a4*oldeta*oldeta*oldeta
								 - kappa*eta_lap);

				//
				// Temperature evolution
				//
	    
				//if (0) // todo ignore temperature for the time being
				//if (Eta_box(amrex::IntVect(i,j))>0.01)
				//{
				amrex::Real temperature_delay = 0.05;
				if (time<temperature_delay) continue;

				AMREX_D_TERM(amrex::Real eta_gradx = (Eta_old_box(m+dx) - Eta_old_box(m-dx))/(2*DX[0]);,
								 amrex::Real eta_grady = (Eta_old_box(m+dy) - Eta_old_box(m-dy))/(2*DX[1]);,
								 amrex::Real eta_gradz = (Eta_old_box(m+dz) - Eta_old_box(m-dz))/(2*DX[2]););
				
				AMREX_D_TERM(amrex::Real T_gradx = (Temp_old_box(m+dx) - Temp_old_box(m-dx))/(2*DX[0]);,
								 amrex::Real T_grady = (Temp_old_box(m+dy) - Temp_old_box(m-dy))/(2*DX[1]);,
								 amrex::Real T_gradz = (Temp_old_box(m+dz) - Temp_old_box(m-dz))/(2*DX[2]););

				amrex::Real eta_grad_mag = sqrt(AMREX_D_TERM(eta_gradx*eta_gradx, + eta_grady*eta_grady, + eta_gradz*eta_gradz));

				amrex::Real T_lap = 
					AMREX_D_TERM((Temp_old_box(m+dx) - 2.*Temp_old_box(m) + Temp_old_box(m-dx))/DX[0]/DX[0],
						     + (Temp_old_box(m+dy) - 2.*Temp_old_box(m) + Temp_old_box(m-dy))/DX[1]/DX[1],
						     + (Temp_old_box(m+dz) - 2.*Temp_old_box(m) + Temp_old_box(m-dz))/DX[2]/DX[2]);
	    
				amrex::Real rho = (rho1-rho0)*Eta_old_box(m) + rho0;
				amrex::Real K   = (k1-k0)*Eta_old_box(m) + k0;
				amrex::Real cp  = (cp1-cp0)*Eta_old_box(m) + cp0;

				Temp_box(m) =
					Temp_old_box(m)
					+ (dt/rho/cp) * ((k1-k0)*(AMREX_D_TERM(eta_gradx*T_gradx,
									       + eta_grady*T_grady,
									       + eta_gradz*T_gradz))
							 + K *T_lap  + (w1 - w0 - qdotburn)*eta_grad_mag);

				if (std::isnan(Temp_box(amrex::IntVect(AMREX_D_DECL(i,j,k)))))
					Util::Abort(INFO, "NaN encountered");
			}
    }
}



void Flame::TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{

	const amrex::Real* dx      = geom[lev].CellSize();

	const amrex::MultiFab& state = *Temp[lev];
	//amrex::Array<int>  itags;
	
	for (amrex::MFIter mfi(state,true); mfi.isValid(); ++mfi)
		{
			const amrex::Box&  bx  = mfi.tilebox();

			amrex::TagBox&     tag  = tags[mfi];
	    
			//amrex::BaseFab<Real> &Eta_box = (*Temp[lev])[mfi];
			amrex::BaseFab<amrex::Real> &Eta_box = (*Eta[lev])[mfi];

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
void Flame::Regrid(int lev, Set::Scalar /* time */)
{
	FlameSpeedFab[lev]->setVal(0.0);
	VoronoiIC->Initialize(lev,FlameSpeedFab);
	Util::Message(INFO,"Regridding on level ", lev);
}
}
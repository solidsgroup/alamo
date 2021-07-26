#include "Flame.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"

namespace Integrator
{

 
Flame::Flame () : Integrator()
{
  IO::ParmParse pp("physics");
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
  //pp.query("std_deviation",std_deviation);
  pp.query("R_min",R_min);
  pp.query("R_max",R_max);
  pp.query("volume_fraction",volume_fraction);
  //pp.query("mean",mean);
  pp.query("fs_min",fs_min);
  pp.query("fs_max",fs_max);
  
  pp.query("fs_ap",fs_ap);
  pp.query("fs_htpb",fs_htpb);
  pp.query("P",P);
  pp.query("r",r);
  pp.query("n",n);
  
  pp.query("refinement_criterion",m_refinement_criterion);
  {

    IO::ParmParse pp("bc");
     TempBC = new BC::Constant(1);
     pp.queryclass("temp", *static_cast<BC::Constant *>(TempBC));
     EtaBC = new BC::Constant(1);
     pp.queryclass("eta", *static_cast<BC::Constant *>(EtaBC));

  }

  VoronoiIC = new IC::Voronoi(geom);
  PackedSpheresIC = new IC::PackedSpheres(geom);
    std::vector<Set::Scalar> fs(fs_number,1);
  //std::vector<Set::Scalar> alpha(fs_number,1);
 // for (int i = 0; i < fs_number; i++)
// fs[i] = 0.5*(1.0 + Util::Random());
 //for (int i = 0; i < fs_number; i++)
 //std::cout<<"f[s]  "<<fs[i];
  VoronoiIC->Define(fs_number,fs,IC::Voronoi::Type::Values);
  //PackedSpheresIC->Define(fs_number,fs,mean,std_deviation,IC::PackedSpheres::Type::Values);
  
  PackedSpheresIC->Define(fs_number,fs,volume_fraction,R_min,R_max,IC::PackedSpheres::Type::Values);

  


  // EtaIC = new IC::Wedge(geom);

  // EtaIC = new IC::Constant(geom,value);

    RegisterNewFab(Temp_mf, TempBC, 1, 1, "Temp", true);
		RegisterNewFab(Temp_old_mf, TempBC, 1, 1, "Temp_old", false);
		RegisterNewFab(Eta_mf, EtaBC, 1, 1, "Eta", true);
		RegisterNewFab(Eta_old_mf, EtaBC, 1, 1, "Eta_old", false);
		//RegisterNewFab(FlameSpeed_mf, EtaBC, 1, 1, "FlameSpeed", true);
		RegisterNewFab(FlameSpeed_mf, EtaBC, 1, 1, "field", true);
}

void Flame::Initialize(int lev)
	{
		Temp_mf[lev]->setVal(1.0);
		Temp_old_mf[lev]->setVal(1.0);

		//EtaIC->Initialize(lev, Eta_mf);
		//EtaIC->Initialize(lev, Eta_old_mf);
		Eta_mf[lev]->setVal(1.0);  //initializing afte removing wedge BC
		Eta_old_mf[lev]->setVal(1.0); // initializing after removing wedge BC

                PackedSpheresIC->Initialize(lev, FlameSpeed_mf);
	}

void Flame::Advance(int lev, amrex::Real time, amrex::Real dt)
	{
		std::swap(Eta_old_mf[lev], Eta_mf[lev]);
		std::swap(Temp_old_mf[lev], Temp_mf[lev]);

		const amrex::Real *DX = geom[lev].CellSize();

		Set::Scalar a0 = w0, a1 = 0.0, a2 = -5 * w1 + 16 * w12 - 11 * a0, a3 = 14 * w1 - 32 * w12 + 18 * a0, a4 = -8 * w1 + 16 * w12 - 8 * a0;


		for (amrex::MFIter mfi(*Temp_mf[lev], true); mfi.isValid(); ++mfi)
		{
			const amrex::Box &bx = mfi.tilebox();
            
			amrex::Array4<Set::Scalar> const &Eta              = (*Eta_mf[lev]).array(mfi);
			amrex::Array4<const Set::Scalar> const &Eta_old    = (*Eta_old_mf[lev]).array(mfi);
			amrex::Array4<Set::Scalar> const &Temp             = (*Temp_mf[lev]).array(mfi);
			amrex::Array4<const Set::Scalar> const &Temp_old   = (*Temp_old_mf[lev]).array(mfi);
			//amrex::Array4<const Set::Scalar> const &FlameSpeed = (*FlameSpeed_mf[lev]).array(mfi);
			amrex::Array4<const Set::Scalar> const &field = (*FlameSpeed_mf[lev]).array(mfi);
			


			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
				//
				// Phase field evolution
				//
                		//Set::Scalar M_dev=field(i,j,k);
				//Set::Scalar M_dev = fs_min + field(i, j, k) *0.1* (fs_max - fs_min);
				Set:: Scalar fs_actual = fs_ap*field(i,j,k) + fs_htpb*(1-field(i,j,k));
				//Util::Message(INFO, "f_actual ",  fs_actual);

				Set::Scalar eta_lap = Numeric::Laplacian(Eta_old, i, j, k, 0, DX);

				Eta(i, j, k) = Eta_old(i, j, k) -
							   (fs_actual*((r*pow(P,n)))) * dt * (a1 + 2 * a2 * Eta_old(i, j, k) + 3 * a3 * Eta_old(i, j, k) * Eta_old(i, j, k) + 4 * a4 * Eta_old(i, j, k) * Eta_old(i, j, k) * Eta_old(i, j, k) - kappa * eta_lap);

				//
				// Temperature evolution
				//

				Set::Scalar temperature_delay = 0.01; // hard coded for now, need to make input
				if (time >= temperature_delay)
				{
					Set::Vector eta_grad = Numeric::Gradient(Eta_old, i, j, k, 0, DX);
					Set::Vector temp_grad = Numeric::Gradient(Temp_old, i, j, k, 0, DX);
					Set::Scalar temp_lap = Numeric::Laplacian(Temp_old, i, j, k, 0, DX);
					Set::Scalar eta_grad_mag = eta_grad.lpNorm<2>();
					amrex::Real rho = (rho1 - rho0) * Eta_old(i,j,k) + rho0;
					amrex::Real K = (k1 - k0) * Eta_old(i,j,k) + k0;
					amrex::Real cp = (cp1 - cp0) * Eta_old(i,j,k) + cp0;
					Set::Vector normvec = eta_grad/Eta_old(i,j,k);
					Set::Scalar test = normvec.dot(temp_grad);				
					Set::Scalar neumbound = 10.0;
				
					if (Eta_old(i,j,k) > 0.001 && Eta_old(i,j,k)<1)
						{
						Temp(i,j,k) = Temp_old(i,j,k) + dt*(K/cp/rho) * (test + temp_lap + eta_grad_mag/Eta_old(i,j,k) * neumbound);
						}	
					else
					{
					
						Temp(i,j,k) = Temp_old(i,j,k)+ dt*(K/cp/rho) * temp_lap;
					
					}
				}
			});
		}
	}

	void Flame::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/)
	{
		const amrex::Real *DX = geom[lev].CellSize();
		Set::Scalar dr  = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

		for (amrex::MFIter mfi(*Temp_mf[lev], true); mfi.isValid(); ++mfi)
		{
			const amrex::Box &bx = mfi.tilebox();
			amrex::Array4<char>              const &tags = a_tags.array(mfi);
			amrex::Array4<const Set::Scalar> const &Eta  = (*Eta_mf[lev]).array(mfi);

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
			{
				Set::Vector gradeta = Numeric::Gradient(Eta,i,j,k,0,DX);
				if (gradeta.lpNorm<2>() * dr > m_refinement_criterion)
					tags(i,j,k) = amrex::TagBox::SET;
				/*Set::Vector temp_grad = Numeric::Gradient(Temp
				if(temp_grad.lpNorm<2>() * dr > m_refinement_criterion)
					tags(i,j,k) = amrex::TagBox::SET;*/
			});
		}
	}
	void Flame::Regrid(int lev, Set::Scalar /* time */)
	{
		if (lev < finest_level) return;
		FlameSpeed_mf[lev]->setVal(0.0);
                PackedSpheresIC->Initialize(lev, FlameSpeed_mf);
		Util::Message(INFO, "Regridding on level ", lev);
	}
} // namespace Integrator

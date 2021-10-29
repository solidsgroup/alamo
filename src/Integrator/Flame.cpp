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
    pp.query("ka",ka);
    pp.query("kh",kh);
    pp.query("k0",k0);
    pp.query("cp1",cp1);
    pp.query("cp0",cp0);
    pp.query("qdotburn",qdotburn);
    pp.query("delA", delA);
    pp.query("delH", delH);
    pp.query("fs_number",fs_number);
    //pp.query("std_deviation",std_deviation);
    pp.query("R_min",R_min);
    pp.query("R_max",R_max);
    pp.query("volume_fraction",volume_fraction);
    pp.query("verbosity",verbosity);
    //pp.query("mean",mean);
    //pp.query("fs_min",fs_min);
    //pp.query("fs_max",fs_max);
    pp.query("R_mean",R_mean);
  
    pp.query("fs_ap",fs_ap);
    pp.query("fs_htpb",fs_htpb);
    pp.query("fs_comb",fs_comb);
    pp.query("P",P);
    pp.query("r_ap",r_ap);
    pp.query("n_ap",n_ap);
    pp.query("r_htpb",r_htpb);
    pp.query("n_htpb",n_htpb);
    pp.query("r_comb",r_comb);
    pp.query("n_comb",n_comb);
  
    pp.query("refinement_criterion",m_refinement_criterion);
    pp.query("temp_refinement", t_refinement_criterion);
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
  
    PackedSpheresIC->Define(fs_number,fs,volume_fraction,R_min,R_max,R_mean,verbosity,IC::PackedSpheres::Type::Values);


  

 // EtaIC = new IC::Wedge(geom);

    // EtaIC = new IC::Constant(geom,value);
    // TODO - replace 2 with nghost that is determined by elastic.on...
        RegisterNewFab(Temp_mf, TempBC, 1, 2, "Temp", true);
        RegisterNewFab(Temp_old_mf, TempBC, 1, 2, "Temp_old", false);
        RegisterNewFab(Eta_mf, EtaBC, 1, 2, "Eta", true);
        RegisterNewFab(Eta_old_mf, EtaBC, 1, 2, "Eta_old", false);
        //RegisterNewFab(FlameSpeed_mf, EtaBC, 1, 1, "FlameSpeed", true);
        RegisterNewFab(FlameSpeed_mf, EtaBC, 1, 2, "field", true);

        // Elastic solver
        RegisterNodalFab(elastic.disp_mf,AMREX_SPACEDIM,2,"disp",true);
        RegisterNodalFab(elastic.rhs_mf,AMREX_SPACEDIM,2,"rhs",true);
        RegisterNodalFab(elastic.res_mf,AMREX_SPACEDIM,2,"res",true);
        RegisterNodalFab(elastic.stress_mf,AMREX_SPACEDIM*AMREX_SPACEDIM, 2, "stress", true);
        RegisterGeneralFab(elastic.model_mf,1,2);
        {
            IO::ParmParse pp("elastic");
            pp.queryclass("ap",elastic.model_ap);
            pp.queryclass("htpb",elastic.model_htpb);
            pp.queryclass("void",elastic.model_void);
            pp.query("interval",elastic.interval);
            pp.queryclass("bc",elastic.bc);
        }
}

void Flame::Initialize(int lev)
    {
        Temp_mf[lev]->setVal(0.0);
        Temp_old_mf[lev]->setVal(0.0);

        //EtaIC->Initialize(lev, Eta_mf);
        //EtaIC->Initialize(lev, Eta_old_mf);
        Eta_mf[lev]->setVal(1.0);  //initializing afte removing wedge BC
        Eta_old_mf[lev]->setVal(1.0); // initializing after removing wedge BC

                                PackedSpheresIC->Initialize(lev, FlameSpeed_mf);
    }

void Flame::TimeStepBegin(Set::Scalar a_time, int a_iter) 
{
    if (!elastic.interval) return;
    if (a_iter%elastic.interval) return;

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        //Eta_mf[lev]->FillBoundary();
        //FlameSpeed_mf[lev]->FillBoundary();
        //Temp_mf[lev]->FillBoundary();
        Util::RealFillBoundary(*Eta_mf[lev],geom[lev]);
        Util::RealFillBoundary(*FlameSpeed_mf[lev],geom[lev]);
        Util::RealFillBoundary(*Temp_mf[lev],geom[lev]);

        elastic.rhs_mf[lev]->setVal(0.0);
        elastic.disp_mf[lev]->setVal(0.0);
        elastic.model_mf[lev]->setVal(elastic.model_ap);

                Set::Vector DX(geom[lev].CellSize());

                for (MFIter mfi(*elastic.model_mf[lev], true); mfi.isValid(); ++mfi)
                {

                        amrex::Box bx = mfi.nodaltilebox();
            bx.grow(1);
                        amrex::Array4<model_type> const &model = elastic.model_mf[lev]->array(mfi);
                        amrex::Array4<const Set::Scalar> const &eta = Eta_mf[lev]->array(mfi);
                        amrex::Array4<const Set::Scalar> const &phi = FlameSpeed_mf[lev]->array(mfi);
                        amrex::Array4<const Set::Scalar> const &temp = Temp_mf[lev]->array(mfi);

                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Scalar phi_avg = Numeric::Interpolate::CellToNodeAverage(phi,i,j,k,0);
                Set::Scalar eta_avg = Numeric::Interpolate::CellToNodeAverage(eta,i,j,k,0);
                Set::Scalar temp_avg = Numeric::Interpolate::CellToNodeAverage(temp,i,j,k,0);

                model_type model_ap = elastic.model_ap;
                model_ap.F0 *= temp_avg;
                model_type model_htpb = elastic.model_htpb;
                model_htpb.F0 *= temp_avg;

                model_type solid = model_ap*phi_avg + model_htpb*(1.-phi_avg);

                model(i,j,k) = solid*eta_avg + elastic.model_void*(1.0-eta_avg);
                        });
                }

        Util::RealFillBoundary(*elastic.model_mf[lev],geom[lev]);
    }
    elastic.bc.Init(elastic.rhs_mf,geom);

        amrex::LPInfo info;
        Operator::Elastic<Model::Solid::Affine::Isotropic::sym> elastic_op(Geom(0,finest_level), grids, DistributionMap(0,finest_level), info);
        elastic_op.SetUniform(false);
        elastic_op.SetBC(&elastic.bc);
    elastic_op.SetAverageDownCoeffs(true);

        Set::Scalar tol_rel = 1E-8, tol_abs = 1E-8;
        IO::ParmParse pp("elastic");
        elastic.solver = new Solver::Nonlocal::Newton<Model::Solid::Affine::Isotropic>(elastic_op);
        pp.queryclass("solver",*elastic.solver);

        elastic.solver->solve(elastic.disp_mf,elastic.rhs_mf,elastic.model_mf,tol_rel,tol_abs);
    elastic.solver->compResidual(elastic.res_mf,elastic.disp_mf,elastic.rhs_mf,elastic.model_mf);

    for (int lev = 0; lev <= elastic.disp_mf.finest_level; lev++)
    {
        const amrex::Real *DX = geom[lev].CellSize();
        for (MFIter mfi(*elastic.disp_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.nodaltilebox();
            amrex::Array4<model_type> const &model = elastic.model_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar> const &stress = elastic.stress_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const &disp = elastic.disp_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                                {
                                    std::array<Numeric::StencilType, AMREX_SPACEDIM>
                                        sten = Numeric::GetStencil(i, j, k, bx);
                                    if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::F)
                                    {
                                        Set::Matrix F = Set::Matrix::Identity() + Numeric::Gradient(disp, i, j, k, DX, sten);
                                        Set::Matrix P = model(i, j, k).DW(F);
                                        Numeric::MatrixToField(stress, i, j, k, P);
                                    }
                                    else
                                    {
                                        Set::Matrix gradu = Numeric::Gradient(disp, i, j, k, DX, sten);
                                        Set::Matrix sigma = model(i, j, k).DW(gradu);
                                        Set::Matrix eps = 0.5 * (gradu + gradu.transpose());
                                        Numeric::MatrixToField(stress, i, j, k, sigma);
                                    }
                                });
        }
    }
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

                                    fmod_ap=fs_ap*(r_ap*pow(P,n_ap));
                                    fmod_htpb=fs_htpb*(r_htpb*pow(P,n_htpb));
                                    //fmod_comb=fs_comb*(r_comb*pow(P,n_comb));
                                    //Util::Message(INFO, "fmod_ap  ", fmod_ap);
                                    //Util::Message(INFO, "fmod_htpb  ", fmod_htpb);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                //
                // Phase field evolution
                //
   

                Set:: Scalar fs_actual;

                /*if (field(i,j,k)==1)
                {
                    fs_actual = fmod_ap;

                }

                else if (field(i,j,k)==0)
                {
                        fs_actual = fmod_htpb;
                }*/


                //else if ((field(i,j,k)>0) && (field(i,j,k)<1))
                //{
                    //fs_actual = fmod_ap*field(i,j,k) + fmod_htpb*(1.0-field(i,j,k))+2*fmod_comb*field(i,j,k)*(1.0-field(i,j,k));
                    fs_actual = fmod_ap*field(i,j,k) + fmod_htpb*(1.0-field(i,j,k));
                    //}

                //Util::Message(INFO, "fs_actual  ", fs_actual);


                Set::Scalar eta_lap = Numeric::Laplacian(Eta_old, i, j, k, 0, DX);

                Eta(i, j, k) = Eta_old(i, j, k) -
                                (fs_actual) * dt * (a1 + 2 * a2 * Eta_old(i, j, k) + 3 * a3 * Eta_old(i, j, k) * Eta_old(i, j, k) + 4 * a4 * Eta_old(i, j, k) * Eta_old(i, j, k) * Eta_old(i, j, k) - kappa * eta_lap);

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
                    Set::Vector normvec = eta_grad/Eta_old(i,j,k);

                    amrex::Real rho = (rho1 - rho0) * Eta_old(i,j,k) + rho0;
                    amrex::Real Ka = (ka - k0) * Eta_old(i,j,k) + k0;
                    amrex::Real Kh = (kh -k0) * Eta_old(i,j,k) + k0;
                    Set:: Scalar K = Ka*field(i,j,k) + Kh*(1-field(i,j,k));

                    amrex::Real cp = (cp1 - cp0) * Eta_old(i,j,k) + cp0;

                    Set::Scalar test = normvec.dot(temp_grad);
                    Set:: Scalar neumbound = delA*field(i,j,k) + delH*(1-field(i,j,k));


                    if (Eta_old(i,j,k) > 0.001 && Eta_old(i,j,k)<1)
                        {
                            Temp(i,j,k) = Temp_old(i,j,k) + dt*(K/cp/rho) * (test + temp_lap + eta_grad_mag/Eta_old(i,j,k) * neumbound);
                        }
                    else
                        {
                            if(Eta_old(i,j,k)<=0.001)
                                {
                                    Temp(i,j,k)= 0;
                                }
                            else
                                {
                                    Temp(i,j,k) = Temp_old(i,j,k)+ dt*(K/cp/rho) * temp_lap;
                                }
                        }
                }
            });
        }
    }

    void Flame::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/)
    {
        const amrex::Real *DX = geom[lev].CellSize();
        Set::Scalar dr  = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));
        //Util::Message(INFO, "dr  ", dr);

        for (amrex::MFIter mfi(*Temp_mf[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.tilebox();
            amrex::Array4<char>              const &tags = a_tags.array(mfi);
            amrex::Array4<const Set::Scalar> const &Eta  = (*Eta_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const &Temp = (*Temp_mf[lev]).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Vector gradeta = Numeric::Gradient(Eta,i,j,k,0,DX);
                if (gradeta.lpNorm<2>() * dr *2 > m_refinement_criterion)
                    tags(i,j,k) = amrex::TagBox::SET;

                Set::Vector tempgrad = Numeric::Gradient(Temp, i,j,k,0,DX);
                if (tempgrad.lpNorm<2>() * dr > t_refinement_criterion)
                    tags(i,j,k) = amrex::TagBox::SET;
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

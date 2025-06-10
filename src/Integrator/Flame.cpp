#include "Flame.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "Model/Regression/Regression.H"
#include "Numeric/Stencil.H"
#include "IC/Laminate.H"
#include "IC/Constant.H"
#include "IC/PSRead.H"
#include "Numeric/Function.H"
#include "IC/Expression.H"
#include "IC/BMP.H"
#include "IC/PNG.H"
#include "Base/Mechanics.H"
#include "Util/Util.H"
#include "Model/Regression/PowerLaw.H"
#include "Model/Regression/Arrhenius.H"
#include "Model/Propellant/Propellant.H"
#include "Model/Propellant/Mesoscale.H"
#include "Model/Propellant/Homogenize.H"

#include <cmath>

namespace Integrator
{

Flame::Flame() : Base::Mechanics<model_type>() {}

Flame::Flame(IO::ParmParse& pp) : Flame()
{
    pp_queryclass(*this);
}


void
Flame::Forbids(IO::ParmParse& pp)
{
    pp.forbid("pressure.P","use chamber.pressure instead");

    pp.forbid("geometry.x_len","This is specified by geometry.prob_lo/hi");
    pp.forbid("geometry.y_len","This is specified by geometry.prob_lo/hi");
    pp.forbid("amr.ghost_cells", "This should not be adjustable ");

    pp.forbid("pf.gamma","implemented in regression rate law objects now"); 

    pp.forbid("pressure.r_ap",   "use regression.powerlaw.r_ap");
    pp.forbid("pressure.r_htpb", "use regression.powerlaw.r_htpb");
    pp.forbid("pressure.r_comb", "use regression.powerlaw.r_comb");
    pp.forbid("pressure.n_ap",   "use regression.powerlaw.n_ap");
    pp.forbid("pressure.n_htpb", "use regression.powerlaw.n_htpb");
    pp.forbid("pressure.n_comb", "use regression.powerlaw.n_comb");

    pp.forbid("thermal.m_ap",   "use regression.arrhenius.m_ap");
    pp.forbid("thermal.m_htpb", "use regression.arrhenius.m_htpb");
    pp.forbid("thermal.E_ap",   "use regression.arrhenius.E_ap");
    pp.forbid("thermal.E_htpb", "use regression.arrhenius.E_htpb");
    pp.forbid("thermal.modeling_ap",   "Old debug variable. Should equal 1 "); 
    pp.forbid("thermal.modeling_htpb", "Old debug variable. Should equal 1"); 

    pp.forbid("pressure.a1", "use propellant.mesoscale.a1 instead"); // Surgate heat flux model paramater - AP
    pp.forbid("pressure.a2", "use propellant.mesoscale.a2 instead"); // Surgate heat flux model paramater - HTPB
    pp.forbid("pressure.a3", "use propellant.mesoscale.a3 instead"); // Surgate heat flux model paramater - Total
    pp.forbid("pressure.b1", "use propellant.mesoscale.b1 instead"); // Surgate heat flux model paramater - AP
    pp.forbid("pressure.b2", "use propellant.mesoscale.b2 instead"); // Surgate heat flux model paramater - HTPB
    pp.forbid("pressure.b3", "use propellant.mesoscale.b3 instead"); // Surgate heat flux model paramater - Total
    pp.forbid("pressure.c1", "use propellant.mesoscale.c1 instead"); // Surgate heat flux model paramater - Total
    pp.forbid("pressure.mob_ap", "no longer used"); // Whether to include pressure to the arrhenius law
    pp.forbid("pressure.dependency", "use propellant.mesoscale.arrhenius_dependency"); // Whether to use pressure to determined the reference Zeta 
    pp.forbid("pressure.h1", "use propellant.mesoscale.h1 instead"); // Surgate heat flux model paramater - Homogenized
    pp.forbid("pressure.h2", "use propellant.mesoscale.h1 instead"); // 
    // AP maximum (or reference) mass flux value - See Meier and Schmidt et al 2024 eq. 16
    pp.forbid("thermal.mlocal_ap", "use propellant.mesoscale.mlocal_ap");
    pp.forbid("thermal.mlocal_comb", "use propellant.mesoscale.mlocal_comb");
    pp.forbid("thermal.mlocal_htpb", "this actually did **nothing** - it was overridden by a hard code using massfraction.");


    pp.forbid("thermal.disperssion1", "use propellant.mesoscale.dispersion1");
    pp.forbid("thermal.disperssion2", "use propellant.mesoscale.dispersion2");
    pp.forbid("thermal.disperssion3", "use propellant.mesoscale.dispersion3"); 


    pp.forbid("thermal.rho_ap", "use propellant.XX.rho_ap ");
    pp.forbid("thermal.rho_htpb","use propellant.XX.rho_htpb ");
    pp.forbid("thermal.k_ap",   "use propellant.XX.k_ap ");
    pp.forbid("thermal.k_htpb", "use propellant.XX.k_htpb ");
    pp.forbid("thermal.cp_ap", "use propellant.XX.cp_ap ");
    pp.forbid("thermal.cp_htpb","use propellant.XX.cp_htpb "); 

}


// [parser]
void
Flame::Parse(Flame& value, IO::ParmParse& pp)
{
    BL_PROFILE("Integrator::Flame::Flame()");

    // pp.pushPrefix("regression");
    // value.regression.Parse(value.regression,pp);
    // pp.popPrefix();
    // pp.select<Model::Regression::PowerLaw,Model::Regression::Arrhenius>
    //     ("regression",value.regression);


    Forbids(pp);

    // Whether to include extra fields (such as mdot, etc) in the plot output
    pp_query_default("plot_field",value.plot_field,true); 
        
    //
    // PHASE FIELD VARIABLES
    //
        
    pp_query_default("pf.eps", value.pf.eps, 0.0); // Burn width thickness
    pp_query_default("pf.kappa", value.pf.kappa, 0.0); // Interface energy param
    pp_query_default("pf.lambda", value.pf.lambda, 0.0); // Chemical potential multiplier
    pp_query_default("pf.w1", value.pf.w1, 0.0); // Unburned rest energy
    pp_query_default("pf.w12", value.pf.w12, 0.0);  // Barrier energy
    pp_query_default("pf.w0", value.pf.w0, 0.0);    // Burned rest energy

    // Boundary conditions for phase field order params
    pp.select<BC::Constant>("pf.eta.bc", value.bc_eta, 1 ); 
    value.RegisterNewFab(value.eta_mf, value.bc_eta, 1, 2, "eta", true);
    value.RegisterNewFab(value.eta_old_mf, value.bc_eta, 1, 2, "eta_old", false);

    // phase field initial condition
    pp.select<IC::Laminate,IC::Constant,IC::Expression,IC::BMP,IC::PNG>("pf.eta.ic",value.ic_eta,value.geom); 

    // regression rate model
    //pp.select<Model::Regression::PowerLaw,Model::Regression::Arrhenius>("regression",value.regression);

    // Whether to use the Thermal Transport Model
    pp_query_default("thermal.on", value.thermal.on, false); 

    // System Initial Temperature. 
    // TFluid and TElastic can be individually assigned, but will otherwise be assigned the value of Tref.
    pp_query_default("thermal.Tref", value.thermal.Tref, 300.0); 


    // Select reduced order model to capture heat feedback
    pp.select<Model::Propellant::PowerLaw, 
              Model::Propellant::Mesoscale,
              Model::Propellant::Homogenize>
        ("propellant",value.propellant);

    if (value.thermal.on) {


        // Used to change heat flux units
        pp_query_default("thermal.hc", value.thermal.hc, 1.0);
        // System AP mass fraction
        //pp_query_default("thermal.massfraction", value.thermal.massfraction, 0.8);

        // Effective fluid temperature
        pp_query_default("thermal.Tfluid", value.thermal.Tfluid, value.thermal.Tref); 


        //Temperature boundary condition
        pp.select_default<BC::Constant>("thermal.temp.bc", value.bc_temp, 1);
            
        value.RegisterNewFab(value.temp_mf, value.bc_temp, 1, 3, "temp", true);
        value.RegisterNewFab(value.temp_old_mf, value.bc_temp, 1, 3, "temp_old", false);
        value.RegisterNewFab(value.temps_mf, value.bc_temp, 1, 3, "temps", false);
        value.RegisterNewFab(value.temps_old_mf, value.bc_temp, 1, 3, "temps_old", false);

        value.RegisterNewFab(value.mdot_mf, value.bc_eta, 1, 3, "mdot", value.plot_field);
        value.RegisterNewFab(value.mob_mf, value.bc_eta, 1, 3, "mob", value.plot_field);
        value.RegisterNewFab(value.alpha_mf, value.bc_temp, 1, 3, "alpha", value.plot_field);
        value.RegisterNewFab(value.heatflux_mf, value.bc_temp, 1, 3, "heatflux", value.plot_field);
        value.RegisterNewFab(value.laser_mf, value.bc_temp, 1, 3, "laser", value.plot_field);

        value.RegisterIntegratedVariable(&value.chamber.volume, "volume");
        value.RegisterIntegratedVariable(&value.chamber.area, "area");
        value.RegisterIntegratedVariable(&value.chamber.massflux, "mass_flux");
        //value.RegisterIntegratedVariable(&value.chamber.pressure, "pressure", false);

        // laser initial condition
        pp.select_default<IC::Constant,IC::Expression>("laser.ic",value.ic_laser, value.geom);

        // thermal initial condition
        pp.select_default<IC::Constant,IC::Expression,IC::BMP,IC::PNG>("temp.ic",value.thermal.ic_temp,value.geom);
    }


    // Constant pressure value
    pp_query_default("chamber.pressure", value.chamber.pressure, 1.0); 


    // Whether to compute the pressure evolution
    pp_query_default("variable_pressure", value.variable_pressure, 0);
    // Whether to initialize Phi with homogenized properties
    pp_query_default("homogeneousSystem", value.homogeneousSystem, 0); 

    // Refinement criterion for eta field   
    pp_query_default("amr.refinement_criterion", value.m_refinement_criterion, 0.001);
    // Refinement criterion for temperature field    
    pp_query_default("amr.refinement_criterion_temp", value.t_refinement_criterion, 0.001);
    // Eta value to restrict the refinament for the temperature field 
    pp_query_default("amr.refinament_restriction", value.t_refinement_restriction, 0.1);
    // Refinement criterion for phi field [infinity]
    pp_query_default("amr.phi_refinement_criterion", value.phi_refinement_criterion, 1.0e100);
    // Lowest value of Eta.
    pp_query_default("small", value.small, 1.0e-8); 

    // Initial condition for $\phi$ field.
    pp.select_default<IC::PSRead,IC::Laminate,IC::Expression,IC::Constant,IC::BMP,IC::PNG>
        ("phi.ic",value.ic_phi,value.geom);


    // // in case we used the laminate IC, we will extract zeta from there.
    // pp_query("phi.ic.laminate.eps", value.zeta);
    // // or in case we used the psread IC, we will extract zeta from there.
    // pp_query("phi.ic.psread.eps", value.zeta); 

    value.RegisterNodalFab(value.phi_mf, 1, value.ghost_count + 1, "phi", true);

    // Whether to use Neo-hookean Elastic model
    pp_query_default("elastic.on", value.elastic.on, 0); 
    // Body force
    pp_query_default("elastic.traction", value.elastic.traction, 0.0); 
    // Phi refinement criteria 
    pp_query_default("elastic.phirefinement", value.elastic.phirefinement, 1); 

    pp.queryclass<Base::Mechanics<model_type>>("elastic",value);

    if (value.m_type != Type::Disable)
    {
        // Reference temperature for thermal expansion 
        // (temperature at which the material is strain-free)
        pp_query_default("Telastic", value.elastic.TElastic, value.thermal.Tref); 
        // elastic model of AP
        pp.queryclass<Model::Solid::Finite::NeoHookeanPredeformed>("model_ap", value.elastic.model_ap);
        // elastic model of HTPB
        pp.queryclass<Model::Solid::Finite::NeoHookeanPredeformed>("model_htpb", value.elastic.model_htpb);

        value.bc_psi = new BC::Nothing();
        value.RegisterNewFab(value.psi_mf, value.bc_psi, 1, value.ghost_count, "psi", value.plot_psi);
        value.psi_on = true;
    }

    bool allow_unused;
    // Set this to true to allow unused inputs without error.
    // (Not recommended.)
    pp.query_default("allow_unused",allow_unused,false);
    if (!allow_unused && pp.AnyUnusedInputs())
    {
        Util::Warning(INFO,"The following inputs were specified but not used:");
        pp.AllUnusedInputs();
        Util::Exception(INFO,"Aborting. Specify 'allow_unused=True` to ignore this error.");
    }
}

void Flame::Initialize(int lev)
{
    BL_PROFILE("Integrator::Flame::Initialize");
    Base::Mechanics<model_type>::Initialize(lev);

    ic_eta->Initialize(lev, eta_mf);
    ic_eta->Initialize(lev, eta_old_mf);
    ic_phi->Initialize(lev, phi_mf);
    //ic_phicell->Initialize(lev, phicell_mf);

    if (elastic.on) {
        psi_mf[lev]->setVal(1.0);
        rhs_mf[lev]->setVal(Set::Vector::Zero());
    }
    if (thermal.on) {
        if (thermal.ic_temp)
        {
            thermal.ic_temp->Initialize(lev,temp_mf);
            thermal.ic_temp->Initialize(lev,temp_old_mf);
            thermal.ic_temp->Initialize(lev,temps_mf);
            thermal.ic_temp->Initialize(lev,temps_old_mf);
        }
        else
        {
            temp_mf[lev]->setVal(thermal.Tref);
            temp_old_mf[lev]->setVal(thermal.Tref);
            temps_mf[lev]->setVal(thermal.Tref);
            temps_old_mf[lev]->setVal(thermal.Tref);
        }
        alpha_mf[lev]->setVal(0.0);
        mob_mf[lev]->setVal(0.0);
        mdot_mf[lev]->setVal(0.0);
        heatflux_mf[lev]->setVal(0.0);
        thermal.w1 = 0.2 * chamber.pressure + 0.9;
        ic_laser->Initialize(lev, laser_mf);
    }
    if (variable_pressure) chamber.pressure = 1.0;
}

void Flame::UpdateModel(int /*a_step*/, Set::Scalar /*a_time*/)
{
    if (m_type == Base::Mechanics<model_type>::Type::Disable) return;

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        amrex::Box domain = this->geom[lev].Domain();
        domain.convert(amrex::IntVect::TheNodeVector());
        const Set::Scalar* DX = geom[lev].CellSize();

        //psi_mf[lev]->setVal(1.0);
        phi_mf[lev]->FillBoundary();
        //phicell_mf[lev]->FillBoundary();
        eta_mf[lev]->FillBoundary();
        temp_mf[lev]->FillBoundary();

        for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box smallbox = mfi.nodaltilebox();
            amrex::Box bx = mfi.grownnodaltilebox() & domain;
            //amrex::Box bx = mfi.nodaltilebox();
            //bx.grow(1);
            Set::Patch<model_type>        model = model_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> phi   = phi_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> eta   = eta_mf.Patch(lev,mfi);
            Set::Patch<Set::Vector>       rhs   = rhs_mf.Patch(lev,mfi);
            // amrex::Array4<const Set::Scalar> const& Pressure = pressure_mf[lev]->array(mfi); // [error]

            if (elastic.on)
            {
                amrex::Array4<const Set::Scalar> const& temp = temp_mf[lev]->array(mfi);
                amrex::ParallelFor(smallbox, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Vector grad_eta = Numeric::CellGradientOnNode(eta, i, j, k, 0, DX);
                    rhs(i, j, k) = elastic.traction * grad_eta;
                });
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar phi_avg = phi(i, j, k, 0);
                    Set::Scalar temp_avg = Numeric::Interpolate::CellToNodeAverage(temp, i, j, k, 0);
                    model_type model_ap = elastic.model_ap;
                    model_ap.F0 -= Set::Matrix::Identity();
                    model_ap.F0 *= (temp_avg - elastic.TElastic);
                    model_ap.F0 += Set::Matrix::Identity();
                    model_type model_htpb = elastic.model_htpb;
                    model_htpb.F0 -= Set::Matrix::Identity();
                    model_htpb.F0 *= (temp_avg - elastic.TElastic);
                    model_htpb.F0 += Set::Matrix::Identity();

                    model(i, j, k) = model_ap * phi_avg + model_htpb * (1. - phi_avg);
                });
            }
            else
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar phi_avg = Numeric::Interpolate::CellToNodeAverage(phi, i, j, k, 0);
                    //phi_avg = phi(i,j,k,0);
                    model_type model_ap = elastic.model_ap;
                    model_ap.F0 *= Set::Matrix::Zero();
                    model_type model_htpb = elastic.model_htpb;
                    model_htpb.F0 *= Set::Matrix::Zero();
                    model(i, j, k) = model_ap * phi_avg + model_htpb * (1. - phi_avg);
                });
            }
        }
        Util::RealFillBoundary(*model_mf[lev], geom[lev]);

        psi_mf[lev]->setVal(1.0);
        amrex::MultiFab::Copy(*psi_mf[lev], *eta_mf[lev], 0, 0, 1, psi_mf[lev]->nGrow());
    }
}

void Flame::TimeStepBegin(Set::Scalar a_time, int a_iter)
{
    BL_PROFILE("Integrator::Flame::TimeStepBegin");
    Base::Mechanics<model_type>::TimeStepBegin(a_time, a_iter);
    if (thermal.on) {
        for (int lev = 0; lev <= finest_level; ++lev)
            ic_laser->Initialize(lev, laser_mf, a_time);
    }
}

void Flame::TimeStepComplete(Set::Scalar /*a_time*/, int /*a_iter*/)
{
    BL_PROFILE("Integrator::Flame::TimeStepComplete");
    if (variable_pressure) {
        Set::Scalar x_len = geom[0].ProbDomain().length(0);
        Set::Scalar y_len = geom[0].ProbDomain().length(1);
        Set::Scalar domain_area = x_len * y_len;
        Util::Message(INFO, "Mass = ", chamber.massflux);
        Util::Message(INFO, "Pressure = ", chamber.pressure);
    }
}

void Flame::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrador::Flame::Advance");
    Base::Mechanics<model_type>::Advance(lev, time, dt);
    const Set::Scalar* DX = geom[lev].CellSize();

    std::swap(eta_old_mf[lev], eta_mf[lev]);

    //
    // Chamber pressure update
    //
    if (variable_pressure) {
        chamber.pressure = exp(0.00075 * chamber.massflux);
        if (chamber.pressure > 10.0) {
            chamber.pressure = 10.0;
        }
        else if (chamber.pressure <= 0.99) {
            chamber.pressure = 0.99;
        }
        elastic.traction = chamber.pressure;
    }


    //
    // Multi-well chemical potential
    //
    Numeric::Function::Polynomial<4> w( pf.w0,
                                        0.0,
                                        -5.0 * pf.w1 + 16.0 * pf.w12 - 11.0 * pf.w0,
                                        14.0 * pf.w1 - 32.0 * pf.w12 + 18.0 * pf.w0,
                                        -8.0 * pf.w1 + 16.0 * pf.w12 - 8.0 * pf.w0);
    Numeric::Function::Polynomial<3> dw = w.D();

    propellant.set_pressure(chamber.pressure);

    for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        // Phase fields
        Set::Patch<Set::Scalar> etanew    = eta_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> eta = eta_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> phi = phi_mf.Patch(lev,mfi);
        // Heat transfer fields
        Set::Patch<const Set::Scalar> temp = temp_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       alpha = alpha_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       laser = laser_mf.Patch(lev,mfi);
        // Diagnostic fields
        Set::Patch<Set::Scalar> mdot     = mdot_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> heatflux = heatflux_mf.Patch(lev,mfi);


        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            //
            // CALCULATE PHI-AVERAGED QUANTITIES
            //
            Set::Scalar phi_avg = Numeric::Interpolate::NodeToCellAverage(phi, i, j, k, 0);
            Set::Scalar T = thermal.on ? temp(i,j,k) : NAN;

            Set::Scalar K = propellant.get_K(phi_avg);
            Set::Scalar rho = propellant.get_rho(phi_avg);
            Set::Scalar cp = propellant.get_cp(phi_avg);

            //
            // CALCULATE MOBILITY
            // 
            Set::Scalar L = propellant.get_L(  phi_avg, T);
            Util::AssertException(INFO,TEST(!isnan(L)));

            // 
            // EVOLVE PHASE FIELD (ETA)
            // 

            Set::Scalar eta_lap = Numeric::Laplacian(eta, i, j, k, 0, DX);
            Set::Scalar df_deta = ((pf.lambda / pf.eps) * dw(eta(i, j, k)) - pf.eps * pf.kappa * eta_lap);
            etanew(i, j, k) = eta(i, j, k) - L * dt * df_deta;
            if (etanew(i, j, k) <= small) etanew(i, j, k) = small;
            Util::AssertException (INFO, TEST(!isnan(etanew(i, j, k))));
            
            if (thermal.on)
            {
                //
                // Calculate thermal diffisivity and store for later gradient
                //

                alpha(i, j, k) = K / rho / cp; 
                if (isnan(alpha(i, j, k))) Util::Exception(INFO); 


                //
                // CALCULATE MASS FLUX BASED ON EVOLVING ETA
                //
            
                mdot(i, j, k) = rho * fabs(eta(i, j, k) - etanew(i, j, k)) / dt; 
                if (isnan(mdot(i, j, k))) Util::Exception(INFO);


                //
                // CALCULATE HEAT FLUX BASED ON THE CALCULATED MASS FLUX
                //

                Set::Scalar q0 = propellant.get_qdot(mdot(i,j,k), phi_avg);
                heatflux(i,j,k) = ( thermal.hc*q0 + laser(i,j,k) ) / K;
                if (isnan(heatflux(i, j, k))) Util::Exception(INFO);
            }
        });

    } // MFi For loop 


    //
    // THERMAL TRANSPORT
    // 
    if (thermal.on)
    {
        std::swap(temp_old_mf[lev], temp_mf[lev]);
        std::swap(temps_old_mf[lev], temps_mf[lev]);

        for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();

            Set::Patch<Set::Scalar>       tempnew = temp_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar>       tempsnew = temps_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> temp = temp_old_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> alpha = alpha_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> temps = temps_old_mf.Patch(lev,mfi);
            // Phase field
            Set::Patch<Set::Scalar>       etanew = (*eta_mf[lev]).array(mfi);
            Set::Patch<const Set::Scalar> eta = (*eta_old_mf[lev]).array(mfi);
            // Diagnostic fields
            Set::Patch<const Set::Scalar> heatflux = heatflux_mf.Patch(lev,mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                auto sten = Numeric::GetStencil(i, j, k, bx);
                Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
                Set::Vector grad_temp = Numeric::Gradient(temp, i, j, k, 0, DX);
                Set::Scalar lap_temp = Numeric::Laplacian(temp, i, j, k, 0, DX);
                Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>();
                Set::Vector grad_alpha = Numeric::Gradient(alpha, i, j, k, 0, DX, sten);
                Set::Scalar dTdt = 0.0;
                dTdt += grad_eta.dot(grad_temp * alpha(i, j, k));
                dTdt += grad_alpha.dot(eta(i, j, k) * grad_temp);
                dTdt += eta(i, j, k) * alpha(i, j, k) * lap_temp;
                dTdt += alpha(i, j, k) * heatflux(i, j, k) * grad_eta_mag;

                Set::Scalar Tsolid = dTdt + temps(i, j, k) * (etanew(i, j, k) - eta(i, j, k)) / dt;
                tempsnew(i, j, k) = temps(i, j, k) + dt * Tsolid;
                tempnew(i, j, k) = etanew(i, j, k) * tempsnew(i, j, k) + (1.0 - etanew(i, j, k)) * thermal.Tfluid;
                if (isnan(tempsnew(i, j, k)) || isnan(temps(i, j, k))) {
                    Util::Message(INFO, tempsnew(i, j, k), "tempsnew contains nan (i=", i, " j= ", j, ")");
                    Util::Message(INFO, temps(i, j, k), "temps contains nan (i=", i, " j= ", j, ")");
                    Util::Abort(INFO);
                }
            });
        }
    }
 
} //Function


void Flame::TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar time, int ngrow)
{
    BL_PROFILE("Integrator::Flame::TagCellsForRefinement");
    Base::Mechanics<model_type>::TagCellsForRefinement(lev, a_tags, time, ngrow);

    const Set::Scalar* DX = geom[lev].CellSize();
    Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

    // Eta criterion for refinement
    for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tags = a_tags.array(mfi);
        amrex::Array4<const Set::Scalar> const& eta = (*eta_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Vector gradeta = Numeric::Gradient(eta, i, j, k, 0, DX);
            if (gradeta.lpNorm<2>() * dr * 2 > m_refinement_criterion && eta(i, j, k) >= t_refinement_restriction)
                tags(i, j, k) = amrex::TagBox::SET;
        });
    }

    // Phi criterion for refinement 
    if (elastic.phirefinement) {
        for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<char> const& tags = a_tags.array(mfi);
            amrex::Array4<const Set::Scalar> const& phi = (*phi_mf[lev]).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Vector gradphi = Numeric::Gradient(phi, i, j, k, 0, DX);
                if (gradphi.lpNorm<2>() * dr >= phi_refinement_criterion)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    }


    // Thermal criterion for refinement 
    if (thermal.on) {
        for (amrex::MFIter mfi(*temp_mf[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<char> const& tags = a_tags.array(mfi);
            amrex::Array4<const Set::Scalar> const& temp = (*temp_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& eta = (*eta_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Vector tempgrad = Numeric::Gradient(temp, i, j, k, 0, DX);
                if (tempgrad.lpNorm<2>() * dr > t_refinement_criterion && eta(i, j, k) >= t_refinement_restriction)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    }
}

void Flame::Regrid(int lev, Set::Scalar time)
{
    BL_PROFILE("Integrator::Flame::Regrid");
    //if (lev < finest_level) return;
    //phi_mf[lev]->setVal(0.0);
    ic_phi->Initialize(lev, phi_mf, time);
    //ic_phicell->Initialize(lev, phi_mf, time);
}

void Flame::Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/,
    const amrex::MFIter& mfi, const amrex::Box& box)
{
    BL_PROFILE("Flame::Integrate");
    const Set::Scalar* DX = geom[amrlev].CellSize();
    Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
    amrex::Array4<amrex::Real> const& eta = (*eta_mf[amrlev]).array(mfi);
    amrex::Array4<amrex::Real> const& mdot = (*mdot_mf[amrlev]).array(mfi);
    if (variable_pressure) {
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            chamber.volume += eta(i, j, k, 0) * dv;
            Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
            Set::Scalar normgrad = grad.lpNorm<2>();
            Set::Scalar da = normgrad * dv;
            chamber.area += da;

            Set::Vector mgrad = Numeric::Gradient(mdot, i, j, k, 0, DX);
            Set::Scalar mnormgrad = mgrad.lpNorm<2>();
            Set::Scalar dm = mnormgrad * dv;
            chamber.massflux += dm;

        });
    }
    else {
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            chamber.volume += eta(i, j, k, 0) * dv;
            Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
            Set::Scalar normgrad = grad.lpNorm<2>();
            Set::Scalar da = normgrad * dv;
            chamber.area += da;
        });
    }
    // time dependent pressure data from experimenta -> p = 0.0954521220950523 * exp(15.289993148880678 * t)
}
} // namespace Integrator



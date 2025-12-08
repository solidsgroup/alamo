#include "Flame.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
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
#include "Model/Propellant/Propellant.H"
#include "Model/Propellant/FullFeedback.H"
#include "Model/Propellant/Homogenize.H"

#include <cmath>

namespace Integrator
{

Flame::Flame() : 
    Base::Mechanics<model_type>(), 
    Hydro(eta_mf, eta_old_mf, true)
{}

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

    pp.forbid("pf.gamma","use propellant.powerlaw.gamma"); 

    pp.forbid("pressure.r_ap",   "use propellant.powerlaw.r_ap");
    pp.forbid("pressure.r_htpb", "use propellant.powerlaw.r_htpb");
    pp.forbid("pressure.r_comb", "use propellant.powerlaw.r_comb");
    pp.forbid("pressure.n_ap",   "use propellant.powerlaw.n_ap");
    pp.forbid("pressure.n_htpb", "use propellant.powerlaw.n_htpb");
    pp.forbid("pressure.n_comb", "use propellant.powerlaw.n_comb");

    pp.forbid("thermal.bound",   "use thermal.Tref");
    pp.forbid("thermal.T_fluid",   "use thermal.Tfluid (or nothing)");
    pp.forbid("thermal.m_ap",   "use propellant.fullfeedback.m_ap");
    pp.forbid("thermal.m_htpb", "use propellant.fullfeedback.m_htpb");
    pp.forbid("thermal.E_ap",   "use propellant.fullfeedback.E_ap");
    pp.forbid("thermal.E_htpb", "use propellant.fullfeedback.E_htpb");
    pp.forbid("thermal.modeling_ap",   "Old debug variable. Should equal 1 "); 
    pp.forbid("thermal.modeling_htpb", "Old debug variable. Should equal 1"); 

    pp.forbid("pressure.a1", "use propellant.fullfeedback.a1 instead");
    pp.forbid("pressure.a2", "use propellant.fullfeedback.a2 instead");
    pp.forbid("pressure.a3", "use propellant.fullfeedback.a3 instead");
    pp.forbid("pressure.b1", "use propellant.fullfeedback.b1 instead");
    pp.forbid("pressure.b2", "use propellant.fullfeedback.b2 instead");
    pp.forbid("pressure.b3", "use propellant.fullfeedback.b3 instead");
    pp.forbid("pressure.c1", "use propellant.fullfeedback.c1 instead");
    pp.forbid("pressure.mob_ap", "no longer used"); 
    pp.forbid("pressure.dependency", "use propellant.fullfeedback.pressure_dependency"); 
    pp.forbid("pressure.h1", "use propellant.homogenize.h1 instead"); 
    pp.forbid("pressure.h2", "use propellant.homogenize.h2 instead"); 
    pp.forbid("thermal.mlocal_ap", "use propellant.homogenize.mlocal_ap");
    pp.forbid("thermal.mlocal_comb", "use propellant.homogenize.mlocal_comb");
    pp.forbid("thermal.mlocal_htpb", "this actually did **nothing** - it was overridden by a hard code using massfraction.");

    pp.forbid("thermal.disperssion1", "use propellant.homogenize.dispersion1");
    pp.forbid("thermal.disperssion2", "use propellant.homogenize.dispersion2");
    pp.forbid("thermal.disperssion3", "use propellant.homogenize.dispersion3"); 

    pp.forbid("thermal.rho_ap", "use propellant.fullfeedback/homogenize.rho_ap ");
    pp.forbid("thermal.rho_htpb","use propellant.fullfeedback/homogenize.rho_htpb ");
    pp.forbid("thermal.k_ap",   "use propellant.fullfeedback/homogenize.k_ap ");
    pp.forbid("thermal.k_htpb", "use propellant.fullfeedback/homogenize.k_htpb ");
    pp.forbid("thermal.cp_ap", "use propellant.fullfeedback/homogenize.cp_ap ");
    pp.forbid("thermal.cp_htpb","use propellant.fullfeedback/homogenize.cp_htpb "); 
}


// [parser]
void
Flame::Parse(Flame& value, IO::ParmParse& pp)
{
    BL_PROFILE("Integrator::Flame::Flame()");

    Forbids(pp);

    // Whether to include extra fields (such as mdot, etc) in the plot output
    pp.query_default("plot_field",value.plot_field,true); 
        
    //
    // PHASE FIELD VARIABLES
    //
        
    // Burn width thickness
    pp.query_default("pf.eps", value.pf.eps, "1.0_m", Unit::Length()); 
    // Interface energy param
    pp.query_default("pf.kappa", value.pf.kappa, "0.0_J/m^2", Unit::Energy() / Unit::Area()); 
    // Chemical potential multiplier
    pp.query_default("pf.lambda", value.pf.lambda, "0.0_J/m^2", Unit::Energy()/Unit::Area()); 
    // Unburned rest energy
    pp.query_default("pf.w1", value.pf.w1, "0.0",Unit::Less()); 
    // Barrier energy
    pp.query_default("pf.w12", value.pf.w12, "0.0", Unit::Less());  
    // Burned rest energy
    pp.query_default("pf.w0", value.pf.w0, "0.0",Unit::Less());    

    // Boundary conditions for phase field order params
    pp.select<BC::Constant>("pf.eta.bc", value.bc_eta, 1 ); 
    value.RegisterNewFab(value.eta_mf, value.bc_eta, 1, 2, "eta", true);
    value.RegisterNewFab(value.eta_old_mf, value.bc_eta, 1, 2, "eta_old", false);

    // phase field initial condition
    pp.select<IC::Laminate,IC::Constant,IC::Expression,IC::BMP,IC::PNG>("pf.eta.ic",value.ic_eta,value.geom); 



    // Select reduced order model to capture heat feedback
    pp.select<  Model::Propellant::PowerLaw, 
                Model::Propellant::FullFeedback,
                Model::Propellant::Homogenize>
        ("propellant",value.propellant);


    // Whether to use the Thermal Transport Model
    pp_query_default("thermal.on", value.thermal.on, false); 

    // Reference temperature
    // Used to set all other reference temperatures by default.
    pp_query_default("thermal.Tref", value.thermal.Tref, "300.0_K",Unit::Temperature()); 

    if (value.thermal.on) {

        // Used to change heat flux units
        pp_query_default("thermal.hc", value.thermal.hc, "1.0", Unit::Power()/Unit::Area());

        // Effective fluid temperature
        pp_query_default("thermal.Tfluid", value.thermal.Tfluid, value.thermal.Tref); 

        //Temperature boundary condition
        pp.select_default<BC::Constant>("thermal.temp.bc", value.bc_temp, 1, Unit::Temperature());
            
        value.RegisterNewFab(value.temp_mf, value.bc_temp, 1, 3, "temp", true);
        value.RegisterNewFab(value.temp_old_mf, value.bc_temp, 1, 3, "temp_old", false);
        value.RegisterNewFab(value.temps_mf, value.bc_temp, 1, 0, "temps", false);

        value.RegisterNewFab(value.mdot_mf, value.bc_temp, 1, 0, "mdot", value.plot_field);
        value.RegisterNewFab(value.alpha_mf, value.bc_temp, 1, 0, "alpha", value.plot_field);
        value.RegisterNewFab(value.heatflux_mf, value.bc_temp, 1, 0, "heatflux", value.plot_field);
        value.RegisterNewFab(value.laser_mf, value.bc_temp, 1, 0, "laser", value.plot_field);
        value.RegisterNewFab(value.rho_htpb_gas_mf, value.bc_temp, 1, 0, "rho_HTPB_gas", value.plot_field); // Density of gaseous Hydroxyl-terminated polybutadiene (HTPB)
        value.RegisterNewFab(value.rho_AP_gas_mf, value.bc_temp, 1, 0, "rho_AP_gas", value.plot_field); // Density of gaseous Ammonium perchlorate (AP)
        value.RegisterNewFab(value.rho_tot_gas_mf, value.bc_temp, 1, 0, "rho_tot_gas", value.plot_field); // Density of total gaseous field

        value.RegisterIntegratedVariable(&value.chamber.volume, "volume");
        value.RegisterIntegratedVariable(&value.chamber.area, "area");
        value.RegisterIntegratedVariable(&value.chamber.massflux, "mass_flux");

        // laser initial condition
        pp.select_default<  IC::Constant,
                            IC::Expression  >
            ("laser.ic",value.ic_laser, value.geom, Unit::Power()/Unit::Area());

        // thermal initial condition
        pp.select_default<  IC::Constant,
                            IC::Expression,
                            IC::BMP,
                            IC::PNG  >
            ("temp.ic",value.thermal.ic_temp,value.geom, Unit::Temperature());
    }


    // Constant pressure value
    pp_query_default("chamber.pressure", value.chamber.pressure, "1.0_MPa", Unit::Pressure()); 

    // Whether to compute the pressure evolution
    pp_query_default("variable_pressure", value.variable_pressure, false);

    // Flag to use Flame with or without Hydro. If Hydro is off the pressure traction from Hydro is not calculated (see UpdateModel function)
    pp_query_default("use_with_Hydro", value.use_with_Hydro, false);

    if (value.use_with_Hydro == 1) {
        //Scalar value to reduce the velocity of the products from regression in interface region being ejected into Hydro to improve Hydro stability
        pp.query_default("velocity_mult", value.velocity_mult, 1.0);
    }

    // Refinement criterion for eta field   
    pp_query_default(   "amr.refinement_criterion", value.m_refinement_criterion, "0.001", 
                        Unit::Less());

    // Refinement criterion for temperature field    
    pp.query_default(   "amr.refinement_criterion_temp", value.t_refinement_criterion, "0.001_K",
                        Unit::Temperature());

    // Eta value to restrict the refinament for the temperature field 
    pp.query_default(   "amr.refinament_restriction", value.t_refinement_restriction, "0.1",
                        Unit::Less());

    // Refinement criterion for phi field [infinity]
    pp_query_default("amr.phi_refinement_criterion", value.phi_refinement_criterion, 1.0e100);

    // Minimum allowable threshold for $\eta$
    pp_query_default("small", value.small, 1.0e-8); 

    // Initial condition for $\phi$ field.
    pp.select_default<IC::PSRead,IC::Laminate,IC::Expression,IC::Constant,IC::BMP,IC::PNG>
        ("phi.ic",value.ic_phi,value.geom);

    value.RegisterNodalFab(value.phi_mf, 1, 2, "phi", true);

    // Whether to use Neo-hookean Elastic model
    pp_query_default("elastic.on", value.elastic.on, 0); 

    // Body force
    pp_query_default("elastic.traction", value.elastic.traction, 0.0); 
    
    // Scalar value to reduce the pressure being applied in the traction boundaay condition. If a linear elastic solve is used, results can be scaled appropriately to recover correct solution
    pp.query_default("elastic.pressure_mult", value.elastic.pressure_mult, 1.0);

    // Phi refinement criteria 
    pp_query_default("elastic.phirefinement", value.elastic.phirefinement, 1); 

    pp.queryclass<Base::Mechanics<model_type>>("elastic",value);

    if (value.m_type != Type::Disable)
    {
        // Reference temperature for thermal expansion 
        // (temperature at which the material is strain-free)
        pp_query_default("Telastic", value.elastic.Telastic, value.thermal.Tref); 
        // elastic model of AP
        pp.queryclass<Model::Solid::Finite::NeoHookeanPredeformed>("model_ap", value.elastic.model_ap);
        // elastic model of HTPB
        pp.queryclass<Model::Solid::Finite::NeoHookeanPredeformed>("model_htpb", value.elastic.model_htpb);

        // Use our current eta field as the psi field for the solver
        value.psi_on = false;
        value.solver.setPsi(value.eta_mf);
    }

    // Determines if hydro is on, by default hydro starts off and turns on a small time into the simulation to aid hydro initilization
    pp.query_default("hydro.on",value.hydro.on,false);
    if (value.hydro.on)
    {   
        // Time value when hydro starts, recommended to be set slightly larger than the start of the simulation to help integrators start correctly
        pp.query_default("hydro.tstart", value.hydro.tstart, 0.0);

        // Density of the gaseous products of Ammonium perchlorate
        pp.query_default("hydro.rho_ap",value.hydro.rho_ap,1.0);

        // Density of the gaseous products of Hydroxyl-terminated Polybutadiene
        pp.query_default("hydro.rho_htpb",value.hydro.rho_htpb,1.0);

        // Velocity source term of gaseous AP products, gets overridden in code by conservation of mass source terms
        pp.query_default("hydro.u0_ap",value.hydro.u0_ap,0.0);

        // Velocity source term of gaseous HTPB products, gets overridden in code by conservation of mass source terms
        pp.query_default("hydro.u0_htpb",value.hydro.u0_htpb,0.0);

        pp.queryclass<Hydro>("hydro", value);
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
    if (hydro.on) Hydro::Initialize(lev);

    ic_eta->Initialize(lev, eta_mf);
    ic_eta->Initialize(lev, eta_old_mf);
    ic_phi->Initialize(lev, phi_mf);
    //ic_phicell->Initialize(lev, phicell_mf);

    if (elastic.on) {
        rhs_mf[lev]->setVal(Set::Vector::Zero());
    }
    if (thermal.on) {
        if (thermal.ic_temp)
        {
            thermal.ic_temp->Initialize(lev,temp_mf);
            thermal.ic_temp->Initialize(lev,temp_old_mf);
            thermal.ic_temp->Initialize(lev,temps_mf);
        }
        else
        {
            temp_mf[lev]->setVal(thermal.Tref);
            temp_old_mf[lev]->setVal(thermal.Tref);
            temps_mf[lev]->setVal(thermal.Tref);
        }
        alpha_mf[lev]->setVal(0.0);
        mdot_mf[lev]->setVal(0.0);
        heatflux_mf[lev]->setVal(0.0);
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

        phi_mf[lev]->FillBoundary();
        eta_mf[lev]->FillBoundary();
        temp_mf[lev]->FillBoundary();

        for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box smallbox = mfi.nodaltilebox();
            amrex::Box bx = mfi.grownnodaltilebox() & domain;
            Set::Patch<model_type>        model = model_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> phi   = phi_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> eta   = eta_mf.Patch(lev,mfi);
            Set::Patch<Set::Vector>       rhs   = rhs_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> pressure = Hydro::pressure_mf.Patch(lev,mfi); // Pressure from Hydro to use as traction force at solid/fluid interface

            if (elastic.on)
            {
                Set::Patch <const Set::Scalar> temp = temp_mf.Patch(lev,mfi);
                amrex::ParallelFor(smallbox, [=] AMREX_GPU_DEVICE(int i, int j, int k)

                {   
                    Set::Vector grad_eta = Numeric::CellGradientOnNode(eta, i, j, k, 0, DX);
                    Set::Vector pres_reg;
                    
                    if (use_with_Hydro) {
                        pres_reg = elastic.pressure_mult*pressure(i,j,k) * grad_eta ; // Add pressure to effect the regression rate
                    } else {
                        // If not using hydro to find the pressure on the regressing surface, set the pressure traction term to 0.0
                        pres_reg = 0.0*grad_eta;
                    }

                    rhs(i, j, k) = 0.0*grad_eta;
                    rhs(i, j, k) = (elastic.traction) * grad_eta - pres_reg;
                });
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar phi_avg = phi(i, j, k, 0);
                    Set::Scalar temp_avg = Numeric::Interpolate::CellToNodeAverage(temp, i, j, k, 0);
                    model_type model_ap = elastic.model_ap;
                    model_ap.F0 -= Set::Matrix::Identity();
                    model_ap.F0 *= (temp_avg - elastic.Telastic);
                    model_ap.F0 += Set::Matrix::Identity();
                    model_type model_htpb = elastic.model_htpb;
                    model_htpb.F0 -= Set::Matrix::Identity();
                    model_htpb.F0 *= (temp_avg - elastic.Telastic);
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

    }
}

void Flame::UpdateFluxes(int lev, Set::Scalar a_time, Set::Scalar dt)
{
    amrex::Box domain = this->geom[lev].Domain();
    domain.convert(amrex::IntVect::TheNodeVector());
    const Set::Scalar* DX = geom[lev].CellSize();

    for (MFIter mfi(*eta_mf[lev], false); mfi.isValid(); ++mfi)
    {
        amrex::Box bx = mfi.tilebox();
        Set::Patch<const Set::Scalar> phi_patch    = phi_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> eta    = eta_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> etaold = eta_old_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> rho_AP_gas = rho_AP_gas_mf.Patch(lev,mfi); // Set scalar value for density of AP
        Set::Patch<Set::Scalar> rho_HTPB_gas = rho_htpb_gas_mf.Patch(lev,mfi); // Set scalar value for density of HTPB
        Set::Patch<Set::Scalar> rho_tot_gas = rho_tot_gas_mf.Patch(lev,mfi); // Set scalar value for total density of the gaseous phase
        Set::Patch<Set::Scalar> pressure = Hydro::pressure_mf.Patch(lev,mfi); // Call the pressure from the Hydro integrator
        Set::Patch<Set::Scalar> u0 = Hydro::u0_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> p = Hydro::pressure_mf.Patch(lev,mfi);

        Real M_AP = 27.645; // Molar mass of mixture after AP undergos pyrolysis (kg/mol)
        Real M_HTPB = 28.0532; // Molar mass of Ethylene, main product of HTPB pyrolysis
        Real R = 8314; // Ideal gas constant (J/kmol-k)
        Real Pref = Hydro::pref; // Find the reference temperature from Hydro
        Real temp_gas = 750; // Set value for temperature of gas phase, this is just an approximation (K)

        Real rho_AP_solid = 1950; // kg/m^3 https://en.wikipedia.org/wiki/Ammonium_perchlorate
        Real rho_HTPB_solid = 920; // kg/m^3 https://www.researchgate.net/publication/279252447_Pocket_Model_for_Aluminum_Agglomeration_Based_on_Propellant_Microstructure

        Set::Patch<Set::Scalar> solidrho  = Hydro::solid.density_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> solidM    = Hydro::solid.momentum_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> m0        = Hydro::m0_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> u0_patch  = Hydro::u0_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {   
            Set::Scalar eta_hydro = 1 - eta(i,j,k) * eta(i,j,k);
            Set::Scalar etaold_hydro = 1 - etaold(i,j,k) * etaold(i,j,k);
            Set::Scalar phi = Numeric::Interpolate::NodeToCellAverage(phi_patch, i, j, k, 0);
            Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
            Set::Vector grad_eta_hydro = -2.0 * eta(i,j,k) * grad_eta;
            Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>();
            Set::Vector N = grad_eta_hydro / (grad_eta_mag + small); // Example of finding the normal vector
            pressure(i,j,k) = pressure(i,j,k) + Pref; // Scale by the reference pressure b/c ideal gas law requires absolute pressure
            rho_AP_gas(i,j,k) = pressure(i,j,k)*M_AP/(R*temp_gas); // Density of AP gaseous products assuming ideal gas
            rho_HTPB_gas(i,j,k) = pressure(i,j,k)*M_HTPB/(R*temp_gas); // Density of HTPB gaseous products assuming ideal gas
            rho_tot_gas(i,j,k) = rho_AP_gas(i,j,k)*phi + rho_HTPB_gas(i,j,k)*(1.0 - phi); // Find the average density of the fluid based on the solid species

            m0(i,j,k) = hydro.rho_ap*phi + hydro.rho_htpb*(1.0 - phi); // example of setting value to m0
            solidrho(i,j,k) = m0(i,j,k);
            
            // Set::Vector u0( hydro.u0_ap*phi + hydro.u0_htpb*(1.0 - phi), 0.0);
            // u0_patch(i,j,k,0) = u0(0); // Take the zeroth component of u0 (x and y) componets of velocity
            // u0_patch(i,j,k,1) = u0(1);

                Set::Vector u0;
            u0(0) = hydro.u0_ap*phi + hydro.u0_htpb*(1.0 - phi);
            u0(1) = 0.0; 
            #if AMREX_SPACEDIM == 3
                u0(2) = 0.0;   // Or whatever the correct z-velocity is
            #endif

            Set::Scalar deta_dt = (eta_hydro - etaold_hydro)/(dt); // time derivate approximation of eta

            if (Hydro::prescribedflowmode == Hydro::PrescribedFlowMode::Relative)
            {
                #if AMREX_SPACEDIM == 2
                    Set::Vector T(N(1), -N(0));
                    u0 = N * u0(0) + T * u0(1);
                #endif

            #if AMREX_SPACEDIM == 3
                Set::Vector T;
                T(0) = N(1);
                T(1) = -N(0);
                T(2) = 0;
                u0 = N*u0(0) + T * u0(1);
                // Might not be physcially accurate, need to find how to extend to 3 dimensions
            #endif
            }

            solidM(i,j,k,0) = solidrho(i,j,k)*u0(0);
            solidM(i,j,k,1) = solidrho(i,j,k)*u0(1);

        #if AMREX_SPACEDIM == 3
        solidM(i,j,k,2) = solidrho(i,j,k)*u0(2);
        #endif

            u0_patch(i,j,k,0) = deta_dt*(rho_AP_solid*phi + rho_HTPB_solid*(1-phi))*N(0)*rho_tot_gas(i,j,k)*velocity_mult; // Update the velocity source term based on conservation of mass
            u0_patch(i,j,k,1) = deta_dt*(rho_AP_solid*phi + rho_HTPB_solid*(1-phi))*N(1)*rho_tot_gas(i,j,k)*velocity_mult;
        #if AMREX_SPACEDIM == 3
        u0_patch(i,j,k,2) = deta_dt*(rho_AP_solid*phi + rho_HTPB_solid*(1-phi))*N(2)*rho_tot_gas(i,j,k)*velocity_mult;
        #endif
        });
    }

    Util::RealFillBoundary(*solid.density_mf[lev],geom[lev]);
    Util::RealFillBoundary(*solid.momentum_mf[lev],geom[lev]);
    Util::RealFillBoundary(*m0_mf[lev],geom[lev]);
    Util::RealFillBoundary(*u0_mf[lev],geom[lev]);
}

void Flame::TimeStepBegin(Set::Scalar a_time, int a_iter)
{
    BL_PROFILE("Integrator::Flame::TimeStepBegin");
    Base::Mechanics<model_type>::TimeStepBegin(a_time, a_iter);
    if (hydro.on && a_time >= hydro.tstart) Hydro::TimeStepBegin(a_time, a_iter);
    if (thermal.on) {
        for (int lev = 0; lev <= finest_level; ++lev)
            ic_laser->Initialize(lev, laser_mf, a_time);
    }
}

void Flame::TimeStepComplete(Set::Scalar a_time, int a_iter)
{
    BL_PROFILE("Integrator::Flame::TimeStepComplete");
    if (hydro.on && a_time >= hydro.tstart) Hydro::TimeStepComplete(a_time,a_iter);
    if (variable_pressure) {
        //Set::Scalar x_len = geom[0].ProbDomain().length(0);
        //Set::Scalar y_len = geom[0].ProbDomain().length(1);
        // Set::Scalar domain_area = x_len * y_len;
        Util::Message(INFO, "Mass = ", chamber.massflux);
        Util::Message(INFO, "Pressure = ", chamber.pressure);
    }
}

void Flame::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrador::Flame::Advance");
    Base::Mechanics<model_type>::Advance(lev, time, dt);
    if (hydro.on && time >= hydro.tstart) 
    {
        Hydro::Advance(lev,time,dt);
    }
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

            // 
            // EVOLVE PHASE FIELD (ETA)
            // 

            Set::Scalar eta_lap = Numeric::Laplacian(eta, i, j, k, 0, DX);
            Set::Scalar df_deta = ((pf.lambda / pf.eps) * dw(eta(i, j, k)) - pf.eps * pf.kappa * eta_lap);
            etanew(i, j, k) = eta(i, j, k) - L * dt * df_deta;
            if (etanew(i, j, k) <= small) etanew(i, j, k) = small;

            if (thermal.on)
            {
                //
                // Calculate thermal diffisivity and store for later gradient
                //

                alpha(i, j, k) = K / rho / cp; 

                //
                // CALCULATE MASS FLUX BASED ON EVOLVING ETA
                //
            
                mdot(i, j, k) = rho * fabs(eta(i, j, k) - etanew(i, j, k)) / dt; 

                //
                // CALCULATE HEAT FLUX BASED ON THE CALCULATED MASS FLUX
                //

                Set::Scalar q0 = propellant.get_qdot(mdot(i,j,k), phi_avg);
                heatflux(i,j,k) = ( thermal.hc*q0 + laser(i,j,k) ) / K;

            }

        });

    } // MFi For loop 


    //
    // THERMAL TRANSPORT
    // 
    if (thermal.on)
    {
        std::swap(temp_old_mf[lev], temp_mf[lev]);

        for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();

            Set::Patch<Set::Scalar>       tempnew = temp_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> temp = temp_old_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> alpha = alpha_mf.Patch(lev,mfi);

            Set::Patch<Set::Scalar>       temps = temps_mf.Patch(lev,mfi);


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
                temps(i, j, k) = temps(i, j, k) + dt * Tsolid;
                tempnew(i, j, k) = etanew(i, j, k) * temps(i, j, k) + (1.0 - etanew(i, j, k)) * thermal.Tfluid;
            });
        }
    }
 
} //Function


void Flame::TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar time, int ngrow)
{
    BL_PROFILE("Integrator::Flame::TagCellsForRefinement");
    Base::Mechanics<model_type>::TagCellsForRefinement(lev, a_tags, time, ngrow);
    if (hydro.on && time >= hydro.tstart) Hydro::TagCellsForRefinement(lev, a_tags, time, ngrow);

    const Set::Scalar* DX = geom[lev].CellSize();
    Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

    // Eta criterion for refinement
    for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tags = a_tags.array(mfi);
        Set::Patch<const Set::Scalar> eta = eta_mf.Patch(lev,mfi);

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
            Set::Patch<const Set::Scalar> phi = phi_mf.Patch(lev,mfi);

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
            Set::Patch<const Set::Scalar> temp = temp_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> eta  = eta_mf.Patch(lev,mfi);
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

void Flame::Integrate(int amrlev, Set::Scalar time, int step,
    const amrex::MFIter& mfi, const amrex::Box& box)
{
    BL_PROFILE("Flame::Integrate");
    
    Base::Mechanics<model_type>::Integrate(amrlev,time,timestep,mfi,box);
    //if (hydro.on && time >= hydro.tstart) Hydro::Integrate(amrlev, time, step, mfi, box);

    const Set::Scalar* DX = geom[amrlev].CellSize();
    Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
    Set::Patch<const Set::Scalar> eta  = eta_mf.Patch(amrlev,mfi);
    Set::Patch<const Set::Scalar> mdot = mdot_mf.Patch(amrlev,mfi);
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



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
#include "Model/Chamber/Ballistic.H"
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
    // Number of pure Allen-Cahn steps run at initialization to relax BMP sharp interface to equilibrium tanh profile
    pp.query_default("pf.relax_steps", value.pf.relax_steps, 0);

    // Boundary conditions for phase field order params
    pp.select<BC::Constant>("pf.eta.bc", value.bc_eta, 1 ); 
    value.RegisterNewFab(value.eta_mf, value.bc_eta, 1, 2, "eta", true);
    value.RegisterNewFab(value.eta_old_mf, value.bc_eta, 1, 2, "eta_old", false);

    // phase field initial condition
    pp.select<IC::Laminate,IC::Constant,IC::Expression,IC::BMP,IC::PNG>("pf.eta.ic",value.ic_eta,value.geom); 



    // Select reduced order model to capture heat feedback
    pp.select<  Model::Propellant::Constant,
                Model::Propellant::PowerLaw, 
                Model::Propellant::FullFeedback, 
                Model::Propellant::Homogenize>
        ("propellant",value.propellant);

    if (value.propellant.get_name() == "homogenize")  value.homogenized = true;

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

        value.RegisterNewFab(value.L_mf, value.bc_temp, 1, 0, "L", value.plot_field);

        value.RegisterIntegratedVariable(&value.chamber.volume, "volume");
        value.RegisterIntegratedVariable(&value.chamber.area, "area");
        value.RegisterIntegratedVariable(&value.chamber.mdot, "mass_flux");
        value.RegisterIntegratedVariable(&value.thermo_max_temp,     "max_temp",      false);
        value.RegisterIntegratedVariable(&value.thermo_mdot_max,     "mdot_max",      false);
        value.RegisterIntegratedVariable(&value.thermo_heatflux_max, "heatflux_max",  false);
        value.RegisterIntegratedVariable(&value.thermo_L_max,        "L_max",         false);
        value.RegisterIntegratedVariable(&value.thermo_eta_min,      "eta_min",       false);

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

    if (value.variable_pressure)
    {
        pp.queryclass<Model::Chamber::Ballistic>("chamber.ballistic",value.chamber.model);
        
        value.RegisterIntegratedVariable(&value.chamber.pressure, "chamber_pressure", false);
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

    if (pf.relax_steps > 0)
    {
        const Set::Scalar* DX = geom[lev].CellSize();
        Set::Scalar dx_min = DX[0];
        for (int d = 1; d < AMREX_SPACEDIM; d++) dx_min = std::min(dx_min, DX[d]);
        const Set::Scalar L_relax = 1.0;
        const Set::Scalar dt_relax = 0.5 * dx_min * dx_min / (4.0 * L_relax * pf.eps * pf.kappa);

        Util::Message(INFO, "pf.relax_steps=", pf.relax_steps, " lev=", lev, " dt_relax=", dt_relax);

        for (int s = 0; s < pf.relax_steps; s++)
        {
            // Use plain FillBoundary (not bc_eta) because bc_eta geometry is not yet defined
            // at the time Initialize(lev) is called — bc_eta->define(geom[lev]) happens AFTER
            // Initialize returns. Using bc_eta->FillBoundary here would apply the Dirichlet BC
            // using a wrong/empty geometry, overwriting ALL valid cells with η=1.
            eta_mf[lev]->FillBoundary(geom[lev].periodicity());
            std::swap(eta_old_mf[lev], eta_mf[lev]);

            for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                Set::Patch<Set::Scalar>       etanew = eta_mf.Patch(lev, mfi);
                Set::Patch<const Set::Scalar> eta    = eta_old_mf.Patch(lev, mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar eta_val = eta(i, j, k);
                    Set::Scalar eta_lap = Numeric::Laplacian(eta, i, j, k, 0, DX);
                    // Use symmetric double-well w_sym = eta^2*(1-eta)^2 for pre-relax.
                    // The simulation's w(eta) has w(0)=0 != w(1)=1 (asymmetric), which would
                    // drive the grain (eta=1) to be consumed by the gas (eta=0) during pre-relax.
                    // w_sym has equal minima at eta=0 and eta=1, so neither phase is preferentially grown.
                    Set::Scalar dw_sym = 2.0 * eta_val * (1.0 - eta_val) * (1.0 - 2.0 * eta_val);
                    Set::Scalar df_deta = (pf.lambda / pf.eps) * dw_sym - pf.eps * pf.kappa * eta_lap;
                    etanew(i, j, k) = eta_val - L_relax * dt_relax * df_deta;
                    if (etanew(i, j, k) < small) etanew(i, j, k) = small;
                    if (etanew(i, j, k) > 1.0)   etanew(i, j, k) = 1.0;
                });
            }
        }
        // Sync eta_old to the relaxed eta so both start consistent
        amrex::MultiFab::Copy(*eta_old_mf[lev], *eta_mf[lev], 0, 0, 1, eta_mf[lev]->nGrow());
    }

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
        L_mf[lev]->setVal(0.0);
        ic_laser->Initialize(lev, laser_mf);
    }
    //if (variable_pressure) chamber.pressure = 1.0;
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

            if (elastic.on)
            {
                Set::Patch <const Set::Scalar> temp = temp_mf.Patch(lev,mfi);
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
                    model_ap.F0 *= (temp_avg - elastic.Telastic);
                    model_ap.F0 += Set::Matrix::Identity();
                    model_type model_htpb = elastic.model_htpb;
                    model_htpb.F0 -= Set::Matrix::Identity();
                    model_htpb.F0 *= (temp_avg - elastic.Telastic);
                    model_htpb.F0 += Set::Matrix::Identity();

                    if (homogenized)
                    {
                        model_type model_void = elastic.model_void;
                        model_void.F0 -= Set::Matrix::Identity();
                        model_void.F0 *= (temp_avg - elastic.Telastic);
                        model_void.F0 += Set::Matrix::Identity();
                        model(i, j, k) = (model_ap * massfraction + model_htpb * (1. - massfraction)) * phi_avg + model_void * (1. - phi_avg);
                    }
                    else
                    {
                        model(i, j, k) = model_ap * phi_avg + model_htpb * (1. - phi_avg);
                    }
                        
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

void Flame::TimeStepBegin(Set::Scalar a_time, int a_iter)
{
    BL_PROFILE("Integrator::Flame::TimeStepBegin");
    Base::Mechanics<model_type>::TimeStepBegin(a_time, a_iter);
    if (thermal.on) {
        for (int lev = 0; lev <= finest_level; ++lev)
            ic_laser->Initialize(lev, laser_mf, a_time);
    }
}

void Flame::TimeStepComplete(Set::Scalar a_time, int /*a_iter*/)
{
    BL_PROFILE("Integrator::Flame::TimeStepComplete");

    Util::Message(INFO, "MONITOR thermal.on=", thermal.on);
    if (thermal.on)
    {
        thermo_max_temp     = 0.0;
        thermo_mdot_max     = 0.0;
        thermo_heatflux_max = 0.0;
        thermo_L_max        = 0.0;
        thermo_eta_min      = 1.0;
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            Set::Scalar max_temp = temp_mf[lev]->max(0);
            Set::Scalar min_temp = temp_mf[lev]->min(0);
            Set::Scalar max_mdot = mdot_mf[lev]->max(0);
            Set::Scalar max_heatflux = heatflux_mf[lev]->max(0);
            Set::Scalar max_L = L_mf[lev]->max(0);
            Set::Scalar min_eta = eta_mf[lev]->min(0);
            thermo_max_temp     = std::max(thermo_max_temp,     max_temp);
            thermo_mdot_max     = std::max(thermo_mdot_max,     max_mdot);
            thermo_heatflux_max = std::max(thermo_heatflux_max, max_heatflux);
            thermo_L_max        = std::max(thermo_L_max,        max_L);
            thermo_eta_min      = std::min(thermo_eta_min,      min_eta);
            Util::Message(INFO, "t=", std::setw(10), std::setprecision(4), std::fixed, a_time,
            " lev=", std::setw(2), lev,
            " T=[", std::setw(8), std::setprecision(2), Unit::Temperature(min_temp),
            ",", std::setw(8), std::setprecision(2), Unit::Temperature(max_temp), "]",
            " mdot_max=", std::setw(10), std::setprecision(4), max_mdot,
            " heatflux_max=", std::setw(10), std::setprecision(4), max_heatflux,
            " L_max=", std::setw(10), std::setprecision(4), max_L,
            " eta_min=", std::setw(10), std::setprecision(4), min_eta);
        }
    }

    if (variable_pressure)
    {
        auto [new_pressure, current_dpdt] = chamber.model.Advance(timestep, chamber.mdot, chamber.volume, chamber.pressure);
        chamber.pressure = new_pressure;
        Util::Message(INFO, "chamber.pressure = ", Unit::Pressure(chamber.pressure));
        Util::Message(INFO, "chamber.mdot = ", Unit::Mass(chamber.mdot) / Unit::Time());
        Util::Message(INFO, "chamber.volume = ", Unit::Volume(chamber.volume));
        Util::Message(INFO, "chamber.dpdt = ", current_dpdt);
    }
}

void Flame::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrador::Flame::Advance");
    Base::Mechanics<model_type>::Advance(lev, time, dt);
    const Set::Scalar* DX = geom[lev].CellSize();

    std::swap(eta_old_mf[lev], eta_mf[lev]);

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
        Set::Patch<Set::Scalar> L_out = L_mf.Patch(lev, mfi);


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
            L_out(i, j, k) = L;
            if (L != L) {
                Util::Message(INFO, "phi = ", phi_avg);
                Util::Message(INFO, "T = ", T);
                Util::Message(INFO, "L = ", L);
                Util::Abort(INFO,"NaN detected in mobility");
            }

            // 
            // EVOLVE PHASE FIELD (ETA)
            // 

            Set::Scalar eta_lap = Numeric::Laplacian(eta, i, j, k, 0, DX);
            Set::Scalar df_deta = ((pf.lambda / pf.eps) * dw(eta(i, j, k)) - pf.eps * pf.kappa * eta_lap);
            etanew(i, j, k) = eta(i, j, k) - L * dt * df_deta;
            if (etanew(i, j, k) <= small) etanew(i, j, k) = small;

            if (etanew(i, j, k) != etanew(i, j, k)) { // Check for NaN
                Util::Message(INFO, "L = ", L);
                Util::Message(INFO, "df_deta = ", df_deta);
                Util::Message(INFO, "eta = ", eta(i,j,k));
                Util::Abort(INFO,"NaN detected in etanew");
                }


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
                if (heatflux(i,j,k) != heatflux(i,j,k)) { // Check for NaN
                    Util::Message(INFO, "hc = ", thermal.hc);
                    Util::Message(INFO, "q0 = ", q0);
                    Util::Message(INFO, "K = ", K);
                    Util::Message(INFO, "laser = ", laser(i,j,k));
                    Util::Abort(INFO,"NaN detected in heatflux");
                }

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

                if (dTdt != dTdt) { // Check for NaN
                    Util::Message(INFO, "grad_eta = ", grad_eta);
                    Util::Message(INFO, "grad_temp = ", grad_temp);
                    Util::Message(INFO, "lap_temp = ", lap_temp);
                    Util::Message(INFO, "grad_alpha = ", grad_alpha);
                    Util::Message(INFO, "alpha = ", alpha(i,j,k));
                    Util::Message(INFO, "heatflux = ", heatflux(i,j,k));
                    Util::Message(INFO, "eta = ", eta(i,j,k));
                    Util::Abort(INFO,"NaN detected in dTdt");
                }


                Set::Scalar Tsolid = dTdt + temps(i, j, k) * (etanew(i, j, k) - eta(i, j, k)) / dt;
                temps(i, j, k) = temps(i, j, k) + dt * Tsolid;
                tempnew(i, j, k) = etanew(i, j, k) * temps(i, j, k) + (1.0 - etanew(i, j, k)) * thermal.Tfluid;

                if (tempnew(i, j, k) != tempnew(i, j, k)) { // Check for NaN
                    Util::Message(INFO, "dTdt = ", dTdt);
                    Util::Message(INFO, "Tsolid = ", Tsolid);
                    Util::Message(INFO, "temps = ", temps(i,j,k));
                    Util::Message(INFO, "temp = ", temp(i,j,k));
                    Util::Abort(INFO,"NaN detected in tempnew");
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

void Flame::Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/,
    const amrex::MFIter& mfi, const amrex::Box& box)
{
    BL_PROFILE("Flame::Integrate");
    const Set::Scalar* DX = geom[amrlev].CellSize();
    Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
    Set::Patch<const Set::Scalar> eta  = eta_mf.Patch(amrlev,mfi);
    Set::Patch<const Set::Scalar> mdot = mdot_mf.Patch(amrlev,mfi);
    if (variable_pressure) {
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            chamber.volume += (1.0 - eta(i, j, k, 0)) * dv;
            Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
            Set::Scalar normgrad = grad.lpNorm<2>();
            Set::Scalar da = normgrad * dv;
            chamber.area += da;

            // Set::Vector mgrad = Numeric::Gradient(mdot, i, j, k, 0, DX);
            // Set::Scalar mnormgrad = mgrad.lpNorm<2>();
            // Set::Scalar dm = mnormgrad * dv;
            //chamber.massflux += dm;

            // new - chamber model
            chamber.mdot += mdot(i, j, k, 0) * dv;
            //chamber.volume += volume;
        });
    }
    else {
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            chamber.volume += (1.0 - eta(i, j, k, 0)) * dv;
            Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
            Set::Scalar normgrad = grad.lpNorm<2>();
            Set::Scalar da = normgrad * dv;
            chamber.area += da;

            chamber.mdot += mdot(i, j, k, 0) * dv;
        });
    }
    // time dependent pressure data from experimenta -> p = 0.0954521220950523 * exp(15.289993148880678 * t)
}
} // namespace Integrator



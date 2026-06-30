#include "IC/StarAftGrain.H"
#include "Flame.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"
#include "IC/Laminate.H"
#include "IC/Constant.H"
#include "IC/PointList.H"
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

Flame::Flame() :
    Base::Mechanics<model_type>() {}

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
    // eta carries 3 ghost cells (not 2): the elastic model blend in UpdateModel
    // does CellToNodeAverage(eta) over the model's grown box, which reaches one
    // cell past a 2-ghost buffer at interior grid edges. eta_old_mf must match
    // because Advance() swaps the current/old handles in place.
    value.RegisterNewFab(value.eta_mf, value.bc_eta, 1, 3, "eta", true);
    value.RegisterNewFab(value.eta_old_mf, value.bc_eta, 1, 3, "eta_old", false);
    value.RegisterNewFab(value.psi_mf, value.bc_eta, 1, 2, "psi", true);

    // Inital value of eta that doesn't evolve and is used during refiment to set the updated values of eta with voids in the domain.
    // Used to fix a bug where duirn refinement, a void won't be updated correctly and would be a square, not a circle
    value.RegisterNewFab(value.eta_0_mf, value.bc_eta, 1, 2, "eta_0", 0);

    // Allen-Cahn mobility L is computed and written into L_mf every Advance
    // regardless of the thermal model, so L_mf must be registered
    // unconditionally. It was previously created only inside the thermal.on
    // block, which left it null and segfaulted Advance (L_out write) when
    // thermal.on=0. nghost=0 so the BC is only nominal; reuse bc_eta.
    value.RegisterNewFab(value.L_mf, value.bc_eta, 1, 0, "L", value.plot_field);

    // phase field initial condition
    pp.select<IC::Laminate,IC::Constant,IC::Expression,IC::BMP,IC::PNG,IC::PSRead,IC::StarAftGrain>("pf.eta.ic",value.ic_eta,value.geom);


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

    // Whether to compute the pressure evolution. Query this before thermal
    // registration because variable-pressure runs need chamber geometry/mass
    // integration even when thermal transport is disabled.
    pp_query_default("variable_pressure", value.variable_pressure, false);

    pp.select_default<BC::Constant>("thermal.temp.bc", value.bc_temp, 1, Unit::Temperature());
    value.RegisterNewFab(value.temp_mf, value.bc_temp, 1, 3, "temp", value.thermal.on && value.plot_field);
    value.RegisterNewFab(value.temp_old_mf, value.bc_temp, 1, 3, "temp_old", false);
    value.RegisterNewFab(value.temps_mf, value.bc_temp, 1, 0, "temps", false);

    value.RegisterNewFab(value.mdot_mf, value.bc_temp, 1, 0, "mdot", value.thermal.on && value.plot_field);
    value.RegisterNewFab(value.alpha_mf, value.bc_temp, 1, 0, "alpha", value.thermal.on && value.plot_field);
    value.RegisterNewFab(value.heatflux_mf, value.bc_temp, 1, 0, "heatflux", value.thermal.on && value.plot_field);
    value.RegisterNewFab(value.laser_mf, value.bc_temp, 1, 0, "laser", value.thermal.on && value.plot_field);
    value.RegisterNewFab(value.thermal.has_exceeded_Tcutoff, value.bc_temp, 1, 2, "exceeded_Tcutoff", false);

    if (value.thermal.on) {

        // Used to change heat flux units
        pp_query_default("thermal.hc", value.thermal.hc, "1.0", Unit::Power()/Unit::Area());

        // Effective fluid temperature, temp of the eta = 0 (fluid) region
        pp_query_default("thermal.Tfluid", value.thermal.Tfluid, value.thermal.Tref);

        // Cutoff value for regression, if T < Tcutoff eta won't evolve/regress
        pp.query_default("thermal.Tcutoff", value.thermal.Tcutoff, "0.0", Unit::Temperature());

        // Switch time of the improved regridding where eta and the temperature field are both used. It is recommended to make this time ~10x the timestep.
        // Before this the refinement is based on the gradient of eta which helps the laser IC start correctly. A regrid is forced when this time is reached.
        pp.query_default("thermal.end_initial_refine_time", value.thermal.end_initial_refine_time, "0.0", Unit::Time());

        // Inital refinement of the phi field based on phi gradient. After time > end_initial_refine_time stops refining at these phi values.
        pp.query_default("thermal.phi_refinement_criterion_inital", value.thermal.phi_refinement_criterion_inital, 1.0e100);

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

    if (value.variable_pressure)
    {
        pp.queryclass<Model::Chamber::Ballistic>("chamber.ballistic",value.chamber.model);

        value.RegisterIntegratedVariable(&value.chamber.pressure, "chamber_pressure", false);
    }

    if (!value.thermal.on && value.variable_pressure)
    {
        value.RegisterIntegratedVariable(&value.chamber.volume, "volume");
        value.RegisterIntegratedVariable(&value.chamber.area, "area");
        value.RegisterIntegratedVariable(&value.chamber.mdot, "mass_flux");
    }

    // Refinement criterion for eta field, if thermal is on, cells will only be tagged for refinement if T>0.9*TCutoff,
    // and the gradient of eta > m_refinement_criterion at each cell
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
    pp.select_default<IC::Laminate,IC::Expression,IC::Constant,IC::BMP,IC::PNG,IC::PSRead,IC::StarAftGrain>
        ("phi.ic",value.ic_phi,value.geom);

    value.RegisterNodalFab(value.phi_mf, 1, 2, "phi", true);

    // Whether to use Neo-hookean Elastic model
    pp_query_default("elastic.on", value.elastic.on, 0);

    // Body force (surface traction applied at the burn interface; pressure units).
    // Bare (unit-less) numbers are still accepted and pass through unchanged.
    pp_query_default("elastic.traction", value.elastic.traction, "0.0", Unit::Pressure());

    // If 1, the interface traction is set each step from the evolving chamber
    // pressure; elastic.traction is then ignored. (Coupling is explicit: the
    // elastic solve at step N uses the pressure from the end of step N-1.)
    pp_query_default("elastic.traction_from_chamber", value.elastic.traction_from_chamber, 0);

    // Phi refinement criteria
    pp_query_default("elastic.phirefinement", value.elastic.phirefinement, 1);

    // Floor for the psi weighting field (0 => psi=eta, the original behavior).
    // A small floor keeps the masked elastic operator non-singular in the gas,
    // which is what lets MLMG converge on thin/complex grain geometries.
    // NOTE: must be read before the Mechanics queryclass below, which strictly
    // validates that every elastic.* key has been consumed.
    pp_query_default("elastic.psi_floor", value.elastic.psi_floor, 0.0);

    // Whether to mask the elastic operator with psi at all (1) or condition the
    // gas purely through the real model_void stiffness with no mask (0). Same
    // strict-validation ordering requirement as psi_floor above.
    pp_query_default("elastic.use_psi", value.elastic.use_psi, 1);

    pp.queryclass<Base::Mechanics<model_type>>("elastic",value);

    if (value.m_type != Type::Disable)
    {
        // Reference temperature for thermal expansion
        // (temperature at which the material is strain-free)
        pp_query_default("Telastic", value.elastic.Telastic, value.thermal.Tref);
        // Elastic model schema follows the propellant resolution:
        //  - homogenized (grain scale): one model_prop + void + casing, blended
        //    by the three-material partition of unity in UpdateModel.
        //  - resolved (mesoscale, e.g. fullfeedback): per-phase model_ap and
        //    model_htpb, blended by the species field phi (rule of mixtures).
        // This keeps both the full-grain and mesoscale AP/HTPB use-cases working
        // from one integrator (resolved was the original behavior; homogenized
        // replaced it -- now both are supported).
        if (value.homogenized)
        {
            // elastic model of the homogenized propellant (single material)
            pp.queryclass<Model::Solid::Finite::NeoHookeanPredeformed>("model_prop", value.elastic.model_prop);
            // elastic model of void (gas phase, lives where phi=1 and eta=0)
            pp.queryclass<Model::Solid::Finite::NeoHookeanPredeformed>("model_void", value.elastic.model_void);
            // elastic model of casing (stiff confinement outside the fuel disk, phi=0)
            pp.queryclass<Model::Solid::Finite::NeoHookeanPredeformed>("model_casing", value.elastic.model_casing);
        }
        else
        {
            // resolved AP/HTPB elastic models (phi=1 -> AP, phi=0 -> HTPB)
            pp.queryclass<Model::Solid::Finite::NeoHookeanPredeformed>("model_ap", value.elastic.model_ap);
            pp.queryclass<Model::Solid::Finite::NeoHookeanPredeformed>("model_htpb", value.elastic.model_htpb);
        }

        // Use (floored) eta as the psi field to weight the elastic operator and zero
        // out the gas region. psi_mf is filled from eta in UpdateModel. When use_psi=0
        // we install no psi at all: psi_avg stays 1 everywhere and the gas is
        // conditioned only by its (soft, unitized) model_void stiffness.
        // Flame owns its psi_mf member separately from Base::Mechanics::psi_mf.
        // Keep the base psi_on flag false so TimeStepBegin does not reattach the
        // empty base field; attach Flame's eta-derived psi directly to the solver.
        value.psi_on = false;
        if (value.elastic.use_psi)
            value.solver.setPsi(value.psi_mf);
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
        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar dx_min = DX(0);
        for (int d = 1; d < AMREX_SPACEDIM; d++) dx_min = std::min(dx_min, DX(d));
        const Set::Scalar L_relax = 1.0;
        const Set::Scalar dt_relax = 0.5 * dx_min * dx_min / (4.0 * L_relax * pf.eps * pf.kappa);

        Util::Message(INFO, "pf.relax_steps=", pf.relax_steps, " lev=", lev, " dt_relax=", dt_relax);

        // GPU: copy POD members into locals so the relax kernel captures by value.
        auto pf = this->pf;
        const Set::Scalar small = this->small;

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
                    Set::Scalar eta_lap = Numeric::Laplacian(eta, i, j, k, 0, DX.data());
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
    else
    {
        temp_mf[lev]->setVal(thermal.Tref);
        temp_old_mf[lev]->setVal(thermal.Tref);
        temps_mf[lev]->setVal(thermal.Tref);
        alpha_mf[lev]->setVal(0.0);
        mdot_mf[lev]->setVal(0.0);
        heatflux_mf[lev]->setVal(0.0);
        L_mf[lev]->setVal(0.0);
        laser_mf[lev]->setVal(0.0);
        thermal.has_exceeded_Tcutoff[lev]->setVal(0.0);
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
        Set::Vector DX(geom[lev].CellSize());

        phi_mf[lev]->FillBoundary();
        eta_mf[lev]->FillBoundary();
        temp_mf[lev]->FillBoundary();

        // Build the (floored) psi weighting field used by the elastic solver:
        //   psi = psi_floor + (1 - psi_floor) * eta
        // With psi_floor=0 this is just eta (original behavior). A small floor keeps
        // the masked operator non-singular in the gas region. eta's ghosts are valid
        // (FillBoundary above), so copy with ghosts to keep psi consistent.
        {
            const int ng = psi_mf[lev]->nGrow();
            amrex::MultiFab::Copy(*psi_mf[lev], *eta_mf[lev], 0, 0, 1, ng);
            psi_mf[lev]->mult(1.0 - elastic.psi_floor, 0, 1, ng);
            psi_mf[lev]->plus(elastic.psi_floor, 0, 1, ng);
        }

        for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box smallbox = mfi.nodaltilebox();
            amrex::Box bx = mfi.grownnodaltilebox() & domain;
            Set::Patch<model_type>        model = model_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> phi   = phi_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> eta   = eta_mf.Patch(lev,mfi);
            Set::Patch<Set::Vector>       rhs   = rhs_mf.Patch(lev,mfi);

            // GPU: copy members into locals so the model-build kernels below
            // capture by value instead of the host `this` pointer.
            auto       elastic     = this->elastic;
            const bool homogenized = this->homogenized;

            if (elastic.on)
            {
                Set::Patch <const Set::Scalar> temp = temp_mf.Patch(lev,mfi);
                const Set::Scalar traction = elastic.traction_from_chamber ? chamber.pressure : elastic.traction;
                amrex::ParallelFor(smallbox, [=] AMREX_GPU_DEVICE(int i, int j, int k)

                {
                    Set::Vector grad_eta = Numeric::CellGradientOnNode(eta, i, j, k, 0, DX.data());
                    // The solver enforces div(sigma) = rhs, so rhs = -b_phys (ALAMO sign
                    // convention). grad_eta points from gas (eta=0) into the solid (eta=1),
                    // i.e. radially inward at a burning surface. A chamber pressure must
                    // COMPRESS the grain, so the physical body force is b_phys = +traction*grad_eta
                    // (inward) and therefore rhs = -traction*grad_eta. The previous sign put
                    // the grain in tension (verified by the static patch test).
                    rhs(i, j, k) = -traction * grad_eta;
                });
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar phi_avg = phi(i, j, k, 0);
                    Set::Scalar temp_avg = Numeric::Interpolate::CellToNodeAverage(temp, i, j, k, 0);

                    if (homogenized)
                    {
                        // Single homogenized propellant model. Apply the thermoelastic
                        // eigenstrain F0 <- I + (F0 - I)*(T - Telastic).
                        model_type model_prop = elastic.model_prop;
                        model_prop.F0 -= Set::Matrix::Identity();
                        model_prop.F0 *= (temp_avg - elastic.Telastic);
                        model_prop.F0 += Set::Matrix::Identity();

                        model_type model_void = elastic.model_void;
                        model_void.F0 -= Set::Matrix::Identity();
                        model_void.F0 *= (temp_avg - elastic.Telastic);
                        model_void.F0 += Set::Matrix::Identity();
                        model_type model_casing = elastic.model_casing;
                        model_casing.F0 -= Set::Matrix::Identity();
                        model_casing.F0 *= (temp_avg - elastic.Telastic);
                        model_casing.F0 += Set::Matrix::Identity();

                        // Three-material partition of unity (sums to 1):
                        //   solid  = phi*eta       -- fuel present AND unburned (propellant)
                        //   void   = phi*(1-eta)   -- fuel region, burned (combustion gas, very soft)
                        //   casing = (1-phi)       -- outside the fuel disk (stiff confinement)
                        // phi (static) masks the cylindrical chamber out of the square domain;
                        // eta (live) splits propellant from gas inside it.
                        // NOTE: this CellToNodeAverage(eta) over the grown box requires eta to carry
                        // 3 ghost cells (like temp) so the read stays in bounds at grid edges.
                        Set::Scalar eta_avg = Numeric::Interpolate::CellToNodeAverage(eta, i, j, k, 0);
                        Set::Scalar w_solid  = phi_avg * eta_avg;
                        Set::Scalar w_void   = phi_avg * (1. - eta_avg);
                        Set::Scalar w_casing = 1. - phi_avg;
                        model(i, j, k) = model_prop * w_solid
                                        + model_void * w_void
                                        + model_casing * w_casing;
                    }
                    else
                    {
                        // Resolved AP/HTPB mesoscale model: blend the two per-phase
                        // elastic models by the species field phi (rule of mixtures),
                        // phi=1 -> AP, phi=0 -> HTPB. Each phase carries the same
                        // thermoelastic eigenstrain F0 <- I + (F0 - I)*(T - Telastic).
                        model_type model_ap = elastic.model_ap;
                        model_ap.F0 -= Set::Matrix::Identity();
                        model_ap.F0 *= (temp_avg - elastic.Telastic);
                        model_ap.F0 += Set::Matrix::Identity();
                        model_type model_htpb = elastic.model_htpb;
                        model_htpb.F0 -= Set::Matrix::Identity();
                        model_htpb.F0 *= (temp_avg - elastic.Telastic);
                        model_htpb.F0 += Set::Matrix::Identity();
                        model(i, j, k) = model_ap * phi_avg
                                        + model_htpb * (1. - phi_avg);
                    }
                });
            }
            else
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    // Elasticity disabled: build the propellant model with no
                    // eigenstrain (F0 = 0), mode-aware so the schema that was
                    // actually parsed is honored (homogenized model_prop vs the
                    // resolved model_ap/model_htpb blend). Not used for forcing.
                    Set::Scalar phi_avg = phi(i, j, k, 0);
                    model_type m = homogenized
                                ? elastic.model_prop
                                : elastic.model_ap * phi_avg + elastic.model_htpb * (1. - phi_avg);
                    for (int ii = 0; ii < AMREX_SPACEDIM; ++ii)
                    {
                        for (int jj = 0; jj < AMREX_SPACEDIM; ++jj)
                        {
                            m.F0(ii, jj) = 0.0;
                        }
                    }
                    model(i, j, k) = m;
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

    if (a_time > thermal.end_initial_refine_time)
    {
        if (!end_initial_refine) {
            for (int lev = 0; lev <= finest_level; ++lev)
                Flame::Regrid(lev, a_time);
            end_initial_refine = 1;
        }

        prev_finest_ba = grids[finest_level];
        prev_finest_level = finest_level;
    }
}

void Flame::TimeStepComplete(Set::Scalar /*a_time*/, int /*a_iter*/)
{
    BL_PROFILE("Integrator::Flame::TimeStepComplete");

    if (thermal.on)
    {
        // Chamber thermo diagnostics written to thermo.dat. The five min/max
        // reductions over every level are accumulated into a single fused device
        // pass (one ReduceData, one host sync) instead of separate per-level
        // MultiFab::min/max calls that each launch a kernel and synchronize. Each
        // rank reduces only its local boxes; the MPI all-reduce below then matches
        // the global result MultiFab::min/max produced. The per-layer debug print
        // has been removed. The 0.0/1.0 floors are preserved so the output stays
        // bit-identical to the previous code (max/min are order-independent).
        amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMax, amrex::ReduceOpMax,
                        amrex::ReduceOpMax, amrex::ReduceOpMin> reduce_op;
        amrex::ReduceData<Set::Scalar, Set::Scalar, Set::Scalar,
                        Set::Scalar, Set::Scalar> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            for (amrex::MFIter mfi(*temp_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& box = mfi.tilebox();
                Set::Patch<const Set::Scalar> temp     = temp_mf.Patch(lev, mfi);
                Set::Patch<const Set::Scalar> mdot     = mdot_mf.Patch(lev, mfi);
                Set::Patch<const Set::Scalar> heatflux = heatflux_mf.Patch(lev, mfi);
                Set::Patch<const Set::Scalar> L        = L_mf.Patch(lev, mfi);
                Set::Patch<const Set::Scalar> eta      = eta_mf.Patch(lev, mfi);
                reduce_op.eval(box, reduce_data,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                    {
                        return { temp(i, j, k, 0), mdot(i, j, k, 0), heatflux(i, j, k, 0),
                                L(i, j, k, 0), eta(i, j, k, 0) };
                    });
            }
        }

        ReduceTuple hv = reduce_data.value(reduce_op);
        thermo_max_temp     = std::max(Set::Scalar(0.0), amrex::get<0>(hv));
        thermo_mdot_max     = std::max(Set::Scalar(0.0), amrex::get<1>(hv));
        thermo_heatflux_max = std::max(Set::Scalar(0.0), amrex::get<2>(hv));
        thermo_L_max        = std::max(Set::Scalar(0.0), amrex::get<3>(hv));
        thermo_eta_min      = std::min(Set::Scalar(1.0), amrex::get<4>(hv));

        // MultiFab::min/max all-reduce internally; replicate that so the value is
        // correct under MPI (the CPU correctness baseline runs np8, GPU runs np1).
        amrex::ParallelDescriptor::ReduceRealMax(thermo_max_temp);
        amrex::ParallelDescriptor::ReduceRealMax(thermo_mdot_max);
        amrex::ParallelDescriptor::ReduceRealMax(thermo_heatflux_max);
        amrex::ParallelDescriptor::ReduceRealMax(thermo_L_max);
        amrex::ParallelDescriptor::ReduceRealMin(thermo_eta_min);
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
    Set::Vector DX(geom[lev].CellSize());

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

    // GPU: device lambdas may not capture the host `this` pointer. Copy the POD
    // phase-field params, the propellant model, and the `small` floor into locals
    // (shadowing the members), and pull the thermal scalars used inside kernels
    // into named locals so the kernels below capture by value, not via `this`.
    auto              propellant      = this->propellant;
    auto              pf              = this->pf;
    const Set::Scalar small           = this->small;
    const bool        thermal_on      = thermal.on;
    const Set::Scalar thermal_hc      = thermal.hc;
    const Set::Scalar thermal_Tcutoff = thermal.Tcutoff;
    const Set::Scalar thermal_Tfluid  = thermal.Tfluid;
    const bool        variable_pressure_on = this->variable_pressure;

    Util::DeviceErrorFlag advance_error;
    int* advance_error_flag = advance_error.dataPtr();

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

        Set::Patch<Set::Scalar> exceeded_Tcutoff = thermal.has_exceeded_Tcutoff.Patch(lev, mfi);
        Set::Scalar Tcutoff = thermal.Tcutoff;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            //
            // CALCULATE PHI-AVERAGED QUANTITIES
            //
            Set::Scalar phi_avg = Numeric::Interpolate::NodeToCellAverage(phi, i, j, k, 0);
            Set::Scalar T = thermal_on ? temp(i,j,k) : NAN;

            Set::Scalar K = propellant.get_K(phi_avg);

            Set::Scalar rho = propellant.get_rho(phi_avg);

            Set::Scalar cp = propellant.get_cp(phi_avg);

            //
            // CALCULATE MOBILITY
            //
            Set::Scalar L = propellant.get_L(  phi_avg, T);
            L_out(i, j, k) = L;
            // L (mobility) is always used by the eta evolution, so validate it
            // unconditionally. K/rho/cp are thermal quantities that are
            // legitimately NAN for burn-rate-only propellant models (e.g.
            // PowerLaw, whose get_K/get_rho/get_cp return NAN) and are only
            // consumed inside the thermal_on block below -- validating them
            // unconditionally spuriously aborts a thermal.on=0 run.
            if (std::isnan(L) || std::isinf(L))
            {
                Util::SetDeviceError(advance_error_flag);
            }
            if ((thermal_on || variable_pressure_on) &&
                (std::isnan(rho) || std::isinf(rho)))
            {
                Util::SetDeviceError(advance_error_flag);
            }
            if (thermal_on &&
                (std::isnan(K) || std::isinf(K) ||
                std::isnan(cp) || std::isinf(cp)))
            {
                Util::SetDeviceError(advance_error_flag);
            }

            //
            // EVOLVE PHASE FIELD (ETA)
            //

            Set::Scalar eta_lap = Numeric::Laplacian(eta, i, j, k, 0, DX.data());
            Set::Scalar df_deta = ((pf.lambda / pf.eps) * dw(eta(i, j, k)) - pf.eps * pf.kappa * eta_lap);

            if (df_deta < 0) {
                // Prevent eta from increasing/healing. A bug was found where if the diffuse thickness was too large compared to a void
                // (region of eta = 0), eta would heal/increase in a non-physcial way, this statement stops that behavior
                df_deta = 0.0;
            }
            if (thermal_on && T < thermal_Tcutoff) {
                // If the temperature is lower then the cutoff temperature don't evolve the eta field
                df_deta = 0.0;
            }
            etanew(i, j, k) = eta(i, j, k) - L * dt * df_deta;
            if (etanew(i, j, k) > eta(i, j, k)) etanew(i, j, k) = eta(i, j, k);
            if (etanew(i, j, k) <= small) etanew(i, j, k) = small;
            if (std::isnan(etanew(i, j, k)) || std::isinf(etanew(i, j, k)))
            {
                Util::SetDeviceError(advance_error_flag);
            }

            if (thermal_on || variable_pressure_on)
            {
                mdot(i, j, k) = rho * fabs(eta(i, j, k) - etanew(i, j, k)) / dt;
                if (std::isnan(mdot(i, j, k)) || std::isinf(mdot(i, j, k)))
                {
                    Util::SetDeviceError(advance_error_flag);
                }
            }
            else
            {
                mdot(i, j, k) = 0.0;
            }

            if (thermal_on)
            {
                //
                // Calculate thermal diffisivity and store for later gradient
                //

                alpha(i, j, k) = K / rho / cp;
                if (std::isnan(alpha(i, j, k)) || std::isinf(alpha(i, j, k)))
                {
                    Util::SetDeviceError(advance_error_flag);
                }

                //
                // CALCULATE HEAT FLUX BASED ON THE CALCULATED MASS FLUX
                //

                Set::Scalar q0 = propellant.get_qdot(mdot(i,j,k), phi_avg);
                heatflux(i,j,k) = ( thermal_hc*q0 + laser(i,j,k) ) / K;
                if (std::isnan(q0) || std::isinf(q0) ||
                    std::isnan(heatflux(i, j, k)) || std::isinf(heatflux(i, j, k)))
                {
                    Util::SetDeviceError(advance_error_flag);
                }

                if (temp(i,j,k) > Tcutoff)
                {
                    exceeded_Tcutoff(i,j,k) = 1;
                }

            }

        });

    } // MFi For loop
    Util::AbortIfDeviceError(advance_error, INFO,
        "non-finite value detected in Flame::Advance phase-field kernel at lev=", lev);


    //
    // THERMAL TRANSPORT
    //
    if (thermal.on)
    {
        std::swap(temp_old_mf[lev], temp_mf[lev]);

        Util::DeviceErrorFlag thermal_error;
        int* thermal_error_flag = thermal_error.dataPtr();

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
                Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX.data());
                Set::Vector grad_temp = Numeric::Gradient(temp, i, j, k, 0, DX.data());
                Set::Scalar lap_temp = Numeric::Laplacian(temp, i, j, k, 0, DX.data());
                Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>();
                Set::Vector grad_alpha = Numeric::Gradient(alpha, i, j, k, 0, DX.data(), sten);
                Set::Scalar dTdt = 0.0;
                dTdt += grad_eta.dot(grad_temp * alpha(i, j, k));
                dTdt += grad_alpha.dot(eta(i, j, k) * grad_temp);
                dTdt += eta(i, j, k) * alpha(i, j, k) * lap_temp;
                dTdt += alpha(i, j, k) * heatflux(i, j, k) * grad_eta_mag;

                Set::Scalar Tsolid = dTdt + temps(i, j, k) * (etanew(i, j, k) - eta(i, j, k)) / dt;
                temps(i, j, k) = temps(i, j, k) + dt * Tsolid;
                tempnew(i, j, k) = etanew(i, j, k) * temps(i, j, k) + (1.0 - etanew(i, j, k)) * thermal_Tfluid;
                if (std::isnan(dTdt) || std::isinf(dTdt) ||
                    std::isnan(temps(i, j, k)) || std::isinf(temps(i, j, k)) ||
                    std::isnan(tempnew(i, j, k)) || std::isinf(tempnew(i, j, k)))
                {
                    Util::SetDeviceError(thermal_error_flag);
                }

            });
        }
        Util::AbortIfDeviceError(thermal_error, INFO,
            "non-finite value detected in Flame::Advance thermal kernel at lev=", lev);
    }

} //Function


void Flame::TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar time, int ngrow)
{
    BL_PROFILE("Integrator::Flame::TagCellsForRefinement");
    Base::Mechanics<model_type>::TagCellsForRefinement(lev, a_tags, time, ngrow);

    Set::Vector DX(geom[lev].CellSize());
    Set::Scalar dr = sqrt(AMREX_D_TERM(DX(0) * DX(0), +DX(1) * DX(1), +DX(2) * DX(2)));

    // GPU: pull refinement-criterion members into locals (shadowing the members)
    // plus the thermal scalars used in kernels, so the tag kernels capture by
    // value instead of the host `this` pointer.
    const Set::Scalar m_refinement_criterion       = this->m_refinement_criterion;
    const Set::Scalar t_refinement_criterion       = this->t_refinement_criterion;
    const Set::Scalar t_refinement_restriction     = this->t_refinement_restriction;
    const Set::Scalar phi_refinement_criterion     = this->phi_refinement_criterion;
    const Set::Scalar thermal_Tcutoff              = thermal.Tcutoff;
    const Set::Scalar thermal_phi_ref_initial      = thermal.phi_refinement_criterion_inital;
    const Set::Scalar thermal_end_initial_refine_t = thermal.end_initial_refine_time;

    const bool thermal_on = thermal.on;
    const bool phi_refinement_on = elastic.phirefinement;

    if (thermal_on) {
        for (amrex::MFIter mfi(*temp_mf[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<char> const& tags = a_tags.array(mfi);
            Set::Patch<const Set::Scalar> eta  = eta_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> phi  = phi_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> temp = temp_mf.Patch(lev,mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Vector gradeta = Numeric::Gradient(eta, i, j, k, 0, DX.data());
                Set::Vector gradphi = Numeric::Gradient(phi, i, j, k, 0, DX.data());
                Set::Vector tempgrad = Numeric::Gradient(temp, i, j, k, 0, DX.data());

                bool tag = false;
                tag = tag ||
                    (gradeta.lpNorm<2>() * dr * 2 > m_refinement_criterion &&
                    eta(i, j, k) >= t_refinement_restriction &&
                    temp(i, j, k) > thermal_Tcutoff * 0.9);
                tag = tag ||
                    (phi_refinement_on &&
                    gradphi.lpNorm<2>() * dr >= phi_refinement_criterion);
                tag = tag ||
                    (tempgrad.lpNorm<2>() * dr > t_refinement_criterion &&
                    eta(i, j, k) >= t_refinement_restriction);
                tag = tag ||
                    ((gradeta.lpNorm<2>() * dr * 2 > m_refinement_criterion ||
                    gradphi.lpNorm<2>() * dr >= thermal_phi_ref_initial) &&
                    time < thermal_end_initial_refine_t);

                if (tag)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    }
    else {
        for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<char> const& tags = a_tags.array(mfi);
            Set::Patch<const Set::Scalar> eta = eta_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> phi = phi_mf.Patch(lev,mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Vector gradeta = Numeric::Gradient(eta, i, j, k, 0, DX.data());
                Set::Vector gradphi = Numeric::Gradient(phi, i, j, k, 0, DX.data());

                bool tag = false;
                tag = tag ||
                    (gradeta.lpNorm<2>() * dr * 2 > m_refinement_criterion &&
                    eta(i, j, k) >= t_refinement_restriction);
                tag = tag ||
                    (phi_refinement_on &&
                    gradphi.lpNorm<2>() * dr >= phi_refinement_criterion);
                tag = tag ||
                    ((gradeta.lpNorm<2>() * dr * 2 > m_refinement_criterion ||
                    gradphi.lpNorm<2>() * dr >= thermal_phi_ref_initial) &&
                    time < thermal_end_initial_refine_t);

                if (tag)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    }
}

void Flame::Regrid(int lev, Set::Scalar time)
{
    BL_PROFILE("Integrator::Flame::Regrid");

    ic_phi->Initialize(lev, phi_mf, time);
    ic_eta->Initialize(lev, eta_0_mf, time);

    if (thermal.on) {
    /*
    This regrid function works by using the "has_exceeded_Tcutoff" field. If the temperature in a cell is greater than Tcutoff,
    eta will change and when regridding won't use the initial eta field. If T < T_cutoff, when regriding happens it applies the inital
    eta field condition. This gives at leat a 4x speed improvement in 2D when doing regression with voids. This is because orgionally
    there was a bug where when regridding, the orgional eta field wouldn't be applied, so there would be "squares" of voids instead of
    circles/spheres when using .xyzr files as the inital condition.
    */
    for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        Set::Patch<Set::Scalar> eta    = eta_mf.Patch(lev, mfi);
        Set::Patch<const Set::Scalar> temp = temp_mf.Patch(lev, mfi);
        Set::Patch<const Set::Scalar> eta_0 = eta_0_mf.Patch(lev, mfi);
        Set::Patch<Set::Scalar> exceeded_Tcutoff = thermal.has_exceeded_Tcutoff.Patch(lev, mfi);

        Set::Scalar Tcutoff = thermal.Tcutoff;

        amrex::BoxList boxes_to_update;
        if (lev == finest_level && prev_finest_level == finest_level)
            boxes_to_update = amrex::complementIn(bx, prev_finest_ba).boxList();
        else
            boxes_to_update.push_back(bx);

        for (const amrex::Box &box : boxes_to_update)
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {

                if (!exceeded_Tcutoff(i,j,k) && temp(i,j,k) < Tcutoff)
                {
                        eta(i, j, k) = eta_0(i, j, k);
                }
            });
    }

    if (lev == finest_level)
    {
        prev_finest_ba    = grids[finest_level];
        prev_finest_level = finest_level;
    }
    }
}

void Flame::Integrate(int amrlev, Set::Scalar time, int /*step*/,
    const amrex::MFIter& mfi, const amrex::Box& box)
{
    BL_PROFILE("Flame::Integrate");

    Base::Mechanics<model_type>::Integrate(amrlev,time,timestep,mfi,box);

    Set::Vector DX(geom[amrlev].CellSize());
    Set::Scalar dv = AMREX_D_TERM(DX(0), *DX(1), *DX(2));
    Set::Patch<const Set::Scalar> eta  = eta_mf.Patch(amrlev,mfi);
    Set::Patch<const Set::Scalar> mdot = mdot_mf.Patch(amrlev,mfi);

    // GPU: accumulating into chamber.* members directly inside a ParallelFor is a
    // parallel race on the device. Use an AMReX sum-reduction over the box (works
    // on CPU and GPU) and add the per-box partial sums to the members on the host.
    // The previous variable_pressure / else branches performed identical
    // accumulation, so they are collapsed here.
    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<Set::Scalar, Set::Scalar, Set::Scalar> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    reduce_op.eval(box, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Set::Scalar dvol  = (1.0 - eta(i, j, k, 0)) * dv;
            Set::Vector grad  = Numeric::Gradient(eta, i, j, k, 0, DX.data());
            Set::Scalar darea = grad.lpNorm<2>() * dv;
            Set::Scalar dmdot = mdot(i, j, k, 0) * dv;
            return {dvol, darea, dmdot};
        });

    ReduceTuple hv = reduce_data.value(reduce_op);
    chamber.volume += amrex::get<0>(hv);
    chamber.area   += amrex::get<1>(hv);
    chamber.mdot   += amrex::get<2>(hv);
    // time dependent pressure data from experimenta -> p = 0.0954521220950523 * exp(15.289993148880678 * t)
}
} // namespace Integrator

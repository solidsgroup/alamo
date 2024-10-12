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
#include "Base/Mechanics.H"

#include <cmath>

namespace Integrator
{

Flame::Flame() : Base::Mechanics<model_type>() {}

Flame::Flame(IO::ParmParse& pp) : Flame()
{
    pp_queryclass(*this);
}

// [parser]
void
Flame::Parse(Flame& value, IO::ParmParse& pp)
{
    BL_PROFILE("Integrator::Flame::Flame()");
    {
        //pp_query("fields_verbose", value.plot_field);
        pp_query("timestep", value.base_time);
        // These are the phase field method parameters
        // that you use to inform the phase field method.
        pp_query("pf.eps", value.pf.eps); // Burn width thickness
        pp_query("pf.kappa", value.pf.kappa); // Interface energy param
        pp_query("pf.gamma", value.pf.gamma); // Scaling factor for mobility
        pp_query("pf.lambda", value.pf.lambda); // Chemical potential multiplier
        pp_query("pf.w1", value.pf.w1); // Unburned rest energy
        pp_query("pf.w12", value.pf.w12);  // Barrier energy
        pp_query("pf.w0", value.pf.w0);    // Burned rest energy
        pp_query("amr.ghost_cells", value.ghost_count); // number of ghost cells in all fields
        pp_query("geometry.x_len", value.x_len); // Domain x length
        pp_query("geometry.y_len", value.y_len); // Domain y length

        value.bc_eta = new BC::Constant(1);
        pp_queryclass("pf.eta.bc", *static_cast<BC::Constant*>(value.bc_eta)); // See :ref:`BC::Constant`
        value.RegisterNewFab(value.eta_mf, value.bc_eta, 1, value.ghost_count, "eta", true);
        value.RegisterNewFab(value.eta_old_mf, value.bc_eta, 1, value.ghost_count, "eta_old", false);

        std::string eta_ic_type = "constant";
        pp_query("eta.ic.type", eta_ic_type); // Eta initial condition [constant, laminate, expression, bmp]
        if (eta_ic_type == "laminate") value.ic_eta = new IC::Laminate(value.geom, pp, "eta.ic.laminate");
        else if (eta_ic_type == "constant") value.ic_eta = new IC::Constant(value.geom, pp, "eta.ic.constant");
        else if (eta_ic_type == "expression") value.ic_eta = new IC::Expression(value.geom, pp, "eta.ic.expression");
        else if (eta_ic_type == "bmp") value.ic_eta = new IC::BMP(value.geom, pp, "eta.ic.bmp");
        else if (eta_ic_type == "png") value.ic_eta = new IC::PNG(value.geom, pp, "eta.ic.png");
        else Util::Abort(INFO, "Invalid eta IC type", eta_ic_type);
    }

    {
        //IO::ParmParse pp("thermal");
        pp_query("thermal.on", value.thermal.on); // Whether to use the Thermal Transport Model
        pp_query("elastic.on", value.elastic.on); // Whether to use Neo-hookean Elastic model
        pp_query("thermal.bound", value.thermal.bound); // System Initial Temperature
        pp_query("elastic.traction", value.elastic.traction); // Body force
        pp_query("elastic.phirefinement", value.elastic.phirefinement); // Phi refinement criteria 



        if (value.thermal.on) {
            pp_query("thermal.rho_ap", value.thermal.rho_ap); // AP Density
            pp_query("thermal.rho_htpb", value.thermal.rho_htpb); // HTPB Density
            pp_query("thermal.k_ap", value.thermal.k_ap); // AP Thermal Conductivity
            pp_query("thermal.k_htpb", value.thermal.k_htpb); // HTPB Thermal Conductivity
            pp_query("thermal.cp_ap", value.thermal.cp_ap); // AP Specific Heat
            pp_query("thermal.cp_htpb", value.thermal.cp_htpb); //HTPB Specific Heat

            pp_query("thermal.q0", value.thermal.q0); // Baseline heat flux

            pp_query("thermal.m_ap", value.thermal.m_ap); // AP Pre-exponential factor for Arrhenius Law
            pp_query("thermal.m_htpb", value.thermal.m_htpb); // HTPB Pre-exponential factor for Arrhenius Law
            pp_query("thermal.E_ap", value.thermal.E_ap); // AP Activation Energy for Arrhenius Law
            pp_query("thermal.E_htpb", value.thermal.E_htpb); // HTPB Activation Energy for Arrhenius Law

            pp_query("thermal.hc", value.thermal.hc); // Used to change heat flux units
            pp_query("thermal.massfraction", value.thermal.massfraction); // Systen AP mass fraction
            pp_query("thermal.mlocal_ap", value.thermal.mlocal_ap); // AP mass flux reference value 
            pp_query("thermal.mlocal_htpb", value.thermal.mlocal_htpb); // HTPB mass flux reference value 
            pp_query("thermal.mlocal_comb", value.thermal.mlocal_comb); // AP/HTPB mass flux reference value 

            pp_query("thermal.T_fluid", value.thermal.T_fluid); // Temperature of the Standin Fluid 

            pp_query("thermal.disperssion1", value.thermal.disperssion1); // K; dispersion variables are use to set the outter field properties for the void grain case.
            pp_query("thermal.disperssion1", value.thermal.disperssion2); // rho; dispersion variables are use to set the outter field properties for the void grain case.
            pp_query("thermal.disperssion1", value.thermal.disperssion3); // cp; dispersion variables are use to set the outter field properties for the void grain case.

            pp_query("thermal.modeling_ap", value.thermal.modeling_ap); // Scaling factor for AP thermal conductivity (default = 1.0)
            pp_query("thermal.modeling_htpb", value.thermal.modeling_htpb); // Scaling factor for HTPB thermal conductivity (default = 1.0)

            value.bc_temp = new BC::Constant(1);
            pp_queryclass("thermal.temp.bc", *static_cast<BC::Constant*>(value.bc_temp));
            value.RegisterNewFab(value.temp_mf, value.bc_temp, 1, value.ghost_count + 1, "temp", true);
            value.RegisterNewFab(value.temp_old_mf, value.bc_temp, 1, value.ghost_count + 1, "temp_old", false);
            value.RegisterNewFab(value.temps_mf, value.bc_temp, 1, value.ghost_count + 1, "temps", false);
            value.RegisterNewFab(value.temps_old_mf, value.bc_temp, 1, value.ghost_count + 1, "temps_old", false);

            value.RegisterNewFab(value.mdot_mf, value.bc_eta, 1, value.ghost_count + 1, "mdot", value.plot_field);
            value.RegisterNewFab(value.mob_mf, value.bc_eta, 1, value.ghost_count + 1, "mob", value.plot_field);
            value.RegisterNewFab(value.alpha_mf, value.bc_temp, 1, value.ghost_count + 1, "alpha", value.plot_field);
            value.RegisterNewFab(value.heatflux_mf, value.bc_temp, 1, value.ghost_count + 1, "heatflux", value.plot_field);
            value.RegisterNewFab(value.laser_mf, value.bc_temp, 1, value.ghost_count + 1, "laser", value.plot_field);

            value.RegisterIntegratedVariable(&value.volume, "total_area");
            value.RegisterIntegratedVariable(&value.area, "Interface_area");
            value.RegisterIntegratedVariable(&value.chamber_area, "chamber_area", false);
            value.RegisterIntegratedVariable(&value.massflux, "mass_flux");
            value.RegisterIntegratedVariable(&value.chamber_pressure, "Pressure", false);

            std::string laser_ic_type = "constant";
            pp_query("laser.ic.type", laser_ic_type); // heat laser initial condition type [constant, expression]
            if (laser_ic_type == "expression") value.ic_laser = new IC::Expression(value.geom, pp, "laser.ic.expression");
            else if (laser_ic_type == "constant") value.ic_laser = new IC::Constant(value.geom, pp, "laser.ic.constant");
            else Util::Abort(INFO, "Invalid eta IC type", laser_ic_type);

            std::string temp_ic_type;
            pp_query_validate("temp.ic.type", temp_ic_type,{"default","constant","expression","bmp","png"}); // Temperature initial condition
            if (temp_ic_type == "constant") value.thermal.ic_temp = new IC::Constant(value.geom, pp, "temp.ic.constant");
            else if (temp_ic_type == "expression") value.thermal.ic_temp = new IC::Expression(value.geom, pp, "temp.ic.expression");
            else if (temp_ic_type == "bmp") value.thermal.ic_temp = new IC::BMP(value.geom, pp, "temp.ic.bmp");
            else if (temp_ic_type == "png") value.thermal.ic_temp = new IC::PNG(value.geom, pp, "temp.ic.png");
            else if (temp_ic_type == "default") value.thermal.ic_temp = nullptr;
        }
    }

    {
        pp_query("pressure.P", value.pressure.P); // Constant pressure value
        if (value.thermal.on)
        {
            pp_query("pressure.a1", value.pressure.arrhenius.a1); // Surgate heat flux model paramater - AP
            pp_query("pressure.a2", value.pressure.arrhenius.a2); // Surgate heat flux model paramater - HTPB
            pp_query("pressure.a3", value.pressure.arrhenius.a3); // Surgate heat flux model paramater - Total
            pp_query("pressure.b1", value.pressure.arrhenius.b1); // Surgate heat flux model paramater - AP
            pp_query("pressure.b2", value.pressure.arrhenius.b2); // Surgate heat flux model paramater - HTPB
            pp_query("pressure.b3", value.pressure.arrhenius.b3); // Surgate heat flux model paramater - Total
            pp_query("pressure.c1", value.pressure.arrhenius.c1); // Surgate heat flux model paramater - Total
            pp_query("pressure.mob_ap", value.pressure.arrhenius.mob_ap); // Whether to include pressure to the arrhenius law
            pp_query("pressure.dependency", value.pressure.arrhenius.dependency); // Whether to use pressure to determined the reference Zeta 
            pp_query("pressure.h1", value.pressure.arrhenius.h1); // Surgate heat flux model paramater - Homogenized
            pp_query("pressure.h2", value.pressure.arrhenius.h2); // Surgate heat flux model paramater - Homogenized

        }
        else
        {
            pp_query("pressure.r_ap", value.pressure.power.r_ap);     // AP power pressure law parameter (r*P^n)
            pp_query("pressure.r_htpb", value.pressure.power.r_htpb); // HTPB power pressure law parameter (r*P^n)
            pp_query("pressure.r_comb", value.pressure.power.r_comb); // AP/HTPB power pressure law parameter (r*P^n)
            pp_query("pressure.n_ap", value.pressure.power.n_ap);     // AP power pressure law parameter (r*P^n)
            pp_query("pressure.n_htpb", value.pressure.power.n_htpb); // HTPB power pressure law parameter (r*P^n)
            pp_query("pressure.n_comb", value.pressure.power.n_comb); // AP/HTPB power pressure law parameter (r*P^n)
        }
        pp_query("variable_pressure", value.variable_pressure); // Whether to compute the pressure evolution
        pp_query("homogeneousSystem", value.homogeneousSystem); // Whether to initialize Phi with homogenized properties
    }

    pp_query("amr.refinement_criterion", value.m_refinement_criterion);// Refinement criterion for eta field   
    pp_query("amr.refinement_criterion_temp", value.t_refinement_criterion);// Refinement criterion for temperature field    
    pp_query("amr.refinament_restriction", value.t_refinement_restriction);// Eta value to restrict the refinament for the temperature field 
    pp_query("amr.phi_refinement_criterion", value.phi_refinement_criterion);// Refinement criterion for phi field [infinity]
    pp_query("small", value.small); // Lowest value of Eta.

    {
        // The material field is referred to as :math:`\phi(\mathbf{x})` and is 
        // specified using these parameters. 
        //IO::ParmParse pp("phi.ic");
        std::string phi_ic_type = "packedspheres";
        pp_query("phi.ic.type", phi_ic_type); // IC type (psread, laminate, constant)
        if (phi_ic_type == "psread") {
            value.ic_phi = new IC::PSRead(value.geom, pp, "phi.ic.psread");
            //value.ic_phicell = new IC::PSRead(value.geom, pp, "phi.ic.psread");
            pp_query("phi.ic.psread.eps", value.zeta); // AP/HTPB interface length
            pp_query("phi.zeta_0", value.zeta_0); // Reference interface length for heat integration
        }
        else if (phi_ic_type == "laminate") {
            value.ic_phi = new IC::Laminate(value.geom, pp, "phi.ic.laminate");
            //value.ic_phicell = new IC::Laminate(value.geom, pp, "phi.ic.laminate");
            pp_query("phi.ic.laminate.eps", value.zeta); // AP/HTPB interface length
            pp_query("phi.zeta_0", value.zeta_0); // Reference interface length for heat integration
        }
        else if (phi_ic_type == "expression") {
            value.ic_phi = new IC::Expression(value.geom, pp, "phi.ic.expression");
            //value.ic_phicell = new IC::Expression(value.geom, pp, "phi.ic.expression");
            pp_query("phi.zeta_0", value.zeta_0); // Reference interface length for heat integration
            pp_query("phi.zeta", value.zeta); // AP/HTPB interface length
        }
        else if (phi_ic_type == "constant") {
            value.ic_phi = new IC::Constant(value.geom, pp, "phi.ic.constant");
            //value.ic_phicell = new IC::Constant(value.geom, pp, "phi.ic.constant");
            pp_query("phi.zeta_0", value.zeta_0); // Reference interface length for heat integration
            pp_query("phi.zeta", value.zeta); // AP/HTPB interface length
        }
        else if (phi_ic_type == "bmp") {
            value.ic_phi = new IC::BMP(value.geom, pp, "phi.ic.bmp");
            //value.ic_phicell = new IC::BMP(value.geom, pp, "phi.ic.bmp");
            pp_query("phi.zeta_0", value.zeta_0); // Reference interface length for heat integration
            pp_query("phi.zeta", value.zeta); // AP/HTPB interface length
        }
        else if (phi_ic_type == "png") {
            value.ic_phi = new IC::PNG(value.geom, pp, "phi.ic.png");
            pp_query("phi.zeta_0", value.zeta_0); // Reference interface length for heat integration
            pp_query("phi.zeta", value.zeta); // AP/HTPB interface length
        }
        else Util::Abort(INFO, "Invalid IC type ", phi_ic_type);
        //value.RegisterNewFab(value.phi_mf, 1, "phi_cell", true);
        value.RegisterNodalFab(value.phi_mf, 1, value.ghost_count + 1, "phi", true);
        //value.RegisterNewFab(value.phicell_mf, value.bc_eta, 1, value.ghost_count + 1, "phi", true);
    }

    pp_queryclass("elastic",static_cast<Base::Mechanics<model_type>&>(value));

    if (value.m_type != Type::Disable)
    {
        value.elastic.Tref = value.thermal.bound;
        pp_query("Tref", value.elastic.Tref); // Initial temperature for thermal expansion computation
        pp_queryclass("model_ap", value.elastic.model_ap);
        pp_queryclass("model_htpb", value.elastic.model_htpb);

        value.bc_psi = new BC::Nothing();
        value.RegisterNewFab(value.psi_mf, value.bc_psi, 1, value.ghost_count, "psi", value.plot_psi);
        value.psi_on = true;
    }
}

void Flame::Initialize(int lev)
{
    BL_PROFILE("Integrator::Flame::Initialize");
    Util::Message(INFO, m_type);
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
            temp_mf[lev]->setVal(thermal.bound);
            temp_old_mf[lev]->setVal(thermal.bound);
            temps_mf[lev]->setVal(thermal.bound);
            temps_old_mf[lev]->setVal(thermal.bound);
        }
        alpha_mf[lev]->setVal(0.0);
        mob_mf[lev]->setVal(0.0);
        mdot_mf[lev]->setVal(0.0);
        heatflux_mf[lev]->setVal(0.0);
        // pressure_mf[lev]->setVal(1.0); // [error]
        thermal.w1 = 0.2 * pressure.P + 0.9;
        thermal.T_fluid = thermal.bound;
        ic_laser->Initialize(lev, laser_mf);
        thermal.mlocal_htpb = 685000.0 - 850e3 * thermal.massfraction;
    }
    if (variable_pressure) pressure.P = 1.0;
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
            amrex::Box bx = mfi.grownnodaltilebox() & domain;
            //amrex::Box bx = mfi.nodaltilebox();
            //bx.grow(1);
            amrex::Array4<model_type>        const& model = model_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& phi = phi_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& eta = eta_mf[lev]->array(mfi);
            amrex::Array4<Set::Vector> const& rhs = rhs_mf[lev]->array(mfi);
            // amrex::Array4<const Set::Scalar> const& Pressure = pressure_mf[lev]->array(mfi); // [error]

            if (elastic.on)
            {
                amrex::Array4<const Set::Scalar> const& temp = temp_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    //auto sten = Numeric::GetStencil(i, j, k, bx);
                    Set::Scalar phi_avg = phi(i, j, k, 0);
                    Set::Scalar temp_avg = Numeric::Interpolate::CellToNodeAverage(temp, i, j, k, 0);
                    Set::Vector grad_eta = Numeric::CellGradientOnNode(eta, i, j, k, 0, DX);
                    rhs(i, j, k) = elastic.traction * grad_eta;

                    model_type model_ap = elastic.model_ap;
                    model_ap.F0 -= Set::Matrix::Identity();
                    model_ap.F0 *= (temp_avg - elastic.Tref);
                    model_ap.F0 += Set::Matrix::Identity();
                    model_type model_htpb = elastic.model_htpb;
                    model_htpb.F0 -= Set::Matrix::Identity();
                    model_htpb.F0 *= (temp_avg - elastic.Tref);
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
        Set::Scalar domain_area = x_len * y_len;
        chamber_pressure = pressure.P;
        chamber_area = domain_area - volume;
        Util::Message(INFO, "Mass = ", massflux);
        Util::Message(INFO, "Pressure = ", pressure.P);
    }
}

void Flame::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrador::Flame::Advance");
    Base::Mechanics<model_type>::Advance(lev, time, dt);
    const Set::Scalar* DX = geom[lev].CellSize();

    if (true) //lev == finest_level) //(true)
    {
        if (thermal.on)
        {
            std::swap(eta_old_mf[lev], eta_mf[lev]);
            std::swap(temp_old_mf[lev], temp_mf[lev]);
            std::swap(temps_old_mf[lev], temps_mf[lev]);

            Numeric::Function::Polynomial<4> w(pf.w0,
                0.0,
                -5.0 * pf.w1 + 16.0 * pf.w12 - 11.0 * pf.w0,
                14.0 * pf.w1 - 32.0 * pf.w12 + 18.0 * pf.w0,
                -8.0 * pf.w1 + 16.0 * pf.w12 - 8.0 * pf.w0);
            Numeric::Function::Polynomial<3> dw = w.D();

            for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                // Phase fields
                amrex::Array4<Set::Scalar> const& etanew = (*eta_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const& eta = (*eta_old_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const& phi = (*phi_mf[lev]).array(mfi);
                //amrex::Array4<const Set::Scalar> const& phicell = (*phicell_mf[lev]).array(mfi);
                // Heat transfer fields
                amrex::Array4<Set::Scalar>       const& tempnew = (*temp_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const& temp = (*temp_old_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar>       const& alpha = (*alpha_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar>       const& tempsnew = (*temps_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar>       const& temps = (*temps_old_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar>       const& laser = (*laser_mf[lev]).array(mfi);
                // Diagnostic fields
                amrex::Array4<Set::Scalar> const& mob = (*mob_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& mdot = (*mdot_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& heatflux = (*heatflux_mf[lev]).array(mfi);

                if (variable_pressure) {
                    pressure.P = exp(0.00075 * massflux);
                    if (pressure.P > 10.0) {
                        pressure.P = 10.0;
                    }
                    else if (pressure.P <= 0.99) {
                        pressure.P = 0.99;
                    }
                    elastic.traction = pressure.P;

                }

                Set::Scalar zeta_2 = 0.000045 - pressure.P * 6.42e-6;
                Set::Scalar zeta_1;
                if (pressure.arrhenius.dependency == 1) zeta_1 = zeta_2;
                else zeta_1 = zeta_0;

                Set::Scalar k1 = pressure.arrhenius.a1 * pressure.P + pressure.arrhenius.b1 - zeta_1 / zeta;
                Set::Scalar k2 = pressure.arrhenius.a2 * pressure.P + pressure.arrhenius.b2 - zeta_1 / zeta;
                Set::Scalar k3 = 4.0 * log((pressure.arrhenius.c1 * pressure.P * pressure.P + pressure.arrhenius.a3 * pressure.P + pressure.arrhenius.b3) - k1 / 2.0 - k2 / 2.0);
                Set::Scalar k4 = pressure.arrhenius.h1 * pressure.P + pressure.arrhenius.h2;

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar phi_avg = Numeric::Interpolate::NodeToCellAverage(phi, i, j, k, 0);
                    Set::Scalar eta_lap = Numeric::Laplacian(eta, i, j, k, 0, DX);
                    Set::Scalar K; // Calculate effective thermal conductivity
                    Set::Scalar rho; // No special interface mixure rule is needed here.
                    Set::Scalar cp;
                    if (homogeneousSystem) {
                        K = (thermal.modeling_ap * thermal.k_ap * thermal.massfraction + thermal.modeling_htpb * thermal.k_htpb * (1.0 - thermal.massfraction)) * phi_avg + thermal.disperssion1 * (1. - phi_avg); // Calculate effective thermal conductivity
                        rho = (thermal.rho_ap * thermal.massfraction + thermal.rho_htpb * (1.0 - thermal.massfraction)) * phi_avg + thermal.disperssion2 * (1. - phi_avg); // No special interface mixure rule is needed here.
                        cp = (thermal.cp_ap * thermal.massfraction + thermal.cp_htpb * (1.0 - thermal.massfraction)) * phi_avg + thermal.disperssion3 * (1. - phi_avg);
                    }
                    else {
                        K = thermal.modeling_ap * thermal.k_ap * phi_avg + thermal.modeling_htpb * thermal.k_htpb * (1.0 - phi_avg); // Calculate effective thermal conductivity
                        rho = thermal.rho_ap * phi_avg + thermal.rho_htpb * (1.0 - phi_avg); // No special interface mixure rule is needed here.
                        cp = thermal.cp_ap * phi_avg + thermal.cp_htpb * (1.0 - phi_avg);
                    }
                    Set::Scalar df_deta = ((pf.lambda / pf.eps) * dw(eta(i, j, k)) - pf.eps * pf.kappa * eta_lap);
                    etanew(i, j, k) = eta(i, j, k) - mob(i, j, k) * dt * df_deta;
                    if (etanew(i, j, k) <= small) etanew(i, j, k) = small;

                    alpha(i, j, k) = K / rho / cp; // Calculate thermal diffusivity and store in fiel
                    mdot(i, j, k) = rho * fabs(eta(i, j, k) - etanew(i, j, k)) / dt; // deta/dt  

                    if (isnan(etanew(i, j, k)) || isnan(alpha(i, j, k)) || isnan(mdot(i, j, k))) {
                        Util::Message(INFO, etanew(i, j, k), "etanew contains nan (i=", i, " j= ", j, ")");
                        Util::Message(INFO, mdot(i, j, k), "mdot contains nan (i=", i, " j= ", j, ")");
                        Util::Message(INFO, alpha(i, j, k), "alpha contains nan (i=", i, " j= ", j, ")");
                        Util::Abort(INFO);
                    }
                });

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar phi_avg = Numeric::Interpolate::NodeToCellAverage(phi, i, j, k, 0);
                    //phi_avg = phicell(i,j,k);
                    Set::Scalar K;
                    Set::Scalar qflux;
                    Set::Scalar mlocal;
                    Set::Scalar mdota;
                    Set::Scalar mbase;

                    if (homogeneousSystem) {
                        qflux = k4 * phi_avg;
                        mlocal = (thermal.mlocal_ap) * thermal.massfraction + (thermal.mlocal_htpb) * (1.0 - thermal.massfraction);
                        mdota = fabs(mdot(i, j, k));
                        mbase = tanh(4.0 * mdota / (mlocal));
                        K = (thermal.modeling_ap * thermal.k_ap * thermal.massfraction + thermal.modeling_htpb * thermal.k_htpb * (1.0 - thermal.massfraction)) * phi_avg + thermal.disperssion1 * (1. - phi_avg);
                        heatflux(i, j, k) = (laser(i, j, k) * phi_avg + thermal.hc * mbase * qflux) / K;
                    }
                    else {
                        qflux = k1 * phi_avg +
                            k2 * (1.0 - phi_avg) +
                            (zeta_1 / zeta) * exp(k3 * phi_avg * (1.0 - phi_avg));
                        mlocal = (thermal.mlocal_ap) * phi_avg + (thermal.mlocal_htpb) * (1.0 - phi_avg) + thermal.mlocal_comb * phi_avg * (1.0 - phi_avg);
                        mdota = fabs(mdot(i, j, k));
                        mbase = tanh(4.0 * mdota / (mlocal));

                        K = thermal.modeling_ap * thermal.k_ap * phi_avg + thermal.modeling_htpb * thermal.k_htpb * (1.0 - phi_avg);
                        heatflux(i, j, k) = (thermal.hc * mbase * qflux + laser(i, j, k)) / K;
                    }

                    if (isnan(heatflux(i, j, k))) {
                        Util::Message(INFO, heatflux(i, j, k), "heat contains nan (i=", i, " j= ", j, ")");
                        Util::Abort(INFO);
                    }
                });

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
                    Set::Scalar Tsolid;
                    Tsolid = dTdt + temps(i, j, k) * (etanew(i, j, k) - eta(i, j, k)) / dt;
                    tempsnew(i, j, k) = temps(i, j, k) + dt * Tsolid;
                    tempnew(i, j, k) = etanew(i, j, k) * tempsnew(i, j, k) + (1.0 - etanew(i, j, k)) * thermal.T_fluid;
                    if (isnan(tempsnew(i, j, k)) || isnan(temps(i, j, k))) {
                        Util::Message(INFO, tempsnew(i, j, k), "tempsnew contains nan (i=", i, " j= ", j, ")");
                        Util::Message(INFO, temps(i, j, k), "temps contains nan (i=", i, " j= ", j, ")");
                        Util::Abort(INFO);
                    }

                });

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar phi_avg = Numeric::Interpolate::NodeToCellAverage(phi, i, j, k, 0);
                    //phi_avg = phicell(i,j,k);

                    Set::Scalar L;
                    if (pressure.arrhenius.mob_ap == 1) L = thermal.m_ap * pressure.P * exp(-thermal.E_ap / tempnew(i, j, k)) * phi_avg;
                    else L = thermal.m_ap * exp(-thermal.E_ap / tempnew(i, j, k)) * phi_avg;
                    L += thermal.m_htpb * exp(-thermal.E_htpb / tempnew(i, j, k)) * (1.0 - phi_avg);
                    //L += thermal.m_comb * (0.5 * tempnew(i, j, k) / thermal.bound) * phi(i, j, k) * (1.0 - phi(i, j, k));
                    if (tempnew(i, j, k) <= thermal.bound) mob(i, j, k) = 0;
                    else mob(i, j, k) = L;
                    //mob(i,j,k) = L;
                    if (isnan(mob(i, j, k))) {
                        Util::Message(INFO, mob(i, j, k), "mob contains nan (i=", i, " j= ", j, ")");
                        Util::Abort(INFO);
                    }
                });
            } // MFi For loop 

        } // thermal IF
        else
        {
            std::swap(eta_old_mf[lev], eta_mf[lev]);
            Numeric::Function::Polynomial<4> w(pf.w0,
                0.0,
                -5.0 * pf.w1 + 16.0 * pf.w12 - 11.0 * pf.w0,
                14.0 * pf.w1 - 32.0 * pf.w12 + 18.0 * pf.w0,
                -8.0 * pf.w1 + 16.0 * pf.w12 - 8.0 * pf.w0
            );
            Numeric::Function::Polynomial<3> dw = w.D();

            for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                //Phase Fields
                amrex::Array4<Set::Scalar> const& etanew = (*eta_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& eta = (*eta_old_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& phi = (*phi_mf[lev]).array(mfi);

                Set::Scalar fmod_ap = pressure.power.r_ap * pow(pressure.P, pressure.power.n_ap);
                Set::Scalar fmod_htpb = pressure.power.r_htpb * pow(pressure.P, pressure.power.n_htpb);
                Set::Scalar fmod_comb = pressure.power.r_comb * pow(pressure.P, pressure.power.n_comb);

                pressure.power.a_fit = -1.16582 * sin(pressure.P) - 0.681788 * cos(pressure.P) + 3.3563;
                pressure.power.b_fit = -0.708225 * sin(pressure.P) + 0.548067 * cos(pressure.P) + 1.55985;
                pressure.power.c_fit = -0.0130849 * sin(pressure.P) - 0.03597 * cos(pressure.P) + 0.00725694;


                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    //Set::Scalar phi_avg = Numeric::Interpolate::NodeToCellAverage(phi, i, j, k, 0);
                    Set::Scalar L;
                    if (homogeneousSystem == 1) {
                        Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
                        Set::Scalar angle = acos(grad_eta[0] / grad_eta.lpNorm<2>()) * 180 / 3.1415;
                        if (angle > 90) angle = angle - 90.0;
                        if (angle > 180) angle = angle - 180.0;
                        if (angle > 270) angle = angle - 270.0;
                        L = pressure.power.a_fit + pressure.power.b_fit * exp(-pressure.power.c_fit * angle);
                    }
                    else {
                        Set::Scalar fs_actual;
                        fs_actual = fmod_ap * phi(i, j, k)
                            + fmod_htpb * (1.0 - phi(i, j, k))
                            + 4.0 * fmod_comb * phi(i, j, k) * (1.0 - phi(i, j, k));
                        L = fs_actual / pf.gamma / (pf.w1 - pf.w0);
                    }
                    Set::Scalar eta_lap = Numeric::Laplacian(eta, i, j, k, 0, DX);
                    Set::Scalar df_deta = ((pf.lambda / pf.eps) * dw(eta(i, j, k)) - pf.eps * pf.kappa * eta_lap);
                    etanew(i, j, k) = eta(i, j, k) - L * dt * df_deta;

                    if (etanew(i, j, k) > eta(i, j, k)) etanew(i, j, k) = eta(i, j, k);
                });
            } // MFI For Loop 
        } // thermal ELSE      
    }// First IF
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
            volume += eta(i, j, k, 0) * dv;
            Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
            Set::Scalar normgrad = grad.lpNorm<2>();
            Set::Scalar da = normgrad * dv;
            area += da;

            Set::Vector mgrad = Numeric::Gradient(mdot, i, j, k, 0, DX);
            Set::Scalar mnormgrad = mgrad.lpNorm<2>();
            Set::Scalar dm = mnormgrad * dv;
            massflux += dm;

        });
    }
    else {
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            volume += eta(i, j, k, 0) * dv;
            Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
            Set::Scalar normgrad = grad.lpNorm<2>();
            Set::Scalar da = normgrad * dv;
            area += da;
        });
    }
    // time dependent pressure data from experimenta -> p = 0.0954521220950523 * exp(15.289993148880678 * t)
}
} // namespace Integrator



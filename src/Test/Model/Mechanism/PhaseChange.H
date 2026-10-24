#ifndef TEST_MODEL_MECHANISM_PHASECHANGE_H
#define TEST_MODEL_MECHANISM_PHASECHANGE_H

#include "AMReX_FArrayBox.H"
#include "Model/Mechanism/PhaseChange.H"

namespace Test::Model::Mechanism
{
inline int PhaseChange()
{
    using Phase = ::Model::Mechanism::PhaseChange;
    int failed = 0;
    auto check = [&](const char* label, bool pass)
    {
        failed += Util::Test::SubMessage(label, !pass);
    };
    auto near = [](double a, double b)
    {
        return std::abs(a - b) <= 2.0e-12 * std::max(1.0, std::abs(b));
    };
    const std::vector<std::string> names = {"binder_gas", "AP_gas", "binder", "AP_solid"};
    const std::vector<double> mw = {20.0, 40.0};
    std::vector<double> densities;
    std::array<double, 4> cp, kappa;
    int case_number = 0;
    auto parse = [&](double fraction, bool homogeneous, bool explicit_flag = true,
                     double ap_conductivity = 0.4186, double pure_heat = 0.0)
    {
        densities = {NAN, NAN, 920.0, 1950.0};
        cp = {NAN, NAN, 2418.29, 1297.90};
        kappa = {NAN, NAN, 0.13, ap_conductivity};
        IO::ParmParse pp("test_phase_" + std::to_string(case_number++));
        pp.add("in", std::string("binder"));
        pp.add("out", std::string("binder_gas"));
        pp.add("rate_multiplier", 1.85e5);
        pp.add("activation_temperature", 3000.0);
        pp.add("temperature_cutoff", 360.0);
        if (pure_heat != 0.0) pp.add("heat_release", pure_heat);
        if (explicit_flag) pp.add("homogeneous", homogeneous ? 1 : 0);
        if (homogeneous)
        {
            pp.add("homogeneous.ap_solid", std::string("AP_solid"));
            if (fraction < 0.0)
            {
                pp.add("homogeneous.total_mass_fraction", 0.85);
                pp.add("homogeneous.resolved_mass_fraction", 0.65);
            }
            else pp.add("homogeneous.mass_fraction", fraction);
            pp.add("homogeneous.ap_rate_multiplier", 2.75e6);
            pp.add("homogeneous.ap_activation_temperature", 6000.0);
            pp.add("homogeneous.binder_heat_release", -1.0e5);
            pp.add("homogeneous.ap_heat_release", 2.0e5);
            pp.add("homogeneous.ap_gas", std::string("AP_gas"));
        }
        pp.add("phase_field.type", std::string("allencahn"));
        pp.add("phase_field.allencahn.mobility", 0.01);
        pp.add("phase_field.allencahn.lambda", 8.0);
        pp.add("phase_field.allencahn.kappa", 1.0e-8);
        pp.add("phase_field.allencahn.w1", 1.0);
        Phase result;
        Phase::Parse(result, pp, names, 2, {2, 3}, densities, mw, 8314.462618,
                     cp.data(), kappa.data());
        return result;
    };
    bool roots = true, volumes = true, capacities = true;
    for (double x : {0.1, 0.999999999, 1.0, 1.000000001, 3.0, 10.0})
        for (double t : {0.0, 0.2, 0.6, 0.999, 1.0})
        {
            const double w = 1950.0*t / (1950.0*t + 920.0*(1.0-t));
            parse(w, true, true, 0.13*x);
            volumes &= near((1.0-w)*densities[2]/920.0+w*densities[2]/1950.0, 1.0);
            capacities &= near(densities[2]*cp[2],
                (1.0-t)*920.0*2418.29+t*1950.0*1297.90);
            const double f = kappa[2] / 0.13;
            roots &= f >= std::min(1.0, x) - 1.0e-14 && f <= std::max(1.0, x) + 1.0e-14;
            roots &= near(f-x, (1.0-t)*(1.0-x)*std::pow(f, 1.0/AMREX_SPACEDIM));
        }
    check("Physical conductivity root and pure/equal-conductivity limits", roots);
    check("Constituent volumes sum to the blended volume", volumes);
    check("Additive volumetric heat capacity", capacities);
    const auto base = parse(0.0, false, false);
    const auto off = parse(0.0, false);
    check("Default disabled", densities[2] == 920.0 && kappa[2] == 0.13 &&
        base.StabilityRate(1.0e-5, 900.0) == off.StabilityRate(1.0e-5, 900.0));
    const auto pure_binder = parse(0.0, true);
    check("Pure binder limit", densities[2] == 920.0 && kappa[2] == 0.13 &&
        pure_binder.StabilityRate(1.0e-5, 900.0) == base.StabilityRate(1.0e-5, 900.0));
    parse(1.0, true);
    check("Pure AP material limit", densities[2] == 1950.0 && kappa[2] == 0.4186);
    const auto blend = parse(-1.0, true);
    const double w = 4.0/7.0, t = 0.38614900314795375;
    const double rho_blend = 1.0 / ((1.0-w)/920.0+w/1950.0);
    check("85/65 whole-propellant mass accounting and blended thermal properties",
        near(densities[2], rho_blend) && near(cp[2], 1778.067142857143));
    const double ratio = std::pow(2.75e6 / 1.85e5, t) * std::exp(-3000.0*t/900.0);
    check("Volume-weighted Arrhenius kinetics", near(
        blend.StabilityRate(1.0e-5, 900.0) / base.StabilityRate(1.0e-5, 900.0), ratio));
    check("Temperature cutoff", blend.StabilityRate(1.0e-5, 360.0) == 0.0);
    const auto direct = parse(w, true);
    check("Equivalent direct blend fraction", near(
        direct.StabilityRate(1.0e-5, 900.0), blend.StabilityRate(1.0e-5, 900.0)));

    amrex::Box box(amrex::IntVect(AMREX_D_DECL(0,0,0)),
                   amrex::IntVect(AMREX_D_DECL(0,0,0)));
    amrex::FArrayBox density(box, 4), eta(box, 1), fractions(box, 2);
    density.setVal(0.0);
    eta.setVal(0.5);
    fractions.setVal(0.0);
    fractions.array()(0,0,0,0) = 0.5;
    const ::Model::Mechanism::State state = {
        density.const_array(), eta.const_array(), fractions.const_array(), 900.0, 3.0e6};
    const double delta_eta = -1.0e-3;
    double heat = 0.0;
    const double dilation = blend.Transfer(density.array(), state, delta_eta, 0,0,0, &heat);
    const auto rho = density.const_array();
    const double released_mass = -rho_blend * delta_eta;
    check("Mass conservation and subgrid AP product split",
        near(rho(0,0,0,0) + rho(0,0,0,1) + rho(0,0,0,2), 0.0) &&
        near(rho(0,0,0,0), released_mass*(1.0-w)) &&
        near(rho(0,0,0,1), released_mass*w) && rho(0,0,0,3) == 0.0);
    check("Product molecular weights enter projection dilation", near(dilation,
        delta_eta + released_mass*8314.462618*900.0/3.0e6*((1.0-w)/20.0+w/40.0)));
    check("Mass-weighted heat release and sign", near(heat,
        released_mass*((1.0-w)*-1.0e5+w*2.0e5)));
    check("Phase-change heat equals additive constituent energy per regressed volume",
        near(heat, -delta_eta*((1.0-t)*920.0*-1.0e5+t*1950.0*2.0e5)));
    heat = 0.0;
    base.Transfer(density.array(), state, delta_eta, 0,0,0, &heat);
    check("Disabled mode adds no heat", heat == 0.0);
    const auto endothermic = parse(0.0, false, true, 0.4186, -418400.0);
    density.setVal(0.0);
    heat = 0.0;
    endothermic.Transfer(density.array(), state, delta_eta, 0,0,0, &heat);
    check("Pure-material decomposition conserves mass and consumes specified heat",
        near(density.const_array()(0,0,0,0), -920.0*delta_eta) &&
        near(density.const_array()(0,0,0,2), 920.0*delta_eta) &&
        near(heat, -920.0*delta_eta*-418400.0));
    check("Local and gradient kinetics use the same blend factor",
        near(blend.Local(state,0,0,0) / base.Local(state,0,0,0), ratio) &&
        near(blend.GradientCoefficient(state,0,0,0) /
             base.GradientCoefficient(state,0,0,0), ratio));
    return failed;
}
}

#endif

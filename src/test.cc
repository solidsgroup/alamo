#include <algorithm>
#include <array>
#include <cmath>
#include <stdlib.h>

#include "Set/Matrix4.H"
#include "Util/Util.H"
#include "IO/FileNameParse.H"

#include "Test/Numeric/Stencil.H"
#include "Test/Set/Matrix4.H"

#include "Operator/Elastic.H"

#include "Numeric/Interpolator/Test.H"
#include "Numeric/Interpolator/Linear.H"

#include "Model/Solid/Linear/Isotropic.H"
#include "Model/Solid/Linear/Cubic.H"
#include "Model/Solid/Linear/Laplacian.H"
#include "Model/Solid/Linear/Transverse.H"
#include "Model/Solid/Affine/Isotropic.H"
#include "Model/Solid/Affine/Cubic.H"
#include "Model/Solid/Finite/NeoHookean.H"
#include "Model/Solid/Finite/NeoHookeanPredeformed.H"
#include "Model/Solid/Finite/PseudoLinear/Cubic.H"
#include "Model/Solid/Finite/PseudoAffine/Cubic.H"
#include "Model/Solid/Linear/Hexagonal.H"
#include "Model/Solid/Affine/Hexagonal.H"
#include "Model/Chemistry/Rocfire.H"
#include "Model/PhaseField/PhaseField.H"

#include "Solver/Local/Riemann/Roe.H"
#include "Solver/Local/ODE/BackwardEuler.H"
#include "Solver/Local/ODE/ForwardEuler.H"
#include "Solver/Local/ODE/ODE.H"

#include "Unit/Test.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc, argv);

    int failed = 0;

    Util::globalprefix = "  │  ";

    Util::Test::Message("Bounded constitutive phase interpolation");
    {
        bool bounded = true;
        constexpr Set::Scalar lo = 0.11920292202211755;
        constexpr Set::Scalar hi = 1.0 - lo;
        for (int i = 0; i <= 10000; ++i)
        {
            for (const Set::Scalar eta :
                 {-0.2 + 1.4 * i / 10000.0,
                  lo + (hi-lo) * (1.0 - i * 1.0e-10)})
            {
                const Set::Scalar h = Model::PhaseField::H(eta);
                bounded &= h >= 0.0 && h <= 1.0;
            }
        }
        failed += Util::Test::SubMessage("Phase weights remain in [0,1], including endpoint roundoff", !bounded);
        failed += Util::Test::SubMessage("Pure-phase limits and interface midpoint preserved",
            Model::PhaseField::H(-1.0) != 0.0 || Model::PhaseField::H(2.0) != 1.0 ||
            std::abs(Model::PhaseField::H(0.5) - 0.5) > 1.0e-14);
    }

    Util::Test::Message("Rocfire phase-change energy ownership test");
    {
        Model::Chemistry::Rocfire chemistry;
        const std::array<Set::Scalar, 4> progress = {2.0, 3.0, 0.0, 0.0};
        const auto legacy = chemistry.source_from_progress(progress);
        chemistry.condensed_heat_in_phase_change = true;
        const auto coupled = chemistry.source_from_progress(progress);
        failed += Util::Test::SubMessage("Gas reaction excludes condensed heat when coupled",
            std::abs(coupled.second - (2.0 * chemistry.qgas[0] + 3.0 * chemistry.qgas[1])) > 1.0e-8);
        failed += Util::Test::SubMessage("Legacy energy and reaction mass sources preserved",
            legacy.first != coupled.first ||
            std::abs(legacy.second - coupled.second -
                     (2.0 * chemistry.qsolid[0] + 3.0 * chemistry.qsolid[1])) > 1.0e-8);
    }

    Util::Test::Message("Rocfire fitted binder gas heat test");
    {
        Model::Chemistry::Rocfire chemistry;
        chemistry.binder_adiabatic_heat_release = true;
        chemistry.condensed_heat_in_phase_change = true;
        bool balance = true, pressure_dependence = true;
        for (const double w : {0.1, 0.52, 0.8})
        {
            chemistry.binder_equivalence_ratio =
                chemistry.equivalence_ratio_from_ap_mass_fraction(w);
            const auto low = chemistry.binder_adiabatic_flame_temperature(101325.0);
            const auto high = chemistry.binder_adiabatic_flame_temperature(100.0 * 101325.0);
            pressure_dependence &= high > low;
            pressure_dependence &= chemistry.binder_adiabatic_flame_temperature(200.0 * 101325.0) ==
                chemistry.binder_adiabatic_flame_temperature(150.0 * 101325.0);
            for (const double solid_heat : {-1.0e6, 0.0, 2.0e5})
            {
                chemistry.binder_solid_heat_release = solid_heat;
                const auto source = chemistry.source_from_progress({0.0, 1.0, 0.0, 0.0}, 100.0 * 101325.0);
                balance &= std::abs(source.second + solid_heat -
                    chemistry.binder_specific_heat * (high - chemistry.binder_initial_temperature)) < 1.0e-8;
            }
        }
        failed += Util::Test::SubMessage("Blended solid heat plus gas heat equals fitted total heat", !balance);
        failed += Util::Test::SubMessage("Pressure-dependent fit and pressure clamp", !pressure_dependence);
        chemistry.binder_adiabatic_heat_release = false;
        failed += Util::Test::SubMessage("Fixed-heat opt-out preserves qgas", chemistry.binder_heat_release(1.0e7) != chemistry.qgas[1]);
    }

    Util::Test::Message("Rocfire homogeneous stoichiometry test");
    {
        Model::Chemistry::Rocfire chemistry;
        bool parsed_composition = true;
        int composition_case = 0;
        for (const auto [w, expected] :
             {std::pair<double, double>{0.4, 4.3844},
              std::pair<double, double>{0.52, 3.30752},
              std::pair<double, double>{0.6, 2.5896}})
        {
            IO::ParmParse pp("test_rocfire_composition_" + std::to_string(composition_case++));
            pp.add("subgrid_ap_mass_fraction", w);
            pp.add("binder_adiabatic_heat_release", 1);
            Model::Chemistry::Rocfire parsed;
            parsed.nspecies = 6;
            Model::Chemistry::Rocfire::Parse(parsed, pp, 6);
            const auto source = parsed.source_from_progress({0.0, 0.0, 1.0, 0.0});
            parsed_composition &= std::abs(parsed.primary_stoichiometry() - expected) < 1.0e-12;
            parsed_composition &= std::abs(source.first[Model::Chemistry::Rocfire::AP] + expected) < 1.0e-12;
            parsed_composition &= std::abs(parsed.binder_equivalence_ratio -
                parsed.equivalence_ratio_from_ap_mass_fraction(w)) < 1.0e-12;
        }
        failed += Util::Test::SubMessage("Composition input determines primary demand and binder fit without beta input", !parsed_composition);
        IO::ParmParse reference_pp("test_rocfire_custom_reference");
        reference_pp.add("reference_beta", 5.0);
        reference_pp.add("stoichiometry_reference_ap_mass_fraction", 0.25);
        reference_pp.add("subgrid_ap_mass_fraction", 0.6);
        Model::Chemistry::Rocfire custom_reference;
        custom_reference.nspecies = 6;
        Model::Chemistry::Rocfire::Parse(custom_reference, reference_pp, 6);
        const auto custom_source = custom_reference.source_from_progress({0.0, 0.0, 1.0, 0.0});
        failed += Util::Test::SubMessage("Configurable reference beta and composition determine the kinetic coefficient",
            std::abs(custom_reference.primary_stoichiometry() - 2.2) > 1.0e-12 ||
            std::abs(custom_source.first[Model::Chemistry::Rocfire::AP] + 2.2) > 1.0e-12);
        IO::ParmParse stoichiometric_reference_pp("test_rocfire_zero_reference_demand");
        stoichiometric_reference_pp.add("reference_beta", 0.0);
        stoichiometric_reference_pp.add("stoichiometry_reference_ap_mass_fraction", 7.974 / 8.974);
        stoichiometric_reference_pp.add("subgrid_ap_mass_fraction", 7.974 / 8.974);
        Model::Chemistry::Rocfire zero_reference;
        Model::Chemistry::Rocfire::Parse(zero_reference, stoichiometric_reference_pp, 6);
        failed += Util::Test::SubMessage("A stoichiometric reference composition accepts zero external AP demand",
            zero_reference.primary_stoichiometry() != 0.0 ||
            zero_reference.source_from_progress({0.0, 0.0, 1.0, 1.0}).second != 0.0);
        IO::ParmParse pure_pp("test_rocfire_pure_htpb_reference");
        pure_pp.add("subgrid_ap_mass_fraction", 0.0);
        Model::Chemistry::Rocfire pure_htpb;
        Model::Chemistry::Rocfire::Parse(pure_htpb, pure_pp, 6);
        failed += Util::Test::SubMessage("Pure HTPB composition uses the derived pure oxidizer demand without enabling the fit",
            std::abs(pure_htpb.primary_stoichiometry() - 7.974) > 1.0e-12);
        bool balance = true, conservation = true;
        for (const double reference : {0.0, 0.52})
        {
            chemistry.stoichiometry_reference_ap_mass_fraction = reference;
            chemistry.homogenized_ap_mass_fraction = reference;
            balance &= std::abs(chemistry.primary_stoichiometry() - chemistry.reference_beta) < 1.0e-12;
            balance &= std::abs(chemistry.secondary_stoichiometry() - chemistry.gamma) < 1.0e-12;
            for (const double w : {0.0, 0.2, 0.52, 0.8})
            {
                chemistry.homogenized_ap_mass_fraction = w;
                balance &= std::abs(chemistry.primary_stoichiometry() - std::max(0.0,
                    (chemistry.reference_beta + reference) / (1.0 - reference) * (1.0 - w) - w)) < 1.0e-12;
                balance &= std::abs(chemistry.secondary_stoichiometry() - std::max(0.0,
                    (chemistry.gamma + reference) / (1.0 - reference) * (1.0 - w) - w)) < 1.0e-12;
                for (const std::array<Set::Scalar, 4> progress :
                     {std::array<Set::Scalar, 4>{0.0, 0.0, 1.0, 0.0},
                      std::array<Set::Scalar, 4>{0.0, 0.0, 0.0, 1.0}})
                {
                    const auto source = chemistry.source_from_progress(progress);
                    conservation &= std::abs(Model::Chemistry::Rocfire::sum_species(source.first)) < 1.0e-12;
                    const int oxidizer = progress[2] > 0.0 ? Model::Chemistry::Rocfire::AP : Model::Chemistry::Rocfire::Mono;
                    const auto coefficient = progress[2] > 0.0 ? chemistry.primary_stoichiometry() :
                        (chemistry.primary_stoichiometry() > 0.0 ? chemistry.secondary_stoichiometry() : 0.0);
                    conservation &= std::abs(source.first[oxidizer] +
                        coefficient) < 1.0e-12;
                }
            }
        }
        failed += Util::Test::SubMessage("Embedded AP plus external demand satisfies reference stoichiometry", !balance);
        failed += Util::Test::SubMessage("Both diffusion reactions use adjusted coefficients and conserve mass", !conservation);

        Model::Chemistry::Rocfire gated;
        gated.gas_constant = Set::Constant::Rg;
        // Positive amounts of every reactant ensure the composition gate,
        // rather than an absent species, is what suppresses the reactions.
        Model::Chemistry::SpeciesArray reactants{};
        for (int n = 0; n < 6; ++n) reactants[n] = 1.0 / 6.0;
        const double stoichiometric_w = 7.974 / 8.974;
        bool shutoff = true;
        gated.gamma = 20.0; // Still positive at primary stoichiometry.
        for (const double w : {stoichiometric_w, 0.95, 1.0})
        {
            gated.homogenized_ap_mass_fraction = w;
            const auto rates = gated.compute_progress_rates(1.0e7, 1500.0, reactants);
            const auto source = gated.source_from_progress({0.0, 0.0, 1.0, 1.0});
            shutoff &= gated.primary_stoichiometry() == 0.0 && rates[2] == 0.0 && rates[3] == 0.0;
            shutoff &= source.second == 0.0;
            for (const auto component : source.first) shutoff &= component == 0.0;
        }
        gated.homogenized_ap_mass_fraction = 0.4;
        const auto active = gated.compute_progress_rates(1.0e7, 1500.0, reactants);
        shutoff &= active[2] > 0.0 && active[3] > 0.0;
        failed += Util::Test::SubMessage("Both diffusion rates and their heat/mass sources vanish at and above stoichiometry", !shutoff);
        gated.gamma = 3.19;
        gated.homogenized_ap_mass_fraction = 0.886;
        const auto secondary_off = gated.compute_progress_rates(1.0e7, 1500.0, reactants);
        failed += Util::Test::SubMessage("Zero secondary demand suppresses conversion rather than creating Final without oxidizer",
            gated.secondary_stoichiometry() != 0.0 || secondary_off[2] <= 0.0 || secondary_off[3] != 0.0 ||
            gated.source_from_progress({0.0, 0.0, 0.0, 1.0}).second != 0.0);
    }

    Util::Test::Message("IO::FileNameParse test");
    {
        int subfailed = 0;
        const char *original = std::getenv("SLURM_FILENAME_PARSE_TEST");
        const bool restore_original = original != nullptr;
        const std::string original_value = restore_original ? original : "";

        setenv("SLURM_FILENAME_PARSE_TEST", "314159", 1);
        std::string available = "output_{SLURM_FILENAME_PARSE_TEST}";
        IO::FileNameParse(available);
        subfailed += Util::Test::SubMessage(
            "Slurm environment substitution", available != "output_314159");

        unsetenv("SLURM_FILENAME_PARSE_TEST");
        std::string unavailable = "output_{SLURM_FILENAME_PARSE_TEST}";
        IO::FileNameParse(unavailable);
        subfailed += Util::Test::SubMessage(
            "Unavailable Slurm variable fallback",
            unavailable != "output_SLURM_FILENAME_PARSE_TEST");

        if (restore_original)
            setenv("SLURM_FILENAME_PARSE_TEST", original_value.c_str(), 1);
        else
            unsetenv("SLURM_FILENAME_PARSE_TEST");

        failed += Util::Test::SubFinalMessage(subfailed);
    }

    #define MODELTEST(TYPE) \
        Util::Test::Message(#TYPE); \
        { \
            int subfailed = 0; \
            subfailed += Util::Test::SubMessage("PODTest",         TYPE::PODTest<TYPE>(true)); \
            subfailed += Util::Test::SubMessage("ArithmeticTest",  TYPE::ArithmeticTest<TYPE>(true)); \
            subfailed += Util::Test::SubMessage("DerivativeTest1", TYPE::DerivativeTest1<TYPE>(true)); \
            subfailed += Util::Test::SubMessage("DerivativeTest2", TYPE::DerivativeTest2<TYPE>(true)); \
            if (TYPE::kinvar == Model::Solid::KinematicVariable::F) \
            { \
                subfailed += Util::Test::SubMessage("MaterialFrameIndifference", TYPE::MaterialFrameIndifference<TYPE>(true)); \
            } \
            failed += Util::Test::SubFinalMessage(subfailed); \
        }
    MODELTEST(Model::Solid::Linear::Isotropic);
    MODELTEST(Model::Solid::Linear::Cubic);
    MODELTEST(Model::Solid::Linear::Laplacian);
    MODELTEST(Model::Solid::Linear::Transverse);
    MODELTEST(Model::Solid::Affine::Isotropic);
    MODELTEST(Model::Solid::Affine::Cubic);
    MODELTEST(Model::Solid::Linear::Hexagonal);
    MODELTEST(Model::Solid::Affine::Hexagonal);
    MODELTEST(Model::Solid::Finite::NeoHookean);
    MODELTEST(Model::Solid::Finite::PseudoLinear::Cubic);
    MODELTEST(Model::Solid::Finite::NeoHookeanPredeformed);
    MODELTEST(Model::Solid::Finite::PseudoAffine::Cubic);


    Test::Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Full>::Test();
    Test::Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic>::Test();
    Test::Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Diagonal>::Test();
    Test::Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Test();
    Test::Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major>::Test();

    Util::Test::Message("Numeric::Interpolator<Linear>");
    {
        int subfailed = 0;
        Numeric::Interpolator::Test<Numeric::Interpolator::Linear<Set::Scalar> > test;
        subfailed += Util::Test::SubMessage("Match",test.Match(0));
        failed += Util::Test::SubFinalMessage(subfailed);
    }

    Util::Test::Message("Numeric::Stencil test");
    {
        int subfailed = 0;
        Test::Numeric::Stencil test;
        test.Define(32);
        // first order
        subfailed += Util::Test::SubMessage("1-0-0",test.Derivative<1,0,0>(0));
        subfailed += Util::Test::SubMessage("0-1-0",test.Derivative<0,1,0>(0));
        // second order
        subfailed += Util::Test::SubMessage("2-0-0",test.Derivative<2,0,0>(0));
        subfailed += Util::Test::SubMessage("0-2-0",test.Derivative<0,2,0>(0));
        subfailed += Util::Test::SubMessage("0-0-1",test.Derivative<0,2,0>(0));
        subfailed += Util::Test::SubMessage("1-1-0",test.Derivative<1,1,0>(0));
        // fourth order
        subfailed += Util::Test::SubMessage("3-1-0",test.Derivative<3,1,0>(0));
        subfailed += Util::Test::SubMessage("1-3-0",test.Derivative<1,3,0>(0));
        subfailed += Util::Test::SubMessage("2-2-0",test.Derivative<2,2,0>(0));
        subfailed += Util::Test::SubMessage("4-0-0",test.Derivative<4,0,0>(0));
        subfailed += Util::Test::SubMessage("0-4-0",test.Derivative<0,4,0>(0));
#if AMREX_SPACEDIM>2
        // first order
        subfailed += Util::Test::SubMessage("0-0-1",test.Derivative<0,0,1>(0));
        // second order
        subfailed += Util::Test::SubMessage("0-0-2",test.Derivative<0,0,2>(0));
        subfailed += Util::Test::SubMessage("1-0-1",test.Derivative<1,0,1>(0));
        subfailed += Util::Test::SubMessage("0-1-1",test.Derivative<0,1,1>(0));
        // fourth order
        subfailed += Util::Test::SubMessage("0-0-4",test.Derivative<0,0,4>(0));
        subfailed += Util::Test::SubMessage("0-1-3",test.Derivative<0,1,3>(0));
        subfailed += Util::Test::SubMessage("0-3-1",test.Derivative<0,3,1>(0));
        subfailed += Util::Test::SubMessage("3-0-1",test.Derivative<3,0,1>(0));
        subfailed += Util::Test::SubMessage("1-0-3",test.Derivative<1,0,3>(0));
        subfailed += Util::Test::SubMessage("0-2-2",test.Derivative<0,2,2>(0));
        subfailed += Util::Test::SubMessage("2-0-2",test.Derivative<2,0,2>(0));
        subfailed += Util::Test::SubMessage("2-1-1",test.Derivative<2,1,1>(0));
        subfailed += Util::Test::SubMessage("1-2-1",test.Derivative<1,2,1>(0));
        subfailed += Util::Test::SubMessage("1-1-2",test.Derivative<1,1,2>(0));
#endif
        failed += Util::Test::SubFinalMessage(subfailed);
    }

    //Util::Test::Message("Solver::Nonlocal::Riemann::Roe test");
    //{
    //    int subfailed = 0;
    //    subfailed += Util::Test::SubMessage("Test",Solver::Local::Riemann::Roe::Test());
    //    failed += Util::Test::SubFinalMessage(subfailed);
    //}

    Util::Test::Message("Solver::Local::ODE test");
    {
        int subfailed = 0;
        using ODESolver = Solver::Local::ODE::ODE<
            Solver::Local::ODE::ForwardEuler,
            Solver::Local::ODE::BackwardEuler>;
        ODESolver solver;
        auto rhs = [] AMREX_GPU_HOST_DEVICE(
            const Set::Scalar* state, Set::Scalar* rate)
        {
            if (!(state[0] >= 0.0)) return false;
            rate[0] = -state[0] * state[0];
            return true;
        };
        auto jacobian = [] AMREX_GPU_HOST_DEVICE(
            const Set::Scalar* state, Set::Scalar* derivative)
        {
            derivative[0] = -2.0 * state[0];
            return true;
        };

        solver.Select<Solver::Local::ODE::ForwardEuler>();
        Set::Scalar explicit_state[1] = {1.0};
        const auto explicit_result = solver.Advance<1>(
            explicit_state, 1, 0.5, rhs, jacobian, false);
        const int explicit_failed = !explicit_result.converged ||
            std::abs(explicit_state[0] - 0.5) > 1.0e-14;
        subfailed += Util::Test::SubMessage("Forward Euler", explicit_failed);

        solver.Select<Solver::Local::ODE::BackwardEuler>();
        solver.Get<Solver::Local::ODE::BackwardEuler>().Configure(
            20, 1.0e-12, 1.0e-14);
        const Set::Scalar exact = std::sqrt(3.0) - 1.0;

        Set::Scalar finite_difference_state[1] = {1.0};
        const auto finite_difference_result = solver.Advance<1>(
            finite_difference_state, 1, 0.5, rhs, jacobian, false);
        const int finite_difference_failed = !finite_difference_result.converged ||
            std::abs(finite_difference_state[0] - exact) > 1.0e-11;
        subfailed += Util::Test::SubMessage(
            "Backward Euler finite-difference Jacobian", finite_difference_failed);

        Set::Scalar analytic_state[1] = {1.0};
        const auto analytic_result = solver.Advance<1>(
            analytic_state, 1, 0.5, rhs, jacobian, true);
        const int analytic_failed = !analytic_result.converged ||
            std::abs(analytic_state[0] - exact) > 1.0e-11;
        subfailed += Util::Test::SubMessage(
            "Backward Euler analytic Jacobian", analytic_failed);
        failed += Util::Test::SubFinalMessage(subfailed);
    }

    Util::Test::Message("Model::Chemistry::Rocfire Jacobian test");
    {
        int subfailed = 0;
        Model::Chemistry::Rocfire chemistry;
        chemistry.homogenized_ap_mass_fraction = 0.8;
        chemistry.condensed_heat_in_phase_change = true;
        chemistry.binder_adiabatic_heat_release = true;
        chemistry.binder_equivalence_ratio = chemistry.equivalence_ratio_from_ap_mass_fraction(0.8);
        chemistry.nspecies = 6;
        chemistry.gas_constant = Set::Constant::Rg;
        constexpr int size = 6;
        std::array<Set::Scalar, size> molecular_weight =
            {{26.0, 28.0, 24.0, 30.0, 22.0, 32.0}};
        std::array<Set::Scalar, size> cp_mass =
            {{800.0, 820.0, 840.0, 860.0, 880.0, 900.0}};
        Model::Gas::Gas::DeviceData gas = {
            size, Set::Constant::Rg, molecular_weight.data(), cp_mass.data(),
            nullptr, nullptr, 0, 0, nullptr, nullptr,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nullptr};
        const int dependent_species = Model::Chemistry::Rocfire::Primary;
        const Set::Scalar pressure = 2.0e6;
        const Set::Scalar mixture_density = 1200.0;
        const Set::Scalar gas_density = 0.6;
        const Set::Scalar temperature_scale = 1350.0;
        Set::Scalar state[size] = {0.18, 0.12, 0.15, 0.10, 0.08, 1.0};
        Set::Scalar analytic[size * size]{};

        auto rhs = [&](const Set::Scalar* z, Set::Scalar* rate)
        {
            Model::Chemistry::SpeciesArray Y{};
            Set::Scalar sum = 0.0;
            for (int column = 0; column < size - 1; ++column)
            {
                const int species = column < dependent_species ?
                    column : column + 1;
                Y[species] = z[column];
                sum += z[column];
            }
            Y[dependent_species] = 1.0 - sum;
            const Set::Scalar temperature = z[size - 1] * temperature_scale;

            Set::Scalar inverse_mw = 0.0;
            Set::Scalar cp = 0.0;
            for (int n = 0; n < size; ++n)
            {
                inverse_mw += Y[n] / molecular_weight[n];
                cp += Y[n] * cp_mass[n];
            }
            const Set::Scalar density = pressure /
                (Set::Constant::Rg * inverse_mw * temperature);
            Model::Chemistry::SpeciesArray rhoY{};
            for (int n = 0; n < size; ++n) rhoY[n] = density * Y[n];
            const auto source = chemistry.ComputeChemistrySources(
                pressure, temperature, rhoY, 0.0, nullptr);

            for (int row = 0; row < size - 1; ++row)
            {
                const int species = row < dependent_species ? row : row + 1;
                rate[row] = source.first[species] / density;
            }
            rate[size - 1] = gas_density * source.second /
                (density * mixture_density * cp * temperature_scale);
        };

        for (const double w : {0.4, 0.8, 0.886, 7.974 / 8.974, 0.95})
        {
            chemistry.homogenized_ap_mass_fraction = w;
            chemistry.binder_equivalence_ratio = chemistry.equivalence_ratio_from_ap_mass_fraction(w);
            const bool evaluated = chemistry.ComputeODEJacobian(
                state, analytic, pressure, mixture_density, gas_density, size,
                dependent_species, temperature_scale, &gas);
            Set::Scalar max_relative_error = 0.0;
            for (int column = 0; column < size && evaluated; ++column)
            {
                const Set::Scalar delta = 1.0e-6;
                Set::Scalar plus[size]{};
                Set::Scalar minus[size]{};
                Set::Scalar rate_plus[size]{};
                Set::Scalar rate_minus[size]{};
                for (int n = 0; n < size; ++n)
                    plus[n] = minus[n] = state[n];
                plus[column] += delta;
                minus[column] -= delta;
                rhs(plus, rate_plus);
                rhs(minus, rate_minus);

                for (int row = 0; row < size; ++row)
                {
                    const Set::Scalar numerical =
                        (rate_plus[row] - rate_minus[row]) / (2.0 * delta);
                    const Set::Scalar scale = std::max(
                        1.0, std::max(std::abs(numerical),
                                      std::abs(analytic[row * size + column])));
                    max_relative_error = std::max(max_relative_error,
                        std::abs(analytic[row * size + column] - numerical) / scale);
                }
            }
            const int jacobian_failed = !evaluated || max_relative_error > 1.0e-6;
            subfailed += Util::Test::SubMessage(
                "Analytic versus centered finite difference", jacobian_failed);
        }
        failed += Util::Test::SubFinalMessage(subfailed);
    }

    Util::Test::Message("Unit test");
    {
        int subfailed = 0;
        subfailed += Util::Test::SubMessage("Equivalence", UnitTest::Equivalence(1));
        failed += subfailed;
    }


    Util::globalprefix = "";
    Util::Message(INFO,failed," tests failed");

    Util::Finalize();
    return failed;
}

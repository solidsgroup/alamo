#include <algorithm>
#include <array>
#include <cmath>
#include <stdlib.h>

#include "AMReX_FArrayBox.H"
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
#include "Model/Chemistry/GrossModel.H"
#include "Model/Mechanism/PhaseChange.H"
#include "Model/Capillarity/MultiphaseFreeEnergy.H"
#include "Model/Capillarity/ConservativeAllenCahn.H"
#include "Model/Capillarity/SinglyDegenerateCahnHilliard.H"

#include "Solver/Local/Riemann/Roe.H"
#include "Solver/Local/ODE/BackwardEuler.H"
#include "Solver/Local/ODE/ForwardEuler.H"
#include "Solver/Local/ODE/ODE.H"

#include "Unit/Test.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc, argv);

    // The documentation builder invokes every executable in input-schema
    // traversal mode.  Unit tests are not input parsers and must not execute
    // against the synthetic traversal state.
    if (IO::ParmParse::InTraversalMode())
    {
        Util::Finalize();
        return 0;
    }

    int failed = 0;

    Util::globalprefix = "  │  ";

    Util::Test::Message("Arrhenius recession speed and implicit decomposition heat");
    {
        int subfailed = 0;
        for (int use_speed = 0; use_speed < 2; ++use_speed)
        {
            IO::ParmParse parameters("arrhenius_units_" + std::to_string(use_speed));
            parameters.add("phase0", std::string("solid"));
            parameters.add("phase1", std::string("binder_gas"));
            parameters.add("kinetics", std::string("arrhenius_surface_flux"));
            parameters.add("reference_pressure", 1.0e5);
            parameters.add("pressure_exponent", 0.0);
            parameters.add("reference_temperature", 1000.0);
            parameters.add("activation_temperature", 5000.0);
            parameters.add("coupled_enthalpy_change", 1.0e5);
            if (use_speed)
                parameters.add("pre_exponential_speed", 2.0);
            else
                parameters.add("reference_mass_flux", 920.0 * 2.0 * std::exp(-5.0));
            Model::Mechanism::PhaseChange mechanism;
            Model::Mechanism::PhaseChange::Parse(mechanism, parameters,
                {"AP_gas", "binder_gas", "solid"}, 2, {2}, {},
                {0.0, 0.0, 920.0}, {26.0, 26.0}, 8314.46261815324);

            const amrex::Box box(amrex::IntVect(AMREX_D_DECL(0,0,0)),
                                 amrex::IntVect(AMREX_D_DECL(2,0,0)));
            amrex::FArrayBox state_storage(box, 5, amrex::The_Managed_Arena());
            const auto state_array = state_storage.array();
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                state_array(i,j,k,0) = 0.0;
                state_array(i,j,k,1) = 1.0;
                state_array(i,j,k,2) = 460.0;
                Set::Scalar temperature = 800.0 + 200.0 * i;
                Set::Scalar heat = 0.0;
                Model::Mechanism::State state{state_array, {}, {},
                    temperature, 1.0e5, 1.0e-4};
                mechanism.ApplyKineticChange(state_array, state,
                    1.0e6, 3.0, temperature, heat, i, j, k);
                state_array(i,j,k,3) = temperature;
                state_array(i,j,k,4) = heat;
            });
            amrex::Gpu::streamSynchronize();
            for (int i = 0; i < 3; ++i)
            {
                const Set::Scalar initial_temperature = 800.0 + 200.0 * i;
                const Set::Scalar final_temperature = state_array(i,0,0,3);
                const Set::Scalar consumed = 460.0 - state_array(i,0,0,2);
                const Set::Scalar expected = 1.0e-4 * 3.0 * 920.0 * 2.0 *
                    std::exp(-5000.0 / final_temperature);
                subfailed += Util::Test::SubMessage(
                    "Dimensional surface flux at coupled temperature",
                    !(consumed > 0.0) || std::abs(consumed - expected) > 2.0e-11);
                subfailed += Util::Test::SubMessage(
                    "All condensed mass becomes binder gas",
                    std::abs(state_array(i,0,0,1) - 1.0 - consumed) > 1.0e-12 ||
                    state_array(i,0,0,0) != 0.0);
                subfailed += Util::Test::SubMessage(
                    "Endothermic heat is applied once",
                    !(final_temperature < initial_temperature) ||
                    std::abs(state_array(i,0,0,4) + 1.0e5 * consumed) > 1.0e-7 ||
                    std::abs(1.0e6 * (final_temperature - initial_temperature) -
                             state_array(i,0,0,4)) > 2.0e-7);
            }
        }
        failed += Util::Test::SubFinalMessage(subfailed);
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

    Util::Test::Message(
        "Model::Capillarity::MultiphaseFreeEnergy surface energy test");
    {
        int subfailed = 0;
        using PhaseModel = Model::Capillarity::MultiphaseFreeEnergy;
        constexpr int liquid = 0;
        constexpr int solid = 1;
        constexpr int gas = 2;
        constexpr int nphase = 3;
        const Set::Scalar sigma_lg = 0.8;
        const Set::Scalar solid_regularization = 0.3;
        const Set::Scalar surface_difference = -0.4;
        const Set::Scalar correction =
            surface_difference - solid_regularization;

        std::vector<Set::Scalar> surface_tension(
            nphase * nphase, 0.0);
        std::vector<Set::Scalar> regularization(
            nphase * nphase, 0.0);
        std::vector<Set::Scalar> difference(
            nphase * nphase, 0.0);
        surface_tension[liquid * nphase + gas] = sigma_lg;
        surface_tension[gas * nphase + liquid] = sigma_lg;
        regularization[liquid * nphase + solid] = solid_regularization;
        regularization[solid * nphase + liquid] = solid_regularization;
        difference[liquid * nphase + solid] = surface_difference;
        difference[solid * nphase + liquid] = surface_difference;

        PhaseModel model;
        const Set::Scalar ell = 0.7;
        const Set::Scalar surface_delta_regularization = 1.0e-12;
        const Set::Scalar regularization_gradient =
            surface_delta_regularization / ell;
        const amrex::GpuArray<int,AMREX_SPACEDIM> nonperiodic{};
        model.Define({"liquid", "solid", "gas"}, 1, true, ell,
                     surface_delta_regularization,
                     surface_tension, regularization, difference);
        subfailed += Util::Test::SubMessage(
            "Surface tension excludes liquid-solid",
            model.SurfaceTension(liquid,solid) != 0.0 ||
            std::abs(model.SurfaceTension(liquid,gas) - sigma_lg) > 1.0e-14);
        subfailed += Util::Test::SubMessage(
            "Signed solid surface correction",
            std::abs(model.SolidSurfaceCorrection(liquid,solid) -
                     correction) > 1.0e-14);
        subfailed += Util::Test::SubMessage(
            "Capillary stiffness includes correction",
            std::abs(model.MaximumCapillaryEnergy() - 1.7) > 1.0e-14);

        const Set::Scalar from_angle =
            PhaseModel::SurfaceEnergyDifferenceFromContactAngle(
                sigma_lg, Set::Constant::Pi / 3.0);
        subfailed += Util::Test::SubMessage(
            "Young contact-angle conversion",
            std::abs(from_angle - surface_difference) > 1.0e-12);

        constexpr int interpolation_points = 10000;
        Set::Scalar normalized_interpolation = 0.0;
        for (int i = 0; i < interpolation_points; ++i)
        {
            const Set::Scalar q =
                (static_cast<Set::Scalar>(i) + 0.5) /
                interpolation_points;
            normalized_interpolation += 2.0 * PhaseModel::Interpolation(q) /
                interpolation_points;
        }
        subfailed += Util::Test::SubMessage(
            "Solid surface delta normalization",
            std::abs(normalized_interpolation - 1.0) > 1.0e-13);
        subfailed += Util::Test::SubMessage(
            "Bulk surface delta is zero",
            PhaseModel::SurfaceDelta(
                0.0, regularization_gradient) != 0.0);
        subfailed += Util::Test::SubMessage(
            "Allen-Cahn pair mobility is interfacial",
            Model::Capillarity::ConservativeAllenCahn::
                    PairMobilityWeight(0.0, 1.0) != 0.0 ||
            Model::Capillarity::ConservativeAllenCahn::
                    PairMobilityWeight(1.0, 0.0) != 0.0 ||
            std::abs(Model::Capillarity::ConservativeAllenCahn::
                    PairMobilityWeight(0.5, 0.5) - 1.0) > 1.0e-14);
        subfailed += Util::Test::SubMessage(
            "Singly-degenerate pair mobility",
            Model::Capillarity::SinglyDegenerateCahnHilliard::
                    PairMobilityWeight(0.0, 1.0) != 0.0 ||
            Model::Capillarity::SinglyDegenerateCahnHilliard::
                    PairMobilityWeight(1.0, 0.0) != 0.0 ||
            std::abs(Model::Capillarity::SinglyDegenerateCahnHilliard::
                    PairMobilityWeight(0.5, 0.5) - 1.0) > 1.0e-14);

        // A diffuse tail can retain a nonzero gradient whose cube underflows
        // even though both phases are exactly absent.  The regularized
        // surface functional and all of its variations must remain finite in
        // that bulk state.
        const amrex::Box tail_domain(
            amrex::IntVect::TheZeroVector(),
            amrex::IntVect(AMREX_D_DECL(4, 4, 4)));
        amrex::FArrayBox tail_phase(tail_domain, nphase);
        tail_phase.setVal(0.0);
        auto tail = tail_phase.array();
        const amrex::Dim3 tail_lo = amrex::lbound(tail_domain);
        const amrex::Dim3 tail_hi = amrex::ubound(tail_domain);
        for (int k = tail_lo.z; k <= tail_hi.z; ++k)
            for (int j = tail_lo.y; j <= tail_hi.y; ++j)
                for (int i = tail_lo.x; i <= tail_hi.x; ++i)
                {
                    tail(i,j,k,solid) = 1.0e-110 * (i + 1.0);
                    tail(i,j,k,gas) = 1.0 - tail(i,j,k,solid);
                }
        const Set::Scalar tail_dx[AMREX_SPACEDIM] =
            {AMREX_D_DECL(1.0, 1.0, 1.0)};
        const auto tail_const = tail_phase.const_array();
        const Set::Scalar tail_mu =
            PhaseModel::SolidSurfaceCorrectionChemicalPotential(
                tail_const, liquid, solid, correction,
                regularization_gradient, 2, 2,
                AMREX_SPACEDIM > 2 ? 2 : 0, tail_dx,
                tail_domain, nonperiodic);
        const Set::Matrix tail_stress =
            PhaseModel::SolidSurfaceCorrectionCapillaryStress(
                tail_const, liquid, solid, correction,
                regularization_gradient, 2, 2,
                AMREX_SPACEDIM > 2 ? 2 : 0, tail_dx,
                Numeric::DefaultType());
        subfailed += Util::Test::SubMessage(
            "Finite liquid-solid bulk-tail variation",
            !std::isfinite(tail_mu) || !std::isfinite(tail_stress.norm()));

        // A disappearing solid can leave one isolated cell with a zero
        // centered gradient and a finite Laplacian.  The continuum-expanded
        // curvature formerly divided that Laplacian by the tiny surface-delta
        // regularization.  The discrete energy variation must instead remain
        // bounded independently of that regularization, agree with a finite
        // difference of the discrete energy, and exert a force that vanishes
        // continuously with the remnant amplitude.
        const amrex::Box remnant_domain(
            amrex::IntVect::TheZeroVector(),
            amrex::IntVect(AMREX_D_DECL(8, 8, 8)));
        amrex::FArrayBox remnant_phase(remnant_domain, nphase);
        remnant_phase.setVal(0.0);
        auto remnant = remnant_phase.array();
        const amrex::Dim3 remnant_lo = amrex::lbound(remnant_domain);
        const amrex::Dim3 remnant_hi = amrex::ubound(remnant_domain);
        for (int k = remnant_lo.z; k <= remnant_hi.z; ++k)
            for (int j = remnant_lo.y; j <= remnant_hi.y; ++j)
                for (int i = remnant_lo.x; i <= remnant_hi.x; ++i)
                {
                    remnant(i,j,k,liquid) = 0.9;
                    remnant(i,j,k,gas) = 0.1;
                }
        const int ri = 4;
        const int rj = AMREX_SPACEDIM > 1 ? 4 : 0;
        const int rk = AMREX_SPACEDIM > 2 ? 4 : 0;
        const Set::Scalar remnant_dx[AMREX_SPACEDIM] =
            {AMREX_D_DECL(1.0, 1.0, 1.0)};
        const Set::Scalar remnant_amplitude = 0.028917;
        remnant(ri,rj,rk,solid) = remnant_amplitude;
        const auto remnant_const = remnant_phase.const_array();
        const Set::Scalar remnant_mu =
            PhaseModel::SolidSurfaceCorrectionChemicalPotential(
                remnant_const, liquid, solid, correction,
                regularization_gradient, ri, rj, rk, remnant_dx,
                remnant_domain, nonperiodic);
        Set::Scalar inverse_spacing_sum = 0.0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            inverse_spacing_sum += 1.0 / remnant_dx[d];
        const Set::Scalar remnant_mu_bound =
            2.0 * std::abs(correction) * inverse_spacing_sum;
        subfailed += Util::Test::SubMessage(
            "Bounded isolated-remnant wetting variation",
            !std::isfinite(remnant_mu) ||
            std::abs(remnant_mu) >
                remnant_mu_bound * (1.0 + 1.0e-12));

        const auto remnant_energy = [&]()
        {
            Set::Scalar energy = 0.0;
            for (int k = remnant_lo.z; k <= remnant_hi.z; ++k)
                for (int j = remnant_lo.y; j <= remnant_hi.y; ++j)
                    for (int i = remnant_lo.x; i <= remnant_hi.x; ++i)
                        energy += PhaseModel::
                            SolidSurfaceCorrectionEnergyDensity(
                                remnant_const, liquid, solid, correction,
                                regularization_gradient, i, j, k,
                                remnant_dx, Numeric::GetStencil(
                                    i, j, k, remnant_domain, nonperiodic));
            return energy;
        };
        const Set::Scalar remnant_perturbation = 1.0e-7;
        remnant(ri,rj,rk,solid) =
            remnant_amplitude + remnant_perturbation;
        const Set::Scalar energy_plus = remnant_energy();
        remnant(ri,rj,rk,solid) =
            remnant_amplitude - remnant_perturbation;
        const Set::Scalar energy_minus = remnant_energy();
        remnant(ri,rj,rk,solid) = remnant_amplitude;
        const Set::Scalar remnant_energy_derivative =
            (energy_plus - energy_minus) /
            (2.0 * remnant_perturbation);
        subfailed += Util::Test::SubMessage(
            "Discrete wetting energy variation",
            std::abs(remnant_mu - remnant_energy_derivative) >
                1.0e-8 * (1.0 + std::abs(remnant_mu)));

        remnant(ri,rj,rk,solid) = 0.5 * remnant_amplitude;
        const Set::Scalar half_remnant_mu =
            PhaseModel::SolidSurfaceCorrectionChemicalPotential(
                remnant_const, liquid, solid, correction,
                regularization_gradient, ri, rj, rk, remnant_dx,
                remnant_domain, nonperiodic);
        subfailed += Util::Test::SubMessage(
            "Vanishing isolated-remnant wetting force",
            std::abs(half_remnant_mu * 0.5 * remnant_amplitude) >
                0.51 * std::abs(remnant_mu * remnant_amplitude));

        const amrex::Box flat_domain(
            amrex::IntVect::TheZeroVector(),
            amrex::IntVect(AMREX_D_DECL(8, 4, 4)));
        amrex::FArrayBox flat_phase(flat_domain, nphase);
        auto flat = flat_phase.array();
        const amrex::Dim3 flat_lo = amrex::lbound(flat_domain);
        const amrex::Dim3 flat_hi = amrex::ubound(flat_domain);
        const Set::Scalar flat_dx = ell / 4.0;
        for (int k = flat_lo.z; k <= flat_hi.z; ++k)
            for (int j = flat_lo.y; j <= flat_hi.y; ++j)
                for (int i = flat_lo.x; i <= flat_hi.x; ++i)
                {
                    const Set::Scalar x = (i - 4) * flat_dx;
                    const Set::Scalar q =
                        0.5 * (1.0 - std::tanh(2.0 * x / ell));
                    flat(i,j,k,liquid) = q;
                    flat(i,j,k,solid) = 1.0 - q;
                    flat(i,j,k,gas) = 0.0;
                }
        const auto flat_const = flat_phase.const_array();
        const int ci = 4;
        const int cj = 2;
        const int ck = AMREX_SPACEDIM > 2 ? 2 : 0;
        const Set::Scalar cell_size[AMREX_SPACEDIM] =
            {AMREX_D_DECL(flat_dx, flat_dx, flat_dx)};
        const auto central = Numeric::DefaultType();
        const Set::Scalar mu_liquid =
            PhaseModel::LiquidSurfaceCorrectionChemicalPotential(
                flat_const, liquid, solid, correction,
                regularization_gradient,
                ci, cj, ck, cell_size, central);
        const Set::Scalar mu_solid =
            PhaseModel::SolidSurfaceCorrectionChemicalPotential(
                flat_const, liquid, solid, correction,
                regularization_gradient,
                ci, cj, ck, cell_size, flat_domain, nonperiodic);
        const Set::Vector grad_liquid = Numeric::Gradient(
            flat_const, ci, cj, ck, liquid, cell_size, central);
        const Set::Vector grad_solid = Numeric::Gradient(
            flat_const, ci, cj, ck, solid, cell_size, central);
        const auto flat_surface_energy = [&]()
        {
            Set::Scalar energy = 0.0;
            for (int k = flat_lo.z; k <= flat_hi.z; ++k)
                for (int j = flat_lo.y; j <= flat_hi.y; ++j)
                    for (int i = flat_lo.x; i <= flat_hi.x; ++i)
                        energy += PhaseModel::
                            SolidSurfaceCorrectionEnergyDensity(
                                flat_const, liquid, solid, correction,
                                regularization_gradient, i, j, k,
                                cell_size, Numeric::GetStencil(
                                    i, j, k, flat_domain, nonperiodic));
            return energy;
        };
        const Set::Scalar flat_perturbation = 1.0e-7;
        const Set::Scalar flat_liquid = flat(ci,cj,ck,liquid);
        flat(ci,cj,ck,liquid) = flat_liquid + flat_perturbation;
        const Set::Scalar flat_liquid_energy_plus = flat_surface_energy();
        flat(ci,cj,ck,liquid) = flat_liquid - flat_perturbation;
        const Set::Scalar flat_liquid_energy_minus = flat_surface_energy();
        flat(ci,cj,ck,liquid) = flat_liquid;
        const Set::Scalar flat_solid = flat(ci,cj,ck,solid);
        flat(ci,cj,ck,solid) = flat_solid + flat_perturbation;
        const Set::Scalar flat_solid_energy_plus = flat_surface_energy();
        flat(ci,cj,ck,solid) = flat_solid - flat_perturbation;
        const Set::Scalar flat_solid_energy_minus = flat_surface_energy();
        flat(ci,cj,ck,solid) = flat_solid;
        const Set::Scalar flat_mu_liquid_fd =
            (flat_liquid_energy_plus - flat_liquid_energy_minus) /
            (2.0 * flat_perturbation);
        const Set::Scalar flat_mu_solid_fd =
            (flat_solid_energy_plus - flat_solid_energy_minus) /
            (2.0 * flat_perturbation);
        subfailed += Util::Test::SubMessage(
            "Flat-interface discrete wetting variation",
            std::abs(mu_liquid - flat_mu_liquid_fd) >
                1.0e-8 * (1.0 + std::abs(mu_liquid)) ||
            std::abs(mu_solid - flat_mu_solid_fd) >
                1.0e-8 * (1.0 + std::abs(mu_solid)));
        const Set::Matrix surface_stress =
            PhaseModel::SolidSurfaceCorrectionCapillaryStress(
                flat_const, liquid, solid, correction,
                regularization_gradient,
                ci, cj, ck, cell_size, central);
        Set::Scalar transverse_surface_stress = 0.0;
        for (int d = 1; d < AMREX_SPACEDIM; ++d)
            transverse_surface_stress = std::max(
                transverse_surface_stress,
                std::abs(surface_stress(d,d)));
        subfailed += Util::Test::SubMessage(
            "Liquid-solid stress uses pressure-reduced gauge",
            transverse_surface_stress > 1.0e-14);
        const Set::Scalar pair_mu_liquid =
            PhaseModel::PairChemicalPotential(
                flat_const, liquid, solid, solid_regularization, ell,
                ci, cj, ck, cell_size, central);
        const Set::Scalar pair_mu_solid =
            PhaseModel::PairChemicalPotential(
                flat_const, solid, liquid, solid_regularization, ell,
                ci, cj, ck, cell_size, central);
        subfailed += Util::Test::SubMessage(
            "Equilibrium binary constrained chemical potential",
            std::abs(pair_mu_liquid - pair_mu_solid) > 1.0e-12);

        constexpr int profile_points = 10000;
        const Set::Scalar dx = 16.0 * ell / profile_points;
        Set::Scalar pair_energy = 0.0;
        Set::Scalar correction_energy = 0.0;
        for (int i = 0; i < profile_points; ++i)
        {
            const Set::Scalar x = -8.0 * ell +
                (static_cast<Set::Scalar>(i) + 0.5) * dx;
            const Set::Scalar q =
                0.5 * (1.0 - std::tanh(2.0 * x / ell));
            const Set::Scalar s = 1.0 - q;
            const Set::Scalar grad_q = -4.0 * q * s / ell;
            const Set::Scalar grad_s = -grad_q;
            const Set::Scalar relative_gradient =
                q * grad_s - s * grad_q;
            pair_energy += solid_regularization *
                (0.75 * ell * relative_gradient * relative_gradient +
                 12.0 / ell * q * q * s * s) * dx;
            correction_energy += 2.0 * correction *
                PhaseModel::Interpolation(q) *
                PhaseModel::SurfaceDelta(
                    std::abs(grad_s), regularization_gradient) * dx;
        }
        subfailed += Util::Test::SubMessage(
            "Flat-interface regularization energy",
            std::abs(pair_energy - solid_regularization) > 2.0e-10);
        subfailed += Util::Test::SubMessage(
            "Flat-interface physical energy",
            std::abs(pair_energy + correction_energy -
                     surface_difference) > 2.0e-10);

        const Set::Matrix flat_pair_stress =
            PhaseModel::PairCapillaryStress(
                flat_const, liquid, solid, solid_regularization, ell,
                ci, cj, ck, cell_size, central);
        const Set::Vector relative_gradient =
            flat_const(ci,cj,ck,liquid) * grad_solid -
            flat_const(ci,cj,ck,solid) * grad_liquid;
        const Set::Matrix expected_pair_stress =
            -1.5 * ell * solid_regularization *
            relative_gradient * relative_gradient.transpose();
        subfailed += Util::Test::SubMessage(
            "Korteweg stress",
            (flat_pair_stress - expected_pair_stress).norm() > 1.0e-14);
        const Set::Matrix flat_stress = flat_pair_stress +
            PhaseModel::SolidSurfaceCorrectionCapillaryStress(
                flat_const, liquid, solid, correction,
                regularization_gradient,
                ci, cj, ck, cell_size, central);
        subfailed += Util::Test::SubMessage(
            "Korteweg stress symmetry",
            (flat_stress - flat_stress.transpose()).norm() > 1.0e-14);

#if AMREX_SPACEDIM == 2
        // Reproduce the finite-volume stress divergence used by LowMach on a
        // deliberately off-center, nonspherical diffuse interface.  Shared
        // face tractions must telescope to zero net force, and a symmetric
        // stress must also give zero first moment without relying on radial
        // symmetry of the phase field.
        const amrex::Box conservation_domain(
            amrex::IntVect(0,0), amrex::IntVect(32,28));
        const amrex::Box phase_box = amrex::grow(conservation_domain, 1);
        amrex::FArrayBox conservation_phase(phase_box, 2);
        auto conservation_eta = conservation_phase.array();
        const amrex::Dim3 phase_lo = amrex::lbound(phase_box);
        const amrex::Dim3 phase_hi = amrex::ubound(phase_box);
        for (int j = phase_lo.y; j <= phase_hi.y; ++j)
            for (int i = phase_lo.x; i <= phase_hi.x; ++i)
            {
                const Set::Scalar x = (i - 12.3) / 6.0;
                const Set::Scalar y = (j - 16.1) / 8.0;
                const Set::Scalar radius = std::sqrt(
                    x*x + y*y + 0.18*x*y);
                Set::Scalar q = 0.0;
                if (radius <= 0.65)
                    q = 1.0;
                else if (radius < 1.35)
                {
                    const Set::Scalar s =
                        (radius - 0.65) / (1.35 - 0.65);
                    q = 1.0 - s*s*(3.0 - 2.0*s);
                }
                conservation_eta(i,j,0,0) = q;
                conservation_eta(i,j,0,1) = 1.0 - q;
            }

        amrex::FArrayBox cell_stress(conservation_domain, 4);
        auto stress = cell_stress.array();
        const auto conservation_const = conservation_phase.const_array();
        const Set::Scalar conservation_dx[2] = {1.0, 1.0};
        const amrex::Dim3 conservation_lo =
            amrex::lbound(conservation_domain);
        const amrex::Dim3 conservation_hi =
            amrex::ubound(conservation_domain);
        for (int j = conservation_lo.y; j <= conservation_hi.y; ++j)
            for (int i = conservation_lo.x; i <= conservation_hi.x; ++i)
            {
                const Set::Matrix value = PhaseModel::PairCapillaryStress(
                    conservation_const, 0, 1, 0.8, 4.0,
                    i, j, 0, conservation_dx, central);
                for (int row = 0; row < 2; ++row)
                    for (int column = 0; column < 2; ++column)
                        stress(i,j,0,2*row+column) = value(row,column);
            }

        amrex::FArrayBox x_traction(
            amrex::surroundingNodes(conservation_domain, 0), 2);
        amrex::FArrayBox y_traction(
            amrex::surroundingNodes(conservation_domain, 1), 2);
        x_traction.setVal(0.0);
        y_traction.setVal(0.0);
        auto tx = x_traction.array();
        auto ty = y_traction.array();
        for (int j = conservation_lo.y; j <= conservation_hi.y; ++j)
            for (int i = conservation_lo.x + 1;
                 i <= conservation_hi.x; ++i)
                for (int component = 0; component < 2; ++component)
                    tx(i,j,0,component) = 0.5 *
                        (stress(i-1,j,0,2*component) +
                         stress(i,j,0,2*component));
        for (int j = conservation_lo.y + 1;
             j <= conservation_hi.y; ++j)
            for (int i = conservation_lo.x; i <= conservation_hi.x; ++i)
                for (int component = 0; component < 2; ++component)
                    ty(i,j,0,component) = 0.5 *
                        (stress(i,j-1,0,2*component+1) +
                         stress(i,j,0,2*component+1));

        Set::Vector total_force = Set::Vector::Zero();
        Set::Scalar total_torque = 0.0;
        Set::Scalar absolute_force = 0.0;
        for (int j = conservation_lo.y; j <= conservation_hi.y; ++j)
            for (int i = conservation_lo.x; i <= conservation_hi.x; ++i)
            {
                Set::Vector force;
                for (int component = 0; component < 2; ++component)
                {
                    force(component) =
                        tx(i+1,j,0,component) - tx(i,j,0,component) +
                        ty(i,j+1,0,component) - ty(i,j,0,component);
                    absolute_force += std::abs(force(component));
                }
                total_force += force;
                total_torque += (i + 0.5 - 12.3) * force(1) -
                                (j + 0.5 - 16.1) * force(0);
            }
        subfailed += Util::Test::SubMessage(
            "Conservative capillary force",
            total_force.norm() > 1.0e-13 * absolute_force);
        subfailed += Util::Test::SubMessage(
            "Conservative capillary torque",
            std::abs(total_torque) >
                1.0e-12 * absolute_force * 33.0);
#endif
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
        const amrex::Box domain(
            amrex::IntVect::TheZeroVector(),
            amrex::IntVect(AMREX_D_DECL(7, 7, 7)));
        const amrex::GpuArray<int,AMREX_SPACEDIM> periodic =
            {AMREX_D_DECL(1, 0, 0)};
        const auto physical_boundary =
            ::Numeric::GetStencil(0, 3, 0, domain);
        const auto periodic_boundary =
            ::Numeric::GetStencil(0, 3, 0, domain, periodic);
        subfailed += Util::Test::SubMessage(
            "physical boundary is one-sided",
            physical_boundary[0] != ::Numeric::StencilType::Hi);
        subfailed += Util::Test::SubMessage(
            "periodic boundary remains centered",
            periodic_boundary[0] != ::Numeric::StencilType::Central);
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

        IO::ParmParse ode_pp;
        ode_pp.add("ode_parse_test.type", "forward_euler");
        ode_pp.add("ode_parse_test.forward_euler.nsubsteps", 2);
        ODESolver parsed_solver;
        ode_pp.select<
            Solver::Local::ODE::ForwardEuler,
            Solver::Local::ODE::BackwardEuler>(
                "ode_parse_test", parsed_solver);
        Set::Scalar parsed_substep_state[1] = {1.0};
        const auto parsed_substep_result = parsed_solver.Advance<1>(
            parsed_substep_state, 1, 0.5, rhs, jacobian, false);
        const int parsed_substep_failed =
            !parsed_substep_result.converged ||
            std::abs(parsed_substep_state[0] - 0.609375) > 1.0e-14;
        subfailed += Util::Test::SubMessage(
            "Parsed ODE substeps", parsed_substep_failed);

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

        solver.Get<Solver::Local::ODE::BackwardEuler>().nsubsteps = 2;
        Set::Scalar implicit_substep_state[1] = {1.0};
        const auto implicit_substep_result = solver.Advance<1>(
            implicit_substep_state, 1, 0.5, rhs, jacobian, true);
        const Set::Scalar first_substep = 2.0 * (std::sqrt(2.0) - 1.0);
        const Set::Scalar implicit_substep_exact =
            2.0 * (std::sqrt(1.0 + first_substep) - 1.0);
        const int implicit_substep_failed =
            !implicit_substep_result.converged ||
            std::abs(implicit_substep_state[0] - implicit_substep_exact) >
                1.0e-11;
        subfailed += Util::Test::SubMessage(
            "Backward Euler substeps", implicit_substep_failed);
        failed += Util::Test::SubFinalMessage(subfailed);
    }

    Util::Test::Message("Model::Chemistry::GrossModel Jacobian test");
    {
        int subfailed = 0;
        Model::Chemistry::GrossModel chemistry;
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
        const int dependent_species = Model::Chemistry::GrossModel::Primary;
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

//
// :code:Model::Chemistry::GrossModel:code: implements a very specific finite rate
// kinetics model for AP/HTPB combustion using a 4 step, 6 species global reaction
// mechanism.  The mechanism utilizes a different form of equations for computing
// reaction rates than typical Cantera/Chemkin style kinetics, so this implementation
// is separate from :code:Model::Chemistry::FiniteRate:code:.
//
// The mechanism is available in public literature.  This implementation is
// named GrossModel after Matt Gross, the original author of the model.
// 
// Reference: Gross, M. L., Hedman, T. D., Son, S. F., Jackson, T. L., &
// Beckstead, M. W. (2013). Coupling micro and meso-scale combustion models of
// AP/HTPB propellants. Combustion and Flame, 160(5), 982–992.
// https://doi.org/10.1016/j.combustflame.2013.01.016
//

#ifndef MODEL_CHEMISTRY_GROSS_MODEL_H
#define MODEL_CHEMISTRY_GROSS_MODEL_H

#include "AMReX_Extension.H"
#include "Model/Chemistry/Chemistry.H"
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace Model
{
namespace Chemistry
{

struct GrossModel
{
    static constexpr const char* name = "gross_model";
    static constexpr bool reactive = true;
    static constexpr bool has_jacobian = true;
    static constexpr int nreactions = 4;
    static constexpr int AP = 0;
    static constexpr int Binder = 1;
    static constexpr int Mono = 2;
    static constexpr int Premixed = 3;
    static constexpr int Primary = 4;
    static constexpr int Final = 5;
    int nspecies = 1;

    enum class IntegrationType { EXPLICIT, SUBSTEP, IMPLICIT };
    IntegrationType integration = IntegrationType::EXPLICIT;

    std::array<Set::Scalar, nreactions> A = {{27000.0, 25300.0, 5000.0, 400.0}};
    std::array<Set::Scalar, nreactions> pressure_exponent = {{1.6, 1.7, 1.7, 1.7}};
    std::array<Set::Scalar, nreactions> activation_energy =
        {{15.9 * 4184.0, 14.5 * 4184.0, 11.3 * 4184.0, 14.5 * 4184.0}};
    // Homogeneous gas-reaction heat only. Condensed-phase energy belongs to
    // the separately configured phase-change mechanisms.
    std::array<Set::Scalar, nreactions> qgas =
        {{430.0 * 4184.0, 447.0 * 4184.0, 8365.0 * 4184.0, 2127.0 * 4184.0}};
    Set::Scalar pressure_reference = 100000.0;
    int nsubsteps = 1;
    int max_iter = 100;
    Set::Scalar reltol = 1.0e-8;
    Set::Scalar abstol = 1.0e-12;
    Set::Scalar beta = 7.974;
    Set::Scalar gamma = 3.19;
    Set::Scalar gas_constant = 0.0;

    GrossModel() = default;

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static Set::Scalar sum_species(const std::array<Set::Scalar, MAX_SPECIES>& values)
    {
        Set::Scalar sum = 0.0;
        for (int n = 0; n < MAX_SPECIES; ++n) sum += values[n];
        return sum;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    std::array<Set::Scalar, nreactions>
    compute_progress_rates(
        const Set::Scalar P,
        const Set::Scalar T,
        const std::array<Set::Scalar, MAX_SPECIES>& rhoY) const
    {
        std::array<Set::Scalar, nreactions> progress = {{0.0, 0.0, 0.0, 0.0}};

        Set::Scalar rho = 0.0;
        SpeciesArray Y{};
        for (int n = 0; n < MAX_SPECIES; ++n)
        {
            rho += rhoY[n];
        }
        if (!(rho > 0.0) || !(P > 0.0) || !(T > 0.0) || !(pressure_reference > 0.0))
        {
            return progress;
        }

        for (int n = 0; n < MAX_SPECIES; ++n)
        {
            Y[n] = std::max(rhoY[n], 0.0) / rho;
        }

        std::array<Set::Scalar, nreactions> k;
        for (int r = 0; r < nreactions; ++r)
        {
            k[r] = A[r] * std::pow(P / pressure_reference, pressure_exponent[r])
                * std::exp(-activation_energy[r] / (gas_constant * T));
        }

        progress[0] = k[0] * Y[AP];
        progress[1] = k[1] * Y[Binder];
        progress[2] = k[2] * Y[AP] * Y[Binder];
        progress[3] = k[3] * Y[Mono] * Y[Premixed];
        return progress;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    std::pair<std::array<Set::Scalar, MAX_SPECIES>, Set::Scalar>
    source_from_progress(const std::array<Set::Scalar, nreactions>& progress) const
    {
        SpeciesArray wdot{};
        for (int n = 0; n < MAX_SPECIES; ++n) wdot[n] = 0.0;

        wdot[AP]       -= progress[0];
        wdot[Mono]     += progress[0];

        wdot[Binder]   -= progress[1];
        wdot[Premixed] += progress[1];

        wdot[AP]       -= beta * progress[2];
        wdot[Binder]   -= progress[2];
        wdot[Primary]  += (beta + 1.0) * progress[2];

        wdot[Mono]     -= gamma * progress[3];
        wdot[Premixed] -= progress[3];
        wdot[Final]    += (gamma + 1.0) * progress[3];

        const Set::Scalar qdot =
            qgas[0] * progress[0]
            + qgas[1] * progress[1]
            + /*(beta + 1.0) * */qgas[2] * progress[2]
            + /*(gamma + 1.0) **/ qgas[3] * progress[3];

        return {wdot, qdot};
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    bool ComputeODEJacobian(
        const Set::Scalar* state,
        Set::Scalar* derivative,
        Set::Scalar pressure,
        Set::Scalar mixture_density,
        Set::Scalar gas_density,
        int ode_size,
        int dependent_species,
        Set::Scalar temperature_scale,
        const Model::Gas::Gas::DeviceData* gas) const
    {
        if (gas == nullptr || ode_size != 6 || nspecies != 6 ||
            dependent_species < 0 || dependent_species >= nspecies ||
            !(pressure > 0.0) || !(mixture_density > 0.0) ||
            !(gas_density > 0.0) || !(temperature_scale > 0.0) ||
            !(pressure_reference > 0.0) ||
            !(Model::Gas::Gas::UniversalGasConstant(*gas) > 0.0))
            return false;

        SpeciesArray Y{};
        int independent_species[5]{};
        Set::Scalar sum = 0.0;
        for (int column = 0; column < nspecies - 1; ++column)
        {
            const int species = column < dependent_species ? column : column + 1;
            independent_species[column] = species;
            Y[species] = state[column];
            if (!(Y[species] >= 0.0) || !std::isfinite(Y[species])) return false;
            sum += Y[species];
        }
        Y[dependent_species] = 1.0 - sum;
        if (!(Y[dependent_species] >= 0.0) ||
            !std::isfinite(Y[dependent_species])) return false;

        const Set::Scalar temperature = state[nspecies - 1] * temperature_scale;
        if (!(temperature > 0.0) || !std::isfinite(temperature)) return false;

        Set::Scalar inverse_mw = 0.0;
        Set::Scalar cp = 0.0;
        Set::Scalar cp_species[6]{};
        Set::Scalar dcp_species_dT[6]{};
        for (int n = 0; n < nspecies; ++n)
        {
            const Set::Scalar molecular_weight =
                Model::Gas::Gas::MolecularWeight(*gas, n);
            if (!(molecular_weight > 0.0)) return false;
            cp_species[n] =
                Model::Gas::Gas::CpMassSpecies(*gas, temperature, n);
            dcp_species_dT[n] =
                Model::Gas::Gas::DCpMassSpeciesDT(*gas, temperature, n);
            inverse_mw += Y[n] / molecular_weight;
            cp += Y[n] * cp_species[n];
        }
        if (!(inverse_mw > 0.0) || !(cp > 0.0)) return false;

        std::array<Set::Scalar, nreactions> rate_constant{};
        std::array<Set::Scalar, nreactions> progress{};
        for (int r = 0; r < nreactions; ++r)
            rate_constant[r] = A[r] *
                std::pow(pressure / pressure_reference, pressure_exponent[r]) *
                std::exp(-activation_energy[r] /
                    (gas_constant * temperature));

        progress[0] = rate_constant[0] * Y[AP];
        progress[1] = rate_constant[1] * Y[Binder];
        progress[2] = rate_constant[2] * Y[AP] * Y[Binder];
        progress[3] = rate_constant[3] * Y[Mono] * Y[Premixed];
        const auto source = source_from_progress(progress);

        const Set::Scalar inverse_density =
            Model::Gas::Gas::UniversalGasConstant(*gas) *
            inverse_mw * temperature / pressure;
        const Set::Scalar temperature_prefactor =
            gas_density / (mixture_density * temperature_scale);

        for (int column = 0; column < nspecies; ++column)
        {
            SpeciesArray dY{};
            Set::Scalar dtemperature = 0.0;
            if (column < nspecies - 1)
            {
                dY[independent_species[column]] = 1.0;
                dY[dependent_species] = -1.0;
            }
            else
                dtemperature = temperature_scale;

            std::array<Set::Scalar, nreactions> dprogress{};
            const Set::Scalar temperature_factor = dtemperature /
                (gas_constant * temperature * temperature);
            dprogress[0] = rate_constant[0] * dY[AP] +
                progress[0] * activation_energy[0] * temperature_factor;
            dprogress[1] = rate_constant[1] * dY[Binder] +
                progress[1] * activation_energy[1] * temperature_factor;
            dprogress[2] = rate_constant[2] *
                (Y[Binder] * dY[AP] + Y[AP] * dY[Binder]) +
                progress[2] * activation_energy[2] * temperature_factor;
            dprogress[3] = rate_constant[3] *
                (Y[Premixed] * dY[Mono] + Y[Mono] * dY[Premixed]) +
                progress[3] * activation_energy[3] * temperature_factor;
            const auto dsource = source_from_progress(dprogress);

            Set::Scalar dinverse_mw = 0.0;
            Set::Scalar dcp = 0.0;
            for (int n = 0; n < nspecies; ++n)
            {
                dinverse_mw += dY[n] /
                    Model::Gas::Gas::MolecularWeight(*gas, n);
                dcp += dY[n] * cp_species[n] +
                    Y[n] * dcp_species_dT[n] * dtemperature;
            }
            const Set::Scalar dinverse_density =
                Model::Gas::Gas::UniversalGasConstant(*gas) / pressure *
                (temperature * dinverse_mw + inverse_mw * dtemperature);

            for (int row = 0; row < nspecies - 1; ++row)
            {
                const int species = independent_species[row];
                derivative[row * nspecies + column] =
                    dinverse_density * source.first[species] +
                    inverse_density * dsource.first[species];
            }
            derivative[(nspecies - 1) * nspecies + column] =
                temperature_prefactor *
                ((dinverse_density * source.second +
                    inverse_density * dsource.second) / cp -
                 inverse_density * source.second * dcp / (cp * cp));
        }

        for (int n = 0; n < nspecies * nspecies; ++n)
            if (!std::isfinite(derivative[n])) return false;
        return true;
    }

    template<typename ODESolver>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Model::Chemistry::Advance::Result
    Advance(
        const Set::Scalar dt,
        const Set::Scalar pressure,
        const Set::Scalar mixture_density,
        SpeciesArray& rhoY,
        Set::Scalar& temperature,
        const Model::Gas::Gas::DeviceData& gas,
        const ODESolver& solver) const
    {
        Model::Chemistry::Advance::Result advance;
        if (!(dt > 0.0)) return advance;
        if (nspecies != 6 || !(pressure > 0.0) ||
            !(mixture_density > 0.0) || !(temperature > 0.0))
        {
            advance.converged = false;
            advance.residual_norm =
                std::numeric_limits<Set::Scalar>::infinity();
            return advance;
        }

        SpeciesArray mass_fraction{};
        Set::Scalar gas_density = 0.0;
        int dependent_species = 0;
        for (int n = 0; n < nspecies; ++n)
        {
            rhoY[n] = rhoY[n] > 0.0 ? rhoY[n] : 0.0;
            gas_density += rhoY[n];
            if (rhoY[n] > rhoY[dependent_species]) dependent_species = n;
        }
        if (!(gas_density > 0.0))
        {
            advance.converged = false;
            advance.residual_norm =
                std::numeric_limits<Set::Scalar>::infinity();
            return advance;
        }

        Set::Scalar independent_mass_fraction = 0.0;
        for (int n = 0; n < nspecies; ++n)
        {
            mass_fraction[n] = rhoY[n] / gas_density;
            if (n != dependent_species)
                independent_mass_fraction += mass_fraction[n];
        }
        mass_fraction[dependent_species] = 1.0 - independent_mass_fraction;

        Set::Scalar inverse_mw = 0.0;
        for (int n = 0; n < nspecies; ++n)
            inverse_mw += mass_fraction[n] /
                Model::Gas::Gas::MolecularWeight(gas, n);
        const Set::Scalar gas_constant =
            Model::Gas::Gas::UniversalGasConstant(gas);
        Set::Scalar rho_eos =
            pressure / (gas_constant * inverse_mw * temperature);
        Set::Scalar remaining = dt;
        Set::Scalar chemistry_dt = dt;
        int rejected_steps = 0;
        const GrossModel model = *this;
        const auto* gas_view = &gas;

        while (remaining > 0.0)
        {
            const SpeciesArray mass_fraction_old = mass_fraction;
            const Set::Scalar temperature_old = temperature;
            const Set::Scalar rho_eos_old = rho_eos;
            const Set::Scalar temperature_scale =
                temperature_old > 1.0 ? temperature_old : 1.0;
            Set::Scalar state[6]{};
            for (int n = 0; n < nspecies - 1; ++n)
            {
                const int species = n < dependent_species ? n : n + 1;
                state[n] = mass_fraction_old[species];
            }
            state[nspecies - 1] = temperature_old / temperature_scale;

            auto rhs = [=] AMREX_GPU_HOST_DEVICE(
                const Set::Scalar* z, Set::Scalar* rate)
            {
                SpeciesArray Y{};
                Set::Scalar sum = 0.0;
                for (int n = 0; n < 5; ++n)
                {
                    const int species =
                        n < dependent_species ? n : n + 1;
                    Y[species] = z[n];
                    if (!(Y[species] >= 0.0) ||
                        !std::isfinite(Y[species])) return false;
                    sum += Y[species];
                }
                Y[dependent_species] = 1.0 - sum;
                if (!(Y[dependent_species] >= 0.0) ||
                    !std::isfinite(Y[dependent_species])) return false;

                const Set::Scalar trial_temperature =
                    z[5] * temperature_scale;
                if (!(trial_temperature > 0.0) ||
                    !std::isfinite(trial_temperature)) return false;

                Set::Scalar trial_inverse_mw = 0.0;
                Set::Scalar cp = 0.0;
                for (int n = 0; n < 6; ++n)
                {
                    trial_inverse_mw += Y[n] /
                        Model::Gas::Gas::MolecularWeight(*gas_view, n);
                    cp += Y[n] * Model::Gas::Gas::CpMassSpecies(
                        *gas_view, trial_temperature, n);
                }
                if (!(trial_inverse_mw > 0.0) || !(cp > 0.0))
                    return false;

                const Set::Scalar rho_state = pressure /
                    (gas_constant * trial_inverse_mw * trial_temperature);
                SpeciesArray intrinsic_rhoY{};
                for (int n = 0; n < 6; ++n)
                    intrinsic_rhoY[n] = rho_state * Y[n];
                const Source source = model.ComputeChemistrySources(
                    pressure, trial_temperature, intrinsic_rhoY, 0.0,
                    nullptr);

                for (int n = 0; n < 5; ++n)
                {
                    const int species =
                        n < dependent_species ? n : n + 1;
                    rate[n] = source.first[species] / rho_state;
                }
                rate[5] = gas_density * source.second /
                    (rho_state * mixture_density * cp * temperature_scale);

                for (int n = 0; n < 6; ++n)
                    if (!std::isfinite(rate[n])) return false;
                return true;
            };

            auto jacobian = [=] AMREX_GPU_HOST_DEVICE(
                const Set::Scalar* z, Set::Scalar* derivative)
            {
                return model.ComputeODEJacobian(
                    z, derivative, pressure, mixture_density, gas_density,
                    6, dependent_species, temperature_scale, gas_view);
            };

            const Solver::Local::ODE::Advance::Result result =
                solver.template Advance<6>(
                    state, 6, chemistry_dt, rhs, jacobian,
                    solver.SupportsJacobian());
            if (!result.converged)
            {
                chemistry_dt *= 0.5;
                ++rejected_steps;
                if (rejected_steps > 40 ||
                    !(chemistry_dt > dt *
                        std::numeric_limits<Set::Scalar>::epsilon()))
                {
                    advance.converged = false;
                    advance.failed_dt = chemistry_dt;
                    advance.residual_norm = result.residual_norm;
                    return advance;
                }
                continue;
            }

            Set::Scalar sum = 0.0;
            for (int n = 0; n < 5; ++n)
            {
                const int species =
                    n < dependent_species ? n : n + 1;
                mass_fraction[species] = state[n];
                sum += state[n];
            }
            mass_fraction[dependent_species] = 1.0 - sum;
            temperature = state[5] * temperature_scale;
            inverse_mw = 0.0;
            for (int n = 0; n < 6; ++n)
                inverse_mw += mass_fraction[n] /
                    Model::Gas::Gas::MolecularWeight(gas, n);
            rho_eos =
                pressure / (gas_constant * inverse_mw * temperature);
            advance.volume_fraction_change +=
                gas_density * (1.0 / rho_eos - 1.0 / rho_eos_old);

            remaining -= chemistry_dt;
            if (remaining <= dt * 1.0e-14) break;
            chemistry_dt =
                2.0 * chemistry_dt < remaining ?
                    2.0 * chemistry_dt : remaining;
            rejected_steps = 0;
        }

        for (int n = 0; n < 6; ++n)
            rhoY[n] = gas_density * mass_fraction[n];
        return advance;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    std::pair<std::array<Set::Scalar, MAX_SPECIES>, Set::Scalar>
    evaluate_source(
        const Set::Scalar P,
        const Set::Scalar T,
        const std::array<Set::Scalar, MAX_SPECIES>& rhoY,
        const Set::Scalar dt) const
    {
        std::array<Set::Scalar, nreactions> progress = compute_progress_rates(P, T, rhoY);
        std::pair<std::array<Set::Scalar, MAX_SPECIES>, Set::Scalar> source = source_from_progress(progress);

        if (dt > 0.0)
        {
            Set::Scalar limiter = 1.0;
            for (int n = 0; n < MAX_SPECIES; ++n)
            {
                if (source.first[n] < 0.0)
                {
                    const Set::Scalar available = std::max(rhoY[n], 0.0);
                    const Set::Scalar requested = -source.first[n] * dt;
                    if (requested > available && requested > 0.0)
                    {
                        limiter = std::min(limiter, available / requested);
                    }
                }
            }

            if (limiter < 1.0)
            {
                for (int n = 0; n < MAX_SPECIES; ++n) source.first[n] *= limiter;
                source.second *= limiter;
            }
        }

        return source;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    std::pair<std::array<Set::Scalar,MAX_SPECIES>,Set::Scalar>
    ComputeChemistrySources(
        const Set::Scalar P,
        const Set::Scalar T,
        const std::array<Set::Scalar, MAX_SPECIES>& rhoY,
        const Set::Scalar dt,
        const Model::Gas::Gas* /*gas*/) const
    {
        SpeciesArray wdot{};
        for (int n = 0; n < MAX_SPECIES; ++n)
        {
            wdot[n] = 0.0;
        }

        if (nspecies != 6)
        {
            Util::Abort(INFO, "GrossModel chemistry requires six gas species.");
            return {wdot, 0.0};
        }

        {
            SpeciesArray rhoY_valid{};
            Set::Scalar density_scale = 0.0;
            for (int n = 0; n < MAX_SPECIES; ++n)
            {
                if (!std::isfinite(rhoY[n]))
                {
                    Util::Abort(INFO, "GrossModel chemistry received non-finite species density for species ",
                        n, ": ", rhoY[n]);
                }
                density_scale += std::abs(rhoY[n]);
            }

            const Set::Scalar density_tol = 1.0e-12 * ((density_scale > 1.0) ? density_scale : 1.0);
            for (int n = 0; n < MAX_SPECIES; ++n)
            {
                rhoY_valid[n] = rhoY[n];
                if (rhoY_valid[n] < 0.0)
                {
                    if (rhoY_valid[n] < -density_tol)
                    {
                        Util::Abort(INFO, "GrossModel chemistry received negative species density for species ",
                            n, ": rhoY=", rhoY_valid[n], ", tolerance=", density_tol);
                    }
                    rhoY_valid[n] = 0.0;
                }
            }

            if (!(sum_species(rhoY_valid) > 0.0)) return {wdot, 0.0};

            if (integration == IntegrationType::SUBSTEP && dt > 0.0 && nsubsteps > 1)
            {
                const Set::Scalar dt_sub = dt / nsubsteps;
                std::array<Set::Scalar, MAX_SPECIES> rhoY_sub = rhoY_valid;
                Set::Scalar qdot_integral = 0.0;

                for (int s = 0; s < nsubsteps; ++s)
                {
                    std::pair<std::array<Set::Scalar, MAX_SPECIES>, Set::Scalar> source =
                        evaluate_source(P, T, rhoY_sub, dt_sub);
                    for (int n = 0; n < MAX_SPECIES; ++n)
                    {
                        rhoY_sub[n] += source.first[n] * dt_sub;
                        if (rhoY_sub[n] < 0.0) rhoY_sub[n] = 0.0;
                    }
                    qdot_integral += source.second * dt_sub;
                }

                for (int n = 0; n < MAX_SPECIES; ++n)
                {
                    wdot[n] = (rhoY_sub[n] - rhoY_valid[n]) / dt;
                }

                return {wdot, qdot_integral / dt};
            }

            return evaluate_source(P, T, rhoY_valid, dt);
        }
    }

    static void Parse(GrossModel& value, IO::ParmParse& pp, int a_nspecies)
    {
        value.nspecies = a_nspecies;
        value.gas_constant = Set::Constant::Rg;
        Util::Message(INFO, pp.getPrefix());
        Util::Assert(INFO, TEST(value.nspecies == 6),
            "GrossModel chemistry requires six gas species.");

        std::vector<Set::Scalar> A;
        std::vector<Set::Scalar> pressure_exponent;
        std::vector<Set::Scalar> activation_energy;
        std::vector<Set::Scalar> qgas;

        std::string integration_type = "explicit";
        pp.query("integration", integration_type);
        if (integration_type == "explicit")
            value.integration = IntegrationType::EXPLICIT;
        else if (integration_type == "substep")
            value.integration = IntegrationType::SUBSTEP;
        else if (integration_type == "implicit")
            value.integration = IntegrationType::IMPLICIT;
        else
            Util::Abort(INFO, "Invalid GrossModel integration type: ", integration_type);

        pp.queryarr_default("A", A,
            {"27.0e3_g/cm^3/s", "25.3e3_g/cm^3/s", "5.0e3_g/cm^3/s", "0.4e3_g/cm^3/s"},
            Unit::Density() / Unit::Time());
        pp.queryarr_default("pressure_exponent", pressure_exponent, {1.6, 1.7, 1.7, 1.7});
        pp.queryarr_default("activation_energy", activation_energy,
            {"15.9_kcal/mol", "14.5_kcal/mol", "11.3_kcal/mol", "14.5_kcal/mol"},
            Unit::Energy() / Unit::Amount());
        pp.queryarr_default("qgas", qgas,
            {"430.0_cal/g", "447.0_cal/g", "8365.0_cal/g", "2127.0_cal/g"},
            Unit::Energy() / Unit::Mass());
        // Reference pressure (default is 1 bar)
        pp.query_default("pressure_reference", value.pressure_reference, "1_bar", Unit::Pressure());
        // Number of explicit chemistry substeps
        pp.query_default("nsubsteps", value.nsubsteps, 1);
        pp.query_default("max_iter", value.max_iter, 100);
        pp.query_default("reltol", value.reltol, 1.0e-8);
        pp.query_default("abstol", value.abstol, 1.0e-12);
        // Stoichiometric coefficient for beta*AP + HTPB -> (beta+1) Primary
        pp.query_default("beta", value.beta, 7.974);
        // Stoichiometric coefficient for gamma*Mono + premixed -> (gamma+1) Final
        pp.query_default("gamma", value.gamma, 3.19);

        Util::Assert(INFO, TEST(A.size() == nreactions), "chemistry.model.gross_model.A must have 4 entries.");
        Util::Assert(INFO, TEST(pressure_exponent.size() == nreactions),
            "chemistry.model.gross_model.pressure_exponent must have 4 entries.");
        Util::Assert(INFO, TEST(activation_energy.size() == nreactions),
            "chemistry.model.gross_model.activation_energy must have 4 entries.");
        Util::Assert(INFO, TEST(qgas.size() == nreactions), "chemistry.model.gross_model.qgas must have 4 entries.");
        Util::Assert(INFO, TEST(value.pressure_reference > 0.0),
            "chemistry.model.gross_model.pressure_reference must be positive.");
        Util::Assert(INFO, TEST(value.nsubsteps > 0), "chemistry.model.gross_model.nsubsteps must be > 0.");
        Util::Assert(INFO, TEST(value.max_iter > 0), "chemistry.model.gross_model.max_iter must be > 0.");
        Util::Assert(INFO, TEST(value.reltol >= 0.0), "chemistry.model.gross_model.reltol must be nonnegative.");
        Util::Assert(INFO, TEST(value.abstol >= 0.0), "chemistry.model.gross_model.abstol must be nonnegative.");
        Util::Assert(INFO, TEST(value.reltol > 0.0 || value.abstol > 0.0),
            "GrossModel implicit integration requires a positive tolerance.");
        Util::Assert(INFO, TEST(value.beta > 0), "chemistry.model.gross_model.beta must be > 0.");
        Util::Assert(INFO, TEST(value.gamma > 0), "chemistry.model.gross_model.gamma must be > 0.");

        std::copy(A.begin(), A.end(), value.A.begin());
        std::copy(pressure_exponent.begin(), pressure_exponent.end(), value.pressure_exponent.begin());
        std::copy(activation_energy.begin(), activation_energy.end(), value.activation_energy.begin());
        std::copy(qgas.begin(), qgas.end(), value.qgas.begin());
    }

};

} // namespace Chemistry
} // namespace Model

#endif

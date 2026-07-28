//
// ``Model::Mechanism::PhaseChange`` transfers one material between a
// condensed phase0 representation and either one condensed phase1
// representation or a fixed-composition bundle of gas species.  Separate
// configured mechanisms form solid-liquid-gas and direct solid-gas paths.
// The sign of the phase0 rate selects the direction, and reversible transfers
// are limited by the locally available phase1 mass.
//

#ifndef MODEL_MECHANISM_PHASECHANGE_H
#define MODEL_MECHANISM_PHASECHANGE_H

#include <array>
#include <cmath>
#include <vector>

#include "Model/Chemistry/Chemistry.H"
#include "Model/Mechanism/Mechanism.H"
#include "Model/PhaseField/AllenCahn.H"
#include "Model/PhaseField/PhaseField.H"
#include "Numeric/Stencil.H"
#include "Set/Set.H"

namespace Model
{
namespace Mechanism
{
class PhaseChange
{
public:
    static constexpr const char* name = "phase_change";
    static constexpr int max_phase1_species = Model::Chemistry::MAX_SPECIES;

    static void Parse(PhaseChange& value, IO::ParmParse& pp,
                    const std::vector<std::string>& species_names,
                    const int ngas_species,
                    const std::vector<int>& rigid_solid_species,
                    const std::vector<Set::Scalar>& reference_density,
                    const std::vector<double>& molecular_weight,
                    const Set::Scalar universal_gas_constant)
    {
        std::string phase0_name;
        std::vector<std::string> phase1_names;
        std::vector<Set::Scalar> phase1_mass_fractions;
        pp.query_required("phase0", phase0_name);
        pp.queryarr_required("phase1", phase1_names);

        for (int n = 0; n < static_cast<int>(species_names.size()); ++n)
            if (species_names[n] == phase0_name) value.phase0 = n;
        if (value.phase0 < 0)
            Util::Exception(INFO, "Unknown species identifier ", phase0_name,
                            " in ", pp.full("phase0"));
        if (value.phase0 < ngas_species)
            Util::Exception(INFO, pp.full("phase0"),
                            " must identify a condensed species");
        for (int n = 0; n < static_cast<int>(rigid_solid_species.size()); ++n)
            if (value.phase0 == rigid_solid_species[n])
            {
                value.phase0_is_rigid = true;
                value.rigid_component = n;
            }
        value.phase0_reference_density = reference_density[value.phase0];
        if (!(value.phase0_reference_density > 0.0))
            Util::Exception(INFO, "No positive reference density is defined for ",
                            phase0_name);

        if (phase1_names.empty() ||
            phase1_names.size() > max_phase1_species)
            Util::Exception(INFO, pp.full("phase1"),
                            " must contain between 1 and ",
                            max_phase1_species, " species");

        for (int n = 0; n < static_cast<int>(phase1_names.size()); ++n)
        {
            value.phase1_species[n] = -1;
            for (int m = 0; m < static_cast<int>(species_names.size()); ++m)
                if (species_names[m] == phase1_names[n])
                    value.phase1_species[n] = m;
            if (value.phase1_species[n] < 0)
                Util::Exception(INFO, "Unknown species identifier ",
                                phase1_names[n], " in ", pp.full("phase1"));
            if (value.phase1_species[n] == value.phase0)
                Util::Exception(INFO, phase0_name,
                                " cannot appear in both phase0 and phase1");
        }

        value.phase1_is_gas = value.phase1_species[0] < ngas_species;
        if (!value.phase1_is_gas && phase1_names.size() != 1)
            Util::Exception(INFO,
                            "A condensed phase1 must contain exactly one species");
        for (int n = 0; n < static_cast<int>(phase1_names.size()); ++n)
            if ((value.phase1_species[n] < ngas_species) != value.phase1_is_gas)
                Util::Exception(INFO,
                                "phase1 cannot mix gas and condensed species");

        if (value.phase1_is_gas)
        {
            if (pp.contains("phase1.mass_fractions"))
                pp.queryarr_required("phase1.mass_fractions",
                                    phase1_mass_fractions);
            else if (phase1_names.size() == 1)
                phase1_mass_fractions.push_back(1.0);
            else
                Util::Exception(INFO, pp.full("phase1.mass_fractions"),
                                " is required for a multi-species gas phase");
            if (phase1_mass_fractions.size() != phase1_names.size())
                Util::Exception(INFO, pp.full("phase1.mass_fractions"),
                                " must contain one value per gas species");

            Set::Scalar sum = 0.0;
            for (int n = 0; n < static_cast<int>(phase1_names.size()); ++n)
            {
                if (phase1_mass_fractions[n] < 0.0)
                    Util::Exception(INFO, pp.full("phase1.mass_fractions"),
                                    " must be nonnegative");
                if (value.phase1_species[n] >=
                    static_cast<int>(molecular_weight.size()) ||
                    !(molecular_weight[value.phase1_species[n]] > 0.0))
                    Util::Exception(INFO, phase1_names[n],
                                    " requires a positive molecular weight");
                value.phase1_mass_fractions[n] = phase1_mass_fractions[n];
                value.phase1_gas_constant += universal_gas_constant *
                    phase1_mass_fractions[n] /
                    molecular_weight[value.phase1_species[n]];
                sum += phase1_mass_fractions[n];
            }
            if (Util::Abs(sum - 1.0) > 1.0e-12)
                Util::Exception(INFO, pp.full("phase1.mass_fractions"),
                                " must sum to one");
        }
        else
        {
            if (pp.contains("phase1.mass_fractions"))
                Util::Exception(INFO, pp.full("phase1.mass_fractions"),
                                " is only valid for a gas phase1");
            const int n = value.phase1_species[0];
            if (!(reference_density[n] > 0.0))
                Util::Exception(INFO, "No positive reference density is defined for ",
                                phase1_names[0]);
            value.phase1_mass_fractions[0] = 1.0;
            value.phase1_inverse_reference_density = 1.0 / reference_density[n];
        }
        value.nphase1_species = static_cast<int>(phase1_names.size());

        pp.query_default("reversible", value.reversible, false);
        pp.query_default("latent_heat", value.latent_heat,
                        "0.0_J/kg", Unit::Energy() / Unit::Mass());
        if (!std::isfinite(value.latent_heat))
            Util::Exception(INFO, "latent_heat must be finite");
        pp.query_default("rate_multiplier", value.rate_multiplier, 1.0);
        if (value.rate_multiplier < 0.0)
            Util::Exception(INFO, "rate_multiplier must be nonnegative");
        pp.query_default("activation_temperature", value.activation_temperature,
                        "0.0_K", Unit::Temperature());
        pp.query_default("temperature_cutoff", value.temperature_cutoff,
                        "0.0_K", Unit::Temperature());
        if (value.activation_temperature < 0.0 ||
            value.temperature_cutoff < 0.0)
            Util::Exception(INFO,
                            "phase-change temperatures must be nonnegative");
        value.use_prescribed_speed = pp.contains("prescribed_speed");
        if (value.use_prescribed_speed)
        {
            pp.query_required("prescribed_speed", value.prescribed_speed,
                            Unit::Velocity());
            if (value.prescribed_speed < 0.0)
                Util::Exception(INFO,
                    "phase-change prescribed_speed must be nonnegative");
        }
        pp.select<Model::PhaseField::AllenCahn>("phase_field", value.model);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar VolumeSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Set::Scalar source = MassSource(state, i, j, k, DX);
        Set::Scalar phase1_specific_volume =
            phase1_inverse_reference_density;
        if (phase1_is_gas)
        {
            if (!(state.pressure > 0.0)) return 0.0;
            phase1_specific_volume = phase1_gas_constant *
                state.temperature / state.pressure;
        }
        return source *
            (1.0 / phase0_reference_density - phase1_specific_volume);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Apply(const amrex::Array4<Set::Scalar>& component_density_rhs,
                const State& state,
                const int i, const int j, const int k,
                const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Set::Scalar source = MassSource(state, i, j, k, DX);
        component_density_rhs(i,j,k,phase0) += source;
        for (int n = 0; n < nphase1_species; ++n)
            component_density_rhs(i,j,k,phase1_species[n]) -=
                phase1_mass_fractions[n] * source;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar HeatSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        return latent_heat * MassSource(state, i, j, k, DX);
    }

    Set::Scalar StabilityRate(const Set::Scalar dx,
                            const Set::Scalar temperature) const
    {
        if (!(temperature > temperature_cutoff)) return 0.0;
        if (use_prescribed_speed) return prescribed_speed / dx;
        const Set::Scalar temperature_factor = activation_temperature > 0.0 ?
            std::exp(-activation_temperature / temperature) : 1.0;
        return rate_multiplier * temperature_factor *
            (phase0_is_rigid ? model.LocalStabilityRate() :
                              model.StabilityRate(dx));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar LocalRate(
        const State& state, const int i, const int j, const int k) const
    {
        if (!phase0_is_rigid || use_prescribed_speed ||
            !(state.temperature > temperature_cutoff))
            return 0.0;

        const Set::Scalar eta =
            Util::Clamp(state.rigid_eta(i,j,k), 0.0, 1.0);
        if (!(eta > 0.0)) return 0.0;
        const Set::Scalar phase0_eta = Util::Clamp(
            state.rigid_species_eta(i,j,k,rigid_component), 0.0, eta);
        const Set::Scalar composition =
            Model::PhaseField::H(phase0_eta / eta);
        Set::Scalar rate = composition *
            model.Local(state.rigid_eta, 0, 1.0, i, j, k);
        if (activation_temperature > 0.0)
            rate *= std::exp(
                -activation_temperature / state.temperature);
        return rate_multiplier * rate;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar GradientCoefficient(
        const State& state, const int i, const int j, const int k) const
    {
        if (!phase0_is_rigid || use_prescribed_speed ||
            !(state.temperature > temperature_cutoff))
            return 0.0;

        const Set::Scalar eta =
            Util::Clamp(state.rigid_eta(i,j,k), 0.0, 1.0);
        if (!(eta > 0.0)) return 0.0;
        const Set::Scalar phase0_eta = Util::Clamp(
            state.rigid_species_eta(i,j,k,rigid_component), 0.0, eta);
        const Set::Scalar composition =
            Model::PhaseField::H(phase0_eta / eta);
        const Set::Scalar temperature_factor =
            activation_temperature > 0.0 ?
                std::exp(-activation_temperature / state.temperature) : 1.0;
        return composition * rate_multiplier * temperature_factor *
            model.GradientCoefficient();
    }

    int RigidComponent() const
    {
        return phase0_is_rigid && !use_prescribed_speed ?
            rigid_component : -1;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar ApplyImplicitChange(
        const amrex::Array4<Set::Scalar>& component_density,
        const State& state, const Set::Scalar eta_change,
        Set::Scalar& integrated_heat,
        const int i, const int j, const int k) const
    {
        Set::Scalar mass_change = phase0_reference_density * eta_change;
        if (!reversible) mass_change = Util::Min(mass_change, 0.0);
        mass_change = Util::Max(mass_change,
            -Util::Max(component_density(i,j,k,phase0), 0.0));
        if (mass_change > 0.0)
        {
            Set::Scalar available = 1.0e300;
            for (int n = 0; n < nphase1_species; ++n)
                if (phase1_mass_fractions[n] > 0.0)
                    available = Util::Min(available,
                        Util::Max(component_density(
                            i,j,k,phase1_species[n]), 0.0) /
                            phase1_mass_fractions[n]);
            mass_change = Util::Min(mass_change, available);
        }

        component_density(i,j,k,phase0) += mass_change;
        for (int n = 0; n < nphase1_species; ++n)
            component_density(i,j,k,phase1_species[n]) -=
                phase1_mass_fractions[n] * mass_change;
        integrated_heat += latent_heat * mass_change;

        Set::Scalar phase1_specific_volume =
            phase1_inverse_reference_density;
        if (phase1_is_gas)
        {
            if (!(state.pressure > 0.0)) return 0.0;
            phase1_specific_volume = phase1_gas_constant *
                state.temperature / state.pressure;
        }
        return mass_change *
            (1.0 / phase0_reference_density - phase1_specific_volume);
    }

private:
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar MassSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if (!(state.temperature > temperature_cutoff)) return 0.0;

        Set::Scalar eta_dot = 0.0;
        if (phase0_is_rigid)
        {
            const Set::Scalar eta =
                Util::Clamp(state.rigid_eta(i,j,k), 0.0, 1.0);
            if (eta > 0.0 && use_prescribed_speed)
            {
                const Set::Scalar phase0_eta = Util::Clamp(
                    state.rigid_species_eta(i,j,k,rigid_component),
                    0.0, eta);
                eta_dot = -phase0_eta / eta * prescribed_speed *
                    Numeric::Gradient(
                        state.rigid_eta, i, j, k, 0, DX).norm();
            }
        }
        else
        {
            if (use_prescribed_speed)
                eta_dot = -prescribed_speed / phase0_reference_density *
                    Numeric::Gradient(
                        state.component_density, i, j, k, phase0, DX).norm();
            else
                eta_dot = model(state.component_density, phase0,
                                1.0 / phase0_reference_density,
                                i, j, k, DX);
        }
        if (activation_temperature > 0.0)
            eta_dot *= std::exp(-activation_temperature / state.temperature);
        Set::Scalar source =
            rate_multiplier * phase0_reference_density * eta_dot;
        if (!reversible) source = Util::Min(source, 0.0);

        if (state.dt > 0.0)
        {
            source = Util::Max(source,
                -Util::Max(state.component_density(i,j,k,phase0), 0.0) /
                    state.dt);
            if (source > 0.0)
            {
                Set::Scalar available = 1.0e300;
                for (int n = 0; n < nphase1_species; ++n)
                    if (phase1_mass_fractions[n] > 0.0)
                        available = Util::Min(available,
                            Util::Max(state.component_density(
                                i,j,k,phase1_species[n]), 0.0) /
                                phase1_mass_fractions[n]);
                source = Util::Min(source, available / state.dt);
            }
        }
        return source;
    }

    Model::PhaseField::PhaseField<Model::PhaseField::AllenCahn> model;
    int phase0 = -1;
    int rigid_component = -1;
    int nphase1_species = 0;
    std::array<int, max_phase1_species> phase1_species{};
    std::array<Set::Scalar, max_phase1_species> phase1_mass_fractions{};
    Set::Scalar phase0_reference_density = NAN;
    Set::Scalar phase1_inverse_reference_density = 0.0;
    Set::Scalar phase1_gas_constant = 0.0;
    Set::Scalar latent_heat = 0.0;
    Set::Scalar rate_multiplier = 1.0;
    Set::Scalar activation_temperature = 0.0;
    Set::Scalar temperature_cutoff = 0.0;
    Set::Scalar prescribed_speed = 0.0;
    bool phase0_is_rigid = false;
    bool phase1_is_gas = false;
    bool reversible = false;
    bool use_prescribed_speed = false;
};
}
}

#endif

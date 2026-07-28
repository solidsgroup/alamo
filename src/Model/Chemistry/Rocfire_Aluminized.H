//
// Rocfire_Aluminized extends the six-species AP/HTPB Rocfire mechanism with
// two inert aluminum-bearing gas species.  The required gas ordering is
//
//   AP, binder, monopropellant, premixed, primary, final, Al_gas, Al2O3_gas.
//
// The first six species retain the original Rocfire reactions and heat
// release.  Al_gas and Al2O3_gas participate in mixture thermodynamics and
// dilution but have zero homogeneous chemistry source until aluminum kinetics
// are added.
//

#ifndef MODEL_CHEMISTRY_ROCFIRE_ALUMINIZED_H
#define MODEL_CHEMISTRY_ROCFIRE_ALUMINIZED_H

#include "Model/Chemistry/Rocfire.H"

namespace Model
{
namespace Chemistry
{

struct Rocfire_Aluminized : public Rocfire
{
    static constexpr const char* name = "rocfire_aluminized";
    static constexpr bool reactive = true;
    static constexpr bool has_jacobian = true;
    static constexpr int aluminum_gas = 6;
    static constexpr int aluminum_oxide_gas = 7;
    static constexpr int required_species = 8;

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    std::pair<SpeciesArray, Set::Scalar>
    ComputeChemistrySources(
        const Set::Scalar pressure,
        const Set::Scalar temperature,
        const SpeciesArray& partial_density,
        const Set::Scalar dt,
        const Model::Gas::Gas* /*gas*/) const
    {
        SpeciesArray source{};
        if (nspecies != required_species)
        {
            Util::Abort(INFO,
                "Rocfire_Aluminized chemistry requires eight gas species.");
            return {source, 0.0};
        }

        SpeciesArray valid_density{};
        Set::Scalar density_scale = 0.0;
        for (int n = 0; n < nspecies; ++n)
        {
            if (!std::isfinite(partial_density[n]))
                Util::Abort(INFO,
                    "Rocfire_Aluminized chemistry received non-finite species "
                    "density for species ", n, ": ", partial_density[n]);
            density_scale += std::abs(partial_density[n]);
        }

        const Set::Scalar density_tolerance = 1.0e-12 *
            ((density_scale > 1.0) ? density_scale : 1.0);
        for (int n = 0; n < nspecies; ++n)
        {
            valid_density[n] = partial_density[n];
            if (valid_density[n] < 0.0)
            {
                if (valid_density[n] < -density_tolerance)
                    Util::Abort(INFO,
                        "Rocfire_Aluminized chemistry received negative species "
                        "density for species ", n, ": rhoY=", valid_density[n],
                        ", tolerance=", density_tolerance);
                valid_density[n] = 0.0;
            }
        }

        if (!(sum_species(valid_density) > 0.0)) return {source, 0.0};

        if (integration == IntegrationType::SUBSTEP && dt > 0.0 && nsubsteps > 1)
        {
            const Set::Scalar substep_dt = dt / nsubsteps;
            SpeciesArray substep_density = valid_density;
            Set::Scalar heat_release_integral = 0.0;
            for (int step = 0; step < nsubsteps; ++step)
            {
                const auto substep_source = evaluate_source(
                    pressure, temperature, substep_density, substep_dt);
                for (int n = 0; n < nspecies; ++n)
                    substep_density[n] = std::max(
                        substep_density[n] +
                            substep_source.first[n] * substep_dt,
                        0.0);
                heat_release_integral += substep_source.second * substep_dt;
            }
            for (int n = 0; n < nspecies; ++n)
                source[n] = (substep_density[n] - valid_density[n]) / dt;
            return {source, heat_release_integral / dt};
        }

        // Rocfire::evaluate_source touches only indices 0--5.  Both aluminum
        // gas sources therefore remain exactly zero while their densities
        // dilute the reacting mass fractions used by the inherited rates.
        return evaluate_source(pressure, temperature, valid_density, dt);
    }

    template<typename Gas>
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
        const Gas* gas) const
    {
        if (gas == nullptr || ode_size != required_species ||
            nspecies != required_species || dependent_species < 0 ||
            dependent_species >= nspecies || !(pressure > 0.0) ||
            !(mixture_density > 0.0) || !(gas_density > 0.0) ||
            !(temperature_scale > 0.0) || !(pressure_reference > 0.0) ||
            !(gas->Rg > 0.0)) return false;

        SpeciesArray mass_fraction{};
        int independent_species[MAX_SPECIES]{};
        Set::Scalar independent_sum = 0.0;
        for (int column = 0; column < nspecies - 1; ++column)
        {
            const int species = column < dependent_species ? column : column + 1;
            independent_species[column] = species;
            mass_fraction[species] = state[column];
            if (!(mass_fraction[species] >= 0.0) ||
                !std::isfinite(mass_fraction[species])) return false;
            independent_sum += mass_fraction[species];
        }
        mass_fraction[dependent_species] = 1.0 - independent_sum;
        if (!(mass_fraction[dependent_species] >= 0.0) ||
            !std::isfinite(mass_fraction[dependent_species])) return false;

        const Set::Scalar temperature =
            state[nspecies - 1] * temperature_scale;
        if (!(temperature > 0.0) || !std::isfinite(temperature)) return false;

        Set::Scalar inverse_molecular_weight = 0.0;
        Set::Scalar heat_capacity = 0.0;
        Set::Scalar species_heat_capacity[MAX_SPECIES]{};
        Set::Scalar species_heat_capacity_derivative[MAX_SPECIES]{};
        for (int n = 0; n < nspecies; ++n)
        {
            if (!(gas->MW[n] > 0.0)) return false;
            species_heat_capacity[n] =
                gas->cp_mol_species(temperature, n) / gas->MW[n];
            species_heat_capacity_derivative[n] =
                gas->dcp_mol_species_dT(temperature, n) / gas->MW[n];
            inverse_molecular_weight += mass_fraction[n] / gas->MW[n];
            heat_capacity += mass_fraction[n] * species_heat_capacity[n];
        }
        if (!(inverse_molecular_weight > 0.0) || !(heat_capacity > 0.0))
            return false;

        std::array<Set::Scalar, nreactions> rate_constant{};
        std::array<Set::Scalar, nreactions> progress{};
        for (int reaction = 0; reaction < nreactions; ++reaction)
            rate_constant[reaction] = A[reaction] *
                std::pow(pressure / pressure_reference,
                    pressure_exponent[reaction]) *
                std::exp(-activation_energy[reaction] /
                    (Set::Constant::Rg * temperature));

        progress[0] = rate_constant[0] * mass_fraction[AP];
        progress[1] = rate_constant[1] * mass_fraction[Binder];
        progress[2] = rate_constant[2] *
            mass_fraction[AP] * mass_fraction[Binder];
        progress[3] = rate_constant[3] *
            mass_fraction[Mono] * mass_fraction[Premixed];
        const auto chemistry_source = source_from_progress(progress);

        const Set::Scalar inverse_density = gas->Rg *
            inverse_molecular_weight * temperature / pressure;
        const Set::Scalar temperature_prefactor =
            gas_density / (mixture_density * temperature_scale);

        for (int column = 0; column < nspecies; ++column)
        {
            SpeciesArray mass_fraction_derivative{};
            Set::Scalar temperature_derivative = 0.0;
            if (column < nspecies - 1)
            {
                mass_fraction_derivative[independent_species[column]] = 1.0;
                mass_fraction_derivative[dependent_species] = -1.0;
            }
            else
                temperature_derivative = temperature_scale;

            std::array<Set::Scalar, nreactions> progress_derivative{};
            const Set::Scalar temperature_factor = temperature_derivative /
                (Set::Constant::Rg * temperature * temperature);
            progress_derivative[0] = rate_constant[0] *
                mass_fraction_derivative[AP] + progress[0] *
                activation_energy[0] * temperature_factor;
            progress_derivative[1] = rate_constant[1] *
                mass_fraction_derivative[Binder] + progress[1] *
                activation_energy[1] * temperature_factor;
            progress_derivative[2] = rate_constant[2] *
                (mass_fraction[Binder] * mass_fraction_derivative[AP] +
                 mass_fraction[AP] * mass_fraction_derivative[Binder]) +
                progress[2] * activation_energy[2] * temperature_factor;
            progress_derivative[3] = rate_constant[3] *
                (mass_fraction[Premixed] * mass_fraction_derivative[Mono] +
                 mass_fraction[Mono] * mass_fraction_derivative[Premixed]) +
                progress[3] * activation_energy[3] * temperature_factor;
            const auto source_derivative =
                source_from_progress(progress_derivative);

            Set::Scalar inverse_molecular_weight_derivative = 0.0;
            Set::Scalar heat_capacity_derivative = 0.0;
            for (int n = 0; n < nspecies; ++n)
            {
                inverse_molecular_weight_derivative +=
                    mass_fraction_derivative[n] / gas->MW[n];
                heat_capacity_derivative +=
                    mass_fraction_derivative[n] * species_heat_capacity[n] +
                    mass_fraction[n] * species_heat_capacity_derivative[n] *
                        temperature_derivative;
            }
            const Set::Scalar inverse_density_derivative = gas->Rg / pressure *
                (temperature * inverse_molecular_weight_derivative +
                 inverse_molecular_weight * temperature_derivative);

            for (int row = 0; row < nspecies - 1; ++row)
            {
                const int species = independent_species[row];
                derivative[row * nspecies + column] =
                    inverse_density_derivative *
                        chemistry_source.first[species] +
                    inverse_density * source_derivative.first[species];
            }
            derivative[(nspecies - 1) * nspecies + column] =
                temperature_prefactor *
                ((inverse_density_derivative * chemistry_source.second +
                  inverse_density * source_derivative.second) /
                     heat_capacity -
                 inverse_density * chemistry_source.second *
                     heat_capacity_derivative /
                     (heat_capacity * heat_capacity));
        }

        for (int n = 0; n < nspecies * nspecies; ++n)
            if (!std::isfinite(derivative[n])) return false;
        return true;
    }

    static void Parse(
        Rocfire_Aluminized& value, IO::ParmParse& pp, int gas_species)
    {
        Util::Assert(INFO, TEST(gas_species == required_species),
            "Rocfire_Aluminized chemistry requires eight gas species ordered "
            "as AP, binder, monopropellant, premixed, primary, final, "
            "Al_gas, Al2O3_gas.");

        // Reuse the original model's parser without weakening its six-species
        // contract.  Parameters are read from the rocfire_aluminized prefix;
        // restore the aluminized species count afterward.
        Rocfire::Parse(static_cast<Rocfire&>(value), pp, 6);
        value.nspecies = gas_species;
    }
};

}
}

#endif

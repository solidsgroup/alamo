//
// ``Model::Mechanism::PhaseChange`` couples a phase-field rate for one
// condensed input species to mass-conserving production of one or more gas
// species. Product mass fractions also determine the ideal-gas dilation sent
// to the low-Mach projection.
//

#ifndef MODEL_MECHANISM_PHASECHANGE_H
#define MODEL_MECHANISM_PHASECHANGE_H

#include <array>
#include <vector>

#include "Model/Mechanism/Mechanism.H"
#include "Model/PhaseField/AllenCahn.H"
#include "Model/PhaseField/PhaseField.H"
#include "Set/Set.H"

namespace Model
{
namespace Mechanism
{
class PhaseChange
{
public:
    static constexpr const char* name = "phase_change";
    static constexpr int max_products = 16;

    static void Parse(PhaseChange& value, IO::ParmParse& pp,
                    const std::vector<std::string>& species_names,
                    const int ngas_species,
                    const std::vector<int>& rigid_solid_species,
                    const std::vector<Set::Scalar>& reference_density,
                    const std::vector<double>& molecular_weight,
                    const Set::Scalar universal_gas_constant)
    {
        std::string input;
        std::vector<std::string> outputs;
        std::vector<Set::Scalar> mass_fractions;
        pp.query_required("in", input);
        pp.queryarr_required("out", outputs);

        for (int n = 0; n < static_cast<int>(species_names.size()); ++n)
            if (species_names[n] == input) value.input = n;
        if (value.input < 0)
            Util::Exception(INFO, "Unknown species identifier ", input,
                            " in ", pp.full("in"));
        if (value.input < ngas_species)
            Util::Exception(INFO, pp.full("in"), " must identify a condensed species");
        for (const int n : rigid_solid_species)
            if (value.input == n) value.input_is_rigid = true;
        value.reference_density = reference_density[value.input];
        if (!(value.reference_density > 0.0))
            Util::Exception(INFO, "No positive reference density is defined for ", input);

        if (outputs.empty() || outputs.size() > max_products)
            Util::Exception(INFO, pp.full("out"), " must contain between 1 and ",
                            max_products, " gas species");
        if (pp.contains("out.mass_fractions"))
            pp.queryarr_required("out.mass_fractions", mass_fractions);
        else if (outputs.size() == 1)
            mass_fractions.push_back(1.0);
        else
            Util::Exception(INFO, pp.full("out.mass_fractions"),
                            " is required for multiple products");
        if (mass_fractions.size() != outputs.size())
            Util::Exception(INFO, pp.full("out.mass_fractions"),
                            " must contain one value per output species");

        value.nproducts = static_cast<int>(outputs.size());
        Set::Scalar sum = 0.0;
        for (int n = 0; n < value.nproducts; ++n)
        {
            value.products[n] = -1;
            for (int m = 0; m < static_cast<int>(species_names.size()); ++m)
                if (species_names[m] == outputs[n]) value.products[n] = m;
            if (value.products[n] < 0)
                Util::Exception(INFO, "Unknown species identifier ", outputs[n],
                                " in ", pp.full("out"));
            if (value.products[n] >= ngas_species ||
                value.products[n] >= static_cast<int>(molecular_weight.size()))
                Util::Exception(INFO, outputs[n], " must identify a gas species");
            if (mass_fractions[n] < 0.0)
                Util::Exception(INFO, pp.full("out.mass_fractions"),
                                " must be nonnegative");
            value.mass_fractions[n] = mass_fractions[n];
            value.product_gas_constant +=
                universal_gas_constant * mass_fractions[n] /
                molecular_weight[value.products[n]];
            sum += mass_fractions[n];
        }
        if (Util::Abs(sum - 1.0) > 1.0e-12)
            Util::Exception(INFO, pp.full("out.mass_fractions"), " must sum to one");

        pp.query_default("interface_only", value.interface_only, true);
        pp.query_default("rate_multiplier", value.rate_multiplier, 1.0);
        if (value.rate_multiplier < 0.0)
            Util::Exception(INFO, "rate_multiplier must be nonnegative");
        pp.select<Model::PhaseField::AllenCahn>("phase_field", value.model);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar VolumeSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Set::Scalar source = MassSource(state, i, j, k, DX);
        const Set::Scalar product_specific_volume =
            product_gas_constant * state.temperature / state.pressure;
        return source * (1.0 / reference_density - product_specific_volume);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Apply(const amrex::Array4<Set::Scalar>& component_density_rhs,
                const State& state,
                const int i, const int j, const int k,
                const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Set::Scalar source = MassSource(state, i, j, k, DX);
        component_density_rhs(i,j,k,input) += source;
        for (int n = 0; n < nproducts; ++n)
            component_density_rhs(i,j,k,products[n]) -= mass_fractions[n] * source;
    }

    Set::Scalar StabilityRate(const Set::Scalar dx) const
    {
        return rate_multiplier * model.StabilityRate(dx);
    }

private:
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar MassSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        Set::Scalar eta_dot = model(
            state.component_density, input, 1.0 / reference_density, i, j, k, DX);
        if (interface_only)
        {
            const Set::Scalar eta = input_is_rigid ?
                Util::Clamp(state.rigid_eta(i,j,k), 0.0, 1.0) :
                Util::Clamp(state.component_density(i,j,k,input) / reference_density,
                            0.0, 1.0);
            eta_dot *= 4.0 * eta * (1.0 - eta);
        }
        return rate_multiplier * reference_density * eta_dot;
    }

    Model::PhaseField::PhaseField<Model::PhaseField::AllenCahn> model;
    int input = -1;
    int nproducts = 0;
    std::array<int, max_products> products{};
    std::array<Set::Scalar, max_products> mass_fractions{};
    Set::Scalar reference_density = NAN;
    Set::Scalar product_gas_constant = 0.0;
    Set::Scalar rate_multiplier = 1.0;
    bool input_is_rigid = false;
    bool interface_only = true;
};
}
}

#endif

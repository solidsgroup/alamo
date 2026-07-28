//
// Surface oxidation of a condensed aluminum phase by effective oxidizer mass
// carried by one or more lumped gas species.  The mechanism is localized to
// the exposed aluminum/gas interface.  Oxygen-bearing donor mass and aluminum
// mass are converted conservatively into liquid and gaseous Al2O3 products.
//

#ifndef MODEL_MECHANISM_SURFACEOXIDATION_H
#define MODEL_MECHANISM_SURFACEOXIDATION_H

#include <array>
#include <vector>

#include "Model/Chemistry/Chemistry.H"
#include "Model/Mechanism/Mechanism.H"
#include "Numeric/Stencil.H"
#include "Unit/Unit.H"

namespace Model
{
namespace Mechanism
{
class SurfaceOxidation
{
public:
    static constexpr const char* name = "surface_oxidation";
    static constexpr int max_oxidizers = 4;

    static void Parse(SurfaceOxidation& value, IO::ParmParse& pp,
                    const std::vector<std::string>& species_names,
                    const int ngas_species,
                    const std::vector<int>& /*rigid_solid_species*/,
                    const std::vector<Set::Scalar>& reference_density,
                    const std::vector<double>& molecular_weight,
                    const Set::Scalar universal_gas_constant)
    {
        auto find_species = [&](const std::string& name)
        {
            for (int n = 0; n < static_cast<int>(species_names.size()); ++n)
                if (species_names[n] == name) return n;
            Util::Exception(INFO, "Unknown species identifier ", name,
                            " in ", pp.full("type"));
            return -1;
        };

        std::string input_name, liquid_name, gas_name;
        std::vector<std::string> oxidizer_names;
        std::vector<Set::Scalar> oxidizer_fractions;
        pp.query_required("in", input_name);
        pp.query_required("out.liquid", liquid_name);
        pp.query_required("out.gas", gas_name);
        pp.queryarr_required("oxidizers", oxidizer_names);
        pp.queryarr_required("oxidizers.available_mass_fractions",
                            oxidizer_fractions);

        value.input = find_species(input_name);
        value.liquid_product = find_species(liquid_name);
        value.gas_product = find_species(gas_name);
        if (value.input < ngas_species)
            Util::Exception(INFO, input_name,
                            " must identify condensed aluminum");
        if (value.liquid_product < ngas_species)
            Util::Exception(INFO, liquid_name,
                            " must identify condensed aluminum oxide");
        if (value.gas_product < 0 || value.gas_product >= ngas_species)
            Util::Exception(INFO, gas_name,
                            " must identify gaseous aluminum oxide");
        if (oxidizer_names.empty() ||
            oxidizer_names.size() > max_oxidizers ||
            oxidizer_fractions.size() != oxidizer_names.size())
            Util::Exception(INFO,
                "surface_oxidation requires one available mass fraction per "
                "oxidizer (maximum ", max_oxidizers, ")");

        value.nspecies = static_cast<int>(species_names.size());
        for (int n = 0; n < value.nspecies; ++n)
            if (n >= ngas_species && reference_density[n] > 0.0)
                value.inverse_reference_density[n] =
                    1.0 / reference_density[n];
        value.input_inverse_reference_density =
            value.inverse_reference_density[value.input];
        value.liquid_product_inverse_reference_density =
            value.inverse_reference_density[value.liquid_product];
        if (!(value.input_inverse_reference_density > 0.0) ||
            !(value.liquid_product_inverse_reference_density > 0.0))
            Util::Exception(INFO,
                "surface_oxidation condensed phases require positive "
                "reference densities");

        value.noxidizers = static_cast<int>(oxidizer_names.size());
        for (int n = 0; n < value.noxidizers; ++n)
        {
            value.oxidizers[n] = find_species(oxidizer_names[n]);
            if (value.oxidizers[n] >= ngas_species)
                Util::Exception(INFO, oxidizer_names[n],
                                " must identify a gas species");
            if (oxidizer_fractions[n] < 0.0 || oxidizer_fractions[n] > 1.0)
                Util::Exception(INFO,
                    "oxidizer available mass fractions must lie in [0,1]");
            value.oxidizer_available_fraction[n] = oxidizer_fractions[n];
            value.oxidizer_gas_constant[n] = universal_gas_constant /
                molecular_weight[value.oxidizers[n]];
        }
        value.product_gas_constant = universal_gas_constant /
            molecular_weight[value.gas_product];

        pp.query_required("aluminum_mass_flux", value.aluminum_mass_flux,
                        Unit::Mass() / Unit::Area() / Unit::Time());
        pp.query_required("gas_product_fraction", value.gas_product_fraction);
        pp.query_default("temperature_cutoff", value.temperature_cutoff,
                        "0.0_K", Unit::Temperature());
        pp.query_default("heat_release.liquid", value.liquid_heat_release,
                        "0.0_J/kg", Unit::Energy() / Unit::Mass());
        pp.query_default("heat_release.gas", value.gas_heat_release,
                        "0.0_J/kg", Unit::Energy() / Unit::Mass());
        if (value.aluminum_mass_flux < 0.0 ||
            value.gas_product_fraction < 0.0 ||
            value.gas_product_fraction > 1.0 ||
            value.temperature_cutoff < 0.0)
            Util::Exception(INFO,
                "surface_oxidation rates must be nonnegative and the gas "
                "product fraction must lie in [0,1]");
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Apply(const amrex::Array4<Set::Scalar>& rhs,
                const State& state, const int i, const int j, const int k,
                const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Source source = ComputeSource(state, i, j, k, DX);
        rhs(i,j,k,input) -= source.aluminum;
        for (int n = 0; n < noxidizers; ++n)
            rhs(i,j,k,oxidizers[n]) -= source.oxidizer[n];
        const Set::Scalar oxide = source.aluminum + source.oxygen;
        rhs(i,j,k,liquid_product) +=
            (1.0 - gas_product_fraction) * oxide;
        rhs(i,j,k,gas_product) += gas_product_fraction * oxide;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar VolumeSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Source source = ComputeSource(state, i, j, k, DX);
        if (!(state.pressure > 0.0)) return 0.0;
        Set::Scalar volume_rate =
            -source.aluminum * input_inverse_reference_density;
        for (int n = 0; n < noxidizers; ++n)
            volume_rate -= source.oxidizer[n] * oxidizer_gas_constant[n] *
                state.temperature / state.pressure;
        const Set::Scalar oxide = source.aluminum + source.oxygen;
        volume_rate += (1.0 - gas_product_fraction) * oxide *
            liquid_product_inverse_reference_density;
        volume_rate += gas_product_fraction * oxide * product_gas_constant *
            state.temperature / state.pressure;
        return volume_rate;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar HeatSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Source source = ComputeSource(state, i, j, k, DX);
        const Set::Scalar oxide = source.aluminum + source.oxygen;
        return oxide * ((1.0 - gas_product_fraction) * liquid_heat_release +
                        gas_product_fraction * gas_heat_release);
    }

    Set::Scalar StabilityRate(const Set::Scalar dx,
                            const Set::Scalar temperature) const
    {
        if (!(temperature > temperature_cutoff) || !(dx > 0.0)) return 0.0;
        return aluminum_mass_flux * input_inverse_reference_density / dx;
    }

    // Surface oxidation is advanced explicitly.  These methods satisfy the
    // common mechanism interface used by the implicit rigid phase-change
    // path and deliberately opt this mechanism out of that solve.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar LocalRate(
        const State& /*state*/, const int /*i*/, const int /*j*/,
        const int /*k*/) const
    {
        return 0.0;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar GradientCoefficient(
        const State& /*state*/, const int /*i*/, const int /*j*/,
        const int /*k*/) const
    {
        return 0.0;
    }

    int RigidComponent() const { return -1; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar ApplyImplicitChange(
        const amrex::Array4<Set::Scalar>& /*component_density*/,
        const State& /*state*/, const Set::Scalar /*eta_change*/,
        const int /*i*/, const int /*j*/, const int /*k*/) const
    {
        return 0.0;
    }

private:
    struct Source
    {
        Set::Scalar aluminum = 0.0;
        Set::Scalar oxygen = 0.0;
        std::array<Set::Scalar, max_oxidizers> oxidizer{};
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Source ComputeSource(const State& state,
                        const int i, const int j, const int k,
                        const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        Source source;
        if (!(state.temperature > temperature_cutoff) ||
            !(aluminum_mass_flux > 0.0)) return source;

        const Set::Scalar aluminum_eta = Util::Clamp(
            state.component_density(i,j,k,input) *
                input_inverse_reference_density, 0.0, 1.0);
        const Set::Scalar oxide_eta = Util::Clamp(
            state.component_density(i,j,k,liquid_product) *
                liquid_product_inverse_reference_density, 0.0, 1.0);
        Set::Scalar condensed_eta = 0.0;
        for (int n = 0; n < nspecies; ++n)
            condensed_eta += Util::Max(state.component_density(i,j,k,n), 0.0) *
                inverse_reference_density[n];
        const Set::Scalar gas_eta = Util::Clamp(1.0 - condensed_eta, 0.0, 1.0);
        const Set::Scalar contacting_phase = gas_eta + oxide_eta;
        const Set::Scalar exposure = contacting_phase > 0.0 ?
            gas_eta / contacting_phase : 0.0;
        if (!(aluminum_eta > 0.0) || !(exposure > 0.0)) return source;

        const Set::Vector grad_aluminum = Numeric::Gradient(
            state.component_density, i, j, k, input, DX) *
            input_inverse_reference_density;
        Set::Scalar aluminum_rate = aluminum_mass_flux *
            grad_aluminum.norm() * exposure;

        Set::Scalar available_oxygen = 0.0;
        for (int n = 0; n < noxidizers; ++n)
            available_oxygen += oxidizer_available_fraction[n] *
                Util::Max(state.component_density(i,j,k,oxidizers[n]), 0.0);

        if (state.dt > 0.0)
        {
            aluminum_rate = Util::Min(aluminum_rate,
                Util::Max(state.component_density(i,j,k,input), 0.0) /
                    state.dt);
            aluminum_rate = Util::Min(aluminum_rate,
                available_oxygen /
                    (oxygen_per_aluminum * state.dt));
        }
        aluminum_rate = Util::Max(aluminum_rate, 0.0);
        source.aluminum = aluminum_rate;
        source.oxygen = oxygen_per_aluminum * aluminum_rate;
        if (available_oxygen > 0.0)
            for (int n = 0; n < noxidizers; ++n)
                source.oxidizer[n] = source.oxygen *
                    oxidizer_available_fraction[n] *
                    Util::Max(state.component_density(i,j,k,oxidizers[n]), 0.0) /
                    available_oxygen;
        return source;
    }

    // 3 O / 2 Al by mass for Al2O3, using standard atomic weights.
    static constexpr Set::Scalar oxygen_per_aluminum =
        (3.0 * 15.999) / (2.0 * 26.9815385);

    int input = -1;
    int liquid_product = -1;
    int gas_product = -1;
    int nspecies = 0;
    int noxidizers = 0;
    std::array<int, max_oxidizers> oxidizers{};
    std::array<Set::Scalar, max_oxidizers> oxidizer_available_fraction{};
    std::array<Set::Scalar, max_oxidizers> oxidizer_gas_constant{};
    std::array<Set::Scalar, Model::Chemistry::MAX_SPECIES>
        inverse_reference_density{};
    Set::Scalar input_inverse_reference_density = 0.0;
    Set::Scalar liquid_product_inverse_reference_density = 0.0;
    Set::Scalar product_gas_constant = 0.0;
    Set::Scalar aluminum_mass_flux = 0.0;
    Set::Scalar gas_product_fraction = 0.0;
    Set::Scalar temperature_cutoff = 0.0;
    Set::Scalar liquid_heat_release = 0.0;
    Set::Scalar gas_heat_release = 0.0;
};
}
}

#endif

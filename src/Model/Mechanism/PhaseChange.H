//
// ``Model::Mechanism::PhaseChange`` couples a phase-field rate for one
// condensed input species to mass-conserving production of one or more gas
// species. Product mass fractions also determine the ideal-gas dilation sent
// to the low-Mach projection. An optional Arrhenius factor permits thermally
// activated regression without coupling this generic mechanism to a specific
// propellant model.
//

#ifndef MODEL_MECHANISM_PHASECHANGE_H
#define MODEL_MECHANISM_PHASECHANGE_H

#include <array>
#include <cmath>
#include <vector>

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
    static constexpr int max_products = 16;

    static void Parse(PhaseChange& value, IO::ParmParse& pp,
                    const std::vector<std::string>& species_names,
                    const int ngas_species,
                    const std::vector<int>& rigid_solid_species,
                    std::vector<Set::Scalar>& reference_density,
                    const std::vector<double>& molecular_weight,
                    const Set::Scalar universal_gas_constant,
                    Set::Scalar* specific_heat = nullptr,
                    Set::Scalar* thermal_conductivity = nullptr,
                    const Set::Scalar default_heat_release = 0.0)
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
        for (int n = 0; n < static_cast<int>(rigid_solid_species.size()); ++n)
            if (value.input == rigid_solid_species[n])
            {
                value.input_is_rigid = true;
                value.rigid_component = n;
            }
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

        pp.query_default("rate_multiplier", value.rate_multiplier, 1.0);
        if (value.rate_multiplier < 0.0)
            Util::Exception(INFO, "rate_multiplier must be nonnegative");
        pp.query_default("activation_temperature", value.activation_temperature,
                        "0.0_K", Unit::Temperature());
        pp.query_default("temperature_cutoff", value.temperature_cutoff,
                        "0.0_K", Unit::Temperature());
        if (value.activation_temperature < 0.0 || value.temperature_cutoff < 0.0)
            Util::Exception(INFO, "phase-change temperatures must be nonnegative");
        value.use_prescribed_speed = pp.contains("prescribed_speed");
        if (value.use_prescribed_speed)
        {
            pp.query_required("prescribed_speed", value.prescribed_speed,
                              Unit::Velocity());
            if (!(value.prescribed_speed >= 0.0) || !std::isfinite(value.prescribed_speed))
                Util::Exception(INFO, "phase-change prescribed_speed must be finite and nonnegative");
        }
        // Energy per unit mass transferred out of the condensed phase.
        // Negative values represent endothermic decomposition. The homogeneous
        // mode below instead computes this from the two constituent heats.
        bool homogeneous;
        pp.query_default("homogeneous", homogeneous, false);
        if (!homogeneous)
        {
            const bool explicit_heat = pp.contains("heat_release");
            pp.query_default("heat_release", value.heat_release,
                             "0.0_J/kg", Unit::Energy() / Unit::Mass());
            if (!explicit_heat) value.heat_release = default_heat_release;
            if (!std::isfinite(value.heat_release))
                Util::Exception(INFO, "phase-change heat_release must be finite");
        }
        if (homogeneous)
        {
            if (pp.contains("heat_release"))
                Util::Exception(INFO, "Use homogeneous.htpb.heat_release and homogeneous.ap.heat_release for a blend");
            if (!value.input_is_rigid || specific_heat == nullptr || thermal_conductivity == nullptr)
                Util::Exception(INFO, "homogeneous phase change requires a rigid solid with thermal properties");
            if (value.nproducts != 1)
                Util::Exception(INFO, "A homogeneous binder must transfer to a single binder gas species");
            // Pure constituents are material inputs, separate from the
            // homogenized binder species and its density initial condition.
            Set::Scalar rho_htpb, rho_ap, cp_htpb, cp_ap, k_htpb, k_ap;
            pp.query_required("homogeneous.htpb.reference_density", rho_htpb, Unit::Density());
            pp.query_required("homogeneous.ap.reference_density", rho_ap, Unit::Density());
            pp.query_required("homogeneous.htpb.specific_heat", cp_htpb, Unit::SpecificHeatCapacity());
            pp.query_required("homogeneous.ap.specific_heat", cp_ap, Unit::SpecificHeatCapacity());
            pp.query_required("homogeneous.htpb.thermal_conductivity", k_htpb, Unit::ThermalConductivity());
            pp.query_required("homogeneous.ap.thermal_conductivity", k_ap, Unit::ThermalConductivity());
            for (const auto property : {rho_htpb, rho_ap, cp_htpb, cp_ap, k_htpb, k_ap})
                if (!(property > 0.0) || !std::isfinite(property))
                    Util::Exception(INFO, "Pure AP/HTPB thermal properties must be positive and finite");

            Set::Scalar w;
            if (pp.contains("homogeneous.mass_fraction"))
            {
                if (pp.contains("homogeneous.total_mass_fraction") ||
                    pp.contains("homogeneous.resolved_mass_fraction"))
                    Util::Exception(INFO, "Specify either blend mass_fraction or total/resolved mass fractions");
                pp.query_required("homogeneous.mass_fraction", w);
            }
            else
            {
                Set::Scalar total, resolved;
                pp.query_required("homogeneous.total_mass_fraction", total);
                pp.query_required("homogeneous.resolved_mass_fraction", resolved);
                if (!(resolved >= 0.0 && resolved < 1.0 && total >= resolved && total <= 1.0))
                    Util::Exception(INFO, "Require 0 <= resolved <= total <= 1 and resolved < 1");
                w = (total - resolved) / (1.0 - resolved);
            }
            if (!(w >= 0.0 && w <= 1.0))
                Util::Exception(INFO, "homogeneous.mass_fraction must lie in [0,1]");
            const int binder = value.input;
            const Set::Scalar t = w * rho_htpb /
                (w * rho_htpb + (1.0 - w) * rho_ap);

            // Add constituent volumes: rho is volume-weighted (equivalently,
            // 1/rho is mass-weighted). Specific heat remains mass-weighted.
            value.reference_density = reference_density[binder] =
                (1.0 - t) * rho_htpb + t * rho_ap;
            specific_heat[binder] = (1.0 - w) * cp_htpb + w * cp_ap;

            // Chen et al. (2002), eqs. (9)/(11): take the unsquared physical
            // branch f-x = (1-t)(1-x) f^(1/d), bounded by the pure values.
            const Set::Scalar kb = k_htpb, kap = k_ap;
            thermal_conductivity[binder] = kb;
            if (t == 1.0) thermal_conductivity[binder] = kap;
            else if (t > 0.0 && kb != kap)
            {
                const Set::Scalar scale = Util::Max(kb, kap);
                const Set::Scalar b = kb / scale, a = kap / scale;
                Set::Scalar lo = Util::Min(a, b), hi = Util::Max(a, b);
                for (int iteration = 0; iteration < 100; ++iteration)
                {
                    const Set::Scalar mid = lo + 0.5 * (hi - lo);
                    if (mid == lo || mid == hi) break;
                    const Set::Scalar residual = mid - a - (1.0 - t) * (b - a) *
                        std::pow(mid / b, 1.0 / AMREX_SPACEDIM);
                    if (residual > 0.0) hi = mid;
                    else lo = mid;
                }
                thermal_conductivity[binder] = scale * (lo + 0.5 * (hi - lo));
            }

            Set::Scalar ap_rate, ap_activation, binder_heat, ap_heat;
            pp.query_required("homogeneous.htpb.rate_multiplier", value.rate_multiplier);
            pp.query_required("homogeneous.htpb.activation_temperature", value.activation_temperature, Unit::Temperature());
            pp.query_required("homogeneous.ap.rate_multiplier", ap_rate);
            pp.query_required("homogeneous.ap.activation_temperature", ap_activation, Unit::Temperature());
            pp.query_required("homogeneous.htpb.heat_release", binder_heat, Unit::Energy() / Unit::Mass());
            pp.query_required("homogeneous.ap.heat_release", ap_heat, Unit::Energy() / Unit::Mass());
            if (!(ap_rate >= 0.0) || !std::isfinite(ap_rate) || !(value.rate_multiplier >= 0.0) || !std::isfinite(value.rate_multiplier) ||
                !(ap_activation >= 0.0) || !std::isfinite(ap_activation) ||
                !(value.activation_temperature >= 0.0) || !std::isfinite(value.activation_temperature) ||
                !std::isfinite(binder_heat) || !std::isfinite(ap_heat))
                Util::Exception(INFO, "Invalid homogeneous kinetic or heat-release parameters");
            if (t == 1.0) value.rate_multiplier = ap_rate;
            else if (t > 0.0)
                value.rate_multiplier = value.rate_multiplier == 0.0 || ap_rate == 0.0 ? 0.0 :
                    std::exp((1.0 - t) * std::log(value.rate_multiplier) + t * std::log(ap_rate));
            value.activation_temperature = (1.0 - t) * value.activation_temperature + t * ap_activation;
            value.heat_release = (1.0 - w) * binder_heat + w * ap_heat;

            // Rocfire's binder gas already represents binder plus subgrid AP
            // decomposition products. Transfer the entire blend to that gas;
            // beta and gamma account for the remaining resolved AP demand.
            // Retain its configured molecular weight for projection expansion.
        }
        pp.select<Model::PhaseField::AllenCahn>(
            "phase_field", value.phase_field);
        if (specific_heat != nullptr)
            value.reactant_specific_heat = specific_heat[value.input];
        if (!(value.reactant_specific_heat > 0.0) || !std::isfinite(value.reactant_specific_heat))
            Util::Exception(INFO, "Phase change requires a positive finite condensed specific heat for enthalpy conservation");
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar operator()(const State& state,
                           const int i, const int j, const int k,
                           const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if (use_prescribed_speed)
        {
            const Set::Scalar factor = KineticFactor(state.temperature);
            if (!(factor > 0.0)) return 0.0;
            if (!input_is_rigid)
                return -factor * prescribed_speed / reference_density *
                    Numeric::Gradient(state.component_density, i, j, k, input, DX).norm();
            const Set::Scalar eta = Util::Clamp(state.rigid_eta(i,j,k), 0.0, 1.0);
            if (!(eta > 0.0)) return 0.0;
            const Set::Scalar input_eta = Util::Clamp(
                state.rigid_species_eta(i,j,k,rigid_component), 0.0, eta);
            return -factor * input_eta / eta * prescribed_speed *
                Numeric::Gradient(state.rigid_eta, i, j, k, 0, DX).norm();
        }
        const Set::Scalar weight = Weight(state, i, j, k);
        return weight * (input_is_rigid ?
            phase_field(state.rigid_eta, 0, 1.0, i, j, k, DX) :
            phase_field(state.component_density, input,
                        1.0 / reference_density, i, j, k, DX));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar VolumeSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if (input_is_rigid && !use_prescribed_speed) return 0.0;
        const Set::Scalar source =
            reference_density * (*this)(state, i, j, k, DX);
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
        if (input_is_rigid && !use_prescribed_speed) return;
        Transfer(component_density_rhs, state,
                 (*this)(state, i, j, k, DX), i, j, k);
    }

    Set::Scalar StabilityRate(const Set::Scalar dx,
                            const Set::Scalar temperature) const
    {
        if (use_prescribed_speed)
            return KineticFactor(temperature) * prescribed_speed / dx;
        return KineticFactor(temperature) * phase_field.StabilityRate(dx);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar Local(
        const State& state, const int i, const int j, const int k) const
    {
        const Set::Scalar weight = Weight(state, i, j, k);
        return weight * (input_is_rigid ?
            phase_field.Local(state.rigid_eta, 0, 1.0, i, j, k) :
            phase_field.Local(state.component_density, input,
                              1.0 / reference_density, i, j, k));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar LocalStabilityRate(
        const State& state, const int i, const int j, const int k) const
    {
        if (input_is_rigid &&
            !(Util::Clamp(state.rigid_eta(i,j,k), 0.0, 1.0) < 1.0))
            return 0.0;
        return Weight(state, i, j, k) *
            phase_field.LocalStabilityRate();
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar GradientCoefficient(
        const State& state, const int i, const int j, const int k) const
    {
        return Weight(state, i, j, k) *
            phase_field.GradientCoefficient();
    }

    int RigidComponent() const
    {
        // Prescribed motion is advanced explicitly through Apply/VolumeSource.
        return use_prescribed_speed ? -1 : rigid_component;
    }

    Set::Scalar HeatRelease() const { return heat_release; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar Transfer(
        const amrex::Array4<Set::Scalar>& component_density,
        const State& state, const Set::Scalar eta_change,
        const int i, const int j, const int k,
        Set::Scalar* heat = nullptr,
        const Model::Gas::Gas::DeviceData* gas = nullptr) const
    {
        const Set::Scalar mass_change = reference_density * eta_change;
        if (heat != nullptr)
        {
            Set::Scalar product_sensible_enthalpy = 0.0;
            if (gas != nullptr)
                for (int n = 0; n < nproducts; ++n)
                    product_sensible_enthalpy += mass_fractions[n] *
                        (Model::Gas::Gas::EnthalpyMassSpecies(*gas, state.temperature, products[n]) -
                         Model::Gas::Gas::EnthalpyMassSpecies(*gas, 300.0, products[n]));
            else if (state.temperature != 300.0)
                Util::Abort(INFO, "Phase-change enthalpy correction requires gas thermodynamic data");
            // Decomposition heat excludes formation enthalpy. Transfer the
            // solid's sensible enthalpy using the same 300 K reference.
            *heat -= mass_change * (heat_release +
                reactant_specific_heat * (state.temperature - 300.0) - product_sensible_enthalpy);
        }
        component_density(i,j,k,input) += mass_change;
        for (int n = 0; n < nproducts; ++n)
            component_density(i,j,k,products[n]) -=
                mass_fractions[n] * mass_change;

        const Set::Scalar product_specific_volume =
            product_gas_constant * state.temperature / state.pressure;
        return mass_change *
            (1.0 / reference_density - product_specific_volume);
    }

private:

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar KineticFactor(const Set::Scalar temperature) const
    {
        if (!(temperature > temperature_cutoff)) return 0.0;
        return rate_multiplier * (activation_temperature > 0.0 ?
            std::exp(-activation_temperature / temperature) : 1.0);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar Weight(const State& state,
                       const int i, const int j, const int k) const
    {
        const Set::Scalar factor = KineticFactor(state.temperature);
        if (!input_is_rigid || !(factor > 0.0)) return factor;

        const Set::Scalar eta =
            Util::Clamp(state.rigid_eta(i,j,k), 0.0, 1.0);
        if (!(eta > 0.0)) return 0.0;
        const Set::Scalar input_eta = Util::Clamp(
            state.rigid_species_eta(i,j,k,rigid_component), 0.0, eta);
        return factor * Model::PhaseField::H(input_eta / eta);
    }

    Model::PhaseField::PhaseField<Model::PhaseField::AllenCahn> phase_field;
    int input = -1;
    int rigid_component = -1;
    int nproducts = 0;
    std::array<int, max_products> products{};
    std::array<Set::Scalar, max_products> mass_fractions{};
    Set::Scalar reference_density = NAN;
    Set::Scalar product_gas_constant = 0.0;
    Set::Scalar rate_multiplier = 1.0;
    Set::Scalar activation_temperature = 0.0;
    Set::Scalar temperature_cutoff = 0.0;
    Set::Scalar prescribed_speed = 0.0;
    bool use_prescribed_speed = false;
    bool input_is_rigid = false;
    Set::Scalar heat_release = 0.0;
    Set::Scalar reactant_specific_heat = NAN;
};
}
}

#endif

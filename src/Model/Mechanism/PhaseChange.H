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
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Model/Mechanism/Mechanism.H"
#include "Model/PhaseField/AllenCahn.H"
#include "Model/PhaseField/PhaseField.H"
#include "Set/Set.H"

namespace Model
{
namespace Mechanism
{
class MillerMobilityTable
{
public:
    struct Entry
    {
        std::array<double, 3> normal;
        double multiplier = 1.0;
    };

    static MillerMobilityTable Read(
        const std::string& filename,
        const std::vector<double>& lattice_constants)
    {
        if (lattice_constants.size() != 3)
            Util::Exception(
                INFO,
                "orientation_mobility.miller_read requires three lattice constants");

        MillerMobilityTable table;
        std::ifstream file(filename);
        if (!file.is_open())
            Util::Exception(
                INFO,
                "Unable to open Miller mobility table ", filename);

        std::string line;
        int line_number = 0;
        while (std::getline(file, line))
        {
            ++line_number;
            const std::size_t comment = line.find('#');
            if (comment != std::string::npos) line.erase(comment);

            std::istringstream row(line);
            double h = 0.0, k = 0.0, l = 0.0, multiplier = 0.0;
            if (!(row >> h >> k >> l >> multiplier)) continue;

            std::array<double, 3> normal = {
                h / lattice_constants[0],
                k / lattice_constants[1],
                l / lattice_constants[2]};
            const double norm = std::sqrt(
                normal[0] * normal[0] +
                normal[1] * normal[1] +
                normal[2] * normal[2]);

            if (!(norm > 0.0) || !std::isfinite(norm))
                Util::Exception(
                    INFO,
                    "Invalid Miller indices in ", filename,
                    " at line ", line_number);
            if (!(multiplier > 0.0) || !std::isfinite(multiplier))
                Util::Exception(
                    INFO,
                    "Miller mobility multiplier must be finite and positive in ",
                    filename, " at line ", line_number);

            for (double& component : normal) component /= norm;
            table.entries.push_back({normal, multiplier});
        }

        if (table.entries.empty())
            Util::Exception(
                INFO,
                "Miller mobility table ", filename, " did not contain any rows");

        return table;
    }

    double Lookup(const std::vector<double>& direction) const
    {
        if (entries.empty())
            Util::Exception(INFO, "Cannot lookup from an empty Miller mobility table");
        if (direction.size() != 3)
            Util::Exception(INFO, "Miller mobility lookup requires a 3D direction");

        std::array<double, 3> normal = {
            direction[0], direction[1], direction[2]};
        const double norm = std::sqrt(
            normal[0] * normal[0] +
            normal[1] * normal[1] +
            normal[2] * normal[2]);
        if (!(norm > 0.0) || !std::isfinite(norm))
            Util::Exception(INFO, "Miller mobility lookup direction must be finite and nonzero");

        for (double& component : normal) component /= norm;

        double best_dot = -1.0;
        double best_multiplier = entries.front().multiplier;
        for (const Entry& entry : entries)
        {
            const double dot = std::abs(
                normal[0] * entry.normal[0] +
                normal[1] * entry.normal[1] +
                normal[2] * entry.normal[2]);
            if (dot > best_dot)
            {
                best_dot = dot;
                best_multiplier = entry.multiplier;
            }
        }

        return best_multiplier;
    }

    std::vector<Entry> entries;
};

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
                    const std::vector<std::vector<double>>& orientation,
                    const std::vector<int>& polycrystal_species,
                    const std::vector<std::vector<std::vector<double>>>&
                        all_grain_orientations,
                    const std::vector<std::vector<double>>& lattice_constants,
                    std::vector<std::vector<double>>&
                        all_grain_mobility_multipliers,
                    const std::vector<double>& molecular_weight,
                    const Set::Scalar universal_gas_constant,
                    Set::Scalar* specific_heat = nullptr,
                    Set::Scalar* thermal_conductivity = nullptr)
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
        pp.select<Model::PhaseField::AllenCahn>(
            "phase_field", value.phase_field);
        for (int polycrystal = 0;
             polycrystal < static_cast<int>(polycrystal_species.size());
             ++polycrystal)
            if (value.input == polycrystal_species[polycrystal])
                value.orientation_mobility_component = polycrystal;
        const std::vector<std::vector<double>> empty_grain_orientations;
        std::vector<double> unused_grain_mobility_multipliers;
        const auto& grain_orientations =
            value.orientation_mobility_component >= 0 ?
                all_grain_orientations[value.orientation_mobility_component] :
                empty_grain_orientations;
        auto& grain_mobility_multipliers =
            value.orientation_mobility_component >= 0 ?
                all_grain_mobility_multipliers[
                    value.orientation_mobility_component] :
                unused_grain_mobility_multipliers;
        const int selected_polycrystal_species =
            value.orientation_mobility_component >= 0 ? value.input : -1;
        const auto& normal = orientation[value.input];
        std::string orientation_mobility_type = "formula";
        pp.query_default(
            "orientation_mobility.type",
            orientation_mobility_type,
            std::string("formula"));

        if (orientation_mobility_type == "formula")
        {
            constexpr Set::Scalar orientation_strength = 0.20;
            constexpr Set::Scalar reference_ab_product = 0.2368;

            const auto formula_multiplier =
                [=](const std::vector<double>& direction)
                {
                    const Set::Scalar na2 = direction[0] * direction[0];
                    const Set::Scalar nb2 = direction[1] * direction[1];
                    const Set::Scalar nc2 = direction[2] * direction[2];
                    return Util::Clamp(
                        1.0 + orientation_strength *
                            (na2 * nb2 / reference_ab_product - nc2),
                        1.0 - orientation_strength,
                        1.0 + orientation_strength);
                };

            value.orientation_mobility_multiplier =
                formula_multiplier(normal);

            if (value.input == selected_polycrystal_species)
            {
                grain_mobility_multipliers.clear();
                grain_mobility_multipliers.reserve(
                    grain_orientations.size());
                for (const auto& grain_orientation : grain_orientations)
                    grain_mobility_multipliers.push_back(
                        formula_multiplier(grain_orientation));
            }
        }
        else if (orientation_mobility_type == "miller_read")
        {
            std::string filename;

            pp.query_file(
                "orientation_mobility.miller_read.filename",
                filename);

            if (!IO::ParmParse::InTraversalMode())
            {
                const MillerMobilityTable mobility_table =
                    MillerMobilityTable::Read(
                        filename, lattice_constants[value.input]);

                value.orientation_mobility_multiplier =
                    mobility_table.Lookup(normal);

                if (value.input == selected_polycrystal_species)
                {
                    grain_mobility_multipliers.clear();
                    grain_mobility_multipliers.reserve(
                        grain_orientations.size());
                    for (const auto& grain_orientation :
                         grain_orientations)
                        grain_mobility_multipliers.push_back(
                            mobility_table.Lookup(grain_orientation));
                }
            }
        }
        else if (orientation_mobility_type == "constant")
        {
            Set::Scalar mobility = 0.0;
            pp.query_required(
                "orientation_mobility.constant.value",
                mobility,
                1.0 / Unit::Pressure() / Unit::Time());

            const Set::Scalar base_mobility = value.phase_field.Mobility();
            if (!(base_mobility > 0.0) || !std::isfinite(base_mobility))
                Util::Exception(
                    INFO,
                    "phase_field.allencahn.mobility must be finite and positive ",
                    "when orientation_mobility.type = constant");
            value.orientation_mobility_multiplier = mobility / base_mobility;

            if (value.input == selected_polycrystal_species)
            {
                grain_mobility_multipliers.clear();
                grain_mobility_multipliers.resize(
                    grain_orientations.size(),
                    value.orientation_mobility_multiplier);
            }
        }
        else
        {
            Util::Exception(
                INFO,
                "Invalid orientation_mobility.type ",
                orientation_mobility_type,
                ". Expected formula, miller_read, or constant.");
        }
        if (!(value.orientation_mobility_multiplier > 0.0) ||
            !std::isfinite(value.orientation_mobility_multiplier))
            Util::Exception(INFO,
                "orientation mobility multiplier must be finite and positive");
        for (const double multiplier : grain_mobility_multipliers)
            if (!(multiplier > 0.0) || !std::isfinite(multiplier))
                Util::Exception(
                    INFO,
                    "Grain mobility multipliers must be finite and positive");
        value.use_local_orientation_mobility =
            value.input == selected_polycrystal_species;
        value.stability_orientation_mobility_multiplier =
            value.orientation_mobility_multiplier;
        if (value.use_local_orientation_mobility)
            for (const double multiplier : grain_mobility_multipliers)
                value.stability_orientation_mobility_multiplier = Util::Max(
                    value.stability_orientation_mobility_multiplier,
                    multiplier);
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
        // Energy per unit mass transferred out of the condensed phase.
        // Negative values represent endothermic decomposition. The homogeneous
        // mode below instead computes this from the two constituent heats.
        bool homogeneous;
        pp.query_default("homogeneous", homogeneous, false);
        if (!homogeneous)
        {
            pp.query_default("heat_release", value.heat_release,
                             "0.0_J/kg", Unit::Energy() / Unit::Mass());
            if (!std::isfinite(value.heat_release))
                Util::Exception(INFO, "phase-change heat_release must be finite");
        }
        if (homogeneous)
        {
            if (pp.contains("heat_release"))
                Util::Exception(INFO, "Use homogeneous.binder_heat_release and homogeneous.ap_heat_release for a blend");
            if (!value.input_is_rigid || specific_heat == nullptr || thermal_conductivity == nullptr)
                Util::Exception(INFO, "homogeneous phase change requires a rigid solid with thermal properties");
            std::string ap_solid, ap_gas;
            pp.query_required("homogeneous.ap_solid", ap_solid);
            pp.query_required("homogeneous.ap_gas", ap_gas);
            int ap = -1, gas = -1;
            for (int n = 0; n < static_cast<int>(species_names.size()); ++n)
            {
                if (n >= ngas_species && species_names[n] == ap_solid) ap = n;
                if (n < ngas_species && species_names[n] == ap_gas) gas = n;
            }
            if (ap < 0 || ap == value.input || gas < 0)
                Util::Exception(INFO, "homogeneous.ap_solid/ap_gas must identify AP species distinct from the binder");

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
            for (const int n : {binder, ap})
                if (!(reference_density[n] > 0.0) || !std::isfinite(reference_density[n]) ||
                    !(specific_heat[n] > 0.0) || !std::isfinite(specific_heat[n]) ||
                    !(thermal_conductivity[n] > 0.0) || !std::isfinite(thermal_conductivity[n]))
                    Util::Exception(INFO, "Homogeneous constituents require positive finite thermal properties");
            const Set::Scalar t = w * reference_density[binder] /
                (w * reference_density[binder] + (1.0 - w) * reference_density[ap]);

            // Add constituent volumes: rho is volume-weighted (equivalently,
            // 1/rho is mass-weighted). Specific heat remains mass-weighted.
            value.reference_density = reference_density[binder] =
                (1.0 - t) * reference_density[binder] + t * reference_density[ap];
            specific_heat[binder] = (1.0 - w) * specific_heat[binder] + w * specific_heat[ap];

            // Chen et al. (2002), eqs. (9)/(11): take the unsquared physical
            // branch f-x = (1-t)(1-x) f^(1/d), bounded by the pure values.
            const Set::Scalar kb = thermal_conductivity[binder], kap = thermal_conductivity[ap];
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
            pp.query_required("homogeneous.ap_rate_multiplier", ap_rate);
            pp.query_required("homogeneous.ap_activation_temperature", ap_activation, Unit::Temperature());
            pp.query_required("homogeneous.binder_heat_release", binder_heat, Unit::Energy() / Unit::Mass());
            pp.query_required("homogeneous.ap_heat_release", ap_heat, Unit::Energy() / Unit::Mass());
            if (!(ap_rate >= 0.0) || !std::isfinite(ap_rate) || !std::isfinite(value.rate_multiplier) ||
                !(ap_activation >= 0.0) || !std::isfinite(ap_activation) ||
                !std::isfinite(value.activation_temperature) ||
                !std::isfinite(binder_heat) || !std::isfinite(ap_heat))
                Util::Exception(INFO, "Invalid homogeneous kinetic or heat-release parameters");
            if (t == 1.0) value.rate_multiplier = ap_rate;
            else if (t > 0.0)
                value.rate_multiplier = value.rate_multiplier == 0.0 || ap_rate == 0.0 ? 0.0 :
                    std::exp((1.0 - t) * std::log(value.rate_multiplier) + t * std::log(ap_rate));
            value.activation_temperature = (1.0 - t) * value.activation_temperature + t * ap_activation;
            value.heat_release = (1.0 - w) * binder_heat + w * ap_heat;

            // The binder species now contains subgrid AP mass as well.
            for (int n = 0; n < value.nproducts; ++n) value.mass_fractions[n] *= 1.0 - w;
            int slot = 0;
            while (slot < value.nproducts && value.products[slot] != gas) ++slot;
            if (slot == max_products) Util::Exception(INFO, "Too many homogeneous gas products");
            if (slot == value.nproducts) value.products[value.nproducts++] = gas;
            value.mass_fractions[slot] += w;
            value.product_gas_constant = 0.0;
            for (int n = 0; n < value.nproducts; ++n)
                value.product_gas_constant += universal_gas_constant * value.mass_fractions[n] /
                    molecular_weight[value.products[n]];
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar operator()(const State& state,
                           const int i, const int j, const int k,
                           const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Set::Scalar weight = Weight(state, i, j, k);
        return OrientationMobility(state, i, j, k) * weight * (input_is_rigid ?
            phase_field(state.rigid_eta, 0, 1.0, i, j, k, DX) :
            phase_field(state.component_density, input,
                        1.0 / reference_density, i, j, k, DX));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar VolumeSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if (input_is_rigid) return 0.0;
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
        if (input_is_rigid) return;
        Transfer(component_density_rhs, state,
                 (*this)(state, i, j, k, DX), i, j, k);
    }

    Set::Scalar StabilityRate(const Set::Scalar dx,
                            const Set::Scalar temperature) const
    {
        return stability_orientation_mobility_multiplier *
            KineticFactor(temperature) * phase_field.StabilityRate(dx);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar Local(
        const State& state, const int i, const int j, const int k) const
    {
        const Set::Scalar weight = Weight(state, i, j, k);
        return OrientationMobility(state, i, j, k) * weight * (input_is_rigid ?
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
        return OrientationMobility(state, i, j, k) *
            Weight(state, i, j, k) *
            phase_field.LocalStabilityRate();
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar GradientCoefficient(
        const State& state, const int i, const int j, const int k) const
    {
        return OrientationMobility(state, i, j, k) *
            Weight(state, i, j, k) *
            phase_field.GradientCoefficient();
    }

    int RigidComponent() const
    {
        return rigid_component;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar Transfer(
        const amrex::Array4<Set::Scalar>& component_density,
        const State& state, const Set::Scalar eta_change,
        const int i, const int j, const int k,
        Set::Scalar* heat = nullptr) const
    {
        const Set::Scalar mass_change = reference_density * eta_change;
        if (heat != nullptr) *heat -= mass_change * heat_release;
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
    Set::Scalar OrientationMobility(
        const State& state, const int i, const int j, const int k) const
    {
        return use_local_orientation_mobility ?
            state.orientation_mobility(
                i, j, k, orientation_mobility_component) :
            orientation_mobility_multiplier;
    }

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
    Set::Scalar orientation_mobility_multiplier = 1.0;
    Set::Scalar stability_orientation_mobility_multiplier = 1.0;
    bool input_is_rigid = false;
    Set::Scalar heat_release = 0.0;
    bool use_local_orientation_mobility = false;
    int orientation_mobility_component = -1;
};
}
}

#endif

#ifndef UTIL_CHEMISTRY_H
#define UTIL_CHEMISTRY_H

#include "AMReX_Extension.H"
#include "IO/CanteraYamlParse.H"
#include <map>
#include <string>
#include <vector>
#include <array>
#include <set>
#include <cstdint>
#include <cmath>

namespace Util
{
namespace Chemistry
{

inline bool IsThirdBodyToken(const std::string& token)
{
    return token == "M" || (!token.empty() && token[0] == '+');
}

inline bool IsThirdBodySpecies(const CanteraYAML::Reaction& reaction, const std::string& token)
{
    if (!reaction.third_body) return false;
    if (IsThirdBodyToken(token)) return true;

    std::map<std::string, double>::const_iterator reactant = reaction.reactants.find(token);
    std::map<std::string, double>::const_iterator product = reaction.products.find(token);
    if (reactant == reaction.reactants.end()) return false;
    if (product == reaction.products.end()) return false;
    return reactant->second == product->second;
}

template <int NS>
struct PreparedKinetics
{
    int NR = 0;
    int NE = 0;
    std::array<int, NS> yaml_species_idx;
    std::map<std::string, int> active_species_idx;
    std::vector<std::array<double, NS> > reactant_nu;
    std::vector<std::array<double, NS> > product_nu;
    std::vector<std::array<double, NS> > nu_net;
    std::vector<std::array<double, NS> > reaction_order;
    std::vector<std::array<double, NS> > efficiency;
    std::vector<std::string> element_names;
    std::vector<std::array<double, NS> > atom_matrix;

    PreparedKinetics()
    {
        for (int n = 0; n < NS; ++n) yaml_species_idx[n] = -1;
    }
};

template <int NS>
inline PreparedKinetics<NS>
PrepareKinetics(const CanteraYAML::CanteraKinetics& kinetics,
                const std::string& context, const int nspecies);

template <int NS>
inline PreparedKinetics<NS>
PrepareKinetics(const CanteraYAML::CanteraKinetics& kinetics, const std::string& context)
{
    Util::Assert(INFO, TEST(!kinetics.phase.empty()),
        context, " YAML must define at least one phase.");
    return PrepareKinetics<NS>(
        kinetics, context, static_cast<int>(kinetics.phase[0].species.size()));
}

template <int NS>
inline PreparedKinetics<NS>
PrepareKinetics(const CanteraYAML::CanteraKinetics& kinetics,
                const std::string& context, const int nspecies)
{
    PreparedKinetics<NS> prepared;

    Util::Assert(INFO, TEST(nspecies > 0 && nspecies <= NS),
        context, " gas species count must be between 1 and ", NS, ".");
    Util::Assert(INFO, TEST(!kinetics.phase.empty()),
        context, " YAML must define at least one phase.");
    Util::Assert(INFO, TEST((int)kinetics.phase[0].species.size() == nspecies),
        context, " active phase species count (", kinetics.phase[0].species.size(),
        ") does not match the gas species count (", nspecies, ").");

    for (int n = 0; n < nspecies; ++n)
    {
        const std::string& simulation_species_name = kinetics.phase[0].species[n];
        Util::Assert(INFO, TEST(prepared.active_species_idx.find(simulation_species_name) == prepared.active_species_idx.end()),
            context, " active phase species '", simulation_species_name, "' appears multiple times.");

        int match = -1;
        for (int j = 0; j < (int)kinetics.species.size(); ++j)
        {
            if (kinetics.species[j].name == simulation_species_name)
            {
                Util::Assert(INFO, TEST(match < 0),
                    "Species '", simulation_species_name, "' appears multiple times in the YAML species list.");
                match = j;
            }
        }

        Util::Assert(INFO, TEST(match >= 0),
            "Species '", simulation_species_name, "' from the active phase is missing from the YAML species table.");

        prepared.active_species_idx[simulation_species_name] = n;
        prepared.yaml_species_idx[n] = match;
    }

    prepared.NR = kinetics.reaction.size();
    prepared.reactant_nu.assign(prepared.NR, std::array<double, NS>());
    prepared.product_nu.assign(prepared.NR, std::array<double, NS>());
    prepared.nu_net.assign(prepared.NR, std::array<double, NS>());
    prepared.reaction_order.assign(prepared.NR, std::array<double, NS>());
    prepared.efficiency.assign(prepared.NR, std::array<double, NS>());

    for (int r = 0; r < prepared.NR; ++r)
    {
        for (int n = 0; n < NS; ++n)
        {
            prepared.reaction_order[r][n] = -1.0;
            prepared.efficiency[r][n] = -1.0;
        }

        const CanteraYAML::Reaction& reaction = kinetics.reaction[r];

        for (const std::pair<const std::string, double>& reactant : reaction.reactants)
        {
            if (IsThirdBodySpecies(reaction, reactant.first)) continue;
            std::map<std::string, int>::const_iterator species_idx = prepared.active_species_idx.find(reactant.first);
            Util::Assert(INFO, TEST(species_idx != prepared.active_species_idx.end()),
                "Reaction '", reaction.equation, "' references reactant species '",
                reactant.first, "' that is not present in the active phase species list.");
            prepared.reactant_nu[r][species_idx->second] = reactant.second;
        }

        for (const std::pair<const std::string, double>& product : reaction.products)
        {
            if (IsThirdBodySpecies(reaction, product.first)) continue;
            std::map<std::string, int>::const_iterator species_idx = prepared.active_species_idx.find(product.first);
            Util::Assert(INFO, TEST(species_idx != prepared.active_species_idx.end()),
                "Reaction '", reaction.equation, "' references product species '",
                product.first, "' that is not present in the active phase species list.");
            prepared.product_nu[r][species_idx->second] = product.second;
        }

        for (const std::pair<const std::string, double>& order : reaction.orders)
        {
            std::map<std::string, int>::const_iterator species_idx = prepared.active_species_idx.find(order.first);
            Util::Assert(INFO, TEST(species_idx != prepared.active_species_idx.end()),
                "Reaction '", reaction.equation, "' defines reaction order for species '",
                order.first, "' that is not present in the active phase species list.");
            prepared.reaction_order[r][species_idx->second] = order.second;
        }

        for (const std::pair<const std::string, double>& collider : reaction.efficiencies)
        {
            std::map<std::string, int>::const_iterator species_idx = prepared.active_species_idx.find(collider.first);
            Util::Assert(INFO, TEST(species_idx != prepared.active_species_idx.end()),
                "Reaction '", reaction.equation, "' defines third-body efficiency for species '",
                collider.first, "' that is not present in the active phase species list.");
            prepared.efficiency[r][species_idx->second] = collider.second;
        }

        std::map<std::string, double> element_balance;
        for (int n = 0; n < nspecies; ++n)
        {
            prepared.nu_net[r][n] = prepared.product_nu[r][n] - prepared.reactant_nu[r][n];
            const CanteraYAML::Species& species = kinetics.species[prepared.yaml_species_idx[n]];
            for (const std::pair<const std::string, double>& composition : species.composition)
            {
                element_balance[composition.first] += prepared.nu_net[r][n] * composition.second;
            }
        }

        for (const std::pair<const std::string, double>& element : element_balance)
        {
            Util::Assert(INFO, TEST(std::abs(element.second) < 1.0e-12),
                "Reaction '", reaction.equation, "' does not conserve element '",
                element.first, "' after species mapping. Residual = ", element.second);
        }
    }

    std::set<std::string> element_set;
    for (int n = 0; n < nspecies; ++n)
    {
        const CanteraYAML::Species& species = kinetics.species[prepared.yaml_species_idx[n]];
        for (const std::pair<const std::string, double>& composition : species.composition)
        {
            element_set.insert(composition.first);
        }
    }

    prepared.element_names.assign(element_set.begin(), element_set.end());
    prepared.NE = prepared.element_names.size();
    prepared.atom_matrix.assign(prepared.NE, std::array<double, NS>());
    for (int e = 0; e < prepared.NE; ++e)
    {
        for (int n = 0; n < NS; ++n) prepared.atom_matrix[e][n] = 0.0;
    }

    for (int n = 0; n < nspecies; ++n)
    {
        const CanteraYAML::Species& species = kinetics.species[prepared.yaml_species_idx[n]];
        for (const std::pair<const std::string, double>& composition : species.composition)
        {
            std::set<std::string>::const_iterator it = element_set.find(composition.first);
            Util::Assert(INFO, TEST(it != element_set.end()),
                "Element '", composition.first, "' for species '", species.name,
                "' was not found while building the atom matrix.");
            int e = std::distance(element_set.begin(), it);
            prepared.atom_matrix[e][n] = composition.second;
        }
    }

    return prepared;
}

} // namespace Chemistry
} // namespace Util

#endif

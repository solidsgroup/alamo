#ifndef MODEL_CHEMISTRY_FINITERATE_H
#define MODEL_CHEMISTRY_FINITERATE_H

#include "AMReX_Extension.H"
#include "Model/Chemistry/Chemistry.H"
#include "IO/YamlParse.H"
#include <vector>
#include <map>
#include <cmath>

namespace Model
{
namespace Chemistry
{

struct FiniteRate
{
    static constexpr const char* name = "finite_rate";

    static bool IsThirdBodyToken(const std::string& token)
    {
        return token == "M" || (!token.empty() && token[0] == '+');
    }

    static bool IsThirdBodySpecies(const CanteraYAML::Reaction& reaction, const std::string& token)
    {
        if (!reaction.third_body) return false;
        if (IsThirdBodyToken(token)) return true;

        std::map<std::string, double>::const_iterator reactant = reaction.reactants.find(token);
        std::map<std::string, double>::const_iterator product = reaction.products.find(token);
        if (reactant == reaction.reactants.end()) return false;
        if (product == reaction.products.end()) return false;
        return reactant->second == product->second;
    }

    FiniteRate() = default;
    FiniteRate(IO::ParmParse &pp, std::string name)
    {
        pp_queryclass(name,*this);
    }

#ifdef ALAMO_YAML
    std::string yaml;
    int NR = 0;
    std::vector<std::array<double, NSPECIES> > reactant_nu;
    std::vector<std::array<double, NSPECIES> > product_nu;
    std::vector<std::array<double, NSPECIES> > reaction_order;
    std::vector<std::array<double, NSPECIES> > efficiency;
    std::vector<std::array<double, NSPECIES> > nu_net;
    std::vector<double> A;
    std::vector<double> b;
    std::vector<double> E;
    std::vector<double> A2;
    std::vector<double> b2;
    std::vector<double> E2;
    std::vector<double> F;
    std::vector<int> reversible;
    std::vector<int> third_body;
    std::vector<int> falloff;
#endif

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    std::pair<std::array<double,NSPECIES>,Set::Scalar> evaluate_source(
        const Set::Scalar P,
        const Set::Scalar T,
        const std::array<double, NSPECIES>& rhoY,
        const Model::Gas::Gas* gas) const
    {
        std::array<double, NSPECIES> wdot;
        std::array<double, NSPECIES> wdot_molar;
        for (int n = 0; n < NSPECIES; ++n)
        {
            wdot[n] = 0.0;
            wdot_molar[n] = 0.0;
        }

#ifdef ALAMO_YAML
        if constexpr (NSPECIES > 1)
        {
            Util::Assert(INFO, TEST(gas != nullptr), "Finite-rate kinetics requires a gas model.");

            const double R = Set::Constant::Rg;
            const double P0 = 101325.0;

            double rho = 0.0;
            double total_concentration = 0.0;
            std::array<double, NSPECIES> concentration;
            std::array<double, NSPECIES> mole_fraction;
            std::array<double, NSPECIES> enthalpy_mol;
            std::array<double, NSPECIES> entropy_mol;

            for (int n = 0; n < NSPECIES; ++n)
            {
                rho += rhoY[n];
                concentration[n] = std::max(rhoY[n], 0.0) / gas->MW[n];
                total_concentration += concentration[n];
                enthalpy_mol[n] = gas->enthalpy_mol_species(T, n);
                entropy_mol[n] = gas->entropy_mol_species(T, n);
            }

            if (!(rho > 0.0) || !(total_concentration > 0.0))
            {
                return {wdot, 0.0};
            }

            for (int n = 0; n < NSPECIES; ++n)
            {
                mole_fraction[n] = concentration[n] / total_concentration;
            }

            for (int r = 0; r < NR; ++r)
            {
                double kf = 0.0;
                if (falloff[r])
                {
                    double k_inf = A[r] * std::pow(T, b[r]) * std::exp(-E[r] / (R * T));
                    double k0 = A2[r] * std::pow(T, b2[r]) * std::exp(-E2[r] / (R * T));
                    double M_concentration = 0.0;
                    for (int n = 0; n < NSPECIES; ++n)
                    {
                        double collider_efficiency = efficiency[r][n];
                        if (collider_efficiency < 0.0) collider_efficiency = 1.0;
                        M_concentration += collider_efficiency * concentration[n];
                    }
                    double Pr = 0.0;
                    if (k_inf > 0.0) Pr = k0 * M_concentration / k_inf;
                    kf = k_inf * Pr / (1.0 + Pr) * F[r];
                }
                else
                {
                    kf = A[r] * std::pow(T, b[r]) * std::exp(-E[r] / (R * T));
                }

                double dnu = 0.0;
                double dH0 = 0.0;
                double dS0 = 0.0;
                for (int n = 0; n < NSPECIES; ++n)
                {
                    dnu += nu_net[r][n];
                    dH0 += nu_net[r][n] * enthalpy_mol[n];
                    dS0 += nu_net[r][n] * (entropy_mol[n] - R * std::log(std::max(P / P0, 1.0e-300)));
                }

                double Kc = 1.0;
                if (reversible[r])
                {
                    double dG0 = dH0 - T * dS0;
                    double Kp = std::exp(-dG0 / (R * T));
                    Kc = Kp * std::pow(P0 / (R * T), dnu);
                }

                double C_react = 1.0;
                double C_prod = 1.0;
                double third_body_factor = 0.0;
                for (int n = 0; n < NSPECIES; ++n)
                {
                    if (!third_body[r])
                    {
                        if (reactant_nu[r][n] > 0.0)
                        {
                            double exponent = reaction_order[r][n] >= 0.0 ? reaction_order[r][n] : reactant_nu[r][n];
                            C_react *= std::pow(concentration[n], exponent);
                        }
                        if (reversible[r] && product_nu[r][n] > 0.0)
                        {
                            C_prod *= std::pow(concentration[n], product_nu[r][n]);
                        }
                        third_body_factor = 1.0;
                    }
                    else
                    {
                        if (reactant_nu[r][n] > 0.0) C_react *= std::pow(concentration[n], reactant_nu[r][n]);
                        if (reversible[r] && product_nu[r][n] > 0.0) C_prod *= std::pow(concentration[n], product_nu[r][n]);

                        double collider_efficiency = efficiency[r][n];
                        if (collider_efficiency < 0.0) collider_efficiency = 1.0;
                        third_body_factor += collider_efficiency * concentration[n];
                    }
                }

                if (falloff[r]) third_body_factor = 1.0;

                double ropf = kf * third_body_factor * C_react;
                double ropr = 0.0;
                if (reversible[r] && std::abs(Kc) > 1.0e-300)
                {
                    ropr = kf / Kc * third_body_factor * C_prod;
                }
                double rop = ropf - ropr;

                for (int n = 0; n < NSPECIES; ++n)
                {
                    wdot_molar[n] += nu_net[r][n] * rop;
                }
            }

            double rho_total = 0.0;
            for (int n = 0; n < NSPECIES; ++n)
            {
                rho_total += std::max(rhoY[n], 0.0);
                wdot[n] = wdot_molar[n] * gas->MW[n];
            }

            double mass_source_total = 0.0;
            for (int n = 0; n < NSPECIES; ++n)
            {
                mass_source_total += wdot[n];
            }

            if (std::abs(mass_source_total) > 0.0 && rho_total > 0.0)
            {
                for (int n = 0; n < NSPECIES; ++n)
                {
                    double Y_n = std::max(rhoY[n], 0.0) / rho_total;
                    wdot[n] -= Y_n * mass_source_total;
                }
            }

            double qdot = 0.0;
            for (int n = 0; n < NSPECIES; ++n)
            {
                double Y_n = 0.0;
                if (rho_total > 0.0) Y_n = std::max(rhoY[n], 0.0) / rho_total;
                qdot -= Y_n * gas->enthalpy_mol_species(T, n) / gas->MW[n] * wdot[n];
            }
            return {wdot, qdot};
        }
#else
        Util::IgnoreUnused(P, T, rhoY, gas);
#endif

        return {wdot, 0.0};
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    void limit_wdot_positivity(
        std::array<double, NSPECIES>& wdot,
        const std::array<double, NSPECIES>& rhoY,
        const Set::Scalar dt) const
    {
        if (!(dt > 0.0)) return;

        double limiter = 1.0;
        for (int n = 0; n < NSPECIES; ++n)
        {
            if (wdot[n] < 0.0)
            {
                double available = std::max(rhoY[n], 0.0);
                double requested = -wdot[n] * dt;
                if (requested > available && requested > 0.0)
                {
                    limiter = std::min(limiter, available / requested);
                }
            }
        }

        if (limiter < 1.0)
        {
            for (int n = 0; n < NSPECIES; ++n) wdot[n] *= limiter;
        }
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    std::pair<std::array<double,NSPECIES>,Set::Scalar>
    compute(
        const Set::Scalar P,
        const Set::Scalar T,
        const std::array<double, NSPECIES>& rhoY,
        const Set::Scalar dt,
        const Model::Gas::Gas* gas) const
    {
        if (dt <= 0.0)
        {
            std::pair<std::array<double, NSPECIES>, Set::Scalar> source = evaluate_source(P, T, rhoY, gas);
            limit_wdot_positivity(source.first, rhoY, dt);
            source.second = 0.0;
            double rho_total = 0.0;
            for (int n = 0; n < NSPECIES; ++n) rho_total += std::max(rhoY[n], 0.0);
            for (int n = 0; n < NSPECIES; ++n)
            {
                double Y_n = 0.0;
                if (rho_total > 0.0) Y_n = std::max(rhoY[n], 0.0) / rho_total;
                source.second -= Y_n * gas->enthalpy_mol_species(T, n) / gas->MW[n] * source.first[n];
            }
            return source;
        }

        return evaluate_source(P, T, rhoY, gas);
    }

    static void Parse(FiniteRate& value, IO::ParmParse& pp)
    {
        Util::Message(INFO,pp.getPrefix());
#ifdef ALAMO_YAML
        if constexpr (NSPECIES > 1)
        {
            pp_query_required("yaml", value.yaml);
            CanteraYAML::CanteraKinetics kinetics = CanteraYAML::ParseYaml(value.yaml);

            Util::Assert(INFO, TEST(!kinetics.phase.empty()),
                "Finite-rate chemistry YAML must define at least one phase.");
            Util::Assert(INFO, TEST((int)kinetics.phase[0].species.size() == NSPECIES),
                "Finite-rate chemistry active phase species count (", kinetics.phase[0].species.size(),
                ") does not match NSPECIES (", NSPECIES, ").");

            value.NR = kinetics.reaction.size();
            value.reactant_nu.assign(value.NR, std::array<double, NSPECIES>());
            value.product_nu.assign(value.NR, std::array<double, NSPECIES>());
            value.reaction_order.assign(value.NR, std::array<double, NSPECIES>());
            value.efficiency.assign(value.NR, std::array<double, NSPECIES>());
            value.nu_net.assign(value.NR, std::array<double, NSPECIES>());
            value.A.assign(value.NR, 0.0);
            value.b.assign(value.NR, 0.0);
            value.E.assign(value.NR, 0.0);
            value.A2.assign(value.NR, 0.0);
            value.b2.assign(value.NR, 0.0);
            value.E2.assign(value.NR, 0.0);
            value.F.assign(value.NR, 1.0);
            value.reversible.assign(value.NR, 0);
            value.third_body.assign(value.NR, 0);
            value.falloff.assign(value.NR, 0);

            std::array<int, NSPECIES> yaml_species_idx;
            for (int n = 0; n < NSPECIES; ++n) yaml_species_idx[n] = -1;
            for (int n = 0; n < NSPECIES; ++n)
            {
                const std::string& simulation_species_name = kinetics.phase[0].species[n];
                for (int j = 0; j < (int)kinetics.species.size(); ++j)
                {
                    if (kinetics.species[j].name == simulation_species_name)
                    {
                        yaml_species_idx[n] = j;
                        break;
                    }
                }
                Util::Assert(INFO, TEST(yaml_species_idx[n] >= 0),
                    "Species '", simulation_species_name, "' from the active phase is missing from the YAML species table.");
            }

            for (int r = 0; r < value.NR; ++r)
            {
                for (int n = 0; n < NSPECIES; ++n)
                {
                    value.reactant_nu[r][n] = 0.0;
                    value.product_nu[r][n] = 0.0;
                    value.reaction_order[r][n] = -1.0;
                    value.efficiency[r][n] = -1.0;
                    value.nu_net[r][n] = 0.0;
                }
            }

            for (int r = 0; r < value.NR; ++r)
            {
                const CanteraYAML::Reaction& reaction = kinetics.reaction[r];
                value.A[r] = reaction.A;
                value.b[r] = reaction.b;
                value.E[r] = reaction.E;
                value.A2[r] = reaction.A2;
                value.b2[r] = reaction.b2;
                value.E2[r] = reaction.E2;
                value.F[r] = reaction.F;
                value.reversible[r] = reaction.reversible ? 1 : 0;
                value.third_body[r] = reaction.third_body ? 1 : 0;
                value.falloff[r] = reaction.falloff ? 1 : 0;

                for (int n = 0; n < NSPECIES; ++n)
                {
                    const std::string& species_name = kinetics.phase[0].species[n];

                    std::map<std::string, double>::const_iterator reactant = reaction.reactants.find(species_name);
                    if (reactant != reaction.reactants.end() && !IsThirdBodySpecies(reaction, species_name))
                    {
                        value.reactant_nu[r][n] = reactant->second;
                    }

                    std::map<std::string, double>::const_iterator product = reaction.products.find(species_name);
                    if (product != reaction.products.end() && !IsThirdBodySpecies(reaction, species_name))
                    {
                        value.product_nu[r][n] = product->second;
                    }

                    std::map<std::string, double>::const_iterator order = reaction.orders.find(species_name);
                    if (order != reaction.orders.end())
                    {
                        value.reaction_order[r][n] = order->second;
                    }

                    std::map<std::string, double>::const_iterator collider = reaction.efficiencies.find(species_name);
                    if (collider != reaction.efficiencies.end())
                    {
                        value.efficiency[r][n] = collider->second;
                    }

                    value.nu_net[r][n] = value.product_nu[r][n] - value.reactant_nu[r][n];
                }

                std::map<std::string, double> element_balance;
                for (int n = 0; n < NSPECIES; ++n)
                {
                    const CanteraYAML::Species& species = kinetics.species[yaml_species_idx[n]];
                    for (const std::pair<const std::string, double>& composition : species.composition)
                    {
                        element_balance[composition.first] += value.nu_net[r][n] * composition.second;
                    }
                }
            for (const std::pair<const std::string, double>& element : element_balance)
            {
                Util::Assert(INFO, TEST(std::abs(element.second) < 1.0e-12),
                    "Finite-rate reaction '", reaction.equation, "' does not conserve element '",
                    element.first, "' after species mapping. Residual = ", element.second);
            }
        }

            Util::Message(INFO, "Finite-rate kinetics initialized from ", value.yaml, " with ", value.NR,
                " reactions using explicit integration.");

        }
#else
        (void)value;
        Util::Abort(INFO, "Alamo was not built with YAML support, but finite-rate chemistry was requested.");
#endif
    }

};

} // namespace Chemistry
} // namespace Model

#endif

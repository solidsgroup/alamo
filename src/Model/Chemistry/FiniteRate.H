#ifndef MODEL_CHEMISTRY_FINITERATE_H
#define MODEL_CHEMISTRY_FINITERATE_H

#include "AMReX_Extension.H"
#include "Model/Chemistry/Chemistry.H"
#include "Util/Chemistry.H"
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

namespace Model
{
namespace Chemistry
{

struct FiniteRate
{
    static constexpr const char* name = "finite_rate";

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static double EvaluateTroeFactor(
        const double T,
        const double Pr,
        const double alpha,
        const double T3,
        const double T1,
        const double T2)
    {
        if (!(Pr > 0.0)) return 1.0;

        double Fcent = 0.0;
        if (T3 > 0.0) Fcent += (1.0 - alpha) * std::exp(-T / T3);
        if (T1 > 0.0) Fcent += alpha * std::exp(-T / T1);
        if (T2 > 0.0) Fcent += std::exp(-T2 / T);

        Fcent = std::max(1.0e-300, std::min(Fcent, 1.0));
        double logFcent = std::log10(Fcent);
        double c = -0.4 - 0.67 * logFcent;
        double n = 0.75 - 1.27 * logFcent;
        double d = 0.14;
        double x = std::log10(std::max(Pr, 1.0e-300)) + c;
        double f1 = x / (n - d * x);
        return std::pow(10.0, logFcent / (1.0 + f1 * f1));
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
    std::vector<double> troe_alpha;
    std::vector<double> troe_T1;
    std::vector<double> troe_T2;
    std::vector<double> troe_T3;
    std::vector<int> reversible;
    std::vector<int> third_body;
    std::vector<int> falloff;
    std::vector<int> troe;
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
            Util::IgnoreUnused(P);

            const double R = Set::Constant::Rg;
            const double P0 = 101325.0;

            double rho = 0.0;
            double total_concentration = 0.0;
            std::array<double, NSPECIES> concentration;
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
                    double falloff_factor = F[r];
                    if (troe[r])
                    {
                        falloff_factor *= EvaluateTroeFactor(
                            T, Pr, troe_alpha[r], troe_T3[r], troe_T1[r], troe_T2[r]);
                    }
                    kf = k_inf * Pr / (1.0 + Pr) * falloff_factor;
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
                    dS0 += nu_net[r][n] * entropy_mol[n];
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
                for (int n = 0; n < NSPECIES; ++n)
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
                }

                double third_body_factor = 1.0;
                if (third_body[r] && !falloff[r])
                {
                    third_body_factor = 0.0;
                    for (int n = 0; n < NSPECIES; ++n)
                    {
                        double collider_efficiency = efficiency[r][n];
                        if (collider_efficiency < 0.0) collider_efficiency = 1.0;
                        third_body_factor += collider_efficiency * concentration[n];
                    }
                }

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

            for (int n = 0; n < NSPECIES; ++n)
            {
                wdot[n] = wdot_molar[n] * gas->MW[n];
            }

            double mass_source_total = 0.0;
            double abs_mass_source_total = 0.0;
            for (int n = 0; n < NSPECIES; ++n)
            {
                mass_source_total += wdot[n];
                abs_mass_source_total += std::abs(wdot[n]);
            }

            double mass_source_tol = 1.0e-12 * std::max(1.0, abs_mass_source_total);
            Util::Assert(INFO, TEST(std::abs(mass_source_total) <= mass_source_tol),
                "Finite-rate kinetics produced a non-conservative mass source. Residual = ",
                mass_source_total, ", tolerance = ", mass_source_tol,
                ". Check species molecular weights and reaction stoichiometry.");

            double qdot = 0.0;
            for (int n = 0; n < NSPECIES; ++n)
            {
                double e_n = (gas->enthalpy_mol_species(T, n) - R * T) / gas->MW[n];
                qdot -= wdot[n] * e_n;
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

    std::pair<std::array<double,NSPECIES>,Set::Scalar>
    compute(
        const Set::Scalar P,
        const Set::Scalar T,
        const std::array<double, NSPECIES>& rhoY,
        const Set::Scalar dt,
        const Model::Gas::Gas* gas) const
    {
        std::pair<std::array<double, NSPECIES>, Set::Scalar> source = evaluate_source(P, T, rhoY, gas);
        if (dt > 0.0)
        {
            limit_wdot_positivity(source.first, rhoY, dt);
            // Re-calculate qdot after limiting wdot to maintain energy consistency
            const double R = Set::Constant::Rg;
            source.second = 0.0;
            for (int n = 0; n < NSPECIES; ++n)
            {
                double e_n = (gas->enthalpy_mol_species(T, n) - R * T) / gas->MW[n];
                source.second -= source.first[n] * e_n;
            }
        }
        return source;
    }

    static void Parse(FiniteRate& value, IO::ParmParse& pp)
    {
        Util::Message(INFO,pp.getPrefix());
#ifdef ALAMO_YAML
        if constexpr (NSPECIES > 1)
        {
            pp_query_required("yaml", value.yaml);
            CanteraYAML::CanteraKinetics kinetics = CanteraYAML::ParseYaml(value.yaml);
            ChemUtil::PreparedKinetics<NSPECIES> prepared =
                ChemUtil::PrepareKinetics<NSPECIES>(kinetics, "Finite-rate chemistry");

            value.NR = prepared.NR;
            value.reactant_nu = prepared.reactant_nu;
            value.product_nu = prepared.product_nu;
            value.reaction_order = prepared.reaction_order;
            value.efficiency = prepared.efficiency;
            value.nu_net = prepared.nu_net;
            value.A.assign(value.NR, 0.0);
            value.b.assign(value.NR, 0.0);
            value.E.assign(value.NR, 0.0);
            value.A2.assign(value.NR, 0.0);
            value.b2.assign(value.NR, 0.0);
            value.E2.assign(value.NR, 0.0);
            value.F.assign(value.NR, 1.0);
            value.troe_alpha.assign(value.NR, 0.0);
            value.troe_T1.assign(value.NR, -1.0);
            value.troe_T2.assign(value.NR, -1.0);
            value.troe_T3.assign(value.NR, -1.0);
            value.reversible.assign(value.NR, 0);
            value.third_body.assign(value.NR, 0);
            value.falloff.assign(value.NR, 0);
            value.troe.assign(value.NR, 0);

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
                value.troe[r] = reaction.Troe.empty() ? 0 : 1;

                Util::Assert(INFO, TEST(!value.troe[r] || value.falloff[r]),
                    "Finite-rate reaction '", reaction.equation,
                    "' defines Troe parameters but is not a falloff reaction.");

                if (value.troe[r])
                {
                    std::map<std::string, double>::const_iterator alpha = reaction.Troe.find("A");
                    if (alpha == reaction.Troe.end()) alpha = reaction.Troe.find("alpha");
                    std::map<std::string, double>::const_iterator T1 = reaction.Troe.find("T1");
                    std::map<std::string, double>::const_iterator T2 = reaction.Troe.find("T2");
                    std::map<std::string, double>::const_iterator T3 = reaction.Troe.find("T3");

                    Util::Assert(INFO, TEST(alpha != reaction.Troe.end()),
                        "Finite-rate Troe reaction '", reaction.equation, "' is missing alpha/A.");
                    Util::Assert(INFO, TEST(T1 != reaction.Troe.end()),
                        "Finite-rate Troe reaction '", reaction.equation, "' is missing T1.");
                    Util::Assert(INFO, TEST(T3 != reaction.Troe.end()),
                        "Finite-rate Troe reaction '", reaction.equation, "' is missing T3.");

                    value.troe_alpha[r] = alpha->second;
                    value.troe_T1[r] = T1->second;
                    value.troe_T3[r] = T3->second;
                    if (T2 != reaction.Troe.end()) value.troe_T2[r] = T2->second;
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

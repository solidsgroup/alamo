#ifndef MODEL_CHEMISTRY_EQUILIBRIUM_H
#define MODEL_CHEMISTRY_EQUILIBRIUM_H

#include "AMReX_Extension.H"
#include "Model/Chemistry/Chemistry.H"
#include "IO/YamlParse.H"
#include "Model/Gas/Gas.H"
#include <eigen3/Eigen/Dense>
#include <set>
#include <algorithm>

namespace Model
{
namespace Chemistry
{

struct Equilibrium
{
    static constexpr const char* name = "equilibrium";

    Equilibrium() = default; 
    Equilibrium(IO::ParmParse &pp, std::string name) 
    { 
        pp_queryclass(name,*this); 
    } 

#ifdef ALAMO_YAML
    std::string yaml;
    CanteraYAML::CanteraKinetics kinetics;
    Eigen::MatrixXd nu; // Stoichiometric matrix (NSPECIES x NREACTIONS)
    Eigen::MatrixXd atom_matrix; // (NELEMENTS x NSPECIES)
    int NE = 0; // Number of elements
#endif

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    std::pair<std::array<double,NSPECIES>,Set::Scalar>
    compute(
        const Set::Scalar /*P_orig*/,
        const Set::Scalar T_orig,
        const std::array<double, NSPECIES>& rhoY,
        const Set::Scalar dt,
        const Model::Gas::Gas* gas) const
    {
        if (dt <= 0) return {};
        const double R = Set::Constant::Rg;
        const double P0 = 101325.0; // Standard pressure

        // 1. Calculate initial element counts and energy
        // n_i is in kmol/m^3 (assuming V=1)
        Eigen::VectorXd n_orig(NSPECIES);
        Eigen::VectorXd b(NE);
        b.setZero();
        double e_total_0 = 0.0;
        
        for (int i = 0; i < NSPECIES; ++i) {
            n_orig(i) = std::max(rhoY[i], 0.0) / gas->MW[i]; // kmol/m^3
            for (int k = 0; k < NE; ++k) b(k) += atom_matrix(k, i) * n_orig(i);
            e_total_0 += n_orig(i) * (gas->enthalpy_mol_species(T_orig, i) - R * T_orig);
        }

        // 2. Newton-Raphson for element potentials and T (UV solve)
        // Variables: gamma_1, ..., gamma_NE, ln(T)
        // gamma_k = lambda_k / RT
        // n_i = (P0 / RT) * exp( sum_k A_ki * gamma_k - g_i^0(T) / RT )

        Eigen::VectorXd gamma = Eigen::VectorXd::Zero(NE);
        // Initial guess for gamma: rough estimate from total moles
        double n_tot_guess = n_orig.sum();
        for (int k=0; k<NE; ++k) gamma(k) = std::log(std::max(b(k), 1e-10) / std::max(n_tot_guess, 1e-10)); 

        double T = T_orig;
        const int max_iter = 100;
        const double tolerance = 1e-10;

        for (int iter = 0; iter < max_iter; ++iter) {
            Eigen::VectorXd n(NSPECIES);
            Eigen::VectorXd mu0(NSPECIES); // Standard state chemical potential (Gibbs energy)
            for (int i = 0; i < NSPECIES; ++i) {
                // mu0_i = g0_i = h0_i - T*s0_i
                mu0(i) = gas->enthalpy_mol_species(T, i) - T * gas->entropy_mol_species(T, i);
                double arg = - mu0(i) / (R * T);
                for (int k = 0; k < NE; ++k) arg += atom_matrix(k, i) * gamma(k);
                n(i) = (P0 / (R * T)) * std::exp(std::clamp(arg, -100.0, 100.0));
            }

            Eigen::VectorXd f(NE + 1);
            Eigen::MatrixXd J(NE + 1, NE + 1);
            f.setZero(); J.setZero();

            // Element conservation residuals
            for (int k = 0; k < NE; ++k) {
                for (int i = 0; i < NSPECIES; ++i) f(k) += atom_matrix(k, i) * n(i);
                f(k) -= b(k);
            }

            // Energy conservation residual (scaled)
            double e_total = 0.0;
            for (int i = 0; i < NSPECIES; ++i) e_total += n(i) * (gas->enthalpy_mol_species(T, i) - R * T);
            f(NE) = (e_total - e_total_0) / (R * T_orig); 

            if (f.norm() < tolerance) break;

            // Jacobian assembly
            // J_kl = df_k / d_gamma_l = sum_i A_ki * n_i * A_li
            for (int k = 0; k < NE; ++k) {
                for (int l = 0; l < NE; ++l) {
                    for (int i = 0; i < NSPECIES; ++i) J(k, l) += atom_matrix(k, i) * n(i) * atom_matrix(l, i);
                }
            }
            // J_k,NE = df_k / dlnT = T * sum_i A_ki * dn_i/dT = sum_i A_ki * n_i * (e_i / RT)
            for (int k = 0; k < NE; ++k) {
                for (int i = 0; i < NSPECIES; ++i) {
                    double e_i = gas->enthalpy_mol_species(T, i) - R * T;
                    J(k, NE) += atom_matrix(k, i) * n(i) * (e_i / (R * T));
                }
            }
            // J_NE,l = df_NE / d_gamma_l = (1/RT_orig) * sum_i e_i * n_i * A_li
            for (int l = 0; l < NE; ++l) {
                for (int i = 0; i < NSPECIES; ++i) {
                    double e_i = gas->enthalpy_mol_species(T, i) - R * T;
                    J(NE, l) += (e_i / (R * T_orig)) * n(i) * atom_matrix(l, i);
                }
            }
            // J_NE,NE = df_NE / dlnT = (T/RT_orig) * dE/dT = (T/RT_orig) * sum(n_i cv_i + e_i dn_i/dT)
            for (int i = 0; i < NSPECIES; ++i) {
                double cv_i = gas->cp_mol_species(T, i) - R;
                double e_i = gas->enthalpy_mol_species(T, i) - R * T;
                J(NE, NE) += (n(i) * cv_i * T + e_i * n(i) * (e_i / (R * T))) / (R * T_orig);
            }

            Eigen::VectorXd step = J.colPivHouseholderQr().solve(-f);
            if (std::isnan(step.norm())) break;

            double alpha = 1.0; 
            // Damping to prevent excessive temperature changes
            if (std::abs(step(NE)) > 0.5) alpha = 0.5 / std::abs(step(NE));
            
            gamma += alpha * step.head(NE);
            T *= std::exp(alpha * step(NE));
            T = std::clamp(T, 200.0, 6000.0);

            if (alpha * step.norm() < 1e-12) break; 
        }

        // 3. Output results
        std::array<double, NSPECIES> wdot;
        Eigen::VectorXd n_eq(NSPECIES);
        for (int i = 0; i < NSPECIES; ++i) {
            double mu0 = gas->enthalpy_mol_species(T, i) - T * gas->entropy_mol_species(T, i);
            double arg = - mu0 / (R * T);
            for (int k = 0; k < NE; ++k) arg += atom_matrix(k, i) * gamma(k);
            n_eq(i) = (P0 / (R * T)) * std::exp(std::clamp(arg, -100.0, 100.0));
            wdot[i] = (n_eq(i) * gas->MW[i] - rhoY[i]) / dt;
        }
        
        Set::Scalar qdot = 0.0;
        for (int i = 0; i < NSPECIES; ++i) {
            double h_i = gas->enthalpy_mol_species(T_orig, i);
            double e_i = (h_i - R * T_orig) / gas->MW[i]; // J/kg
            qdot -= wdot[i] * e_i;
        }

        return {wdot, qdot};
    }

    static void Parse(Equilibrium& value, IO::ParmParse& pp)
    {
        Util::Message(INFO,pp.getPrefix());
#ifdef ALAMO_YAML
        pp_query_required("yaml", value.yaml);
        value.kinetics = CanteraYAML::ParseYaml(value.yaml);
        
        int NS = NSPECIES;
        int NR = value.kinetics.reaction.size();
        value.nu.resize(NS, NR);
        value.nu.setZero();
        
        // Build mapping from YAML species names to simulation species indices
        std::map<std::string, int> species_map;
        for (int i = 0; i < (int)value.kinetics.species.size(); ++i) {
            species_map[value.kinetics.species[i].name] = i;
        }
        
        // Populate stoichiometric matrix
        for (int r = 0; r < NR; ++r) {
            for (auto const& [name, coeff] : value.kinetics.reaction[r].reactants) {
                if (species_map.count(name) && species_map[name] < NS) value.nu(species_map[name], r) -= coeff;
            }
            for (auto const& [name, coeff] : value.kinetics.reaction[r].products) {
                if (species_map.count(name) && species_map[name] < NS) value.nu(species_map[name], r) += coeff;
            }
        }

        // Build element list and atom matrix
        std::set<std::string> element_set;
        for (int i = 0; i < NS; ++i) {
            for (auto const& [el, count] : value.kinetics.species[i].composition) {
                element_set.insert(el);
            }
        }
        value.NE = element_set.size();
        std::vector<std::string> elements(element_set.begin(), element_set.end());
        value.atom_matrix.resize(value.NE, NS);
        value.atom_matrix.setZero();
        for (int i = 0; i < NS; ++i) {
            for (auto const& [el, count] : value.kinetics.species[i].composition) {
                auto it = std::find(elements.begin(), elements.end(), el);
                int k = std::distance(elements.begin(), it);
                value.atom_matrix(k, i) = count;
            }
        }

        Util::Message(INFO, "Equilibrium solver initialized from ", value.yaml, " with ", NR, " reactions and ", value.NE, " elements.");
#else
        (void)value;
        Util::Abort(INFO, "Alamo was not built with YAML support, but Equilibrium chemistry was requested.");
#endif
    }

};

} // namespace Chemistry
} // namespace Model

#endif

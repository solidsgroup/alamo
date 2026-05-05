#ifndef MODEL_CHEMISTRY_EQUILIBRIUM_H
#define MODEL_CHEMISTRY_EQUILIBRIUM_H

#include "AMReX_Extension.H"
#include "Model/Chemistry/Chemistry.H"
#include "Util/Chemistry.H"
#include "Model/Gas/Gas.H"
#include <eigen3/Eigen/Dense>
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
        if (dt <= 0)
        {
            std::array<double, NSPECIES> wdot;
            for (int n = 0; n < NSPECIES; ++n) wdot[n] = 0.0;
            return {wdot, 0.0};
        }
        
        if constexpr (NSPECIES <= 1)
        {
            Util::IgnoreUnused(T_orig, rhoY, gas);
            std::array<double, NSPECIES> wdot;
            for (int n = 0; n < NSPECIES; ++n) wdot[n] = 0.0;
            return {wdot, 0.0};
        }
#ifdef ALAMO_YAML
        else
        {
            const double R = Set::Constant::Rg;
            const double P0 = 101325.0; // Standard pressure
            std::array<double, NSPECIES> wdot;
            for (int n = 0; n < NSPECIES; ++n) wdot[n] = 0.0;

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
            {
                int n_active_species = 0;
                for (int i = 0; i < NSPECIES; ++i)
                {
                    if (n_orig(i) > 1.0e-30) ++n_active_species;
                }

                if (n_active_species >= NE)
                {
                    Eigen::MatrixXd A_ls(n_active_species, NE);
                    Eigen::VectorXd rhs_ls(n_active_species);
                    int row = 0;
                    for (int i = 0; i < NSPECIES; ++i)
                    {
                        if (n_orig(i) <= 1.0e-30) continue;
                        double mu0_i = gas->enthalpy_mol_species(T_orig, i) - T_orig * gas->entropy_mol_species(T_orig, i);
                        rhs_ls(row) = std::log(std::max(n_orig(i) * R * T_orig / P0, 1.0e-300)) + mu0_i / (R * T_orig);
                        for (int k = 0; k < NE; ++k) A_ls(row, k) = atom_matrix(k, i);
                        ++row;
                    }
                    gamma = A_ls.colPivHouseholderQr().solve(rhs_ls);
                    if (!std::isfinite(gamma.norm())) gamma.setZero();
                }
                else
                {
                    double n_tot_guess = n_orig.sum();
                    for (int k = 0; k < NE; ++k)
                        gamma(k) = std::log(std::max(b(k), 1e-10) / std::max(n_tot_guess, 1e-10));
                }
            }

            double T = T_orig;
            const int max_iter = 200;
            const double tolerance = 1e-10;
            bool finite_state = true;
            bool converged = false;
            bool best_state_valid = false;
            double initial_residual_norm = -1.0;
            double best_residual_norm = std::numeric_limits<double>::infinity();
            Eigen::VectorXd best_gamma = gamma;
            double best_T = T;

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

                double residual_norm = f.norm();
                if (initial_residual_norm < 0.0) initial_residual_norm = residual_norm;
                if (std::isfinite(residual_norm) && residual_norm < best_residual_norm)
                {
                    best_residual_norm = residual_norm;
                    best_gamma = gamma;
                    best_T = T;
                    best_state_valid = true;
                }
                if (residual_norm < tolerance)
                {
                    converged = true;
                    break;
                }

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
                if (!std::isfinite(step.norm()))
                {
                    finite_state = false;
                    break;
                }

                double alpha = 1.0;
                if (std::abs(step(NE)) > 0.5) alpha = 0.5 / std::abs(step(NE));

                bool accepted = false;
                Eigen::VectorXd gamma_trial(NE);
                Eigen::VectorXd n_trial(NSPECIES);
                Eigen::VectorXd f_trial(NE + 1);
                double T_trial = T;

                for (int line_iter = 0; line_iter < 12; ++line_iter)
                {
                    gamma_trial = gamma + alpha * step.head(NE);
                    T_trial = std::clamp(T * std::exp(alpha * step(NE)), 200.0, 6000.0);

                    if (!std::isfinite(gamma_trial.norm()) || !std::isfinite(T_trial))
                    {
                        alpha *= 0.5;
                        continue;
                    }

                    for (int i = 0; i < NSPECIES; ++i)
                    {
                        double mu0_i = gas->enthalpy_mol_species(T_trial, i) - T_trial * gas->entropy_mol_species(T_trial, i);
                        double arg = - mu0_i / (R * T_trial);
                        for (int k = 0; k < NE; ++k) arg += atom_matrix(k, i) * gamma_trial(k);
                        n_trial(i) = (P0 / (R * T_trial)) * std::exp(std::clamp(arg, -100.0, 100.0));
                    }

                    f_trial.setZero();
                    for (int k = 0; k < NE; ++k)
                    {
                        for (int i = 0; i < NSPECIES; ++i) f_trial(k) += atom_matrix(k, i) * n_trial(i);
                        f_trial(k) -= b(k);
                    }
                    double e_total_trial = 0.0;
                    for (int i = 0; i < NSPECIES; ++i)
                        e_total_trial += n_trial(i) * (gas->enthalpy_mol_species(T_trial, i) - R * T_trial);
                    f_trial(NE) = (e_total_trial - e_total_0) / (R * T_orig);

                    double trial_residual_norm = f_trial.norm();
                    if (std::isfinite(trial_residual_norm) && trial_residual_norm < best_residual_norm)
                    {
                        best_residual_norm = trial_residual_norm;
                        best_gamma = gamma_trial;
                        best_T = T_trial;
                        best_state_valid = true;
                    }

                    if (std::isfinite(trial_residual_norm) && trial_residual_norm < residual_norm)
                    {
                        accepted = true;
                        break;
                    }

                    alpha *= 0.5;
                }

                if (!accepted)
                {
                    finite_state = false;
                    break;
                }

                gamma = gamma_trial;
                T = T_trial;

                if (alpha * step.norm() < 1e-12) break;
            }

            if (!converged)
            {
                double acceptable_residual_norm = std::max(1.0e-8, 1.0e-6 * std::max(initial_residual_norm, 1.0));
                if (!best_state_valid || !std::isfinite(best_residual_norm) || best_residual_norm > acceptable_residual_norm)
                {
                    return {wdot, 0.0};
                }
                gamma = best_gamma;
                T = best_T;
            }
            else if (!finite_state)
            {
                return {wdot, 0.0};
            }

            // 3. Output results
            Eigen::VectorXd n_eq(NSPECIES);
            for (int i = 0; i < NSPECIES; ++i) {
                double mu0 = gas->enthalpy_mol_species(T, i) - T * gas->entropy_mol_species(T, i);
                double arg = - mu0 / (R * T);
                for (int k = 0; k < NE; ++k) arg += atom_matrix(k, i) * gamma(k);
                if (!std::isfinite(arg)) return {wdot, 0.0};
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
#else
        else
        {
            Util::IgnoreUnused(T_orig, rhoY, gas);
            std::array<double, NSPECIES> wdot;
            for (int n = 0; n < NSPECIES; ++n) wdot[n] = 0.0;
            return {wdot, 0.0};
        }
#endif
    }

    static void Parse(Equilibrium& value, IO::ParmParse& pp)
    {
        Util::Message(INFO,pp.getPrefix());
#ifdef ALAMO_YAML
        if constexpr (NSPECIES > 1)
        {
            pp_query_required("yaml", value.yaml);
            value.kinetics = CanteraYAML::ParseYaml(value.yaml);
            ChemUtil::PreparedKinetics<NSPECIES> prepared =
                ChemUtil::PrepareKinetics<NSPECIES>(value.kinetics, "Equilibrium chemistry");

            int NS = NSPECIES;
            int NR = prepared.NR;
            value.nu.resize(NS, NR);
            value.nu.setZero();

            for (int r = 0; r < NR; ++r)
            {
                for (int i = 0; i < NS; ++i)
                {
                    value.nu(i, r) = prepared.nu_net[r][i];
                }
            }

            value.NE = prepared.NE;
            Util::Assert(INFO, TEST(value.NE > 0),
                "Equilibrium chemistry YAML did not define any elemental composition for the active phase species.");
            value.atom_matrix.resize(value.NE, NS);
            value.atom_matrix.setZero();

            for (int k = 0; k < value.NE; ++k)
            {
                for (int i = 0; i < NS; ++i)
                {
                    value.atom_matrix(k, i) = prepared.atom_matrix[k][i];
                }
            }

            Util::Message(INFO, "Equilibrium solver initialized from ", value.yaml, " with ", NR, " reactions and ", value.NE, " elements.");
        }
#else
        (void)value;
        Util::Abort(INFO, "Alamo was not built with YAML support, but Equilibrium chemistry was requested.");
#endif
    }

};

} // namespace Chemistry
} // namespace Model

#endif

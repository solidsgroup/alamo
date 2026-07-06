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

    enum class IntegrationType { EXPLICIT, SUBSTEP, IMPLICIT };
    IntegrationType integration = IntegrationType::EXPLICIT;
    int nsubsteps = 0;
    int max_iter = 0;
    double reltol = 0.0;
    double abstol = 0.0;

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

    AMREX_FORCE_INLINE
    static double norm2(const std::array<double, NSPECIES>& values)
    {
        double norm = 0.0;
        for (int i = 0; i < NSPECIES; ++i) norm += values[i] * values[i];
        return std::sqrt(norm);
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static double sum_species(const std::array<double, NSPECIES>& values)
    {
        double sum = 0.0;
        for (int i = 0; i < NSPECIES; ++i) sum += values[i];
        return sum;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static bool sanitize_species_densities(
        const std::array<double, NSPECIES>& input,
        std::array<double, NSPECIES>& output,
        const char* context)
    {
        double density_scale = 0.0;
        for (int i = 0; i < NSPECIES; ++i)
        {
            if (!std::isfinite(input[i]))
                Util::Abort(INFO, context, " received non-finite species density for species ", i, ": ", input[i]);
            density_scale += std::abs(input[i]);
        }

        const double density_tol = 1.0e-12 * ((density_scale > 1.0) ? density_scale : 1.0);
        double clipped_density = 0.0;
        for (int i = 0; i < NSPECIES; ++i)
        {
            output[i] = input[i];
            if (output[i] < 0.0)
            {
                if (output[i] < -density_tol)
                    Util::Abort(INFO, context, " received negative species density for species ", i,
                        ": rhoY=", output[i], ", tolerance=", density_tol);
                output[i] = 0.0;
            }
            clipped_density += output[i];
        }

        return clipped_density > 0.0;
    }

    AMREX_FORCE_INLINE
    static bool project_nonnegative_mass(std::array<double, NSPECIES>& values, const double target_mass)
    {
        if (target_mass < 0.0) return false;

        double positive_sum = 0.0;
        for (int i = 0; i < NSPECIES; ++i)
        {
            if (values[i] < 0.0) values[i] = 0.0;
            positive_sum += values[i];
        }

        if (target_mass == 0.0)
        {
            for (int i = 0; i < NSPECIES; ++i) values[i] = 0.0;
            return true;
        }
        if (!(positive_sum > 0.0)) return false;

        const double scale = target_mass / positive_sum;
        for (int i = 0; i < NSPECIES; ++i) values[i] *= scale;
        return true;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar compute_qdot_from_wdot(
        const Set::Scalar T,
        const std::array<double, NSPECIES>& wdot,
        const Model::Gas::Gas* gas) const
    {
        const double R = Set::Constant::Rg;
        Set::Scalar qdot = 0.0;
        for (int n = 0; n < NSPECIES; ++n)
        {
            double e_n = (gas->enthalpy_mol_species(T, n) - R * T) / gas->MW[n];
            qdot -= wdot[n] * e_n;
        }
        return qdot;
    }

    AMREX_FORCE_INLINE
    double compute_backward_euler_residual(
        const Set::Scalar P,
        const Set::Scalar T,
        const std::array<double, NSPECIES>& rhoY_old,
        const std::array<double, NSPECIES>& rhoY,
        const Set::Scalar dt,
        const Model::Gas::Gas* gas,
        std::array<double, NSPECIES>& G,
        std::array<double, NSPECIES>& wdot) const
    {
        auto source = evaluate_source(P, T, rhoY, gas);
        wdot = source.first;

        double residual_norm = 0.0;
        for (int i = 0; i < NSPECIES; ++i)
        {
            G[i] = rhoY[i] - rhoY_old[i] - dt * wdot[i];
            residual_norm += G[i] * G[i];
        }
        return std::sqrt(residual_norm);
    }

    AMREX_FORCE_INLINE
    static bool solve_linear_system(double A[NSPECIES][NSPECIES], double b[NSPECIES])
    {
        for (int i = 0; i < NSPECIES; ++i)
        {
            int pivot = i;
            for (int j = i + 1; j < NSPECIES; ++j)
                if (std::abs(A[j][i]) > std::abs(A[pivot][i])) pivot = j;
            
            for (int j = 0; j < NSPECIES; ++j) { double tmp = A[i][j]; A[i][j] = A[pivot][j]; A[pivot][j] = tmp; }
            { double tmp = b[i]; b[i] = b[pivot]; b[pivot] = tmp; }

            if (std::abs(A[i][i]) < 1.0e-30) return false;

            for (int j = i + 1; j < NSPECIES; ++j)
            {
                double factor = A[j][i] / A[i][i];
                A[j][i] = 0.0;
                for (int k = i + 1; k < NSPECIES; ++k) A[j][k] -= factor * A[i][k];
                b[j] -= factor * b[i];
            }
        }

        for (int i = NSPECIES - 1; i >= 0; --i)
        {
            for (int j = i + 1; j < NSPECIES; ++j) b[i] -= A[i][j] * b[j];
            if (std::abs(A[i][i]) > 1.0e-30) b[i] /= A[i][i];
            else return false;
        }
        return true;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    std::pair<std::array<double,NSPECIES>,Set::Scalar>
    compute_substep(
        const Set::Scalar P,
        const Set::Scalar T,
        const std::array<double, NSPECIES>& rhoY_old,
        const Set::Scalar dt,
        const Model::Gas::Gas* gas) const
    {
        if (!(dt > 0.0)) return evaluate_source(P, T, rhoY_old, gas);

        Util::Assert(INFO, TEST(nsubsteps > 0),
            "Finite-rate substep integration requires nsubsteps > 0.");

        std::array<double, NSPECIES> rhoY = rhoY_old;
        Set::Scalar dt_sub = dt / (Set::Scalar)nsubsteps;
        const double mass_old = sum_species(rhoY_old);

        for (int s = 0; s < nsubsteps; ++s)
        {
            auto source = evaluate_source(P, T, rhoY, gas);
            limit_wdot_positivity(source.first, rhoY, dt_sub);
            for (int n = 0; n < NSPECIES; ++n)
            {
                rhoY[n] += source.first[n] * dt_sub;
            }
            Util::Assert(INFO, TEST(project_nonnegative_mass(rhoY, mass_old)),
                "Finite-rate substep integration failed to project species densities.");
        }

        std::array<double, NSPECIES> wdot_eff;
        for (int n = 0; n < NSPECIES; ++n)
        {
            wdot_eff[n] = (rhoY[n] - rhoY_old[n]) / dt;
        }

        const double mass_change = sum_species(rhoY) - mass_old;
        const double mass_tol = 1.0e-10 * std::max(1.0, std::abs(mass_old));
        Util::Assert(INFO, TEST(std::abs(mass_change) <= mass_tol),
            "Finite-rate substep integration produced a non-conservative update. Residual = ",
            mass_change, ", tolerance = ", mass_tol, ".");

        return {wdot_eff, compute_qdot_from_wdot(T, wdot_eff, gas)};
    }

    AMREX_FORCE_INLINE
    std::pair<std::array<double,NSPECIES>,Set::Scalar>
    compute_implicit(
        const Set::Scalar P,
        const Set::Scalar T,
        const std::array<double, NSPECIES>& rhoY_old,
        const Set::Scalar dt,
        const Model::Gas::Gas* gas) const
    {
        if (!(dt > 0.0)) return evaluate_source(P, T, rhoY_old, gas);

        Util::Assert(INFO, TEST(max_iter > 0),
            "Finite-rate implicit integration requires max_iter > 0.");

        std::array<double, NSPECIES> rhoY = rhoY_old;
        std::array<double, NSPECIES> G;
        std::array<double, NSPECIES> wdot;
        const double state_scale = std::max(1.0, norm2(rhoY_old));
        const double mass_old = sum_species(rhoY_old);
        double res_norm = compute_backward_euler_residual(P, T, rhoY_old, rhoY, dt, gas, G, wdot);
        const double residual_scale = std::max(state_scale, res_norm);
        const double conv_tol = abstol + reltol * residual_scale;
        bool converged = res_norm <= conv_tol;
        
        for (int iter = 0; iter < max_iter; ++iter)
        {
            if (converged) break;

            double J[NSPECIES][NSPECIES];
            double b[NSPECIES];
            for (int j = 0; j < NSPECIES; ++j)
            {
                double eps = std::sqrt(1.0e-16) * std::max(std::abs(rhoY[j]), 1.0e-8 * state_scale);
                std::array<double, NSPECIES> rhoY_pert = rhoY;
                rhoY_pert[j] += eps;
                auto source_pert = evaluate_source(P, T, rhoY_pert, gas);
                for (int i = 0; i < NSPECIES; ++i)
                {
                    double dwdot_drhoY = (source_pert.first[i] - wdot[i]) / eps;
                    J[i][j] = (i == j ? 1.0 : 0.0) - dt * dwdot_drhoY;
                }
                b[j] = -G[j];
            }

            Util::Assert(INFO, TEST(solve_linear_system(J, b)),
                "Finite-rate implicit integration encountered a singular Newton system.");

            bool accepted = false;
            double alpha = 1.0;
            for (int ls = 0; ls < 20; ++ls)
            {
                std::array<double, NSPECIES> rhoY_trial = rhoY;
                for (int i = 0; i < NSPECIES; ++i)
                {
                    rhoY_trial[i] += alpha * b[i];
                }

                const bool projected = project_nonnegative_mass(rhoY_trial, mass_old);
                std::array<double, NSPECIES> G_trial;
                std::array<double, NSPECIES> wdot_trial;
                const double trial_res_norm =
                    compute_backward_euler_residual(P, T, rhoY_old, rhoY_trial, dt, gas, G_trial, wdot_trial);

                if (projected && (trial_res_norm <= conv_tol || trial_res_norm < res_norm))
                {
                    rhoY = rhoY_trial;
                    G = G_trial;
                    wdot = wdot_trial;
                    res_norm = trial_res_norm;
                    converged = res_norm <= conv_tol;
                    accepted = true;
                    break;
                }

                alpha *= 0.5;
            }

            Util::Assert(INFO, TEST(accepted),
                "Finite-rate implicit integration failed to find a positivity-preserving Newton update. ",
                "Residual = ", res_norm, ", tolerance = ", conv_tol, ".");
        }

        Util::Assert(INFO, TEST(converged),
            "Finite-rate implicit integration failed to converge after ", max_iter,
            " Newton iterations. Residual = ", res_norm, ", tolerance = ", conv_tol, ".");

        const double mass_change = sum_species(rhoY) - mass_old;
        const double mass_tol =
            std::sqrt((double)NSPECIES) * conv_tol + 1.0e-12 * std::max(1.0, std::abs(mass_old));
        Util::Assert(INFO, TEST(std::abs(mass_change) <= mass_tol),
            "Finite-rate implicit integration produced a non-conservative update. Residual = ",
            mass_change, ", tolerance = ", mass_tol, ".");

        std::array<double, NSPECIES> wdot_eff;
        for (int n = 0; n < NSPECIES; ++n)
        {
            wdot_eff[n] = (rhoY[n] - rhoY_old[n]) / dt;
        }

        return {wdot_eff, compute_qdot_from_wdot(T, wdot_eff, gas)};
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

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    std::pair<std::array<double,NSPECIES>,Set::Scalar>
    compute(
        const Set::Scalar P,
        const Set::Scalar T,
        const std::array<double, NSPECIES>& rhoY,
        const Set::Scalar dt,
        const Model::Gas::Gas* gas) const
    {
        std::array<double, NSPECIES> rhoY_valid;
        if (!sanitize_species_densities(rhoY, rhoY_valid, "Finite-rate chemistry"))
        {
            std::array<double, NSPECIES> wdot;
            for (int n = 0; n < NSPECIES; ++n) wdot[n] = 0.0;
            return {wdot, 0.0};
        }

        if (integration == IntegrationType::SUBSTEP)
        {
            return compute_substep(P, T, rhoY_valid, dt, gas);
        }
        else if (integration == IntegrationType::IMPLICIT)
        {
#ifdef AMREX_USE_GPU
            Util::Abort(INFO,
                "Finite-rate implicit integration is not supported inside GPU chemistry kernels. ",
                "Use substep or explicit integration.");
            return evaluate_source(P, T, rhoY_valid, gas);
#else
            return compute_implicit(P, T, rhoY_valid, dt, gas);
#endif
        }
        else
        {
            std::pair<std::array<double, NSPECIES>, Set::Scalar> source = evaluate_source(P, T, rhoY_valid, gas);
            if (dt > 0.0)
            {
                limit_wdot_positivity(source.first, rhoY_valid, dt);
                source.second = compute_qdot_from_wdot(T, source.first, gas);
            }
            return source;
        }
    }

    static void Parse(FiniteRate& value, IO::ParmParse& pp)
    {
        Util::Message(INFO,pp.getPrefix());

        std::string integration_type_str = "explicit";
        pp.query("integration", integration_type_str);
        if (integration_type_str == "explicit") value.integration = IntegrationType::EXPLICIT;
        else if (integration_type_str == "substep") value.integration = IntegrationType::SUBSTEP;
        else if (integration_type_str == "implicit") value.integration = IntegrationType::IMPLICIT;
        else Util::Abort(INFO, "Invalid integration type: ", integration_type_str);

        pp.query_default("nsubsteps", value.nsubsteps, 10);
        pp.query_default("max_iter", value.max_iter, 100);
        pp.query_default("reltol", value.reltol, 1e-6);
        pp.query_default("abstol", value.abstol, 1e-12);

        Util::Assert(INFO, TEST(value.nsubsteps > 0),
            "chemistry.finite_rate.nsubsteps must be > 0.");
        Util::Assert(INFO, TEST(value.max_iter > 0),
            "chemistry.finite_rate.max_iter must be > 0.");
        Util::Assert(INFO, TEST(value.reltol >= 0.0),
            "chemistry.finite_rate.reltol must be >= 0.");
        Util::Assert(INFO, TEST(value.abstol >= 0.0),
            "chemistry.finite_rate.abstol must be >= 0.");
        Util::Assert(INFO, TEST(value.reltol > 0.0 || value.abstol > 0.0),
            "At least one of chemistry.finite_rate.reltol or chemistry.finite_rate.abstol must be > 0.");

#ifdef AMREX_USE_GPU
        Util::Assert(INFO, TEST(value.integration != IntegrationType::IMPLICIT),
            "Finite-rate implicit integration currently uses a dense per-cell finite-difference Jacobian ",
            "and is not supported on GPU builds. Use substep or explicit integration.");
#endif

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
                " reactions using ", integration_type_str, " integration.");

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

#ifndef MODEL_GAS_EOS_TPG_H_
#define MODEL_GAS_EOS_TPG_H_

#include <cmath>
#include <memory>
#include "IO/ParmParse.H"
#include "Model/Gas/Gas.H"
#include "Model/Gas/EOS/EOS.H"
#include "Util/Util.H"

namespace Model {
namespace Gas {

namespace EOS {

namespace
{
    constexpr double sensible_reference_temperature = 298.15;
}

class TPG
{
public:
    static constexpr const char* name = "tpg";

public:
    static void Parse(TPG & /*value*/, IO::ParmParse & /*pp*/) {}

    template<typename GAS>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    double ComputeT(const GAS &gas, double density, double momentumx, double momentumy, double E, double Tguess,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol) const
    {
        if (!std::isfinite(density) || density <= 0.0)
            Util::Abort(INFO, "Invalid density in TPG::ComputeT: ", density);
        if (!std::isfinite(E) || !std::isfinite(momentumx) || !std::isfinite(momentumy))
            Util::Abort(INFO, "Invalid conserved state in TPG::ComputeT.");
        if (!std::isfinite(rtol) || rtol <= 0.0) rtol = 1e-12;
    
        const double inv_density = 1.0 / density;
        const double e_int = E * inv_density
            - 0.5 * (momentumx*momentumx + momentumy*momentumy) * inv_density * inv_density;
        const double mw = gas.GetMW(X, i, j, k);
        if (!std::isfinite(mw) || mw <= 0.0)
            Util::Abort(INFO, "Invalid mixture molecular weight in TPG::ComputeT: ", mw);
    
        const double inv_mw = 1.0 / mw;
        const double e_ref = (gas.enthalpy_mol(sensible_reference_temperature, X, i, j, k)
            - Set::Constant::Rg * sensible_reference_temperature) * inv_mw;
        const double T_floor = 1.0e-6;
        const double T_ceiling = 1.0e7;
    
        auto eval = [&](double T, double &f, double &cv) {
            const double h_mol = gas.enthalpy_mol(T, X, i, j, k);
            const double cp_mol = gas.cp_mol(T, X, i, j, k);
            cv = (cp_mol - Set::Constant::Rg) * inv_mw;
            f = (h_mol - Set::Constant::Rg * T) * inv_mw - e_ref - e_int;
            return std::isfinite(f) && std::isfinite(cv) && cv > 0.0;
        };
    
        auto residual_converged = [&](double f, double cv, double T) {
            const double T_scale = (T > 1.0) ? T : 1.0;
            double e_scale = std::fabs(cv) * T_scale;
            if (e_scale < 1.0) e_scale = 1.0;
            return std::fabs(f) <= rtol * e_scale;
        };
    
        auto converged = [&](double dT, double f, double cv, double T) {
            const double T_scale = (T > 1.0) ? T : 1.0;
            return std::fabs(dT) <= rtol * T_scale || residual_converged(f, cv, T);
        };
    
        double T = (std::isfinite(Tguess) && Tguess > T_floor) ? Tguess : sensible_reference_temperature;
        if (T > T_ceiling) T = T_ceiling;
    
        double f = 0.0, cv = 0.0;
        if (!eval(T, f, cv))
        {
            T = sensible_reference_temperature;
            if (!eval(T, f, cv))
                Util::Abort(INFO, "Unable to evaluate thermodynamic state in TPG::ComputeT at T=", T);
        }
        if (residual_converged(f, cv, T)) return T;
    
        double T_lo = T, T_hi = T, f_lo = f, f_hi = f;
        if (f > 0.0)
        {
            T_hi = T;
            f_hi = f;
            for (int n = 0; n < 64 && f_lo > 0.0 && T_lo > T_floor; ++n)
            {
                T_lo *= 0.5;
                if (T_lo < T_floor) T_lo = T_floor;
                double cv_lo = 0.0;
                if (!eval(T_lo, f_lo, cv_lo))
                    Util::Abort(INFO, "Unable to bracket temperature in TPG::ComputeT at T=", T_lo);
            }
            if (f_lo > 0.0) return T_floor;
        }
        else
        {
            T_lo = T;
            f_lo = f;
            for (int n = 0; n < 64 && f_hi < 0.0 && T_hi < T_ceiling; ++n)
            {
                T_hi *= 2.0;
                if (T_hi > T_ceiling) T_hi = T_ceiling;
                double cv_hi = 0.0;
                if (!eval(T_hi, f_hi, cv_hi))
                    Util::Abort(INFO, "Unable to bracket temperature in TPG::ComputeT at T=", T_hi);
            }
            if (f_hi < 0.0)
                Util::Abort(INFO, "Unable to bracket temperature in TPG::ComputeT. T_hi=", T_hi, " residual=", f_hi);
        }
    
        for (int counter = 0; counter < 64; ++counter)
        {
            if (f > 0.0)
            {
                T_hi = T;
                f_hi = f;
            }
            else
            {
                T_lo = T;
                f_lo = f;
            }
    
            double T_new = T - f / cv;
            if (!std::isfinite(T_new) || T_new <= T_lo || T_new >= T_hi)
                T_new = 0.5 * (T_lo + T_hi);
    
            double f_new = 0.0, cv_new = 0.0;
            if (!eval(T_new, f_new, cv_new))
            {
                T_new = 0.5 * (T_lo + T_hi);
                if (!eval(T_new, f_new, cv_new))
                    Util::Abort(INFO, "Unable to evaluate thermodynamic state in TPG::ComputeT at T=", T_new);
            }
    
            if (converged(T_new - T, f_new, cv_new, T_new)) return T_new;
    
            T = T_new;
            f = f_new;
            cv = cv_new;
        }
    
        Util::Abort(INFO, "Temperature didn't converge after 64 iterations. T_lo: ",
            T_lo, " T_hi: ", T_hi, " f_lo: ", f_lo, " f_hi: ", f_hi);
        return 0.5 * (T_lo + T_hi);
    }

    template<typename GAS>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    double ComputeT(const GAS &gas, double density, double momentumx, double momentumy, double momentumz, double E, double Tguess,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol) const
    {
        if (!std::isfinite(density) || density <= 0.0)
            Util::Abort(INFO, "Invalid density in TPG::ComputeT: ", density);
        if (!std::isfinite(E) || !std::isfinite(momentumx) || !std::isfinite(momentumy) || !std::isfinite(momentumz))
            Util::Abort(INFO, "Invalid conserved state in TPG::ComputeT.");
        if (!std::isfinite(rtol) || rtol <= 0.0) rtol = 1e-12;
    
        const double inv_density = 1.0 / density;
        const double e_int = E * inv_density
            - 0.5 * (momentumx*momentumx + momentumy*momentumy + momentumz*momentumz) * inv_density * inv_density;
        const double mw = gas.GetMW(X, i, j, k);
        if (!std::isfinite(mw) || mw <= 0.0)
            Util::Abort(INFO, "Invalid mixture molecular weight in TPG::ComputeT: ", mw);
    
        const double inv_mw = 1.0 / mw;
        const double e_ref = (gas.enthalpy_mol(sensible_reference_temperature, X, i, j, k)
            - Set::Constant::Rg * sensible_reference_temperature) * inv_mw;
        const double T_floor = 1.0e-6;
        const double T_ceiling = 1.0e7;
    
        auto eval = [&](double T, double &f, double &cv) {
            const double h_mol = gas.enthalpy_mol(T, X, i, j, k);
            const double cp_mol = gas.cp_mol(T, X, i, j, k);
            cv = (cp_mol - Set::Constant::Rg) * inv_mw;
            f = (h_mol - Set::Constant::Rg * T) * inv_mw - e_ref - e_int;
            return std::isfinite(f) && std::isfinite(cv) && cv > 0.0;
        };
    
        auto residual_converged = [&](double f, double cv, double T) {
            const double T_scale = (T > 1.0) ? T : 1.0;
            double e_scale = std::fabs(cv) * T_scale;
            if (e_scale < 1.0) e_scale = 1.0;
            return std::fabs(f) <= rtol * e_scale;
        };
    
        auto converged = [&](double dT, double f, double cv, double T) {
            const double T_scale = (T > 1.0) ? T : 1.0;
            return std::fabs(dT) <= rtol * T_scale || residual_converged(f, cv, T);
        };
    
        double T = (std::isfinite(Tguess) && Tguess > T_floor) ? Tguess : sensible_reference_temperature;
        if (T > T_ceiling) T = T_ceiling;
    
        double f = 0.0, cv = 0.0;
        if (!eval(T, f, cv))
        {
            T = sensible_reference_temperature;
            if (!eval(T, f, cv))
                Util::Abort(INFO, "Unable to evaluate thermodynamic state in TPG::ComputeT at T=", T);
        }
        if (residual_converged(f, cv, T)) return T;
    
        double T_lo = T, T_hi = T, f_lo = f, f_hi = f;
        if (f > 0.0)
        {
            T_hi = T;
            f_hi = f;
            for (int n = 0; n < 64 && f_lo > 0.0 && T_lo > T_floor; ++n)
            {
                T_lo *= 0.5;
                if (T_lo < T_floor) T_lo = T_floor;
                double cv_lo = 0.0;
                if (!eval(T_lo, f_lo, cv_lo))
                    Util::Abort(INFO, "Unable to bracket temperature in TPG::ComputeT at T=", T_lo);
            }
            if (f_lo > 0.0) return T_floor;
        }
        else
        {
            T_lo = T;
            f_lo = f;
            for (int n = 0; n < 64 && f_hi < 0.0 && T_hi < T_ceiling; ++n)
            {
                T_hi *= 2.0;
                if (T_hi > T_ceiling) T_hi = T_ceiling;
                double cv_hi = 0.0;
                if (!eval(T_hi, f_hi, cv_hi))
                    Util::Abort(INFO, "Unable to bracket temperature in TPG::ComputeT at T=", T_hi);
            }
            if (f_hi < 0.0)
                Util::Abort(INFO, "Unable to bracket temperature in TPG::ComputeT. T_hi=", T_hi, " residual=", f_hi);
        }
    
        for (int counter = 0; counter < 64; ++counter)
        {
            if (f > 0.0)
            {
                T_hi = T;
                f_hi = f;
            }
            else
            {
                T_lo = T;
                f_lo = f;
            }
    
            double T_new = T - f / cv;
            if (!std::isfinite(T_new) || T_new <= T_lo || T_new >= T_hi)
                T_new = 0.5 * (T_lo + T_hi);
    
            double f_new = 0.0, cv_new = 0.0;
            if (!eval(T_new, f_new, cv_new))
            {
                T_new = 0.5 * (T_lo + T_hi);
                if (!eval(T_new, f_new, cv_new))
                    Util::Abort(INFO, "Unable to evaluate thermodynamic state in TPG::ComputeT at T=", T_new);
            }
    
            if (converged(T_new - T, f_new, cv_new, T_new)) return T_new;
    
            T = T_new;
            f = f_new;
            cv = cv_new;
        }
    
        Util::Abort(INFO, "Temperature didn't converge after 64 iterations. T_lo: ",
            T_lo, " T_hi: ", T_hi, " f_lo: ", f_lo, " f_hi: ", f_hi);
        return 0.5 * (T_lo + T_hi);
    }

    template<typename GAS>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    double ComputeE(const GAS &gas, double density, double momentumx, double momentumy, double T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        double h_ref = gas.enthalpy_mass(sensible_reference_temperature, X, i, j, k);
        double e_ref = h_ref - gas.R(X, i, j, k) * sensible_reference_temperature;
        double h = gas.enthalpy_mass(T, X, i, j, k);
        double e = h - gas.R(X, i, j, k) * T - e_ref;
        double rhoE = density * e;
        double E = rhoE + 0.5*(momentumx*momentumx + momentumy*momentumy)/density;
        return E;
    }

    template<typename GAS>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    double ComputeE(const GAS &gas, double density, double momentumx, double momentumy, double momentumz, double T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        double h_ref = gas.enthalpy_mass(sensible_reference_temperature, X, i, j, k);
        double e_ref = h_ref - gas.R(X, i, j, k) * sensible_reference_temperature;
        double h = gas.enthalpy_mass(T, X, i, j, k);
        double e = h - gas.R(X, i, j, k) * T - e_ref;
        double rhoE = density * e;
        double E = rhoE + 0.5*(momentumx*momentumx + momentumy*momentumy + momentumz*momentumz)/density;
        return E;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    double ComputeT_from_primitives(double pressure, double density,
                    double R) const
    {
        return pressure / density / R;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    double ComputeP(double density, double T,
                    double R) const
    {
        return density * R * T;
    }

}; // class TPG

} // namespace EOS
} // namespace Gas
} // namespace Model

#endif

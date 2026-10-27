//
// :code:Model::Gas::EOS::TPG:code: implements a thermally perfect ideal-gas equation
// of state.
//
// The pressure relation remains :code:P = rho R T:code:, while conserved energy uses
// temperature-dependent thermodynamic properties from the gas model. Temperature
// recovery solves the nonlinear internal-energy equation with safeguarded Newton
// iteration and bracketing.
//

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
    constexpr Set::Scalar sensible_reference_temperature = 298.15;
}

class TPG
{
public:
    static constexpr const char* name = "tpg";

public:
    static void Parse(TPG & /*value*/, IO::ParmParse & /*pp*/) {}

    template<typename GAS, class Composition>
    AMREX_FORCE_INLINE
    Set::Scalar ComputeT(const GAS &gas, Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar E, Set::Scalar Tguess,
                    const Composition& X, int i, int j, int k, Set::Scalar rtol) const
    {
        if (!std::isfinite(density) || density <= 0.0)
            Util::Abort(INFO, "Invalid density in TPG::ComputeT: ", density);
        if (!std::isfinite(E) || !std::isfinite(momentumx) || !std::isfinite(momentumy))
            Util::Abort(INFO, "Invalid conserved state in TPG::ComputeT.");
        if (!std::isfinite(rtol) || rtol <= 0.0) rtol = 1e-12;

        const Set::Scalar inv_density = 1.0 / density;
        const Set::Scalar e_int = E * inv_density
            - 0.5 * (momentumx*momentumx + momentumy*momentumy) * inv_density * inv_density;
        const Set::Scalar mw = gas.GetMW(X, i, j, k);
        if (!std::isfinite(mw) || mw <= 0.0)
            Util::Abort(INFO, "Invalid mixture molecular weight in TPG::ComputeT: ", mw);

        const Set::Scalar inv_mw = 1.0 / mw;
        const Set::Scalar e_ref = (gas.enthalpy_mol(sensible_reference_temperature, X, i, j, k)
            - Set::Constant::Rg * sensible_reference_temperature) * inv_mw;
        const Set::Scalar T_floor = 1.0e-6;
        const Set::Scalar T_ceiling = 1.0e7;

        auto eval = [&](Set::Scalar T, Set::Scalar &f, Set::Scalar &cv) {
            const Set::Scalar h_mol = gas.enthalpy_mol(T, X, i, j, k);
            const Set::Scalar cp_mol = gas.cp_mol(T, X, i, j, k);
            cv = (cp_mol - Set::Constant::Rg) * inv_mw;
            f = (h_mol - Set::Constant::Rg * T) * inv_mw - e_ref - e_int;
            return std::isfinite(f) && std::isfinite(cv) && cv > 0.0;
        };

        auto residual_converged = [&](Set::Scalar f, Set::Scalar cv, Set::Scalar T) {
            const Set::Scalar T_scale = (T > 1.0) ? T : 1.0;
            Set::Scalar e_scale = std::fabs(cv) * T_scale;
            if (e_scale < 1.0) e_scale = 1.0;
            return std::fabs(f) <= rtol * e_scale;
        };

        auto converged = [&](Set::Scalar dT, Set::Scalar f, Set::Scalar cv, Set::Scalar T) {
            const Set::Scalar T_scale = (T > 1.0) ? T : 1.0;
            return std::fabs(dT) <= rtol * T_scale || residual_converged(f, cv, T);
        };

        Set::Scalar T = (std::isfinite(Tguess) && Tguess > T_floor) ? Tguess : sensible_reference_temperature;
        if (T > T_ceiling) T = T_ceiling;

        Set::Scalar f = 0.0, cv = 0.0;
        if (!eval(T, f, cv))
        {
            T = sensible_reference_temperature;
            if (!eval(T, f, cv))
                Util::Abort(INFO, "Unable to evaluate thermodynamic state in TPG::ComputeT at T=", T);
        }
        if (residual_converged(f, cv, T)) return T;

        Set::Scalar T_lo = T, T_hi = T, f_lo = f, f_hi = f;
        if (f > 0.0)
        {
            T_hi = T;
            f_hi = f;
            for (int n = 0; n < 64 && f_lo > 0.0 && T_lo > T_floor; ++n)
            {
                T_lo *= 0.5;
                if (T_lo < T_floor) T_lo = T_floor;
                Set::Scalar cv_lo = 0.0;
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
                Set::Scalar cv_hi = 0.0;
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

            Set::Scalar T_new = T - f / cv;
            if (!std::isfinite(T_new) || T_new <= T_lo || T_new >= T_hi)
                T_new = 0.5 * (T_lo + T_hi);

            Set::Scalar f_new = 0.0, cv_new = 0.0;
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

    template<typename GAS, class Composition>
    AMREX_FORCE_INLINE
    Set::Scalar ComputeT(const GAS &gas, Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz, Set::Scalar E, Set::Scalar Tguess,
                    const Composition& X, int i, int j, int k, Set::Scalar rtol) const
    {
        if (!std::isfinite(density) || density <= 0.0)
            Util::Abort(INFO, "Invalid density in TPG::ComputeT: ", density);
        if (!std::isfinite(E) || !std::isfinite(momentumx) || !std::isfinite(momentumy) || !std::isfinite(momentumz))
            Util::Abort(INFO, "Invalid conserved state in TPG::ComputeT.");
        if (!std::isfinite(rtol) || rtol <= 0.0) rtol = 1e-12;

        const Set::Scalar inv_density = 1.0 / density;
        const Set::Scalar e_int = E * inv_density
            - 0.5 * (momentumx*momentumx + momentumy*momentumy + momentumz*momentumz) * inv_density * inv_density;
        const Set::Scalar mw = gas.GetMW(X, i, j, k);
        if (!std::isfinite(mw) || mw <= 0.0)
            Util::Abort(INFO, "Invalid mixture molecular weight in TPG::ComputeT: ", mw);

        const Set::Scalar inv_mw = 1.0 / mw;
        const Set::Scalar e_ref = (gas.enthalpy_mol(sensible_reference_temperature, X, i, j, k)
            - Set::Constant::Rg * sensible_reference_temperature) * inv_mw;
        const Set::Scalar T_floor = 1.0e-6;
        const Set::Scalar T_ceiling = 1.0e7;

        auto eval = [&](Set::Scalar T, Set::Scalar &f, Set::Scalar &cv) {
            const Set::Scalar h_mol = gas.enthalpy_mol(T, X, i, j, k);
            const Set::Scalar cp_mol = gas.cp_mol(T, X, i, j, k);
            cv = (cp_mol - Set::Constant::Rg) * inv_mw;
            f = (h_mol - Set::Constant::Rg * T) * inv_mw - e_ref - e_int;
            return std::isfinite(f) && std::isfinite(cv) && cv > 0.0;
        };

        auto residual_converged = [&](Set::Scalar f, Set::Scalar cv, Set::Scalar T) {
            const Set::Scalar T_scale = (T > 1.0) ? T : 1.0;
            Set::Scalar e_scale = std::fabs(cv) * T_scale;
            if (e_scale < 1.0) e_scale = 1.0;
            return std::fabs(f) <= rtol * e_scale;
        };

        auto converged = [&](Set::Scalar dT, Set::Scalar f, Set::Scalar cv, Set::Scalar T) {
            const Set::Scalar T_scale = (T > 1.0) ? T : 1.0;
            return std::fabs(dT) <= rtol * T_scale || residual_converged(f, cv, T);
        };

        Set::Scalar T = (std::isfinite(Tguess) && Tguess > T_floor) ? Tguess : sensible_reference_temperature;
        if (T > T_ceiling) T = T_ceiling;

        Set::Scalar f = 0.0, cv = 0.0;
        if (!eval(T, f, cv))
        {
            T = sensible_reference_temperature;
            if (!eval(T, f, cv))
                Util::Abort(INFO, "Unable to evaluate thermodynamic state in TPG::ComputeT at T=", T);
        }
        if (residual_converged(f, cv, T)) return T;

        Set::Scalar T_lo = T, T_hi = T, f_lo = f, f_hi = f;
        if (f > 0.0)
        {
            T_hi = T;
            f_hi = f;
            for (int n = 0; n < 64 && f_lo > 0.0 && T_lo > T_floor; ++n)
            {
                T_lo *= 0.5;
                if (T_lo < T_floor) T_lo = T_floor;
                Set::Scalar cv_lo = 0.0;
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
                Set::Scalar cv_hi = 0.0;
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

            Set::Scalar T_new = T - f / cv;
            if (!std::isfinite(T_new) || T_new <= T_lo || T_new >= T_hi)
                T_new = 0.5 * (T_lo + T_hi);

            Set::Scalar f_new = 0.0, cv_new = 0.0;
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

    template<typename GAS, class Composition>
    AMREX_FORCE_INLINE
    Set::Scalar ComputeE(const GAS &gas, Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar T,
                    const Composition& X, int i, int j, int k) const
    {
        Set::Scalar h_ref = gas.enthalpy_mass(sensible_reference_temperature, X, i, j, k);
        Set::Scalar e_ref = h_ref - gas.R(X, i, j, k) * sensible_reference_temperature;
        Set::Scalar h = gas.enthalpy_mass(T, X, i, j, k);
        Set::Scalar e = h - gas.R(X, i, j, k) * T - e_ref;
        Set::Scalar rhoE = density * e;
        Set::Scalar E = rhoE + 0.5*(momentumx*momentumx + momentumy*momentumy)/density;
        return E;
    }

    template<typename GAS, class Composition>
    AMREX_FORCE_INLINE
    Set::Scalar ComputeE(const GAS &gas, Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz, Set::Scalar T,
                    const Composition& X, int i, int j, int k) const
    {
        Set::Scalar h_ref = gas.enthalpy_mass(sensible_reference_temperature, X, i, j, k);
        Set::Scalar e_ref = h_ref - gas.R(X, i, j, k) * sensible_reference_temperature;
        Set::Scalar h = gas.enthalpy_mass(T, X, i, j, k);
        Set::Scalar e = h - gas.R(X, i, j, k) * T - e_ref;
        Set::Scalar rhoE = density * e;
        Set::Scalar E = rhoE + 0.5*(momentumx*momentumx + momentumy*momentumy + momentumz*momentumz)/density;
        return E;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeT_from_primitives(Set::Scalar pressure, Set::Scalar density,
                    Set::Scalar R) const
    {
        return pressure / density / R;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeP(Set::Scalar density, Set::Scalar T,
                    Set::Scalar R) const
    {
        return density * R * T;
    }

}; // class TPG

} // namespace EOS
} // namespace Gas
} // namespace Model

#endif

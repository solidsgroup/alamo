#ifndef MODEL_GAS_EOS_CPG_H_
#define MODEL_GAS_EOS_CPG_H_

#include "IO/ParmParse.H"
#include "Model/Gas/EOS/EOS.H"
#include <memory>

namespace Model {
namespace Gas {

namespace EOS {

class CPG
{
public:
    static constexpr const char* name = "cpg";

public:
    static void Parse(CPG & /*value*/, IO::ParmParse & /*pp*/) {}

public:
    template<typename GAS>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeT(const GAS &gas, Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar E, Set::Scalar Tguess,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k, Set::Scalar /*rtol=1e-12*/) const
    {
        Set::Scalar P = (E - 0.5*(momentumx*momentumx + momentumy*momentumy)/density) * (gas.gamma(Tguess, X, i, j, k) - 1.0);
        Set::Scalar T = P / density / gas.R(X, i, j, k);
        return T;
    }

    template<typename GAS>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeT(const GAS &gas, Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz, Set::Scalar E, Set::Scalar Tguess,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k, Set::Scalar /*rtol=1e-12*/) const
    {
        Set::Scalar P = (E - 0.5*(momentumx*momentumx + momentumy*momentumy + momentumz*momentumz)/density) * (gas.gamma(Tguess, X, i, j, k) - 1.0);
        Set::Scalar T = P / density / gas.R(X, i, j, k);
        return T;
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

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeE(Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar T,
                    Set::Scalar R, Set::Scalar gamma) const
    {
        Set::Scalar P = density * R * T;
        Set::Scalar rhoE = P / (gamma - 1.0);
        Set::Scalar E = rhoE + 0.5*(momentumx*momentumx + momentumy*momentumy)/density;
        return E;
    }

    template<typename GAS>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeE(const GAS &gas, Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        return ComputeE(density, momentumx, momentumy, T, gas.R(X, i, j, k), gas.gamma(T, X, i, j, k));
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeE(Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz, Set::Scalar T,
                    Set::Scalar R, Set::Scalar gamma) const
    {
        Set::Scalar P = density * R * T;
        Set::Scalar rhoE = P / (gamma - 1.0);
        Set::Scalar E = rhoE + 0.5*(momentumx*momentumx + momentumy*momentumy + momentumz*momentumz)/density;
        return E;
    }

    template<typename GAS>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeE(const GAS &gas, Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz, Set::Scalar T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        return ComputeE(density, momentumx, momentumy, momentumz, T, gas.R(X, i, j, k), gas.gamma(T, X, i, j, k));
    }

}; // class CPG

} // namespace EOS
} // namespace Gas
} // namespace Model

#endif


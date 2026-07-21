//
// :code:Model::Gas::EOS::Rocfire:code: implements the reduced Rocfire gas equation
// of state.  Rocfire assumes a single molecular weight and heat capacity for all
// species, so the mixture gas constant and heat-capacity ratio are independent of
// composition.
//

#ifndef MODEL_GAS_EOS_ROCFIRE_H_
#define MODEL_GAS_EOS_ROCFIRE_H_

#include <cmath>
#include "IO/ParmParse.H"
#include "Model/Gas/EOS/EOS.H"
#include "Util/Util.H"

namespace Model {
namespace Gas {

namespace EOS {

class Rocfire
{
public:
    static constexpr const char* name = "rocfire";

public:
    static void Parse(Rocfire & /*value*/, IO::ParmParse & /*pp*/) {}

public:
    template<typename GAS>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeT(const GAS &gas, Set::Scalar density,
            Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar E,
            Set::Scalar Tguess, Set::Patch<const Set::Scalar>& X,
            int i, int j, int k, Set::Scalar /*rtol=1e-12*/) const
    {
        const Set::Scalar kinetic = 0.5 * (momentumx*momentumx + momentumy*momentumy) / density;
        const Set::Scalar P = (E - kinetic) * (gas.gamma(Tguess, X, i, j, k) - 1.0);
        return P / density / gas.R(X, i, j, k);
    }

    template<typename GAS>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeT(const GAS &gas, Set::Scalar density,
            Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz,
            Set::Scalar E, Set::Scalar Tguess, Set::Patch<const Set::Scalar>& X,
            int i, int j, int k, Set::Scalar /*rtol=1e-12*/) const
    {
        const Set::Scalar kinetic =
            0.5 * (momentumx*momentumx + momentumy*momentumy + momentumz*momentumz) / density;
        const Set::Scalar P = (E - kinetic) * (gas.gamma(Tguess, X, i, j, k) - 1.0);
        return P / density / gas.R(X, i, j, k);
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeT_from_primitives(Set::Scalar pressure, Set::Scalar density,
            Set::Scalar R) const
    {
        return pressure / density / R;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeP(Set::Scalar density, Set::Scalar T, Set::Scalar R) const
    {
        return density * R * T;
    }

    template<typename GAS>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeE(const GAS &gas, Set::Scalar density,
            Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar T,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        const Set::Scalar P = density * gas.R(X, i, j, k) * T;
        return P / (gas.gamma(T, X, i, j, k) - 1.0) +
            0.5 * (momentumx*momentumx + momentumy*momentumy) / density;
    }

    template<typename GAS>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeE(const GAS &gas, Set::Scalar density,
            Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz,
            Set::Scalar T, Set::Patch<const Set::Scalar>& X,
            int i, int j, int k) const
    {
        const Set::Scalar P = density * gas.R(X, i, j, k) * T;
        return P / (gas.gamma(T, X, i, j, k) - 1.0) +
            0.5 * (momentumx*momentumx + momentumy*momentumy + momentumz*momentumz) / density;
    }
};

} // namespace EOS
} // namespace Gas
} // namespace Model

#endif

//
// :code:Model::Gas::EOS::Rocfire:code: implements the reduced Rocfire gas equation
// of state. Rocfire assumes a single molecular weight and heat capacity for all
// species, so the mixture gas constant and heat-capacity ratio are independent of
// composition.
//

#ifndef MODEL_GAS_EOS_ROCFIRE_H_
#define MODEL_GAS_EOS_ROCFIRE_H_

#include "IO/ParmParse.H"
#include "Model/Gas/EOS/EOS.H"

namespace Model {
namespace Gas {
namespace EOS {

class Rocfire
{
public:
    static constexpr const char* name = "rocfire";

    static void Parse(Rocfire & /*value*/, IO::ParmParse & /*pp*/) {}

    template<typename GAS, class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeT(const GAS &gas, Set::Scalar density,
            Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar E,
            Set::Scalar Tguess, const Composition& X,
            int i, int j, int k, Set::Scalar /*rtol=1e-12*/) const
    {
        const Set::Scalar kinetic =
            0.5 * (momentumx*momentumx + momentumy*momentumy) / density;
        const Set::Scalar pressure =
            (E - kinetic) * (gas.gamma(Tguess, X, i, j, k) - 1.0);
        return pressure / density / gas.R(X, i, j, k);
    }

    template<typename GAS, class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeT(const GAS &gas, Set::Scalar density,
            Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz,
            Set::Scalar E, Set::Scalar Tguess, const Composition& X,
            int i, int j, int k, Set::Scalar /*rtol=1e-12*/) const
    {
        const Set::Scalar kinetic = 0.5 *
            (momentumx*momentumx + momentumy*momentumy + momentumz*momentumz) /
            density;
        const Set::Scalar pressure =
            (E - kinetic) * (gas.gamma(Tguess, X, i, j, k) - 1.0);
        return pressure / density / gas.R(X, i, j, k);
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeT_from_primitives(Set::Scalar pressure,
            Set::Scalar density, Set::Scalar R) const
    {
        return pressure / density / R;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeP(Set::Scalar density, Set::Scalar T, Set::Scalar R) const
    {
        return density * R * T;
    }

    template<typename GAS, class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeE(const GAS &gas, Set::Scalar density,
            Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar T,
            const Composition& X, int i, int j, int k) const
    {
        const Set::Scalar pressure = density * gas.R(X, i, j, k) * T;
        return pressure / (gas.gamma(T, X, i, j, k) - 1.0) +
            0.5 * (momentumx*momentumx + momentumy*momentumy) / density;
    }

    template<typename GAS, class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeE(const GAS &gas, Set::Scalar density,
            Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz,
            Set::Scalar T, const Composition& X,
            int i, int j, int k) const
    {
        const Set::Scalar pressure = density * gas.R(X, i, j, k) * T;
        return pressure / (gas.gamma(T, X, i, j, k) - 1.0) + 0.5 *
            (momentumx*momentumx + momentumy*momentumy + momentumz*momentumz) /
            density;
    }
};

} // namespace EOS
} // namespace Gas
} // namespace Model

#endif

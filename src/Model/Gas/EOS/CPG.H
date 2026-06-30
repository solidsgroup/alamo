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
    double ComputeT(const GAS &gas, double density, double momentumx, double momentumy, double E, double Tguess,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k, double /*rtol=1e-12*/) const
    {
        double P = (E - 0.5*(momentumx*momentumx + momentumy*momentumy)/density) * (gas.gamma(Tguess, X, i, j, k) - 1.0);
        double T = P / density / gas.R(X, i, j, k);
        return T;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    double ComputeT(double pressure, double density,
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

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    double ComputeE(double density, double momentumx, double momentumy, double T,
                    double R, double gamma) const
    {
        double P = density * R * T;
        double rhoE = P / (gamma - 1.0);
        double E = rhoE + 0.5*(momentumx*momentumx + momentumy*momentumy)/density;
        return E;
    }

}; // class CPG

} // namespace EOS
} // namespace Gas
} // namespace Model

#endif


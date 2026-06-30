#ifndef MODEL_GAS_TRANSPORT_ROCFIRE_H_
#define MODEL_GAS_TRANSPORT_ROCFIRE_H_

#include <array>
#include <vector>
#include <memory>
#include "IO/ParmParse.H"
#include "Model/Gas/Transport/Transport.H"

namespace Model {
namespace Gas {
namespace Transport {

class Rocfire
{
public:
    static constexpr const char* name = "rocfire";

private:
    double Rg;
    std::array<double, NSPECIES> MW;
    Thermo::Thermo<Thermo::CpConstant,Thermo::NASA7,Thermo::Rocfire> thermo;
    Set::Scalar mw, lewis, prandtl, lambda_a, lambda_b;

public:
    Rocfire()

    static void Parse(Rocfire & value, IO::ParmParse & pp)
    {
        // molecular weight
        pp.query_default("mw", value.mw, "26.0_g/mol", Unit::Mass()/Unit::Amount());
        // prandtl number
        pp.query_default("prandtl", value.prandtl, "0.69", Unit::Less());
        // lewis number
        pp.query_default("lewis", value.lewis, "1.0", Unit::Less());
        // thermal conductivity coefficient, lambda_a
        pp.query_default("lambda_a", value.lambda_a, "2.13037e-7_cal/cm/s/K^2", Unit::Energy()/Unit::Length()/Unit::Time()/Unit::Temperature()/Unit::Temperature());
        // thermal conductivity coefficient, lambda_b
        pp.query_default("lambda_b", value.lambda_b, "5.32743e-6_cal/cm/s/K", Unit::Energy()/Unit::Length()/Unit::Time()/Unit::Temperature());

        Util::Assert(INFO, TEST(value.mw > 0.0), "gas.transport.rocfire.mw must be positive.");
        Util::Assert(INFO, TEST(value.prandtl > 0.0), "gas.transport.rocfire.prandtl must be positive.");
        Util::Assert(INFO, TEST(value.lewis > 0.0), "gas.transport.rocfire.lewis must be positive.");
    }

public:
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    double mixture_mw(Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        double mw_mix = 0.0;
        for (int n = 0; n < NSPECIES; ++n) mw_mix += X(i,j,k,n) * MW[n];
        return mw_mix;
    }

    double dynamic_viscosity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        const double cp_mass = thermo->cp_mol(T, X, i, j, k) / mw;
        return prandtl * thermal_conductivity(T, X, i, j, k) / cp_mass;
    }

    double thermal_conductivity(double T, Set::Patch<const Set::Scalar>& /*X*/, int /*i*/, int /*j*/, int /*k*/) const
    {
        return lambda_a * T + lambda_b;
    }

    void diffusion_coeffs(
            Set::Patch<Set::Scalar>& DKM, double T, double P,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k)
    {
        const double kappa = thermal_conductivity(T, X, i, j, k);
        const double cp_mass = thermo->cp_mol(T, X, i, j, k) / mw;
        const double rho = (P > 0.0 && T > 0.0) ? P * mixture_mw(X, i, j, k) / (Set::Constant::Rg * T) : 0.0;
        for (int a=0; a<NSPECIES; ++a)
        {
            DKM(i, j, k, a) = (rho > 0.0) ? kappa / rho / cp_mass / lewis : 0.0;
        }
    }

}; // class Rocfire

} // namespace Transport
} // namespace Gas
} // namespace Model

#endif

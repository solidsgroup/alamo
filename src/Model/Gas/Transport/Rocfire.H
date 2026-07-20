//
// :code:Model::Gas::Transport::Rocfire:code: provides a reduced transport closure used
// with Rocfire-style propellant gas models.
//
// Thermal conductivity is modeled as a linear function of temperature. Dynamic viscosity
// is inferred from the Prandtl number and heat capacity, and species diffusion is inferred
// from the Lewis number.
//
// Reference: Gross, M. L., Hedman, T. D., Son, S. F., Jackson, T. L., &
// Beckstead, M. W. (2013). Coupling micro and meso-scale combustion models of
// AP/HTPB propellants. Combustion and Flame, 160(5), 982-992.
// https://doi.org/10.1016/j.combustflame.2013.01.016
//

#ifndef MODEL_GAS_TRANSPORT_ROCFIRE_H_
#define MODEL_GAS_TRANSPORT_ROCFIRE_H_

#include <array>
#include <vector>
#include <memory>
#include <cmath>
#include "IO/ParmParse.H"
#include "Model/Gas/Thermo/CpConstant.H"
#include "Model/Gas/Thermo/NASA7.H"
#include "Model/Gas/Thermo/Rocfire.H"
#include "Model/Gas/Transport/Transport.H"

namespace Model {
namespace Gas {
namespace Transport {

class Rocfire
{
public:
    static constexpr const char* name = "rocfire";

private:
    Set::Scalar mw, cp_mass, lewis, prandtl, lambda_a, lambda_b;

public:
    Rocfire() = default;

    static void Parse(  Rocfire & value, IO::ParmParse & pp,
                        const std::array<Set::Scalar,NSPECIES> &a_MW,
                        const Thermo::Thermo<Thermo::CpConstant,Thermo::Rocfire,Thermo::NASA7> &a_thermo)
    {
        // prandtl number
        pp.query_default("prandtl", value.prandtl, "0.69", Unit::Less());
        // lewis number
        pp.query_default("lewis", value.lewis, "1.0", Unit::Less());
        // thermal conductivity coefficient, lambda_a
        pp.query_default("lambda_a", value.lambda_a, "2.13037e-7_cal/cm/s/K^2", Unit::Energy()/Unit::Length()/Unit::Time()/Unit::Temperature()/Unit::Temperature());
        // thermal conductivity coefficient, lambda_b
        pp.query_default("lambda_b", value.lambda_b, "5.32743e-6_cal/cm/s/K", Unit::Energy()/Unit::Length()/Unit::Time()/Unit::Temperature());

        value.mw = a_MW[0];
        const Set::Scalar mw_tol = 1.0e-12 * ((std::fabs(value.mw) > 1.0) ? std::fabs(value.mw) : 1.0);
        for (int n=0; n<NSPECIES; ++n)
        {
            Util::Assert(INFO, TEST(std::isfinite(a_MW[n]) && a_MW[n] > 0.0),
                "Invalid molecular weight for species ", n, ": ", a_MW[n]);
            if (std::fabs(a_MW[n] - value.mw) > mw_tol)
            {
                Util::ParmParseException(INFO, pp.getPrefix(),
                    "transport.type=rocfire requires all gas mw entries to be equal.");
            }
        }
        value.cp_mass = a_thermo.cp_mol_species(0.0, 0) / value.mw;

        Util::Assert(INFO, TEST(value.mw > 0.0), "gas mw entries must be positive.");
        Util::Assert(INFO, TEST(value.prandtl > 0.0), "gas.transport.rocfire.prandtl must be positive.");
        Util::Assert(INFO, TEST(value.lewis > 0.0), "gas.transport.rocfire.lewis must be positive.");
        Util::Assert(INFO, TEST(value.cp_mass > 0.0),
            "gas.transport.rocfire cp_mass from thermo must be positive.");
    }

public:
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar mixture_mw(Set::Patch<const Set::Scalar>& /*X*/, int /*i*/, int /*j*/, int /*k*/) const
    {
        return mw;
    }

    AMREX_FORCE_INLINE
    Set::Scalar dynamic_viscosity(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        return prandtl * thermal_conductivity(T, X, i, j, k) / cp_mass;
    }

    AMREX_FORCE_INLINE
    Set::Scalar thermal_conductivity(Set::Scalar T, Set::Patch<const Set::Scalar>& /*X*/, int /*i*/, int /*j*/, int /*k*/) const
    {
        return lambda_a * T + lambda_b;
    }

    AMREX_FORCE_INLINE
    void diffusion_coeffs(
            Set::Patch<Set::Scalar>& DKM, Set::Scalar T, Set::Scalar P,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        const Set::Scalar kappa = thermal_conductivity(T, X, i, j, k);
        const Set::Scalar rho = (P > 0.0 && T > 0.0) ?
            P * mw / (Set::Constant::Rg * T) : 0.0;
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

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

#include <vector>
#include <memory>
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
    int nspecies = 1;
    Set::Scalar Rg;
    std::vector<double> MW;
    Thermo::Thermo<Thermo::CpConstant,Thermo::NASA7,Thermo::Rocfire> thermo;
    Set::Scalar mw, lewis, prandtl, lambda_a, lambda_b;

public:
    Rocfire()
    {
        Rg = Set::Constant::Rg;
    }

    static void Parse(  Rocfire & value, IO::ParmParse & pp,
                        int a_nspecies, const std::vector<double> &a_MW,
                        const Thermo::Thermo<Thermo::CpConstant,Thermo::NASA7,Thermo::Rocfire> &a_thermo)
    {
        value.nspecies = a_nspecies;
        value.MW = a_MW;
        value.thermo = a_thermo;
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
    template<class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar mixture_mw(const Composition& X, int i, int j, int k) const
    {
        Set::Scalar mw_mix = 0.0;
        for (int n = 0; n < nspecies; ++n) mw_mix += X(i,j,k,n) * MW[n];
        return mw_mix;
    }

    template<class Composition>
    Set::Scalar dynamic_viscosity(Set::Scalar T, const Composition& X, int i, int j, int k) const
    {
        const Set::Scalar cp_mass = thermo.cp_mol(T, X, i, j, k) / mw;
        return prandtl * thermal_conductivity(T, X, i, j, k) / cp_mass;
    }

    template<class Composition>
    Set::Scalar thermal_conductivity(Set::Scalar T, const Composition& /*X*/, int /*i*/, int /*j*/, int /*k*/) const
    {
        return lambda_a * T + lambda_b;
    }

    template<class Composition>
    Set::Scalar diffusion_coefficient(Set::Scalar T, Set::Scalar P,
                                    const Composition& X, int i, int j, int k,
                                    int /*species*/) const
    {
        const Set::Scalar kappa = thermal_conductivity(T, X, i, j, k);
        const Set::Scalar cp_mass = thermo.cp_mol(T, X, i, j, k) / mw;
        const Set::Scalar rho = (P > 0.0 && T > 0.0) ? P * mixture_mw(X, i, j, k) / (Set::Constant::Rg * T) : 0.0;
        return (rho > 0.0) ? kappa / rho / cp_mass / lewis : 0.0;
    }

    template<class Composition>
    void diffusion_coeffs(Set::Patch<Set::Scalar>& DKM, Set::Scalar T,
                        Set::Scalar P, const Composition& X,
                        int i, int j, int k) const
    {
        for (int a=0; a<nspecies; ++a)
            DKM(i,j,k,a) = diffusion_coefficient(T, P, X, i, j, k, a);
    }

}; // class Rocfire

} // namespace Transport
} // namespace Gas
} // namespace Model

#endif

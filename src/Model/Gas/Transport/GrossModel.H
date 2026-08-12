//
// :code:Model::Gas::Transport::GrossModel:code: provides a reduced transport closure used
// with GrossModel-style propellant gas models.
//
// Thermal conductivity is modeled as a linear function of temperature. Dynamic viscosity
// is inferred from the Prandtl number and heat capacity, and species diffusion is inferred
// from Lewis numbers.  gas.transport.gross_model.lewis accepts either one
// nondimensional value, broadcast to every gas species, or one value per gas
// species in gas.mw order.  Uniform values retain the common-diffusivity
// implicit solve; nonuniform values use conservative corrected species fluxes.
//
// Reference: Gross, M. L., Hedman, T. D., Son, S. F., Jackson, T. L., &
// Beckstead, M. W. (2013). Coupling micro and meso-scale combustion models of
// AP/HTPB propellants. Combustion and Flame, 160(5), 982-992.
// https://doi.org/10.1016/j.combustflame.2013.01.016
//

#ifndef MODEL_GAS_TRANSPORT_GROSS_MODEL_H_
#define MODEL_GAS_TRANSPORT_GROSS_MODEL_H_

#include <array>
#include <vector>
#include <memory>
#include <cmath>
#include "IO/ParmParse.H"
#include "Model/Gas/Thermo/CpConstant.H"
#include "Model/Gas/Thermo/NASA7.H"
#include "Model/Gas/Thermo/GrossModel.H"
#include "Model/Gas/Transport/Transport.H"

namespace Model {
namespace Gas {
namespace Transport {

class GrossModel
{
public:
    static constexpr const char* name = "gross_model";
    static constexpr bool implicit_diffusivity = true;

private:
    static constexpr int max_species = 64;
    int nspecies = 1;
    std::array<Set::Scalar, max_species> molecular_weight{};
    std::array<Set::Scalar, max_species> lewis_numbers{};
    Set::Scalar cp_mass_value, prandtl, lambda_a, lambda_b;
    bool common_lewis_number = true;

public:
    GrossModel() = default;

    static void Parse(  GrossModel & value, IO::ParmParse & pp,
                        int a_nspecies, const std::vector<double> &a_MW,
                        const Thermo::Thermo<Thermo::CpConstant,Thermo::NASA7,Thermo::GrossModel> &a_thermo)
    {
        // prandtl number
        pp.query_default("prandtl", value.prandtl, "0.69", Unit::Less());
        // One Lewis number broadcasts to all gas species; otherwise the input
        // must contain exactly one nondimensional value per gas species.
        std::vector<Set::Scalar> parsed_lewis;
        pp.queryarr_default("lewis", parsed_lewis, {1.0});
        // thermal conductivity coefficient, lambda_a
        pp.query_default("lambda_a", value.lambda_a, "2.13037e-7_cal/cm/s/K^2", Unit::Energy()/Unit::Length()/Unit::Time()/Unit::Temperature()/Unit::Temperature());
        // thermal conductivity coefficient, lambda_b
        pp.query_default("lambda_b", value.lambda_b, "5.32743e-6_cal/cm/s/K", Unit::Energy()/Unit::Length()/Unit::Time()/Unit::Temperature());

        // Schema traversal records the model inputs without constructing valid
        // parent gas and thermo data.  Avoid validating or dereferencing those
        // placeholders until parsing a real simulation input.
        if (IO::ParmParse::InTraversalMode()) return;

        value.nspecies = a_nspecies;
        Util::Assert(INFO, TEST(static_cast<int>(a_MW.size()) == value.nspecies),
            "Transport::GrossModel molecular-weight count does not match its species count.");
        Util::Assert(INFO, TEST(!a_MW.empty()),
            "gas.mw must contain at least one molecular weight.");
        Util::Assert(INFO, TEST(a_MW.size() <= max_species),
            "gas.transport.gross_model supports at most ", max_species, " species.");

        for (int n = 0; n < value.nspecies; ++n)
        {
            Util::Assert(INFO,
                TEST(std::isfinite(a_MW[n]) && a_MW[n] > 0.0),
                "Invalid molecular weight for species ", n, ": ", a_MW[n]);
            value.molecular_weight[n] = a_MW[n];
        }
        if (parsed_lewis.size() != 1 &&
            static_cast<int>(parsed_lewis.size()) != value.nspecies)
            Util::Exception(INFO,
                "gas.transport.gross_model.lewis requires either one value "
                "or exactly one value per gas species (expected ",
                value.nspecies, ", got ", parsed_lewis.size(), ")");
        value.common_lewis_number = true;
        for (int n = 0; n < value.nspecies; ++n)
        {
            value.lewis_numbers[n] = parsed_lewis.size() == 1 ?
                parsed_lewis[0] : parsed_lewis[n];
            Util::Assert(INFO,
                TEST(std::isfinite(value.lewis_numbers[n]) &&
                     value.lewis_numbers[n] > 0.0),
                "gas.transport.gross_model.lewis value for gas species ", n,
                " must be finite and positive (got ",
                value.lewis_numbers[n], ")");
            if (n > 0 && value.lewis_numbers[n] != value.lewis_numbers[0])
                value.common_lewis_number = false;
        }
        value.cp_mass_value =
            a_thermo.cp_mol_species(0.0, 0) / value.molecular_weight[0];

        Util::Assert(INFO, TEST(value.prandtl > 0.0), "gas.transport.gross_model.prandtl must be positive.");
        Util::Assert(INFO, TEST(value.cp_mass_value > 0.0),
            "gas.transport.gross_model requires positive cp_mass from the thermo model.");
    }

public:
    auto DeviceData() const
    {
        return std::tie(lewis_numbers, prandtl, lambda_a, lambda_b);
    }

    bool UsesCommonDiffusivity() const
    {
        return common_lewis_number;
    }

    template<class Composition>
    Set::Scalar dynamic_viscosity(Set::Scalar T, const Composition& X, int i, int j, int k) const
    {
        return prandtl * thermal_conductivity(T, X, i, j, k) /
            cp_mass_value;
    }

    template<class Composition>
    Set::Scalar thermal_conductivity(Set::Scalar T, const Composition& /*X*/, int /*i*/, int /*j*/, int /*k*/) const
    {
        return lambda_a * T + lambda_b;
    }

    template<class Composition>
    Set::Scalar diffusion_coefficient(Set::Scalar T, Set::Scalar P,
                                    const Composition& X, int i, int j, int k,
                                    int species) const
    {
        const Set::Scalar kappa = thermal_conductivity(T, X, i, j, k);
        Set::Scalar mixture_molecular_weight = 0.0;
        for (int n = 0; n < nspecies; ++n)
            mixture_molecular_weight += X(i,j,k,n) * molecular_weight[n];
        const Set::Scalar rho = (P > 0.0 && T > 0.0) ?
            P * mixture_molecular_weight / (Set::Constant::Rg * T) : 0.0;
        return (rho > 0.0) ?
            kappa / rho / cp_mass_value / lewis_numbers[species] : 0.0;
    }

    template<class Composition>
    void diffusion_coeffs(Set::Patch<Set::Scalar> DKM, Set::Scalar T,
                        Set::Scalar P, const Composition& X,
                        int i, int j, int k) const
    {
        for (int a=0; a<nspecies; ++a)
            DKM(i,j,k,a) = diffusion_coefficient(T, P, X, i, j, k, a);
    }

}; // class GrossModel

} // namespace Transport
} // namespace Gas
} // namespace Model

#endif

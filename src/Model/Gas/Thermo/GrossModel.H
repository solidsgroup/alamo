//
// :code:Model::Gas::Thermo::GrossModel:code: provides a reduced thermodynamic model used
// with GrossModel-style propellant gas closures.
//
// The model assumes a common constant mass heat capacity. Species enthalpy is
// linear in temperature relative to :code:Tref:code:, entropy is zero, and
// mixture molecular weight is evaluated from the configured species values.
//
// Reference: Gross, M. L., Hedman, T. D., Son, S. F., Jackson, T. L., &
// Beckstead, M. W. (2013). Coupling micro and meso-scale combustion models of
// AP/HTPB propellants. Combustion and Flame, 160(5), 982-992.
// https://doi.org/10.1016/j.combustflame.2013.01.016
//

#ifndef MODEL_GAS_THERMO_GROSS_MODEL_H_
#define MODEL_GAS_THERMO_GROSS_MODEL_H_

#include <array>
#include <vector>
#include <memory>
#include <cmath>
#include "IO/ParmParse.H"
#include "Model/Gas/Thermo/Thermo.H"

namespace Model {
namespace Gas {
namespace Thermo {

class GrossModel
{
public:
    static constexpr const char* name = "gross_model";

private:
    static constexpr int max_species = 64;
    int nspecies = 0;
    std::array<Set::Scalar, max_species> molecular_weight{};
    Set::Scalar cp_mass_value;
    Set::Scalar Tref;

public:
    GrossModel() {} ;

    static void Parse(GrossModel & value, IO::ParmParse & pp,
            const std::vector<double>& molecular_weight)
    {
        Util::Assert(INFO, TEST(!molecular_weight.empty()),
            "gas.mw must contain at least one molecular weight.");
        Util::Assert(INFO, TEST(molecular_weight.size() <= max_species),
            "gas.thermo.gross_model supports at most ", max_species, " species.");
        value.nspecies = static_cast<int>(molecular_weight.size());
        // gravimetric specific heat
        pp.query_default("cp_mass", value.cp_mass_value, "0.3_cal/g/K", Unit::SpecificHeatCapacity());
        // reference temperature
        pp.query_default("Tref", value.Tref, "300.0_K", Unit::Temperature());

        for (std::size_t n = 0; n < molecular_weight.size(); ++n)
        {
            Util::Assert(INFO,
                TEST(std::isfinite(molecular_weight[n]) && molecular_weight[n] > 0.0),
                "Invalid molecular weight for species ", n, ": ", molecular_weight[n]);
            value.molecular_weight[n] = molecular_weight[n];
        }

        Util::Assert(INFO, TEST(value.cp_mass_value > 0.0), "gas.thermo.gross_model.cp_mass must be positive.");
        Util::Assert(INFO, TEST(value.Tref >= 0.0), "gas.thermo.gross_model.Tref must be greater than or equal to zero.");
    }

public:
    auto DeviceData() const
    {
        return std::tie(cp_mass_value, Tref);
    }

    template<class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar cp_mol(Set::Scalar /*T*/, const Composition& X,
            int i, int j, int k) const
    {
        return cp_mass_value * GetMW(X, i, j, k);
    }
    template<class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar enthalpy_mol(Set::Scalar T, const Composition& X,
            int i, int j, int k) const
    {
        return cp_mass_value * GetMW(X, i, j, k) * (T - Tref);
    }
    template<class Composition>
    Set::Scalar entropy_mol(Set::Scalar /*T*/, const Composition& /*X*/, int /*i*/, int /*j*/, int /*k*/) const
    {
        return 0.0;
    }
    Set::Scalar cp_mol_species(Set::Scalar /*T*/, int n) const
    {
        return cp_mass_value * molecular_weight[n];
    }
    Set::Scalar dcp_mol_species_dT(Set::Scalar /*T*/, int /*n*/) const
    {
        return 0.0;
    }
    Set::Scalar enthalpy_mol_species(Set::Scalar T, int n) const
    {
        return cp_mass_value * molecular_weight[n] * (T - Tref);
    }
    Set::Scalar entropy_mol_species(Set::Scalar /*T*/, int /*n*/) const
    {
        return 0.0;
    }

    template<class Fraction>
    AMREX_FORCE_INLINE
    Fraction MoleFractions(
            const amrex::Array4<const Set::Scalar>& partial_density,
            int a_nspecies, int i, int j, int k) const
    {
        return Fraction(partial_density, molecular_weight.data(),
                        a_nspecies, i, j, k);
    }

    template<class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar GetMW(const Composition& X,
            int i, int j, int k) const
    {
        Set::Scalar value = 0.0;
        for (int n = 0; n < nspecies; ++n)
            value += X(i,j,k,n) * molecular_weight[n];
        return value;
    }

}; // class GrossModel

} // namespace Thermo
} // namespace Gas
} // namespace Model

#endif

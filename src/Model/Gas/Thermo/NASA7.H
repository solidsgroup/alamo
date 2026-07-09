#ifndef MODEL_GAS_THERMO_NASA7_H_
#define MODEL_GAS_THERMO_NASA7_H_

#include <array>
#include <vector>
#include <cmath>
#include <memory>
#include "Model/Gas/Thermo/Thermo.H"
#include "IO/ParmParse.H"
#include "IO/CanteraYamlParse.H"
#include "Set/Set.H"

namespace Model {
namespace Gas {
namespace Thermo {

class NASA7
{
public:
    static constexpr const char* name = "nasa7";

private:
    std::vector<std::array<Set::Scalar, 3>> Tranges;
    std::vector<std::array<std::array<Set::Scalar, 7>, 2>> coeffs;

public:
    NASA7() {};

    static void Parse(NASA7 & value, IO::ParmParse & pp)
    {
        value.Tranges.resize(NSPECIES);
        value.coeffs.resize(NSPECIES);

        std::string type;
        // declare if parameters are defined in a yaml file [yaml] or the alamo input file [input]
        pp.query_validate("type", type, {"yaml", "input"});

        if (type == "yaml")
        {
#ifdef ALAMO_YAML
            std::string yaml;
            // Cantera yaml file
            pp_query_required("yaml", yaml);
            CanteraYAML::CanteraKinetics kinetics = CanteraYAML::ParseYaml(yaml);

            Util::Assert(INFO, TEST(!kinetics.phase.empty()),
                "NASA7 YAML input must define at least one phase.");
            Util::Assert(INFO, TEST((int)kinetics.phase[0].species.size() == NSPECIES),
                "NASA7 YAML active phase species count (", kinetics.phase[0].species.size(),
                ") does not match NSPECIES (", NSPECIES, ").");

            std::array<int, NSPECIES> yaml_species_idx;
            for (int n = 0; n < NSPECIES; ++n) yaml_species_idx[n] = -1;

            for (int n = 0; n < NSPECIES; ++n)
            {
                const std::string& simulation_species_name = kinetics.phase[0].species[n];
                int match = -1;
                for (int j = 0; j < (int)kinetics.species.size(); ++j)
                {
                    if (kinetics.species[j].name == simulation_species_name)
                    {
                        Util::Assert(INFO, TEST(match < 0),
                            "Species '", simulation_species_name, "' appears multiple times in the YAML species list.");
                        match = j;
                    }
                }
                Util::Assert(INFO, TEST(match >= 0),
                    "Species '", simulation_species_name, "' from the active phase is missing from the YAML species table.");
                yaml_species_idx[n] = match;
            }

            for (int n = 0; n < NSPECIES; ++n)
            {
                const CanteraYAML::Species& species = kinetics.species[yaml_species_idx[n]];
                Util::Assert(INFO, TEST(species.thermoModel == "NASA7"),
                    "Species '", species.name, "' in ", yaml, " does not use the NASA7 thermo model.");
                Util::Assert(INFO, TEST(species.thermoTemp.size() == 3),
                    "Species '", species.name, "' must provide exactly 3 NASA7 temperature ranges.");
                Util::Assert(INFO, TEST(species.thermoData.size() == 2),
                    "Species '", species.name, "' must provide exactly 2 NASA7 coefficient sets.");
                Util::Assert(INFO, TEST(species.thermoData[0].size() == 7),
                    "Species '", species.name, "' low-temperature NASA7 coefficients must have length 7.");
                Util::Assert(INFO, TEST(species.thermoData[1].size() == 7),
                    "Species '", species.name, "' high-temperature NASA7 coefficients must have length 7.");

                for (int i = 0; i < 3; ++i) value.Tranges[n][i] = species.thermoTemp[i];
                for (int i = 0; i < 7; ++i)
                {
                    value.coeffs[n][0][i] = species.thermoData[0][i];
                    value.coeffs[n][1][i] = species.thermoData[1][i];
                }
            }
#else
            Util::Abort(INFO, "Alamo was not built with YAML support, but NASA7.yaml was requested.");
#endif
        }
        else if (type == "input")
        {
            std::vector<Set::Scalar> flat_Tranges;
            std::vector<Set::Scalar> flat_coeffs_lo;
            std::vector<Set::Scalar> flat_coeffs_hi;

            // Array of temperature ranges for NASA7 polynomial curve fits, must be size 3 or 3 x NSPECIES
            // T_lo, T_mid, T_high
            pp.queryarr_required("temperature_ranges", flat_Tranges, Unit::Temperature());
            // NASA7 polynomial coefficients for lower temperature range, T_lo:T_mid
            pp.queryarr_required("coeffs_lo", flat_coeffs_lo);
            // NASA7 polynomial coefficients for upper temperature range, T_mid:T_high
            pp.queryarr_required("coeffs_hi", flat_coeffs_hi);

            Util::Assert(INFO, TEST((int)flat_Tranges.size() == 3 || (int)flat_Tranges.size() == 3 * NSPECIES),
                "NASA7 temperature_ranges must contain either 3 shared values or 3 values per species (expected 3 or ",
                3 * NSPECIES, ", got ", flat_Tranges.size(), ").");
            Util::Assert(INFO, TEST((int)flat_coeffs_lo.size() == 7 * NSPECIES),
                "NASA7 coeffs_lo must contain 7 values per species (expected ",
                7 * NSPECIES, ", got ", flat_coeffs_lo.size(), ").");
            Util::Assert(INFO, TEST((int)flat_coeffs_hi.size() == 7 * NSPECIES),
                "NASA7 coeffs_hi must contain 7 values per species (expected ",
                7 * NSPECIES, ", got ", flat_coeffs_hi.size(), ").");

            for (int n = 0; n < NSPECIES; ++n)
            {
                for (int i = 0; i < 3; ++i)
                {
                    if ((int)flat_Tranges.size() == 3) value.Tranges[n][i] = flat_Tranges[i];
                    else                               value.Tranges[n][i] = flat_Tranges[3 * n + i];
                }
                for (int i = 0; i < 7; ++i)
                {
                    value.coeffs[n][0][i] = flat_coeffs_lo[7 * n + i];
                    value.coeffs[n][1][i] = flat_coeffs_hi[7 * n + i];
                }
            }
        }
        else
        {
            Util::Abort(INFO, "Unrecognized type for thermo model ", type);
        }
    }


private:
    // Compute nondimensional quantities (cp/R, h/(RT), s/R) from NASA7 coefficients
    Set::Scalar cp_R(Set::Scalar T, const std::array<Set::Scalar, 3>& Trange, const std::array<std::array<Set::Scalar, 7>, 2>& coeff) const
    {
        int hilo = (T <= Trange[1]) ? 0 : 1;
        const std::array<Set::Scalar, 7>& a = coeff[hilo];
        return a[0] + a[1]*T + a[2]*T*T + a[3]*T*T*T + a[4]*T*T*T*T;
    }

    Set::Scalar h_RT(Set::Scalar T, const std::array<Set::Scalar, 3>& Trange, const std::array<std::array<Set::Scalar, 7>, 2>& coeff) const
    {
        int hilo = (T <= Trange[1]) ? 0 : 1;
        const std::array<Set::Scalar, 7>& a = coeff[hilo];
        return a[0] + 0.5*a[1]*T + (1.0/3.0)*a[2]*T*T + 0.25*a[3]*T*T*T + 0.2*a[4]*T*T*T*T + a[5]/T;
    }

    Set::Scalar s_R(Set::Scalar T, const std::array<Set::Scalar, 3>& Trange, const std::array<std::array<Set::Scalar, 7>, 2>& coeff) const
    {
        int hilo = (T <= Trange[1]) ? 0 : 1;
        const std::array<Set::Scalar, 7>& a = coeff[hilo];
        return a[0]*log(T) + a[1]*T + 0.5*a[2]*T*T + (1.0/3.0)*a[3]*T*T*T + 0.25*a[4]*T*T*T*T + a[6];
    }

public:
    Set::Scalar cp_mol(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        // Specific heat (constant pressure), J/(kmol-K)
        Set::Scalar cp = 0.0;
        for (int n=0; n<NSPECIES; ++n) cp += X(i,j,k,n) * cp_R(T, Tranges[n], coeffs[n]) * Set::Constant::Rg;
        return cp;
    }
    Set::Scalar enthalpy_mol(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        // Specific enthalpy, J/kmol)
        Set::Scalar h = 0.0;
        for (int n=0; n<NSPECIES; ++n) h += X(i,j,k,n) * h_RT(T, Tranges[n], coeffs[n]) * Set::Constant::Rg * T;
        return h;
    }
    Set::Scalar entropy_mol(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        // specific entropy, J/(kmol-K)
        Set::Scalar s = 0.0;
        for (int n=0; n<NSPECIES; ++n) s += X(i,j,k,n) * s_R(T, Tranges[n], coeffs[n]) * Set::Constant::Rg;
        return s;
    }
    Set::Scalar cp_mol_species(Set::Scalar T, int n) const
    {
        return cp_R(T, Tranges[n], coeffs[n]) * Set::Constant::Rg;
    }
    Set::Scalar enthalpy_mol_species(Set::Scalar T, int n) const
    {
        return h_RT(T, Tranges[n], coeffs[n]) * Set::Constant::Rg * T;
    }
    Set::Scalar entropy_mol_species(Set::Scalar T, int n) const
    {
        return s_R(T, Tranges[n], coeffs[n]) * Set::Constant::Rg;
    }

}; // class NASA7

} // namespace Thermo
} // namespace Gas
} // namespace Model

#endif

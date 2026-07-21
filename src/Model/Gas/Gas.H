//
// This header defines the :code:Model::Gas::Gas:code: class, which provides the
// gas-property interface used by the hydrodynamics and chemistry models.
//
// The :code:Gas:code: class owns species molecular weights and selected thermodynamic,
// transport, and equation-of-state models. These submodels are selected from input
// using static polymorphism, then exposed through a single interface for mixture
// thermodynamic properties, transport coefficients, equation-of-state conversions, and
// local mass- and mole-fraction updates.
//

#ifndef MODEL_GAS_H_
#define MODEL_GAS_H_

#include <array>
#include <vector>
#include <cmath>
#include <memory>
#include <string>
#include "Model/Gas/Thermo/CpConstant.H"
#include "Model/Gas/Thermo/Rocfire.H"
#include "Model/Gas/Thermo/NASA7.H"
#include "Model/Gas/Transport/Mixture_Averaged.H"
#include "Model/Gas/Transport/Rocfire.H"
#include "Model/Gas/EOS/CPG.H"
#include "Model/Gas/EOS/TPG.H"
#include "Model/Gas/EOS/Rocfire.H"
#include "Util/Util.H"
#include "Set/Base.H"
#include "Set/Set.H"
#include "IO/ParmParse.H"

namespace Model {
namespace Gas {

class Gas {
protected:
public:
    Set::Scalar Rg;
    std::array<Set::Scalar, NSPECIES> MW;      // Species molecular weights, kg/kmol
    Thermo::Thermo<Thermo::CpConstant,Thermo::Rocfire,Thermo::NASA7> thermo;
    Transport::Transport<Transport::Mixture_Averaged,Transport::Rocfire> transport;
    EOS::EOS<EOS::CPG,EOS::TPG,EOS::Rocfire> eos;
    bool rocfire = false;

public:
    Gas()
    {
        this->Rg = Set::Constant::Rg;
        this->MW.fill(NAN);
    }
    Gas(IO::ParmParse& pp, std::string name) : Gas()
    {
        pp.queryclass(name, *this);
    }

    static void Parse(Gas & value, IO::ParmParse & pp)
    {
        // number of species
        pp.queryarr_required("mw",value.MW,Unit::Mass()/Unit::Amount());
        Util::Assert(INFO, TEST(value.MW.size() == NSPECIES), "Error: Runtime species count from input file (", value.MW.size(), ") doesn't match the compile-time NSPECIES (", NSPECIES, ")");
        for (int n=0; n<NSPECIES; ++n)
        {
            Util::Assert(INFO, TEST(std::isfinite(value.MW[n]) && value.MW[n] > 0.0),
                "Invalid molecular weight for species ", n, ": ", value.MW[n]);
        }

        // thermo model
        pp.select<Thermo::CpConstant, Thermo::NASA7, Thermo::Rocfire>("thermo", value.thermo);

        // transport model
        pp.select<Transport::Mixture_Averaged, Transport::Rocfire>("transport", value.transport, value.MW, value.thermo);

        // eos model
        pp.select<EOS::CPG, EOS::TPG, EOS::Rocfire>("eos", value.eos, value.thermo);

        value.rocfire = std::string(value.eos.model_name()) == "rocfire";
    }

public:
    // Thermodynamic quantities
    Set::Scalar cp_mol(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    Set::Scalar enthalpy_mol(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    Set::Scalar entropy_mol(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    Set::Scalar cp_mol_species(Set::Scalar T, int n) const;
    Set::Scalar enthalpy_mol_species(Set::Scalar T, int n) const;
    Set::Scalar entropy_mol_species(Set::Scalar T, int n) const;

    // Transport quantities
    Set::Scalar dynamic_viscosity(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    Set::Scalar thermal_conductivity(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    void diffusion_coeffs(Set::Patch<Set::Scalar>& DKM, Set::Scalar T, Set::Scalar P, Set::Patch<const Set::Scalar>& X, int i, int j, int k);

    // EOS
    #if AMREX_SPACEDIM == 2
    Set::Scalar ComputeT(
            Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar E, Set::Scalar Tguess,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k, Set::Scalar rtol=1e-12) const;
    Set::Scalar ComputeE(
            Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar T,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    #elif AMREX_SPACEDIM == 3
    Set::Scalar ComputeT(
            Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz, Set::Scalar E, Set::Scalar Tguess,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k, Set::Scalar rtol=1e-12) const;
    Set::Scalar ComputeE(
            Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz, Set::Scalar T,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    #endif
    Set::Scalar ComputeT_from_primitives(
            Set::Scalar pressure, Set::Scalar density,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    Set::Scalar ComputeP(Set::Scalar density, Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;

// Universal methods inherent to Gas
public:
    Set::Scalar ComputeLocalFractions(
            Set::Patch<const Set::Scalar>& density_mf,
            Set::Patch<Set::Scalar>& mass_fraction_mf,
            Set::Patch<Set::Scalar>& mole_fraction_mf,
            const int i, const int j, const int k)
    {
        // In place update for mass and mole fractions
        // This needs to be called before any other functions if cell composition
        // has changed, since they are dependent on mole fraction
        if (rocfire)
        {
            Set::Scalar raw_density = 0.0;
            Set::Scalar density_scale = 0.0;
            Set::Scalar rho_clipped[NSPECIES];
            for (int n=0; n<NSPECIES; ++n)
            {
                const Set::Scalar rhoYn = density_mf(i,j,k,n);
                if (!std::isfinite(rhoYn))
                {
                    Util::Abort(INFO, "Non-finite species density in Rocfire local fractions at (",
                        i, ",", j, ",", k, "), species ", n, ": ", rhoYn);
                }
                raw_density += rhoYn;
                density_scale += std::fabs(rhoYn);
            }

            const Set::Scalar density_tol = 1.0e-12 * ((density_scale > 1.0) ? density_scale : 1.0);
            Set::Scalar density = 0.0;
            for (int n=0; n<NSPECIES; ++n)
            {
                Set::Scalar rhoYn = density_mf(i,j,k,n);
                if (rhoYn < 0.0)
                {
                    if (rhoYn < -density_tol)
                    {
                        Util::Abort(INFO, "Negative species density in Rocfire local fractions at (",
                            i, ",", j, ",", k, "), species ", n, ": rhoY=", rhoYn,
                            ", raw density=", raw_density, ", tolerance=", density_tol);
                    }
                    rhoYn = 0.0;
                }
                rho_clipped[n] = rhoYn;
                density += rhoYn;
            }
            if (!std::isfinite(density) || density <= 0.0)
            {
                Util::Abort(INFO, "Invalid Rocfire mixture state at (",
                    i, ",", j, ",", k, "): raw density=", raw_density,
                    ", clipped density=", density, ", tolerance=", density_tol);
            }

            const Set::Scalar inv_density = 1.0 / density;
            for (int n=0; n<NSPECIES; ++n)
            {
                const Set::Scalar fraction = rho_clipped[n] * inv_density;
                mass_fraction_mf(i,j,k,n) = fraction;
                mole_fraction_mf(i,j,k,n) = fraction;
            }
            return density;
        }

        Set::Scalar raw_density = 0.0;
        Set::Scalar density_scale = 0.0;
        for (int n=0; n<NSPECIES; ++n)
        {
            Util::Assert(INFO, TEST(std::isfinite(MW[n]) && MW[n] > 0.0),
                "Invalid molecular weight for species ", n, ": ", MW[n]);
            raw_density += density_mf(i,j,k,n);
            density_scale += std::fabs(density_mf(i,j,k,n));
        }

        const Set::Scalar density_tol = 1.0e-12 * ((density_scale > 1.0) ? density_scale : 1.0);
        Set::Scalar density = 0.0;
        Set::Scalar moles = 0.0;
        Set::Scalar rho_clipped[NSPECIES];
        for (int n=0; n<NSPECIES; ++n)
        {
            rho_clipped[n] = density_mf(i,j,k,n);
            if (!std::isfinite(rho_clipped[n]))
            {
                Util::Abort(INFO, "Non-finite species density in ComputeLocalFractions at (",
                    i, ",", j, ",", k, "), species ", n, ": ", rho_clipped[n]);
            }
            if (rho_clipped[n] < 0.0)
            {
                if (rho_clipped[n] < -density_tol)
                {
                    Util::Abort(INFO, "Negative species density in ComputeLocalFractions at (",
                        i, ",", j, ",", k, "), species ", n, ": rhoY=", rho_clipped[n],
                        ", raw density=", raw_density, ", tolerance=", density_tol);
                }
                rho_clipped[n] = 0.0;
            }
            density += rho_clipped[n];
            moles += rho_clipped[n] / MW[n];
        }
        if (!std::isfinite(density) || !std::isfinite(moles) || density <= 0.0 || moles <= 0.0)
        {
            Util::Abort(INFO, "Invalid mixture state in ComputeLocalFractions at (",
                i, ",", j, ",", k, "): raw density=", raw_density,
                ", clipped density=", density, ", clipped moles=", moles,
                ", tolerance=", density_tol);
        }
        Set::Scalar mass_fraction_sum = 0.0;
        Set::Scalar mole_fraction_sum = 0.0;
        for (int n=0; n<NSPECIES; ++n) {
            mass_fraction_mf(i,j,k,n) = rho_clipped[n] / density;
            mole_fraction_mf(i,j,k,n) = rho_clipped[n] / MW[n] / moles;
            mass_fraction_sum += mass_fraction_mf(i,j,k,n);
            mole_fraction_sum += mole_fraction_mf(i,j,k,n);
            //Util::Message(INFO, "Mass fraction[n]: ", n, " ", mass_fraction_mf(i,j,k,n));
            Util::Assert(INFO, TEST(mass_fraction_mf(i,j,k,n) >= -1e-12),
                "Error: Mass fraction less than 0.0 at (", i, ",", j, ",", k,
                "), species ", n, ": ", mass_fraction_mf(i,j,k,n));
            Util::Assert(INFO, TEST(mass_fraction_mf(i,j,k,n) <= 1.0 + 1e-12),
                "Error: Mass fraction greater than 1.0 at (", i, ",", j, ",", k,
                "), species ", n, ": ", mass_fraction_mf(i,j,k,n));
            Util::Assert(INFO, TEST(mole_fraction_mf(i,j,k,n) >= -1e-12),
                "Error: Mole fraction less than 0.0 at (", i, ",", j, ",", k,
                "), species ", n, ": ", mole_fraction_mf(i,j,k,n));
            Util::Assert(INFO, TEST(mole_fraction_mf(i,j,k,n) <= 1.0 + 1e-12),
                "Error: Mole fraction greater than 1.0 at (", i, ",", j, ",", k,
                "), species ", n, ": ", mole_fraction_mf(i,j,k,n));
        }
        Util::Assert(INFO, TEST(std::fabs(mass_fraction_sum - 1.0) <= 1e-12),
            "Error: Mass fractions do not sum to 1.0 at (", i, ",", j, ",", k,
            "): sum=", mass_fraction_sum);
        Util::Assert(INFO, TEST(std::fabs(mole_fraction_sum - 1.0) <= 1e-12),
            "Error: Mole fractions do not sum to 1.0 at (", i, ",", j, ",", k,
            "): sum=", mole_fraction_sum);

        return raw_density;
    }

    Set::Scalar GetMW(Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Molecular weight, kg/kmol
        if (rocfire) return MW[0];
        Set::Scalar mw = 0.0;
        for (int n=0; n<NSPECIES; ++n) mw += X(i,j,k,n) * MW[n];
        return mw;
    }
    Set::Scalar ComputeD(Set::Patch<const Set::Scalar>& rhoY, int i, int j, int k) const {
        // Mixture density kg/m^3
        Set::Scalar rho = 0.0;
        for (int n=0; n<NSPECIES; ++n) rho += rhoY(i,j,k,n);
        return rho;
    }
    Set::Scalar R(Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific gas constant, J/(kg-K)
        return Rg/GetMW(X, i, j, k);
    }
    Set::Scalar cp_mass(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific heat (constant pressure), J/(kg-K)
        if (rocfire) return cp_mol_species(T, 0) / MW[0];
        return cp_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }
    Set::Scalar cv_mol(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific heat (constant volume), J/(kmol-K)
        return cp_mol(T, X, i, j, k) - Rg;
    }
    Set::Scalar cv_mass(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific heat (constant volume), J/(kg-K)
        return cv_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }
    Set::Scalar gamma(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const  {
        // Specific heat ratio
        if (rocfire)
        {
            const Set::Scalar cp_mol = cp_mol_species(T, 0);
            return cp_mol / (cp_mol - Rg);
        }
        return cp_mol(T, X, i, j, k) / cv_mol(T, X, i, j, k);
    }
    Set::Scalar enthalpy_mass(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific enthalpy, J/kg
        return enthalpy_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }
    Set::Scalar entropy_mass(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific entropy, J/(kg-K)
        return entropy_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }

}; // class Gas

} // namespace Gas
} // namespace Model

#endif

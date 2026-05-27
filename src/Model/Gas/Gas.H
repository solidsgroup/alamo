#ifndef MODEL_GAS_H_
#define MODEL_GAS_H_

#include <vector>
#include <cmath>
#include <memory>
#include "Model/Gas/Thermo/CpConstant.H"
#include "Model/Gas/Thermo/Rocfire.H"
#include "Model/Gas/Thermo/NASA7.H"
#include "Model/Gas/Transport/Mixture_Averaged.H"
#include "Model/Gas/Transport/Rocfire.H"
#include "Model/Gas/EOS/CPG.H"
#include "Model/Gas/EOS/TPG.H"
#include "Util/Util.H"
#include "Set/Base.H"
#include "Set/Set.H"
#include "IO/ParmParse.H"

namespace Model {
namespace Gas {
class Gas {
protected:
public:
    std::vector<double> MW;      // Species molecular weights, kg/kmol
    Thermo::Thermo *thermo = nullptr;
    Transport::Transport *transport = nullptr;
    EOS::EOS *eos = nullptr;

public:
    virtual ~Gas() {
        delete thermo;
        delete transport;
        delete eos;
    }

    static void Parse(Gas & value, IO::ParmParse & pp)
    {
        // number of species
        pp.queryarr_required("mw",value.MW,Unit::Mass()/Unit::Amount());
        Util::Assert(INFO, TEST(value.MW.size() == NSPECIES), "Error: Runtime species count from input file (", value.MW.size(), ") doesn't match the compile-time NSPECIES (", NSPECIES, ")");

        delete value.thermo;
        value.thermo = nullptr;
        // thermo model
        pp.select<Thermo::CpConstant, Thermo::NASA7, Thermo::Rocfire>("thermo", value.thermo);

        delete value.transport;
        value.transport = nullptr;
        // transport model
        pp.select<Transport::Mixture_Averaged, Transport::Rocfire>("transport", value.transport, value.MW, value.thermo);

        delete value.eos;
        value.eos = nullptr;
        // eos model
        pp.select<EOS::CPG, EOS::TPG>("eos", value.eos, &value);
    }

public:
    // Thermodynamic quantities
    double cp_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    double enthalpy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    double entropy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    double cp_mol_species(double T, int n) const;
    double enthalpy_mol_species(double T, int n) const;
    double entropy_mol_species(double T, int n) const;

    // Transport quantities
    double dynamic_viscosity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    double thermal_conductivity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    void diffusion_coeffs(Set::Patch<Set::Scalar>& DKM, double T, double P, Set::Patch<const Set::Scalar>& X, int i, int j, int k);

    // EOS
    #if AMREX_SPACEDIM == 2
    double ComputeT(
            double density, double momentumx, double momentumy, double E, double Tguess,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol=1e-12) const;
    double ComputeE(    
            double density, double momentumx, double momentumy, double T,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    #elif AMREX_SPACEDIM == 3
    double ComputeT(
            double density, double momentumx, double momentumy, double momentumz, double E, double Tguess,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol=1e-12) const;
    double ComputeE(    
            double density, double momentumx, double momentumy, double momentumz, double T,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    #endif
    double ComputeT_from_primitives(
            double pressure, double density,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    double ComputeP(double density, double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;

// Universal methods inherent to Gas
public:
    void ComputeLocalFractions(
            Set::Patch<const Set::Scalar>& density_mf,
            Set::Patch<Set::Scalar>& mass_fraction_mf,
            Set::Patch<Set::Scalar>& mole_fraction_mf,
            const int i, const int j, const int k)
    {
        // In place update for mass and mole fractions
        // This needs to be called before any other functions if cell composition 
        // has changed, since they are dependent on mole fraction
        double raw_density = 0.0;
        double density_scale = 0.0;
        for (int n=0; n<NSPECIES; ++n)
        {
            Util::Assert(INFO, TEST(std::isfinite(MW[n]) && MW[n] > 0.0),
                "Invalid molecular weight for species ", n, ": ", MW[n]);
            raw_density += density_mf(i,j,k,n);
            density_scale += std::fabs(density_mf(i,j,k,n));
        }

        const double density_tol = 1.0e-12 * ((density_scale > 1.0) ? density_scale : 1.0);
        double density = 0.0;
        double moles = 0.0;
        double rho_clipped[NSPECIES];
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
        double mass_fraction_sum = 0.0;
        double mole_fraction_sum = 0.0;
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
    }
    double GetMW(Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Molecular weight, kg/kmol
        double mw = 0.0;
        for (int n=0; n<NSPECIES; ++n) mw += X(i,j,k,n) * MW[n];
        return mw;
    }
    double ComputeD(Set::Patch<const Set::Scalar>& rhoY, int i, int j, int k) const {
        // Mixture density kg/m^3
        double rho = 0.0;
        for (int n=0; n<NSPECIES; ++n) rho += rhoY(i,j,k,n);
        return rho;
    }
    double R(Set::Patch<const Set::Scalar>& X, int i, int j, int k) const { 
        // Specific gas constant, J/(kg-K)
        return Set::Constant::Rg/GetMW(X, i, j, k);
    }
    double cp_mass(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific heat (constant pressure), J/(kg-K)
        return cp_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }
    double cv_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific heat (constant volume), J/(kmol-K)
        return cp_mol(T, X, i, j, k) - Set::Constant::Rg;
    }
    double cv_mass(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific heat (constant volume), J/(kg-K)
        return cv_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }
    double gamma(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const  {
        // Specific heat ratio
        return cp_mol(T, X, i, j, k) / cv_mol(T, X, i, j, k);
    }
    double enthalpy_mass(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific enthalpy, J/kg
        return enthalpy_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }
    double entropy_mass(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific entropy, J/(kg-K)
        return entropy_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }

}; // class Gas

} // namespace Gas
} // namespace Model

#endif

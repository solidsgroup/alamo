#ifndef MODEL_GAS_H_
#define MODEL_GAS_H_

#include <vector>
#include <cmath>
#include <memory>
#include "Model/Gas/Thermo/CpConstant.H"
#include "Model/Gas/Transport/Mixture_Averaged.H"
#include "Model/Gas/EOS/CPG.H"
#include "Util/Util.H"
#include "Set/Base.H"
#include "Set/Set.H"
#include "IO/ParmParse.H"

namespace Model {
namespace Gas {

class Gas {
protected:
public:
    int nspecies;
    double Rg;
    std::vector<double> MW;      // Species molecular weights, kg/kmol
    Thermo::Thermo<Thermo::CpConstant> thermo;
    Transport::Transport<Transport::Mixture_Averaged> transport;
    // [replace!] EOS::EOS *eos = nullptr;

public:
    Gas() {} ;
    Gas(IO::ParmParse& pp, std::string name) : Gas()
    {
        pp.queryclass(name, *this);
    }

    static void Parse(Gas & value, IO::ParmParse & pp)
    {
        // number of species
        pp.queryarr_required("mw",value.MW,Unit::Mass()/Unit::Amount());
        value.nspecies = value.MW.size();
        value.Rg = Set::Constant::Rg;

        // thermo model
        // [replace!] pp.select<Thermo::CpConstant>("thermo", value.thermo, value.nspecies);

        // transport model
        // [replace!] pp.select<Transport::Mixture_Averaged>("transport", value.transport, value.nspecies, value.MW, value.thermo);

        // eos model
        // [replace!] pp.select<EOS::CPG>("eos", value.eos, &value);
    }
    
public:
    // Thermodynamic quantities
    AMREX_GPU_HOST_DEVICE
    double cp_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    AMREX_GPU_HOST_DEVICE
    double enthalpy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    AMREX_GPU_HOST_DEVICE
    double entropy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;

    // Transport quantities
    AMREX_GPU_HOST_DEVICE
    double dynamic_viscosity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    AMREX_GPU_HOST_DEVICE
    double thermal_conductivity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    AMREX_GPU_HOST_DEVICE
    void diffusion_coeffs(Set::Patch<Set::Scalar>& DKM, double T, double P, Set::Patch<const Set::Scalar>& X, int i, int j, int k);

    // EOS
    AMREX_GPU_HOST_DEVICE
    double ComputeT(
            double density, double momentumx, double momentumy, double E, double Tguess,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol=1e-12) const;
    AMREX_GPU_HOST_DEVICE
    double ComputeT(
            double pressure, double density,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    AMREX_GPU_HOST_DEVICE
    double ComputeP(double density, double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    AMREX_GPU_HOST_DEVICE
    double ComputeE(    
            double density, double momentumx, double momentumy, double T,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;

// Universal methods inherent to Gas
public:
    AMREX_GPU_HOST_DEVICE
    void ComputeLocalFractions(
            Set::Patch<const Set::Scalar>& density_mf,
            Set::Patch<Set::Scalar>& mass_fraction_mf,
            Set::Patch<Set::Scalar>& mole_fraction_mf,
            const int i, const int j, const int k) const
    {
        // In place update for mass and mole fractions
        // This needs to be called before any other functions if cell composition 
        // has changed, since they are dependent on mole fraction
        double density = 0.0;
        double moles = 0.0;
        for (int n=0; n<nspecies; ++n) {
            density += density_mf(i,j,k,n);
            moles += density_mf(i,j,k,n) / MW[n];
        }
        for (int n=0; n<nspecies; ++n) {
            mass_fraction_mf(i,j,k,n) = density_mf(i,j,k,n) / density;
            mole_fraction_mf(i,j,k,n) = density_mf(i,j,k,n) / MW[n] / moles;
        }
    }
    AMREX_GPU_HOST_DEVICE
    double GetMW(Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Molecular weight, kg/kmol
        double mw = 0.0;
        for (int n=0; n<nspecies; ++n) mw += X(i,j,k,n) * MW[n];
        return mw;
    }
    AMREX_GPU_HOST_DEVICE
    double ComputeD(Set::Patch<const Set::Scalar>& rhoY, int i, int j, int k) const {
        // Mixture density kg/m^3
        double rho = 0.0;
        for (int n=0; n<nspecies; ++n) rho += rhoY(i,j,k,n);
        return rho;
    }
    AMREX_GPU_HOST_DEVICE
    double R(Set::Patch<const Set::Scalar>& X, int i, int j, int k) const { 
        // Specific gas constant, J/(kg-K)
        return Rg/GetMW(X, i, j, k);
    }
    AMREX_GPU_HOST_DEVICE
    double cp_mass(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific heat (constant pressure), J/(kg-K)
        return cp_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }
    AMREX_GPU_HOST_DEVICE
    double cv_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific heat (constant volume), J/(kmol-K)
        return cp_mol(T, X, i, j, k) - Rg;
    }
    AMREX_GPU_HOST_DEVICE
    double cv_mass(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific heat (constant volume), J/(kg-K)
        return cv_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }
    AMREX_GPU_HOST_DEVICE
    double gamma(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const  {
        // Specific heat ratio
        return cp_mol(T, X, i, j, k) / cv_mol(T, X, i, j, k);
    }
    AMREX_GPU_HOST_DEVICE
    double enthalpy_mass(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific enthalpy, J/kg
        return enthalpy_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }
    AMREX_GPU_HOST_DEVICE
    double entropy_mass(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific entropy, J/(kg-K)
        return entropy_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }

}; // class Gas

} // namespace Gas
} // namespace Model

#endif

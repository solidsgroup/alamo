#ifndef MODEL_GAS_H_
#define MODEL_GAS_H_

#include <vector>
#include <cmath>
#include <memory>
#include "Model/Gas/Thermo/CpConstant.H"
#include "Model/Gas/Thermo/NASA7.H"
#include "Model/Gas/Thermo/Rocfire.H"
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

/// Cell-local mole fractions reconstructed lazily from gas partial densities.
class MoleFraction
{
public:
    MoleFraction(const amrex::Array4<const Set::Scalar>& a_partial_density,
                const double* a_molecular_weight,
                const int a_nspecies,
                const int i, const int j, const int k)
        : partial_density(a_partial_density),
        molecular_weight(a_molecular_weight),
        moles(0.0)
    {
        for (int n = 0; n < a_nspecies; ++n)
            moles += partial_density(i,j,k,n) / molecular_weight[n];
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar operator()(const int i, const int j, const int k, const int n) const
    {
        return moles > 0.0 ?
            (partial_density(i,j,k,n) / molecular_weight[n]) / moles :
            (n == 0 ? 1.0 : 0.0);
    }

private:
    amrex::Array4<const Set::Scalar> partial_density;
    const double* molecular_weight;
    Set::Scalar moles;
};

class Gas {
protected:
public:
    int nspecies;
    double Rg;
    std::vector<double> MW;      // Species molecular weights, kg/kmol
    Thermo::Thermo<Thermo::CpConstant,Thermo::NASA7,Thermo::Rocfire> thermo;
    Transport::Transport<Transport::Mixture_Averaged,Transport::Rocfire> transport;
    EOS::EOS<EOS::CPG,EOS::TPG> eos;

public:
    Gas() 
    {
        this->Rg = Set::Constant::Rg;
    }
    Gas(IO::ParmParse& pp, std::string name) : Gas()
    {
        pp.queryclass(name, *this);
    }

    static void Parse(Gas & value, IO::ParmParse & pp)
    {
        // number of species
        pp.queryarr_required("mw",value.MW,Unit::Mass()/Unit::Amount());
        value.nspecies = value.MW.size();

        // thermo model
        pp.select<Thermo::CpConstant,Thermo::NASA7,Thermo::Rocfire>("thermo", value.thermo, pp.forward_args(value.nspecies));

        // transport model
        pp.select<Transport::Mixture_Averaged,Transport::Rocfire>("transport", value.transport, pp.forward_args(value.nspecies, value.MW, value.thermo));

        // eos model
        pp.select<EOS::CPG,EOS::TPG>("eos", value.eos);
    }
    
public:
    // Thermodynamic quantities
    double cp_mol(double T, Set::Patch<const Set::Scalar> X, int i, int j, int k) const;
    double enthalpy_mol(double T, Set::Patch<const Set::Scalar> X, int i, int j, int k) const;
    double entropy_mol(double T, Set::Patch<const Set::Scalar> X, int i, int j, int k) const;
    double cp_mol_species(double T, int n) const;
    double dcp_mol_species_dT(double T, int n) const;
    double enthalpy_mol_species(double T, int n) const;
    double entropy_mol_species(double T, int n) const;

    // Transport quantities
    double dynamic_viscosity(double T, Set::Patch<const Set::Scalar> X, int i, int j, int k) const;
    double thermal_conductivity(double T, Set::Patch<const Set::Scalar> X, int i, int j, int k) const;
    void diffusion_coeffs(Set::Patch<Set::Scalar> DKM, double T, double P, Set::Patch<const Set::Scalar> X, int i, int j, int k);

    // EOS
    double ComputeT(
            double density, double momentumx, double momentumy, double E, double Tguess,
            Set::Patch<const Set::Scalar> X, int i, int j, int k, double rtol=1e-12) const;
    double ComputeT(
            double pressure, double density,
            Set::Patch<const Set::Scalar> X, int i, int j, int k) const;
    double ComputeP(double density, double T, Set::Patch<const Set::Scalar> X, int i, int j, int k) const;
    double ComputeE(    
            double density, double momentumx, double momentumy, double T,
            Set::Patch<const Set::Scalar> X, int i, int j, int k) const;

    MoleFraction MoleFractions(const amrex::Array4<const Set::Scalar>& partial_density,
                                const int i, const int j, const int k) const
    {
        return MoleFraction(partial_density, MW.data(), nspecies, i, j, k);
    }

    template<class Composition>
    AMREX_FORCE_INLINE
    double cp_mol(double T, const Composition& X, int i, int j, int k) const
    {
        return thermo.cp_mol(T, X, i, j, k);
    }

    template<class Composition>
    AMREX_FORCE_INLINE
    double enthalpy_mol(double T, const Composition& X, int i, int j, int k) const
    {
        return thermo.enthalpy_mol(T, X, i, j, k);
    }

    template<class Composition>
    AMREX_FORCE_INLINE
    double entropy_mol(double T, const Composition& X, int i, int j, int k) const
    {
        return thermo.entropy_mol(T, X, i, j, k);
    }

    template<class Composition>
    AMREX_FORCE_INLINE
    double dynamic_viscosity(double T, const Composition& X, int i, int j, int k) const
    {
        return transport.dynamic_viscosity(T, X, i, j, k);
    }

    template<class Composition>
    AMREX_FORCE_INLINE
    double thermal_conductivity(double T, const Composition& X, int i, int j, int k) const
    {
        return transport.thermal_conductivity(T, X, i, j, k);
    }

    template<class Composition>
    AMREX_FORCE_INLINE
    double diffusion_coefficient(double T, double P, const Composition& X,
                                int i, int j, int k, int species) const
    {
        return transport.diffusion_coefficient(T, P, X, i, j, k, species);
    }

    template<class Composition>
    AMREX_FORCE_INLINE
    double ComputeE(double density, double momentumx, double momentumy, double T,
                    const Composition& X, int i, int j, int k) const
    {
        return eos.ComputeE(*this, density, momentumx, momentumy, T, X, i, j, k);
    }

// Universal methods inherent to Gas
public:
    void ComputeLocalFractions(
            Set::Patch<const Set::Scalar> density_mf,
            Set::Patch<Set::Scalar> mass_fraction_mf,
            Set::Patch<Set::Scalar> mole_fraction_mf,
            const int i, const int j, const int k)
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
    template<class Composition>
    double GetMW(const Composition& X, int i, int j, int k) const {
        // Molecular weight, kg/kmol
        double mw = 0.0;
        for (int n=0; n<nspecies; ++n) mw += X(i,j,k,n) * MW[n];
        return mw;
    }
    double ComputeD(Set::Patch<const Set::Scalar> rhoY, int i, int j, int k) const {
        // Mixture density kg/m^3
        double rho = 0.0;
        for (int n=0; n<nspecies; ++n) rho += rhoY(i,j,k,n);
        return rho;
    }
    template<class Composition>
    double R(const Composition& X, int i, int j, int k) const {
        // Specific gas constant, J/(kg-K)
        return Rg/GetMW(X, i, j, k);
    }
    template<class Composition>
    double cp_mass(double T, const Composition& X, int i, int j, int k) const {
        // Specific heat (constant pressure), J/(kg-K)
        return cp_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }
    template<class Composition>
    double cv_mol(double T, const Composition& X, int i, int j, int k) const {
        // Specific heat (constant volume), J/(kmol-K)
        return cp_mol(T, X, i, j, k) - Rg;
    }
    template<class Composition>
    double cv_mass(double T, const Composition& X, int i, int j, int k) const {
        // Specific heat (constant volume), J/(kg-K)
        return cv_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }
    template<class Composition>
    double gamma(double T, const Composition& X, int i, int j, int k) const  {
        // Specific heat ratio
        return cp_mol(T, X, i, j, k) / cv_mol(T, X, i, j, k);
    }
    double enthalpy_mass(double T, Set::Patch<const Set::Scalar> X, int i, int j, int k) const {
        // Specific enthalpy, J/kg
        return enthalpy_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }
    template<class Composition>
    double enthalpy_mass(double T, const Composition& X, int i, int j, int k) const {
        return enthalpy_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }
    double entropy_mass(double T, Set::Patch<const Set::Scalar> X, int i, int j, int k) const {
        // Specific entropy, J/(kg-K)
        return entropy_mol(T, X, i, j, k) / GetMW(X, i, j, k);
    }

}; // class Gas

} // namespace Gas
} // namespace Model

#endif

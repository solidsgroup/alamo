#ifndef MODEL_GAS_H_
#define MODEL_GAS_H_

#include <vector>
#include <cmath>
#include <memory>
#include "Util/Util.H"
#include "Set/Base.H"
#include "Set/Set.H"

namespace Model {
namespace Gas {

// Forward declarations
namespace Thermo { class Thermo; }
namespace Transport { class Transport; }
namespace EOS { class EOS; }

struct GasData {
    int nspecies;                // Number of species
    std::vector<double> MW;      // Species molecular weights, kg/kmol
    static constexpr double Ru = 8314.46261815324;     // Universal gas constant, J/(kmol-K)
    static constexpr double Boltzmann = 1.380649e-23;  // Boltzmann constant, J/K
    static constexpr double Avogadro = 6.02214076e26;  // Avogradro's number, particles/kmol
};

class Gas {
protected:
    std::shared_ptr<GasData> data;
    std::shared_ptr<Thermo::Thermo> thermo;
    std::shared_ptr<Transport::Transport> transport;
    std::shared_ptr<EOS::EOS> eos;

public:
    Gas(int nspecies_in, const std::vector<double>& MW_in) {
        data = std::make_shared<GasData>();
        data->nspecies = nspecies_in;
        data->MW = MW_in;
    }
    virtual ~Gas() = default;

    // Attach subsystems
    void SetThermo(std::shared_ptr<Thermo::Thermo> t) { thermo = std::move(t); }
    void SetTransport(std::shared_ptr<Transport::Transport> tr) { transport = std::move(tr); }
    void SetEOS(std::shared_ptr<EOS::EOS> e) { eos = std::move(e); }

    const std::shared_ptr<GasData>& dataPtr() const { return data; }
    //const GasData& data() const { return *data; }
    
// Methods that need to be defined by inherited class
public:
    // Thermodynamic quantities
    double cp_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    double enthalpy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    double entropy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;

    // Transport quantities
    double dynamic_viscosity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    double thermal_conductivity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    void diffusion_coeffs(Set::Patch<Set::Scalar>& DKM, double T, double P, Set::Patch<const Set::Scalar>& X, int i, int j, int k);

    // EOS
    double ComputeT(
            double density, double momentumx, double momentumy, double E, double Tguess,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol=1e-12) const;
    double ComputeT(
            double pressure, double density,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    double ComputeP(double density, double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;
    double ComputeE(    
            double density, double momentumx, double momentumy, double T,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k) const;

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
        double density = 0.0;
        double moles = 0.0;
        for (int n=0; n<data->nspecies; ++n) {
            density += density_mf(i,j,k,n);
            moles += density_mf(i,j,k,n) / data->MW[n];
        }
        for (int n=0; n<data->nspecies; ++n) {
            mass_fraction_mf(i,j,k,n) = density_mf(i,j,k,n) / density;
            mole_fraction_mf(i,j,k,n) = density_mf(i,j,k,n) / data->MW[n] / moles;
        }
    }
    double MW(Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Molecular weight, kg/kmol
        double mw = 0.0;
        for (int n=0; n<data->nspecies; ++n) mw += X(i,j,k,n) * data->MW[n];
        return mw;
    }
    double ComputeD(Set::Patch<const Set::Scalar>& rhoY, int i, int j, int k) const {
        // Mixture density kg/m^3
        double rho = 0.0;
        for (int n=0; n<data->nspecies; ++n) rho += rhoY(i,j,k,n);
        return rho;
    }
    double R(Set::Patch<const Set::Scalar>& X, int i, int j, int k) const { 
        // Specific gas constant, J/(kg-K)
        return GasData::Ru/MW(X, i, j, k);
    }
    double cp_mass(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific heat (constant pressure), J/(kg-K)
        return cp_mol(T, X, i, j, k) / MW(X, i, j, k);
    }
    double cv_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific heat (constant volume), J/(kmol-K)
        return cp_mol(T, X, i, j, k) - GasData::Ru;
    }
    double cv_mass(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific heat (constant volume), J/(kg-K)
        return cv_mol(T, X, i, j, k) / MW(X, i, j, k);
    }
    double gamma(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const  {
        // Specific heat ratio
        return cp_mol(T, X, i, j, k) / cv_mol(T, X, i, j, k);
    }
    double enthalpy_mass(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific enthalpy, J/kg
        return enthalpy_mol(T, X, i, j, k) / MW(X, i, j, k);
    }
    double entropy_mass(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific entropy, J/(kg-K)
        return entropy_mol(T, X, i, j, k) / MW(X, i, j, k);
    }

}; // class Gas

} // namespace Gas
} // namespace Model

#endif

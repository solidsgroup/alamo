#ifndef MODEL_GAS_H_
#define MODEL_GAS_H_

#include <vector>
#include <cmath>
#include <memory>
#include "Util/Util.H"

namespace Model {
namespace Gas {

// Forward declarations
namespace Thermo { class Thermo; }
namespace Transport { class Transport; }
namespace EOS { class EOS; }

struct GasData {
    int nspecies;                               // Number of species
    std::vector<double> MW;                     // Species molecular weights, kg/kmol
    static constexpr Ru = 8314.46261815324;     // Universal gas constant, J/(kmol-K)
    static constexpr Boltzmann = 1.380649e-23;  // Boltzmann constant, J/K
    static constexpr Avogadro = 6.02214076e26;  // Avogradro's number, particles/kmol
};

class Gas {
protected:
    std::shared_ptr<GasData> data_;
    std::unique_ptr<Thermo::Thermo> thermo_;
    std::unique_ptr<Transport::Transport> transport_;
    std::unique_ptr<EOS::EOS> eos_;

public:
    Gas(int nspecies, std::vector<double>& MW);
    virtual ~Gas() = default;

    // Attach subsystems
    void SetThermo(std::unique_ptr<Thermo::Thermo> t) { thermo_ = std::move(t); }
    void SetTransport(std::unique_ptr<Transport::Transport> tr) { transport_ = std::move(tr); }
    void SetEOS(std::unique_ptr<EOS::EOS> e) { eos_ = std::move(e); }

    const std::shared_ptr<GasData>& dataPtr() const { return data_; }
    const GasData& data() const { return *data_; }
    
// Methods that need to be defined by inherited class
public:
    // Thermodynamic quantities
    double cp_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific heat (constant pressure), J/(kmol-K)
        if (!thermo_) Util::Abort(INFO, "[Gas::cp_mol] No Thermo model attached.");
        return thermo_-> cp_mol(T, X, i , j, k);
    }
    double enthalpy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific enthalpy, J/kmol)
        if (!thermo_) Util::Abort(INFO, "[Gas::enthalpy_mol] No Thermo model attached.");
        return thermo_-> enthalpy_mol(T, X, i , j, k);
    }
    double entropy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // specific entropy, J/(kmol-K)
        if (!thermo_) Util::Abort(INFO, "[Gas::entropy_mol] No Thermo model attached.");
        return thermo_-> entropy_mol(T, X, i , j, k);
    }

    // Transport quantities
    double dynamic_viscosity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Dynamic viscosity, Pa-s
        if (!transport_) Util::Abort(INFO, "[Gas::dynamic_viscosity] No Transport model attached.");
        return transport_-> dynamic_viscosity(T, X, i , j, k);
    }
    double thermal_conductivity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Thermal conductivity coefficient, W/(m-K)
        if (!transport_) Util::Abort(INFO, "[Gas::thermal_conductivity] No Transport model attached.");
        return transport_-> thermal_conductivity(T, X, i , j, k);
    }
    void diffusion_coeffs(Set::Patch<Set::Scalar>& DKM, double T, double P, Set::Patch<const Set::Scalar>& X, int i, int j, int k) {
        // Species diffusion coefficients, m^2/s
        if (!transport_) Util::Abort(INFO, "[Gas::diffusion_coeffs] No Transport model attached.");
        return transport_-> diffusion_coeffs(DKM, T, P, X, i , j, k);
    }

    // EOS
    double ComputeT(
            double density, double momentumx, double momentumy, double E, double Tguess,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol) const {
        // Temperature, K
        if (!eos_) Util::Abort(INFO, "[Gas::ComputeT] No EOS model attached.");
        return eos_-> ComputeT(density, momentumx, momentumy, E, Tguess, X, i , j, k, rtol);
    }
    double ComputeP(double density, double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Pressure, Pa
        if (!eos_) Util::Abort(INFO, "[Gas::ComputeP] No EOS model attached.");
        return eos_-> ComputeP(density, T, X, i , j, k);
    }
    double ComputeE(    
            double density, double momentumx, double momentumy, double T,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Energy, J/m^3
        if (!eos_) Util::Abort(INFO, "[Gas::ComputeE] No EOS model attached.");
        return eos_-> ComputeE(density, momentumx, momentumy, T, X, i , j, k);
    }

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
        for (int n=0; n<data_->nspecies; ++n) {
            density += density_mf(i,j,k,n);
            moles += density_mf(i,j,k,n) / data_->MW[n];
        }
        for (int n=0; n<data_->nspecies; ++n) {
            mass_fraction_mf(i,j,k,n) = density_mf(i,j,k,n) / density;
            mole_fraction_mf(i,j,k,n) = density_mf(i,j,k,n / data_->MW[n] / moles;
        }
    }
    double MW(Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Molecular weight, kg/kmol
        double mw = 0.0;
        for (int n=0; n<data_->nspecies; ++n) mw += X(i,j,k,n) * data_->MW[n];
        return mw;
    }
    double ComputeD(Set::Patch<const Set::Scalar>& rhoY, int i, int j, int k) const {
        // Mixture density kg/m^3
        double rho = 0.0;
        for (int n=0; n<data_->nspecies; ++n) rho += rhoY(i,j,k,n);
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


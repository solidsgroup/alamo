#ifndef MODEL_GAS_H_
#define MODEL_GAS_H_

namespace Model {
namespace Gas {

class Gas {
public:
    const int nspecies;
    std::vector<double> _MW;
    Gas(int nsp, std::vector<double>& nMW) : nspecies(nsp), _MW(nMW) {};
    virtual ~Gas() = default;

    static constexpr Ru = 8314.46261815324;     // Universal gas constant, J/(kmol-K)
    static constexpr Boltzmann = 1.380649e-23;  // Boltzmann constant, J/K
    static constexpr Avogadro = 6.02214076e26;  // Avogradro's number, particles/kmol
    
// Methods that need to be defined by inherited class
public:
    // Thermodynamic quantities
    virtual double cp_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;               // Specific heat (constant pressure), J/(kmol-K)
    virtual double enthalpy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;         // Specific enthalpy, J/kmol)
    virtual double entropy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;          // specific entropy, J/(kmol-K)
    // Transport quantities
    virtual double dynamic_viscosity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;    // Dynamic viscosity, Pa-s
    virtual double thermal_conductivity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0; // Thermal conductivity coefficient, W/(m-K)
    virtual void diffusion_coeffs(                                                                                  // Species diffusion coefficients, m^2/s
        Set::Patch<Set::Scalar>& DKM, double T, double P,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) = 0;
    // EOS
    virtual double ComputeT(    // Temperature, K
        double density, double momentumx, double momentumy, double E, double Tguess,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol) const = 0;
    virtual double ComputeP(    // Pressure, Pa
        double density, double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;
    virtual double ComputeE(    // Energy, J/m^3
        double density, double momentumx, double momentumy, double T,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;

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
        for (int n=0; n<nspecies; ++n) {
            density += density_mf(i,j,k,n);
            moles += density_mf(i,j,k,n) / _MW[n];
        }
        for (int n=0; n<nspecies; ++n) {
            mass_fraction_mf(i,j,k,n) = density_mf(i,j,k,n) / density;
            mole_fraction_mf(i,j,k,n) = density_mf(i,j,k,n / _MW[n] / moles;
        }
    }
    double MW(Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Molecular weight, kg/kmol
        double mw = 0.0;
        for (int n=0; n<nspecies; ++n) mw += X(i,j,k,n) * _MW[n];
        return mw;
    }
    double ComputeD(Set::Patch<const Set::Scalar>& rhoY, int i, int j, int k) const {
        // Mixture density kg/m^3
        double rho = 0.0;
        for (int n=0; n<nspecies; ++n) rho += rhoY(i,j,k,n);
        return rho;
    }
    double R(Set::Patch<const Set::Scalar>& X, int i, int j, int k) const { 
        // Specific gas constant, J/(kg-K)
        return Ru/MW(X, i, j, k);
    }
    double cp_mass(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific heat (constant pressure), J/(kg-K)
        return cp_mol(T, X, i, j, k) / MW(X, i, j, k);
    }
    double cv_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
        // Specific heat (constant volume), J/(kmol-K)
        return cp_mol(T, X, i, j, k) - Ru;
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


#ifndef MODEL_GAS_SPECIES_H
#define MODEL_GAS_SPECIES_H

#include <AMReX_Array.H>
#include <cmath>
#include <string>

namespace Model {
namespace Gas {

constexpr double Ru = 8314.46;    // J/(kmol-K)
constexpr double boltzmann = 1.38066e-23; // J/K

// ===================================================================================
// 1. Thermo model: NASA-7 polynomial coefficients
// ===================================================================================
struct ThermoNASA7
{
    amrex::GpuArray<double, 3> Trange;   // [T_low, T_mid, T_high]
    amrex::GpuArray<double, 14> coeffs;  // 2 x 7 coefficients (flattened)

    AMREX_GPU_HOST_DEVICE
    ThermoNASA7(const double Trange_in[3],
                const double coeffs_in[2][7])
    {
        for (int i=0; i<3; ++i) Trange[i] = Trange_in[i];
        for (int i=0; i<2; ++i)
            for (int j=0; j<7; ++j)
                coeffs[i*7 + j] = coeffs_in[i][j];
    }

    AMREX_GPU_HOST_DEVICE
    double cp_R(double T) const noexcept {
        int hilo = (T <= Trange[1]) ? 0 : 1;
        const double* a = coeffs.data() + 7*hilo;
        return a[0] + a[1]*T + a[2]*pow(T,2.0)
             + a[3]*pow(T,3.0) + a[4]*pow(T,4.0);
    }

    AMREX_GPU_HOST_DEVICE
    double h_RT(double T) const noexcept {
        int hilo = (T <= Trange[1]) ? 0 : 1;
        const double* a = coeffs.data() + 7*hilo;
        return a[0] + 0.5*a[1]*T + (1.0/3.0)*a[2]*pow(T,2.0)
             + 0.25*a[3]*pow(T,3.0) + 0.2*a[4]*pow(T,4.0) + a[5]/T;
    }

    AMREX_GPU_HOST_DEVICE
    double s_R(double T) const noexcept {
        int hilo = (T <= Trange[1]) ? 0 : 1;
        const double* a = coeffs.data() + 7*hilo;
        return a[0]*log(T) + a[1]*T + 0.5*a[2]*pow(T,2.0)
             + (1.0/3.0)*a[3]*pow(T,3.0) + 0.25*a[4]*pow(T,4.0) + a[6];
    }
};

// ===================================================================================
// 2. Thermo model: Constant cp, with h(T) = cp*(T - Tref) + h0
// ===================================================================================
struct ThermoConst
{
    double cp_const;   // J/(kmol-K)
    double h0;         // J/kmol  (reference enthalpy)
    double s0;         // J/(kmol-K)  (reference entropy)
    double Tref;       // reference temperature (K)

    AMREX_GPU_HOST_DEVICE
    ThermoConst(double cp_in, double h0_in, double s0_in, double Tref_in = 298.15)
        : cp_const(cp_in), h0(h0_in), s0(s0_in), Tref(Tref_in) {}

    AMREX_GPU_HOST_DEVICE
    double cp_R(double /*T*/) const noexcept { return cp_const / Ru; }

    AMREX_GPU_HOST_DEVICE
    double h_RT(double T) const noexcept {
        double hT = cp_const * (T - Tref) + h0;
        return hT / (Ru * T);
    }

    AMREX_GPU_HOST_DEVICE
    double s_R(double T) const noexcept {
        double sT = cp_const * log(T / Tref) + s0;
        return sT / Ru;
    }
};

// ===================================================================================
// 3. Main templated Species struct
// ===================================================================================
template <typename ThermoModel>
struct Composition
{
    amrex::GpuArray<char, 12> species;  // Symbol (e.g. "H2O")
    double _MW;              // kg/kmol
    double _LJ_welldepth;    // K
    double _LJ_diameter;     // Å
    ThermoModel thermo;     // embedded thermodynamic model

    // Host-side constructor
    template <typename... Args>
    Composition(const char* name,
                double MW_in, double LJw_in, double LJd_in,
                Args&&... args)
        : _MW(MW_in),
          _LJ_welldepth(LJw_in),
          _LJ_diameter(LJd_in),
          thermo(std::forward<Args>(args)...)
    {
        for (int i=0; i<12; ++i) species[i] = (name[i] != '\0') ? name[i] : '\0';
    }

    // ---------- Thermodynamic properties ----------
    AMREX_GPU_HOST_DEVICE
    double MW(const amrex::Array4<double>& /* mole fraction */,
              const int /* i */, const int /* j */, const int /* k */,
              const int /* nspecies */ ) const noexcept {
        return _MW;
    }

    AMREX_GPU_HOST_DEVICE
    double R(const amrex::Array4<double>& /* mole fraction */,
             const int /* i */, const int /* j */, const int /* k */,
             const int /* nspecies */ ) const noexcept {
        return Ru/_MW;
    }

    AMREX_GPU_HOST_DEVICE
    double cp_mol(const double T, const amrex::Array4<double>& /* mole fraction */,
                            const int /* i */, const int /* j */, const int /* k */,
                            const int /* nspecies */ ) const noexcept {
        return thermo.cp_R(T) * Ru;
    }

    AMREX_GPU_HOST_DEVICE
    double cp_mass(const double T, const amrex::Array4<double>& /* mole fraction */,
                   const int /* i */, const int /* j */, const int /* k */,
                   const int /* nspecies */ ) const noexcept {
        return cp_mol(T, nullptr, 0, 0, 0, 0) / _MW;
    }

    AMREX_GPU_HOST_DEVICE
    double cv_mol(const double T, const amrex::Array4<double>& /* mole fraction */,
                  const int /* i */, const int /* j */, const int /* k */,
                  const int /* nspecies */ ) const noexcept {
        return cp_mol(T, nullptr, 0, 0, 0, 0) - Ru;
    }

    AMREX_GPU_HOST_DEVICE
    double cv_mass(const double T, const amrex::Array4<double>& /* mole fraction */,
                   const int /* i */, const int /* j */, const int /* k */,
                   const int /* nspecies */ ) const noexcept {
        return cv_mol(T, nullptr, 0, 0, 0, 0) / _MW;
    }

    AMREX_GPU_HOST_DEVICE
    double gamma(const double T, const amrex::Array4<double>& /* mole fraction */,
                 const int /* i */, const int /* j */, const int /* k */,
                 const int /* nspecies */ ) const noexcept {
        return cp_mol(T, nullptr, 0, 0, 0, 0) / cv_mol(T, nullptr, 0, 0, 0, 0);
    }

    AMREX_GPU_HOST_DEVICE
    double enthalpy_mol(const double T, const amrex::Array4<double>& /* mole fraction */,
                        const int /* i */, const int /* j */, const int /* k */,
                        const int /* nspecies */ ) const noexcept {
        return thermo.h_RT(T) * Ru * T;
    }

    AMREX_GPU_HOST_DEVICE
    double enthalpy_mass(const double T, const amrex::Array4<double>& /* mole fraction */,
                         const int /* i */, const int /* j */, const int /* k */,
                         const int /* nspecies */ ) const noexcept {
        return enthalpy_mol(T, nullptr, 0, 0, 0, 0) / _MW;
    }

    AMREX_GPU_HOST_DEVICE
    double entropy_mol(const double T, const amrex::Array4<double>& /* mole fraction */,
                       const int /* i */, const int /* j */, const int /* k */,
                       const int /* nspecies */ ) const noexcept {
        return thermo.s_R(T) * Ru;
    }

    AMREX_GPU_HOST_DEVICE
    double entropy_mass(const double T, const amrex::Array4<double>& /* mole fraction */,
                        const int /* i */, const int /* j */, const int /* k */,
                        const int /* nspecies */ ) const noexcept {
        return entropy_mol(T, nullptr, 0, 0, 0, 0) / _MW;
    }

    // ---------- Transport properties ----------
    AMREX_GPU_HOST_DEVICE
    double dynamic_viscosity(const double T, const amrex::Array4<double>& /* mole fraction */,
                             const int /* i */, const int /* j */, const int /* k */,
                             const int /* nspecies */ ) const noexcept {
        double Tstar = T / _LJ_welldepth;
        double sigma_star =
            1.16145 * pow(Tstar, -0.14874)
          + 0.52487 / exp(0.77320 * Tstar)
          + 2.16178 / exp(2.43787 * Tstar);
        return 2.6693e-5 * sqrt(_MW * T) /
               (pow(_LJ_diameter, 2.0) * sigma_star);
    }

    AMREX_GPU_HOST_DEVICE
    double thermal_conductivity(const double T, const amrex::Array4<double>& /* mole fraction */,
                                const int /* i */, const int /* j */, const int /* k */,
                                const int /* nspecies */ ) const noexcept {
        return gamma(T, nullptr, 0, 0, 0, 0) * dynamic_viscosity(T, nullptr, 0, 0, 0, 0) * R(nullptr, 0, 0, 0, 0);
    }

    // ---------- Composition properties ----------
    AMREX_GPU_HOST_DEVICE
    void ComputeLocalFractions(amrex::Array4<double>& X, amrex::Array4<double>& Y,
                               const amrex::Array4<double>& /* rhoY */, 
                               const int i, const int j, const int k, const int /* nspecies */) {
        AMREX_GPU_ASSERT(nspecies == 1);
        X(i,j,k,0) = 1.0;
        Y(i,j,k,0) = 1.0;
    }
};

} // namespace Gas
} // namespace Model

#endif

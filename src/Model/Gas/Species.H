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

    // Host-side constructor
    ThermoNASA7(const amrex::Array1D<double,0,2>& Trange_in,
                const amrex::Array2D<double,0,1,0,6>& coeffs_in)
    {
        for (int i=0; i<3; ++i) Trange[i] = Trange_in(i);
        for (int i=0; i<2; ++i)
            for (int j=0; j<7; ++j)
                coeffs[i*7 + j] = coeffs_in(i,j);
    }

    AMREX_GPU_HOST_DEVICE
    double cp_R(double Tval) const noexcept {
        int hilo = (Tval <= Trange[1]) ? 0 : 1;
        const double* a = coeffs.data() + 7*hilo;
        return a[0] + a[1]*Tval + a[2]*pow(Tval,2.0)
             + a[3]*pow(Tval,3.0) + a[4]*pow(Tval,4.0);
    }

    AMREX_GPU_HOST_DEVICE
    double h_RT(double Tval) const noexcept {
        int hilo = (Tval <= Trange[1]) ? 0 : 1;
        const double* a = coeffs.data() + 7*hilo;
        return a[0] + 0.5*a[1]*Tval + (1.0/3.0)*a[2]*pow(Tval,2.0)
             + 0.25*a[3]*pow(Tval,3.0) + 0.2*a[4]*pow(Tval,4.0) + a[5]/Tval;
    }

    AMREX_GPU_HOST_DEVICE
    double s_R(double Tval) const noexcept {
        int hilo = (Tval <= Trange[1]) ? 0 : 1;
        const double* a = coeffs.data() + 7*hilo;
        return a[0]*log(Tval) + a[1]*Tval + 0.5*a[2]*pow(Tval,2.0)
             + (1.0/3.0)*a[3]*pow(Tval,3.0) + 0.25*a[4]*pow(Tval,4.0) + a[6];
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
struct Species
{
    amrex::GpuArray<char, 12> species;  // Symbol (e.g. "H2O")
    double MW;              // kg/kmol
    double LJ_welldepth;    // K
    double LJ_diameter;     // Å
    double R;               // J/(kg-K)
    double T;               // current temperature
    ThermoModel thermo;     // embedded thermodynamic model

    // Host-side constructor
    template <typename... Args>
    Species(const std::string& name,
            double MW_in, double LJw_in, double LJd_in, double T_in=300.0,
            Args&&... args)
        : MW(MW_in),
          LJ_welldepth(LJw_in),
          LJ_diameter(LJd_in),
          R(Ru/MW_in),
          T(T_in),
          thermo(std::forward<Args>(args)...)
    {
        for (int i=0; i<12; ++i)
            species[i] = (i < static_cast<int>(name.size())) ? name[i] : '\0';
    }

    // ---------- Thermodynamic properties ----------
    AMREX_GPU_HOST_DEVICE double cp_mol(double Tval) const noexcept {
        return thermo.cp_R(Tval) * Ru;
    }

    AMREX_GPU_HOST_DEVICE double cp_mass(double Tval) const noexcept {
        return cp_mol(Tval) / MW;
    }

    AMREX_GPU_HOST_DEVICE double cv_mol(double Tval) const noexcept {
        return cp_mol(Tval) - Ru;
    }

    AMREX_GPU_HOST_DEVICE double cv_mass(double Tval) const noexcept {
        return cv_mol(Tval) / MW;
    }

    AMREX_GPU_HOST_DEVICE double gamma(double Tval) const noexcept {
        return cp_mol(Tval) / cv_mol(Tval);
    }

    AMREX_GPU_HOST_DEVICE double enthalpy_mol(double Tval) const noexcept {
        return thermo.h_RT(Tval) * Ru * Tval;
    }

    AMREX_GPU_HOST_DEVICE double enthalpy_mass(double Tval) const noexcept {
        return enthalpy_mol(Tval) / MW;
    }

    AMREX_GPU_HOST_DEVICE double entropy_mol(double Tval) const noexcept {
        return thermo.s_R(Tval) * Ru;
    }

    AMREX_GPU_HOST_DEVICE double entropy_mass(double Tval) const noexcept {
        return entropy_mol(Tval) / MW;
    }

    // ---------- Transport properties ----------
    AMREX_GPU_HOST_DEVICE double dynamic_viscosity(double Tval) const noexcept {
        double Tstar = Tval / LJ_welldepth;
        double sigma_star =
            1.16145 * pow(Tstar, -0.14874)
          + 0.52487 / exp(0.77320 * Tstar)
          + 2.16178 / exp(2.43787 * Tstar);
        return 2.6693e-5 * sqrt(MW * Tval) /
               (pow(LJ_diameter, 2.0) * sigma_star);
    }

    AMREX_GPU_HOST_DEVICE double thermal_conductivity(double Tval) const noexcept {
        return gamma(Tval) * dynamic_viscosity(Tval) * R;
    }

    AMREX_GPU_HOST_DEVICE void UpdateTemperature(double Tnew) noexcept { T = Tnew; }
};

} // namespace Gas
} // namespace Model

#endif

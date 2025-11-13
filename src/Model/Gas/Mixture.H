#ifndef MODEL_GAS_MIXTURE_H
#define MODEL_GAS_MIXTURE_H

#include <AMReX_GpuQualifiers.H>
#include <cmath>
#include "Model/Gas/Species.H"

namespace Model {
namespace Gas {

// ===================================================================================
// GPU-safe Mixture Composition class (preallocated device arrays for Y, X, DKM)
// ===================================================================================
template <typename ThermoModel>
struct CompositionMixture
{
    // --- Constant pointers (device accessible)
    const Composition<ThermoModel>* species_data;  // array of species

    // --- Scalars
    int nspecies;
    double _rho_mix;    // mixture density [kg/m³]
    double _mu_mix;     // dynamic viscosity [Pa-s]
    double _k_mix;      // thermal conductivity [W/m-K]

    // --- Constructor (GPU-safe)
    AMREX_GPU_HOST_DEVICE
    CompositionMixture(const Composition<ThermoModel>* sdata)
        : species_data(sdata),
          _rho_mix(0.0),
          _mu_mix(0.0),
          _k_mix(0.0)
    {}

    // ===================================================================================
    // --- Thermodynamic properties
    // ===================================================================================
    AMREX_GPU_HOST_DEVICE
    double MW(const amrex::Array4<double>& X,
              const int i, const int j, const int k,
              const int nspecies) const noexcept {
        double _MW = 0.0;
        for (int n = 0; n < nspecies; ++n) _MW += X(i,j,k,n) * species_data[n].MW(nullptr, 0, 0, 0, 0);
        return _MW;
    }

    AMREX_GPU_HOST_DEVICE
    double R(const amrex::Array4<double>& X, 
             const int i, const int j, const int k,
             const int nspecies) const noexcept { return Ru / MW(X, i, j, k, nspecies); }

    AMREX_GPU_HOST_DEVICE
    double cp_mol(const double T, const amrex::Array4<double>& X, 
                  const int i, const int j, const int k,
                  const int nspecies) const noexcept {
        double cp_mix = 0.0, sumx = 0.0;
        double MWmix = MW();
        for (int n = 0; n < nspecies; ++n) {
            double x = X(i,j,k,n);
            cp_mix += x * species_data[n].cp_mol(T);
            sumx += x;
        }
        return cp_mix / sumx;
    }

    AMREX_GPU_HOST_DEVICE
    double cp_mass(double T) const noexcept { return cp_mol(T) / MW(); }

    AMREX_GPU_HOST_DEVICE
    double cv_mol(double T) const noexcept { return cp_mol(T) - Ru; }

    AMREX_GPU_HOST_DEVICE
    double cv_mass(double T) const noexcept { return cv_mol(T) / MW(); }

    AMREX_GPU_HOST_DEVICE
    double gamma(double T) const noexcept { return cp_mol(T) / cv_mol(T); }

    AMREX_GPU_HOST_DEVICE
    double enthalpy_mol(double T) const noexcept {
        double h_mix = 0.0, sumx = 0.0;
        for (int n = 0; n < nspecies; ++n) {
            double x = X[n];
            h_mix += x * species_data[n].enthalpy_mol(T);
            sumx += x;
        }
        return h_mix / sumx;
    }

    AMREX_GPU_HOST_DEVICE
    double enthalpy_mass(double T) const noexcept { return enthalpy_mol(T) / MW(); }

    AMREX_GPU_HOST_DEVICE
    double entropy_mol(double T) const noexcept {
        double s_mix = 0.0, sumx = 0.0;
        for (int n = 0; n < nspecies; ++n) {
            double x = X[n];
            s_mix += x * species_data[n].entropy_mol(T);
            sumx += x;
        }
        return s_mix / sumx;
    }

    AMREX_GPU_HOST_DEVICE
    double entropy_mass(double T) const noexcept { return entropy_mol(T) / MW(); }

    // ===================================================================================
    // --- Compute transport properties (μ, k, DKM)
    // ===================================================================================
    AMREX_GPU_HOST_DEVICE
    void compute_transport_properties(double T, double P) noexcept
    {
        mu_mix = 0.0;
        k_mix  = 0.0;

        // Wilke's mixture averaging
        for (int a = 0; a < nspecies; ++a) {
            double mua = species_data[a].dynamic_viscosity(T);
            double Mwa = species_data[a].MW();
            double phi = 0.0;

            for (int b = 0; b < nspecies; ++b) {
                double mub = species_data[b].dynamic_viscosity(T);
                double Mwb = species_data[b].MW();
                double ratio = pow(Mwb / Mwa, 0.25);
                double mu_ratio = sqrt(mua / mub);
                double term = sqrt(1.0 + Mwa / Mwb);
                phi += X[b] * pow(1.0 + mu_ratio * ratio, 2.0) / (sqrt(8.0) * term);
            }

            mu_mix += X[a] * mua / phi;
            k_mix  += X[a] * species_data[a].thermal_conductivity(T) / phi;
        }

        // Mixture diffusion coefficients (Chapman–Enskog)
        for (int a = 0; a < nspecies; ++a) {
            double sumD = 0.0;
            for (int b = 0; b < nspecies; ++b) {
                if (a == b) continue;

                double sigmaAB = 0.5 * (species_data[a]._LJ_diameter + species_data[b]._LJ_diameter);
                double epsAB   = sqrt(species_data[a]._LJ_welldepth * species_data[b]._LJ_welldepth);
                double Tstar   = T / epsAB;

                double omegaAB =
                    1.06036 / pow(Tstar, 0.15610) +
                    0.19300 / exp(0.47635 * Tstar) +
                    1.03587 / exp(1.52996 * Tstar) +
                    1.76474 / exp(3.89411 * Tstar);

                double DAB = 0.0018583 * sqrt(pow(T,3.0) * (1.0 / species_data[a].MW() + 1.0 / species_data[b].MW()))
                            / ( (P / 101325.0) * pow(sigmaAB,2.0) * omegaAB );
                DAB /= 1.0e4;

                sumD += X[b] / DAB;
            }

            DKM[a] = (1.0 - X[a]) / sumD;
            if (!(DKM[a] == DKM[a])) DKM[a] = 0.0; // NaN guard
        }
    }

    // ===================================================================================
    // --- Access transport properties
    // ===================================================================================
    AMREX_GPU_HOST_DEVICE double dynamic_viscosity() const noexcept { return mu_mix; }
    AMREX_GPU_HOST_DEVICE double thermal_conductivity() const noexcept { return k_mix; }
    AMREX_GPU_HOST_DEVICE double diffusion_coeff(int i) const noexcept { return DKM[i]; }
};

} // namespace Gas
} // namespace Model

#endif


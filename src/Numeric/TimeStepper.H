#ifndef NUMERIC_TIMESTEPPER_H_
#define NUMERIC_TIMESTEPPER_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Box.H>
#include <AMReX_MFIter.H>
#include <AMReX_Array4.H>
#include "Set/Set.H"
#include "Integrator/Integrator.H"

namespace Numeric {
// Template-based base class for time-stepping schemes
template <typename IntegratorDerivedClass>
class TimeStepper {
public:
    virtual void ComputeSubStep(int lev, Set::Scalar dt, int stage, IntegratorDerivedClass* solver) = 0;
    virtual ~TimeStepper() = default;

    struct EulerForwardConstants {
        static constexpr int stages = 1;
        static constexpr Set::Scalar a[] = {1.0};
        static constexpr Set::Scalar b[] = {0.0};
        static constexpr Set::Scalar c[] = {1.0};
    };

    struct RK3Constants {
        static constexpr int stages = 3;
        static constexpr Set::Scalar a[] = {1.0, 0.25, 2.0 / 3.0};
        static constexpr Set::Scalar b[] = {0.0, 0.75, 1.0 / 3.0};
        static constexpr Set::Scalar c[] = {1.0, 0.25, 2.0 / 3.0};
    };
};

// Euler Forward time-stepping scheme
template <typename IntegratorDerivedClass>
class EulerForwardScheme : public TimeStepper<IntegratorDerivedClass> {
public:
    void ComputeSubStep(int lev, Set::Scalar dt, int stage, IntegratorDerivedClass* solver) override {
        if (!solver) {
            amrex::Abort("Solver pointer is null.");
            return;
        }

        const Set::Scalar* DX = solver->geom[lev].CellSize();

        for (amrex::MFIter mfi(*solver->QVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();
            auto const& q_arr = solver->QVec_mf.Patch(lev, mfi);
            auto const& q_old_arr = solver->QVec_old_mf.Patch(lev, mfi);
            auto const& xflux_arr = solver->XFlux_mf.Patch(lev, mfi);
            auto const& yflux_arr = solver->YFlux_mf.Patch(lev, mfi);
#if (AMREX_SPACEDIM == 3)
            auto const& zflux_arr = solver->ZFlux_mf.Patch(lev, mfi);
#endif

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                for (int n = 0; n < solver->number_of_components; ++n) {
                    Set::Scalar flux_div = 0.0;
                    if (solver->XFlux_mf[lev]) {
                        flux_div += (xflux_arr(i, j, k, n) - xflux_arr(i - 1, j, k, n)) / DX[0];
                    }
                    if (solver->YFlux_mf[lev]) {
                        flux_div += (yflux_arr(i, j, k, n) - yflux_arr(i, j - 1, k, n)) / DX[1];
                    }
#if (AMREX_SPACEDIM == 3)
                    if (solver->ZFlux_mf[lev]) {
                        flux_div += (zflux_arr(i, j, k, n) - zflux_arr(i, j, k - 1, n)) / DX[2];
                    }
#endif
                    q_arr(i, j, k, n) = q_old_arr(i, j, k, n) - dt * flux_div;
                }
            });
        }
    }
};

// RK3 time-stepping scheme
template <typename IntegratorDerivedClass>
class RK3Scheme : public TimeStepper<IntegratorDerivedClass> {
public:
    void ComputeSubStep(int lev, Set::Scalar dt, int stage, IntegratorDerivedClass* solver) override {
        if (!solver) {
            amrex::Abort("Solver pointer is null.");
            return;
        }

        const Set::Scalar* DX = solver->geom[lev].CellSize();

        for (amrex::MFIter mfi(*solver->QVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();
            auto const& q_arr = solver->QVec_mf.Patch(lev, mfi);
            auto const& q_old_arr = solver->QVec_old_mf.Patch(lev, mfi);
            auto const& xflux_arr = solver->XFlux_mf.Patch(lev, mfi);
            auto const& yflux_arr = solver->YFlux_mf.Patch(lev, mfi);
#if (AMREX_SPACEDIM == 3)
            auto const& zflux_arr = solver->ZFlux_mf.Patch(lev, mfi);
#endif

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                for (int n = 0; n < solver->number_of_components; ++n) {
                    Set::Scalar flux_div = 0.0;
                    if (solver->XFlux_mf[lev]) {
                        flux_div += (xflux_arr(i, j, k, n) - xflux_arr(i - 1, j, k, n)) / DX[0];
                    }
                    if (solver->YFlux_mf[lev]) {
                        flux_div += (yflux_arr(i, j, k, n) - yflux_arr(i, j - 1, k, n)) / DX[1];
                    }
#if (AMREX_SPACEDIM == 3)
                    if (solver->ZFlux_mf[lev]) {
                        flux_div += (zflux_arr(i, j, k, n) - zflux_arr(i, j, k - 1, n)) / DX[2];
                    }
#endif
                    q_arr(i, j, k, n) = TimeStepper<IntegratorDerivedClass>::RK3Constants::a[stage] * q_old_arr(i, j, k, n) +
                                        TimeStepper<IntegratorDerivedClass>::RK3Constants::b[stage] * q_arr(i, j, k, n) -
                                        TimeStepper<IntegratorDerivedClass>::RK3Constants::c[stage] * dt * flux_div;
                }
            });
        }
    }
};

}  // namespace Numeric

#endif  // NUMERIC_TIMESTEPPER_H_


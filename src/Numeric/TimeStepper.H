#ifndef NUMERIC_TIMESTEPPER_H_
#define NUMERIC_TIMESTEPPER_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Box.H>
#include <AMReX_MFIter.H>
#include <AMReX_Array4.H>
#include "Set/Set.H"
#include "Integrator/Integrator.H"
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"
#include <memory>
#include <string>

namespace Numeric {

// Abstract interface for time-stepping schemes
template <typename IntegratorDerivedClass>
class TimeSteppingScheme {
public:
    virtual void ComputeSubStep(int lev, Set::Scalar dt, int stage, IntegratorDerivedClass* solver) = 0;
    virtual int GetNumberOfStages() const = 0;
    virtual ~TimeSteppingScheme() = default;
};

// TimeStepper base class (like FluxHandler)
template <typename IntegratorDerivedClass>
class TimeStepper {
private:
    std::shared_ptr<TimeSteppingScheme<IntegratorDerivedClass>> scheme;  // Shared pointer to time-stepping scheme

public:
    TimeStepper() = default;

    // Set the time-stepping scheme
    void SetTimeSteppingScheme(const std::shared_ptr<TimeSteppingScheme<IntegratorDerivedClass>>& newScheme) {
        scheme = newScheme;
    }

    // Call the scheme's ComputeSubStep directly
    void ComputeSubStep(int lev, Set::Scalar dt, int stage, IntegratorDerivedClass* solver) {
        if (!scheme) {
            Util::Abort(__FILE__, __func__, __LINE__, "TimeSteppingScheme is not set.");
        }
        scheme->ComputeSubStep(lev, dt, stage, solver);
    }

    int GetNumberOfStages() const {
        if (!scheme) {
            Util::Abort(__FILE__, __func__, __LINE__, "TimeSteppingScheme is not set.");
        }
        return scheme->GetNumberOfStages();
    }
};

// EulerForwardScheme
template <typename IntegratorDerivedClass>
class EulerForwardScheme : public TimeSteppingScheme<IntegratorDerivedClass> {
public:
    void ComputeSubStep(int lev, Set::Scalar dt, int stage, IntegratorDerivedClass* solver) override {
        if (stage != 0) {
            Util::Abort(__FILE__, __func__, __LINE__, "Invalid stage for EulerForward: only one stage exists.");
        }

        Util::ScimitarX_Util::Debug debug;

        //debug.SetTargetDebugLocation(0.5, 0.2, 0.0, solver->geom[lev]);
        debug.SetTargetDebugLocationIndices(50, 0, 0);

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
                    Set::Scalar xflux_right = 0.0, xflux_left = 0.0;
                    Set::Scalar yflux_top = 0.0, yflux_bottom = 0.0;
                    Set::Scalar zflux_front = 0.0, zflux_back = 0.0;
            
                    if (solver->XFlux_mf[lev]) {
                        xflux_right = xflux_arr(i, j, k, n);  // Flux at right face
                        xflux_left = xflux_arr(i - 1, j, k, n);  // Flux at left face
                        flux_div += (xflux_right - xflux_left) / DX[0];
                    }
                    if (solver->YFlux_mf[lev]) {
                        yflux_top = yflux_arr(i, j, k, n);  // Flux at top face
                        yflux_bottom = yflux_arr(i, j - 1, k, n);  // Flux at bottom face
                        flux_div += (yflux_top - yflux_bottom) / DX[1];
                    }
#if (AMREX_SPACEDIM == 3)
                    if (solver->ZFlux_mf[lev]) {
                        zflux_front = zflux_arr(i, j, k, n);  // Flux at front face
                        zflux_back = zflux_arr(i, j, k - 1, n);  // Flux at back face
                        flux_div += (zflux_front - zflux_back) / DX[2];
                    }
#endif
            
                    Set::Scalar q_old = q_old_arr(i, j, k, n);
                    Set::Scalar dt_val = dt;
            
                    // Call the debug function
                    debug.DebugEulerForwardStep(i, j, k, xflux_right, xflux_left, yflux_top, yflux_bottom, 
                                                            zflux_front, zflux_back, q_old, dt_val, "Euler Forward Substep", true);
            
                    // Update solution
                    q_arr(i, j, k, n) = q_old - dt * flux_div;
                }
            });

        }
    }

    int GetNumberOfStages() const override {
        return 1;  // EulerForward has 1 stage
    }
};

// RK3 Scheme
template <typename IntegratorDerivedClass>
class RK3Scheme : public TimeSteppingScheme<IntegratorDerivedClass> {
public:
    struct RK3Constants {
        static constexpr int stages = 3;
        static constexpr Set::Scalar a[] = {1.0, 0.25, 2.0 / 3.0};
        static constexpr Set::Scalar b[] = {0.0, 0.75, 1.0 / 3.0};
        static constexpr Set::Scalar c[] = {1.0, 0.25, 2.0 / 3.0};
    };

    void ComputeSubStep(int lev, Set::Scalar dt, int stage, IntegratorDerivedClass* solver) override {
        if (stage >= RK3Constants::stages) {
            Util::Abort(__FILE__, __func__, __LINE__, "Invalid stage for RK3.");
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
                    q_arr(i, j, k, n) = RK3Constants::a[stage] * q_old_arr(i, j, k, n) +
                                        RK3Constants::b[stage] * q_arr(i, j, k, n) -
                                        RK3Constants::c[stage] * dt * flux_div;
                }
            });
        }
    }

    int GetNumberOfStages() const override {
        return RK3Constants::stages;
    }
};

}  // namespace Numeric

#endif  // NUMERIC_TIMESTEPPER_H_


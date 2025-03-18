#ifndef NUMERIC_TIMESTEPPER_H_
#define NUMERIC_TIMESTEPPER_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Box.H>
#include <AMReX_MFIter.H>
#include <AMReX_Array4.H>
#include "Set/Set.H"
#include "Integrator/Integrator.H"
#include "Integrator/NarrowBandLevelset.H"
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"
#include <memory>
#include <string>

// Include foward declare for ScimitarX.H
namespace Integrator {
    class ScimitarX;
    class NarrowBandLevelset;
}

namespace Numeric {

/*// Helper function to compute face indices for a control volume's flux calculation
AMREX_GPU_DEVICE
inline void GetFaceIndices(int i, int j, int k, amrex::IntVect& west_face, amrex::IntVect& east_face,
                            amrex::IntVect& south_face, amrex::IntVect& north_face,
                            amrex::IntVect& bottom_face, amrex::IntVect& top_face) {
    // West and east face indices for the x-direction
    west_face = amrex::IntVect(AMREX_D_DECL(i, j, k));
    east_face = amrex::IntVect(AMREX_D_DECL(i + 1, j, k));

    // South and north face indices for the y-direction
    south_face = amrex::IntVect(AMREX_D_DECL(i, j, k));
    north_face = amrex::IntVect(AMREX_D_DECL(i, j + 1, k));

#if (AMREX_SPACEDIM == 3)
    // Bottom and top face indices for the z-direction
    bottom_face = amrex::IntVect(AMREX_D_DECL(i, j, k));
    top_face = amrex::IntVect(AMREX_D_DECL(i, j, k + 1));
#else
    bottom_face = amrex::IntVect(AMREX_D_DECL(0, 0, 0));  // Placeholder for 2D
    top_face = amrex::IntVect(AMREX_D_DECL(0, 0, 0));     // Placeholder for 2D
#endif
}*/

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

        //Util::ScimitarX_Util::Debug debug;
        //debug.SetTargetDebugLocationIndices(101, 0, 0);
        //debug.SetTargetDebugLocation(0.49, 0.12, 0.0, solver->geom[lev]);

        const Set::Scalar* DX = solver->geom[lev].CellSize();  // Grid spacing

        //for (amrex::MFIter mfi(*solver->QVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        for (amrex::MFIter mfi(*solver->QVec_mf[lev], false); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();  // Valid region of the box
            auto const& q_arr = solver->QVec_mf.Patch(lev, mfi);  // Updated conserved variables
            auto const& q_old_arr = solver->QVec_old_mf.Patch(lev, mfi);  // Old conserved variables
            auto const& xflux_arr = solver->XFlux_mf.Patch(lev, mfi);  // X-direction flux
            auto const& yflux_arr = solver->YFlux_mf.Patch(lev, mfi);  // Y-direction flux
#if (AMREX_SPACEDIM == 3)
            auto const& zflux_arr = solver->ZFlux_mf.Patch(lev, mfi);  // Z-direction flux
#endif

            // Parallel loop over cells in the valid box
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                for (int n = 0; n < solver->number_of_components; ++n) {
                    // **Get Face Indices for the Control Volume**
                    // amrex::IntVect west_face, east_face, south_face, north_face, bottom_face, top_face;
                    // GetFaceIndices(i, j, k, west_face, east_face, south_face, north_face, bottom_face, top_face);

                    // **Flux Divergence Computation**
                    Set::Scalar flux_div = 0.0;

                    // X-direction flux divergence
                    flux_div += (xflux_arr(i, j, k, n) - xflux_arr(i-1, j, k, n)) / DX[0];

                    // Y-direction flux divergence
                    flux_div += (yflux_arr(i, j, k, n) - yflux_arr(i, j-1, k, n)) / DX[1];

#if (AMREX_SPACEDIM == 3)
                    // Z-direction flux divergence
                    flux_div += (zflux_arr(i, j, k, n) - zflux_arr(i, j, k-1, n)) / DX[2];
#endif

                    // Update conserved variable using Euler Forward step
                    Set::Scalar q_old = q_old_arr(i, j, k, n);
                    q_arr(i, j, k, n) = q_old - dt * flux_div;

                    //Debug output
                    /* debug.DebugEulerForwardStep(i, j, k, xflux_arr(i,j,k,n), xflux_arr(i-1,j,k,n),
                                                yflux_arr(i,j,k,n), yflux_arr(i,j-1,k,n),
                                                0.0,
                                                0.0,
                                                q_old, dt, "Euler Forward Substep", true); */               
                }
            });
        }
    }

    int GetNumberOfStages() const override {
        return 1;  // Euler Forward has only one stage
    }
};

template <>
class EulerForwardScheme<Integrator::NarrowBandLevelset> : public TimeSteppingScheme<Integrator::NarrowBandLevelset> {
public:
    void ComputeSubStep(int lev, Set::Scalar dt, int stage, Integrator::NarrowBandLevelset* solver) override {
        if (stage != 0) {
            Util::Abort(__FILE__, __func__, __LINE__, "Invalid stage for EulerForward: only one stage exists.");
        }
        const Set::Scalar* DX = solver->geom[lev].CellSize();  // Grid spacing

        for (amrex::MFIter mfi(*solver->ls_mf[lev], false); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();  // Valid region of the box
            auto const& ls_arr = solver->ls_mf.Patch(lev, mfi);
            auto const& ls_old_arr = solver->ls_old_mf.Patch(lev, mfi);

            // Get face-centered flux fields
            auto const& xflux_arr = solver->XFlux_mf.Patch(lev, mfi);;  // X-direction flux
#if AMREX_SPACEDIM >= 2
            auto const& yflux_arr = solver->YFlux_mf.Patch(lev, mfi);  // Y-direction flux
#endif
#if AMREX_SPACEDIM == 3
            auto const& zflux_arr = solver->ZFlux_mf.Patch(lev, mfi);  // Z-direction flux
#endif

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                Set::Scalar flux_div = 0.0;

                // Compute flux divergence in X-direction
                flux_div += (xflux_arr(i, j, k, 0) - xflux_arr(i - 1, j, k, 0)) / DX[0];

#if AMREX_SPACEDIM >= 2
                // Compute flux divergence in Y-direction
                flux_div += (yflux_arr(i, j, k, 0) - yflux_arr(i, j - 1, k, 0)) / DX[1];
#endif

#if AMREX_SPACEDIM == 3
                // Compute flux divergence in Z-direction
                flux_div += (zflux_arr(i, j, k, 0) - zflux_arr(i, j, k - 1, 0)) / DX[2];
#endif

                // Euler forward update for level set function
                ls_arr(i, j, k, 0) = ls_old_arr(i, j, k, 0) - dt * flux_div;
            });
        }
    }

    int GetNumberOfStages() const override {
        return 1;  // Euler Forward has only one stage
    }
};

/*// EulerForwardScheme
template <>
class EulerForwardScheme<Integrator::NarrowBandLevelset> : public TimeSteppingScheme<Integrator::NarrowBandLevelset> {
public:
    void ComputeSubStep(int lev, Set::Scalar dt, int stage, Integrator::NarrowBandLevelset* solver) override {
        if (stage != 0) {
            Util::Abort(__FILE__, __func__, __LINE__, "Invalid stage for EulerForward: only one stage exists.");
        }
        const Set::Scalar* DX = solver->geom[lev].CellSize();  // Grid spacing

        for (amrex::MFIter mfi(*solver->ls_mf[lev], false); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();  // Valid region of the box
            auto const& ls_arr = solver->ls_mf[lev]->array(mfi);  // Current level set field
            auto const& ls_old_arr = solver->ls_old_mf[lev]->array(mfi);  // Old level set field
            auto const& flux_arr = solver->flux_mf[lev]->array(mfi);  // Computed fluxes 

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Compute the flux divergence (advection equation: ∂ϕ/∂t + ∇⋅F = 0)
                Set::Scalar flux_div = 0.0;

                // Loop over spatial dimensions
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    flux_div += (flux_arr(i, j, k, d) - flux_arr(i - (d == 0), j - (d == 1), k - (d == 2), d)) / DX[d];
                }

                // Euler forward update for level set function
                ls_arr(i, j, k) = ls_old_arr(i, j, k) - dt * flux_div;
            });
        }
    }

    int GetNumberOfStages() const override {
        return 1;  // Euler Forward has only one stage
    }
};*/

// RK3 Scheme
template <typename IntegratorDerivedClass>
class RK3Scheme : public TimeSteppingScheme<IntegratorDerivedClass> {
public:
    struct RK3Constants {
        static constexpr int stages = 3;  // RK3 has 3 stages
        static constexpr Set::Scalar a[] = {1.0, 0.25, 2.0 / 3.0};  // Coefficients
        static constexpr Set::Scalar b[] = {0.0, 0.75, 1.0 / 3.0};
        static constexpr Set::Scalar c[] = {1.0, 0.25, 2.0 / 3.0};
    };

    void ComputeSubStep(int lev, Set::Scalar dt, int stage, IntegratorDerivedClass* solver) override {
        if (stage >= RK3Constants::stages) {
            Util::Abort(__FILE__, __func__, __LINE__, "Invalid stage for RK3.");
        }

        const Set::Scalar* DX = solver->geom[lev].CellSize();  // Grid spacing

        for (amrex::MFIter mfi(*solver->QVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();
            auto const& q_arr = solver->QVec_mf.Patch(lev, mfi);  // Updated conserved variables
            auto const& q_old_arr = solver->QVec_old_mf.Patch(lev, mfi);  // Old conserved variables
            auto const& xflux_arr = solver->XFlux_mf.Patch(lev, mfi);  // X-direction flux
            auto const& yflux_arr = solver->YFlux_mf.Patch(lev, mfi);  // Y-direction flux
#if (AMREX_SPACEDIM == 3)
            auto const& zflux_arr = solver->ZFlux_mf.Patch(lev, mfi);  // Z-direction flux
#endif

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                for (int n = 0; n < solver->number_of_components; ++n) {
                    // **Flux Divergence Computation**
                    Set::Scalar flux_div = 0.0;
                    Set::Scalar xflux_right = 0.0, xflux_left = 0.0;
                    Set::Scalar yflux_top = 0.0, yflux_bottom = 0.0;
                    Set::Scalar zflux_front = 0.0, zflux_back = 0.0;

                    // X-direction flux divergence
                    if (solver->XFlux_mf[lev]) {
                        xflux_right = xflux_arr(i, j, k, n);  // Right face flux
                        xflux_left = xflux_arr(i - 1, j, k, n);  // Left face flux
                        flux_div += (xflux_right - xflux_left) / DX[0];
                    }

                    // Y-direction flux divergence (if 2D or 3D)
                    if (solver->YFlux_mf[lev]) {
                        yflux_top = yflux_arr(i, j, k, n);  // Top face flux
                        yflux_bottom = yflux_arr(i, j - 1, k, n);  // Bottom face flux
                        flux_div += (yflux_top - yflux_bottom) / DX[1];
                    }

                    // Z-direction flux divergence (if 3D)
#if (AMREX_SPACEDIM == 3)
                    if (solver->ZFlux_mf[lev]) {
                        zflux_front = zflux_arr(i, j, k, n);  // Front face flux
                        zflux_back = zflux_arr(i, j, k - 1, n);  // Back face flux
                        flux_div += (zflux_front - zflux_back) / DX[2];
                    }
#endif

                    // **Runge-Kutta 3rd order update formula**
                    q_arr(i, j, k, n) = RK3Constants::a[stage] * q_old_arr(i, j, k, n) +
                                        RK3Constants::b[stage] * q_arr(i, j, k, n) -
                                        RK3Constants::c[stage] * dt * flux_div;
                }
            });
        }
    }

    int GetNumberOfStages() const override {
        return RK3Constants::stages;  // RK3 has 3 stages
    }
};

}  // namespace Numeric

#endif  // NUMERIC_TIMESTEPPER_H_

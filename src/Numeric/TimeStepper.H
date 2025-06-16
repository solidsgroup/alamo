#ifndef NUMERIC_TIMESTEPPER_H_
#define NUMERIC_TIMESTEPPER_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Box.H>
#include <AMReX_MFIter.H>
#include <AMReX_Array4.H>
#include "Set/Set.H"
#include "Integrator/Integrator.H"
#include "Integrator/NarrowBandLevelset.H"  // or whatever the correct path is
#include "Numeric/NumericFactory.H"
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
                    flux_div += (xflux_arr(i, j, k, n) - xflux_arr(i - 1, j, k, n)) / DX[0];

                    // Y-direction flux divergence
                    flux_div += (yflux_arr(i, j, k, n) - yflux_arr(i, j - 1, k, n)) / DX[1];

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
class EulerForwardScheme<Integrator::NarrowBandLevelset>
    : public TimeSteppingScheme<Integrator::NarrowBandLevelset>
{
public:
    void ComputeSubStep(int lev, Set::Scalar dt, int stage, Integrator::NarrowBandLevelset* solver) override {
        if (stage != 0) {
            Util::Abort(__FILE__, __func__, __LINE__,
                "Invalid stage for EulerForward: only one stage exists.");
        }

        const Set::Scalar* DX = solver->geom[lev].CellSize();  // Grid spacing
        const int ncomp = solver->number_of_components;
        //const int OuterGhost= solver->GetOuterGhost();

        for (int ils = 0; ils < ncomp; ++ils) {
            ComputeEulerUpdateForComponent(lev, ils, dt, solver, DX, 2);
        }
    }

    int GetNumberOfStages() const override {
        return 1;  // Euler Forward has only one stage
    }

private:
    void ComputeEulerUpdateForComponent(
        int lev, int ils, Set::Scalar dt,
        Integrator::NarrowBandLevelset* solver,
        const Set::Scalar* DX,
        int OuterGhost)
    {
        const auto& ls_data = solver->level_sets[lev][ils];
        const auto& ba = ls_data.narrowband_ba;
        const auto& dm = ls_data.narrowband_dm;

        amrex::iMultiFab narrowband(ba, dm, 1, 0);
        amrex::MultiFab LS(ba, dm, 1, 0);
        amrex::MultiFab LS_old(ba, dm, 1, 0);

        narrowband.ParallelCopy(*ls_data.Tube, 0, 0, 1, 0, 0);
        LS.ParallelCopy(*solver->ls_mf[lev], ils, 0, 1, 0, 0);
        LS_old.ParallelCopy(*solver->ls_old_mf[lev], ils, 0, 1, 0, 0);

        // Define face-centered multifabs from solver
        amrex::IntVect face_type_x = amrex::IntVect::TheCellVector();  // All 0s: cell-centered
        face_type_x[0] = 1;  // Set x-direction to node-centered (i.e., face-centered in x)

        amrex::BoxArray xflux_ba(ba);
        xflux_ba.convert(amrex::IndexType(face_type_x));
        amrex::MultiFab xflux(xflux_ba, dm, 1, 1);
        xflux.ParallelCopy(*solver->XFlux_mf[lev], ils, 0, 1, 1, 1);

#if AMREX_SPACEDIM >= 2
        amrex::IntVect face_type_y = amrex::IntVect::TheCellVector();  // All 0s: cell-centered
        face_type_y[1] = 1;  

        amrex::BoxArray yflux_ba(ba);
        yflux_ba.convert(amrex::IndexType(face_type_y));
        amrex::MultiFab yflux(yflux_ba, dm, 1, 1);
        yflux.ParallelCopy(*solver->YFlux_mf[lev], ils, 0, 1, 1, 1);
#endif

#if AMREX_SPACEDIM == 3
        amrex::IntVect face_type_z = amrex::IntVect::TheCellVector();  // All 0s: cell-centered
        face_type_z[2] = 1; 

        amrex::BoxArray zflux_ba(ba);
        zflux_ba.convert(amrex::IndexType(face_type_z));
        amrex::MultiFab zflux(zflux_ba, dm, 1, 1);
        zflux.ParallelCopy(*solver->ZFlux_mf[lev], ils, 0, 1, 1, 1);
#endif
        /*for (amrex::MFIter mfi(LS); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            const amrex::FArrayBox& ls_fab     = (LS)[mfi];
            const amrex::FArrayBox& ls_old_fab = (LS_old)[mfi];
    
            if (bx.contains(amrex::IntVect(20, 49))) {
                amrex::Print() << "Before TimeStepper:\n";
                amrex::Print() << "LS(20,49): " << ls_fab(amrex::IntVect(20, 49), 0) << "\n";
                amrex::Print() << "LS_OLD(20,49): " << ls_old_fab(amrex::IntVect(20, 49), 0) << "\n";
            }
        }*/

        for (amrex::MFIter mfi(LS); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.validbox();
            auto const& nb_arr = narrowband.const_array(mfi);
            auto const& ls_arr = LS.array(mfi);
            auto const& ls_old_arr = LS_old.const_array(mfi);
            auto const& xflux_arr = xflux.const_array(mfi);
#if AMREX_SPACEDIM >= 2
            auto const& yflux_arr = yflux.const_array(mfi);
#endif
#if AMREX_SPACEDIM == 3
            auto const& zflux_arr = zflux.const_array(mfi);
#endif

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::IntVect coord(AMREX_D_DECL(i,j,k));
                if (std::abs(nb_arr(i, j, k)) >= OuterGhost) return;

                Set::Scalar flux_div = (xflux_arr(i, j, k) - xflux_arr(i - 1, j, k)) / DX[0];
                
#if AMREX_SPACEDIM >= 2
                flux_div += (yflux_arr(i, j, k) - yflux_arr(i, j - 1, k)) / DX[1];
#endif
#if AMREX_SPACEDIM == 3
                flux_div += (zflux_arr(i, j, k) - zflux_arr(i, j, k - 1)) / DX[2];
#endif
                Set::Scalar LS_old = ls_arr(coord);
                ls_arr(coord) = LS_old - dt * flux_div;
                if ((i==46 &&j==46)){
                    amrex::Print() << "LS_old" << coord << ": " << LS_old << std::endl;
                    amrex::Print() << "xflux" << coord << ": " << xflux_arr(coord) << std::endl;
                    amrex::Print() << "xlfux" << coord << "(-1): " << xflux_arr(i-1, j, k) << std::endl;
                    amrex::Print() << "LS" << coord << ": " << ls_arr(coord) << std::endl;
                }
            });
        }

        solver->ls_mf[lev]->ParallelCopy(LS, 0, ils, 1, 0, 0);
    }
};

/*template <>
class EulerForwardScheme<Integrator::NarrowBandLevelset> : public TimeSteppingScheme<Integrator::NarrowBandLevelset> {
public:
    void ComputeSubStep(int lev, Set::Scalar dt, int stage, Integrator::NarrowBandLevelset* solver) override {
        if (stage != 0) {
            Util::Abort(__FILE__, __func__, __LINE__, "Invalid stage for EulerForward: only one stage exists.");
        }
        // Define constants
        const Set::Scalar* DX = solver->geom[lev].CellSize();  // Grid spacing
        const int ncomp = solver->number_of_components;

        // Define outside of narrowband
        const int OutsideNarrowband = solver->GetOutsideNarrowBandPos();

        // Loop through all components
        for (int ils=0; ils < ncomp; ils++){
            // Create all temporary multifabs here so no need to perform multiple copies
            auto const& ls_data = solver->level_sets[ils];

            // Get proper distribution mapping for narrowband
            const amrex::BoxArray& narrowband_ba = ls_data.narrowband_ba;
            const amrex::DistributionMapping& narrowband_dm = ls_data.narrowband_dm;

            // Define cell centered multifabs from solver
            amrex::iMultiFab narrowband(narrowband_ba, narrowband_dm, 1, 0);
            amrex::MultiFab LS(narrowband_ba, narrowband_dm, 1, 0);
            amrex::MultiFab LS_old(narrowband_ba, narrowband_dm, 1, 0);

            // ParallelCopy
            auto const& Tube = ls_data.Tube_imf;
            narrowband.ParallelCopy(*Tube, 0, 0, 1, 0, 0);
            LS.ParallelCopy(*solver->ls_mf[lev], ils, 0, 1, 0, 0);
            LS_old.ParallelCopy(*solver->ls_old_mf[lev], ils, 0, 1, 0, 0);

            // Define face-centered multifabs from solver
            amrex::IntVect face_type_x = amrex::IntVect::TheCellVector();  // All 0s: cell-centered
            face_type_x[0] = 1;  // Set x-direction to node-centered (i.e., face-centered in x)

            amrex::BoxArray xflux_ba(narrowband_ba);
            xflux_ba.convert(amrex::IndexType(face_type_x));
            amrex::MultiFab xflux(xflux_ba, narrowband_dm, 1, 1);
            xflux.ParallelCopy(*solver->XFlux_mf[lev], ils, 0, 1, 1, 1);

#if AMREX_SPACEDIM >= 2
            amrex::IntVect face_type_y = amrex::IntVect::TheCellVector();
            face_type_y[1] = 1;

            amrex::BoxArray yflux_ba(narrowband_ba);
            yflux_ba.convert(amrex::IndexType(face_type_y));
            amrex::MultiFab yflux(yflux_ba, narrowband_dm, 1, 1);
            yflux.ParallelCopy(*solver->YFlux_mf[lev], ils, 0, 1, 1, 1);
#endif

#if AMREX_SPACEDIM == 3
            amrex::IntVect face_type_z = amrex::IntVect::TheCellVector();
            face_type_z[2] = 1;

            amrex::BoxArray zflux_ba(narrowband_ba);
            zflux_ba.convert(amrex::IndexType(face_type_z));
            amrex::MultiFab zflux(zflux_ba, narrowband_dm, 1, 1);
            zflux.ParallelCopy(*solver->ZFlux_mf[lev], ils, 0, 1, 1, 1);
#endif

            // Loop through all valid boxes
            for (amrex::MFIter mfi(LS, false); mfi.isValid(); ++mfi) {
                // Define the valid box
                const amrex::Box& bx = mfi.validbox();  // Valid region of the box

                // Define array views for cell center variables
                auto const& nb_arr = narrowband.const_array(mfi);
                auto const& ls_arr = LS.array(mfi);
                auto const& ls_old_arr = LS_old.const_array(mfi);

                // Get face-centered flux fields
                auto const& xflux_arr = xflux.const_array(mfi);  // X-direction flux
#if AMREX_SPACEDIM >= 2
                auto const& yflux_arr = yflux.const_array(mfi);  // Y-direction flux
#endif
#if AMREX_SPACEDIM == 3
                auto const& zflux_arr = zflux.const_array(mfi);  // Z-direction flux
#endif

                // Loop through individual cells
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::IntVect coord(AMREX_D_DECL(i, j, k));

                    // Exit early if not narrowband cell
                    Set::Scalar nb_val = nb_arr(coord);
                    if (std::abs(nb_val) == OutsideNarrowband) return;

                    // Initialize flux to 0
                    Set::Scalar flux_div = 0.0;

                    // Compute flux divergence in X-direction
                    flux_div += (xflux_arr(coord, 0) - xflux_arr(i - 1, j, k, 0)) / DX[0];

#if AMREX_SPACEDIM >= 2
                    // Compute flux divergence in Y-direction
                    flux_div += (yflux_arr(coord, 0) - yflux_arr(i, j - 1, k, 0)) / DX[1];
#endif

#if AMREX_SPACEDIM == 3
                    // Compute flux divergence in Z-direction
                    flux_div += (zflux_arr(i, j, k, 0) - zflux_arr(i, j, k - 1, 0)) / DX[2];
#endif

                    // Euler forward update for level set function
                    //if (coord == debug) amrex::Print() << "LS: " << ls_arr(coord, ils) << std::endl;
                    ls_arr(coord, 0) = ls_old_arr(coord, 0) - dt * flux_div;
                });
            }

            // Parallel copy back
            solver->ls_mf[lev]->ParallelCopy(LS, 0, ils, 1, 0, 0);
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
        static constexpr Set::Scalar a[] = {1.0, 0.75, 1.0 / 3.0};  // Coefficients
        static constexpr Set::Scalar b[] = {0.0, 0.25, 2.0 / 3.0};
        static constexpr Set::Scalar c[] = {1.0, 0.25, 2.0 / 3.0};
    };

    void ComputeSubStep(int lev, Set::Scalar dt, int stage, IntegratorDerivedClass* solver) override {
        if (stage >= RK3Constants::stages) {
            Util::Abort(__FILE__, __func__, __LINE__, "Invalid stage for RK3.");
        }

        const Set::Scalar* DX = solver->geom[lev].CellSize();  // Grid spacing

        //for (amrex::MFIter mfi(*solver->QVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        for (amrex::MFIter mfi(*solver->QVec_mf[lev], false); mfi.isValid(); ++mfi) {
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
#if (AMREX_SPACEDIM == 3)
                    Set::Scalar zflux_front = 0.0, zflux_back = 0.0;
#endif

                    // X-direction flux divergence
                        xflux_right = xflux_arr(i, j, k, n);  // Right face flux
                        xflux_left = xflux_arr(i-1, j, k, n);  // Left face flux
                        flux_div += (xflux_right - xflux_left) / DX[0];

                    // Y-direction flux divergence (if 2D or 3D)
                        yflux_top = yflux_arr(i, j, k, n);  // Top face flux
                        yflux_bottom = yflux_arr(i, j-1, k, n);  // Bottom face flux
                        flux_div += (yflux_top - yflux_bottom) / DX[1];

                    // Z-direction flux divergence (if 3D)
#if (AMREX_SPACEDIM == 3)
                        zflux_front = zflux_arr(i, j, k, n);  // Front face flux
                        zflux_back = zflux_arr(i, j, k-1, n);  // Back face flux
                        flux_div += (zflux_front - zflux_back) / DX[2];
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

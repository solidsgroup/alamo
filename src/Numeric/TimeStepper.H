#ifndef NUMERIC_TIMESTEPPER_H_
#define NUMERIC_TIMESTEPPER_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Box.H>
#include <AMReX_MFIter.H>
#include <AMReX_Array4.H>
#include "Set/Set.H"
#include "Integrator/Integrator.H"
#include "Numeric/NumericFactory.H"
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"
#include <memory>
#include <string>

namespace Numeric {

/*// Helper function to compute face indices for a control volume's flux calculation
AMREX_GPU_DEVICE
inline void GetFaceIndices(int i, int j, int k, amrex::IntVect& west_face, amrex::IntVect& east_face,
                            amrex::IntVect& south_face, amrex::IntVect& north_face,
                            amrex::IntVect& bottom_face, amrex::IntVect& top_face) {
    // West and east face indices for the x-direction
    west_face = amrex::IntVect(AMREX_D_DECL(i - 1, j, k));
    east_face = amrex::IntVect(AMREX_D_DECL(i, j, k));

    // South and north face indices for the y-direction
    south_face = amrex::IntVect(AMREX_D_DECL(i, j - 1, k));
    north_face = amrex::IntVect(AMREX_D_DECL(i, j, k));

#if (AMREX_SPACEDIM == 3)
    // Bottom and top face indices for the z-direction
    bottom_face = amrex::IntVect(AMREX_D_DECL(i, j, k - 1));
    top_face = amrex::IntVect(AMREX_D_DECL(i, j, k));
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

        Util::ScimitarX_Util::Debug debug;
        debug.SetTargetDebugLocationIndices(32, 7, 7);
        //debug.SetTargetDebugLocation(0.49, 0.12, 0.0, solver->geom[lev]);

        const Set::Scalar* DX = solver->geom[lev].CellSize();  // Grid spacing
        

        //for (amrex::MFIter mfi(*solver->QVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        for (amrex::MFIter mfi(*solver->QVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();  // Valid region of the box

            auto const& q_arr = solver->QVec_mf.Patch(lev, mfi);  // Updated conserved variables
            auto const& q_old_arr = solver->QVec_old_mf.Patch(lev, mfi);  // Old conserved variables
            amrex::MultiFab const& XFLUX = *solver->XFlux_mf[lev];
            amrex::MultiFab const& YFLUX = *solver->YFlux_mf[lev];
            #if (AMREX_SPACEDIM == 3)
            amrex::MultiFab const& ZFLUX = *solver->ZFlux_mf[lev];
            #endif

            auto const& xflux_arr = XFLUX.const_array(mfi);
            auto const& yflux_arr = YFLUX.const_array(mfi);
            #if (AMREX_SPACEDIM == 3)
            auto const& zflux_arr = ZFLUX.const_array(mfi);
            #endif


            
//             auto const& xflux_arr = solver->XFlux_mf.Patch(lev, mfi);  // X-direction flux
//             auto const& yflux_arr = solver->YFlux_mf.Patch(lev, mfi);  // Y-direction flux
// #if (AMREX_SPACEDIM == 3)
//             auto const& zflux_arr = solver->ZFlux_mf.Patch(lev, mfi);  // Z-direction flux
// #endif
            auto const& sourceTerm_arr = solver->SourceTermVec.Patch(lev, mfi);  // Source-term derivatives from Viscous stress

            // Parallel loop over cells in the valid box
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                for (int n = 0; n < solver->number_of_components; ++n) {
                    // **Get Face Indices for the Control Volume**
                    //amrex::IntVect west_face, east_face, south_face, north_face, bottom_face, top_face;
                    //GetFaceIndices(i, j, k, west_face, east_face, south_face, north_face, bottom_face, top_face);

                    // **Flux Divergence Computation**
                    Set::Scalar flux_div = 0.0;

                    // X-direction flux divergence
                    flux_div += (xflux_arr(i+1, j, k, n) - xflux_arr(i, j, k, n)) / DX[0];

                    // Y-direction flux divergence
                    flux_div += (yflux_arr(i, j+1, k, n) - yflux_arr(i, j, k, n)) / DX[1];

#if (AMREX_SPACEDIM == 3)
                    // Z-direction flux divergence
                    flux_div += (zflux_arr(i, j, k+1, n) - zflux_arr(i, j, k, n)) / DX[2];
#endif

                    // Update conserved variable using Euler Forward step
                    Set::Scalar q_old = q_old_arr(i, j, k, n);
                    q_arr(i, j, k, n) = q_old - dt * (flux_div - sourceTerm_arr(i, j, k, n));

                    // //Debug output
                    // debug.DebugEulerForwardStep(i, j, k, xflux_arr(i+1,j,k,n), xflux_arr(i,j,k,n),
                    //                             yflux_arr(i,j+1,k,n), yflux_arr(i,j,k,n),
                    //                             0.0,
                    //                             0.0,
                    //                             q_old, dt, "Euler Forward Substep", true); 

                                    // Only print at target cell
// if (Util::ScimitarX_Util::Debug::IsTargetLocation(i,j,k))
// {
//     // ===================== comp 0 =====================
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, xflux_arr(i+1,j,k,0), "XFlux(0)@i+1",   "EULER: flux faces", false, 0, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, xflux_arr(i  ,j,k,0), "XFlux(0)@i", "EULER: flux faces", false, 0, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, yflux_arr(i,j+1,k,0), "YFlux(0)@j+1",   "EULER: flux faces", false, 0, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, yflux_arr(i,j  ,k,0), "YFlux(0)@j", "EULER: flux faces", false, 0, true);
// #if (AMREX_SPACEDIM == 3)
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, zflux_arr(i,j,k+1,0), "ZFlux(0)@k+1",   "EULER: flux faces", false, 0, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, zflux_arr(i,j,k  ,0), "ZFlux(0)@k", "EULER: flux faces", false, 0, true);
// #endif
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, flux_div, "flux_div(0)", "EULER: divergence", false, 0, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, q_arr(i,j,k,0), "Qnew(0)", "EULER: update", false, 0, true);

//     // ===================== comp 1 =====================
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, xflux_arr(i+1,j,k,1), "XFlux(1)@i+1",   "EULER: flux faces", false, 1, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, xflux_arr(i  ,j,k,1), "XFlux(1)@i", "EULER: flux faces", false, 1, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, yflux_arr(i,j+1,k,1), "YFlux(1)@j+1",   "EULER: flux faces", false, 1, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, yflux_arr(i,j  ,k,1), "YFlux(1)@j", "EULER: flux faces", false, 1, true);
// #if (AMREX_SPACEDIM == 3)
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, zflux_arr(i,j,k+1,1), "ZFlux(1)@k+1",   "EULER: flux faces", false, 1, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, zflux_arr(i,j,k  ,1), "ZFlux(1)@k", "EULER: flux faces", false, 1, true);
// #endif
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, flux_div, "flux_div(1)", "EULER: divergence", false, 1, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, q_arr(i,j,k,1), "Qnew(1)", "EULER: update", false, 1, true);

//     // ===================== comp 2 =====================
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, xflux_arr(i+1,j,k,2), "XFlux(2)@i+1",   "EULER: flux faces", false, 2, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, xflux_arr(i  ,j,k,2), "XFlux(2)@i", "EULER: flux faces", false, 2, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, yflux_arr(i,j+1,k,2), "YFlux(2)@j+1",   "EULER: flux faces", false, 2, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, yflux_arr(i,j  ,k,2), "YFlux(2)@j", "EULER: flux faces", false, 2, true);
// #if (AMREX_SPACEDIM == 3)
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, zflux_arr(i,j,k+1,2), "ZFlux(2)@k+1",   "EULER: flux faces", false, 2, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, zflux_arr(i,j,k  ,2), "ZFlux(2)@k", "EULER: flux faces", false, 2, true);
// #endif
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, flux_div, "flux_div(2)", "EULER: divergence", false, 2, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, q_arr(i,j,k,2), "Qnew(2)", "EULER: update", false, 2, true);

//     // ===================== comp 3 =====================
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, xflux_arr(i+1,j,k,3), "XFlux(3)@i+1",   "EULER: flux faces", false, 3, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, xflux_arr(i  ,j,k,3), "XFlux(3)@i", "EULER: flux faces", false, 3, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, yflux_arr(i,j+1,k,3), "YFlux(3)@j+1",   "EULER: flux faces", false, 3, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, yflux_arr(i,j  ,k,3), "YFlux(3)@j", "EULER: flux faces", false, 3, true);
// #if (AMREX_SPACEDIM == 3)
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, zflux_arr(i,j,k+1,3), "ZFlux(3)@k+1",   "EULER: flux faces", false, 3, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, zflux_arr(i,j,k  ,3), "ZFlux(3)@k", "EULER: flux faces", false, 3, true);
// #endif
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, flux_div, "flux_div(3)", "EULER: divergence", false, 3, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, q_arr(i,j,k,3), "Qnew(3)", "EULER: update", false, 3, true);

//     // ===================== comp 4 =====================
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, xflux_arr(i+1,j,k,4), "XFlux(4)@i+1",   "EULER: flux faces", false, 4, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, xflux_arr(i  ,j,k,4), "XFlux(4)@i", "EULER: flux faces", false, 4, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, yflux_arr(i,j+1,k,4), "YFlux(4)@j+1",   "EULER: flux faces", false, 4, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, yflux_arr(i,j  ,k,4), "YFlux(4)@j", "EULER: flux faces", false, 4, true);
// #if (AMREX_SPACEDIM == 3)
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, zflux_arr(i,j,k+1,4), "ZFlux(4)@k+1",   "EULER: flux faces", false, 4, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, zflux_arr(i,j,k  ,4), "ZFlux(4)@k", "EULER: flux faces", false, 4, true);
// #endif
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, flux_div, "flux_div(4)", "EULER: divergence", false, 4, true);
//     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
//         i,j,k, q_arr(i,j,k,4), "Qnew(4)", "EULER: update", false, 4, true);
// }

                }
            });
        }
    }

    int GetNumberOfStages() const override {
        return 1;  // Euler Forward has only one stage
    }
};

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
            auto const& sourceTerm_arr = solver->SourceTermVec.Patch(lev, mfi);  // Source-term derivatives from Viscous stress

            {
                const int ic=32, jc=7, kc=7;
                const int nprobe = 3;  

                const amrex::IntVect iv(ic,jc,kc);
                const amrex::Box& tbx = mfi.tilebox();

                
                if (tbx.contains(iv)) {
                
                    const auto FxR = xflux_arr(ic+1, jc, kc, nprobe);
                    const auto FxL = xflux_arr(ic  , jc, kc, nprobe);

                    const auto FyT = yflux_arr(ic, jc+1, kc, nprobe);
                    const auto FyB = yflux_arr(ic, jc  , kc, nprobe);

            #if (AMREX_SPACEDIM == 3)
                    const auto FzF = zflux_arr(ic, jc, kc+1, nprobe);
                    const auto FzB = zflux_arr(ic, jc, kc  , nprobe);
            #endif

                    amrex::AllPrint()
                        << "\n[TS PROBE pre] proc=" << amrex::ParallelDescriptor::MyProc()
                        << " stage=" << stage
                        << " tile=" << tbx
                        << " cell=("<<ic<<","<<jc<<","<<kc<<") n="<<nprobe
                        << "\n  q_old=" << q_old_arr(ic,jc,kc,nprobe)
                        << "\n  q_cur=" << q_arr(ic,jc,kc,nprobe)
                        << "\n  FxR=" << FxR << " FxL=" << FxL << " dFx=" << (FxR-FxL)
                        << "\n  FyT=" << FyT << " FyB=" << FyB << " dFy=" << (FyT-FyB)
            #if (AMREX_SPACEDIM == 3)
                        << "\n  FzF=" << FzF << " FzB=" << FzB << " dFz=" << (FzF-FzB)
            #endif
                        << "\n";
                }
            }
                  

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
                        xflux_right = xflux_arr(i+1, j, k, n);  // Right face flux
                        xflux_left  = xflux_arr(i  , j, k, n);  // Left face flux
                        flux_div += (xflux_right - xflux_left) / DX[0];

                    // Y-direction flux divergence (if 2D or 3D)
                        yflux_top    = yflux_arr(i, j+1, k, n);  // Top face flux
                        yflux_bottom = yflux_arr(i, j  , k, n);  // Bottom face flux
                        flux_div += (yflux_top - yflux_bottom) / DX[1];

                        

                    // Z-direction flux divergence (if 3D)
#if (AMREX_SPACEDIM == 3)
                        zflux_front = zflux_arr(i, j, k+1, n);  // Front face flux
                        zflux_back  = zflux_arr(i, j, k  , n);  // Back face flux
                        flux_div += (zflux_front - zflux_back) / DX[2];
#endif

                    // **Runge-Kutta 3rd order update formula**
                    q_arr(i, j, k, n) = RK3Constants::a[stage] * q_old_arr(i, j, k, n) +
                                        RK3Constants::b[stage] * q_arr(i, j, k, n) -
                                        RK3Constants::c[stage] * dt * (flux_div - sourceTerm_arr(i, j, k, n));
                                        

                                        // Update conserved variable
                    
                    // ---- PROBE: after timestepper update (device-safe) ----
                    const int rho1 = solver->variableIndex.DENS;

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        i,j,k, q_old_arr(i,j,k,rho1),
                        "Qold(rho)", "TS: before update", false, rho1, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        i,j,k, q_arr(i,j,k,rho1),
                        "Qnew(rho)", "TS: after update", false, rho1, true);
                    

                    // Util::ScimitarX_Util::Debug::SetTargetDebugLocationIndices(32, 7, 7, true);

                    
                    // if (Util::ScimitarX_Util::Debug::IsTargetLocation(i,j,k) && n == rho)
                    // {
                    //     // raw fluxes at the faces used by this cell
                    //     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                    //         i,j,k, xflux_arr(i,  j,k,n), "XFlux(i)",   "RK3: flux faces", false, n, true);
                    //     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                    //         i,j,k, xflux_arr(i-1,j,k,n), "XFlux(i-1)", "RK3: flux faces", false, n, true);
                    
                    //     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                    //         i,j,k, yflux_arr(i,j,  k,n), "YFlux(j)",   "RK3: flux faces", false, n, true);
                    //     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                    //         i,j,k, yflux_arr(i,j-1,k,n), "YFlux(j-1)", "RK3: flux faces", false, n, true);
                    
                    // #if (AMREX_SPACEDIM == 3)
                    //     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                    //         i,j,k, zflux_arr(i,j,k,  n), "ZFlux(k)",   "RK3: flux faces", false, n, true);
                    //     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                    //         i,j,k, zflux_arr(i,j,k-1,n), "ZFlux(k-1)", "RK3: flux faces", false, n, true);
                    // #endif
                    
                    //     // the exact scalars you used in flux_div (nice sanity check)
                    //     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                    //         i,j,k, xflux_right, "xflux_right", "RK3: derived", false, n, true);
                    //     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                    //         i,j,k, xflux_left,  "xflux_left",  "RK3: derived", false, n, true);
                    
                    //     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                    //         i,j,k, yflux_top,    "yflux_top",    "RK3: derived", false, n, true);
                    //     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                    //         i,j,k, yflux_bottom, "yflux_bottom", "RK3: derived", false, n, true);
                    
                    // #if (AMREX_SPACEDIM == 3)
                    //     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                    //         i,j,k, zflux_front, "zflux_front", "RK3: derived", false, n, true);
                    //     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                    //         i,j,k, zflux_back,  "zflux_back",  "RK3: derived", false, n, true);
                    // #endif
                    
                    //     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                    //         i,j,k, flux_div, "flux_div", "RK3: flux div", false, n, true);
                    
                    //     Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                    //         i,j,k, sourceTerm_arr(i,j,k,n), "Src", "RK3: source", false, n, true);
                    // }
                    

                }
            });
            
        
            const int rho = solver->variableIndex.DENS;
            const amrex::IntVect p7(32, 7, 7);
            const amrex::IntVect p8(32, 8, 7);
            const amrex::Box& fabx = (*solver->QVec_mf[lev])[mfi].box();

            if (fabx.contains(p7) || fabx.contains(p8)) {

                auto print_q = [&](const amrex::IntVect& p, const char* tag) {
                    if (!fabx.contains(p)) return;

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        32,7,7, q_arr(p[0],p[1],p[2],rho),
                        std::string("q_new(rho)@") + tag, "TS: POST-UPDATE", false, rho, true);
                };

                print_q(p7, "j7");
                print_q(p8, "j8");
}

            
        }
            
        
        {   
            
            const int ic1 = 32, jc1 = 7, kc1 = 7;
            const int ic2 = 32, jc2 = 8, kc2 = 7;
            const int rho = solver->variableIndex.DENS;

        
            for (amrex::MFIter mfi(*solver->QVec_mf[lev]); mfi.isValid(); ++mfi) {

                const amrex::Box& fbx = (*solver->QVec_mf[lev])[mfi].box();
                if (fbx.contains(amrex::IntVect(ic1,jc1,kc1)) ||
                    fbx.contains(amrex::IntVect(ic2,jc2,kc2))) {

                    auto const& q = (*solver->QVec_mf[lev]).const_array(mfi);

                    if (fbx.contains(amrex::IntVect(ic1,jc1,kc1))) {
                        amrex::AllPrint()
                            << "\n[TS SIMPLE] stage=" << stage
                            << " DENS(32,7,7)=" << q(ic1,jc1,kc1,rho);
                    }

                    if (fbx.contains(amrex::IntVect(ic2,jc2,kc2))) {
                        amrex::AllPrint()
                            << "  DENS(32,8,7)=" << q(ic2,jc2,kc2,rho)
                            << "\n";
                    }
                }
            }
        }
// ===== END SIMPLE PROBE =====

        // solver->XFlux_mf.FillBoundary();
        // solver->YFlux_mf.FillBoundary();
        // solver->ZFlux_mf.FillBoundary();
    }
    
    int GetNumberOfStages() const override {
        return RK3Constants::stages;  // RK3 has 3 stages
    }
};

}  // namespace Numeric

#endif  // NUMERIC_TIMESTEPPER_H_

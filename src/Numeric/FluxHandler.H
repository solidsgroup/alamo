#ifndef NUMERIC_FLUXHANDLER_H
#define NUMERIC_FLUXHANDLER_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Set/Set.H"
#include "Integrator/ScimitarX.H"
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"
#include "Model/Fluid/Fluid.H"  // Include fluid model helper functions
#include "Numeric/Stencil.H"                                
#include <memory>
#include <cmath>

namespace Numeric {

// Base class for flux reconstruction
template <typename IntegratorDerivedClass>
class FluxReconstruction {
public:
    virtual ~FluxReconstruction() = default;

    /**
     * @brief Perform flux reconstruction for the specified direction.
     */
    virtual void ReconstructFlux(int direction, int lev, IntegratorDerivedClass* solver, amrex::MultiFab& QL_stencil, amrex::MultiFab& QR_stencil) const = 0;
};

// Base class for flux method (e.g., Local Lax-Friedrichs, HLLC)
template <typename IntegratorDerivedClass>
class FluxMethod {
public:
    virtual ~FluxMethod() = default;

    /**
     * @brief Compute the flux for the specified direction.
     */
    virtual void ComputeFlux(int direction, int lev, IntegratorDerivedClass* solver, const amrex::MultiFab& QL_stencil, const amrex::MultiFab& QR_stencil) const = 0;
};

// FluxHandler class for managing reconstruction and flux computation
template <typename IntegratorDerivedClass>
class FluxHandler {
public:
    enum Directions { Xdir = 0, Ydir = 1, Zdir = 2 };

    FluxHandler() {}

    void SetReconstruction(std::shared_ptr<FluxReconstruction<IntegratorDerivedClass>> recon) {
        reconstruction = std::move(recon);
    }

    void SetFluxMethod(std::shared_ptr<FluxMethod<IntegratorDerivedClass>> method) {
        flux_method = std::move(method);
    }

    void ConstructFluxes(int lev, IntegratorDerivedClass* solver) const {
        if (!reconstruction || !flux_method) {
            amrex::Abort("FluxHandler: Reconstruction or flux method is not set.");
        }

        Util::Message(INFO, "Constructing fluxes for level " + std::to_string(lev));

        // Create temporary QL and QR stencils for the level
        amrex::MultiFab QL_stencil(solver->grids[lev], solver->dmap[lev], solver->number_of_components, solver->number_of_ghost_cells);
        amrex::MultiFab QR_stencil(solver->grids[lev], solver->dmap[lev], solver->number_of_components, solver->number_of_ghost_cells);

        reconstruction->ReconstructFlux(Xdir, lev, solver, QL_stencil, QR_stencil);
        flux_method->ComputeFlux(Xdir, lev, solver, QL_stencil, QR_stencil);

#if (AMREX_SPACEDIM >= 2)
        reconstruction->ReconstructFlux(Ydir, lev, solver, QL_stencil, QR_stencil);
        flux_method->ComputeFlux(Ydir, lev, solver, QL_stencil, QR_stencil);
#endif
#if (AMREX_SPACEDIM == 3)
        reconstruction->ReconstructFlux(Zdir, lev, solver, QL_stencil, QR_stencil);
        flux_method->ComputeFlux(Zdir, lev, solver, QL_stencil, QR_stencil);
#endif
    }

private:
    std::shared_ptr<FluxReconstruction<IntegratorDerivedClass>> reconstruction;  // Reconstruction method
    std::shared_ptr<FluxMethod<IntegratorDerivedClass>> flux_method;             // Flux method

};

// First-order upwind reconstruction
template <typename IntegratorDerivedClass>
class FirstOrderReconstruction : public FluxReconstruction<IntegratorDerivedClass> {
public:
    void ReconstructFlux(int direction, int lev, IntegratorDerivedClass* solver, amrex::MultiFab& QL_stencil, amrex::MultiFab& QR_stencil) const override {
        Util::Message(INFO, "Performing first-order upwind flux reconstruction for direction: " + std::to_string(direction));
        Util::ScimitarX_Util::Debug debug;

        Model::Fluid::Fluid fluid_model;

        debug.SetTargetDebugLocationIndices(50, 0, 0);  // Debug target location

        const amrex::Box& domain_box = solver->geom[lev].Domain();  // Get physical domain for boundary detection

        // **Direction to velocity component mapping**
        int velocity_component[3] = {
            solver->variableIndex.UVEL,  // X-direction -> UVEL
            solver->variableIndex.VVEL,  // Y-direction -> VVEL
            solver->variableIndex.WVEL   // Z-direction -> WVEL
        };



        for (amrex::MFIter mfi(*solver->QVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();  // Get the valid computational box for this patch
            auto const& QL_arr = QL_stencil.array(mfi);  // Left reconstructed state
            auto const& QR_arr = QR_stencil.array(mfi);  // Right reconstructed state
            auto const& p_arr = solver->PVec_mf.Patch(lev, mfi);  // Access primitive variables
            auto const& pressure_arr = solver->Pressure_mf.Patch(lev, mfi);  // Access pressure field

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    
                Set::Scalar vel[3] = {0.0, 0.0, 0.0};  

                // **Primitive variables**
                Set::Scalar rho = p_arr(i, j, k, solver->variableIndex.DENS);  // Density
                            vel[0] = p_arr(i, j, k, velocity_component[direction]);  // Primary velocity
                            vel[1] = (AMREX_SPACEDIM >= 2) ? p_arr(i, j, k, velocity_component[(direction + 1) % AMREX_SPACEDIM]) : 0.0;  // First transverse velocity
                            vel[2] = (AMREX_SPACEDIM == 3) ? p_arr(i, j, k, velocity_component[(direction + 2) % AMREX_SPACEDIM]) : 0.0;  // Second transverse velocity
                Set::Scalar p = pressure_arr(i, j, k);  // Pressure
                Set::Scalar gamma = 1.4;  // Ratio of specific heats

                // **Compute speed of sound and characteristic wave speed**
                Set::Scalar c = fluid_model.ComputeWaveSpeed(rho, p, gamma);  // Speed of sound
                Set::Scalar lambda = (vel[0] >= 0) ? (vel[0] + c) : (vel[0] - c);  // Characteristic velocity
                int wave_sign = (lambda > 0) ? 1 : -1;  // Determine wave direction

                // Determine offsets for left and right states
                int offset_i = (direction == FluxHandler<IntegratorDerivedClass>::Xdir) ? wave_sign : 0;
                int offset_j = (direction == FluxHandler<IntegratorDerivedClass>::Ydir) ? wave_sign : 0;
                int offset_k = (direction == FluxHandler<IntegratorDerivedClass>::Zdir) ? wave_sign : 0;

                // **Get stencil type (Lo, Hi, Central) for each direction**
                std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil = Numeric::GetStencil(i, j, k, domain_box);

                debug.DebugValuesIfTarget(i, j, k, rho, "Density", "Reconstruction Step", false);

                for (int n = 0; n < solver->number_of_components; ++n) {
                    if (stencil[direction] == Numeric::StencilType::Lo) {
                        // Near lower boundary (backward difference for left state)
                        QL_arr(i, j, k, n) = p_arr(i, j, k, n);  // Zero-gradient (Neumann) or inflow value
                        QR_arr(i, j, k, n) = p_arr(i + offset_i, j + offset_j, k + offset_k, n);  // Use neighbor on the "right"
                    } else if (stencil[direction] == Numeric::StencilType::Hi) {
                        // Near upper boundary (forward difference for right state)
                        QL_arr(i, j, k, n) = p_arr(i - offset_i, j - offset_j, k - offset_k, n);  // Use neighbor on the "left"
                        QR_arr(i, j, k, n) = p_arr(i, j, k, n);  // Zero-gradient (Neumann) or outflow value
                    } else {
                        // Interior point (Upwind reconstruction based on wave direction)
                        QL_arr(i, j, k, n) = (wave_sign > 0) ? p_arr(i, j, k, n) : p_arr(i - offset_i, j - offset_j, k - offset_k, n);
                        QR_arr(i, j, k, n) = (wave_sign > 0) ? p_arr(i + offset_i, j + offset_j, k + offset_k, n) : p_arr(i, j, k, n);
                    }

                    // **Debug output for QL and QR arrays**
                    debug.DebugValuesIfTarget(i, j, k, QL_arr(i, j, k, n), "QL_arr", "Upwind QL", false, n);
                    debug.DebugValuesIfTarget(i, j, k, QR_arr(i, j, k, n), "QR_arr", "Upwind QR", false, n);
                }
            });
        }
    }
};


// Local Lax-Friedrichs flux method
template <typename IntegratorDerivedClass>
class LocalLaxFriedrichsMethod : public FluxMethod<IntegratorDerivedClass> {
public:
    void ComputeFlux(int direction, int lev, IntegratorDerivedClass* solver, const amrex::MultiFab& QL_stencil, const amrex::MultiFab& QR_stencil) const override {
        Util::Message(INFO, "Computing Local Lax-Friedrichs flux for direction: " + std::to_string(direction));

        const Set::Scalar gamma = 1.4;  // Ratio of specific heats
        Util::ScimitarX_Util::Debug debug;
        debug.SetTargetDebugLocationIndices(50, 0, 0);  // Debug target location

        for (amrex::MFIter mfi(*solver->QVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();
            auto const& flux_arr = (direction == FluxHandler<IntegratorDerivedClass>::Xdir) ? solver->XFlux_mf.Patch(lev, mfi) :
                                   (direction == FluxHandler<IntegratorDerivedClass>::Ydir) ? solver->YFlux_mf.Patch(lev, mfi) :
                                   solver->ZFlux_mf.Patch(lev, mfi);
            auto const& QL_arr = QL_stencil.array(mfi);
            auto const& QR_arr = QR_stencil.array(mfi);

            // **Direction to velocity mapping**
            int velocity_component[3] = {
                solver->variableIndex.UVEL,  // X-direction -> UVEL
                solver->variableIndex.VVEL,  // Y-direction -> VVEL
                solver->variableIndex.WVEL   // Z-direction -> WVEL
            };


            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // **Left state (QL) primitive variables**
                Set::Scalar rho_L = QL_arr(i, j, k, solver->variableIndex.DENS);  // Density
                Set::Scalar vel_L[3] = {0.0, 0.0, 0.0};  // Velocity array for left state
                vel_L[0] = QL_arr(i, j, k, velocity_component[direction]);  // Primary velocity
                vel_L[1] = (AMREX_SPACEDIM >= 2) ? QL_arr(i, j, k, velocity_component[(direction + 1) % AMREX_SPACEDIM]) : 0.0;  // Transverse 1
                vel_L[2] = (AMREX_SPACEDIM == 3) ? QL_arr(i, j, k, velocity_component[(direction + 2) % AMREX_SPACEDIM]) : 0.0;  // Transverse 2
                Set::Scalar e_L = QL_arr(i, j, k, solver->variableIndex.IE);  // Specific internal energy
                Set::Scalar KE_L = 0.5 * (vel_L[0] * vel_L[0] + vel_L[1] * vel_L[1] + vel_L[2] * vel_L[2]);  // Kinetic energy
                Set::Scalar E_L = e_L + KE_L;  // Total energy
                Set::Scalar pressure_L = (gamma - 1.0) * (e_L);  // Pressure

                // **Right state (QR) primitive variables**
                Set::Scalar rho_R = QR_arr(i, j, k, solver->variableIndex.DENS);  // Density
                Set::Scalar vel_R[3] = {0.0, 0.0, 0.0};  // Velocity array for right state
                vel_R[0] = QR_arr(i, j, k, velocity_component[direction]);  // Primary velocity
                vel_R[1] = (AMREX_SPACEDIM >= 2) ? QR_arr(i, j, k, velocity_component[(direction + 1) % AMREX_SPACEDIM]) : 0.0;  // Transverse 1
                vel_R[2] = (AMREX_SPACEDIM == 3) ? QR_arr(i, j, k, velocity_component[(direction + 2) % AMREX_SPACEDIM]) : 0.0;  // Transverse 2
                Set::Scalar e_R = QR_arr(i, j, k, solver->variableIndex.IE);  // Specific internal energy
                Set::Scalar KE_R = 0.5 * (vel_R[0] * vel_R[0] + vel_R[1] * vel_R[1] + vel_R[2] * vel_R[2]);  // Kinetic energy
                Set::Scalar E_R = e_R + KE_R;  // Total energy
                Set::Scalar pressure_R = (gamma - 1.0) * (e_R);  // Pressure

                // **Speed of sound**
                Set::Scalar c_L = std::sqrt(gamma * std::max(pressure_L, 1e-6) / rho_L);  // Left speed of sound
                Set::Scalar c_R = std::sqrt(gamma * std::max(pressure_R, 1e-6) / rho_R);  // Right speed of sound

                // **Maximum wave speed (lambda_max)**
                Set::Scalar lambda_max = std::max(std::abs(vel_L[0]) + c_L, std::abs(vel_R[0]) + c_R);

                // **Fluxes (F_L and F_R)**
                Set::Scalar F_L[5], F_R[5];
                F_L[0] = rho_L * vel_L[0];                                 // Mass flux
                F_R[0] = rho_R * vel_R[0];
                F_L[1] = rho_L * vel_L[0] * vel_L[0] + pressure_L;         // Primary momentum flux
                F_R[1] = rho_R * vel_R[0] * vel_R[0] + pressure_R;
                F_L[2] = rho_L * vel_L[0] * vel_L[1];                      // Cross-momentum flux 1
                F_R[2] = rho_R * vel_R[0] * vel_R[1];
                F_L[3] = (AMREX_SPACEDIM == 3) ? rho_L * vel_L[0] * vel_L[2] : 0.0;  // Cross-momentum flux 2 (3D only)
                F_R[3] = (AMREX_SPACEDIM == 3) ? rho_R * vel_R[0] * vel_R[2] : 0.0;
                F_L[4] = (rho_L * E_L + pressure_L) * vel_L[0];            // Energy flux
                F_R[4] = (rho_R * E_R + pressure_R) * vel_R[0];

                // **Local Lax-Friedrichs flux calculation**
                for (int n = 0; n < solver->number_of_components; ++n) {
                    flux_arr(i, j, k, n) = 0.5 * (F_L[n] + F_R[n]) - 0.5 * lambda_max * (QR_arr(i, j, k, n) - QL_arr(i, j, k, n));
                }

                // Debug
                debug.DebugLocalLaxFriedrichsFlux(i, j, k,
                                                  rho_L, vel_L[0], vel_L[1], vel_L[2], e_L, pressure_L,
                                                  rho_R, vel_R[0], vel_R[1], vel_R[2], e_R, pressure_R,
                                                  c_L, c_R, lambda_max,
                                                  F_L, F_R,
                                                  "ComputeFlux", true);
                debug.CheckNaNAndAbort(i, j, k, lambda_max, "Lambda Max", "wave speed", true);
            });
        }

        Util::Message(INFO, "Local Lax-Friedrichs flux computed for direction: " + std::to_string(direction));
    }
};

}  // namespace Numeric

#endif  // NUMERIC_FLUXHANDLER_H


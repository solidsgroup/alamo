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
#include <algorithm>

namespace Numeric {

// ClampIndices: Clamps indices to valid bounds (branch-free logic)
AMREX_GPU_DEVICE
inline void ClampIndices(int indices[3], const int lower_bounds[3], const int upper_bounds[3]) {
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        indices[d] = lower_bounds[d] + std::min(std::max(indices[d] - lower_bounds[d], 0), upper_bounds[d] - lower_bounds[d]);
    }
}

// ShiftAndClampIndices: Shifts and clamps indices (branch-free logic)
AMREX_GPU_DEVICE
inline void ShiftAndClampIndices(int indices[3], const int offsets[3],
                                 const int lower_bounds[3], const int upper_bounds[3], int direction) {
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        int shifted_index = indices[d] + (d == direction ? offsets[d] : 0);
        indices[d] = lower_bounds[d] + std::min(std::max(shifted_index - lower_bounds[d], 0), upper_bounds[d] - lower_bounds[d]);
    }
}

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

      //  Util::Message(INFO, "Constructing fluxes for level " + std::to_string(lev));

   for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        amrex::BoxArray face_ba = solver->grids[lev];
        amrex::IntVect face_type = amrex::IntVect::TheZeroVector();  // `{0, 0, 0}`
        face_type[dir] = 1;  // Make it face-centered in `dir`
        face_ba.convert(face_type);

        // Create temporary QL and QR MultiFabs at face-centered positions
        amrex::MultiFab QL_stencil(face_ba, solver->dmap[lev], solver->number_of_components, solver->number_of_ghost_cells);
        amrex::MultiFab QR_stencil(face_ba, solver->dmap[lev], solver->number_of_components, solver->number_of_ghost_cells);

        reconstruction->ReconstructFlux(dir, lev, solver, QL_stencil, QR_stencil);
        flux_method->ComputeFlux(dir, lev, solver, QL_stencil, QR_stencil);
    }

    }

private:
    std::shared_ptr<FluxReconstruction<IntegratorDerivedClass>> reconstruction;  // Reconstruction method
    std::shared_ptr<FluxMethod<IntegratorDerivedClass>> flux_method;             // Flux method

};

// First-order upwind reconstruction
template <typename IntegratorDerivedClass>
class FirstOrderReconstruction : public FluxReconstruction<IntegratorDerivedClass> {
public:
// ReconstructFlux: First-order upwind flux reconstruction with branch-free framework
void ReconstructFlux(int direction, int lev, IntegratorDerivedClass* solver, amrex::MultiFab& QL_stencil, amrex::MultiFab& QR_stencil) const override {
    //Util::Message(INFO, "Performing first-order upwind flux reconstruction for direction: " + std::to_string(direction));
    //Util::ScimitarX_Util::Debug debug;
    //debug.SetTargetDebugLocationIndices(101, 0, 0);  // Debug target location
    //debug.SetTargetDebugLocation(0.49, 0.12, 0.0, solver->geom[lev]);  // Debug target location
    const int nghosts = solver->number_of_ghost_cells; 
// Loop through face-centered indices for flux reconstruction
//for (amrex::MFIter mfi(QL_stencil, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
    const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);

    auto const& QL_arr = QL_stencil.array(mfi);  // Left reconstructed state
    auto const& QR_arr = QR_stencil.array(mfi);  // Right reconstructed state
    auto const& p_arr = solver->PVec_mf.Patch(lev, mfi);  // Conservative variables

    const int num_components = solver->number_of_components;  // Number of components (e.g., rho, u, v, E)

    // Define velocity components based on direction
    int velocity_component[3] = {
        solver->variableIndex.UVEL,  // X-direction -> UVEL
        solver->variableIndex.VVEL,  // Y-direction -> VVEL
        solver->variableIndex.WVEL   // Z-direction -> WVEL
    };

    amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {

   // Construct index arrays with correct number of components
    int index[3];
    int left_index[3];
    int right_index[3];
    int lower_bounds[3];
    int upper_bounds[3];

#if AMREX_SPACEDIM == 2
    index[0] = iface; index[1] = jface; index[2] = kface;
    left_index[0] = iface; left_index[1] = jface; left_index[2] = kface;
    right_index[0] = iface; right_index[1] = jface; right_index[2] = kface;

    lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); lower_bounds[2] = 0;
    upper_bounds[0] = face_bx_with_ghosts.bigEnd(0)-1; upper_bounds[1] = face_bx_with_ghosts.bigEnd(1)-1; upper_bounds[2] = 0;
#else
    index[0] = iface; index[1] = jface; index[2] = kface;
    left_index[0] = iface; left_index[1] = jface; left_index[2] = kface;
    right_index[0] = iface; right_index[1] = jface; right_index[2] = kface;

    lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); lower_bounds[2] = face_bx_with_ghosts.smallEnd(2);
    upper_bounds[0] = face_bx_with_ghosts.bigEnd(0)-1; upper_bounds[1] = face_bx_with_ghosts.bigEnd(1)-1; upper_bounds[2] = face_bx_with_ghosts.bigEnd(2)-1;
#endif

    // Clamp the main index
    ClampIndices(index, lower_bounds, upper_bounds);

    // Shift and clamp for left and right indices
    int left_offset[3] = {0, 0, 0};
    int right_offset[3] = {0, 0, 0};
    left_offset[direction] = -1;  // Shift left
    right_offset[direction] = 1;  // Shift right

    ShiftAndClampIndices(left_index, left_offset, lower_bounds, upper_bounds, direction);
    ShiftAndClampIndices(right_index, right_offset, lower_bounds, upper_bounds, direction);
    
    //debug.DebugIndicesAndBounds(iface, jface, kface, index, left_index, right_index, lower_bounds, upper_bounds, "QL and QR Construction");

        const Set::Scalar vel = p_arr(index[0], index[1], index[2], velocity_component[direction]);  // Normal velocity
        const Set::Scalar sigma = std::copysign(1.0, vel);  // 1 if vel > 0, 0 if vel <= 0

        for (int n = 0; n < num_components; ++n) {
            // Branch-free computation of QL and QR
            const Set::Scalar QL = 0.5 * ((1 + sigma)* p_arr(index[0], index[1], index[2], n) + (1.0 - sigma) * p_arr(left_index[0], left_index[1], left_index[2], n));
            const Set::Scalar QR = 0.5 * ((1 + sigma)* p_arr(right_index[0], right_index[1], right_index[2], n) + (1.0 - sigma)*p_arr(index[0], index[1], index[2],n));  // Since we are using first-order upwinding  
            
            // Assign values to the stencil arrays
            QL_arr(index[0], index[1], index[2], n) = QL;
            QR_arr(index[0], index[1], index[2], n) = QR;

            // Debugging output
            //debug.DebugValuesIfTarget(iface, jface, kface, QL, "QL_arr", "Reconstructed Left State", false, n);
            //debug.DebugValuesIfTarget(iface, jface, kface, QR, "QR_arr", "Reconstructed Right State", false, n);
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
        //Util::Message(INFO, "Computing Local Lax-Friedrichs flux for direction: " + std::to_string(direction));

        const Set::Scalar gamma = 1.4;  // Ratio of specific heats
        //Util::ScimitarX_Util::Debug debug;
        //debug.SetTargetDebugLocationIndices(101, 0, 0);  // Debug target location
        //debug.SetTargetDebugLocation(0.49, 0.12, 0.0,solver->geom[lev]);  // Debug target location

        const int nghosts = solver->number_of_ghost_cells;

        //for (amrex::MFIter mfi(QL_stencil, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {

            const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts); //Includes two ghost layers

            auto const& flux_arr = (direction == FluxHandler<IntegratorDerivedClass>::Xdir) ? solver->XFlux_mf.Patch(lev, mfi) :
                                   (direction == FluxHandler<IntegratorDerivedClass>::Ydir) ? solver->YFlux_mf.Patch(lev, mfi) :
                                   solver->ZFlux_mf.Patch(lev, mfi);
            auto const& QL_arr = QL_stencil.array(mfi);
            auto const& QR_arr = QR_stencil.array(mfi);

            // **Direction to velocity and energy mapping**
            int normal = direction;                 // Normal velocity component
            int trans1 = (direction + 1) % AMREX_SPACEDIM;  // First transverse velocity
#if AMREX_SPACEDIM == 3                                                              
            int trans2 = (direction + 2) % AMREX_SPACEDIM;  // Second transverse velocity (3D only)
#endif
            // **Direction to velocity mapping**
            int velocity_component[3] = {
                solver->variableIndex.UVEL,  // X-direction -> UVEL
                solver->variableIndex.VVEL,  // Y-direction -> VVEL
                solver->variableIndex.WVEL   // Z-direction -> WVEL
            };

        amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {


           // Construct index arrays with correct number of components
            int index[3];
           // int left_index[3];
           // int right_index[3];
            int lower_bounds[3];
            int upper_bounds[3];
        
#if AMREX_SPACEDIM == 2
            index[0] = iface; index[1] = jface; index[2] = kface;
            //left_index[0] = iface; left_index[1] = jface; left_index[2] = kface;
            //right_index[0] = iface; right_index[1] = jface; right_index[2] = kface;

            lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); lower_bounds[2] = 0;
            upper_bounds[0] = face_bx_with_ghosts.bigEnd(0)-1; upper_bounds[1] = face_bx_with_ghosts.bigEnd(1)-1; upper_bounds[2] = 0;
#else
            index[0] = iface; index[1] = jface; index[2] = kface;
            //left_index[0] = iface; left_index[1] = jface; left_index[2] = kface;
            //right_index[0] = iface; right_index[1] = jface; right_index[2] = kface;

            lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); lower_bounds[2] = face_bx_with_ghosts.smallEnd(2);
            upper_bounds[0] = face_bx_with_ghosts.bigEnd(0)-1; upper_bounds[1] = face_bx_with_ghosts.bigEnd(1)-1; upper_bounds[2] = face_bx_with_ghosts.bigEnd(2);
#endif

            // Clamp the main index
            ClampIndices(index, lower_bounds, upper_bounds);
                
            // Shift and clamp for left and right indices
            //int left_offset[3] = {0, 0, 0};
            //int right_offset[3] = {0, 0, 0};
            //left_offset[direction] = 0;  // No Left offset. For Flux left has the same indices as the cell center indices.
            //right_offset[direction] = 1;  // Shift right

            //ShiftAndClampIndices(left_index, left_offset, lower_bounds, upper_bounds, direction);
            //ShiftAndClampIndices(right_index, right_offset, lower_bounds, upper_bounds, direction);

            // **Left and right primitive variables (QL and QR)**
            Set::Scalar rho_L = QL_arr(index[0], index[1], index[2], solver->variableIndex.DENS);
            Set::Scalar rho_R = QR_arr(index[0], index[1], index[2], solver->variableIndex.DENS);

            rho_L = std::max(rho_L, 1e-8);  // Avoid zero density
            rho_R = std::max(rho_R, 1e-8);  // Avoid zero density
        
            Set::Scalar e_L = std::max(QL_arr(index[0], index[1], index[2], solver->variableIndex.IE), 1e-8);
            Set::Scalar e_R = std::max(QR_arr(index[0], index[1], index[2], solver->variableIndex.IE), 1e-8);
        
            Set::Scalar pressure_L = (gamma - 1.0) * e_L;
            Set::Scalar pressure_R = (gamma - 1.0) * e_R;
        
            Set::Scalar vel_L[3] = {0.0, 0.0, 0.0};
            Set::Scalar vel_R[3] = {0.0, 0.0, 0.0};

            vel_L[normal] = QL_arr(index[0], index[1], index[2], velocity_component[normal]);
            vel_R[normal] = QR_arr(index[0], index[1], index[2], velocity_component[normal]);
        
#if AMREX_SPACEDIM >= 2    
                vel_L[trans1] = QL_arr(index[0], index[1], index[2], velocity_component[trans1]);
                vel_R[trans1] = QR_arr(index[0], index[1], index[2], velocity_component[trans1]);
#endif
#if AMREX_SPACEDIM == 3    
                vel_L[trans2] = QL_arr(index[0], index[1], index[2], velocity_component[trans2]);
                vel_R[trans2] = QR_arr(index[0], index[1], index[2], velocity_component[trans2]);
#endif            
        

            // **Compute Kinetic Energy for Left and Right States**
            Set::Scalar KE_L = 0.5 * (vel_L[0] * vel_L[0] + vel_L[1] * vel_L[1] + vel_L[2] * vel_L[2]);  // Kinetic energy for left state
            Set::Scalar KE_R = 0.5 * (vel_R[0] * vel_R[0] + vel_R[1] * vel_R[1] + vel_R[2] * vel_R[2]);  // Kinetic energy for right state

            // **Total Energy (E = e + KE)**
            Set::Scalar E_L = e_L + KE_L;  // Total energy for left state
            Set::Scalar E_R = e_R + KE_R;  // Total energy for right state

            // **Speed of sound**
            Set::Scalar c_L = std::sqrt(gamma * std::max(pressure_L, 1e-8) / rho_L);
            Set::Scalar c_R = std::sqrt(gamma * std::max(pressure_R, 1e-8) / rho_R);
        
            // **Maximum wave speed (lambda_max)**
            Set::Scalar lambda_max = std::max(std::abs(vel_L[normal]) + c_L, std::abs(vel_R[normal]) + c_R);
        

                //debug.CheckNegativePressureAndAbort(iface, jface, kface, pressure_L, "Pressure Field Left", "LLF Step", true);
                //debug.CheckNegativePressureAndAbort(iface, jface, kface, pressure_R, "Pressure Field Right", "LLF Step", true);

            // **Flux arrays**
            int num_components = solver->number_of_components;
            std::vector<Set::Scalar> F_L(num_components, 0.0);
            std::vector<Set::Scalar> F_R(num_components, 0.0);
        
            F_L[solver->variableIndex.DENS] = rho_L * vel_L[direction];  // Mass flux
            F_R[solver->variableIndex.DENS] = rho_R * vel_R[direction];
        
            F_L[solver->variableIndex.UVEL + normal] = rho_L * vel_L[normal] * vel_L[normal] + pressure_L;
            F_R[solver->variableIndex.UVEL + normal] = rho_R * vel_R[normal] * vel_R[normal] + pressure_R;

#if AMREX_SPACEDIM >= 2
                F_L[solver->variableIndex.UVEL + trans1] = rho_L * vel_L[normal] * vel_L[trans1];
                F_R[solver->variableIndex.UVEL + trans1] = rho_R * vel_R[normal] * vel_R[trans1];
#endif
#if AMREX_SPACEDIM == 3
                F_L[solver->variableIndex.UVEL + trans2] = rho_L * vel_L[normal] * vel_L[trans2];
                F_R[solver->variableIndex.UVEL + trans2] = rho_R * vel_R[normal] * vel_R[trans2];
#endif
                    
            F_L[solver->variableIndex.IE] = (rho_L * E_L + pressure_L) * vel_L[normal];
            F_R[solver->variableIndex.IE] = (rho_R * E_R + pressure_R) * vel_R[normal];
        

            // **Conservative state vectors (left and right)**
            std::vector<Set::Scalar> QL_cons(num_components, 0.0);
            std::vector<Set::Scalar> QR_cons(num_components, 0.0);
            
            // **Convert primitive variables to conservative variables**
            QL_cons[solver->variableIndex.DENS] = rho_L;  // Density
            QR_cons[solver->variableIndex.DENS] = rho_R;
            
            QL_cons[solver->variableIndex.UVEL + normal] = rho_L * vel_L[normal];  // Momentum in normal direction
            QR_cons[solver->variableIndex.UVEL + normal] = rho_R * vel_R[normal];
            
#if AMREX_SPACEDIM >= 2
            QL_cons[solver->variableIndex.UVEL + trans1] = rho_L * vel_L[trans1];  // Momentum in first transverse direction
            QR_cons[solver->variableIndex.UVEL + trans1] = rho_R * vel_R[trans1];
#endif
            
#if AMREX_SPACEDIM == 3
            QL_cons[solver->variableIndex.UVEL + trans2] = rho_L * vel_L[trans2];  // Momentum in second transverse direction
            QR_cons[solver->variableIndex.UVEL + trans2] = rho_R * vel_R[trans2];
#endif
            
            QL_cons[solver->variableIndex.IE] = rho_L * E_L;  // Total energy
            QR_cons[solver->variableIndex.IE] = rho_R * E_R;
            

            // **Local Lax-Friedrichs flux calculation**
            for (int n = 0; n < num_components; ++n) {

                flux_arr(index[0], index[1], index[2], n) =  0.5 * (F_L[n] + F_R[n]) - 0.5 * lambda_max * (QR_cons[n] - QL_cons[n]);

            }
        
            // Debug check for NaNs
               // Debug
                /*debug.DebugLocalLaxFriedrichsFlux(iface, jface, kface,num_components,
                                                  rho_L, vel_L[0], vel_L[1], vel_L[2], e_L, pressure_L,
                                                  rho_R, vel_R[0], vel_R[1], vel_R[2], e_R, pressure_R,
                                                  c_L, c_R, lambda_max,
                                                  F_L.data(), F_R.data(),
                                                  "ComputeFlux", true);

            debug.CheckNaNAndAbort(iface, jface, kface, lambda_max, "Lambda Max", "wave speed", true);*/
        });
        
        
        }
        
        //Util::Message(INFO, "Local Lax-Friedrichs flux computed for direction: " + std::to_string(direction));
    }
};

}  // namespace Numeric

#endif  // NUMERIC_FLUXHANDLER_H


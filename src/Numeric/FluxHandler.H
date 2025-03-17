#ifndef NUMERIC_FLUXHANDLER_H
#define NUMERIC_FLUXHANDLER_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Set/Set.H"
#include "Integrator/ScimitarX.H"
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"
#include "Model/Fluid/Fluid.H" 
#include "Numeric/Stencil.H" 
#include "Numeric/NumericTypes.H"
#include "Numeric/SolverCapabilities.H"
#include "Numeric/IntegratorVariableAccessLayer.H"
#include <memory>
#include <cmath>
#include <algorithm>

//forward declaration for ScimitarX
namespace Integrator {
    class ScimitarX;

}

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
    virtual void ReconstructFluxVariables(int direction, int lev, IntegratorDerivedClass* solver, amrex::MultiFab& working_buffer, amrex::MultiFab& QL_stencil, amrex::MultiFab& QR_stencil) const = 0;
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

    // Constructor accepts a variable accessor
    FluxHandler(std::shared_ptr<GenericVariableAccessor> accessor) 
        : variable_accessor(accessor) {
        if (!accessor) {
            Util::Warning(INFO, "FluxHandler created with null variable accessor");
        }
    }

    // Set reconstruction method 
    void SetReconstruction(std::shared_ptr<FluxReconstruction<IntegratorDerivedClass>> recon) {
        reconstruction = std::move(recon);
    }

    // Set flux computation method
    void SetFluxMethod(std::shared_ptr<FluxMethod<IntegratorDerivedClass>> method) {
        flux_method = std::move(method);
    }

    // ConstructFluxes method with improved validation and error handling
    void ConstructFluxes(int lev, IntegratorDerivedClass* solver) const {
        // Validate reconstruction and flux method are set
        if (!reconstruction || !flux_method) {
            Util::Abort(INFO, "FluxHandler: Reconstruction or flux method is not set.");
        }

        // Check if variable accessor is valid
        if (!variable_accessor) {
            Util::Abort(INFO, "FluxHandler: Variable accessor is not set.");
        }

        // Get solver capabilities for validation
        auto solverCapabilities = variable_accessor->getSolverCapabilities();
        if (!solverCapabilities) {
            Util::Abort(INFO, "FluxHandler: Unable to retrieve solver capabilities.");
        }

        // Get current reconstruction mode
        ReconstructionMode mode = variable_accessor->getCurrentReconstructionMode();

        // Validate method combinations through solver capabilities
        auto validationResult = solverCapabilities->validateMethodCombination(
            FluxReconstructionType::WENO,  // Default assumption, could be parameterized
            FluxScheme::LocalLaxFriedrichs, // Default assumption, could be parameterized
            TimeSteppingSchemeType::RK3,    // Default assumption, could be parameterized
            mode
        );

        // Process for each spatial direction
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            // Create face-centered BoxArray
            amrex::BoxArray face_ba = solver->grids[lev];
            amrex::IntVect face_type = amrex::IntVect::TheZeroVector();
            face_type[dir] = 1;
            face_ba.convert(face_type);

            // Determine required ghost cells based on reconstruction method and mode
            int ghost_cells = solver->number_of_ghost_cells;

            // Create stencil MultiFabs for left and right states
            amrex::MultiFab QL_stencil(face_ba, solver->dmap[lev], solver->number_of_components, ghost_cells);
            amrex::MultiFab QR_stencil(face_ba, solver->dmap[lev], solver->number_of_components, ghost_cells);

            // Create working buffer using enhanced variable accessor method
            amrex::MultiFab working_buffer = variable_accessor->CreateWorkingBuffer(
                solver->grids[lev], 
                solver->dmap[lev], 
                solver->number_of_components, 
                ghost_cells
            );

            // Transform variables using enhanced method with validation
            variable_accessor->TransformVariables(
                lev, solver, working_buffer, mode, validationResult, dir
            );

            // Perform flux reconstruction
            reconstruction->ReconstructFluxVariables(dir, lev, solver, working_buffer, QL_stencil, QR_stencil);

            // Reverse transform with validation
            variable_accessor->ReverseTransform(
                lev, solver, QL_stencil, QR_stencil, mode, validationResult, dir
            );

            // Compute numerical flux
            flux_method->ComputeFlux(dir, lev, solver, QL_stencil, QR_stencil);
        }
    }

    // Get the current variable accessor
    std::shared_ptr<GenericVariableAccessor> GetVariableAccessor() const {
        return variable_accessor;
    }

private:
    // Reconstruction and flux method pointers
    std::shared_ptr<FluxReconstruction<IntegratorDerivedClass>> reconstruction;
    std::shared_ptr<FluxMethod<IntegratorDerivedClass>> flux_method;
    std::shared_ptr<GenericVariableAccessor> variable_accessor;
};

// First-order upwind reconstruction
template <typename IntegratorDerivedClass>
class FirstOrderReconstruction : public FluxReconstruction<IntegratorDerivedClass> {
public:
// ReconstructFlux: First-order upwind flux reconstruction with branch-free framework
void ReconstructFluxVariables(int direction, int lev, IntegratorDerivedClass* solver,amrex::MultiFab& working_buffer, amrex::MultiFab& QL_stencil, amrex::MultiFab& QR_stencil) const override {
    
    const int nghosts = solver->number_of_ghost_cells; 
// Loop through face-centered indices for flux reconstruction
//for (amrex::MFIter mfi(QL_stencil, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
    const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);

    auto const& QL_arr = QL_stencil.array(mfi);  // Left reconstructed state
    auto const& QR_arr = QR_stencil.array(mfi);  // Right reconstructed state
    auto const& W_arr = working_buffer.array(mfi);  // Primitive variables

    const int num_components = solver->number_of_components;  // Number of components (e.g., rho, u, v, E)


    amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {

    // Construct index arrays with correct number of components
    int index[3];
    int left_index[3];
    int right_index[3];
    int lower_bounds[3];
    int upper_bounds[3];

    // Set up indices
    index[0] = iface; index[1] = jface; index[2] = kface;
    left_index[0] = iface; left_index[1] = jface; left_index[2] = kface;
    right_index[0] = iface; right_index[1] = jface; right_index[2] = kface;

    // Define bounds including ghost cells
    lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); 
    lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); 
    lower_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.smallEnd(2) : 0;
    upper_bounds[0] = face_bx_with_ghosts.bigEnd(0) -1;
    upper_bounds[1] = face_bx_with_ghosts.bigEnd(1) -1;
    upper_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.bigEnd(2) - 1 : 0;

    // Clamp the main index
    ClampIndices(index, lower_bounds, upper_bounds);

    // Shift and clamp for left and right indices
    int left_offset[3] = {0, 0, 0};
    int right_offset[3] = {0, 0, 0};
    left_offset[direction] = 0;  // Shift left
    right_offset[direction] = 1;  // Shift right

    ShiftAndClampIndices(left_index, left_offset, lower_bounds, upper_bounds, direction);
    ShiftAndClampIndices(right_index, right_offset, lower_bounds, upper_bounds, direction);
    

        for (int n = 0; n < num_components; ++n) {
            
            const Set::Scalar QL = W_arr(left_index[0], left_index[1], left_index[2], n);
            const Set::Scalar QR = W_arr(right_index[0], right_index[1], right_index[2], n);  // Since we are using first-order upwinding  
            // Assign values to the stencil arrays
            QL_arr(index[0], index[1], index[2], n) = QL;
            QR_arr(index[0], index[1], index[2], n) = QR;
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

        const Set::Scalar gamma = 1.4;  // Ratio of specific heats

        const int nghosts = solver->number_of_ghost_cells;

        //for (amrex::MFIter mfi(QL_stencil, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {

            const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts); //Includes two ghost layers

            auto const& flux_arr = (direction == FluxHandler<IntegratorDerivedClass>::Xdir) ? solver->XFlux_mf.Patch(lev, mfi) :
                                    (direction == FluxHandler<IntegratorDerivedClass>::Ydir) ? solver->YFlux_mf.Patch(lev, mfi) :
                                    solver->ZFlux_mf.Patch(lev, mfi);
            auto const& QL_arr = QL_stencil.array(mfi);
            auto const& QR_arr = QR_stencil.array(mfi);
            auto const& p_arr  = solver->PVec_mf.Patch(lev,mfi);

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
            int left_index[3];
            int right_index[3];
            int lower_bounds[3];
            int upper_bounds[3];
        
            // Set up indices
            index[0] = iface; index[1] = jface; index[2] = kface;
            left_index[0] = iface; left_index[1] = jface; left_index[2] = kface;
            right_index[0] = iface; right_index[1] = jface; right_index[2] = kface;
        
            // Define bounds including ghost cells
            lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); 
            lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); 
            lower_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.smallEnd(2) : 0;
            upper_bounds[0] = face_bx_with_ghosts.bigEnd(0) -1;
            upper_bounds[1] = face_bx_with_ghosts.bigEnd(1) -1;
            upper_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.bigEnd(2) - 1 : 0;
        
            // Clamp the main index
            ClampIndices(index, lower_bounds, upper_bounds);
        
            // Shift and clamp for left and right indices
            int left_offset[3] = {0, 0, 0};
            int right_offset[3] = {0, 0, 0};
            left_offset[direction] = 0;  // Shift left
            right_offset[direction] = 1;  // Shift right
        
            ShiftAndClampIndices(left_index, left_offset, lower_bounds, upper_bounds, direction);
            ShiftAndClampIndices(right_index, right_offset, lower_bounds, upper_bounds, direction);
                        
            // **Left and right primitive variables (QL and QR)**
            Set::Scalar rho_L = p_arr(left_index[0], left_index[1], left_index[2], solver->variableIndex.DENS);
            Set::Scalar rho_R = p_arr(right_index[0], right_index[1], right_index[2], solver->variableIndex.DENS);

            rho_L = std::max(rho_L, 1e-8);  // Avoid zero density
            rho_R = std::max(rho_R, 1e-8);  // Avoid zero density
        
            Set::Scalar e_L = std::max(p_arr(left_index[0], left_index[1], left_index[2], solver->variableIndex.IE), 1e-8);
            Set::Scalar e_R = std::max(p_arr(right_index[0], right_index[1], right_index[2], solver->variableIndex.IE), 1e-8);
        
            Set::Scalar pressure_L = (gamma - 1.0) * rho_L * e_L;
            Set::Scalar pressure_R = (gamma - 1.0) * rho_R * e_R;
        
            Set::Scalar vel_L[3] = {0.0, 0.0, 0.0};
            Set::Scalar vel_R[3] = {0.0, 0.0, 0.0};

            vel_L[normal] = p_arr(left_index[0], left_index[1], left_index[2], velocity_component[normal]);
            vel_R[normal] = p_arr(right_index[0], right_index[1], right_index[2], velocity_component[normal]);
        
#if AMREX_SPACEDIM >= 2    
                vel_L[trans1] = p_arr(left_index[0], left_index[1], left_index[2], velocity_component[trans1]);
                vel_R[trans1] = p_arr(right_index[0], right_index[1], right_index[2], velocity_component[trans1]);
#endif
#if AMREX_SPACEDIM == 3    
                vel_L[trans2] = p_arr(left_index[0], left_index[1], left_index[2], velocity_component[trans2]);
                vel_R[trans2] = p_arr(right_index[0], right_index[1], right_index[2], velocity_component[trans2]);
#endif            
        

            // **Compute Kinetic Energy for Left and Right States**
            Set::Scalar KE_L = 0.5 * (vel_L[0] * vel_L[0] + vel_L[1] * vel_L[1] + vel_L[2] * vel_L[2]);  // Kinetic energy for left state
            Set::Scalar KE_R = 0.5 * (vel_R[0] * vel_R[0] + vel_R[1] * vel_R[1] + vel_R[2] * vel_R[2]);  // Kinetic energy for right state

            // **Total Internal Energy (E = e + KE)**
            Set::Scalar E_L = e_L + KE_L;  // Total internal energy for left state
            Set::Scalar E_R = e_R + KE_R;  // Total internal energy for right state

            // **Speed of sound**
            Set::Scalar c_L = std::sqrt(gamma * std::max(pressure_L, 1e-8) / rho_L);
            Set::Scalar c_R = std::sqrt(gamma * std::max(pressure_R, 1e-8) / rho_R);
        
            // **Maximum wave speed (lambda_max)**
            Set::Scalar lambda_max = std::max(std::abs(vel_L[normal]) + c_L, std::abs(vel_R[normal]) + c_R);
        
            // **Flux arrays**
            int num_components = solver->number_of_components;
            std::vector<Set::Scalar> F_L(num_components, 0.0);
            std::vector<Set::Scalar> F_R(num_components, 0.0);

            F_L[solver->variableIndex.DENS] = QL_arr(index[0], index[1], index[2], solver->variableIndex.DENS);  // Mass flux
            F_R[solver->variableIndex.DENS] = QR_arr(index[0], index[1], index[2], solver->variableIndex.DENS);
        
            F_L[solver->variableIndex.UVEL + normal] = QL_arr(index[0], index[1], index[2], solver->variableIndex.UVEL + normal);
            F_R[solver->variableIndex.UVEL + normal] = QR_arr(index[0], index[1], index[2], solver->variableIndex.UVEL + normal);

#if AMREX_SPACEDIM >= 2
            F_L[solver->variableIndex.UVEL + trans1] = QL_arr(index[0], index[1], index[2], solver->variableIndex.UVEL + trans1);
            F_R[solver->variableIndex.UVEL + trans1] = QR_arr(index[0], index[1], index[2], solver->variableIndex.UVEL + trans1);
#endif
#if AMREX_SPACEDIM == 3
            F_L[solver->variableIndex.UVEL + trans2] = QL_arr(index[0], index[1], index[2], solver->variableIndex.UVEL + trans2);
            F_R[solver->variableIndex.UVEL + trans2] = QR_arr(index[0], index[1], index[2], solver->variableIndex.UVEL + trans2);
#endif
                    
            F_L[solver->variableIndex.IE] = QL_arr(index[0], index[1], index[2], solver->variableIndex.IE);
            F_R[solver->variableIndex.IE] = QR_arr(index[0], index[1], index[2], solver->variableIndex.IE);
        

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
        
        });
        
        
        }
        
    }
};

// HLLC flux method
template <typename IntegratorDerivedClass>
class HLLCMethod : public FluxMethod<IntegratorDerivedClass> {
public:
    void ComputeFlux(int direction, int lev, IntegratorDerivedClass* solver, 
                    const amrex::MultiFab& QL_stencil, 
                    const amrex::MultiFab& QR_stencil) const override {

        const Set::Scalar gamma = 1.4;  // Ratio of specific heats

        const int nghosts = solver->number_of_ghost_cells;

        // Loop through the grid
        for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
            const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);

            // Get appropriate flux array based on direction
            auto const& flux_arr = (direction == FluxHandler<IntegratorDerivedClass>::Xdir) ? solver->XFlux_mf.Patch(lev, mfi) :
                                  (direction == FluxHandler<IntegratorDerivedClass>::Ydir) ? solver->YFlux_mf.Patch(lev, mfi) :
                                   solver->ZFlux_mf.Patch(lev, mfi);
            
            auto const& QL_arr = QL_stencil.array(mfi);
            auto const& QR_arr = QR_stencil.array(mfi);
            auto const& p_arr = solver->PVec_mf.Patch(lev, mfi);

            // Direction-dependent indices
            int normal = direction;  // Normal velocity component
            int trans1 = (direction + 1) % AMREX_SPACEDIM;  // First transverse velocity
#if AMREX_SPACEDIM == 3
            int trans2 = (direction + 2) % AMREX_SPACEDIM;  // Second transverse velocity (3D only)
#endif
            // Map direction to velocity components
            int velocity_component[3] = {
                solver->variableIndex.UVEL,  // X-direction -> UVEL
                solver->variableIndex.VVEL,  // Y-direction -> VVEL
                solver->variableIndex.WVEL   // Z-direction -> WVEL
            };

            int num_components = 5;

            // Parallel loop over cells
            amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {
                // Setup indices
                int index[3];
                int left_index[3];
                int right_index[3];
                int lower_bounds[3];
                int upper_bounds[3];
                
                // Initialize index arrays
                index[0] = iface; index[1] = jface; index[2] = kface;
                left_index[0] = iface; left_index[1] = jface; left_index[2] = kface;
                right_index[0] = iface; right_index[1] = jface; right_index[2] = kface;
                
                // Define bounds including ghost cells
                lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); 
                lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); 
                lower_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.smallEnd(2) : 0;
                upper_bounds[0] = face_bx_with_ghosts.bigEnd(0) - 1;
                upper_bounds[1] = face_bx_with_ghosts.bigEnd(1) - 1;
                upper_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.bigEnd(2) - 1 : 0;
                
                // Clamp indices to ensure they're within bounds
                ClampIndices(index, lower_bounds, upper_bounds);
                
                // Shift and clamp for left and right indices
                int left_offset[3] = {0, 0, 0};
                int right_offset[3] = {0, 0, 0};
                right_offset[direction] = 1;  // Shift right
                
                ShiftAndClampIndices(left_index, left_offset, lower_bounds, upper_bounds, direction);
                ShiftAndClampIndices(right_index, right_offset, lower_bounds, upper_bounds, direction);
                
                // Get primitive variables for left and right states
                // Left state
                Set::Scalar rho_L = p_arr(left_index[0], left_index[1], left_index[2], solver->variableIndex.DENS);
                rho_L = std::max(rho_L, 1e-8);  // Avoid zero density
                
                Set::Scalar e_L = std::max(p_arr(left_index[0], left_index[1], left_index[2], solver->variableIndex.IE), 1e-8);
                Set::Scalar pressure_L = (gamma - 1.0) * rho_L * e_L;
                
                Set::Scalar vel_L[3] = {0.0, 0.0, 0.0};
                vel_L[normal] = p_arr(left_index[0], left_index[1], left_index[2], velocity_component[normal]);
#if AMREX_SPACEDIM >= 2    
                vel_L[trans1] = p_arr(left_index[0], left_index[1], left_index[2], velocity_component[trans1]);
#endif
#if AMREX_SPACEDIM == 3    
                vel_L[trans2] = p_arr(left_index[0], left_index[1], left_index[2], velocity_component[trans2]);
#endif
                
                // Right state
                Set::Scalar rho_R = p_arr(right_index[0], right_index[1], right_index[2], solver->variableIndex.DENS);
                rho_R = std::max(rho_R, 1e-8);  // Avoid zero density
                
                Set::Scalar e_R = std::max(p_arr(right_index[0], right_index[1], right_index[2], solver->variableIndex.IE), 1e-8);
                Set::Scalar pressure_R = (gamma - 1.0) * rho_R * e_R;
                
                Set::Scalar vel_R[3] = {0.0, 0.0, 0.0};
                vel_R[normal] = p_arr(right_index[0], right_index[1], right_index[2], velocity_component[normal]);
#if AMREX_SPACEDIM >= 2    
                vel_R[trans1] = p_arr(right_index[0], right_index[1], right_index[2], velocity_component[trans1]);
#endif
#if AMREX_SPACEDIM == 3    
                vel_R[trans2] = p_arr(right_index[0], right_index[1], right_index[2], velocity_component[trans2]);
#endif
                
                // Compute kinetic energy
                Set::Scalar KE_L = 0.5 * (vel_L[0] * vel_L[0] + vel_L[1] * vel_L[1] + vel_L[2] * vel_L[2]);
                Set::Scalar KE_R = 0.5 * (vel_R[0] * vel_R[0] + vel_R[1] * vel_R[1] + vel_R[2] * vel_R[2]);
                
                // Total energy
                Set::Scalar E_L = e_L + KE_L;
                Set::Scalar E_R = e_R + KE_R;
                
                // Sound speeds
                Set::Scalar c_L = std::sqrt(gamma * pressure_L / rho_L);
                Set::Scalar c_R = std::sqrt(gamma * pressure_R / rho_R);
                
                // Create vectors for Roe averaging
                Set::MultiVector WL(num_components);  // num_components is a local variable and is set to 5
                Set::MultiVector WR(num_components);
                WL.setZero();
                WR.setZero();
                
                // Fill left state vector
                WL(0) = rho_L;
                WL(1) = vel_L[0];
                WL(2) = vel_L[1];
                WL(3) = vel_L[2] * (AMREX_SPACEDIM == 3 ? 1.0 : 0.0);
                WL(4) = e_L;
                
                // Fill right state vector
                WR(0) = rho_R;
                WR(1) = vel_R[0];
                WR(2) = vel_R[1];
                WR(3) = vel_R[2] * (AMREX_SPACEDIM == 3 ? 1.0 : 0.0);
                WR(4) = e_R;
                
                // Compute Roe averages using the function from IntegratorVariableAccessLayer
                Set::MultiVector Wavg = Numeric::CompressibleEuler::CompressibleEulerVariableAccessor::ComputeRoeAverages(WL, WR, num_components);
                
                // Extract Roe-averaged values
                Set::Scalar u_Roe = Wavg(1);
                Set::Scalar c_Roe = std::sqrt(gamma * (gamma - 1.0) * (Wavg(4) - 0.5 * 
                                   (Wavg(1)*Wavg(1) + Wavg(2)*Wavg(2) + Wavg(3)*Wavg(3))));
                
                // Compute wave speeds using Batten approach
                Set::Scalar SL = std::min(vel_L[normal] - c_L, u_Roe - c_Roe);
                Set::Scalar SR = std::max(vel_R[normal] + c_R, u_Roe + c_Roe);
                
                // Calculate middle wave speed (S*)
                Set::Scalar rho_L_inv = 1.0 / rho_L;
                Set::Scalar rho_R_inv = 1.0 / rho_R;
                Set::Scalar S_star = (pressure_R - pressure_L + 
                                    rho_L * vel_L[normal] * (SL - vel_L[normal]) - 
                                    rho_R * vel_R[normal] * (SR - vel_R[normal])) / 
                                   (rho_L * (SL - vel_L[normal]) - rho_R * (SR - vel_R[normal]));
                
                // Conservative variables for left and right states
                std::vector<Set::Scalar> QL_cons(solver->number_of_components, 0.0);
                std::vector<Set::Scalar> QR_cons(solver->number_of_components, 0.0);
                
                // Convert primitive to conservative variables
                QL_cons[solver->variableIndex.DENS] = rho_L;  // Density
                QR_cons[solver->variableIndex.DENS] = rho_R;
                
                QL_cons[solver->variableIndex.UVEL + normal] = rho_L * vel_L[normal];  // Normal momentum
                QR_cons[solver->variableIndex.UVEL + normal] = rho_R * vel_R[normal];
                
#if AMREX_SPACEDIM >= 2
                QL_cons[solver->variableIndex.UVEL + trans1] = rho_L * vel_L[trans1];  // Transverse momentum 1
                QR_cons[solver->variableIndex.UVEL + trans1] = rho_R * vel_R[trans1];
#endif
                
#if AMREX_SPACEDIM == 3
                QL_cons[solver->variableIndex.UVEL + trans2] = rho_L * vel_L[trans2];  // Transverse momentum 2
                QR_cons[solver->variableIndex.UVEL + trans2] = rho_R * vel_R[trans2];
#endif
                
                QL_cons[solver->variableIndex.IE] = rho_L * E_L;  // Total energy
                QR_cons[solver->variableIndex.IE] = rho_R * E_R;
                
                // Compute physical fluxes for left and right states
                std::vector<Set::Scalar> FL(solver->number_of_components, 0.0);
                std::vector<Set::Scalar> FR(solver->number_of_components, 0.0);
                
                // Mass flux
                FL[solver->variableIndex.DENS] = rho_L * vel_L[normal];
                FR[solver->variableIndex.DENS] = rho_R * vel_R[normal];
                
                // Momentum fluxes
                FL[solver->variableIndex.UVEL + normal] = rho_L * vel_L[normal] * vel_L[normal] + pressure_L;
                FR[solver->variableIndex.UVEL + normal] = rho_R * vel_R[normal] * vel_R[normal] + pressure_R;
                
#if AMREX_SPACEDIM >= 2
                FL[solver->variableIndex.UVEL + trans1] = rho_L * vel_L[normal] * vel_L[trans1];
                FR[solver->variableIndex.UVEL + trans1] = rho_R * vel_R[normal] * vel_R[trans1];
#endif
                
#if AMREX_SPACEDIM == 3
                FL[solver->variableIndex.UVEL + trans2] = rho_L * vel_L[normal] * vel_L[trans2];
                FR[solver->variableIndex.UVEL + trans2] = rho_R * vel_R[normal] * vel_R[trans2];
#endif
                
                // Energy flux
                FL[solver->variableIndex.IE] = (rho_L * E_L + pressure_L) * vel_L[normal];
                FR[solver->variableIndex.IE] = (rho_R * E_R + pressure_R) * vel_R[normal];
                
                // Compute star state variables
                std::vector<Set::Scalar> QL_star(solver->number_of_components, 0.0);
                std::vector<Set::Scalar> QR_star(solver->number_of_components, 0.0);
                
                // Density ratios for star states
                Set::Scalar rho_L_star_factor = (SL - vel_L[normal]) / (SL - S_star);
                Set::Scalar rho_R_star_factor = (SR - vel_R[normal]) / (SR - S_star);
                
                // Left star state
                QL_star[solver->variableIndex.DENS] = rho_L * rho_L_star_factor;
                QL_star[solver->variableIndex.UVEL + normal] = QL_star[solver->variableIndex.DENS] * S_star;
#if AMREX_SPACEDIM >= 2
                QL_star[solver->variableIndex.UVEL + trans1] = QL_star[solver->variableIndex.DENS] * vel_L[trans1];
#endif
#if AMREX_SPACEDIM == 3
                QL_star[solver->variableIndex.UVEL + trans2] = QL_star[solver->variableIndex.DENS] * vel_L[trans2];
#endif
                QL_star[solver->variableIndex.IE] = QL_star[solver->variableIndex.DENS] * 
                                                   (E_L + (S_star - vel_L[normal]) * 
                                                   (S_star + pressure_L / (rho_L * (SL - vel_L[normal]))));
                
                // Right star state
                QR_star[solver->variableIndex.DENS] = rho_R * rho_R_star_factor;
                QR_star[solver->variableIndex.UVEL + normal] = QR_star[solver->variableIndex.DENS] * S_star;
#if AMREX_SPACEDIM >= 2
                QR_star[solver->variableIndex.UVEL + trans1] = QR_star[solver->variableIndex.DENS] * vel_R[trans1];
#endif
#if AMREX_SPACEDIM == 3
                QR_star[solver->variableIndex.UVEL + trans2] = QR_star[solver->variableIndex.DENS] * vel_R[trans2];
#endif
                QR_star[solver->variableIndex.IE] = QR_star[solver->variableIndex.DENS] * 
                                                   (E_R + (S_star - vel_R[normal]) * 
                                                   (S_star + pressure_R / (rho_R * (SR - vel_R[normal]))));
                
                // Compute fluxes for star states
                std::vector<Set::Scalar> FL_star(solver->number_of_components, 0.0);
                std::vector<Set::Scalar> FR_star(solver->number_of_components, 0.0);
                
                for (int n = 0; n < solver->number_of_components; ++n) {
                    FL_star[n] = FL[n] + SL * (QL_star[n] - QL_cons[n]);
                    FR_star[n] = FR[n] + SR * (QR_star[n] - QR_cons[n]);
                }
                
                // Apply HLLC flux based on wave speeds
                std::vector<Set::Scalar> HLLC_flux(solver->number_of_components, 0.0);
                
                for (int n = 0; n < solver->number_of_components; ++n) {
                    if (SL >= 0.0) {
                        HLLC_flux[n] = FL[n];
                    } else if (S_star >= 0.0) {
                        HLLC_flux[n] = FL_star[n];
                    } else if (SR >= 0.0) {
                        HLLC_flux[n] = FR_star[n];
                    } else {
                        HLLC_flux[n] = FR[n];
                    }
                    
                    // Store computed flux
                    flux_arr(index[0], index[1], index[2], n) = HLLC_flux[n];
                }
            });
        }
    }
};

}  // namespace Numeric

#endif  // NUMERIC_FLUXHANDLER_H

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
    // Enum to identify what type of entities we're reconstructing
    enum class EntityType {
        Variables,  // Reconstruct variables (for HLLC-like methods)
        Fluxes      // Reconstruct fluxes (for LLF-like methods)
    };

    virtual ~FluxReconstruction() = default;

    // Original method for flux-based reconstruction (for LLF-like methods)
    virtual void ReconstructFluxes(int direction, int lev, IntegratorDerivedClass* solver,
                                  const amrex::MultiFab& Flux_Negative, 
                                  const amrex::MultiFab& Flux_Positive,
                                  amrex::MultiFab& SummedFlux) const = 0;
    
    // New method for variable-based reconstruction (for HLLC-like methods)
    virtual void ReconstructVariables(int direction, int lev, IntegratorDerivedClass* solver,
                                     const amrex::MultiFab& VariableBuffer,
                                     amrex::MultiFab& LeftStates, 
                                     amrex::MultiFab& RightStates) const = 0;
};

// Base class for flux methods (e.g., Local Lax-Friedrichs, HLLC)
template <typename IntegratorDerivedClass>
class FluxMethod {
public:
    // Enum to determine reconstruction approach based on method requirements
    enum class ReconstructionTarget {
        Variables,  // Reconstruct variables first (for HLLC-like methods)
        Fluxes      // Reconstruct fluxes first (for LLF-like methods)
    };

    virtual ~FluxMethod() = default;

    // Return whether this method requires reconstructing variables or fluxes first
    virtual ReconstructionTarget getReconstructionTarget() const = 0;

    // Method-specific flux computation implementations
    virtual void SplitFluxes(int direction, int lev, IntegratorDerivedClass* solver, 
                             amrex::MultiFab& VariableBuffer, amrex::MultiFab& CellFluxBuffer, 
                             amrex::MultiFab& Flux_Negative, amrex::MultiFab& Flux_Positive) const = 0;
    
    // New method for directly computing fluxes from reconstructed states (for HLLC-like methods)
    virtual void ComputeFluxes(int direction, int lev, IntegratorDerivedClass* solver,
                              const amrex::MultiFab& LeftStates, const amrex::MultiFab& RightStates, 
                              amrex::MultiFab& SummedFlux) const = 0;
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

    // Main method for constructing fluxes
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

        // Create VariableBuffer for reuse across all directions
        int ghost_cells = solver->number_of_ghost_cells;  
        amrex::MultiFab VariableBuffer = variable_accessor->CreateWorkingBuffer(
            solver->grids[lev], solver->dmap[lev], solver->number_of_components, ghost_cells
        );

        // Process each direction
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            // Choose processing path based on the flux method's reconstruction target
            if (flux_method->getReconstructionTarget() == FluxMethod<IntegratorDerivedClass>::ReconstructionTarget::Variables) {
                // For methods like HLLC that reconstruct variables first
                ProcessVariableBasedReconstruction(dir, lev, solver, VariableBuffer);
            } else {
                // For methods like LLF that reconstruct fluxes first
                ProcessFluxBasedReconstruction(dir, lev, solver, VariableBuffer);
            }
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

    // For variable-based reconstruction methods like HLLC
    void ProcessVariableBasedReconstruction(
        int direction, int lev, IntegratorDerivedClass* solver, 
        amrex::MultiFab& VariableBuffer) const;

    // For flux-based reconstruction methods like LLF
    void ProcessFluxBasedReconstruction(
        int direction, int lev, IntegratorDerivedClass* solver, 
        amrex::MultiFab& VariableBuffer) const;
};

// For variable-based reconstruction methods like HLLC
template <typename IntegratorDerivedClass>
void FluxHandler<IntegratorDerivedClass>::ProcessVariableBasedReconstruction(
    int direction, int lev, IntegratorDerivedClass* solver, amrex::MultiFab& VariableBuffer) const {
    
    // Validate inputs and components
    if (!reconstruction || !flux_method || !variable_accessor) {
        Util::Abort(INFO, "FluxHandler: Required components not initialized");
    }
    
    // Get solver capabilities for validation
    auto solverCapabilities = variable_accessor->getSolverCapabilities();
    if (!solverCapabilities) {
        Util::Abort(INFO, "FluxHandler: Unable to retrieve solver capabilities");
    }
    
    // Get current reconstruction mode
    ReconstructionMode mode = variable_accessor->getCurrentReconstructionMode();
    
    // Validate method combinations
    auto validationResult = solverCapabilities->validateMethodCombination(
        FluxReconstructionType::WENO,  // Assume WENO, could be parameterized
        FluxScheme::HLLC,             // For this path, we're using HLLC
        TimeSteppingSchemeType::RK3,  // Assume RK3, could be parameterized
        mode
    );
    
    // 1. Transform variables if needed (based on reconstruction mode)
    variable_accessor->TransformVariables(direction, lev, solver, VariableBuffer, 
                                         mode, validationResult);
    
    // 2. Create intermediate buffers for reconstructed left/right states
    int ghost_cells = solver->number_of_ghost_cells;


    // 3. Create a face-centered BoxArray for the flux
    amrex::BoxArray face_ba = solver->grids[lev];
    amrex::IntVect face_type = amrex::IntVect::TheZeroVector();
    face_type[direction] = 1;  // Set face-centered in the specified direction
    face_ba.convert(face_type);


    amrex::MultiFab LeftStates = variable_accessor->CreateWorkingBuffer(
        face_ba, solver->dmap[lev], solver->number_of_components, ghost_cells);
    
    amrex::MultiFab RightStates = variable_accessor->CreateWorkingBuffer(
        face_ba, solver->dmap[lev], solver->number_of_components, ghost_cells);
    
    // 4. Perform reconstruction on variables to get left/right states
    reconstruction->ReconstructVariables(direction, lev, solver, 
                                       VariableBuffer, LeftStates, RightStates);
    
    // 5. Transform reconstructed states back if needed (for characteristic mode)
    variable_accessor->ReverseTransformVariables(direction, lev, solver, 
                                                  LeftStates, RightStates, 
                                                  mode, validationResult);
    
    
    // 6. Create MultiFab for the computed flux
    amrex::MultiFab SummedFlux(face_ba, solver->dmap[lev], 
                              solver->number_of_components, ghost_cells);
    
    // 7. Compute fluxes from reconstructed states
    flux_method->ComputeFluxes(direction, lev, solver, 
                             LeftStates, RightStates, SummedFlux);
    
    // 8. Store fluxes in the appropriate directional flux arrays
    variable_accessor->StoreDirectionalFlux(direction, lev, solver, 
                                          SummedFlux, mode, validationResult);
}

// For flux-based reconstruction methods like LLF
template <typename IntegratorDerivedClass>
void FluxHandler<IntegratorDerivedClass>::ProcessFluxBasedReconstruction(
    int direction, int lev, IntegratorDerivedClass* solver, amrex::MultiFab& VariableBuffer) const {
    
    // Validate inputs and components
    if (!reconstruction || !flux_method || !variable_accessor) {
        Util::Abort(INFO, "FluxHandler: Required components not initialized");
    }
    
    // Get solver capabilities for validation
    auto solverCapabilities = variable_accessor->getSolverCapabilities();
    if (!solverCapabilities) {
        Util::Abort(INFO, "FluxHandler: Unable to retrieve solver capabilities");
    }
    
    // Get current reconstruction mode
    ReconstructionMode mode = variable_accessor->getCurrentReconstructionMode();
    
    // Validate method combinations
    auto validationResult = solverCapabilities->validateMethodCombination(
        FluxReconstructionType::WENO,  // Assume WENO, could be parameterized
        FluxScheme::LocalLaxFriedrichs,  // For this path, we're using LLF
        TimeSteppingSchemeType::RK3,   // Assume RK3, could be parameterized
        mode
    );
    
    // 1. Transform variables if needed (based on reconstruction mode)
    variable_accessor->TransformVariables(direction, lev, solver, VariableBuffer, 
                                         mode, validationResult);
    
    // 2. Create a cell flux buffer for storing computed fluxes
    amrex::MultiFab CellFluxBuffer = variable_accessor->CreateWorkingBuffer(
        solver->grids[lev], solver->dmap[lev], solver->number_of_components, 
        solver->number_of_ghost_cells);
    
    // 3. Compute and transform fluxes if needed
    variable_accessor->TransformFluxes(direction, lev, solver, CellFluxBuffer, 
                                      mode, validationResult);
    
    // 4. Create buffers for split fluxes (positive and negative)
    int ghost_cells = solver->number_of_ghost_cells;
    amrex::MultiFab Flux_Negative = variable_accessor->CreateWorkingBuffer(
        solver->grids[lev], solver->dmap[lev], solver->number_of_components, ghost_cells);
    
    amrex::MultiFab Flux_Positive = variable_accessor->CreateWorkingBuffer(
        solver->grids[lev], solver->dmap[lev], solver->number_of_components, ghost_cells);
    
    // 5. Split fluxes into positive and negative components
    flux_method->SplitFluxes(direction, lev, solver, VariableBuffer, 
                           CellFluxBuffer, Flux_Negative, Flux_Positive);
    
    // 6. Create a face-centered BoxArray for the flux
    amrex::BoxArray face_ba = solver->grids[lev];
    amrex::IntVect face_type = amrex::IntVect::TheZeroVector();
    face_type[direction] = 1;  // Set face-centered in the specified direction
    face_ba.convert(face_type);
    
    // 7. Create MultiFab for the reconstructed summed flux
    amrex::MultiFab SummedFlux(face_ba, solver->dmap[lev], 
                              solver->number_of_components, ghost_cells);
    
    // 8. Perform reconstruction on fluxes
    reconstruction->ReconstructFluxes(direction, lev, solver, 
                                    Flux_Negative, Flux_Positive, SummedFlux);
    
    // 9. Transform fluxes back if needed
    variable_accessor->ReverseTransformFluxes(direction, lev, solver, SummedFlux, 
                                           mode, validationResult);
    
    // 10. Store fluxes in the appropriate directional flux arrays
    variable_accessor->StoreDirectionalFlux(direction, lev, solver, 
                                          SummedFlux, mode, validationResult);
}

// First-order upwind reconstruction
template <typename IntegratorDerivedClass>
class FirstOrderReconstruction : public FluxReconstruction<IntegratorDerivedClass> {
public:
    // Implement flux-based reconstruction (for LLF-like methods)
    void ReconstructFluxes(int direction, int lev, IntegratorDerivedClass* solver,
                         const amrex::MultiFab& Flux_Negative, const amrex::MultiFab& Flux_Positive,
                         amrex::MultiFab& SummedFlux) const override {
        
        const int nghosts = solver->number_of_ghost_cells; 

        // Loop through face-centered indices for flux reconstruction
        for (amrex::MFIter mfi(SummedFlux, false); mfi.isValid(); ++mfi) {
            const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);

            auto const& FluxNeg_arr = Flux_Negative.array(mfi);  // Left reconstructed state
            auto const& FluxPos_arr = Flux_Positive.array(mfi);  // Right reconstructed state
            auto const& TotalFlux_arr = SummedFlux.array(mfi);  // Summed flux

            const int num_components = solver->number_of_components;  // Number of components

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
                upper_bounds[0] = face_bx_with_ghosts.bigEnd(0) - 1;
                upper_bounds[1] = face_bx_with_ghosts.bigEnd(1) - 1;
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
                    const Set::Scalar FL = FluxNeg_arr(left_index[0], left_index[1], left_index[2], n);
                    const Set::Scalar FR = FluxPos_arr(right_index[0], right_index[1], right_index[2], n);  
                    // Assign values to the stencil arrays
                    TotalFlux_arr(index[0], index[1], index[2], n) = FL + FR;
                }
            });
        }
    }
    
    // Implement variable-based reconstruction (for HLLC-like methods)
    void ReconstructVariables(int direction, int lev, IntegratorDerivedClass* solver,
                            const amrex::MultiFab& VariableBuffer,
                            amrex::MultiFab& LeftStates, amrex::MultiFab& RightStates) const override {
        
        const int nghosts = solver->number_of_ghost_cells;
        
        // Loop through grid cells to reconstruct left/right states at interfaces
        for (amrex::MFIter mfi(VariableBuffer, false); mfi.isValid(); ++mfi) {
            const amrex::Box& bx_with_ghosts = mfi.growntilebox();
            
            auto const& Var_arr = VariableBuffer.array(mfi);
            auto const& LeftStates_arr = LeftStates.array(mfi);
            auto const& RightStates_arr = RightStates.array(mfi);
            
            const int num_components = solver->number_of_components;
            
            amrex::ParallelFor(bx_with_ghosts, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Construct index arrays
                int index[3];
                int left_cell[3];
                int right_cell[3];
                int lower_bounds[3];
                int upper_bounds[3];

                // Set up indices
                index[0] = i; index[1] = j; index[2] = k;
                left_cell[0] = i; left_cell[1] = j; left_cell[2] = k;
                right_cell[0] = i; right_cell[1] = j; right_cell[2] = k;
                
                // Define bounds including ghost cells
                lower_bounds[0] = bx_with_ghosts.smallEnd(0); 
                lower_bounds[1] = bx_with_ghosts.smallEnd(1); 
                lower_bounds[2] = AMREX_SPACEDIM == 3 ? bx_with_ghosts.smallEnd(2) : 0;
                upper_bounds[0] = bx_with_ghosts.bigEnd(0);
                upper_bounds[1] = bx_with_ghosts.bigEnd(1);
                upper_bounds[2] = AMREX_SPACEDIM == 3 ? bx_with_ghosts.bigEnd(2) : 0;
                
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
                    // First-order reconstruction simply uses the cell values directly
                    LeftStates_arr(index[0], index[1], index[2], n) = 
                        Var_arr(left_cell[0], left_cell[1], left_cell[2], n);
                    
                    RightStates_arr(index[0], index[1], index[2], n) = 
                        Var_arr(right_cell[0], right_cell[1], right_cell[2], n);
                }
            });
        }
    }
};

// Local Lax-Friedrichs flux method
template <typename IntegratorDerivedClass>
class LocalLaxFriedrichsMethod : public FluxMethod<IntegratorDerivedClass> {
public:
   using typename FluxMethod<IntegratorDerivedClass>::ReconstructionTarget;

      // Indicate that LLF works with flux-based reconstruction
    ReconstructionTarget getReconstructionTarget() const override {
        return ReconstructionTarget::Fluxes;
    }

    void SplitFluxes(int direction, int lev, IntegratorDerivedClass* solver, amrex::MultiFab& VariableBuffer, amrex::MultiFab& CellFluxBuffer, amrex::MultiFab& Flux_Negative, amrex::MultiFab& Flux_Positive) const override {

        const Set::Scalar gamma = 1.4;  // Ratio of specific heats

        const int nghosts = solver->number_of_ghost_cells;

        //for (amrex::MFIter mfi(QL_stencil, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        for (amrex::MFIter mfi(Flux_Negative, false); mfi.isValid(); ++mfi) {

            const amrex::Box& bx_with_ghosts = mfi.growntilebox();

           // auto const& flux_arr = (direction == FluxHandler<IntegratorDerivedClass>::Xdir) ? solver->XFlux_mf.Patch(lev, mfi) :
           //                         (direction == FluxHandler<IntegratorDerivedClass>::Ydir) ? solver->YFlux_mf.Patch(lev, mfi) :
           //                         solver->ZFlux_mf.Patch(lev, mfi);
            auto const& FluxNeg_arr = Flux_Negative.array(mfi);
            auto const& FluxPos_arr = Flux_Positive.array(mfi);
            auto const& p_arr  = solver->PVec_mf.Patch(lev,mfi);
            auto const& pres_arr  = solver->Pressure_mf.Patch(lev,mfi);
            auto const& q_arr  = VariableBuffer.array(mfi);
            auto const& flux_arr  = CellFluxBuffer.array(mfi);

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

        amrex::ParallelFor(bx_with_ghosts, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {


            // Construct index arrays with correct number of components
            int index[3];
            int left_index[3];
            int right_index[3];
            int lower_bounds[3];
            int upper_bounds[3];
        
            // Set up indices
            index[0] = i; index[1] = j; index[2] = k;
            left_index[0] = i; left_index[1] = j; left_index[2] = k;
            right_index[0] = i; right_index[1] = j; right_index[2] = k;
        
            // Define bounds including ghost cells
            lower_bounds[0] = bx_with_ghosts.smallEnd(0); 
            lower_bounds[1] = bx_with_ghosts.smallEnd(1); 
            lower_bounds[2] = AMREX_SPACEDIM == 3 ? bx_with_ghosts.smallEnd(2) : 0;
            upper_bounds[0] = bx_with_ghosts.bigEnd(0);
            upper_bounds[1] = bx_with_ghosts.bigEnd(1);
            upper_bounds[2] = AMREX_SPACEDIM == 3 ? bx_with_ghosts.bigEnd(2) : 0;
        
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
        
            // **Left and right primitive variables (QL and QR)**
            Set::Scalar pressure_L = pres_arr(left_index[0], left_index[1], left_index[2]);
            Set::Scalar pressure_R = pres_arr(right_index[0], right_index[1], right_index[2]);

            Set::Scalar vel_L[3] = {0.0, 0.0, 0.0};
            Set::Scalar vel_R[3] = {0.0, 0.0, 0.0};

            vel_L[normal] = p_arr(left_index[0], left_index[1], left_index[2], velocity_component[normal]);
            vel_R[normal] = p_arr(right_index[0], right_index[1], right_index[2], velocity_component[normal]);
        
            // **Speed of sound**
            Set::Scalar c_L = std::sqrt(gamma * std::max(pressure_L, 1e-8) / std::max(rho_L, 1e-8));
            Set::Scalar c_R = std::sqrt(gamma * std::max(pressure_R, 1e-8) / std::max(rho_R, 1e-8));
        
            // **Maximum wave speed (lambda_max)**
            Set::Scalar lambda_max = std::max(std::abs(vel_L[normal]) + c_L, std::abs(vel_R[normal]) + c_R);
        
            // **Flux arrays**
            int num_components = solver->number_of_components;
            std::vector<Set::Scalar> FL(num_components, 0.0);
            std::vector<Set::Scalar> FR(num_components, 0.0);

            // Mass flux
            FL[solver->variableIndex.DENS] = flux_arr(left_index[0], left_index[1], left_index[2], solver->variableIndex.DENS);
            FR[solver->variableIndex.DENS] = flux_arr(right_index[0], right_index[1], right_index[2], solver->variableIndex.DENS);
        
            // Momentum fluxes
            FL[solver->variableIndex.UVEL + normal] = flux_arr(left_index[0], left_index[1], left_index[2], solver->variableIndex.UVEL + normal);
            FR[solver->variableIndex.UVEL + normal] = flux_arr(right_index[0], right_index[1], right_index[2], solver->variableIndex.UVEL + normal);

#if AMREX_SPACEDIM >= 2
            FL[solver->variableIndex.UVEL + trans1] = flux_arr(left_index[0], left_index[1], left_index[2], solver->variableIndex.UVEL + trans1);
            FR[solver->variableIndex.UVEL + trans1] = flux_arr(right_index[0], right_index[1], right_index[2], solver->variableIndex.UVEL + trans1);
#endif
#if AMREX_SPACEDIM == 3
            FL[solver->variableIndex.UVEL + trans2] = flux_arr(left_index[0], left_index[1], left_index[2], solver->variableIndex.UVEL + trans2);
            FR[solver->variableIndex.UVEL + trans2] = flux_arr(right_index[0], right_index[1], right_index[2], solver->variableIndex.UVEL + trans2);
#endif
                    
            // Energy flux
            FL[solver->variableIndex.IE] = flux_arr(left_index[0], left_index[1], left_index[2], solver->variableIndex.IE);
            FR[solver->variableIndex.IE] = flux_arr(right_index[0], right_index[1], right_index[2], solver->variableIndex.IE);
        

            // **Conservative state vectors (left and right)**
            std::vector<Set::Scalar> QL_cons(num_components, 0.0);
            std::vector<Set::Scalar> QR_cons(num_components, 0.0);
            
            // **Convert primitive variables to conservative variables**
            QL_cons[solver->variableIndex.DENS] = q_arr(left_index[0], left_index[1], left_index[2], solver->variableIndex.DENS); //rho_L;  // Density
            QR_cons[solver->variableIndex.DENS] = q_arr(right_index[0], right_index[1], right_index[2], solver->variableIndex.DENS); //rho_R;
            
            QL_cons[solver->variableIndex.UVEL + normal] = q_arr(left_index[0], left_index[1], left_index[2], solver->variableIndex.UVEL + normal);//rho_L * vel_L[normal];  // Momentum in normal direction
            QR_cons[solver->variableIndex.UVEL + normal] = q_arr(right_index[0], right_index[1], right_index[2], solver->variableIndex.UVEL + normal); //rho_R * vel_R[normal];
            
#if AMREX_SPACEDIM >= 2
            QL_cons[solver->variableIndex.UVEL + trans1] = q_arr(left_index[0], left_index[1], left_index[2], solver->variableIndex.UVEL + trans1);//rho_L * vel_L[trans1];  // Momentum in first transverse direction
            QR_cons[solver->variableIndex.UVEL + trans1] = q_arr(right_index[0], right_index[1], right_index[2], solver->variableIndex.UVEL + trans1); //rho_R * vel_R[trans1];
#endif
            
#if AMREX_SPACEDIM == 3
            QL_cons[solver->variableIndex.UVEL + trans2] = q_arr(left_index[0], left_index[1], left_index[2], solver->variableIndex.UVEL + trans2);//rho_L * vel_L[trans2];  // Momentum in second transverse direction
            QR_cons[solver->variableIndex.UVEL + trans2] = q_arr(right_index[0], right_index[1], right_index[2], solver->variableIndex.UVEL + trans2); //rho_R * vel_R[trans2];
#endif
            
            QL_cons[solver->variableIndex.IE] = q_arr(left_index[0], left_index[1], left_index[2], solver->variableIndex.IE);//rho_L * E_L;  // Total energy
            QR_cons[solver->variableIndex.IE] = q_arr(right_index[0], right_index[1], right_index[2], solver->variableIndex.IE); //rho_R * E_R;
            

            // **Local Lax-Friedrichs flux calculation**
            for (int n = 0; n < num_components; ++n) {
              FluxNeg_arr(index[0], index[1], index[2], n) =  0.5 * FL[n] + 0.5 * lambda_max * QL_cons[n];
              FluxPos_arr(index[0], index[1], index[2], n) =  0.5 * FR[n] - 0.5 * lambda_max * QR_cons[n];
            }
        
        });
        
        
        }
        
    }

    // Implement ComputeFluxes for completeness (though it won't be used for LLF)
    void ComputeFluxes(int direction, int lev, IntegratorDerivedClass* solver,
                      const amrex::MultiFab& LeftStates, const amrex::MultiFab& RightStates,
                      amrex::MultiFab& SummedFlux) const override {
        // This is not the primary pathway for LLF, but provide an implementation
        // for completeness or potential future use
        Util::Warning(INFO, "LLF ComputeFluxes called directly - not optimal path");
        
        // If needed, a direct implementation could be added here
    }    

};

// HLLC flux method
template <typename IntegratorDerivedClass>
class HLLCMethod : public FluxMethod<IntegratorDerivedClass> {
public:
    
   using typename FluxMethod<IntegratorDerivedClass>::ReconstructionTarget;

   ReconstructionTarget getReconstructionTarget() const override {
       return ReconstructionTarget::Variables;
   }
   
   void SplitFluxes(int direction, int lev, IntegratorDerivedClass* solver, 
                   amrex::MultiFab& VariableBuffer, amrex::MultiFab& CellFluxBuffer,
                   amrex::MultiFab& Flux_Negative, amrex::MultiFab& Flux_Positive) const override {
       // Not the primary pathway for HLLC
       Util::Warning(INFO, "HLLC SplitFluxes called - not optimal path");
   }
   
   void ComputeFluxes(int direction, int lev, IntegratorDerivedClass* solver,
                  const amrex::MultiFab& LeftStates, const amrex::MultiFab& RightStates,
                  amrex::MultiFab& SummedFlux) const override {
    
        const Set::Scalar gamma = 1.4;  // Ratio of specific heats
        const int nghosts = solver->number_of_ghost_cells;

        // Loop through the grid
        for (amrex::MFIter mfi(SummedFlux, false); mfi.isValid(); ++mfi) {
            const amrex::Box& bx_with_ghosts = mfi.growntilebox();


            // Updated array references for LeftStates, RightStates, and SummedFlux
            auto const& Left_arr = LeftStates.array(mfi);
            auto const& Right_arr = RightStates.array(mfi);
            auto const& TotalFlux_arr = SummedFlux.array(mfi);
            auto const& pres_arr = solver->Pressure_mf.Patch(lev, mfi);
                  
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

            int num_components = solver->number_of_components;

            // Parallel loop over cells
            amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Setup indices
                int index[3];
                int lower_bounds[3];
                int upper_bounds[3];
                
                // Initialize index arrays
                index[0] = i; index[1] = j; index[2] = k;
                
                // Define bounds including ghost cells
                lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); 
                lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); 
                lower_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.smallEnd(2) : 0;
                upper_bounds[0] = face_bx_with_ghosts.bigEnd(0) - 1;
                upper_bounds[1] = face_bx_with_ghosts.bigEnd(1) - 1;
                upper_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.bigEnd(2) - 1 : 0;
                
                // Clamp indices to ensure they're within bounds
                ClampIndices(index, lower_bounds, upper_bounds);
                                
                // Get primitive variables for left and right states
                // Left state
                Set::Scalar rho_L = Left_arr(index[0], index[1], index[2], solver->variableIndex.DENS);
                rho_L = std::max(rho_L, 1e-8);  // Avoid zero density
                
                Set::Scalar e_L = Left_arr(index[0], index[1], index[2], solver->variableIndex.IE);
                Set::Scalar pressure_L = pres_arr(index[0], index[1], index[2]);
                pressure_L = std::max(pressure_L, 1e-8);  // Avoid zero pressure
                
                Set::Scalar vel_L[3] = {0.0, 0.0, 0.0};
                vel_L[normal] = Left_arr(index[0], index[1], index[2], velocity_component[normal]);
#if AMREX_SPACEDIM >= 2    
                vel_L[trans1] = Left_arr(index[0], index[1], index[2], velocity_component[trans1]);
#endif
#if AMREX_SPACEDIM == 3    
                vel_L[trans2] = Left_arr(index[0], index[1], index[2], velocity_component[trans2]);
#endif
                
                // Right state
                Set::Scalar rho_R = Right_arr(index[0], index[1], index[2], solver->variableIndex.DENS);
                rho_R = std::max(rho_R, 1e-8);  // Avoid zero density
                
                Set::Scalar e_R = Right_arr(index[0], index[1], index[2], solver->variableIndex.IE);
                Set::Scalar pressure_R = pres_arr(index[0], index[1], index[2]);
                pressure_R = std::max(pressure_R, 1e-8);  // Avoid zero pressure
                
                Set::Scalar vel_R[3] = {0.0, 0.0, 0.0};
                vel_R[normal] = Right_arr(index[0], index[1], index[2], velocity_component[normal]);
#if AMREX_SPACEDIM >= 2    
                vel_R[trans1] = Right_arr(index[0], index[1], index[2], velocity_component[trans1]);
#endif
#if AMREX_SPACEDIM == 3    
                vel_R[trans2] = Right_arr(index[0], index[1], index[2], velocity_component[trans2]);
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
                
                // Conservative variables for left and right states
                std::vector<Set::Scalar> QL_cons(num_components, 0.0);
                std::vector<Set::Scalar> QR_cons(num_components, 0.0);
                
                // Convert primitive variables to conservative variables
                QL_cons[solver->variableIndex.DENS] = rho_L;
                QL_cons[solver->variableIndex.UVEL + normal] = rho_L * vel_L[normal];
#if AMREX_SPACEDIM >= 2
                QL_cons[solver->variableIndex.UVEL + trans1] = rho_L * vel_L[trans1];
#endif
#if AMREX_SPACEDIM == 3
                QL_cons[solver->variableIndex.UVEL + trans2] = rho_L * vel_L[trans2];
#endif
                QL_cons[solver->variableIndex.IE] = rho_L * E_L;
                
                QR_cons[solver->variableIndex.DENS] = rho_R;
                QR_cons[solver->variableIndex.UVEL + normal] = rho_R * vel_R[normal];
#if AMREX_SPACEDIM >= 2
                QR_cons[solver->variableIndex.UVEL + trans1] = rho_R * vel_R[trans1];
#endif
#if AMREX_SPACEDIM == 3
                QR_cons[solver->variableIndex.UVEL + trans2] = rho_R * vel_R[trans2];
#endif
                QR_cons[solver->variableIndex.IE] = rho_R * E_R;
                
                // Compute physical fluxes for left and right states
                std::vector<Set::Scalar> FL(num_components, 0.0);
                std::vector<Set::Scalar> FR(num_components, 0.0);
                
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
                
                // Wave speeds using Davis estimates (more robust than Roe averages)
                Set::Scalar SL = std::min(vel_L[normal] - c_L, vel_R[normal] - c_R);
                Set::Scalar SR = std::max(vel_L[normal] + c_L, vel_R[normal] + c_R);
                
                Set::Scalar epsilon = 0.1 * std::min(c_L, c_R); 
                if (std::abs(SR - SL) < epsilon) {
                    //Apply entropy fix by adjusting wave speeds
                    SR += epsilon;
                    SL -= epsilon;    
                }   

                // Calculate the middle wave speed (S*)
                Set::Scalar numerator = pressure_R - pressure_L + 
                                     rho_L * vel_L[normal] * (SL - vel_L[normal]) - 
                                     rho_R * vel_R[normal] * (SR - vel_R[normal]);
                                     
                Set::Scalar denominator = rho_L * (SL - vel_L[normal]) - rho_R * (SR - vel_R[normal]);
                
                // Prevent division by zero with a safeguard
                if (std::abs(denominator) < 1e-10) {
                    denominator = (denominator >= 0) ? 1e-10 : -1e-10;
                }
                
                Set::Scalar S_star = numerator / denominator;
                
                // Compute star state variables
                std::vector<Set::Scalar> QL_star(num_components, 0.0);
                std::vector<Set::Scalar> QR_star(num_components, 0.0);
                
                // Density factors for star regions
                Set::Scalar rho_L_factor = (SL - vel_L[normal]) / (SL - S_star);
                Set::Scalar rho_R_factor = (SR - vel_R[normal]) / (SR - S_star);
                
                // Left star state
                QL_star[solver->variableIndex.DENS] = rho_L * rho_L_factor;
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
                QR_star[solver->variableIndex.DENS] = rho_R * rho_R_factor;
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
                std::vector<Set::Scalar> FL_star(num_components, 0.0);
                std::vector<Set::Scalar> FR_star(num_components, 0.0);
                
                for (int n = 0; n < num_components; ++n) {
                    FL_star[n] = FL[n] + SL * (QL_star[n] - QL_cons[n]);
                    FR_star[n] = FR[n] + SR * (QR_star[n] - QR_cons[n]);
                }
                
                // Compute branch-free HLLC flux
                Set::Scalar HL = 0.5 * (1.0 + std::copysign(1.0, SL));  // H(SL)
                Set::Scalar HS = 0.5 * (1.0 + std::copysign(1.0, S_star));  // H(S*)
                Set::Scalar HR = 0.5 * (1.0 + std::copysign(1.0, SR));  // H(SR)
                

                for (int n = 0; n < num_components; ++n) {
                    // Branch-free HLLC flux calculation
                    Set::Scalar HLLC_flux = 
                        HL * FL[n] + 
                        (1.0 - HL) * HS * FL_star[n] + 
                        (1.0 - HS) * HR * FR_star[n] + 
                        (1.0 - HR) * FR[n];

                    TotalFlux_arr(index[0], index[1], index[2], n) = HLLC_flux;                    
                }
            });
        }
    }
};

}  // namespace Numeric

#endif  // NUMERIC_FLUXHANDLER_H

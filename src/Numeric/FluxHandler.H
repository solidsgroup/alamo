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
#include "Numeric/SymmetryPreservingRoeAveragingOperations.H"
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

    // Original method for flux-based reconstruction (for LLF-like methods)
    virtual void ReconstructFluxes(int direction, int lev, IntegratorDerivedClass* solver,
                                const amrex::MultiFab& Flux_Negative, 
                                const amrex::MultiFab& Flux_Positive,
                                amrex::MultiFab& SummedFlux,
                                const std::shared_ptr<GenericVariableAccessor>& variable_accessor,
                                const amrex::MultiFab* AverageStateBuffer = nullptr) const = 0;

    // New method for variable-based reconstruction (for HLLC-like methods)
    virtual void ReconstructVariables(int direction, int lev, IntegratorDerivedClass* solver,
                                    const amrex::MultiFab& VariableBuffer,
                                    amrex::MultiFab& LeftStates, 
                                    amrex::MultiFab& RightStates,
                                    const std::shared_ptr<GenericVariableAccessor>& variable_accessor,
                                    const amrex::MultiFab* AverageStateBuffer = nullptr) const = 0;
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

    static inline void ProbeCellCC(
        char const* tag,
        amrex::MultiFab const& mf,
        int lev, int i, int j, int k,
        int comp)
    {
        amrex::Gpu::streamSynchronize();
        amrex::IntVect iv(AMREX_D_DECL(i,j,k));
    
        for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {      
            const amrex::Box& vb = mfi.validbox();
            if (!vb.contains(iv)) continue;
    
            auto const& a = mf.const_array(mfi);
            amrex::Print() << tag
                           << " lev="<<lev
                           << " iv=("<<i<<","<<j<<","<<k<<")"
                           << " valid=1 vb="<<vb
                           << " comp="<<comp
                           << " val="<<a(i,j,k,comp)
                           << "\n";
            return;
        }
    
        amrex::Print() << tag << " lev="<<lev
                       << " iv=("<<i<<","<<j<<","<<k<<") valid=0\n";
    }
    
    static inline void ProbeFace(
        char const* tag,
        amrex::MultiFab const& mf_face,
        int lev, int dir, int i, int j, int k,
        int comp)
    {
        amrex::Gpu::streamSynchronize();
        amrex::IntVect iv(AMREX_D_DECL(i,j,k));
    
        for (amrex::MFIter mfi(mf_face); mfi.isValid(); ++mfi) { // NOTE: no tiling here
            const amrex::Box& vb = mfi.validbox();
            if (!vb.contains(iv)) continue;
    
            auto const& a = mf_face.const_array(mfi);
            amrex::Print() << tag
                           << " lev="<<lev
                           << " dir="<<dir
                           << " face=("<<i<<","<<j<<","<<k<<")"
                           << " valid=1 vb="<<vb
                           << " comp="<<comp
                           << " val="<<a(i,j,k,comp)
                           << "\n";
            return;
        }
    
        amrex::Print() << tag
                       << " lev="<<lev
                       << " dir="<<dir
                       << " face=("<<i<<","<<j<<","<<k<<") valid=0\n";
    }
    

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


        // Process each direction
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            // Choose processing path based on the flux method's reconstruction target
            if (flux_method->getReconstructionTarget() == FluxMethod<IntegratorDerivedClass>::ReconstructionTarget::Variables) {
                // For methods like HLLC that reconstruct variables first
                ProcessVariableBasedReconstruction(dir, lev, solver);
            } else {
                // For methods like LLF that reconstruct fluxes first
                ProcessFluxBasedReconstruction(dir, lev, solver);
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
        int direction, int lev, IntegratorDerivedClass* solver) const;

    // For flux-based reconstruction methods like LLF
    void ProcessFluxBasedReconstruction(
        int direction, int lev, IntegratorDerivedClass* solver) const;
};

// For variable-based reconstruction methods like HLLC
template <typename IntegratorDerivedClass>
void FluxHandler<IntegratorDerivedClass>::ProcessVariableBasedReconstruction(
    int direction, int lev, IntegratorDerivedClass* solver) const {
    Util::ScimitarX_Util::Debug::SetTargetDebugLocationIndices(32, 7, 7, "true");

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

    int ghost_cells = solver->number_of_ghost_cells;


    // 1. Create intermediate buffers for storing cell centered variables
    amrex::MultiFab VariableBuffer = variable_accessor->CreateWorkingBuffer(
        solver->grids[lev], solver->dmap[lev], solver->number_of_components, ghost_cells
    );


    // 2. Transform variables if needed (based on reconstruction mode)
    variable_accessor->CopyVariables(direction, lev, solver, VariableBuffer, 
                                        mode, validationResult);

    // if (direction == 0) {
    // const int rho = solver->variableIndex.DENS;
    // ProbeCellCC("[PROBE] PRE-RECON VBUF rho @ (31,32,8)", VariableBuffer, lev, 31,32,8, rho);
    // ProbeCellCC("[PROBE] PRE-RECON VBUF rho @ (32,32,8)", VariableBuffer, lev, 32,32,8, rho);
    // ProbeCellCC("[PROBE] PRE-RECON VBUF rho @ (33,32,8)", VariableBuffer, lev, 33,32,8, rho);
    // }



    // if (direction == 0) {
    //     const int i=32, j=32, k=8;
    //     const int dens = solver->variableIndex.DENS;


    //     for (amrex::MFIter mfi(VariableBuffer); mfi.isValid(); ++mfi) {
    //         const amrex::Box& bx = VariableBuffer[mfi].box();   // use FAB box, not tilebox
    //         if (!bx.contains(amrex::IntVect(i-1,j,k)) ||
    //             !bx.contains(amrex::IntVect(i,j,k)) ||
    //             !bx.contains(amrex::IntVect(i+1,j,k))) continue;

    //         auto const& V = VariableBuffer.const_array(mfi);

    //         amrex::AllPrint()
    //             << "\n[PRE-RECON BUF] stage=" 
    //             << " cell(31,32,8) DENS=" << V(i-1,j,k,dens)
    //             << " cell(32,32,8) DENS=" << V(i  ,j,k,dens)
    //             << " cell(33,32,8) DENS=" << V(i+1,j,k,dens)
    //             << "\n";
    //         break;
    //     }
    // }



    
    // 4. Create a face-centered BoxArray for the flux
    amrex::BoxArray face_ba = solver->grids[lev];
    amrex::IntVect face_type = amrex::IntVect::TheZeroVector();
    face_type[direction] = 1;  // Set face-centered in the specified direction
    face_ba.convert(face_type);


    amrex::MultiFab AverageStateBuffer;

    // 3. If Characteristic Reconstruction, Populate Average States in a buffer
    if ( mode == ReconstructionMode::Characteristic) {
    // 3a. Create intermediate buffers for storing cell centered variables
        AverageStateBuffer = variable_accessor->CreateWorkingBuffer(
        face_ba, solver->dmap[lev], solver->number_of_components, 0
    );

    // 3b. Compute average state buffers
    variable_accessor->PopulateAverageStates(direction, lev, solver, AverageStateBuffer);  

    }

    amrex::MultiFab LeftStates = variable_accessor->CreateWorkingBuffer(
        face_ba, solver->dmap[lev], solver->number_of_components, 0);

    amrex::MultiFab RightStates = variable_accessor->CreateWorkingBuffer(
        face_ba, solver->dmap[lev], solver->number_of_components, 0);

    // 5. Perform reconstruction on variables to get left/right states
    reconstruction->ReconstructVariables(
            direction, lev, solver, 
            VariableBuffer, LeftStates, RightStates, 
            variable_accessor,
            (mode == ReconstructionMode::Characteristic) ? &AverageStateBuffer : nullptr);


    // if (direction == 0) {
    // const int rho = solver->variableIndex.DENS;
    // ProbeFace("[PROBE] POST-RECON L rho @ face(32,32,8)", LeftStates,  lev, 0, 32,32,8, rho);
    // ProbeFace("[PROBE] POST-RECON R rho @ face(32,32,8)", RightStates, lev, 0, 32,32,8, rho);
    // }


    // if (direction == 0) {
    //     const int pi = 32, pj = 32, pk = 8;
    //     for (amrex::MFIter mfi(LeftStates, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
    //         const amrex::Box& bx = mfi.tilebox();

    //         if (bx.contains(amrex::IntVect(pi, pj, pk))) {
    //             auto L = LeftStates.array(mfi);
    //             auto R = RightStates.array(mfi);

    //             amrex::Print() << "PROBE1 Left/Right DENS @ (32,32,8): "
    //                         << L(pi,pj,pk, solver->variableIndex.DENS) << "  "
    //                         << R(pi,pj,pk, solver->variableIndex.DENS) << "\n";

    //             // amrex::Print() << "PROBE1 Left/Right DENS @ (33,32,8): "
    //             //             << L(pi+1,pj,pk, solver->variableIndex.DENS) << "  "
    //             //             << R(pi+1,pj,pk, solver->variableIndex.DENS) << "\n";
    //             break;
    //         }
    //     }
    // }

    
    //     const int ii = 32, jj = 32, kk = 8;
    //     if (direction == 0) {
    //         for (amrex::MFIter mfi(LeftStates, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    //             const amrex::Box& tbx = mfi.tilebox();
    //             if (!tbx.contains(amrex::IntVect(ii,jj,kk))) continue;
        
    //             auto const& LS = LeftStates.const_array(mfi);
    //             auto const& RS = RightStates.const_array(mfi);

    //                 amrex::Print() << "\n[POST-RECON] dir=" << direction
    //                             << " face=(" << ii << "," << jj << "," << kk << ")"
    //                             << " tilebox=" << tbx << "\n";


    //                 for (int n = 0; n < solver->number_of_components; ++n) {
    //                     amrex::Print() << "  n=" << n
    //                                 << " Left State=" << LS(ii,jj,kk,n)
    //                                 << " Right State=" << RS(ii,jj,kk,n) << "\n";
    //                 }
    //         }

    //     }

    // 6. Create MultiFab for the computed flux
    amrex::MultiFab SummedFlux(face_ba, solver->dmap[lev], 
                            solver->number_of_components, ghost_cells);

    // 7. Compute fluxes from reconstructed states
    flux_method->ComputeFluxes(direction, lev, solver, 
                            LeftStates, RightStates, SummedFlux);

    // if (direction == 0) {
    // const int comp_to_watch = 0; // pick the comp you know is nonzero (your log shows n=3 is 0.1)
    // ProbeFace("[PROBE] POST-FLUX SF comp0 @ face(32,32,8)", SummedFlux, lev, 0, 32,32,8, comp_to_watch);
    // }




                        
    // const int iii=32, jjj=32, kkk=8;

    // // 1) Point probe
    // if (direction == 0) {
    //     for (amrex::MFIter mfi(SummedFlux); mfi.isValid(); ++mfi) {
    //         const amrex::Box& vbx = mfi.validbox();
    //         if (!vbx.contains(amrex::IntVect(AMREX_D_DECL(iii,jjj,kkk)))) continue;

    //         auto const& F = SummedFlux.const_array(mfi);

    //         amrex::Print()
    //         << "\n[POST-FLUX POINT] dir=" << direction
    //         << " F("<<iii<<","<<jjj<<","<<kkk<<",:)\n";
    //         for (int n=0; n<solver->number_of_components; ++n) {
    //             amrex::Print() << "  n="<<n<<" F="<<F(iii,jjj,kkk,n) << "\n";
    //         }
    //         break;
    //     }
    // }


    // 7. Store fluxes in the appropriate directional flux arrays
    variable_accessor->StoreDirectionalFlux(direction, lev, solver, 
                                        SummedFlux);


    // if (direction == 0) {
    // const int comp_to_watch = 0;
    // ProbeFace("[PROBE] POST-STORE SF comp0 @ face(32,32,8)", SummedFlux, lev, 0, 32,32,8, comp_to_watch);
    // ProbeFace("[PROBE] POST-STORE XF comp0 @ face(32,32,8)", *solver->XFlux_mf[lev], lev, 0, 32,32,8, comp_to_watch);
    // }


    

    // if (direction == 0)
    // {
    //     const int ii=32, jj=32, kk=8;
    //     const int n = 0;

    //     amrex::MultiFab& YF = *solver->YFlux_mf[lev]; 
    //     amrex::MultiFab& XF = *solver->XFlux_mf[lev]; 
    //     amrex::MultiFab& ZF = *solver->ZFlux_mf[lev];   // unique_ptr -> deref


    //     for (amrex::MFIter mfi(SummedFlux); mfi.isValid(); ++mfi) {
    //         const amrex::Box& bx = SummedFlux[mfi].box();
    //         if (!bx.contains(amrex::IntVect(AMREX_D_DECL(ii,jj,kk)))) continue;

    //         auto const& SF = SummedFlux.const_array(mfi);
    //         auto const& DF2 = YF.const_array(mfi);
    //         auto const& DF1 = XF.const_array(mfi);
    //         auto const& DF3 = ZF.const_array(mfi);

    //         amrex::Real a = SF(ii,jj,kk,n);
    //         amrex::Real b = DF1(ii,jj,kk,n);
    //         amrex::Real c = DF2(ii,jj,kk,n);
    //         amrex::Real d = DF3(ii,jj,kk,n);

    //         amrex::Print() << "[POST-STORE] stage Y="
    //                     << " face=("<<ii<<","<<jj<<","<<kk<<") comp="<<n
    //                     << " SF="<<a<<" XF="<<b<<" YF="<<c << "ZF="<<d<<"\n";
    //         break;
    //     }
    // }
  


                                        
}

// For flux-based reconstruction methods like LLF
template <typename IntegratorDerivedClass>
void FluxHandler<IntegratorDerivedClass>::ProcessFluxBasedReconstruction(
    int direction, int lev, IntegratorDerivedClass* solver) const {

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


    int ghost_cells = solver->number_of_ghost_cells;


    // 1. Create intermediate buffers for storing cell centered variables
    amrex::MultiFab VariableBuffer = variable_accessor->CreateWorkingBuffer(
        solver->grids[lev], solver->dmap[lev], solver->number_of_components, ghost_cells
    );

    // 2. Transform variables if needed (based on reconstruction mode)
    variable_accessor->CopyVariables(direction, lev, solver, VariableBuffer, 
                                        mode, validationResult);


    // 4. Create a cell flux buffer for storing computed fluxes
    amrex::MultiFab CellFluxBuffer = variable_accessor->CreateWorkingBuffer(
        solver->grids[lev], solver->dmap[lev], solver->number_of_components, 
        solver->number_of_ghost_cells);

    // 5. Compute and transform fluxes if needed
    variable_accessor->CopyFluxes(direction, lev, solver, CellFluxBuffer, 
                                    mode, validationResult);

    // 6. Create buffers for split fluxes (positive and negative)
    amrex::MultiFab Flux_Negative = variable_accessor->CreateWorkingBuffer(
        solver->grids[lev], solver->dmap[lev], solver->number_of_components, ghost_cells);

    amrex::MultiFab Flux_Positive = variable_accessor->CreateWorkingBuffer(
        solver->grids[lev], solver->dmap[lev], solver->number_of_components, ghost_cells);

    // 7. Split fluxes into positive and negative components
    flux_method->SplitFluxes(direction, lev, solver, VariableBuffer, 
                            CellFluxBuffer, Flux_Negative, Flux_Positive);

    // 8. Create a face-centered BoxArray for the flux
    amrex::BoxArray face_ba = solver->grids[lev];
    amrex::IntVect face_type = amrex::IntVect::TheZeroVector();
    face_type[direction] = 1;  // Set face-centered in the specified direction
    face_ba.convert(face_type);


    amrex::MultiFab AverageStateBuffer;

    // 3. If Characteristic Reconstruction, Populate Average States in a buffer
    if (mode == ReconstructionMode::Characteristic) {
    // 3a. Create intermediate buffers for storing cell centered variables
    AverageStateBuffer = variable_accessor->CreateWorkingBuffer(
        face_ba, solver->dmap[lev], solver->number_of_components, ghost_cells);

    // 3b. Compute average state buffers
    variable_accessor->PopulateAverageStates(direction, lev, solver, AverageStateBuffer);  

    }

    // 9. Create MultiFab for the reconstructed summed flux
    amrex::MultiFab SummedFlux(face_ba, solver->dmap[lev], 
                            solver->number_of_components, ghost_cells);

    // 10. Perform reconstruction on fluxes
    reconstruction->ReconstructFluxes(
            direction, lev, solver, 
            Flux_Negative, Flux_Positive, SummedFlux, 
            variable_accessor,
            (mode == ReconstructionMode::Characteristic) ? &AverageStateBuffer : nullptr);

    // 11. Store fluxes in the appropriate directional flux arrays
    variable_accessor->StoreDirectionalFlux(direction, lev, solver, 
                                        SummedFlux);

}

// First-order upwind reconstruction
template <typename IntegratorDerivedClass>
class FirstOrderReconstruction : public FluxReconstruction<IntegratorDerivedClass> {
public:
    // Implement flux-based reconstruction (for LLF-like methods)
    void ReconstructFluxes(int direction, int /*lev*/, IntegratorDerivedClass* solver,
                        const amrex::MultiFab& Flux_Negative, const amrex::MultiFab& Flux_Positive,
                        amrex::MultiFab& SummedFlux,
                        [[maybe_unused]] const std::shared_ptr<GenericVariableAccessor>& variable_accessor,
                        [[maybe_unused]] const amrex::MultiFab* AverageStateBuffer = nullptr) const override {

        const int nghosts = solver->number_of_ghost_cells; 

        // Loop through face-centered indices for flux reconstruction
        for (amrex::MFIter mfi(SummedFlux, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
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
    void ReconstructVariables(int direction, int /*lev*/, IntegratorDerivedClass* solver,
                            const amrex::MultiFab& VariableBuffer,
                            amrex::MultiFab& LeftStates, amrex::MultiFab& RightStates,
                        [[maybe_unused]] const std::shared_ptr<GenericVariableAccessor>& variable_accessor,
                        [[maybe_unused]] const amrex::MultiFab* AverageStateBuffer = nullptr) const override {

        // const int nghosts = solver->number_of_ghost_cells;

        // Loop through grid cells to reconstruct left/right states at interfaces
        for (amrex::MFIter mfi(LeftStates, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& face_bx = mfi.tilebox();
            // const amrex::Box& vbx = VariableBuffer[mfi].box();

            auto const& Var_arr = VariableBuffer.array(mfi);
            auto const& LeftStates_arr = LeftStates.array(mfi);
            auto const& RightStates_arr = RightStates.array(mfi);

            const int num_components = solver->number_of_components;
            
            // if (amrex::ParallelDescriptor::MyProc() == 0 && direction == 2) {

            //     const amrex::Box tbx  = mfi.tilebox();
            //     const amrex::Box ntbx = mfi.nodaltilebox(direction);
            //     const amrex::Box gnbx = mfi.grownnodaltilebox(direction, nghosts);

            //     amrex::Print()
            //     << "TBX  " << bx_with_ghosts << " size=(" << bx_with_ghosts.length(0)  << "," << bx_with_ghosts.length(1)  << "," << bx_with_ghosts.length(2)  << ")\n"
            //     << "ValidBox = " << vbx   << " size=("<< vbx.length(0) << "," << vbx.length(1) << "," << vbx.length(2) << ")\n"
            //     // << "NTBX " << ntbx << " size=(" << ntbx.length(0) << "," << ntbx.length(1) << "," << ntbx.length(2) << ")\n"
            //     // << "GNBX " << gnbx << " size=(" << gnbx.length(0) << "," << gnbx.length(1) << "," << gnbx.length(2) << ")\n"
            //     << "----\n";
            // }

            amrex::ParallelFor(face_bx, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {
                // Construct index arrays
                
                // int index[3];
                // int left_index[3];
                // int right_index[3];
                // int lower_bounds[3];
                // int upper_bounds[3];

                // // Set up indices
                // index[0] = i; index[1] = j; index[2] = k;
                // left_index[0] = i; left_index[1] = j; left_index[2] = k;
                // right_index[0] = i; right_index[1] = j; right_index[2] = k;

                // // Define bounds including ghost cells
                // lower_bounds[0] = bx_with_ghosts.smallEnd(0); 
                // lower_bounds[1] = bx_with_ghosts.smallEnd(1); 
                // lower_bounds[2] = AMREX_SPACEDIM == 3 ? bx_with_ghosts.smallEnd(2) : 0;
                // upper_bounds[0] = bx_with_ghosts.bigEnd(0) - 1;
                // upper_bounds[1] = bx_with_ghosts.bigEnd(1) - 1;
                // upper_bounds[2] = AMREX_SPACEDIM == 3 ? bx_with_ghosts.bigEnd(2) - 1 : 0;

                // lower_bounds[0] = vbx.smallEnd(0); 
                // lower_bounds[1] = vbx.smallEnd(1); 
                // lower_bounds[2] = AMREX_SPACEDIM == 3 ? vbx.smallEnd(2): 0;
                // upper_bounds[0] = vbx.bigEnd(0);
                // upper_bounds[1] = vbx.bigEnd(1);
                // upper_bounds[2] = AMREX_SPACEDIM == 3 ? vbx.bigEnd(2) : 0;


                // // Clamp the main index
                // ClampIndices(index, lower_bounds, upper_bounds);

                // // Shift and clamp for left and right indices
                // int left_offset[3] = {0, 0, 0};
                // int right_offset[3] = {0, 0, 0};
                // left_offset[direction] = 0;  // Shift left
                // right_offset[direction] = 1;  // Shift right

                //     // Save before shift+clamp (so we can see what changed)
                // int left_index_pre[3]  = { left_index[0],  left_index[1],  left_index[2]  };
                // int right_index_pre[3] = { right_index[0], right_index[1], right_index[2] };


                // ShiftAndClampIndices(left_index, left_offset, lower_bounds, upper_bounds, direction);
                // ShiftAndClampIndices(right_index, right_offset, lower_bounds, upper_bounds, direction);

                // Create branchless masks for left (i-1) and right (i) indices for each face based on direction
                int dx = (direction == 0);
                int dy = (direction == 1);
                int dz = (direction == 2);

                // Compute left/right indices with no branching
                int iL = iface - dx;
                int iR = iface;

                int jL = jface - dy;
                int jR = jface;

                int kL = kface - dz;
                int kR = kface;

                if (Util::ScimitarX_Util::Debug::IsTargetLocation(iface, jface, kface))
                {
                    // 1) Print indices before/after
                    // Util::ScimitarX_Util::Debug::DebugIndicesAndBounds(
                    //     index[0], index[1], index[2],
                    //     index,
                    //     left_index_pre, right_index_pre,   
                    //     lower_bounds, upper_bounds,
                    //     "RECON: indices PRE shift+clamp", true);

                    // Util::ScimitarX_Util::Debug::DebugIndicesAndBounds(
                    //     index[0], index[1], index[2],
                    //     index,
                    //     left_index, right_index,           
                    //     lower_bounds, upper_bounds,
                    //     "RECON: indices POST shift+clamp", true);

                    const int rho = solver->variableIndex.DENS;
                    const int u   = solver->variableIndex.UVEL;   
                    const int v   = solver->variableIndex.VVEL;
                #if (AMREX_SPACEDIM == 3)
                    const int w   = solver->variableIndex.WVEL;
                #endif
                      

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        Var_arr(iL, jL, kL, rho),
                        "PVec[rho]@left_index", "RECON: source left", false, rho, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        Var_arr(iR, jR, kR, rho),
                        "PVec[rho]@right_index", "RECON: source right", false, rho, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        Var_arr(iL, jL, kL, u),
                        "PVec[u]@left_index", "RECON: source left", false, u, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        Var_arr(iR, jR, kR, u),
                        "PVec[u]@right_index", "RECON: source right", false, u, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        Var_arr(iL, jL, kL, v),
                        "PVec[v]@left_index", "RECON: source left", false, v, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        Var_arr(iR, jR, kR, v),
                        "PVec[v]@right_index", "RECON: source right", false, v, true);

                #if (AMREX_SPACEDIM == 3)
                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        Var_arr(iL, jL, kL, w),
                        "PVec[w]@left_index", "RECON: source left", false, w, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iL, jL, kL,
                        Var_arr(iR, jR, kR, w),
                        "PVec[w]@right_index", "RECON: source right", false, w, true);
                #endif

                   
                }


            
            
                for (int n = 0; n < num_components; ++n) {
                    LeftStates_arr(iface, jface, kface, n) =
                        Var_arr(iL, jL, kL, n);
                
                    RightStates_arr(iface, jface, kface, n) =
                        Var_arr(iR, jR, kR, n);
                }
                // === PROBE: after reconstruction write (face states) ===
                if (Util::ScimitarX_Util::Debug::IsTargetLocation(iface, jface, kface))
                {
                    // Util::ScimitarX_Util::Debug::DebugIndicesAndBounds(
                    //     index[0], index[1], index[2],
                    //     index, left_index, right_index,
                    //     lower_bounds, upper_bounds,
                    //     "RECON: indices/bounds (after clamp/shift)", true);

                    // -------- component 0 --------
                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        Var_arr(iL, jL, kL, 0),
                        "VarBuffer(0)@left_index", "RECON: source left", false, 0, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        Var_arr(iR, jR, kR, 0),
                        "VarBuffer(0)@right_index", "RECON: source right", false, 0, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        LeftStates_arr(iface, jface, kface, 0),
                        "LeftStates(0)", "RECON: wrote LeftStates", false, 0, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        RightStates_arr(iface, jface, kface, 0),
                        "RightStates(0)", "RECON: wrote RightStates", false, 0, true);

                    // -------- component 1 --------
                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        Var_arr(iL, jL, kL, 1),
                        "VarBuffer(1)@left_index", "RECON: source left", false, 1, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        Var_arr(iR, jR, kR, 1),
                        "VarBuffer(1)@right_index", "RECON: source right", false, 1, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        LeftStates_arr(iface, jface, kface, 1),
                        "LeftStates(1)", "RECON: wrote LeftStates", false, 1, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iL, jL, kL,
                        RightStates_arr(iface, jface, kface, 1),
                        "RightStates(1)", "RECON: wrote RightStates", false, 1, true);

                    // -------- component 2 --------
                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        Var_arr(iL, jL, kL, 2),
                        "VarBuffer(2)@left_index", "RECON: source left", false, 2, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        Var_arr(iR, jR, kR, 2),
                        "VarBuffer(2)@right_index", "RECON: source right", false, 2, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        LeftStates_arr(iface, jface, kface, 2),
                        "LeftStates(2)", "RECON: wrote LeftStates", false, 2, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        RightStates_arr(iface, jface, kface, 2),
                        "RightStates(2)", "RECON: wrote RightStates", false, 2, true);

                    // -------- component 3 --------
                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        Var_arr(iL, jL, kL, 3),
                        "VarBuffer(3)@left_index", "RECON: source left", false, 3, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        Var_arr(iR, jR, kR, 3),
                        "VarBuffer(3)@right_index", "RECON: source right", false, 3, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        LeftStates_arr(iface, jface, kface, 3),
                        "LeftStates(3)", "RECON: wrote LeftStates", false, 3, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        RightStates_arr(iface, jface, kface, 3),
                        "RightStates(3)", "RECON: wrote RightStates", false, 3, true);

                    // -------- component 4 --------
                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        Var_arr(iL, jL, kL, 4),
                        "VarBuffer(4)@left_index", "RECON: source left", false, 4, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        Var_arr(iR, jR, kR, 4),
                        "VarBuffer(4)@right_index", "RECON: source right", false, 4, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        LeftStates_arr(iface, jface, kface, 4),
                        "LeftStates(4)", "RECON: wrote LeftStates", false, 4, true);

                    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                        iface, jface, kface,
                        RightStates_arr(iface, jface, kface, 4),
                        "RightStates(4)", "RECON: wrote RightStates", false, 4, true);
                }

                // if (i == 32 && j == 25 && (k == 7 || k == 8)) {
                //     const int rho = solver->variableIndex.DENS;
                //     amrex::Print() << "[RECON WRITE] face=(" << i << "," << j << "," << k << ")"
                //                    << " LS(rho)=" << LeftStates_arr(index[0],index[1],index[2],rho)
                //                    << " RS(rho)=" << RightStates_arr(index[0],index[1],index[2],rho)
                //                    << "\n";
                // }
                
            });
        }
        // LeftStates.FillBoundary();
        // RightStates.FillBoundary();
    }
};

// Local Lax-Friedrichs flux method
template <typename IntegratorDerivedClass>
class LocalLaxFriedrichsMethod : public FluxMethod
<IntegratorDerivedClass> {
public:
    using typename FluxMethod<IntegratorDerivedClass>::ReconstructionTarget;

    // Indicate that LLF works with flux-based reconstruction
    ReconstructionTarget getReconstructionTarget() const override {
        return ReconstructionTarget::Fluxes;
    }

    void SplitFluxes(int direction, int lev, IntegratorDerivedClass* solver, amrex::MultiFab& VariableBuffer, amrex::MultiFab& CellFluxBuffer, amrex::MultiFab& Flux_Negative, amrex::MultiFab& Flux_Positive) const override {

        const Set::Scalar gamma = 1.4;  // Ratio of specific heats

        //for (amrex::MFIter mfi(QL_stencil, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        for (amrex::MFIter mfi(Flux_Negative, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

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
            FluxPos_arr(index[0], index[1], index[2], n) =  0.5 * FL[n] + 0.5 * lambda_max * QL_cons[n];
            FluxNeg_arr(index[0], index[1], index[2], n) =  0.5 * FR[n] - 0.5 * lambda_max * QR_cons[n];
            }

        });


        }

    }

    // Implement ComputeFluxes for completeness (though it won't be used for LLF)
    void ComputeFluxes(int /*direction*/, int /*lev*/, IntegratorDerivedClass* /*solver*/,
                    const amrex::MultiFab& /*LeftStates*/, const amrex::MultiFab& /*RightStates*/,
                    amrex::MultiFab& /*SummedFlux*/) const override {
        // This is not the primary pathway for LLF, but provide an implementation
        // for completeness or potential future use
        Util::Warning(INFO, "LLF ComputeFluxes called directly - not optimal path");

        // If needed, a direct implementation could be added here
    }    

};

// HLLC flux method with symmetry-preserving operations
template <typename IntegratorDerivedClass>
class HLLCMethod : public FluxMethod<IntegratorDerivedClass> {
public:

    using typename FluxMethod<IntegratorDerivedClass>::ReconstructionTarget;

    ReconstructionTarget getReconstructionTarget() const override {
        return ReconstructionTarget::Variables;
    }

    void SplitFluxes(int /*direction*/, int /*lev*/, IntegratorDerivedClass* /*solver*/, 
                    amrex::MultiFab& /*VariableBuffer*/, amrex::MultiFab& /*CellFluxBuffer*/,
                    amrex::MultiFab& /*Flux_Negative*/, amrex::MultiFab& /*Flux_Positive*/) const override {
        // Not the primary pathway for HLLC
        Util::Warning(INFO, "HLLC SplitFluxes called - not optimal path");
    }

void ComputeFluxes(int direction, int lev [[maybe_unused]], IntegratorDerivedClass* solver,
            const amrex::MultiFab& LeftStates, const amrex::MultiFab& RightStates,
            amrex::MultiFab& SummedFlux) const override {

    const Set::Scalar gamma = 1.4;  // Ratio of specific heats
    // const int nghosts = solver->number_of_ghost_cells;

    // Loop through the grid
    for (amrex::MFIter mfi(SummedFlux, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);

        const amrex::Box& face_bx = mfi.tilebox();
        // const amrex::Box& vbx = mfi.validbox();

        // Get arrays for left/right reconstructed states and flux output
        auto const& Left_arr = LeftStates.array(mfi);
        auto const& Right_arr = RightStates.array(mfi);
        auto const& TotalFlux_arr = SummedFlux.array(mfi);

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
        amrex::ParallelFor(face_bx, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {
            // // Setup indices
            // int index[3];
            // int lower_bounds[3];
            // int upper_bounds[3];

            // // Initialize index arrays
            // index[0] = i; index[1] = j; index[2] = k;

            // // Define bounds including ghost cells
            // lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); 
            // lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); 
            // lower_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.smallEnd(2) : 0;
            // upper_bounds[0] = face_bx_with_ghosts.bigEnd(0) - 1;
            // upper_bounds[1] = face_bx_with_ghosts.bigEnd(1) - 1;
            // upper_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.bigEnd(2) - 1 : 0;

            // lower_bounds[0] = vbx.smallEnd(0); 
            // lower_bounds[1] = vbx.smallEnd(1); 
            // lower_bounds[2] = AMREX_SPACEDIM == 3 ? vbx.smallEnd(2): 0;
            // upper_bounds[0] = vbx.bigEnd(0);
            // upper_bounds[1] = vbx.bigEnd(1);
            // upper_bounds[2] = AMREX_SPACEDIM == 3 ? vbx.bigEnd(2) : 0;

            // // Clamp indices to ensure they're within bounds
            // ClampIndices(index, lower_bounds, upper_bounds);

            // Get LEFT conservative variables from reconstructed state
            Set::Scalar rho_L = Left_arr(iface, jface, kface, solver->variableIndex.DENS);
            rho_L = std::max(rho_L, 1e-8);  // Avoid zero density

            Set::Scalar mom_L_normal = Left_arr(iface, jface, kface, velocity_component[normal]);
            Set::Scalar vel_L_normal = mom_L_normal / rho_L;

            Set::Scalar vel_L[3] = {0.0, 0.0, 0.0};
            vel_L[normal] = vel_L_normal;

#if AMREX_SPACEDIM >= 2    
            Set::Scalar mom_L_trans1 = Left_arr(iface, jface, kface, velocity_component[trans1]);
            vel_L[trans1] = mom_L_trans1 / rho_L;
#endif
#if AMREX_SPACEDIM == 3    
            Set::Scalar mom_L_trans2 = Left_arr(iface, jface, kface, velocity_component[trans2]);
            vel_L[trans2] = mom_L_trans2 / rho_L;
#endif

            Set::Scalar rhoE_L = Left_arr(iface, jface, kface, solver->variableIndex.IE);
            Set::Scalar E_L = rhoE_L / rho_L;

            // Calculate kinetic energy for left state using symmetry-preserving operation
            Set::Scalar KE_L = Numeric::SymmetryPreserving::ConsistentKineticEnergy(
                rho_L, vel_L[0], vel_L[1], vel_L[2]);

            // Calculate internal energy and pressure for left state
            Set::Scalar e_L = E_L - KE_L / rho_L;
            Set::Scalar pressure_L = (gamma - 1.0) * rho_L * e_L;
            pressure_L = std::max(pressure_L, 1e-8);  // Ensure positive pressure

            // Get RIGHT conservative variables from reconstructed state
            Set::Scalar rho_R = Right_arr(iface, jface, kface, solver->variableIndex.DENS);
            rho_R = std::max(rho_R, 1e-8);  // Avoid zero density

            Set::Scalar mom_R_normal = Right_arr(iface, jface, kface, velocity_component[normal]);
            Set::Scalar vel_R_normal = mom_R_normal / rho_R;

            Set::Scalar vel_R[3] = {0.0, 0.0, 0.0};
            vel_R[normal] = vel_R_normal;

#if AMREX_SPACEDIM >= 2    
            Set::Scalar mom_R_trans1 = Right_arr(iface, jface, kface, velocity_component[trans1]);
            vel_R[trans1] = mom_R_trans1 / rho_R;
#endif
#if AMREX_SPACEDIM == 3    
            Set::Scalar mom_R_trans2 = Right_arr(iface, jface, kface, velocity_component[trans2]);
            vel_R[trans2] = mom_R_trans2 / rho_R;
#endif

            Set::Scalar rhoE_R = Right_arr(iface, jface, kface, solver->variableIndex.IE);
            Set::Scalar E_R = rhoE_R / rho_R;

            // Calculate kinetic energy for right state using symmetry-preserving operation
            Set::Scalar KE_R = Numeric::SymmetryPreserving::ConsistentKineticEnergy(
                rho_R, vel_R[0], vel_R[1], vel_R[2]);

            // Calculate internal energy and pressure for right state
            Set::Scalar e_R = E_R - KE_R / rho_R;
            Set::Scalar pressure_R = (gamma - 1.0) * rho_R * e_R;
            pressure_R = std::max(pressure_R, 1e-8);  // Ensure positive pressure

            // Sound speeds
            Set::Scalar c_L = std::sqrt(gamma * pressure_L / rho_L);
            Set::Scalar c_R = std::sqrt(gamma * pressure_R / rho_R);

            // Store conservative variables for reconstructed left and right states
            std::vector<Set::Scalar> QL_cons(num_components, 0.0);
            std::vector<Set::Scalar> QR_cons(num_components, 0.0);

            // These are already in conservative form from the reconstruction
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

            // Momentum fluxes - use consistent operations
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

            // Energy flux - use consistent operations
            FL[solver->variableIndex.IE] = (rho_L * E_L + pressure_L) * vel_L[normal];
            FR[solver->variableIndex.IE] = (rho_R * E_R + pressure_R) * vel_R[normal];

            // Wave speeds using Davis estimates (more robust than Roe averages)
            Set::Scalar SL = std::min(vel_L[normal] - c_L, vel_R[normal] - c_R);
            Set::Scalar SR = std::max(vel_L[normal] + c_L, vel_R[normal] + c_R);

            // Apply entropy fix for transonic rarefactions
            Set::Scalar epsilon = 0.1 * std::min(c_L, c_R); 
            if (std::abs(SR - SL) < epsilon) {
                SR += epsilon;
                SL -= epsilon;    
            }   

            // Calculate the middle wave speed (S*) using the symmetry-preserving operation
            Set::Scalar S_star = Numeric::SymmetryPreserving::ConsistentContactWaveSpeed(
                pressure_R, pressure_L, 
                rho_L, vel_L[normal], SL, 
                rho_R, vel_R[normal], SR);

            // Compute star state variables
            std::vector<Set::Scalar> QL_star(num_components, 0.0);
            std::vector<Set::Scalar> QR_star(num_components, 0.0);

            // Density factors for star regions
            Set::Scalar rho_L_star = rho_L * (SL - vel_L[normal]) / (SL - S_star);
            Set::Scalar rho_R_star = rho_R * (SR - vel_R[normal]) / (SR - S_star);

            // Left star state
            QL_star[solver->variableIndex.DENS] = rho_L_star;
            QL_star[solver->variableIndex.UVEL + normal] = rho_L_star * S_star;
#if AMREX_SPACEDIM >= 2
            QL_star[solver->variableIndex.UVEL + trans1] = rho_L_star * vel_L[trans1];
#endif
#if AMREX_SPACEDIM == 3
            QL_star[solver->variableIndex.UVEL + trans2] = rho_L_star * vel_L[trans2];
#endif
            // Use consistent operations for the star state energy
            Set::Scalar pressureTerm_L = pressure_L / (rho_L * (SL - vel_L[normal]));
            Set::Scalar velocityDiff_L = S_star - vel_L[normal];

            QL_star[solver->variableIndex.IE] = rho_L_star * (E_L + velocityDiff_L * (S_star + pressureTerm_L));

            // Right star state
            QR_star[solver->variableIndex.DENS] = rho_R_star;
            QR_star[solver->variableIndex.UVEL + normal] = rho_R_star * S_star;
#if AMREX_SPACEDIM >= 2
            QR_star[solver->variableIndex.UVEL + trans1] = rho_R_star * vel_R[trans1];
#endif
#if AMREX_SPACEDIM == 3
            QR_star[solver->variableIndex.UVEL + trans2] = rho_R_star * vel_R[trans2];
#endif
            // Use consistent operations for the star state energy
            Set::Scalar pressureTerm_R = pressure_R / (rho_R * (SR - vel_R[normal]));
            Set::Scalar velocityDiff_R = S_star - vel_R[normal];

            QR_star[solver->variableIndex.IE] = rho_R_star * (E_R + velocityDiff_R * (S_star + pressureTerm_R));

            // Use a centralized formulation for computing consistent fluxes
            // Precompute star region fluxes
            std::vector<Set::Scalar> FL_star(num_components, 0.0);
            std::vector<Set::Scalar> FR_star(num_components, 0.0);

            for (int n = 0; n < num_components; ++n) {
                FL_star[n] = FL[n] + SL * (QL_star[n] - QL_cons[n]);
                FR_star[n] = FR[n] + SR * (QR_star[n] - QR_cons[n]);
            }

            // Define Heaviside-like functions for the wave speeds - outside the loop since they're the same for all components
            Set::Scalar HL = 0.5 * (1.0 + std::copysign(1.0, SL));     // H(SL): 1 if SL >= 0, 0 if SL < 0
            Set::Scalar HStar = 0.5 * (1.0 + std::copysign(1.0, S_star)); // H(S*): 1 if S* >= 0, 0 if S* < 0
            Set::Scalar HR = 0.5 * (1.0 + std::copysign(1.0, SR));     // H(SR): 1 if SR >= 0, 0 if SR < 0

            // Branch-free HLLC flux calculation
            for (int n = 0; n < num_components; ++n) {
                // Compute full HLLC flux using branch-free formulation
                Set::Scalar HLLC_flux = 
                    HL * FL[n] +                                  // Case: SL >= 0
                    (1.0 - HL) * HStar * FL_star[n] +            // Case: SL < 0, S* >= 0
                    (1.0 - HStar) * HR * FR_star[n] +            // Case: S* < 0, SR >= 0
                    (1.0 - HR) * FR[n];                          // Case: SR < 0

                // Store the computed flux
                TotalFlux_arr(iface, jface, kface, n) = HLLC_flux;

                
                if (Util::ScimitarX_Util::Debug::IsTargetLocation(iface, jface, kface))
{
    // -------- component 0 --------

    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
        iface, jface, kface,
        Left_arr(iface, jface, kface, 0),
        "LeftStates(0)", "FLUX: input left", false, 0, true);

    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
        iface, jface, kface,
        Right_arr(iface, jface, kface, 0),
        "RightStates(0)", "FLUX: input right", false, 0, true);

    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
        iface, jface, kface,
        TotalFlux_arr(iface, jface, kface, 0),
        "SummedFlux(0)", "FLUX: output", false, 0, true);

    // -------- component 1 --------
    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
        iface, jface, kface,
        Left_arr(iface, jface, kface, 1),
        "LeftStates(1)", "FLUX: input left", false, 1, true);

    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
        iface, jface, kface,
        Right_arr(iface, jface, kface, 1),
        "RightStates(1)", "FLUX: input right", false, 1, true);

    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
        iface, jface, kface,
        TotalFlux_arr(iface, jface, kface, 1),
        "SummedFlux(1)", "FLUX: output", false, 1, true);

    // -------- component 2 --------
    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
        iface, jface, kface,
        Left_arr(iface, jface, kface, 2),
        "LeftStates(2)", "FLUX: input left", false, 2, true);

    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
        iface, jface, kface,
        Right_arr(iface, jface, kface, 2),
        "RightStates(2)", "FLUX: input right", false, 2, true);

    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
        iface, jface, kface,
        TotalFlux_arr(iface, jface, kface, 2),
        "SummedFlux(2)", "FLUX: output", false, 2, true);

    // -------- component 3 --------
    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
        iface, jface, kface,
        Left_arr(iface, jface, kface, 3),
        "LeftStates(3)", "FLUX: input left", false, 3, true);

    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
        iface, jface, kface,
        Right_arr(iface, jface, kface, 3),
        "RightStates(3)", "FLUX: input right", false, 3, true);

    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
        iface, jface, kface,
        TotalFlux_arr(iface, jface, kface, 3),
        "SummedFlux(3)", "FLUX: output", false, 3, true);

    // -------- component 4 --------
    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
        iface, jface, kface,
        Left_arr(iface, jface, kface, 4),
        "LeftStates(4)", "FLUX: input left", false, 4, true);

    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
        iface, jface, kface,
        Right_arr(iface, jface, kface, 4),
        "RightStates(4)", "FLUX: input right", false, 4, true);

    Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
        iface, jface, kface,
        TotalFlux_arr(iface, jface, kface, 4),
        "SummedFlux(4)", "FLUX: output", false, 4, true);
}



            }

        });
    }
    // SummedFlux.FillBoundary();
}

};


/**
* @class HLLCLowMachMethod
* @brief HLLC flux with reduced dissipation for low Mach regions
* 
* This is an implementation of the HLLC-LM method as described in:
* - Fleischmann, N., Adami, S., & Adams, N. A. (2020). A shock-stable modification 
*   of the HLLC Riemann solver with reduced numerical dissipation. Journal of 
*   Computational Physics, Volume 423, 109762.
* 
* The method reduces acoustic dissipation in low Mach number regions to prevent
* grid-aligned shock instabilities while maintaining accuracy.
*/

// HLLCLowMach flux method with reduced dissipation for improved shock stability
/*template <typename IntegratorDerivedClass>
class HLLCLowMachMethod : public FluxMethod<IntegratorDerivedClass> {
public:
    using typename FluxMethod<IntegratorDerivedClass>::ReconstructionTarget;

    // Use variable-based reconstruction like the base HLLC
    ReconstructionTarget getReconstructionTarget() const override {
        return ReconstructionTarget::Variables;
    }

    void SplitFluxes(int direction, int lev, IntegratorDerivedClass* solver, 
                    amrex::MultiFab& VariableBuffer, amrex::MultiFab& CellFluxBuffer,
                    amrex::MultiFab& Flux_Negative, amrex::MultiFab& Flux_Positive) const override {
        // Not the primary pathway for HLLCLowMach
        Util::Warning(INFO, "HLLCLowMach SplitFluxes called - not optimal path");
    }

    void ComputeFluxes(int direction, int lev, IntegratorDerivedClass* solver,
                    const amrex::MultiFab& LeftStates, const amrex::MultiFab& RightStates,
                    amrex::MultiFab& SummedFlux) const override {

        const Set::Scalar gamma = 1.4;  // Ratio of specific heats
        const int nghosts = solver->number_of_ghost_cells;
        const Set::Scalar Ma_limit = 0.1; // Mach number limit for applying dissipation reduction

        // Loop through the grid
        for (amrex::MFIter mfi(SummedFlux, false); mfi.isValid(); ++mfi) {
            const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);

            // Get arrays for left/right reconstructed states and flux output
            auto const& Left_arr = LeftStates.array(mfi);
            auto const& Right_arr = RightStates.array(mfi);
            auto const& TotalFlux_arr = SummedFlux.array(mfi);

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

                // Get LEFT conservative variables from reconstructed state
                Set::Scalar rho_L = Left_arr(index[0], index[1], index[2], solver->variableIndex.DENS);
                rho_L = std::max(rho_L, 1e-8);  // Avoid zero density

                Set::Scalar mom_L_normal = Left_arr(index[0], index[1], index[2], velocity_component[normal]);
                Set::Scalar vel_L_normal = mom_L_normal / rho_L;

                Set::Scalar vel_L[3] = {0.0, 0.0, 0.0};
                vel_L[normal] = vel_L_normal;

#if AMREX_SPACEDIM >= 2    
                Set::Scalar mom_L_trans1 = Left_arr(index[0], index[1], index[2], velocity_component[trans1]);
                vel_L[trans1] = mom_L_trans1 / rho_L;
#endif
#if AMREX_SPACEDIM == 3    
                Set::Scalar mom_L_trans2 = Left_arr(index[0], index[1], index[2], velocity_component[trans2]);
                vel_L[trans2] = mom_L_trans2 / rho_L;
#endif

                Set::Scalar rhoE_L = Left_arr(index[0], index[1], index[2], solver->variableIndex.IE);
                Set::Scalar E_L = rhoE_L / rho_L;

                // Calculate kinetic energy for left state using symmetry-preserving operation
                Set::Scalar KE_L = Numeric::SymmetryPreserving::ConsistentKineticEnergy(
                    rho_L, vel_L[0], vel_L[1], vel_L[2]);

                // Calculate internal energy and pressure for left state
                Set::Scalar e_L = E_L - KE_L / rho_L;
                Set::Scalar pressure_L = (gamma - 1.0) * rho_L * e_L;
                pressure_L = std::max(pressure_L, 1e-8);  // Ensure positive pressure

                // Get RIGHT conservative variables from reconstructed state
                Set::Scalar rho_R = Right_arr(index[0], index[1], index[2], solver->variableIndex.DENS);
                rho_R = std::max(rho_R, 1e-8);  // Avoid zero density

                Set::Scalar mom_R_normal = Right_arr(index[0], index[1], index[2], velocity_component[normal]);
                Set::Scalar vel_R_normal = mom_R_normal / rho_R;

                Set::Scalar vel_R[3] = {0.0, 0.0, 0.0};
                vel_R[normal] = vel_R_normal;

#if AMREX_SPACEDIM >= 2    
                Set::Scalar mom_R_trans1 = Right_arr(index[0], index[1], index[2], velocity_component[trans1]);
                vel_R[trans1] = mom_R_trans1 / rho_R;
#endif
#if AMREX_SPACEDIM == 3    
                Set::Scalar mom_R_trans2 = Right_arr(index[0], index[1], index[2], velocity_component[trans2]);
                vel_R[trans2] = mom_R_trans2 / rho_R;
#endif

                Set::Scalar rhoE_R = Right_arr(index[0], index[1], index[2], solver->variableIndex.IE);
                Set::Scalar E_R = rhoE_R / rho_R;

                // Calculate kinetic energy for right state using symmetry-preserving operation
                Set::Scalar KE_R = Numeric::SymmetryPreserving::ConsistentKineticEnergy(
                    rho_R, vel_R[0], vel_R[1], vel_R[2]);

                // Calculate internal energy and pressure for right state
                Set::Scalar e_R = E_R - KE_R / rho_R;
                Set::Scalar pressure_R = (gamma - 1.0) * rho_R * e_R;
                pressure_R = std::max(pressure_R, 1e-8);  // Ensure positive pressure

                // Sound speeds
                Set::Scalar c_L = std::sqrt(gamma * pressure_L / rho_L);
                Set::Scalar c_R = std::sqrt(gamma * pressure_R / rho_R);

                // Store conservative variables for reconstructed left and right states
                std::vector<Set::Scalar> QL_cons(num_components, 0.0);
                std::vector<Set::Scalar> QR_cons(num_components, 0.0);

                // These are already in conservative form from the reconstruction
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

                // Non-modified signal speeds (original HLLC)
                Set::Scalar SL_orig = std::min(vel_L[normal] - c_L, vel_R[normal] - c_R);
                Set::Scalar SR_orig = std::max(vel_L[normal] + c_L, vel_R[normal] + c_R);

                // Local Mach number calculation - use the maximum of left and right
                Set::Scalar Ma_local = std::max(std::abs(vel_L[normal])/c_L, std::abs(vel_R[normal])/c_R);

                // Calculate phi reduction factor based on local Mach number
                // Using the sine function as described in the paper
                Set::Scalar phi = std::sin(std::min(1.0, Ma_local/Ma_limit) * (3.14159265358979323846/2.0));

                // Apply the low Mach number modification to signal speeds as per equation (23)
                Set::Scalar SL = phi * SL_orig;
                Set::Scalar SR = phi * SR_orig;

                // Apply entropy fix for transonic rarefactions
                Set::Scalar epsilon = 0.1 * std::min(c_L, c_R); 
                if (std::abs(SR - SL) < epsilon) {
                    SR += epsilon;
                    SL -= epsilon;    
                }   

                // Calculate the middle wave speed (S*) using symmetry-preserving operation
                Set::Scalar S_star = Numeric::SymmetryPreserving::ConsistentContactWaveSpeed(
                    pressure_R, pressure_L, rho_L, vel_L[normal], SL_orig, rho_R, vel_R[normal], SR_orig);

                // Compute star state variables
                std::vector<Set::Scalar> QL_star(num_components, 0.0);
                std::vector<Set::Scalar> QR_star(num_components, 0.0);

                // Density factors for star regions - use original (non-modified) signal speeds for star state computation
                Set::Scalar rho_L_star = rho_L * (SL_orig - vel_L[normal]) / (SL_orig - S_star);
                Set::Scalar rho_R_star = rho_R * (SR_orig - vel_R[normal]) / (SR_orig - S_star);

                // Left star state
                QL_star[solver->variableIndex.DENS] = rho_L_star;
                QL_star[solver->variableIndex.UVEL + normal] = rho_L_star * S_star;
#if AMREX_SPACEDIM >= 2
                QL_star[solver->variableIndex.UVEL + trans1] = rho_L_star * vel_L[trans1];
#endif
#if AMREX_SPACEDIM == 3
                QL_star[solver->variableIndex.UVEL + trans2] = rho_L_star * vel_L[trans2];
#endif
                QL_star[solver->variableIndex.IE] = rho_L_star * (E_L + (S_star - vel_L[normal]) * 
                                                (S_star + pressure_L / (rho_L * (SL_orig - vel_L[normal]))));

                // Right star state
                QR_star[solver->variableIndex.DENS] = rho_R_star;
                QR_star[solver->variableIndex.UVEL + normal] = rho_R_star * S_star;
#if AMREX_SPACEDIM >= 2
                QR_star[solver->variableIndex.UVEL + trans1] = rho_R_star * vel_R[trans1];
#endif
#if AMREX_SPACEDIM == 3
                QR_star[solver->variableIndex.UVEL + trans2] = rho_R_star * vel_R[trans2];
#endif
                QR_star[solver->variableIndex.IE] = rho_R_star * (E_R + (S_star - vel_R[normal]) * 
                                                (S_star + pressure_R / (rho_R * (SR_orig - vel_R[normal]))));

                // Precompute star region fluxes (using centralized formulation from equation 7)
                std::vector<Set::Scalar> FL_star(num_components, 0.0);
                std::vector<Set::Scalar> FR_star(num_components, 0.0);

                for (int n = 0; n < num_components; ++n) {
                    FL_star[n] = FL[n] + SL * (QL_star[n] - QL_cons[n]);
                    FR_star[n] = FR[n] + SR * (QR_star[n] - QR_cons[n]);
                }

                // Centralized HLLC-LM flux computation as per equation (19)
                std::vector<Set::Scalar> HLLC_flux(num_components, 0.0);
                // Branch-free centralized HLLC-LM flux computation
                for (int n = 0; n < num_components; ++n) {
                    // Create Heaviside-like functions for SL and SR using sign bits
                    Set::Scalar HSL = 0.5 * (1.0 + std::copysign(1.0, SL));  // 1 if SL >= a0, 0 otherwise
                    Set::Scalar HSR = 0.5 * (1.0 + std::copysign(1.0, SR));  // 1 if SR >= 0, 0 otherwise
                    Set::Scalar HStar = 0.5 * (1.0 + std::copysign(1.0, S_star));  // 1 if S* >= 0, 0 otherwise

                    // Calculate star region dissipation contributions
                    Set::Scalar term1 = SL * (QL_star[n] - QL_cons[n]);
                    Set::Scalar term2 = std::abs(S_star) * (2.0 * HStar - 1.0) * (QL_star[n] - QR_star[n]);
                    Set::Scalar term3 = SR * (QR_star[n] - QR_cons[n]);

                    // Use consistent sum for the three dissipation terms
                    Set::Scalar subsonic_dissipation = Numeric::SymmetryPreserving::ConsistentSum3(term1, term2, term3);

                    // Combine fluxes without branching:
                    // If SL >= 0: FL[n]
                    // If SR <= 0: FR[n]
                    // Otherwise: 0.5*(FL[n]+FR[n]) + 0.5*dissipation
                    HLLC_flux[n] = HSL * FL[n] + 
                                    (1.0 - HSR) * FR[n] + 
                                    (1.0 - HSL) * HSR * (0.5 * (FL[n] + FR[n]) + 0.5 * subsonic_dissipation);

                    TotalFlux_arr(index[0], index[1], index[2], n) = HLLC_flux[n];
                }

            });
        }
    }
};*/

// AUSM+up flux method with symmetry-preserving operations

template <typename IntegratorDerivedClass>
class AUSMupMethod : public FluxMethod<IntegratorDerivedClass> {
public:
    
    using typename FluxMethod<IntegratorDerivedClass>::ReconstructionTarget;

    ReconstructionTarget getReconstructionTarget() const override {
        return ReconstructionTarget::Variables;
    }
   
    void SplitFluxes(int /*direction*/, int /*lev*/, IntegratorDerivedClass* /*solver*/, 
                    amrex::MultiFab& /*VariableBuffer*/, amrex::MultiFab& /*CellFluxBuffer*/,
                    amrex::MultiFab& /*Flux_Negative*/, amrex::MultiFab& /*Flux_Positive*/) const override {
        // Not the primary pathway for AUSM+up
        Util::Warning(INFO, "AUSM+up SplitFluxes called - not optimal path");
    }

Set::Scalar M_infinity;
Set::Scalar Mp_cut;
Set::Scalar Kp;
Set::Scalar Ku;
Set::Scalar sigma;
 
static void Parse(AUSMupMethod& value, IO::ParmParse& pp)
{
pp.query_default("freeStreamMachNumber", value.M_infinity, 1.0);
pp.query_default("Mp_cut", value.Mp_cut, 0.1);
pp.query_default("Kp", value.Kp, 0.25);
pp.query_default("Ku", value.Ku, 0.75);
pp.query_default("sigma", value.sigma, 1.0);
}
  
 
void ComputeFluxes(int direction, int /*lev*/, IntegratorDerivedClass* solver,
            const amrex::MultiFab& LeftStates, const amrex::MultiFab& RightStates,
            amrex::MultiFab& SummedFlux) const override {

    const Set::Scalar gamma = 1.4;  // Ratio of specific heats
    const int nghosts = solver->number_of_ghost_cells;

    IO::ParmParse pp;  
    AUSMupMethod::Parse(*const_cast<AUSMupMethod*>(this), pp);
    const Set::Scalar M_inf = this->M_infinity; //SHOULD NOT BE ZERO

//    Util::ScimitarX_Util::Debug::SetTargetDebugLocationIndices(50, 10, 0, "true");

    // Loop through the grid
    for (amrex::MFIter mfi(SummedFlux, false); mfi.isValid(); ++mfi) {
        const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);

        // Get arrays for left/right reconstructed states and flux output
        auto const& Left_arr = LeftStates.array(mfi);
        auto const& Right_arr = RightStates.array(mfi);
        auto const& TotalFlux_arr = SummedFlux.array(mfi);
        
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


            // Get LEFT conservative variables from reconstructed state
            Set::Scalar rho_L = Left_arr(index[0], index[1], index[2], solver->variableIndex.DENS);
            rho_L = std::max(rho_L, 1.0e-8);  // Avoid zero density
            
            Set::Scalar mom_L_normal = Left_arr(index[0], index[1], index[2], velocity_component[normal]);
            Set::Scalar vel_L_normal = mom_L_normal / rho_L;
            
            Set::Scalar vel_L[3] = {0.0, 0.0, 0.0};
            vel_L[normal] = vel_L_normal;
            
#if AMREX_SPACEDIM >= 2    
            Set::Scalar mom_L_trans1 = Left_arr(index[0], index[1], index[2], velocity_component[trans1]);
            vel_L[trans1] = mom_L_trans1 / rho_L;
#endif
#if AMREX_SPACEDIM == 3    
            Set::Scalar mom_L_trans2 = Left_arr(index[0], index[1], index[2], velocity_component[trans2]);
            vel_L[trans2] = mom_L_trans2 / rho_L;
#endif
            
            Set::Scalar rhoE_L = Left_arr(index[0], index[1], index[2], solver->variableIndex.IE);
            Set::Scalar E_L = std::max(rhoE_L / rho_L, 1.0e-8);
            
            // Calculate kinetic energy for left state using symmetry-preserving operation
            Set::Scalar KE_L = Numeric::SymmetryPreserving::ConsistentKineticEnergy(
                rho_L, vel_L[0], vel_L[1], vel_L[2]);
            
            // Calculate internal energy and pressure for left state
            Set::Scalar e_L = std::max(E_L - KE_L / rho_L, 1.0e-8);
            Set::Scalar pressure_L = (gamma - 1.0) * rho_L * e_L;
            pressure_L = std::max(pressure_L, 1.0e-8);  // Ensure positive pressure
            
            // Get RIGHT conservative variables from reconstructed state
            Set::Scalar rho_R = Right_arr(index[0], index[1], index[2], solver->variableIndex.DENS);
            rho_R = std::max(rho_R, 1.0e-8);  // Avoid zero density
            
            Set::Scalar mom_R_normal = Right_arr(index[0], index[1], index[2], velocity_component[normal]);
            Set::Scalar vel_R_normal = mom_R_normal / rho_R;
            
            Set::Scalar vel_R[3] = {0.0, 0.0, 0.0};
            vel_R[normal] = vel_R_normal;
            
#if AMREX_SPACEDIM >= 2    
            Set::Scalar mom_R_trans1 = Right_arr(index[0], index[1], index[2], velocity_component[trans1]);
            vel_R[trans1] = mom_R_trans1 / rho_R;
#endif
#if AMREX_SPACEDIM == 3    
            Set::Scalar mom_R_trans2 = Right_arr(index[0], index[1], index[2], velocity_component[trans2]);
            vel_R[trans2] = mom_R_trans2 / rho_R;
#endif
            
            Set::Scalar rhoE_R = Right_arr(index[0], index[1], index[2], solver->variableIndex.IE);
            Set::Scalar E_R = std::max(rhoE_R / rho_R, 1.0e-8);
            
            // Calculate kinetic energy for right state using symmetry-preserving operation
            Set::Scalar KE_R = Numeric::SymmetryPreserving::ConsistentKineticEnergy(
                rho_R, vel_R[0], vel_R[1], vel_R[2]);
            
            // Calculate internal energy and pressure for right state
            Set::Scalar e_R = std::max(E_R - KE_R / rho_R, 1.0e-8);
            Set::Scalar pressure_R = (gamma - 1.0) * rho_R * e_R;
            pressure_R = std::max(pressure_R, 1.0e-8);  // Ensure positive pressure
            
            // Sound speeds
            Set::Scalar c_L = std::sqrt(std::max(1.0e-8, gamma * pressure_L / rho_L));
            Set::Scalar c_R = std::sqrt(std::max(1.0e-8, gamma * pressure_R / rho_R));

        // ######################################
            // Define specific total Enthalpy H_total 
        Set::Scalar H_total_L = std::max(1.0e-08,E_L + pressure_L/rho_L);  // = Cv.T + (Vel^2)/2 + P/rho = specific Total Enthalpy
        Set::Scalar H_total_R = std::max(1.0e-08, E_R + pressure_R/rho_R);

        // Define critical speed of sound c_star
        Set::Scalar c_star_L = std::sqrt(2*(gamma-1)*H_total_L/(gamma+1));
        Set::Scalar c_star_R = std::sqrt(2*(gamma-1)*H_total_R/(gamma+1));
                    
        // Define c_L_hat and c_R_hat
        Set::Scalar c_L_hat = std::pow(c_star_L,2)/nonZeroValue(std::max(c_star_L,std::abs(vel_L[normal])));
        Set::Scalar c_R_hat = std::pow(c_star_R,2)/nonZeroValue(std::max(c_star_R,std::abs(vel_R[normal])));

        // Define interface acoustic speed c_interface
        Set::Scalar c_interface = std::min((0.5*(c_L + c_R)), std::min(c_L_hat,c_R_hat)); 
        // using total enthalpy and critical acoustic speed and simple mean

        // Define M_L and M_R
        Set::Scalar M_L = vel_L[normal]/nonZeroValue(c_interface);
        Set::Scalar M_R = vel_R[normal]/nonZeroValue(c_interface);

        // Define M_localMean and M0
        Set::Scalar M_localMean = std::sqrt(0.5*(M_L*M_L + M_R*M_R));
        Set::Scalar M0 = std::sqrt(std::min(1.0,std::max(M_localMean*M_localMean, M_inf*M_inf)));
        Set::Scalar M0_p = std::min(1.0,std::max(M_localMean, Mp_cut));

        // Define scaling function fa, interface density rho_interface, sigma and Kp constant terms
        //Set::Scalar fa = std::sqrt(std::pow((1.0-M0*M0)*M_localMean , 2) + 4*M0*M0)/(1+M0*M0);
        Set::Scalar fa = M0*(2.0-M0); // scaling factor for velocity diffusion term in pressure flux
        Set::Scalar fa_p = M0_p*(2.0-M0_p); // scaling factor for pressure diffusion term in interface Mach number

        Set::Scalar rho_interface = 0.5*(rho_L + rho_R);                

        Set::Scalar alpha = 3.0/16.0*(-4.0+5.0*fa*fa);
        Set::Scalar beta = 1.0/8.0;

        Set::Scalar M1plus_L = 0.5*(M_L + std::abs(M_L));
        Set::Scalar M1minus_R = 0.5*(M_R - std::abs(M_R));

        Set::Scalar M2plus_L = 0.25*(M_L + 1)*(M_L + 1);
        Set::Scalar M2minus_L = -0.25*(M_L - 1)*(M_L - 1);
        Set::Scalar M2plus_R = 0.25*(M_R + 1)*(M_R + 1);
        Set::Scalar M2minus_R = -0.25*(M_R - 1)*(M_R - 1);

        Set::Scalar A_M_L = 0.5*(1.0 + std::copysign(1.0, std::abs(M_L)-1.0));
        Set::Scalar B_M_L = std::abs(1.0 - A_M_L);
        Set::Scalar A_M_R = 0.5*(1.0 + std::copysign(1.0, std::abs(M_R)-1.0));
        Set::Scalar B_M_R = std::abs(1.0 - A_M_R);

        Set::Scalar M4plus_L = A_M_L*(M1plus_L) + B_M_L*(M2plus_L*(1.0 - 16.0*beta*M2minus_L));
        Set::Scalar M4minus_R = A_M_R*(M1minus_R) + B_M_R*(M2minus_R*(1.0 + 16.0*beta*M2plus_R));

        Set::Scalar P5plus_L = A_M_L*M1plus_L/(M_L + 1.0e-10) + B_M_L*M2plus_L*(2.0 - M_L - 16.0*alpha*M_L*M2minus_L);
        Set::Scalar P5minus_R = A_M_R*M1minus_R/(M_R + 1.0e-10) + B_M_R*M2minus_R*(-2.0 - M_R + 16.0*alpha*M_R*M2plus_R);
        
        // Calculate interface Mach number, M_interface
        Set::Scalar M_pressureDiffusion = (-Kp/(fa_p + 1.0e-10))*std::max((1.0-sigma*M_localMean*M_localMean),0.0)*
                                          (pressure_R-pressure_L)/(rho_interface*c_interface*c_interface + 1.0e-10);

        Set::Scalar M_interface = M4plus_L + M4minus_R + M_pressureDiffusion;

        // Calculate the Mass flux...
        Set::Scalar A = 0.5*(1.0+std::copysign(1.0,(M_interface - 1.0e-10)));
        Set::Scalar B = std::abs(1.0-A);  // when A=1: B=0 and vice versa/ then can only be 0 or 1
        Set::Scalar massFlux = c_interface*M_interface*(A * rho_L + B * rho_R);

        // Pressure flux calculations...        

        Set::Scalar P_velocityDiffusion = -Ku*(P5plus_L)*(P5minus_R)*(2*rho_interface)*(fa*c_interface)*(vel_R[normal]-vel_L[normal]);

        Set::Scalar pressureFlux = (P5plus_L)*pressure_L + (P5minus_R)*pressure_R + P_velocityDiffusion;

        // Debugging for variables and Fluxes
        /*Util::ScimitarX_Util::Debug::DebugAUSMupFlux(i, j, k,
                                                rho_L, vel_L[0], vel_L[1], vel_L[2], E_L, pressure_L,
                                                rho_R, vel_R[0], vel_R[1], vel_R[2], E_R, pressure_R,
                                                c_L, c_R, c_interface, fa, M_L, M_R, M_inf, M_pressureDiffusion, P_velocityDiffusion,
                                                M4plus_L,M4minus_R,P5plus_L,P5minus_R,
                                                M_interface, massFlux, pressureFlux,
                                                "true"); */
        
            // Assigning Fluxes
        std::vector<Set::Scalar> AUSMupFluxes(num_components, 0.0);

        AUSMupFluxes[solver->variableIndex.DENS] = massFlux;
        AUSMupFluxes[solver->variableIndex.UVEL + normal] = massFlux * (A*vel_L[normal] + B*vel_R[normal]) +  pressureFlux;
#if AMREX_SPACEDIM >= 2
        AUSMupFluxes[solver->variableIndex.UVEL + trans1] = massFlux * (A*vel_L[trans1] + B*vel_R[trans1]);
#endif
#if AMREX_SPACEDIM == 3
        AUSMupFluxes[solver->variableIndex.UVEL + trans2] = massFlux * (A*vel_L[trans2] + B*vel_R[trans2]);
#endif
        AUSMupFluxes[solver->variableIndex.IE] = massFlux * (A*E_L + B*E_R) + pressureFlux*(c_interface*M_interface);
         
/*
        Set::Scalar vel_interface = c_interface*M_interface;
        Set::Scalar magVel_interface = std::abs(c_interface*M_interface);

        AUSMupFluxes[solver->variableIndex.DENS] = 0.5*(vel_interface*(rho_L + rho_R) - magVel_interface*(rho_L - rho_R));

        AUSMupFluxes[solver->variableIndex.UVEL + normal] = 0.5*(vel_interface*(mom_L_normal  + mom_R_normal) - 
                                                            magVel_interface*(mom_L_normal - mom_R_normal)) +  
                                                            pressureFlux;
#if AMREX_SPACEDIM >= 2
        AUSMupFluxes[solver->variableIndex.UVEL + trans1] = 0.5*(vel_interface*(mom_L_trans1  + mom_R_trans1) -
                                                            magVel_interface*(mom_L_trans1 - mom_R_trans1));
#endif
#if AMREX_SPACEDIM == 3
        AUSMupFluxes[solver->variableIndex.UVEL + trans2] = 0.5*(vel_interface*(mom_L_trans2  + mom_R_trans2) -
                                                            magVel_interface*(mom_L_trans2 - mom_R_trans2));
#endif
        AUSMupFluxes[solver->variableIndex.IE] = 0.5*(vel_interface*(rhoE_L + rhoE_R) - magVel_interface*(rhoE_L - rhoE_R)) + pressureFlux*vel_interface;
*/

        for (int n = 0; n < num_components; ++n) {
                TotalFlux_arr(index[0], index[1], index[2], n) = AUSMupFluxes[n];
        }
        // ################################################### 
        });
    }
}

private:

// polynomial functions for Mass and Pressure flux calculation:
#include <iostream>

// for mass flux:
Set::Scalar massFluxPolynomialFunction(Set::Scalar M, Set::Scalar sign) const {

    // Define switch variables A and B
    Set::Scalar A = 0.5 * (1.0 + std::copysign(1.0, std::abs(M) - 1.0));
    Set::Scalar B = std::abs(A-1.0);
    // Hence if A is 1, B will be 0, and vice versa
    
    const Set::Scalar beta = 0.125;
    Set::Scalar result = A*(0.5*(M+sign*std::abs(M))) + sign*B*(0.25*std::pow((M+sign),2) + beta*std::pow((M*M-1.0),2));
    return result;
}

// for pressure flux:
Set::Scalar pressureFluxPolynomialFunction(Set::Scalar M, Set::Scalar sign, Set::Scalar fa) const {

    // Define switch variables A and B
    Set::Scalar A = 0.5 * (1.0 + std::copysign(1.0, std::abs(M)-1.0));
    Set::Scalar B = std::abs(A - 1.0);
    // Hence if A is 1, B will be 0, and vice versa

    const Set::Scalar alpha = 3.0*(5.0*fa*fa - 4.0)/16.0;
    Set::Scalar result = A*(0.5*(1.0 + sign * std::copysign(1.0, M))) + B*(0.5*std::pow((M + sign),2) - sign*0.25*M*std::pow((M + sign),2)  + sign*alpha*M*std::pow((M*M - 1),2));
    return result;
}

// to ensure non-zero value
Set::Scalar nonZeroValue(Set::Scalar value) const {
        Set::Scalar result = value + 1.0e-10;
        return result;
}

};


}  // namespace Numeric

#endif  // NUMERIC_FLUXHANDLER_H

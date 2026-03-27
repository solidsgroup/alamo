#include "Numeric/IntegratorVariableAccessLayer.H"
#include "Numeric/SymmetryPreservingRoeAveragingOperations.H"
#include "Integrator/ScimitarX.H"
#include "Numeric/Stencil.H"
#include "Model/Fluid/ElastoPlastic.H"


namespace Numeric {

// Define the static member outside of the class definition
GenericVariableAccessor::VariableIndices GenericVariableAccessor::variableIndices;

namespace CompressibleEuler {    

CompressibleEulerVariableAccessor::CompressibleEulerVariableAccessor(
    int total_ghosts,
    ReconstructionMode mode
) : current_mode(mode), total_ghost_cells(total_ghosts) {}

std::shared_ptr<SolverCapabilities> 
CompressibleEulerVariableAccessor::getSolverCapabilities() const {
    // Create a detailed solver capabilities object for Compressible Euler solver
    return std::make_shared<CompressibleEuler::CompressibleEulerCapabilities>();
}

amrex::MultiFab 
CompressibleEulerVariableAccessor::CreateWorkingBuffer(
    const amrex::BoxArray& baseGrids, 
    const amrex::DistributionMapping& dm, 
    const int num_components, 
    const int ghost_cells
) const {

    return amrex::MultiFab(
        baseGrids, 
        dm, 
        num_components, 
        ghost_cells
    );
}

// Implementations of other methods remain similar to previous version
// with added validation and more flexible handling

int 
CompressibleEulerVariableAccessor::getRequiredGhostCells(
    ReconstructionMode mode [[maybe_unused]], 
    FluxReconstructionType reconstructionType
) const {
    switch(reconstructionType) {
        case FluxReconstructionType::FirstOrder:
            return 1;
        case FluxReconstructionType::WENO:
            return 4;  // WENO5 requires 3 ghost cells
        default:
            return total_ghost_cells;
    }
}
    

// Add these implementations to IntegratorVariableAccessLayer.cpp

void 
CompressibleEulerVariableAccessor::CopyVariables(
    int direction [[maybe_unused]],
    int lev, 
    void* solver_void,
    amrex::MultiFab& VariableBuffer,
    ReconstructionMode mode, 
    const SolverCapabilities::MethodValidationResult& validationResult [[maybe_unused]]
) const {
    // Cast the void pointer to ScimitarX pointer
    auto solver = static_cast<Integrator::ScimitarX*>(solver_void);
    
    // Handle different reconstruction modes
    switch(mode) {
        case ReconstructionMode::Primitive:
            // Copy primitive variables and fluxes to the working buffer
            CopyPrimitiveVariablesToVariableBuffer(lev, solver, VariableBuffer);
            break;
            
        case ReconstructionMode::Conservative:
        case ReconstructionMode::Characteristic:
            // Copy conservative variables to the working buffer
            CopyConservativeVariablesToVariableBuffer(lev, solver, VariableBuffer);
            break;
                        
        default:
            Util::Abort(__FILE__, __func__, __LINE__, "Unsupported reconstruction mode");
    }
}

void 
CompressibleEulerVariableAccessor::CopyFluxes(
    int direction,
    int lev, 
    void* solver_void,
    amrex::MultiFab& CellFluxBuffer,
    ReconstructionMode mode, 
    const SolverCapabilities::MethodValidationResult& validationResult [[maybe_unused]]
) const {
    // Cast the void pointer to ScimitarX pointer
    auto solver = static_cast<Integrator::ScimitarX*>(solver_void);
    
    // Handle different reconstruction modes
    switch(mode) {
        case ReconstructionMode::Primitive:
            // Copy primitive variables and fluxes to the working buffer
            //CopyPrimitiveFluxesToCellFluxBuffer(direction, lev, solver, CellFluxBuffer);
            break;
            
        case ReconstructionMode::Conservative:
            // Copy conservative variables to the working buffer
            CopyConservativeFluxesToCellFluxBuffer(direction, lev, solver, CellFluxBuffer);
            break;

        case ReconstructionMode::Characteristic:
            // Copy conservative variables to the working buffer
            CopyConservativeFluxesToCellFluxBuffer(direction, lev, solver, CellFluxBuffer);
            break;

        default:
            Util::Abort(__FILE__, __func__, __LINE__, "Unsupported reconstruction mode");
    }
}

/*void CompressibleEulerVariableAccessor::CopyPrimitiveVariablesToVariableBuffer(
    int lev,
    Integrator::ScimitarX* solver,
    amrex::MultiFab& VariableBuffer
) const {
    // Copy from PVec_mf to working buffer
    amrex::MultiFab::Copy(
        VariableBuffer, 
        *(solver->PVec_mf[lev]),
        0, 0, 
        solver->number_of_components, 
        VariableBuffer.nGrow()
    );

    VariableBuffer.FillBoundary();
}*/

void 
CompressibleEulerVariableAccessor::CopyPrimitiveVariablesToVariableBuffer(
    int lev,
    Integrator::ScimitarX* solver,
    amrex::MultiFab& VariableBuffer
) const {
    // Loop through all valid regions of the MultiFab
    for (amrex::MFIter mfi(VariableBuffer); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        
        // Get data arrays for the primitive variables, pressure, and the destination buffer
        auto const& p_arr = solver->PVec_mf.Patch(lev, mfi);   // Primitive variables
        auto const& press_arr = solver->Pressure_mf.Patch(lev, mfi); // Pressure
        auto var_arr = VariableBuffer.array(mfi);             // Destination buffer
        
        // Get number of components
        int num_components = solver->number_of_components;
        
        // Parallel loop over the valid box
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Copy primitive variables (density, velocity components)
            for (int n = 0; n < num_components - 1; ++n) {
                var_arr(i, j, k, n) = p_arr(i, j, k, n);
            }
            
            // Copy pressure as the last component
            var_arr(i, j, k, num_components - 1) = press_arr(i, j, k);
        });
    }

    VariableBuffer.FillBoundary();
}


void 
CompressibleEulerVariableAccessor::CopyConservativeVariablesToVariableBuffer(
    int lev,
    Integrator::ScimitarX* solver,
    amrex::MultiFab& VariableBuffer
) const {
    // Copy from PVec_mf to working buffer
    amrex::MultiFab::Copy(
        VariableBuffer, 
        *(solver->QVec_mf[lev]),
        0, 0, 
        solver->number_of_components, 
        VariableBuffer.nGrow()
    );

    VariableBuffer.FillBoundary();
}


/*void CompressibleEulerVariableAccessor::CopyPrimitiveFluxVariablesToWorkingBuffer(
    int direction,    
    int lev,
    Integrator::ScimitarX* solver,
    amrex::MultiFab& CellFluxBuffer
) const {

    //Requires Primitive Flux Variable Implementation Here   
    CellFluxBuffer.FillBoundary();
}*/

void 
CompressibleEulerVariableAccessor::CopyConservativeFluxesToCellFluxBuffer(
    int direction,    
    int lev,
    Integrator::ScimitarX* solver,
    amrex::MultiFab& CellFluxBuffer
) const {

    for (amrex::MFIter mfi(CellFluxBuffer, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.growntilebox();
        auto const& p_arr = solver->PVec_mf.Patch(lev, mfi);
        // auto const& pres_arr = solver->Pressure_mf.Patch(lev, mfi);
        auto const& W_arr  = CellFluxBuffer.array(mfi);
        
        const Set::Scalar gamma = 1.4;  // Ratio of specific heats

        // Retrieve variable indices for better GPU performance
        int num_components = solver->number_of_components;   

        int rho_idx = solver->variableIndex.DENS;   // density
        int ie_idx  = solver->variableIndex.IE;     // energy
        int u_idx   = solver->variableIndex.UVEL;   // x-momentum

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


        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

            Set::Scalar rho = p_arr(i, j, k, rho_idx);
            
            Set::Scalar vel[3] = {0.0, 0.0, 0.0};
            vel[normal] = p_arr(i, j, k, velocity_component[normal]);
            vel[trans1] = p_arr(i, j, k, velocity_component[trans1]);
#if AMREX_SPACEDIM == 3                                                              
            vel[trans2] = p_arr(i, j, k, velocity_component[trans2]);
#endif

            Set::Scalar KE = 0.5 * (vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);

            Set::Scalar IE = p_arr(i, j, k, ie_idx);

            Set::Scalar E = IE + KE;

            Set::Scalar pressure = (gamma - 1.0) * rho * IE;


            Set::MultiVector Flux_Vector(num_components);
            Flux_Vector.setZero();

            Flux_Vector(rho_idx) = rho * vel[normal];
            Flux_Vector(u_idx + normal) = rho * vel[normal] * vel[normal] + pressure;
#if AMREX_SPACEDIM >= 2                
            Flux_Vector(u_idx + trans1) = rho * vel[normal] * vel[trans1];
#endif
#if AMREX_SPACEDIM == 3
            Flux_Vector(u_idx + trans2) = rho * vel[normal] * vel[trans2]; 
#endif            
            Flux_Vector(ie_idx) = (rho * E + pressure) * vel[normal];    
                         
            for (int n = 0; n < num_components; ++n) {

                W_arr(i, j, k, n) = Flux_Vector(n);

            }
                        
        });
    }

    
    CellFluxBuffer.FillBoundary();
}

void 
CompressibleEulerVariableAccessor::PopulateAverageStates(
    int direction, 
    int lev, 
    void* solver_void, 
    amrex::MultiFab& AverageStateBuffer
) const {

    // Cast the void pointer to ScimitarX pointer
    auto solver = static_cast<Integrator::ScimitarX*>(solver_void);

    const int nghosts = solver->number_of_ghost_cells; 

    for (amrex::MFIter mfi(AverageStateBuffer, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);
        auto const& p_arr = solver->PVec_mf.Patch(lev, mfi);
        auto const& W_arr  = AverageStateBuffer.array(mfi);
        
        // Retrieve variable indices for better GPU performance
        int num_components = 5;
        int solver_components = solver->number_of_components;   

        int rho_idx = solver->variableIndex.DENS;   // density
        int ie_idx  = solver->variableIndex.IE;     // energy
        int u_idx   = solver->variableIndex.UVEL;   // x-momentum
        int v_idx   = solver->variableIndex.VVEL;   // y-momentum
#if AMREX_SPACEDIM == 3
        int w_idx   = solver->variableIndex.WVEL;   // z-momentum
#endif
        
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
                
            // Extract original left state primitive variables
            Set::MultiVector UL = Numeric::FieldToMultiVector(p_arr, left_index[0], left_index[1], left_index[2], solver_components);
                        
            // Extract original right state primitive variables
            Set::MultiVector UR = Numeric::FieldToMultiVector(p_arr, right_index[0], right_index[1], right_index[2], solver_components);

            Set::MultiVector WL(num_components);
            Set::MultiVector WR(num_components);            
            Set::MultiVector Wavg(num_components);
            WL.setZero();
            WR.setZero();
            Wavg.setZero();

            WL(0) = UL(rho_idx);
            WL(1) = UL(u_idx);
            WL(2) = UL(v_idx);
#if AMREX_SPACEDIM == 3
            WL(3) = UL(w_idx);  
#else
            WL(3) = 0.0;
#endif
            WL(4) = UL(ie_idx);
            
            WR(0) = UR(rho_idx);
            WR(1) = UR(u_idx);
            WR(2) = UR(v_idx);
#if AMREX_SPACEDIM == 3
            WR(3) = UR(w_idx);  
#else
            WR(3) = 0.0;
#endif
            WR(4) = UR(ie_idx);
                
            Wavg = CompressibleEulerVariableAccessor::ComputeRoeAverages(WL, WR, num_components);

            W_arr(index[0], index[1], index[2], rho_idx) = Wavg(0);  // density 
            W_arr(index[0], index[1], index[2], u_idx)   = Wavg(1);  // x-momentum
            W_arr(index[0], index[1], index[2], v_idx)   = Wavg(2);  // y-momentum
#if AMREX_SPACEDIM == 3
            W_arr(index[0], index[1], index[2], w_idx)   = Wavg(3);  // z-momentum
#endif
            W_arr(index[0], index[1], index[2],ie_idx)   = Wavg(4);  // Total Enthalpy 
                        
        });
    }


    AverageStateBuffer.FillBoundary();
}

Set::MultiMatrix
CompressibleEulerVariableAccessor::TransformStencilToCharacteristic(
    int i, int j, int k,    
    const Set::MultiMatrix& stencil_matrix,
    int direction,
    const Set::MultiVector& avg_state
) const {
    
        int total_stencilpoints = stencil_matrix.rows();
        int num_components = stencil_matrix.cols();
        
        // Create output matrix of same dimensions
        Set::MultiMatrix char_stencil_matrix(total_stencilpoints, num_components);
        char_stencil_matrix.setZero();
        
        // Use the static variable indices
        const int rho_idx = variableIndices.DENS;
        const int u_idx = variableIndices.UVEL;
        const int v_idx = variableIndices.VVEL;
        const int w_idx = variableIndices.WVEL;
        const int ie_idx = variableIndices.IE;
        
        // Create standardized vector for the eigenvector calculation
        Set::MultiVector std_avg_state(5); // Always 5 components for eigenvector calculation
        std_avg_state(0) = avg_state(rho_idx);
        std_avg_state(1) = avg_state(u_idx);
        std_avg_state(2) = avg_state(v_idx);
        std_avg_state(3) = (w_idx >= 0) ? avg_state(w_idx) : 0.0;
        std_avg_state(4) = avg_state(ie_idx);

            // Debug the average state first
        Util::ScimitarX_Util::Debug::DebugAverageState(
        i, j, k,  // These will be ignored if not the target location
        Set::MultiVector::Zero(avg_state.size()),  // We don't have WL/WR here
        Set::MultiVector::Zero(avg_state.size()),  // We don't have WL/WR here
        avg_state,
        "Transform to Characteristic",
        false,    // Don't abort, just warn
        false);    // Enable debug
        
        // Get left eigenvector matrix for the transformation
        Set::MultiMatrix L_n = ComputeLeftEigenvectorMatrix(i, j, k, std_avg_state, direction, 5);
        
        // Transform each stencil point to characteristic variables
        for (int s = 0; s < total_stencilpoints; ++s) {
            // Map stencil values to standardized vector
            Set::MultiVector std_state_vec(5);
            std_state_vec(0) = stencil_matrix(s, rho_idx);
            std_state_vec(1) = stencil_matrix(s, u_idx);
            std_state_vec(2) = stencil_matrix(s, v_idx);
            std_state_vec(3) = (w_idx >= 0) ? stencil_matrix(s, w_idx) : 0.0;
            std_state_vec(4) = stencil_matrix(s, ie_idx);
            
            // Transform to characteristic space
            Set::MultiVector char_vec = Numeric::SymmetryPreserving::ConsistentMatrixVectorMultiply(
                L_n, std_state_vec);
            
            // Map back to solver's indexing scheme
            char_stencil_matrix(s, rho_idx) = char_vec(0);
            char_stencil_matrix(s, u_idx) = char_vec(1);
            char_stencil_matrix(s, v_idx) = char_vec(2);
            if (w_idx >= 0) char_stencil_matrix(s, w_idx) = char_vec(3);
            char_stencil_matrix(s, ie_idx) = char_vec(4);
       
            // Debug each transformation
            Util::ScimitarX_Util::Debug::DebugCharacteristicTransformation(
            i, j, k,  // These will be ignored if not the target location
            stencil_matrix, char_stencil_matrix, L_n, s,
            "Transform Stencil Point " + std::to_string(s), 
            false,    // Don't abort, just warn
            false);    // Enable debug       
       
        }
        
        return char_stencil_matrix;
}

Set::MultiVector
CompressibleEulerVariableAccessor::TransformFromCharacteristic(
    const Set::MultiVector& char_vec,
    int direction,
    const Set::MultiVector& avg_state
) const {
        // Create output vector of the same type and size
        Set::MultiVector phys_vec(char_vec.size());
         
        // Use the static variable indices
        const int rho_idx = variableIndices.DENS;
        const int u_idx = variableIndices.UVEL;
        const int v_idx = variableIndices.VVEL;
        const int w_idx = variableIndices.WVEL;
        const int ie_idx = variableIndices.IE;
         
        // Create standardized vectors
        Set::MultiVector std_avg_state(5);
        std_avg_state(0) = avg_state(rho_idx);
        std_avg_state(1) = avg_state(u_idx);
        std_avg_state(2) = avg_state(v_idx);
        std_avg_state(3) = (w_idx >= 0) ? avg_state(w_idx) : 0.0;
        std_avg_state(4) = avg_state(ie_idx);
         
        Set::MultiVector std_char_vec(5);
        std_char_vec(0) = char_vec(rho_idx);
        std_char_vec(1) = char_vec(u_idx);
        std_char_vec(2) = char_vec(v_idx);
        std_char_vec(3) = (w_idx >= 0) ? char_vec(w_idx) : 0.0;
        std_char_vec(4) = char_vec(ie_idx);
         
        // Get right eigenvector matrix for the inverse transformation
        Set::MultiMatrix R_n = ComputeRightEigenvectorMatrix(std_avg_state, direction, 5);
         
        // Transform back to physical space
        Set::MultiVector std_phys_vec = Numeric::SymmetryPreserving::ConsistentMatrixVectorMultiply(
            R_n, std_char_vec);
         
        // Map back to solver's indexing scheme
        phys_vec(rho_idx) = std_phys_vec(0);
        phys_vec(u_idx) = std_phys_vec(1);
        phys_vec(v_idx) = std_phys_vec(2);
        if (w_idx >= 0) phys_vec(w_idx) = std_phys_vec(3);
        phys_vec(ie_idx) = std_phys_vec(4);
         
        return phys_vec;
}


void 
CompressibleEulerVariableAccessor::StoreDirectionalFlux(
    int direction, 
    int lev, 
    void* solver_void,
    amrex::MultiFab& SummedFlux
) const {

    // Cast the void pointer to ScimitarX pointer
    auto solver = static_cast<Integrator::ScimitarX*>(solver_void);

    const int nghosts = solver->number_of_ghost_cells;
    const int num_components = solver->number_of_components; 

    for (amrex::MFIter mfi(SummedFlux, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.grownnodaltilebox(direction, nghosts);
        auto const& TotalFlux_arr = SummedFlux.array(mfi);
        auto const& flux_arr = (direction == Directions::Xdir) ? solver->XFlux_mf.Patch(lev, mfi) :
                                (direction == Directions::Ydir) ? solver->YFlux_mf.Patch(lev, mfi) :
                                solver->ZFlux_mf.Patch(lev, mfi);
        
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {

    
            // Construct index arrays with correct number of components
            int index[3];
            int lower_bounds[3];
            int upper_bounds[3];

            // Set up indices
            index[0] = iface; index[1] = jface; index[2] = kface;

            // Define bounds including ghost cells
            lower_bounds[0] = bx.smallEnd(0); 
            lower_bounds[1] = bx.smallEnd(1); 
            lower_bounds[2] = AMREX_SPACEDIM == 3 ? bx.smallEnd(2) : 0;
            upper_bounds[0] = bx.bigEnd(0) - 1;
            upper_bounds[1] = bx.bigEnd(1) - 1;
            upper_bounds[2] = AMREX_SPACEDIM == 3 ? bx.bigEnd(2) - 1 : 0;
                  
      
            // Clamp the main index
            ClampIndices(index, lower_bounds, upper_bounds);

            for (int n = 0; n < num_components; ++n) { 
            flux_arr(index[0], index[1], index[2], n) = TotalFlux_arr(index[0], index[1], index[2], n);
            }

        });

    } 
}


Set::MultiMatrix
CompressibleEulerVariableAccessor::ComputeRightEigenvectorMatrix(
    const Set::MultiVector& W, int dir, int num_components) {
    // Constants
    const Set::Scalar gamma = 1.4;
    const Set::Scalar g = gamma - 1.0;  // g = gamma - 1.0

    // Extract primitive variables from W with consistent ordering
    const Set::Scalar u = W(1);
    const Set::Scalar v = W(2);
    const Set::Scalar w = (num_components > 3 && AMREX_SPACEDIM == 3) ? W(3) : 0.0;
    const Set::Scalar H = W(4);  // Total enthalpy

    // Derived quantities with consistent ordering
    // Use symmetry-preserving operations for velocity magnitude calculations
    const Set::Scalar Va2 = Numeric::SymmetryPreserving::ConsistentVelocityMagnitudeSquared(u, v, w);

    // Sound speed with consistent ordering
    const Set::Scalar a2 = g * (H - 0.5*Va2);  // Sound speed squared
    const Set::Scalar a = std::sqrt(a2);       // Sound speed

    // Initialize right eigenvector matrix
    Set::MultiMatrix R_n = Set::MultiMatrix::Zero(num_components, num_components);

    // Fill matrix based on direction with consistent ordering of operations
    switch(dir) {
    case 0:  // X-direction
        // Column 1: First acoustic wave
        R_n(0, 0) = 1.0;
        R_n(1, 0) = u - a;
        R_n(2, 0) = v;
        R_n(3, 0) = w;
        R_n(4, 0) = H - u*a;

        // Column 2: Entropy wave
        R_n(0, 1) = 1.0;
        R_n(1, 1) = u;
        R_n(2, 1) = v;
        R_n(3, 1) = w;
        R_n(4, 1) = 0.5*Va2;

        // Column 3: First shear wave
        R_n(0, 2) = 0.0;
        R_n(1, 2) = 0.0;
        R_n(2, 2) = 1.0;
        R_n(3, 2) = 0.0;
        R_n(4, 2) = v;

        // Column 4: Second shear wave
        R_n(0, 3) = 0.0;
        R_n(1, 3) = 0.0;
        R_n(2, 3) = 0.0;
        R_n(3, 3) = 1.0;
        R_n(4, 3) = w;

        // Column 5: Second acoustic wave
        R_n(0, 4) = 1.0;
        R_n(1, 4) = u + a;
        R_n(2, 4) = v;
        R_n(3, 4) = w;
        R_n(4, 4) = H + u*a;
        break;

    case 1:  // Y-direction
        // Column 1: First acoustic wave
        R_n(0, 0) = 1.0;
        R_n(1, 0) = u;
        R_n(2, 0) = v - a;
        R_n(3, 0) = w;
        R_n(4, 0) = H - v*a;

        // Column 2: First shear wave
        R_n(0, 1) = 0.0;
        R_n(1, 1) = 1.0;
        R_n(2, 1) = 0.0;
        R_n(3, 1) = 0.0;
        R_n(4, 1) = u;

        // Column 3: Entropy wave
        R_n(0, 2) = 1.0;
        R_n(1, 2) = u;
        R_n(2, 2) = v;
        R_n(3, 2) = w;
        R_n(4, 2) = 0.5*Va2;

        // Column 4: Second shear wave
        R_n(0, 3) = 0.0;
        R_n(1, 3) = 0.0;
        R_n(2, 3) = 0.0;
        R_n(3, 3) = 1.0;
        R_n(4, 3) = w;

        // Column 5: Second acoustic wave
        R_n(0, 4) = 1.0;
        R_n(1, 4) = u;
        R_n(2, 4) = v + a;
        R_n(3, 4) = w;
        R_n(4, 4) = H + v*a;
        break;

    case 2:  // Z-direction
        // Column 1: First acoustic wave
        R_n(0, 0) = 1.0;
        R_n(1, 0) = u;
        R_n(2, 0) = v;
        R_n(3, 0) = w - a;
        R_n(4, 0) = H - w*a;

        // Column 2: First shear wave
        R_n(0, 1) = 0.0;
        R_n(1, 1) = 1.0;
        R_n(2, 1) = 0.0;
        R_n(3, 1) = 0.0;
        R_n(4, 1) = u;

        // Column 3: Second shear wave
        R_n(0, 2) = 0.0;
        R_n(1, 2) = 0.0;
        R_n(2, 2) = 1.0;
        R_n(3, 2) = 0.0;
        R_n(4, 2) = v;

        // Column 4: Entropy wave
        R_n(0, 3) = 1.0;
        R_n(1, 3) = u;
        R_n(2, 3) = v;
        R_n(3, 3) = w;
        R_n(4, 3) = 0.5*Va2;

        // Column 5: Second acoustic wave
        R_n(0, 4) = 1.0;
        R_n(1, 4) = u;
        R_n(2, 4) = v;
        R_n(3, 4) = w + a;
        R_n(4, 4) = H + w*a;
        break;
    }

    return R_n;
}

Set::MultiMatrix 
CompressibleEulerVariableAccessor::ComputeLeftEigenvectorMatrix(
    int i, int j, int k,    
    const Set::MultiVector& W, int dir, int num_components) {
    // Constants
    const Set::Scalar gamma = 1.4;
    const Set::Scalar g = gamma - 1.0;  // g = gamma - 1.0
    
    // Extract primitive variables from W with consistent ordering
    const Set::Scalar u = W(1);
    const Set::Scalar v = W(2);
    const Set::Scalar w = (num_components > 3 && AMREX_SPACEDIM == 3) ? W(3) : 0.0;
    const Set::Scalar H = W(4);  // Total enthalpy

    // Derived quantities with consistent ordering using symmetry-preserving operations
    // Use the consistent sum function for velocity magnitude calculation
    const Set::Scalar Va2 = Numeric::SymmetryPreserving::ConsistentVelocityMagnitudeSquared(u, v, w);
    
    // Sound speed with consistent ordering
    const Set::Scalar a2 = g * (H - 0.5*Va2);  // Sound speed squared
    const Set::Scalar a = std::sqrt(a2);       // Sound speed
    
    // Add debug check
    Util::ScimitarX_Util::Debug::DebugSoundSpeedComputation(
        i, j, k,  // These will be ignored if not the target location
        H, Va2, gamma,
        a2, a,    // Pass by reference to allow correction
        "Left Eigenvector Computation", 
        false,    // Don't abort, just warn
        false);    // Enable this debug

    // Scaling factor used in eigenvector computation
    const Set::Scalar M = (0.5*g)/(a2);  // Scaling factor 
    
    // Initialize left eigenvector matrix
    Set::MultiMatrix L_n = Set::MultiMatrix::Zero(num_components, num_components);
    
    // Fill matrix based on direction with consistent ordering of operations
    switch(dir) {
    case 0:  // X-direction
        // Row 1: First acoustic wave - λ₁
        L_n(0, 0) = M*(H + a*(u-a)/g);
        L_n(0, 1) = -M*(u + a/g);
        L_n(0, 2) = -M*v;
        L_n(0, 3) = -M*w;
        L_n(0, 4) = M;
        
        // Row 2: Entropy wave - λ₂
        L_n(1, 0) = M*(-2.0*H + 4.0*a2/g);
        L_n(1, 1) = M*2.0*u;
        L_n(1, 2) = M*2.0*v;
        L_n(1, 3) = M*2.0*w;
        L_n(1, 4) = -2.0*M;
        
        // Row 3: First shear wave - λ₃
        L_n(2, 0) = -M*2.0*v*a2/g;
        L_n(2, 1) = 0.0;
        L_n(2, 2) = M*2.0*a2/g;
        L_n(2, 3) = 0.0;
        L_n(2, 4) = 0.0;
        
        // Row 4: Second shear wave - λ₄
        L_n(3, 0) = -M*2.0*w*a2/g;
        L_n(3, 1) = 0.0;
        L_n(3, 2) = 0.0;
        L_n(3, 3) = M*2.0*a2/g;
        L_n(3, 4) = 0.0;
        
        // Row 5: Second acoustic wave - λ₅
        L_n(4, 0) = M*(H - a*(u+a)/g);
        L_n(4, 1) = M*(-u + a/g);
        L_n(4, 2) = -M*v;
        L_n(4, 3) = -M*w;
        L_n(4, 4) = M;
        break;
        
    case 1:  // Y-direction
        // Row 1: First acoustic wave - λ₁
        L_n(0, 0) = M*(H + a*(v-a)/g);
        L_n(0, 1) = -M*u;
        L_n(0, 2) = -M*(v + a/g);
        L_n(0, 3) = -M*w;
        L_n(0, 4) = M;
        
        // Row 2: First shear wave - λ₂
        L_n(1, 0) = -M*2.0*u*a2/g;
        L_n(1, 1) = M*2.0*a2/g;
        L_n(1, 2) = 0.0;
        L_n(1, 3) = 0.0;
        L_n(1, 4) = 0.0;
        
        // Row 3: Entropy wave - λ₃
        L_n(2, 0) = M*(-2.0*H + 4.0*a2/g);
        L_n(2, 1) = M*2.0*u;
        L_n(2, 2) = M*2.0*v;
        L_n(2, 3) = M*2.0*w;
        L_n(2, 4) = -2.0*M;
        
        // Row 4: Second shear wave - λ₄
        L_n(3, 0) = -M*2.0*w*a2/g;
        L_n(3, 1) = 0.0;
        L_n(3, 2) = 0.0;
        L_n(3, 3) = M*2.0*a2/g;
        L_n(3, 4) = 0.0;
        
        // Row 5: Second acoustic wave - λ₅
        L_n(4, 0) = M*(H - a*(v+a)/g);
        L_n(4, 1) = -M*u;
        L_n(4, 2) = M*(-v + a/g);
        L_n(4, 3) = -M*w;
        L_n(4, 4) = M;
        break;
        
    case 2:  // Z-direction
        // Row 1: First acoustic wave - λ₁
        L_n(0, 0) = M*(H + a*(w-a)/g);
        L_n(0, 1) = -M*u;
        L_n(0, 2) = -M*v;
        L_n(0, 3) = -M*(w + a/g);
        L_n(0, 4) = M;
        
        // Row 2: First shear wave - λ₂
        L_n(1, 0) = -M*2.0*u*a2/g;
        L_n(1, 1) = M*2.0*a2/g;
        L_n(1, 2) = 0.0;
        L_n(1, 3) = 0.0;
        L_n(1, 4) = 0.0;
        
        // Row 3: Second shear wave - λ₃
        L_n(2, 0) = -M*2.0*v*a2/g;
        L_n(2, 1) = 0.0;
        L_n(2, 2) = M*2.0*a2/g;
        L_n(2, 3) = 0.0;
        L_n(2, 4) = 0.0;
        
        // Row 4: Entropy wave - λ₄
        L_n(3, 0) = M*(-2.0*H + 4.0*a2/g);
        L_n(3, 1) = M*2.0*u;
        L_n(3, 2) = M*2.0*v;
        L_n(3, 3) = M*2.0*w;
        L_n(3, 4) = -2.0*M;
        
        // Row 5: Second acoustic wave - λ₅
        L_n(4, 0) = M*(H - a*(w+a)/g);
        L_n(4, 1) = -M*u;
        L_n(4, 2) = -M*v;
        L_n(4, 3) = M*(-w + a/g);
        L_n(4, 4) = M;
        break;
    }

    // Debug check on the final matrix
    Util::ScimitarX_Util::Debug::DebugEigenvectorMatrix(
        i, j, k,  // These will be ignored if not the target location
        L_n, W, 
        "Left Eigenvector Computation", 
        false,    // Don't abort, just warn
        false);    // Enable this debug

    return L_n;
}

    
Set::MultiVector
CompressibleEulerVariableAccessor::ComputeRoeAverages(
    const Set::MultiVector& WL, 
    const Set::MultiVector& WR,
    int num_components)
{
    // Simply delegate to the symmetry-preserving templated implementation
    return Numeric::SymmetryPreserving::ConsistentRoeAverage(WL, WR, num_components);
}        

} // namespace CompressibleEuler


namespace ElastoPlastic {    

ElastoPlasticVariableAccessor::ElastoPlasticVariableAccessor(
    int total_ghosts,
    ReconstructionMode mode
) : current_mode(mode), total_ghost_cells(total_ghosts) {}

std::shared_ptr<SolverCapabilities> 
ElastoPlasticVariableAccessor::getSolverCapabilities() const {
    // Create a detailed solver capabilities object for ElastoPlastic solver
    return std::make_shared<ElastoPlastic::ElastoPlasticCapabilities>();
}

amrex::MultiFab 
ElastoPlasticVariableAccessor::CreateWorkingBuffer(
    const amrex::BoxArray& baseGrids, 
    const amrex::DistributionMapping& dm, 
    const int num_components, 
    const int ghost_cells
) const {

    return amrex::MultiFab(
        baseGrids, 
        dm, 
        num_components, 
        ghost_cells
    );
}

// Implementations of other methods remain similar to previous version
// with added validation and more flexible handling

int 
ElastoPlasticVariableAccessor::getRequiredGhostCells(
    ReconstructionMode mode [[maybe_unused]], 
    FluxReconstructionType reconstructionType
) const {
    switch(reconstructionType) {
        case FluxReconstructionType::FirstOrder:
            return 1;
        case FluxReconstructionType::WENO:
            return 4;  // WENO5 requires 3 ghost cells
        default:
            return total_ghost_cells;
    }
}
    

// Add these implementations to IntegratorVariableAccessLayer.cpp

void 
ElastoPlasticVariableAccessor::CopyVariables(
    int direction [[maybe_unused]],
    int lev, 
    void* solver_void,
    amrex::MultiFab& VariableBuffer,
    ReconstructionMode mode, 
    const SolverCapabilities::MethodValidationResult& validationResult [[maybe_unused]]
) const {
    // Cast the void pointer to ScimitarX pointer
    auto solver = static_cast<Integrator::ScimitarX*>(solver_void);
    
    // Handle different reconstruction modes
    switch(mode) {
        case ReconstructionMode::Primitive:
            // Copy primitive variables and fluxes to the working buffer
            CopyPrimitiveVariablesToVariableBuffer(lev, solver, VariableBuffer);
            break;
            
        case ReconstructionMode::Conservative:
        case ReconstructionMode::Characteristic:
            // Copy conservative variables to the working buffer
            CopyConservativeVariablesToVariableBuffer(lev, solver, VariableBuffer);
            break;
                        
        default:
            Util::Abort(__FILE__, __func__, __LINE__, "Unsupported reconstruction mode");
    }
}

void 
ElastoPlasticVariableAccessor::CopyFluxes(
    int direction,
    int lev, 
    void* solver_void,
    amrex::MultiFab& CellFluxBuffer,
    ReconstructionMode mode, 
    const SolverCapabilities::MethodValidationResult& validationResult [[maybe_unused]]
) const {
    // Cast the void pointer to ScimitarX pointer
    auto solver = static_cast<Integrator::ScimitarX*>(solver_void);
    
    // Handle different reconstruction modes
    switch(mode) {
        case ReconstructionMode::Primitive:
            // Copy primitive variables and fluxes to the working buffer
            //CopyPrimitiveFluxesToCellFluxBuffer(direction, lev, solver, CellFluxBuffer);
            break;
            
        case ReconstructionMode::Conservative:
            // Copy conservative variables to the working buffer
            CopyConservativeFluxesToCellFluxBuffer(direction, lev, solver, CellFluxBuffer);
            break;

        case ReconstructionMode::Characteristic:
            // Copy conservative variables to the working buffer
            CopyConservativeFluxesToCellFluxBuffer(direction, lev, solver, CellFluxBuffer);
            break;

        default:
            Util::Abort(__FILE__, __func__, __LINE__, "Unsupported reconstruction mode");
    }
}

/*void CompressibleEulerVariableAccessor::CopyPrimitiveVariablesToVariableBuffer(
    int lev,
    Integrator::ScimitarX* solver,
    amrex::MultiFab& VariableBuffer
) const {
    // Copy from PVec_mf to working buffer
    amrex::MultiFab::Copy(
        VariableBuffer, 
        *(solver->PVec_mf[lev]),
        0, 0, 
        solver->number_of_components, 
        VariableBuffer.nGrow()
    );

    VariableBuffer.FillBoundary();
}*/

void 
ElastoPlasticVariableAccessor::CopyPrimitiveVariablesToVariableBuffer(
    int lev,
    Integrator::ScimitarX* solver,
    amrex::MultiFab& VariableBuffer
) const {
    // Loop through all valid regions of the MultiFab
    for (amrex::MFIter mfi(VariableBuffer); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        
        // Get data arrays for the primitive variables, pressure, and the destination buffer
        auto const& p_arr = solver->PVec_mf.Patch(lev, mfi);   // Primitive variables
        auto const& press_arr = solver->Pressure_mf.Patch(lev, mfi); // Pressure
        auto var_arr = VariableBuffer.array(mfi);             // Destination buffer
        
        // Get number of components
        int num_components = solver->number_of_components;
        
        // Parallel loop over the valid box
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Copy primitive variables (density, velocity components)
            for (int n = 0; n < num_components - 1; ++n) {
                var_arr(i, j, k, n) = p_arr(i, j, k, n);
            }
            
            // Copy pressure as the last component
            var_arr(i, j, k, num_components - 1) = press_arr(i, j, k);
        });
    }

    VariableBuffer.FillBoundary();
}


void 
ElastoPlasticVariableAccessor::CopyConservativeVariablesToVariableBuffer(
    int lev,
    Integrator::ScimitarX* solver,
    amrex::MultiFab& VariableBuffer
) const {
    // Copy from PVec_mf to working buffer
    amrex::MultiFab::Copy(
        VariableBuffer, 
        *(solver->QVec_mf[lev]),
        0, 0, 
        solver->number_of_components, 
        VariableBuffer.nGrow()
    );

    VariableBuffer.FillBoundary();
}


/*void CompressibleEulerVariableAccessor::CopyPrimitiveFluxVariablesToWorkingBuffer(
    int direction,    
    int lev,
    Integrator::ScimitarX* solver,
    amrex::MultiFab& CellFluxBuffer
) const {

    //Requires Primitive Flux Variable Implementation Here   
    CellFluxBuffer.FillBoundary();
}*/

void 
ElastoPlasticVariableAccessor::CopyConservativeFluxesToCellFluxBuffer(
    int direction,    
    int lev,
    Integrator::ScimitarX* solver,
    amrex::MultiFab& CellFluxBuffer
) const {

//        Util::ScimitarX_Util::Debug debug;

//        debug.SetTargetDebugLocationIndices(3, 3, 0);
    for (amrex::MFIter mfi(CellFluxBuffer, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.growntilebox();
        auto const& p_arr = solver->PVec_mf.Patch(lev, mfi);
        // auto const& pres_arr = solver->Pressure_mf.Patch(lev, mfi);
        auto const& W_arr  = CellFluxBuffer.array(mfi);
        
        // Retrieve variable indices for better GPU performance
        int num_components = solver->number_of_components;   

        int rho_idx = solver->variableIndex.DENS;   // density
        int ie_idx  = solver->variableIndex.IE;     // energy
        int u_idx   = solver->variableIndex.UVEL;   // x-momentum

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

	int sxx_idx = solver->variableIndex.SXX;
	int syy_idx = solver->variableIndex.SYY;
	int szz_idx = solver->variableIndex.SZZ;
	int sxy_idx = solver->variableIndex.SXY;
#if AMREX_SPACEDIM == 3
	int sxz_idx = solver->variableIndex.SXZ;
	int syz_idx = solver->variableIndex.SYZ;
#endif
	int epsbar_idx = solver->variableIndex.EPSBAR;
	int ie_elastic_idx = solver->variableIndex.IE_ELASTIC;


        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
			
	       Set::Scalar rho = p_arr(i, j, k, rho_idx);
	       Set::Scalar vel[3] = {0.0, 0.0, 0.0};
	       vel[normal] = p_arr(i, j, k, velocity_component[normal]);
	       vel[trans1] = p_arr(i, j, k, velocity_component[trans1]);
#if AMREX_SPACEDIM == 3
	       vel[trans2] = p_arr(i, j, k, velocity_component[trans2]);
#endif

            Set::Scalar KE = 0.5 * (vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);

            Set::Scalar IE = p_arr(i, j, k, ie_idx);

            Set::Scalar E = IE + KE;

	    Set::Scalar pressure = Model::ElastoPlastic::ElastoPlastic().ComputePressureFromDensityAndInternalEnergy(rho, IE);


	    Set::Scalar sxx = p_arr(i, j, k, sxx_idx);
	    Set::Scalar syy = p_arr(i, j, k, syy_idx);
	    Set::Scalar szz = p_arr(i, j, k, szz_idx);
	    Set::Scalar sxy = p_arr(i, j, k, sxy_idx);
#if AMREX_SPACEDIM == 3
	    Set::Scalar sxz = p_arr(i, j, k, sxz_idx);
	    Set::Scalar syz = p_arr(i, j, k, syz_idx);
#endif
	    Set::Scalar epsbar = p_arr(i, j, k, epsbar_idx);
	    Set::Scalar ie_elastic = p_arr(i, j, k, ie_elastic_idx);

/*
                debug.DebugValuesIfTarget(i, j, k, rho, "Density", "Variable access", false, -1, true);
                debug.DebugValuesIfTarget(i, j, k, vel[normal], "UVEL", "Variable access", false, -1, true);
                debug.DebugValuesIfTarget(i, j, k, vel[trans1], "UVEL", "Variable access", false, -1, true);
                debug.DebugValuesIfTarget(i, j, k, pressure, "Pressure", "UpdateSolutions", false, -1, true);
                debug.DebugValuesIfTarget(i, j, k, IE, "internal_energy", "UpdateSolutions", false, -1, true);
                debug.DebugValuesIfTarget(i, j, k, sxx, "SXX", "UpdateSolutions", false, -1, true);
                debug.DebugValuesIfTarget(i, j, k, syy, "SYY", "UpdateSolutions", false, -1, true);
                debug.DebugValuesIfTarget(i, j, k, sxy, "SXY", "UpdateSolutions", false, -1, true);
*/
            Set::MultiVector Flux_Vector(num_components);
            Flux_Vector.setZero();

            Flux_Vector(rho_idx) = rho * vel[normal];
            Flux_Vector(u_idx + normal) = rho * vel[normal] * vel[normal] + pressure;
#if AMREX_SPACEDIM >= 2                
            Flux_Vector(u_idx + trans1) = rho * vel[normal] * vel[trans1];
#endif
#if AMREX_SPACEDIM == 3
            Flux_Vector(u_idx + trans2) = rho * vel[normal] * vel[trans2]; 
#endif            
            Flux_Vector(ie_idx) = (rho * E + pressure) * vel[normal];

	    //Elastoplastic
	    Flux_Vector(sxx_idx) = rho * vel[normal] * sxx;
	    Flux_Vector(syy_idx) = rho * vel[normal] * syy;
	    Flux_Vector(sxy_idx) = rho * vel[normal] * sxy;
	    Flux_Vector(szz_idx) = rho * vel[normal] * szz;
#if AMREX_SPACEDIM == 3
	    Flux_Vector(sxz_idx) = rho * vel[normal] * sxz;
	    Flux_Vector(syz_idx) = rho * vel[normal] * syz;
#endif
	    Flux_Vector(epsbar_idx) = 0.0;
	    Flux_Vector(ie_elastic_idx) = 0.0;

            for (int n = 0; n < num_components; ++n) {

                W_arr(i, j, k, n) = Flux_Vector(n);

            }
                        
        });
    }

    
    CellFluxBuffer.FillBoundary();
}


void 
ElastoPlasticVariableAccessor::PopulateAverageStates(
    int direction, 
    int lev, 
    void* solver_void, 
    amrex::MultiFab& AverageStateBuffer
) const {

}


Set::MultiMatrix
ElastoPlasticVariableAccessor::TransformStencilToCharacteristic(
    int i, int j, int k,    
    const Set::MultiMatrix& stencil_matrix,
    int direction,
    const Set::MultiVector& avg_state
) const {
}


Set::MultiVector
ElastoPlasticVariableAccessor::TransformFromCharacteristic(
    const Set::MultiVector& char_vec,
    int direction,
    const Set::MultiVector& avg_state
) const {
}

void 
ElastoPlasticVariableAccessor::StoreDirectionalFlux(
    int direction, 
    int lev, 
    void* solver_void,
    amrex::MultiFab& SummedFlux
) const {

    // Cast the void pointer to ScimitarX pointer
    auto solver = static_cast<Integrator::ScimitarX*>(solver_void);

    const int nghosts = solver->number_of_ghost_cells;
    const int num_components = solver->number_of_components; 

    for (amrex::MFIter mfi(SummedFlux, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.grownnodaltilebox(direction, nghosts);
        auto const& TotalFlux_arr = SummedFlux.array(mfi);
        auto const& flux_arr = (direction == Directions::Xdir) ? solver->XFlux_mf.Patch(lev, mfi) :
                                (direction == Directions::Ydir) ? solver->YFlux_mf.Patch(lev, mfi) :
                                solver->ZFlux_mf.Patch(lev, mfi);
        
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {

    
            // Construct index arrays with correct number of components
            int index[3];
            int lower_bounds[3];
            int upper_bounds[3];

            // Set up indices
            index[0] = iface; index[1] = jface; index[2] = kface;

            // Define bounds including ghost cells
            lower_bounds[0] = bx.smallEnd(0); 
            lower_bounds[1] = bx.smallEnd(1); 
            lower_bounds[2] = AMREX_SPACEDIM == 3 ? bx.smallEnd(2) : 0;
            upper_bounds[0] = bx.bigEnd(0) - 1;
            upper_bounds[1] = bx.bigEnd(1) - 1;
            upper_bounds[2] = AMREX_SPACEDIM == 3 ? bx.bigEnd(2) - 1 : 0;
                  
      
            // Clamp the main index
            ClampIndices(index, lower_bounds, upper_bounds);

            for (int n = 0; n < num_components; ++n) { 
            flux_arr(index[0], index[1], index[2], n) = TotalFlux_arr(index[0], index[1], index[2], n);
            }

        });

    } 
}

} // namespace ElastoPlastic

}  // namespace Numeric 

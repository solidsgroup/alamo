#include "Numeric/IntegratorVariableAccessLayer.H"
#include "Integrator/ScimitarX.H"
#include "Numeric/Stencil.H"

namespace Numeric {

namespace CompressibleEuler {    

CompressibleEulerVariableAccessor::CompressibleEulerVariableAccessor(
    int total_ghosts,
    ReconstructionMode mode
) : current_mode(mode), total_ghost_cells(total_ghosts) {}

std::shared_ptr<SolverCapabilities> 
CompressibleEulerVariableAccessor::getSolverCapabilities() const {
    // Create a detailed solver capabilities object for Compressible Euler solver
    class CompressibleEulerCapabilities : public SolverCapabilities {
    public:
        std::string getIdentifier() const override { 
            return "CompressibleEuler"; 
        }

        std::string getDescription() const override { 
            return "Compressible Euler Equations Solver"; 
        }

        MethodSupport supportsFluxReconstruction(FluxReconstructionType method) const override {
            switch(method) {
                case FluxReconstructionType::FirstOrder:
                    return MethodSupport::Supported();
                case FluxReconstructionType::WENO:
                    return MethodSupport::Supported();
                default:
                    return MethodSupport::Unsupported(
                        {"Unsupported reconstruction method"},
                        {"Use FirstOrder or WENO reconstruction"}
                    );
            }
        }

        MethodSupport supportsFluxScheme(FluxScheme scheme) const override {
            switch(scheme) {
                case FluxScheme::LocalLaxFriedrichs:
                    return MethodSupport::Supported();
                case FluxScheme::HLLC:
                    return MethodSupport::Supported();
                default:
                    return MethodSupport::Unsupported(
                        {"Unknown flux scheme"},
                        {"Use standard flux methods"}
                    );
            }
        }

        MethodSupport supportsTimeSteppingScheme(TimeSteppingSchemeType scheme) const override {
            switch(scheme) {
                case TimeSteppingSchemeType::ForwardEuler:
                case TimeSteppingSchemeType::RK3:
                    return MethodSupport::Supported();
                default:
                    return MethodSupport::Unsupported(
                        {"Unsupported time stepping scheme"},
                        {"Use ForwardEuler or RK3"}
                    );
            }
        }

        MethodSupport supportsReconstructionMode(ReconstructionMode mode) const override {
            switch(mode) {
                case ReconstructionMode::Primitive:
                case ReconstructionMode::Conservative:
                case ReconstructionMode::Characteristic:
                    return MethodSupport::Supported();
                default:
                    return MethodSupport::Unsupported(
                        {"Unsupported reconstruction mode"},
                        {"Use Primitive, Conservative, or Characteristic"}
                    );
            }
        }

        MethodSupport supportsWenoVariant(WenoVariant variant) const override {
            switch(variant) {
                case WenoVariant::WENOJS5:
                case WenoVariant::WENOZ5:
                    return MethodSupport::Supported();
                default:
                    return MethodSupport::Unsupported(
                        {"Unsupported WENO variant"},
                        {"Use WENOJS5 or WENOZ5"}
                    );
            }
        }

        MethodValidationResult validateMethodCombination(
            FluxReconstructionType fluxReconstruction,
            FluxScheme fluxScheme,
            TimeSteppingSchemeType timeSteppingScheme,
            ReconstructionMode reconstructionMode,
            WenoVariant wenoVariant = WenoVariant::WENOJS5
        ) const override {
            MethodValidationResult result;
            result.isValid = true;

            // Add specific validation logic here
            if (fluxReconstruction == FluxReconstructionType::WENO && 
                fluxScheme == FluxScheme::HLLC) {
                result.isValid = false;
                result.warnings.push_back(
                    "WENO reconstruction with HLLC flux might be numerically unstable"
                );
            }

            return result;
        }

        DefaultConfiguration getDefaultConfiguration() const override {
            return {
                FluxReconstructionType::WENO,
                FluxScheme::LocalLaxFriedrichs,
                TimeSteppingSchemeType::RK3,
                ReconstructionMode::Primitive,
                WenoVariant::WENOJS5
            };
        }

        std::shared_ptr<GenericVariableAccessor> createVariableAccessor(
            ReconstructionMode mode,
            int numGhostCells = 2
        ) const override {
            return std::make_shared<CompressibleEulerVariableAccessor>(
                numGhostCells, mode);
        }
    };

    return std::make_shared<CompressibleEulerCapabilities>();
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
    ReconstructionMode mode, 
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

void CompressibleEulerVariableAccessor::TransformVariables(
    int direction,
    int lev, 
    void* solver_void,
    amrex::MultiFab& VariableBuffer,
    ReconstructionMode mode, 
    const SolverCapabilities::MethodValidationResult& validationResult
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
            // Copy conservative variables to the working buffer
            CopyConservativeVariablesToVariableBuffer(lev, solver, VariableBuffer);
            break;
            
        case ReconstructionMode::Characteristic:
            // transform to characteristic from primitive space
            CopyConservativeVariablesToVariableBuffer(lev, solver, VariableBuffer);
            ToCharacteristic(direction, lev, solver, VariableBuffer);
            break;
            
        default:
            Util::Abort(__FILE__, __func__, __LINE__, "Unsupported reconstruction mode");
    }
}

void CompressibleEulerVariableAccessor::TransformFluxes(
    int direction,
    int lev, 
    void* solver_void,
    amrex::MultiFab& CellFluxBuffer,
    ReconstructionMode mode, 
    const SolverCapabilities::MethodValidationResult& validationResult
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
            // transform to characteristic from primitive space
            CopyConservativeFluxesToCellFluxBuffer(direction, lev, solver, CellFluxBuffer);
            ToCharacteristic(direction, lev, solver, CellFluxBuffer);
            break;
            
        default:
            Util::Abort(__FILE__, __func__, __LINE__, "Unsupported reconstruction mode");
    }
}

void CompressibleEulerVariableAccessor::ReverseTransformFluxes(
    int direction,
    int lev,
    void* solver_void,
    amrex::MultiFab& SummedFlux,
    ReconstructionMode mode, 
    const SolverCapabilities::MethodValidationResult& validationResult
) const {
    // Cast the void pointer to ScimitarX pointer
    auto solver = static_cast<Integrator::ScimitarX*>(solver_void);
    
    // Handle different reconstruction modes
    switch(mode) {
        case ReconstructionMode::Primitive:
            // No transformation needed - the reconstructed values are already primitive
            break;
            
        case ReconstructionMode::Conservative:
            // Transform from Conservative to Primitive
            FromConservative(direction, lev, solver, SummedFlux);
            break;
            
        case ReconstructionMode::Characteristic:
            // Transform from characteristic variables back to primitive/conservative
            FromCharacteristic(direction, lev, solver, SummedFlux);
            break;
            
        default:
            Util::Abort(__FILE__, __func__, __LINE__, "Unsupported reconstruction mode");
    }
}

void CompressibleEulerVariableAccessor::CopyPrimitiveVariablesToVariableBuffer(
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
}


void CompressibleEulerVariableAccessor::CopyConservativeVariablesToVariableBuffer(
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

void CompressibleEulerVariableAccessor::CopyConservativeFluxesToCellFluxBuffer(
    int direction,    
    int lev,
    Integrator::ScimitarX* solver,
    amrex::MultiFab& CellFluxBuffer
) const {

    for (amrex::MFIter mfi(CellFluxBuffer, false); mfi.isValid(); ++mfi) {
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

}

void CompressibleEulerVariableAccessor::ToCharacteristic(
    int direction, 
    int lev, 
    const Integrator::ScimitarX* solver, 
    amrex::MultiFab& input_mf
) const {
    for (amrex::MFIter mfi(input_mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.growntilebox();
        auto const& p_arr = solver->PVec_mf.Patch(lev, mfi);
        auto const& W_arr  = input_mf.array(mfi);
        
        // Retrieve variable indices for better GPU performance
        int num_components = 5;
        int solver_components = solver->number_of_components;   

        int rho_idx = solver->variableIndex.DENS;   // density
        int ie_idx  = solver->variableIndex.IE;     // energy
        int u_idx   = solver->variableIndex.UVEL;   // x-momentum
        int v_idx   = solver->variableIndex.VVEL;   // y-momentum
#if AMREX_SPACEDIM == 3
        int w_idx   = solver->variableIndex.WVEL;   // z-momentum
#else
        int w_idx   = -1;                          // placeholder for 2D
#endif
        
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

            // Construct index arrays with correct number of components
            int index[3];
            int left_index[3];
            int right_index[3];
            int lower_bounds[3];
            int upper_bounds[3];
          
#if AMREX_SPACEDIM == 2
            index[0] = i; index[1] = j; index[2] = k;
            left_index[0] = i; left_index[1] = j; left_index[2] = k;
            right_index[0] = i; right_index[1] = j; right_index[2] = k;
          
            lower_bounds[0] = box.smallEnd(0); lower_bounds[1] = box.smallEnd(1); lower_bounds[2] = 0;
            upper_bounds[0] = box.bigEnd(0); upper_bounds[1] = box.bigEnd(1); upper_bounds[2] = 0;
#else
            index[0] = i; index[1] = j; index[2] = k;
            left_index[0] = i; left_index[1] = j; left_index[2] = k;
            right_index[0] = i; right_index[1] = j; right_index[2] = k;
          
            lower_bounds[0] = box.smallEnd(0); lower_bounds[1] = box.smallEnd(1); lower_bounds[2] = box.smallEnd(2);
            upper_bounds[0] = box.bigEnd(0); upper_bounds[1] = box.bigEnd(1); upper_bounds[2] = box.bigEnd(2);
#endif
      
            // Clamp the main index
            ClampIndices(index, lower_bounds, upper_bounds);
      
            // Shift and clamp for left and right indices
            int left_offset[3] = {0, 0, 0};
            int right_offset[3] = {0, 0, 0};
            left_offset[direction] = 0;  // Shift left
            right_offset[direction] = 1;  // Shift right
      
            ShiftAndClampIndices(left_index, left_offset, lower_bounds, upper_bounds, direction);
            ShiftAndClampIndices(right_index, right_offset, lower_bounds, upper_bounds, direction);
                
            // Extract original Current Index primitive variables
            Set::MultiVector U_orig = Numeric::FieldToMultiVector(W_arr, index[0], index[1], index[2], solver_components);
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

            // Compute right eigenvector matrix
            Set::MultiMatrix L_n = ComputeLeftEigenvectorMatrix(Wavg, direction, num_components);
                        
            // Create reordered conservative variables to match eigenvector matrix ordering
            Set::MultiVector U_reordered(num_components);
            U_reordered.setZero();
            
            // Reordering
            U_reordered(0) = U_orig(rho_idx);  // density 
            U_reordered(1) = U_orig(u_idx);  // x-momentum
            U_reordered(2) = U_orig(v_idx);  // y-momentum
#if AMREX_SPACEDIM == 3
            U_reordered(3) = U_orig(w_idx);  // z-momentum (3D)
#else
            U_reordered(3) = 0.0;                // z-momentum (2D)
#endif
            U_reordered(4) = U_orig(ie_idx);  // energy
                        
            Set::MultiVector W_vec_reordered = L_n * U_reordered;
                                    
            // Unreordering
            W_arr(index[0], index[1], index[2], rho_idx) = W_vec_reordered(0);  // density 
            W_arr(index[0], index[1], index[2], u_idx) = W_vec_reordered(1);  // x-momentum
            W_arr(index[0], index[1], index[2], v_idx) = W_vec_reordered(2);  // y-momentum
#if AMREX_SPACEDIM == 3
            W_arr(index[0], index[1], index[2], w_idx) = W_vec_reordered(3);  // z-momentum
#endif
            W_arr(index[0], index[1], index[2],ie_idx) = W_vec_reordered(4);  // energy
                        
        });
    }
}


void CompressibleEulerVariableAccessor::FromConservative(
    int direction, 
    int lev, 
    const Integrator::ScimitarX* solver,
    amrex::MultiFab& SummedFlux
) const {

    const int nghosts = solver->number_of_ghost_cells;
    const int num_components = solver->number_of_components; 

    for (amrex::MFIter mfi(SummedFlux, false); mfi.isValid(); ++mfi) {
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

void CompressibleEulerVariableAccessor::FromCharacteristic(
    int direction, 
    int lev, 
    const Integrator::ScimitarX* solver,
    amrex::MultiFab& SummedFlux
) const {
     
    const int nghosts = solver->number_of_ghost_cells;


    for (amrex::MFIter mfi(SummedFlux, false); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.grownnodaltilebox(direction, nghosts);
        auto const& TotalFlux_arr = SummedFlux.array(mfi);
        auto const& flux_arr = (direction == Directions::Xdir) ? solver->XFlux_mf.Patch(lev, mfi) :
                                 (direction == Directions::Ydir) ? solver->YFlux_mf.Patch(lev, mfi) :
                                 solver->ZFlux_mf.Patch(lev, mfi);
        auto const& p_arr = solver->PVec_mf.Patch(lev, mfi);
        //Characteristic Reconstruction is Performed through Fixed 5X5 Matrix
        //This works for both 2D and 3D. 
        //The matrix is re-arranged such that for 2D, last column, last row are zeros.
        int num_components = 5; 
        int solver_components = solver->number_of_components; 
        // Retrieve variable indices for better GPU performance
        int rho_idx = solver->variableIndex.DENS;   // density
        int ie_idx = solver->variableIndex.IE;     // energy
        int u_idx = solver->variableIndex.UVEL;   // x-momentum
        int v_idx = solver->variableIndex.VVEL;   // y-momentum
#if AMREX_SPACEDIM == 3
        int w_idx = solver->variableIndex.WVEL;   // z-momentum
#else
        int w_idx = -1;                          // placeholder for 2D
#endif
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {

    
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
            lower_bounds[0] = bx.smallEnd(0); 
            lower_bounds[1] = bx.smallEnd(1); 
            lower_bounds[2] = AMREX_SPACEDIM == 3 ? bx.smallEnd(2) : 0;
            upper_bounds[0] = bx.bigEnd(0)-1;
            upper_bounds[1] = bx.bigEnd(1)-1;
            upper_bounds[2] = AMREX_SPACEDIM == 3 ? bx.bigEnd(2) - 1 : 0;
                  
      
            // Clamp the main index
            ClampIndices(index, lower_bounds, upper_bounds);
      
            // Shift and clamp for left and right indices
            int left_offset[3] = {0, 0, 0};
            int right_offset[3] = {0, 0, 0};
            left_offset[direction] = 0;  // Shift left
            right_offset[direction] = 1;  // Shift right
      
            ShiftAndClampIndices(left_index, left_offset, lower_bounds, upper_bounds, direction);
            ShiftAndClampIndices(right_index, right_offset, lower_bounds, upper_bounds, direction);
    
            // Extract characteristic variables in original ordering
            Set::MultiVector W_orig = Numeric::FieldToMultiVector(TotalFlux_arr, index[0], index[1], index[2], solver_components);


            // Extract original Current Index primitive variables
            Set::MultiVector U_orig = Numeric::FieldToMultiVector(p_arr, index[0], index[1], index[2], solver_components);
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
                        
            // Compute right eigenvector matrices
            Set::MultiMatrix R_n = ComputeRightEigenvectorMatrix(Wavg, direction, num_components);

            // Create reordered characteristic variables
            Set::MultiVector W_reordered(num_components);
            W_reordered.setZero();

            // Reordering
            W_reordered(0) = W_orig(rho_idx);
            W_reordered(1) = W_orig(u_idx);
            W_reordered(2) = W_orig(v_idx);
#if AMREX_SPACEDIM == 3
            W_reordered(3) = W_orig(w_idx);
#else
            W_reordered(3) = 0.0;
#endif
            W_reordered(4) = W_orig(ie_idx);
            
                                    
            // Convert to conservative variables (in reordered space)
            Set::MultiVector U_reordered = R_n * W_reordered;
                        
            // Unreordering and storing back in arrays
            flux_arr(index[0], index[1], index[2], rho_idx) = U_reordered(0);
            flux_arr(index[0], index[1], index[2], u_idx)   = U_reordered(1);
            flux_arr(index[0], index[1], index[2], v_idx)   = U_reordered(2);
#if AMREX_SPACEDIM == 3
            flux_arr(index[0], index[1], index[2], w_idx)   = U_reordered(3);
#endif
            flux_arr(index[0], index[1], index[2], ie_idx)  = U_reordered(4);
            
        });
    }
}

Set::MultiMatrix 
CompressibleEulerVariableAccessor::ComputeRightEigenvectorMatrix(
    const Set::MultiVector& W, int dir, int num_components) {
    // Constants
    const Set::Scalar gamma = 1.4;
    const Set::Scalar g = gamma - 1.0;  // g = gamma - 1.0
    
    // Extract primitive variables from W    
    const Set::Scalar u = W(1);
    const Set::Scalar v = W(2);
    const Set::Scalar w = W(3) * (AMREX_SPACEDIM == 3 ? 1.0 : 0.0);
    const Set::Scalar H = W(4);  // Total enthalpy

    // Derived quantities
    const Set::Scalar Va = std::sqrt(u*u + v*v + w*w);  // Velocity magnitude
    const Set::Scalar a = std::sqrt(g * (H - 0.5*Va*Va));  // Sound speed
    
    // Initialize right eigenvector matrix
    Set::MultiMatrix R_n = Set::MultiMatrix::Zero(num_components, num_components);
    
    // Use switch for better GPU performance than if-else chain
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
        R_n(4, 1) = 0.5*Va*Va;
        
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
        R_n(4, 2) = 0.5*Va*Va;
        
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
        R_n(4, 3) = 0.5*Va*Va;
        
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
    const Set::MultiVector& W, int dir, int num_components) {
    // Constants
    const Set::Scalar gamma = 1.4;
    const Set::Scalar g = gamma - 1.0;  // g = gamma - 1.0
    
    // Extract primitive variables from W    
    const Set::Scalar u = W(1);
    const Set::Scalar v = W(2);
    const Set::Scalar w = W(3) * (AMREX_SPACEDIM == 3 ? 1.0 : 0.0);
    const Set::Scalar H = W(4);  // Total enthalpy

    // Derived quantities
    const Set::Scalar Va = std::sqrt(u*u + v*v + w*w);  // Velocity magnitude
    const Set::Scalar a = std::sqrt(g * (H - 0.5*Va*Va));  // Sound speed
    const Set::Scalar a2 = a * a;
    const Set::Scalar M = (0.5*g)/(a2);  // Scaling factor from Fortran code
    
    // Initialize left eigenvector matrix
    Set::MultiMatrix L_n = Set::MultiMatrix::Zero(num_components, num_components);
    
    // Use switch for better GPU performance than if-else chain
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
    
    return L_n;
}
    
Set::MultiVector
CompressibleEulerVariableAccessor::ComputeRoeAverages(
    const Set::MultiVector& WL, 
    const Set::MultiVector& WR,
    int num_components)
{
        
    const Set::Scalar gamma = 1.4;  // Ratio of specific heats (hardcoded for now)
    Set::Scalar Va_L = std::sqrt(WL(1)*WL(1) + WL(2)*WL(2) + WL(3)*WL(3));
    Set::Scalar p_L = (gamma - 1.0) * WL(0) * WL(4);
    Set::Scalar c_L = std::sqrt(gamma * p_L / WL(0));  // Sound speed
    Set::Scalar h_L = 0.5*Va_L*Va_L + c_L*c_L/(gamma - 1.0);  // Specific enthalpy
        
    Set::Scalar Va_R = std::sqrt(WR(1)*WR(1) + WR(2)*WR(2) + WR(3)*WR(3));
    Set::Scalar p_R = (gamma - 1.0) * WR(0) * WR(4);
    Set::Scalar c_R = std::sqrt(gamma * p_R / WR(0));  // Sound speed
    Set::Scalar h_R = 0.5*Va_R*Va_R + c_R * c_R / (gamma - 1.0);  // Specific enthalpy
    
    // Compute Roe weights
    Set::Scalar sqrt_rho_L = std::sqrt(WL(0));
    Set::Scalar sqrt_rho_R = std::sqrt(WR(0));
    Set::Scalar inv_sqrt_sum = 1.0 / (sqrt_rho_L + sqrt_rho_R);

    Set::MultiVector W_avg(num_components);
    W_avg.setZero();
    
    // Compute Roe-averaged values using weighting formula
    W_avg(0) = sqrt_rho_L * sqrt_rho_R;  // Roe-averaged density (geometric mean)
    W_avg(1) = (sqrt_rho_L * WL(1) + sqrt_rho_R * WR(1)) * inv_sqrt_sum;  // X-velocity
    W_avg(2) = (sqrt_rho_L * WL(2) + sqrt_rho_R * WR(2)) * inv_sqrt_sum;  // Y-velocity
    W_avg(3) = (sqrt_rho_L * WL(3) + sqrt_rho_R * WR(3)) * inv_sqrt_sum;  // Z-velocity
    
    // Compute Roe-averaged enthalpy and use it to derive sound speed and internal energy
    Set::Scalar h_avg = (sqrt_rho_L * h_L + sqrt_rho_R * h_R) * inv_sqrt_sum;
    
    // Store total enthalpy
    W_avg(4) = h_avg;
    
    return W_avg;
}

} // namespace CompressibleEuler

}  // namespace Numeric 

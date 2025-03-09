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
                    return MethodSupport::Unsupported(
                        {"HLLC not fully implemented"},
                        {"Use LocalLaxFriedrichs for now"}
                    );
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
            return 3;  // WENO5 requires 3 ghost cells
        default:
            return total_ghost_cells;
    }
}
    

// Add these implementations to IntegratorVariableAccessLayer.cpp

void CompressibleEulerVariableAccessor::TransformVariables(
    int lev, 
    void* solver_void,
    amrex::MultiFab& working_buffer,
    ReconstructionMode mode, 
    const SolverCapabilities::MethodValidationResult& validationResult,
    int direction
) const {
    // Cast the void pointer to ScimitarX pointer
    auto solver = static_cast<Integrator::ScimitarX*>(solver_void);
    
    // Handle different reconstruction modes
    switch(mode) {
        case ReconstructionMode::Primitive:
            // Copy primitive variables to the working buffer
            CopyPrimitiveVariablesToWorkingBuffer(lev, solver, working_buffer);
            break;
            
        case ReconstructionMode::Conservative:
            // Copy conservative variables to the working buffer
            CopyConservativeVariablesToWorkingBuffer(lev, solver, working_buffer);
            break;
            
        case ReconstructionMode::Characteristic:
            // First copy primitive variables, then transform to characteristic
            CopyConservativeVariablesToWorkingBuffer(lev, solver, working_buffer);
            ToCharacteristic(direction, lev, solver, working_buffer);
            break;
            
        default:
            Util::Abort(__FILE__, __func__, __LINE__, "Unsupported reconstruction mode");
    }
}

void CompressibleEulerVariableAccessor::ReverseTransform(
    int lev,
    void* solver_void,
    amrex::MultiFab& QL_stencil,
    amrex::MultiFab& QR_stencil, 
    ReconstructionMode mode, 
    const SolverCapabilities::MethodValidationResult& validationResult,
    int direction
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
            FromConservative(direction, lev, solver, QL_stencil, QR_stencil);
            break;
            
        case ReconstructionMode::Characteristic:
            // Transform from characteristic variables back to primitive/conservative
            FromCharacteristic(direction, lev, solver, QL_stencil, QR_stencil);
            break;
            
        default:
            Util::Abort(__FILE__, __func__, __LINE__, "Unsupported reconstruction mode");
    }
}

void CompressibleEulerVariableAccessor::CopyPrimitiveVariablesToWorkingBuffer(
    int lev,
    Integrator::ScimitarX* solver,
    amrex::MultiFab& working_buffer
) const {
    // Copy from PVec_mf to working buffer
    amrex::MultiFab::Copy(
        working_buffer, 
        *(solver->PVec_mf[lev]),
        0, 0, 
        solver->number_of_components, 
        working_buffer.nGrow()
    );
    working_buffer.FillBoundary(solver->geom[lev].periodicity());
}

void CompressibleEulerVariableAccessor::CopyConservativeVariablesToWorkingBuffer(
    int lev,
    Integrator::ScimitarX* solver,
    amrex::MultiFab& working_buffer
) const {
    // Copy from QVec_mf to working buffer
    amrex::MultiFab::Copy(
        working_buffer, 
        *(solver->QVec_mf[lev]),
        0, 0, 
        solver->number_of_components, 
        working_buffer.nGrow()
    );
    working_buffer.FillBoundary(solver->geom[lev].periodicity());
}

void CompressibleEulerVariableAccessor::ToCharacteristic(
    int direction, 
    int lev, 
    const Integrator::ScimitarX* solver, 
    amrex::MultiFab& input_mf
) const {
    for (amrex::MFIter mfi(input_mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.growntilebox();
        auto const& Q_arr = solver->QVec_mf.Patch(lev, mfi);
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

            // Extract original conservative variables
            Set::MultiVector U_orig = Numeric::FieldToMultiVector(Q_arr, i, j, k, solver_components);
            
            // Create primitive variables for eigenvector computation
            Set::MultiVector W(num_components);
            W.setZero();

            W(0) = U_orig(rho_idx);                         // density
            W(1) = U_orig(ie_idx) / U_orig(rho_idx);    // internal energy
            W(2) = U_orig(u_idx) / U_orig(rho_idx);    // x-velocity
            W(3) = U_orig(v_idx) / U_orig(rho_idx);    // y-velocity
#if AMREX_SPACEDIM == 3
            W(4) = U_orig(w_idx) / U_orig(rho_idx);    // z-velocity (3D)
#else
            W(4) = 0.0;                                        // z-velocity (2D)
#endif
            
            // Compute right eigenvector matrix
            Set::MultiMatrix R_n = ComputeRightEigenvectorMatrix(W, direction, num_components);
            
            // Create reordered conservative variables to match eigenvector matrix ordering
            Set::MultiVector U_reordered(num_components);
            U_reordered.setZero();
            
            // Reordering
            U_reordered(0) = U_orig(rho_idx);  // density 
            U_reordered(1) = U_orig(ie_idx);  // energy
            U_reordered(2) = U_orig(u_idx);  // x-momentum
            U_reordered(3) = U_orig(v_idx);  // y-momentum
#if AMREX_SPACEDIM == 3
            U_reordered(4) = U_orig(w_idx);  // z-momentum (3D)
#else
            U_reordered(4) = 0.0;                // z-momentum (2D)
#endif
            
            // Compute characteristic variables
            Set::MultiMatrix L_n = R_n.inverse();  // Left eigenvector matrix
            Set::MultiVector W_vec_reordered = L_n * U_reordered;
            
            // Create vector for original-ordered characteristic variables
            Set::MultiVector W_vec_original(solver_components);
            W_vec_original.setZero();
            
            // Unreordering
            W_vec_original(rho_idx) = W_vec_reordered(0);  // density 
            W_vec_original(ie_idx) = W_vec_reordered(1);  // energy
            W_vec_original(u_idx) = W_vec_reordered(2);  // x-momentum
            W_vec_original(v_idx) = W_vec_reordered(3);  // y-momentum
#if AMREX_SPACEDIM == 3
            W_vec_original(w_idx) = W_vec_reordered(4);  // z-momentum
#endif
            
            // Store transformed characteristic variables in original ordering
            Numeric::MultiVectorToField(W_arr, i, j, k, W_vec_original, solver_components);
        });
    }
}

void CompressibleEulerVariableAccessor::FromConservative(
    int direction, 
    int lev, 
    const Integrator::ScimitarX* solver,
    amrex::MultiFab& QL_stencil, 
    amrex::MultiFab& QR_stencil
) const {

    const int nghosts = solver->number_of_ghost_cells;

    for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
        const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);
        auto const& QL_arr = QL_stencil.array(mfi);
        auto const& QR_arr = QR_stencil.array(mfi);
        
        int num_components = solver->number_of_components;   
        
        int rho_idx = solver->variableIndex.DENS;   // density
        int ie_idx  = solver->variableIndex.IE;     // energy
        int u_idx   = solver->variableIndex.UVEL;   // x-momentum
        int v_idx   = solver->variableIndex.VVEL;   // y-momentum
#if AMREX_SPACEDIM == 3
        int w_idx = solver->variableIndex.WVEL;   // z-momentum
#endif


        amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Extract conservative variables in original ordering
            Set::MultiVector QCons_L = Numeric::FieldToMultiVector(QL_arr, i, j, k, num_components);
            Set::MultiVector QCons_R = Numeric::FieldToMultiVector(QR_arr, i, j, k, num_components);
            
            // Create vector for storing Primitive QR and QL variables
            Set::MultiVector QPrim_L(num_components);
            QPrim_L.setZero();
            Set::MultiVector QPrim_R(num_components);
            QPrim_R.setZero();
            
            // Conversion of Left States
            QPrim_L[rho_idx] = QCons_L(rho_idx);
            QPrim_L[u_idx]   = QCons_L(u_idx)/QCons_L(rho_idx);
            QPrim_L[v_idx]   = QCons_L(v_idx)/QCons_L(rho_idx);
#if AMREX_SPACEDIM == 3
            QPrim_L[w_idx]   = QCons_L(w_idx)/QCons_L(rho_idx);
#endif
            Set::Scalar KE_L = 0.5 * (QPrim_L[u_idx]*QPrim_L[u_idx] + QPrim_L[v_idx]*QPrim_L[v_idx]
#if AMREX_SPACEDIM == 3
                                    + QPrim_L[w_idx]*QPrim_L[w_idx]
#endif                    
                    );
            Set::Scalar TE_L =  QCons_L(ie_idx)/QCons_L(rho_idx);

            QPrim_L[ie_idx]   = TE_L - KE_L;
            
            // Conversion of Right States
            QPrim_R[rho_idx] = QCons_R(rho_idx);
            QPrim_R[u_idx]   = QCons_R(u_idx)/QCons_R(rho_idx);
            QPrim_R[v_idx]   = QCons_R(v_idx)/QCons_R(rho_idx);
#if AMREX_SPACEDIM == 3
            QPrim_R[w_idx]   = QCons_R(w_idx)/QCons_R(rho_idx);
#endif
            Set::Scalar KE_R = 0.5 * (QPrim_R[u_idx]*QPrim_R[u_idx] + QPrim_R[v_idx]*QPrim_R[v_idx]
#if AMREX_SPACEDIM == 3
                                    + QPrim_R[w_idx]*QPrim_R[w_idx]
#endif                    
                    );
            Set::Scalar TE_R =  QCons_R(ie_idx)/QCons_R(rho_idx);

            QPrim_R[ie_idx]   = TE_R - KE_R;

            Numeric::MultiVectorToField(QL_arr, i, j, k, QPrim_L, num_components); 
            Numeric::MultiVectorToField(QR_arr, i, j, k, QPrim_R, num_components); 

        });
    }
}


void CompressibleEulerVariableAccessor::FromCharacteristic(
    int direction, 
    int lev, 
    const Integrator::ScimitarX* solver,
    amrex::MultiFab& QL_stencil, 
    amrex::MultiFab& QR_stencil
) const {
     
    const int nghosts = solver->number_of_ghost_cells;

    for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
        const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);
        auto const& QL_arr = QL_stencil.array(mfi);
        auto const& QR_arr = QR_stencil.array(mfi);
        
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
        
        amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Extract characteristic variables in original ordering
            Set::MultiVector W_L_orig = Numeric::FieldToMultiVector(QL_arr, i, j, k, solver_components);
            Set::MultiVector W_R_orig = Numeric::FieldToMultiVector(QR_arr, i, j, k, solver_components);
            
            // Create reordered characteristic variables
            Set::MultiVector W_L_reordered(num_components);
            W_L_reordered.setZero();
            Set::MultiVector W_R_reordered(num_components);
            W_R_reordered.setZero();
            
            // Reordering
            W_L_reordered(0) = W_L_orig(rho_idx);
            W_L_reordered(1) = W_L_orig(ie_idx);
            W_L_reordered(2) = W_L_orig(u_idx);
            W_L_reordered(3) = W_L_orig(v_idx);
#if AMREX_SPACEDIM == 3
            W_L_reordered(4) = W_L_orig(w_idx);
#else
            W_L_reordered(4) = 0.0;
#endif
            
            W_R_reordered(0) = W_R_orig(rho_idx);
            W_R_reordered(1) = W_R_orig(ie_idx);
            W_R_reordered(2) = W_R_orig(u_idx);
            W_R_reordered(3) = W_R_orig(v_idx);
#if AMREX_SPACEDIM == 3
            W_R_reordered(4) = W_R_orig(w_idx);
#else
            W_R_reordered(4) = 0.0;
#endif
            
            // Convert to conservative variables (in reordered space)
            Set::MultiMatrix R_n_L = ComputeRightEigenvectorMatrix(W_L_reordered, direction, num_components);
            Set::MultiMatrix R_n_R = ComputeRightEigenvectorMatrix(W_R_reordered, direction, num_components);
            
            Set::MultiVector U_L_reordered = R_n_L * W_L_reordered;
            Set::MultiVector U_R_reordered = R_n_R * W_R_reordered;
            
            // Unreordering and storing back in arrays
            QL_arr(i, j, k, rho_idx) = U_L_reordered(0);
            QL_arr(i, j, k, ie_idx)  = U_L_reordered(1);
            QL_arr(i, j, k, u_idx)   = U_L_reordered(2);
            QL_arr(i, j, k, v_idx)   = U_L_reordered(3);
#if AMREX_SPACEDIM == 3
            QL_arr(i, j, k, w_idx)   = U_L_reordered(4);
#endif
            
            QR_arr(i, j, k, rho_idx) = U_R_reordered(0);
            QR_arr(i, j, k, ie_idx)  = U_R_reordered(1);
            QR_arr(i, j, k, u_idx)   = U_R_reordered(2);
            QR_arr(i, j, k, v_idx)   = U_R_reordered(3);
#if AMREX_SPACEDIM == 3
            QR_arr(i, j, k, w_idx)   = U_R_reordered(4);
#endif
        });
    }
}

Set::MultiMatrix 
CompressibleEulerVariableAccessor::ComputeRightEigenvectorMatrix(const Set::MultiVector& W, int dir, int num_components) {
    // Constants
    const Set::Scalar gamma = 1.4;
    
    // Extract primitive variables from W    
    const Set::Scalar rho = W(0);
    const Set::Scalar e = W(1);
    const Set::Scalar u = W(2);
    const Set::Scalar v = W(3);
    const Set::Scalar w = W(4) * (AMREX_SPACEDIM == 3 ? 1.0 : 0.0);

    // Derived quantities
    const Set::Scalar q2 = u*u + v*v + w*w;
    const Set::Scalar p = (gamma - 1.0) * rho * (e - 0.5*q2);
    const Set::Scalar c = std::sqrt(std::max(gamma * p / rho, 1e-10)); // Prevent division by zero
    const Set::Scalar H = e + p/rho;

  // Create orthonormal basis using arrays instead of Set::Vector
    Set::Scalar n[3] = {0.0, 0.0, 0.0};
    Set::Scalar m[3] = {0.0, 0.0, 0.0};
    Set::Scalar l[3] = {0.0, 0.0, 0.0};
    
    // Set normal vector components based on direction
    n[0] = (dir == 0) ? 1.0 : 0.0;
    n[1] = (dir == 1) ? 1.0 : 0.0;
    n[2] = (dir == 2) ? 1.0 : 0.0;
    
    // Circular permutation for the other two orthogonal vectors
    m[0] = (dir == 1) ? 1.0 : ((dir == 2) ? 1.0 : 0.0);
    m[1] = (dir == 2) ? 1.0 : ((dir == 0) ? 1.0 : 0.0);
    m[2] = (dir == 0) ? 1.0 : ((dir == 1) ? 1.0 : 0.0);
    
    l[0] = (dir == 2) ? 1.0 : ((dir == 1) ? 1.0 : 0.0);
    l[1] = (dir == 0) ? 1.0 : ((dir == 2) ? 1.0 : 0.0);
    l[2] = (dir == 1) ? 1.0 : ((dir == 0) ? 1.0 : 0.0);


    // Velocity projections
    const Set::Scalar qn = u*n[0] + v*n[1] + w*n[2];
    const Set::Scalar ql = u*l[0] + v*l[1] + w*l[2];
    const Set::Scalar qm = u*m[0] + v*m[1] + w*m[2];

    // Initialize right eigenvector matrix
    Set::MultiMatrix R_n = Set::MultiMatrix::Zero(num_components, num_components);

    // Row 1: Density (ρ)
    R_n.row(0) << 1.0, 1.0, 1.0, 0.0, 0.0;

    // Row 2: Energy (ρE)
    R_n.row(1) << H - c*qn, 0.5*q2, H + c*qn, ql, qm;

    // Rows 3-5: Momentum components
    R_n.row(2) << u - c*n[0], u, u + c*n[0], l[0], m[0];
    R_n.row(3) << v - c*n[1], v, v + c*n[1], l[1], m[1];
    R_n.row(4) << w - c*n[2], w, w + c*n[2], l[2], m[2];

    return R_n;
}

} // namespace CompressibleEuler

}  // namespace Numeric 

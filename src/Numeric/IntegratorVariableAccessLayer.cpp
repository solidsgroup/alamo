#include "Numeric/IntegratorVariableAccessLayer.H"
#include "Integrator/ScimitarX.H"
#include "Numeric/Stencil.H"

namespace Numeric {

void CompressibleEuler::CompressibleEulerVariableAccessor::CopyPrimitiveVariablesToWorkingBuffer(
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

void CompressibleEuler::CompressibleEulerVariableAccessor::CopyConservativeVariablesToWorkingBuffer(
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

void CompressibleEuler::CompressibleEulerVariableAccessor::ToCharacteristic(
    int direction, 
    int lev, 
    const Integrator::ScimitarX* solver, 
    amrex::MultiFab& input_mf
) const {
    for (amrex::MFIter mfi(input_mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        auto const& Q_arr = solver->QVec_mf.Patch(lev, mfi);
        auto const& W_arr = input_mf.array(mfi);
        
        // Retrieve variable indices for better GPU performance
        int indices[5];
        indices[0] = solver->variableIndex.DENS;   // density
        indices[1] = solver->variableIndex.IE;     // energy
        indices[2] = solver->variableIndex.UVEL;   // x-momentum
        indices[3] = solver->variableIndex.VVEL;   // y-momentum
#if AMREX_SPACEDIM == 3
        indices[4] = solver->variableIndex.WVEL;   // z-momentum
#else
        indices[4] = -1;                          // placeholder for 2D
#endif
        
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Extract original conservative variables
            Set::Vector U_orig = Numeric::FieldToVector(Q_arr, i, j, k);
            
            // Create primitive variables for eigenvector computation
            Set::Vector W;
            W(0) = U_orig(indices[0]);                         // density
            W(1) = U_orig(indices[1]) / U_orig(indices[0]);    // internal energy
            W(2) = U_orig(indices[2]) / U_orig(indices[0]);    // x-velocity
            W(3) = U_orig(indices[3]) / U_orig(indices[0]);    // y-velocity
#if AMREX_SPACEDIM == 3
            W(4) = U_orig(indices[4]) / U_orig(indices[0]);    // z-velocity (3D)
#else
            W(4) = 0.0;                                        // z-velocity (2D)
#endif
            
            // Compute right eigenvector matrix
            Set::Matrix R_n = ComputeRightEigenvectorMatrix(W, direction);
            
            // Create reordered conservative variables to match eigenvector matrix ordering
            Set::Vector U_reordered;
            
            // Reordering
            U_reordered(0) = U_orig(indices[0]);  // density 
            U_reordered(1) = U_orig(indices[1]);  // energy
            U_reordered(2) = U_orig(indices[2]);  // x-momentum
            U_reordered(3) = U_orig(indices[3]);  // y-momentum
#if AMREX_SPACEDIM == 3
            U_reordered(4) = U_orig(indices[4]);  // z-momentum (3D)
#else
            U_reordered(4) = 0.0;                // z-momentum (2D)
#endif
            
            // Compute characteristic variables
            Set::Matrix L_n = R_n.inverse();  // Left eigenvector matrix
            Set::Vector W_vec_reordered = L_n * U_reordered;
            
            // Create vector for original-ordered characteristic variables
            Set::Vector W_vec_original = Set::Vector::Zero();
            
            // Unreordering
            W_vec_original(indices[0]) = W_vec_reordered(0);  // density 
            W_vec_original(indices[1]) = W_vec_reordered(1);  // energy
            W_vec_original(indices[2]) = W_vec_reordered(2);  // x-momentum
            W_vec_original(indices[3]) = W_vec_reordered(3);  // y-momentum
#if AMREX_SPACEDIM == 3
            W_vec_original(indices[4]) = W_vec_reordered(4);  // z-momentum
#endif
            
            // Store transformed characteristic variables in original ordering
            Numeric::VectorToField(W_arr, i, j, k, W_vec_original);
        });
    }
}

void CompressibleEuler::CompressibleEulerVariableAccessor::FromCharacteristic(
    int direction, 
    int lev, 
    const Integrator::ScimitarX* solver,
    amrex::MultiFab& QL_stencil, 
    amrex::MultiFab& QR_stencil
) const {
    for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        auto const& QL_arr = QL_stencil.array(mfi);
        auto const& QR_arr = QR_stencil.array(mfi);
        
        // Retrieve variable indices for better GPU performance
        int indices[5];
        indices[0] = solver->variableIndex.DENS;   // density
        indices[1] = solver->variableIndex.IE;     // energy
        indices[2] = solver->variableIndex.UVEL;   // x-momentum
        indices[3] = solver->variableIndex.VVEL;   // y-momentum
#if AMREX_SPACEDIM == 3
        indices[4] = solver->variableIndex.WVEL;   // z-momentum
#else
        indices[4] = -1;                          // placeholder for 2D
#endif
        
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Extract characteristic variables in original ordering
            Set::Vector W_L_orig = Numeric::FieldToVector(QL_arr, i, j, k);
            Set::Vector W_R_orig = Numeric::FieldToVector(QR_arr, i, j, k);
            
            // Create reordered characteristic variables
            Set::Vector W_L_reordered, W_R_reordered;
            
            // Reordering
            W_L_reordered(0) = W_L_orig(indices[0]);
            W_L_reordered(1) = W_L_orig(indices[1]);
            W_L_reordered(2) = W_L_orig(indices[2]);
            W_L_reordered(3) = W_L_orig(indices[3]);
#if AMREX_SPACEDIM == 3
            W_L_reordered(4) = W_L_orig(indices[4]);
#else
            W_L_reordered(4) = 0.0;
#endif
            
            W_R_reordered(0) = W_R_orig(indices[0]);
            W_R_reordered(1) = W_R_orig(indices[1]);
            W_R_reordered(2) = W_R_orig(indices[2]);
            W_R_reordered(3) = W_R_orig(indices[3]);
#if AMREX_SPACEDIM == 3
            W_R_reordered(4) = W_R_orig(indices[4]);
#else
            W_R_reordered(4) = 0.0;
#endif
            
            // Convert to conservative variables (in reordered space)
            Set::Matrix R_n_L = ComputeRightEigenvectorMatrix(W_L_reordered, direction);
            Set::Matrix R_n_R = ComputeRightEigenvectorMatrix(W_R_reordered, direction);
            
            Set::Vector U_L_reordered = R_n_L * W_L_reordered;
            Set::Vector U_R_reordered = R_n_R * W_R_reordered;
            
            // Unreordering and storing back in arrays
            QL_arr(i, j, k, indices[0]) = U_L_reordered(0);
            QL_arr(i, j, k, indices[1]) = U_L_reordered(1);
            QL_arr(i, j, k, indices[2]) = U_L_reordered(2);
            QL_arr(i, j, k, indices[3]) = U_L_reordered(3);
#if AMREX_SPACEDIM == 3
            QL_arr(i, j, k, indices[4]) = U_L_reordered(4);
#endif
            
            QR_arr(i, j, k, indices[0]) = U_R_reordered(0);
            QR_arr(i, j, k, indices[1]) = U_R_reordered(1);
            QR_arr(i, j, k, indices[2]) = U_R_reordered(2);
            QR_arr(i, j, k, indices[3]) = U_R_reordered(3);
#if AMREX_SPACEDIM == 3
            QR_arr(i, j, k, indices[4]) = U_R_reordered(4);
#endif
        });
    }
}

Set::Matrix 
CompressibleEuler::CompressibleEulerVariableAccessor::ComputeRightEigenvectorMatrix(const Set::Vector& W, int dir) {
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
    Set::Matrix R_n = Set::Matrix::Zero();

    // Row 1: Density (ρ)
    R_n.row(0) << 1.0, 1.0, 0.0, 0.0, 1.0;

    // Row 2: Energy (ρE)
    R_n.row(1) << H - c*qn, 0.5*q2, ql, qm, H + c*qn;

    // Rows 3-5: Momentum components
    R_n.row(2) << u - c*n[0], u, l[0], m[0], u + c*n[0];
    R_n.row(3) << v - c*n[1], v, l[1], m[1], v + c*n[1];
    R_n.row(4) << w - c*n[2], w, l[2], m[2], w + c*n[2];

    return R_n;
}

}  // namespace Numeric 

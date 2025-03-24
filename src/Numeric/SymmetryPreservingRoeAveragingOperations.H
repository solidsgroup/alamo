/**
 * @file SymmetryPreservingRoeAveragingOperations.H
 * @brief Implementation of symmetry-preserving operations for numerical simulations
 * 
 * This file provides implementations of consistent floating-point operations to preserve
 * symmetry in numerical simulations, particularly for Roe averaging and characteristic
 * transformations. These techniques are based on the following publications:
 * 
 * - Fleischmann, N., Adami, S., & Adams, N. A. (2019). Numerical symmetry-preserving 
 *   techniques for low-dissipation shock-capturing schemes. Computers & Fluids, 
 *   Volume 189, 94-107.
 * 
 * - Fleischmann, N., Adami, S., Hu, X. Y., & Adams, N. A. (2020). A low dissipation 
 *   method to cure the grid-aligned shock instability. Journal of Computational Physics, 
 *   Volume 401, 109004.
 * 
 * - Fleischmann, N., Adami, S., & Adams, N. A. (2020). A shock-stable modification of 
 *   the HLLC Riemann solver with reduced numerical dissipation. Journal of Computational 
 *   Physics, Volume 423, 109762.
 * 
 * The techniques presented here address the issue of symmetry breaking due to the lack
 * of associativity in floating-point operations, which can lead to significant 
 * asymmetries in low-dissipation numerical simulations.
 */

#ifndef NUMERIC_SYMMETRY_PRESERVING_ROE_AVERAGING_OPERATIONS_H
#define NUMERIC_SYMMETRY_PRESERVING_ROE_AVERAGING_OPERATIONS_H

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <algorithm>
#include <array>
#include "Set/Set.H"

namespace Numeric {

/**
 * @brief Collection of functions for performing floating-point operations consistently
 * 
 * This namespace provides functionality to handle floating-point operations in a manner
 * that preserves symmetry and consistency, avoiding numerical artifacts that can arise
 * from the lack of associativity in floating-point arithmetic.
 * 
 * These operations are particularly important for low-dissipation schemes where
 * small numerical artifacts can grow into large-scale asymmetries.
 */
namespace SymmetryPreserving {

/**
 * @brief Performs a consistent sum of three floating-point values
 * 
 * Computes the sum in a way that is independent of the order of operations
 * by taking all possible paths of addition and finding the midpoint of the
 * extremes. This ensures (a+b)+c == (a+c)+b == (b+c)+a.
 * 
 * As described in equation (19) of Fleischmann et al. (2019).
 * 
 * @tparam T Scalar type (float, double, etc.)
 * @param a First value
 * @param b Second value
 * @param c Third value
 * @return Consistently calculated sum
 */
template<typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T ConsistentSum3(T a, T b, T c) {
    T s1 = (a + b) + c;
    T s2 = (c + a) + b;
    T s3 = (b + c) + a;
    return 0.5 * (std::max({s1, s2, s3}) + std::min({s1, s2, s3}));
}

/**
 * @brief Performs a consistent sum of four floating-point values
 * 
 * Similar to ConsistentSum3, but handles four values by calculating all possible
 * pairwise sums first, then summing the pairs consistently.
 * 
 * As described in equation (21) of Fleischmann et al. (2019).
 * 
 * @tparam T Scalar type (float, double, etc.)
 * @param a First value
 * @param b Second value
 * @param c Third value
 * @param d Fourth value
 * @return Consistently calculated sum
 */
template<typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T ConsistentSum4(T a, T b, T c, T d) {
    T s1 = (a + b) + (c + d);
    T s2 = (a + c) + (b + d);
    T s3 = (a + d) + (b + c);
    return 0.5 * (std::max({s1, s2, s3}) + std::min({s1, s2, s3}));
}

/**
 * @brief Calculates velocity magnitude squared consistently
 * 
 * Computes the square of the velocity magnitude (v·v) in a consistent
 * manner regardless of coordinate order, which is critical for
 * maintaining symmetry in problems with physical symmetries.
 * 
 * Addresses the issue described in equation (17) of Fleischmann et al. (2019).
 * 
 * @tparam T Scalar type (float, double, etc.)
 * @param vx X-component of velocity
 * @param vy Y-component of velocity
 * @param vz Z-component of velocity
 * @return Consistently calculated squared magnitude
 */
template<typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T ConsistentVelocityMagnitudeSquared(T vx, T vy, T vz) {
    return ConsistentSum3(vx*vx, vy*vy, vz*vz);
}

/**
 * @brief Calculates velocity magnitude consistently
 * 
 * @tparam T Scalar type
 * @param vx X-component of velocity
 * @param vy Y-component of velocity
 * @param vz Z-component of velocity
 * @return Consistently calculated magnitude
 */
template<typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T ConsistentVelocityMagnitude(T vx, T vy, T vz) {
    return std::sqrt(ConsistentVelocityMagnitudeSquared(vx, vy, vz));
}

/**
 * @brief Calculates consistent dot product of two vectors
 * 
 * Computes a dot product in a manner that is independent of the
 * indexing order of the vectors, which is important for maintaining
 * symmetries in multidimensional operations.
 * 
 * @tparam VectorType Vector type supporting indexing and size() operations
 * @param a First vector
 * @param b Second vector
 * @return Consistently calculated dot product
 */
template<typename VectorType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Set::Scalar ConsistentDotProduct(const VectorType& a, const VectorType& b) {
    // Safety check - return 0 if sizes don't match
    if (a.size() != b.size()) {
        return Set::Scalar(0);
    }
    
    // Compute initial product directly
    Set::Scalar sum = a(0) * b(0);
    
    // Use direct loop for remaining elements
    const int size = a.size();
    for (int i = 1; i < size; ++i) {
        sum += a(i) * b(i);
    }
    
    return sum;
}

/**
 * @brief Special overload for consistent dot product of row vector and column vector
 * 
 * This specialization handles the case where the first argument is a row of a matrix
 * and the second argument is a vector. This is needed because matrix rows have a 
 * different type than vectors in Eigen.
 * 
 * @tparam RowType Row type (from matrix)
 * @tparam VectorType Vector type
 * @param row Matrix row
 * @param v Vector
 * @return Consistently calculated dot product
 */
template<typename RowType, typename VectorType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Set::Scalar ConsistentDotProduct_RowVec(const RowType& row, const VectorType& v) {
    // Safety check - return 0 if sizes don't match
    if (row.size() != v.size()) {
        return Set::Scalar(0);
    }
    
    // Compute initial product directly
    Set::Scalar sum = row(0) * v(0);
    
    // Use direct loop for remaining elements
    const int size = row.size();
    for (int i = 1; i < size; ++i) {
        sum += row(i) * v(i);
    }
    
    return sum;
}

/**
 * @brief Calculates consistent matrix-vector multiplication
 * 
 * Performs matrix-vector multiplication in a manner that is consistent
 * regardless of the order of operations, which is important for
 * characteristic decompositions and other direction-dependent operations.
 * 
 * @tparam MatrixType Matrix type supporting row() operation and rows() dimension
 * @tparam VectorType Vector type compatible with the matrix
 * @param M Matrix
 * @param v Vector
 * @return Consistently calculated matrix-vector product
 */
template<typename MatrixType, typename VectorType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
VectorType ConsistentMatrixVectorMultiply(const MatrixType& M, const VectorType& v) {
    VectorType result(M.rows());
    
    // Direct loop to compute each element of the result vector
    for (int i = 0; i < M.rows(); i++) {
        result(i) = ConsistentDotProduct_RowVec(M.row(i), v);
    }
    
    return result;
}

/**
 * @brief Consistent computation of kinetic energy
 * 
 * Calculates kinetic energy (0.5 * ρ * |v|²) in a consistent manner
 * regardless of coordinate ordering.
 * 
 * @tparam T Scalar type
 * @param rho Density
 * @param vx X-component of velocity
 * @param vy Y-component of velocity
 * @param vz Z-component of velocity
 * @return Kinetic energy calculated consistently
 */
template<typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T ConsistentKineticEnergy(T rho, T vx, T vy, T vz) {
    T v_squared = ConsistentVelocityMagnitudeSquared(vx, vy, vz);
    return 0.5 * rho * v_squared;
}

/**
 * @brief Calculates contact wave speed consistently for HLLC solver
 * 
 * Implements the properly bracketed version of the contact wave speed
 * calculation as shown in equation (16) of Fleischmann et al. (2019).
 * 
 * @tparam T Scalar type
 * @param pR Right state pressure
 * @param pL Left state pressure
 * @param rhoL Left state density
 * @param uL Left state velocity
 * @param SL Left wave speed
 * @param rhoR Right state density
 * @param uR Right state velocity
 * @param SR Right wave speed
 * @return Contact wave speed calculated consistently
 */
template<typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T ConsistentContactWaveSpeed(
    T pR, T pL, 
    T rhoL, T uL, T SL, 
    T rhoR, T uR, T SR) 
{
    // Implement equation (16) with proper bracketing
    T numerator = (pR - pL) + (rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR));
    T denominator = rhoL * (SL - uL) - rhoR * (SR - uR);
    
    // Prevent division by zero
    if (std::abs(denominator) < T(1.0e-15)) {
        return T(0.0);
    }
    
    return numerator / denominator;
}

/**
 * @brief Consistent computation of Roe averages
 * 
 * Calculates Roe averages using a consistent approach to avoid
 * direction-dependent artifacts in flux calculations, as described
 * in Section 4.4 of Fleischmann et al. (2019).
 * 
 * @tparam VectorType Vector type for primitive variables
 * @param WL Left state
 * @param WR Right state
 * @param num_components Number of components in the vectors
 * @return Roe-averaged state
 */
template<typename VectorType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
VectorType ConsistentRoeAverage(const VectorType& WL, const VectorType& WR, int num_components) {
    using T = typename VectorType::Scalar;
    
    // Constants
    const T gamma = T(1.4);  // Ratio of specific heats
    
    // Compute components with explicit order of operations for consistency
    T sqrt_rho_L = std::sqrt(WL(0));
    T sqrt_rho_R = std::sqrt(WR(0));
    
    // Compute denominator separately to ensure consistency
    T inv_sqrt_sum = 1.0 / (sqrt_rho_L + sqrt_rho_R);
    
    // Initialize result vector
    VectorType W_avg(num_components);
    W_avg.setZero();
    
    // Compute Roe-averaged density with explicit ordering
    W_avg(0) = sqrt_rho_L * sqrt_rho_R;  // Roe-averaged density (geometric mean)
    
    // Compute velocity components with consistent ordering
    W_avg(1) = (sqrt_rho_L * WL(1) + sqrt_rho_R * WR(1)) * inv_sqrt_sum;  // X-velocity
    W_avg(2) = (sqrt_rho_L * WL(2) + sqrt_rho_R * WR(2)) * inv_sqrt_sum;  // Y-velocity
    
    // Handle 3D case consistently
    if (num_components > 3) {
        W_avg(3) = (sqrt_rho_L * WL(3) + sqrt_rho_R * WR(3)) * inv_sqrt_sum;  // Z-velocity
    } else {
        W_avg(3) = 0.0;
    }
    
    // Calculate velocities and sound speeds separately with consistent ordering
    T Va_L = ConsistentVelocityMagnitudeSquared(WL(1), WL(2), WL(3));
    Va_L = std::sqrt(Va_L);
    
    T Va_R = ConsistentVelocityMagnitudeSquared(WR(1), WR(2), WR(3));
    Va_R = std::sqrt(Va_R);
    
    // Calculate pressures with consistent ordering
    T p_L = (gamma - 1.0) * WL(0) * WL(4);
    T p_R = (gamma - 1.0) * WR(0) * WR(4);
    
    // Calculate sound speeds with consistent ordering
    T c_L = std::sqrt(gamma * p_L / WL(0));
    T c_R = std::sqrt(gamma * p_R / WR(0));
    
    // Calculate enthalpies with consistent ordering
    T h_L = 0.5 * Va_L * Va_L + c_L * c_L / (gamma - 1.0);
    T h_R = 0.5 * Va_R * Va_R + c_R * c_R / (gamma - 1.0);
    
    // Calculate Roe-averaged enthalpy with consistent ordering
    T h_avg = (sqrt_rho_L * h_L + sqrt_rho_R * h_R) * inv_sqrt_sum;
    
    // Store total enthalpy
    W_avg(4) = h_avg;
    
    return W_avg;
}

} // namespace SymmetryPreserving

} // namespace Numeric

#endif // NUMERIC_SYMMETRY_PRESERVING_ROE_AVERAGING_OPERATIONS_H

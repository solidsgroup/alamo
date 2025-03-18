/**
 * @file WENOReconstruction.H
 * @brief Implementation of Weighted Essentially Non-Oscillatory (WENO) schemes
 * 
 * This file provides a template-based, traits-driven implementation of WENO
 * reconstruction schemes for hyperbolic conservation laws. It supports multiple
 * WENO variants (JS, Z) and orders (3rd, 5th) through specialized traits classes.
 * 
 * References:
 * - Jiang, G. S., & Shu, C. W. (1996). Efficient implementation of weighted ENO schemes.
 *   Journal of computational physics, 126(1), 202-228.
 * - Borges, R., Carmona, M., Costa, B., & Don, W. S. (2008). An improved weighted 
 *   essentially non-oscillatory scheme for hyperbolic conservation laws. 
 *   Journal of Computational Physics, 227(6), 3191-3211.
 */

#ifndef NUMERIC_WENO_RECONSTRUCTION_H
#define NUMERIC_WENO_RECONSTRUCTION_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <vector>
#include <map>
#include <memory>
#include <array>
#include "Set/Set.H"
#include "Integrator/Integrator.H"
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"
#include "Numeric/NumericTypes.H"
#include "Numeric/IntegratorVariableAccessLayer.H"
#include "Numeric/FluxHandler.H"
#include "Numeric/NumericFactory.H"

namespace Numeric {

// Forward declarations
template<typename VariantTraits> class WENOReconstructionPolicy;

/**
 * @class DifferenceOperator
 * @brief Base trait for difference operators used in smoothness indicators
 * 
 * This structure defines the coefficients for computing finite differences
 * used in WENO smoothness indicators, represented as matrices.
 * 
 * @tparam NumDifferenceTerms Number of difference terms in smoothness indicator
 * @tparam PointsPerStencil Number of points in each candidate stencil
 * @tparam NumCandidateStencils Number of candidate stencils
 */
template<int NumDifferenceTerms, int PointsPerStencil, int NumCandidateStencils>
struct DifferenceOperator {
    // Term weights for the smoothness indicator calculation
    static constexpr std::array<Set::Scalar, NumDifferenceTerms> TermWeights = {};
    
    // Function to get coefficient matrices
    static std::vector<Set::MultiMatrix> GetCoefficients() {
        // Default implementation returns empty vector
        return std::vector<Set::MultiMatrix>(NumCandidateStencils);
    }
};

//==============================================================================
// WENO-JS3 Traits - 3rd order Jiang-Shu WENO scheme
//==============================================================================

/**
 * @class WENOJS3Traits
 * @brief Traits for 3rd order WENO-JS scheme
 * 
 * Implements the 3rd order WENO scheme from Jiang & Shu (1996).
 * This scheme uses 2 candidate stencils with 2 points each.
 */
struct WENOJS3Traits {
    static constexpr int OrderValue = 3;
    static constexpr Set::Scalar Epsilon = 1.0e-6;  // Small parameter to avoid division by zero
    static constexpr int NumCandidateStencils = 2;  // Number of candidate stencils
    static constexpr int PointsPerStencil = 2;      // Points per candidate stencil
    static constexpr int TotalPoints = 3;           // Total points in the full stencil
    static constexpr int NumDifferenceTerms = 1;    // Number of difference terms in smoothness indicator
    static constexpr int Power = 2;                 // Power parameter in weights formula
    
    // Optimal (linear) weights for WENO-JS3 as a MultiVector
    static Set::MultiVector GetOptimalWeights() {
        Set::MultiVector weights(NumCandidateStencils);
        weights << 1.0/3.0, 2.0/3.0;
        return weights;
    }
    
    // Coefficients for candidate stencils reconstruction as a MultiMatrix
    static Set::MultiMatrix GetCandidateStencilCoefficients() {
        Set::MultiMatrix coeffs(NumCandidateStencils, PointsPerStencil);
        coeffs << -1.0/2.0, 3.0/2.0,   // Left stencil S₀ coefficients
                  1.0/2.0, 1.0/2.0;    // Right stencil S₁ coefficients
        return coeffs;
    }
    
    // Stencil offsets define which points to use for each candidate stencil
    static Set::MultiMatrix GetSubStencilOffsets() {
        Set::MultiMatrix offsets(NumCandidateStencils, PointsPerStencil);
        offsets << -1, 0,   // Left stencil S₀ uses points (i-1, i)
                   0, 1;    // Right stencil S₁ uses points (i, i+1)
        return offsets;
    }
    
    /**
     * @struct WENOJS3DifferenceOperator
     * @brief Difference operators for WENO-JS3 smoothness indicators
     */
    struct WENOJS3DifferenceOperator : 
        public DifferenceOperator<NumDifferenceTerms, PointsPerStencil, NumCandidateStencils> {
        // Define term weights for the smoothness indicator calculation
        static constexpr std::array<Set::Scalar, NumDifferenceTerms> TermWeights = {
            1.0  // Weight for the first difference term
        };
        
        // Function to get coefficient matrices
        static std::vector<Set::MultiMatrix> GetCoefficients() {
            std::vector<Set::MultiMatrix> coeffs(NumCandidateStencils);
            
            // Stencil 0 coefficients
            coeffs[0] = (Set::MultiMatrix(1, 2) << 
                1.0, -1.0    // 1st diff for stencil 0: f[i-1] - f[i]
            ).finished();
            
            // Stencil 1 coefficients
            coeffs[1] = (Set::MultiMatrix(1, 2) << 
                1.0, -1.0    // 1st diff for stencil 1: f[i] - f[i+1]
            ).finished();
            
            return coeffs;
        }
    };
    
    using DifferenceOperators = WENOJS3DifferenceOperator;
};

//==============================================================================
// WENO-JS5 Traits - 5th order Jiang-Shu WENO scheme
//==============================================================================

/**
 * @class WENOJS5Traits 
 * @brief Traits for 5th order WENO-JS scheme
 * 
 * Implements the 5th order WENO scheme from Jiang & Shu (1996).
 * This scheme uses 3 candidate stencils with 3 points each.
 */
struct WENOJS5Traits {
    static constexpr int OrderValue = 5;
    static constexpr Set::Scalar Epsilon = 1.0e-40;  // Small parameter to avoid division by zero
    static constexpr int NumCandidateStencils = 3;  // Number of candidate stencils
    static constexpr int PointsPerStencil = 3;      // Points per candidate stencil
    static constexpr int TotalPoints = 5;           // Total points in the full stencil
    static constexpr int NumDifferenceTerms = 2;    // Number of difference terms in smoothness indicator
    static constexpr int Power = 2;                 // Power parameter in weights formula
    
    // Optimal (linear) weights for WENO-JS5 as a MultiVector
    static Set::MultiVector GetOptimalWeights() {
        Set::MultiVector weights(NumCandidateStencils);
        weights << 0.1, 0.6, 0.3;
        return weights;
    }
    
    // Coefficients for candidate stencils reconstruction as a MultiMatrix
    static Set::MultiMatrix GetCandidateStencilCoefficients() {
        Set::MultiMatrix coeffs(NumCandidateStencils, PointsPerStencil);
        coeffs << 1.0/3.0, -7.0/6.0, 11.0/6.0,   // Left stencil S₀ coefficients
                 -1.0/6.0, 5.0/6.0, 1.0/3.0,     // Center stencil S₁ coefficients
                  1.0/3.0, 5.0/6.0, -1.0/6.0;    // Right stencil S₂ coefficients
        return coeffs;
    }
    
    // Stencil offsets define which points to use for each candidate stencil
    static Set::MultiMatrix GetSubStencilOffsets() {
        Set::MultiMatrix offsets(NumCandidateStencils, PointsPerStencil);
        offsets << -2, -1, 0,   // Left stencil S₀ uses points (i-2, i-1, i)
                  -1, 0, 1,     // Center stencil S₁ uses points (i-1, i, i+1)
                   0, 1, 2;     // Right stencil S₂ uses points (i, i+1, i+2)
        return offsets;
    }
    
    /**
     * @struct WENOJS5DifferenceOperator
     * @brief Difference operators for WENO-JS5 smoothness indicators
     */
    struct WENOJS5DifferenceOperator : 
        public DifferenceOperator<NumDifferenceTerms, PointsPerStencil, NumCandidateStencils> {
        
        // Define term weights for the smoothness indicator calculation
        static constexpr std::array<Set::Scalar, NumDifferenceTerms> TermWeights = {
            13.0/12.0,  // Weight for the first difference term
            1.0/4.0     // Weight for the second difference term
        };
        
        // Function to get coefficient matrices
        static std::vector<Set::MultiMatrix> GetCoefficients() {
            std::vector<Set::MultiMatrix> coeffs(NumCandidateStencils);
            
            // Stencil 0 coefficients
            coeffs[0] = (Set::MultiMatrix(2, 3) << 
                1.0, -2.0, 1.0,   // 1st diff for stencil 0: f[i-2] - 2*f[i-1] + f[i]
                1.0, -4.0, 3.0    // 2nd diff for stencil 0: f[i-2] - 4*f[i-1] + 3*f[i]
            ).finished();
            
            // Stencil 1 coefficients
            coeffs[1] = (Set::MultiMatrix(2, 3) << 
                1.0, -2.0, 1.0,   // 1st diff for stencil 1: f[i-1] - 2*f[i] + f[i+1]
                1.0, 0.0, -1.0    // 2nd diff for stencil 1: f[i-1] - f[i+1]
            ).finished();
            
            // Stencil 2 coefficients
            coeffs[2] = (Set::MultiMatrix(2, 3) << 
                1.0, -2.0, 1.0,   // 1st diff for stencil 2: f[i] - 2*f[i+1] + f[i+2]
                3.0, -4.0, 1.0    // 2nd diff for stencil 2: 3*f[i] - 4*f[i+1] + f[i+2]
            ).finished();
            
            return coeffs;
        }
    };
    
    using DifferenceOperators = WENOJS5DifferenceOperator;
};

//==============================================================================
// WENO-Z5 Traits - 5th order WENO-Z scheme
//==============================================================================

/**
 * @class WENOZ5Traits
 * @brief Traits for 5th order WENO-Z scheme from Borges et al. (2008)
 * 
 * Extends the 5th order WENO-JS scheme with a global smoothness measure
 * to improve accuracy near critical points.
 */
struct WENOZ5Traits : public WENOJS5Traits {
    // Inherits everything from WENOJS5Traits
    // Weight calculation is overridden in the WENOReconstructionPolicy
};

//==============================================================================
// Generic WENO Reconstruction Policy
//==============================================================================

/**
 * @class WENOReconstructionPolicy
 * @brief Implementation of WENO reconstruction algorithm
 * 
 * This template class implements the core WENO reconstruction algorithm, 
 * using the traits class to customize behavior for different variants.
 * 
 * @tparam VariantTraits Traits class defining WENO variant parameters
 */
template<typename VariantTraits>
class WENOReconstructionPolicy {
public:
    using Traits = VariantTraits;
    
    /**
     * Compute smoothness indicators for a given stencil
     * 
     * @param stencil Array of stencil values
     * @param i, j, k Cell indices for debugging
     * @param enableDebug Enable detailed debugging output
     * @return Array of smoothness indicators for each candidate stencil
     */
    template <typename StencilArray>
    AMREX_GPU_DEVICE
    static auto ComputeSmoothnessIndicators(
        const StencilArray& stencil, 
        int i, int j, int k,  // Cell indices for debug output
        bool enableDebug = false
    ) {
        std::array<Set::Scalar, Traits::NumCandidateStencils> indicators;
        
        // Get difference operator coefficients
        auto stencilCoefficients = Traits::DifferenceOperators::GetCoefficients();
        auto termWeights = Traits::DifferenceOperators::TermWeights;
        auto subStencilOffsets = Traits::GetSubStencilOffsets();
        
        // Loop over each candidate stencil
        for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
            indicators[r] = 0.0;
            
            // Extract stencil values for this candidate stencil
            Set::MultiVector stencil_values(Traits::PointsPerStencil);
            for (int p = 0; p < Traits::PointsPerStencil; ++p) {
                int stencil_idx = subStencilOffsets(r, p) + Traits::TotalPoints/2;
                stencil_values(p) = stencil[stencil_idx];
            }
            
            // Loop over each difference term in the smoothness indicator
            for (int term = 0; term < Traits::NumDifferenceTerms; ++term) {
                // Apply the difference operator for this term using dot product
                Set::Scalar diff_value = 
                    stencilCoefficients[r].row(term).dot(stencil_values);
                
                // Add the weighted squared difference to the indicator
                indicators[r] += termWeights[term] * diff_value * diff_value;
            }
        }
        
        // Debug call - separated from the main algorithm logic
        Util::ScimitarX_Util::Debug::DebugWENOSmoothnessIndicators(
            i, j, k, stencil, indicators, 
            stencilCoefficients,
            termWeights,
            subStencilOffsets, 
            enableDebug
        );
        
        return indicators;
    }
    
    /**
     * Compute candidate fluxes for each stencil
     * 
     * @param stencil Array of stencil values
     * @param i, j, k Cell indices for debugging
     * @param enableDebug Enable detailed debugging output
     * @return Array of candidate reconstructed values for each stencil
     */
    template <typename StencilArray>
    AMREX_GPU_DEVICE
    static auto ComputeCandidateFluxes(
        const StencilArray& stencil,
        int i, int j, int k,  // Cell indices for debug output
        bool enableDebug = false
    ) {
        std::array<Set::Scalar, Traits::NumCandidateStencils> fluxes;
        
        // Get coefficients and offsets
        auto candidateStencilCoefficients = Traits::GetCandidateStencilCoefficients();
        auto subStencilOffsets = Traits::GetSubStencilOffsets();
        
        // Compute reconstruction for each candidate stencil
        for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
            // Create vector for stencil values
            Set::MultiVector stencil_values(Traits::PointsPerStencil);
            
            // Fill vector with stencil values
            for (int p = 0; p < Traits::PointsPerStencil; ++p) {
                int stencil_idx = subStencilOffsets(r, p) + Traits::TotalPoints/2;
                stencil_values(p) = stencil[stencil_idx];
            }
            
            // Extract coefficients for this stencil
            Set::MultiVector coeffs = candidateStencilCoefficients.row(r);
            
            // Compute flux using dot product
            fluxes[r] = coeffs.dot(stencil_values);
        }
        
        // Debug call for candidate fluxes
        Util::ScimitarX_Util::Debug::DebugWENOCandidateFluxes(
            i, j, k, stencil, fluxes,
            candidateStencilCoefficients,
            subStencilOffsets,
            enableDebug
        );
        
        return fluxes;
    }
    
    /**
     * Compute non-linear weights for WENO-JS scheme
     * 
     * @param indicators Array of smoothness indicators
     * @param i, j, k Cell indices for debugging
     * @param enableDebug Enable detailed debugging output
     * @return Array of normalized non-linear weights
     */
    template <typename IndicatorArray>
    AMREX_GPU_DEVICE
    static auto ComputeNonlinearWeights(
        const IndicatorArray& indicators,
        int i, int j, int k,  // Cell indices for debug output
        bool enableDebug = false
    ) {
        std::array<Set::Scalar, Traits::NumCandidateStencils> weights;
        
        // Get optimal weights
        auto optimalWeights = Traits::GetOptimalWeights();
        
        // Calculate alpha values and their sum
        Set::Scalar alpha_sum = 0.0;
        std::array<Set::Scalar, Traits::NumCandidateStencils> alphas;
        
        for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
            // Safe indicator value to prevent division by zero
            Set::Scalar safe_indicator = std::max(indicators[r], 1.0e-15);
            
            // WENO-JS weight computation: α₍ᵣ₎ = d₍ᵣ₎/(ε + β₍ᵣ₎)ᵖ
            alphas[r] = optimalWeights(r) / 
                std::pow(safe_indicator + Traits::Epsilon, Traits::Power);
            alpha_sum += alphas[r];
        }
        
        // Normalize weights
        Set::Scalar inv_alpha_sum = (alpha_sum > 1.0e-15) ? (1.0 / alpha_sum) : 1.0e15;
        
        for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
            weights[r] = alphas[r] * inv_alpha_sum;
        }
        
        // Debug call - separated from the main algorithm logic
        Util::ScimitarX_Util::Debug::DebugWENOWeights(
            i, j, k, indicators, weights, 
            optimalWeights, Traits::Epsilon, Traits::Power, 
            enableDebug
        );
        
        return weights;
    }
    
    /**
     * Perform WENO reconstruction for a single point
     * 
     * @param stencil Array of stencil values
     * @param i, j, k Cell indices for debugging
     * @param enableDebug Enable detailed debugging output
     * @return Reconstructed value at the interface
     */
    template <typename StencilArray>
    AMREX_GPU_DEVICE
    static Set::Scalar ReconstructValue(
        const StencilArray& stencil,
        int i, int j, int k,  // Cell indices for debugging
        bool enableDebug = false
    ) {
        // 1. Compute smoothness indicators
        auto indicators = ComputeSmoothnessIndicators(stencil, i, j, k, enableDebug);
        
        // 2. Compute candidate reconstructions
        auto candidate_fluxes = ComputeCandidateFluxes(stencil, i, j, k, enableDebug);
        
        // 3. Compute non-linear weights
        auto weights = ComputeNonlinearWeights(indicators, i, j, k, enableDebug);
        
        // 4. Form convex combination for final reconstruction
        Set::Scalar reconstructed_value = 0.0;
        for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
            reconstructed_value += weights[r] * candidate_fluxes[r];
        }
        
        // Debug call for the final reconstruction
        Util::ScimitarX_Util::Debug::DebugWENOReconstruction(
            i, j, k, candidate_fluxes, weights, 
            Traits::GetOptimalWeights(), reconstructed_value,
            enableDebug
        );
        
        return reconstructed_value;
    }
};

/**
 * @class WENOReconstructionPolicy (Specialization for WENOZ5Traits)
 * @brief Specialized implementation for WENO-Z scheme
 * 
 * This specialization implements the WENO-Z weight calculation from Borges et al. (2008)
 * which uses a global smoothness indicator to improve accuracy near critical points.
 */
template<>
class WENOReconstructionPolicy<WENOZ5Traits> {
public:
    using Traits = WENOZ5Traits;
    
    // Use the same implementation for ComputeSmoothnessIndicators and ComputeCandidateFluxes
    template <typename StencilArray>
    AMREX_GPU_DEVICE
    static auto ComputeSmoothnessIndicators(
        const StencilArray& stencil, 
        int i, int j, int k,
        bool enableDebug = false
    ) {
        return WENOReconstructionPolicy<WENOJS5Traits>::ComputeSmoothnessIndicators(
            stencil, i, j, k, enableDebug);
    }
    
    template <typename StencilArray>
    AMREX_GPU_DEVICE
    static auto ComputeCandidateFluxes(
        const StencilArray& stencil,
        int i, int j, int k,
        bool enableDebug = false
    ) {
        return WENOReconstructionPolicy<WENOJS5Traits>::ComputeCandidateFluxes(
            stencil, i, j, k, enableDebug);
    }
    
    /**
     * Compute non-linear weights for WENO-Z scheme with enhanced robustness
     * 
     * The Z-variant uses a global smoothness indicator τ₅ = |β₀ - β₂|
     * to improve accuracy near smooth critical points.
     * 
     * @param indicators Smoothness indicators for each candidate stencil
     * @param i, j, k Cell indices for debugging
     * @param enableDebug Enable detailed debugging output
     * @return Array of normalized non-linear weights
     */
    template <typename IndicatorArray>
    AMREX_GPU_DEVICE
    static auto ComputeNonlinearWeights(
        const IndicatorArray& indicators,
        int i, int j, int k,
        bool enableDebug = false
    ) {
        std::array<Set::Scalar, Traits::NumCandidateStencils> weights;
        
        // Get optimal weights
        auto optimalWeights = Traits::GetOptimalWeights();
        
        // Compute global smoothness indicator tau = |β₀ - β₂|
        // This detects higher regularity in the solution
        Set::Scalar tau = std::abs(indicators[0] - indicators[2]);
        
        // Calculate alpha values and their sum
        Set::Scalar alpha_sum = 0.0;
        std::array<Set::Scalar, Traits::NumCandidateStencils> alphas;
        
        for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
            // Safe indicator value to prevent numerical issues
            Set::Scalar safe_indicator = std::max(indicators[r], 1.0e-15);
            
            // WENO-Z weight computation: α₍ᵣ₎ = d₍ᵣ₎(1 + τ/(β₍ᵣ₎ + ε))^p
            alphas[r] = optimalWeights(r) * 
                std::pow(1.0 + tau / (safe_indicator + Traits::Epsilon), Traits::Power);
            
            alpha_sum += alphas[r];
        }
        
        // Normalize weights with protection against division by zero
        Set::Scalar inv_alpha_sum = (alpha_sum > 1.0e-15) ? (1.0 / alpha_sum) : 1.0e15;
        
        for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
            weights[r] = alphas[r] * inv_alpha_sum;
        }
        
        // Debug call - separated from the main algorithm logic
        Util::ScimitarX_Util::Debug::DebugWENOWeights(
            i, j, k, indicators, weights, 
            optimalWeights, Traits::Epsilon, Traits::Power, 
            enableDebug
        );
        
        return weights;
    }
    
    /**
     * Perform WENO-Z reconstruction for a single point
     * 
     * @param stencil Array of stencil values
     * @param i, j, k Cell indices for debugging
     * @param enableDebug Enable detailed debugging output
     * @return Reconstructed value at the interface
     */
    template <typename StencilArray>
    AMREX_GPU_DEVICE
    static Set::Scalar ReconstructValue(
        const StencilArray& stencil,
        int i, int j, int k,
        bool enableDebug = false
    ) {
        // 1. Compute smoothness indicators
        auto indicators = ComputeSmoothnessIndicators(stencil, i, j, k, enableDebug);
        
        // 2. Compute candidate reconstructions
        auto candidate_fluxes = ComputeCandidateFluxes(stencil, i, j, k, enableDebug);
        
        // 3. Compute non-linear weights with WENO-Z method
        auto weights = ComputeNonlinearWeights(indicators, i, j, k, enableDebug);
        
        // 4. Form convex combination for final reconstruction
        Set::Scalar reconstructed_value = 0.0;
        for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
            reconstructed_value += weights[r] * candidate_fluxes[r];
        }
        
        // Debug call for the final reconstruction
        Util::ScimitarX_Util::Debug::DebugWENOReconstruction(
            i, j, k, candidate_fluxes, weights, 
            Traits::GetOptimalWeights(), reconstructed_value,
            enableDebug
        );
        
        return reconstructed_value;
    }
};

//==============================================================================
// WENO Reconstruction class (implements FluxReconstruction interface)
//==============================================================================

/**
 * @class WENOReconstruction
 * @brief WENO flux reconstruction implementation
 * 
 * This class implements the FluxReconstruction interface using WENO schemes.
 * It uses the WENOReconstructionPolicy with the specified traits to 
 * perform reconstruction at cell interfaces.
 * 
 * @tparam IntegratorDerivedClass The integrator class (e.g., ScimitarX)
 * @tparam VariantTraits Traits defining the WENO variant to use
 */
template<typename IntegratorDerivedClass, typename VariantTraits>
class WENOReconstruction : public FluxReconstruction<IntegratorDerivedClass> {
public:
    using Traits = VariantTraits;
    using Policy = WENOReconstructionPolicy<Traits>;
    
    WENOReconstruction() {}

    /**
     * Reconstruct flux at cell interfaces using WENO scheme
     * 
     * @param direction Spatial direction (0=x, 1=y, 2=z)
     * @param lev AMR refinement level
     * @param solver Pointer to the solver object
     * @param working_buffer MultiFab containing the values to reconstruct
     * @param QL_stencil MultiFab to store left-biased reconstructed values
     * @param QR_stencil MultiFab to store right-biased reconstructed values
     */
    void ReconstructFluxVariables(
        int direction, 
        int /*lev*/, 
        IntegratorDerivedClass* solver,
        amrex::MultiFab& working_buffer,
        amrex::MultiFab& QL_stencil, 
        amrex::MultiFab& QR_stencil
    ) const override {
        const int nghosts = solver->number_of_ghost_cells;
        const bool enableDebug = false;  // Enable detailed debugging

        // Process each grid patch
        for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
            const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);
            
            auto const& in_arr = working_buffer.const_array(mfi);
            auto const& ql_arr = QL_stencil.array(mfi);
            
            // Process each cell in the valid box
            amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {
                // Construct index arrays with correct number of components
                int index[3];
                int right_index[3];
                int lower_bounds[3];
                int upper_bounds[3];
                
                // Set up indices and bounds
                index[0] = iface; index[1] = jface; index[2] = kface;
                right_index[0] = iface; right_index[1] = jface; right_index[2] = kface;

                // Define bounds including ghost cells
                lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); 
                lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); 
                lower_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.smallEnd(2) : 0;
                upper_bounds[0] = face_bx_with_ghosts.bigEnd(0) -1;
                upper_bounds[1] = face_bx_with_ghosts.bigEnd(1) -1;
                upper_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.bigEnd(2) - 1 : 0;
                
                // Clamp the main index to ensure it's within bounds
                ClampIndices(index, lower_bounds, upper_bounds);

                int right_offset[3] = {0, 0, 0};
                right_offset[direction] = 1;  // Shift right

                ShiftAndClampIndices(right_index, right_offset, lower_bounds, upper_bounds, direction);
                
                // Process each component
                for (int n = 0; n < solver->number_of_components; ++n) {
                    // Extract stencil values based on direction
                    std::array<Set::Scalar, Traits::TotalPoints> stencil;
                    
                    // Base stencil center point - the left interface point
                    int center_idx[3] = {index[0], index[1], index[2]};
                    
                    // Extract the full stencil
                    for (int s = 0; s < Traits::TotalPoints; ++s) {
                        // Calculate stencil point index with proper offset
                        int stencil_idx[3] = {center_idx[0], center_idx[1], center_idx[2]};
                        
                        // Apply offset in the reconstruction direction
                        // Map s to range [-(TotalPoints/2), +(TotalPoints/2)]
                        int offset = s - Traits::TotalPoints/2;
                        stencil_idx[direction] += offset;
                        
                        // Clamp indices to ensure they're within bounds
                        ClampIndices(stencil_idx, lower_bounds, upper_bounds);
                        
                        // Extract data value for this stencil point
                        stencil[s] = in_arr(stencil_idx[0], stencil_idx[1], stencil_idx[2], n);
                    }
                    
                    // Reconstruct left state (QL) at interface using WENO
                    ql_arr(index[0], index[1], index[2], n) = 
                        Policy::ReconstructValue(stencil, iface, jface, kface, enableDebug);
                }
            });
        }

        // Process each grid patch to set QR_stencil
        for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
            const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);
            
            auto const& ql_arr = QL_stencil.array(mfi);
            auto const& qr_arr = QR_stencil.array(mfi);
            
            // Process each cell in the valid box
            amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {
                // Construct index arrays with correct number of components
                int index[3];
                int right_index[3];
                int lower_bounds[3];
                int upper_bounds[3];
                
                // Set up indices and bounds
                index[0] = iface; index[1] = jface; index[2] = kface;
                right_index[0] = iface; right_index[1] = jface; right_index[2] = kface;

                // Define bounds including ghost cells
                lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); 
                lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); 
                lower_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.smallEnd(2) : 0;
                upper_bounds[0] = face_bx_with_ghosts.bigEnd(0) -1;
                upper_bounds[1] = face_bx_with_ghosts.bigEnd(1) -1;
                upper_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.bigEnd(2) - 1 : 0;
                
                // Clamp the main index to ensure it's within bounds
                ClampIndices(index, lower_bounds, upper_bounds);

                int right_offset[3] = {0, 0, 0};
                right_offset[direction] = 1;  // Shift right

                ShiftAndClampIndices(right_index, right_offset, lower_bounds, upper_bounds, direction);
                
                // Process each component
                for (int n = 0; n < solver->number_of_components; ++n) {
                    // Set QR state for the neighboring cell's left interface (which is our right interface)
                    qr_arr(index[0], index[1], index[2], n) = 
                        ql_arr(right_index[0], right_index[1], right_index[2], n);
                }
            });
        }
    }    
};

//==============================================================================
// Convenience type aliases for common WENO variants
//==============================================================================

/**
 * @typedef WENOJS3
 * @brief 3rd order WENO-JS scheme
 */
template<typename IntegratorDerivedClass>
using WENOJS3 = WENOReconstruction<IntegratorDerivedClass, WENOJS3Traits>;

/**
 * @typedef WENOJS5
 * @brief 5th order WENO-JS scheme
 */
template<typename IntegratorDerivedClass>
using WENOJS5 = WENOReconstruction<IntegratorDerivedClass, WENOJS5Traits>;

/**
 * @typedef WENOZ5
 * @brief 5th order WENO-Z scheme
 */
template<typename IntegratorDerivedClass>
using WENOZ5 = WENOReconstruction<IntegratorDerivedClass, WENOZ5Traits>;

} // namespace Numeric

#endif // NUMERIC_WENO_RECONSTRUCTION_H

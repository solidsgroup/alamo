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
 * This structure defines the coefficients used to compute finite differences
 * for smoothness indicators in WENO schemes.
 * 
 * @tparam NumDifferenceTerms Number of difference terms in smoothness indicator
 * @tparam PointsPerStencil Number of points in each candidate stencil
 * @tparam NumCandidateStencils Number of candidate stencils
 */
template<int NumDifferenceTerms, int PointsPerStencil, int NumCandidateStencils>
struct DifferenceOperator {
    /**
     * Coefficients array with the following structure:
     * - First index: Candidate stencil (0 to NumCandidateStencils-1)
     * - Second index: Difference term (0 to NumDifferenceTerms-1)
     * - Third index: Stencil point coefficient (0 to PointsPerStencil-1)
     */
    static constexpr std::array<
        std::array<
            std::array<Set::Scalar, PointsPerStencil>, 
            NumDifferenceTerms
        >, 
        NumCandidateStencils
    > Coefficients = {};
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
    
    /**
     * Optimal (linear) weights for WENO-JS3
     * These weights would give 3rd order accuracy in smooth regions
     */
    static constexpr std::array<Set::Scalar, NumCandidateStencils> OptimalWeights = {
        1.0/3.0, 2.0/3.0
    };
    
    /**
     * Coefficients for candidate stencils reconstruction
     * First index: Candidate stencil (0 to NumCandidateStencils-1)
     * Second index: Coefficient for each point in the stencil
     */
    static constexpr std::array<std::array<Set::Scalar, PointsPerStencil>, NumCandidateStencils> 
    CandidateStencilCoefficients = {{
        {-1.0/2.0, 3.0/2.0},  // Left stencil S₀ coefficients
        {1.0/2.0, 1.0/2.0}    // Right stencil S₁ coefficients
    }};
    
    /**
     * Stencil offsets define which points to use for each candidate stencil
     * These are offsets relative to the interface being reconstructed
     */
    static constexpr std::array<std::array<int, PointsPerStencil>, NumCandidateStencils> 
    SubStencilOffsets = {{
        {-1, 0},  // Left stencil S₀ uses points (i-1, i)
        {0, 1}    // Right stencil S₁ uses points (i, i+1)
    }};
    
    /**
     * @struct WENOJS3DifferenceOperator
     * @brief Difference operators for WENO-JS3 smoothness indicators
     */
    struct WENOJS3DifferenceOperator : 
        public DifferenceOperator<NumDifferenceTerms, PointsPerStencil, NumCandidateStencils> {
        static constexpr std::array<
            std::array<
                std::array<Set::Scalar, PointsPerStencil>, 
                NumDifferenceTerms
            >, 
            NumCandidateStencils
        > Coefficients = {{
            {{ {1.0, -1.0} }},  // 1st difference for left stencil: f[i-1] - f[i]
            {{ {1.0, -1.0} }}   // 1st difference for right stencil: f[i] - f[i+1]
        }};
    };
    
    using DifferenceOperators = WENOJS3DifferenceOperator;
    
    /**
     * Calculate non-linear weights for WENO-JS3
     * 
     * @param indicators Smoothness indicators for each candidate stencil
     * @return Array of normalized non-linear weights
     */
    template <typename IndicatorArray>
    AMREX_GPU_DEVICE
    static auto ComputeNonlinearWeights(const IndicatorArray& indicators) {
        std::array<Set::Scalar, NumCandidateStencils> weights;
        Set::Scalar alpha_sum = 0.0;
        
        // WENO-JS weight computation: α₍ᵣ₎ = d₍ᵣ₎/(ε + β₍ᵣ₎)ᵖ
        for (int r = 0; r < NumCandidateStencils; ++r) {
            Set::Scalar alpha = OptimalWeights[r] / 
                std::pow(indicators[r] + Epsilon, Power);
            
            weights[r] = alpha;
            alpha_sum += alpha;
        }
        
        // Normalize weights to ensure they sum to 1
        const Set::Scalar inv_alpha_sum = 1.0 / alpha_sum;
        for (int r = 0; r < NumCandidateStencils; ++r) {
            weights[r] *= inv_alpha_sum;
        }
        
        return weights;
    }
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
    
    /**
     * Optimal (linear) weights for WENO-JS5
     * These weights would give 5th order accuracy in smooth regions
     */
    static constexpr std::array<Set::Scalar, NumCandidateStencils> OptimalWeights = {
        0.1, 0.6, 0.3
    };
    
    /**
     * Coefficients for candidate stencils reconstruction
     * First index: Candidate stencil (0 to NumCandidateStencils-1)
     * Second index: Coefficient for each point in the stencil
     */
    static constexpr std::array<std::array<Set::Scalar, PointsPerStencil>, NumCandidateStencils> 
    CandidateStencilCoefficients = {{
        {1.0/3.0, -7.0/6.0, 11.0/6.0},   // Left stencil S₀ coefficients
        {-1.0/6.0, 5.0/6.0, 1.0/3.0},    // Center stencil S₁ coefficients
        {1.0/3.0, 5.0/6.0, -1.0/6.0}     // Right stencil S₂ coefficients
    }};
    
    /**
     * Stencil offsets define which points to use for each candidate stencil
     * These are offsets relative to the interface being reconstructed
     */
    static constexpr std::array<std::array<int, PointsPerStencil>, NumCandidateStencils> 
    SubStencilOffsets = {{
        {-2, -1, 0},   // Left stencil S₀ uses points (i-2, i-1, i)
        {-1, 0, 1},    // Center stencil S₁ uses points (i-1, i, i+1)
        {0, 1, 2}      // Right stencil S₂ uses points (i, i+1, i+2)
    }};
    
    /**
     * @struct WENOJS5DifferenceOperator
     * @brief Difference operators for WENO-JS5 smoothness indicators
     * 
     * These operators implement the smoothness indicators from Jiang & Shu (1996):
     * β₀ = (13/12)(f₍ᵢ₋₂₎ - 2f₍ᵢ₋₁₎ + fᵢ)² + (1/4)(f₍ᵢ₋₂₎ - 4f₍ᵢ₋₁₎ + 3fᵢ)²
     * β₁ = (13/12)(f₍ᵢ₋₁₎ - 2fᵢ + f₍ᵢ₊₁₎)² + (1/4)(f₍ᵢ₋₁₎ - f₍ᵢ₊₁₎)²
     * β₂ = (13/12)(fᵢ - 2f₍ᵢ₊₁₎ + f₍ᵢ₊₂₎)² + (1/4)(3fᵢ - 4f₍ᵢ₊₁₎ + f₍ᵢ₊₂₎)²
     */
    struct WENOJS5DifferenceOperator : 
        public DifferenceOperator<NumDifferenceTerms, PointsPerStencil, NumCandidateStencils> {
        
        // Weights for the different squared terms
        static constexpr std::array<Set::Scalar, NumDifferenceTerms> TermWeights = {
            13.0/12.0,  // Weight for the first difference term
            1.0/4.0     // Weight for the second difference term
        };
        
        // Coefficients for each difference term in the smoothness indicators
        static constexpr std::array<
            std::array<
                std::array<Set::Scalar, PointsPerStencil>, 
                NumDifferenceTerms
            >, 
            NumCandidateStencils
        > Coefficients = {{
            {{ 
                {1.0, -2.0, 1.0},    // 1st diff for stencil 0: f[i-2] - 2*f[i-1] + f[i]
                {1.0, -4.0, 3.0}     // 2nd diff for stencil 0: f[i-2] - 4*f[i-1] + 3*f[i]
            }},
            {{ 
                {1.0, -2.0, 1.0},    // 1st diff for stencil 1: f[i-1] - 2*f[i] + f[i+1]
                {1.0, 0.0, -1.0}     // 2nd diff for stencil 1: f[i-1] - f[i+1]
            }},
            {{ 
                {1.0, -2.0, 1.0},    // 1st diff for stencil 2: f[i] - 2*f[i+1] + f[i+2]
                {3.0, -4.0, 1.0}     // 2nd diff for stencil 2: 3*f[i] - 4*f[i+1] + f[i+2]
            }}
        }};
    };
    
    using DifferenceOperators = WENOJS5DifferenceOperator;
    
    /**
     * Calculate non-linear weights for WENO-JS5
     * 
     * @param indicators Smoothness indicators for each candidate stencil
     * @return Array of normalized non-linear weights
     */
    template <typename IndicatorArray>
    AMREX_GPU_DEVICE
    static auto ComputeNonlinearWeights(const IndicatorArray& indicators) {
        std::array<Set::Scalar, NumCandidateStencils> weights;
        Set::Scalar alpha_sum = 0.0;
        
        // WENO-JS weight computation: α₍ᵣ₎ = d₍ᵣ₎/(ε + β₍ᵣ₎)ᵖ
        for (int r = 0; r < NumCandidateStencils; ++r) {
            Set::Scalar alpha = OptimalWeights[r] / 
                std::pow(indicators[r] + Epsilon, Power);
            
            weights[r] = alpha;
            alpha_sum += alpha;
        }
        
        // Normalize weights to ensure they sum to 1
        const Set::Scalar inv_alpha_sum = 1.0 / alpha_sum;
        for (int r = 0; r < NumCandidateStencils; ++r) {
            weights[r] *= inv_alpha_sum;
        }
        
        return weights;
    }
};

//==============================================================================
// WENO-Z5 Traits - 5th order WENO-Z scheme
//==============================================================================

/**
 * @class WENOZ5Traits
 * @brief Traits for 5th order WENO-Z scheme from Borges et al. (2008)
 * 
 * Extends the 5th order WENO-JS scheme with a global smoothness measure
 * to improve accuracy near critical points. Inherits most properties from
 * WENOJS5Traits but overrides the weight calculation.
 */
struct WENOZ5Traits : public WENOJS5Traits {
    /**
     * Calculate non-linear weights for WENO-Z5 scheme
     * 
     * The Z-variant uses a global smoothness indicator τ₅ = |β₀ - β₂|
     * to improve accuracy near smooth critical points.
     * 
     * @param indicators Smoothness indicators for each candidate stencil
     * @return Array of normalized non-linear weights
     */
    template <typename IndicatorArray>
    AMREX_GPU_DEVICE
    static auto ComputeNonlinearWeights(const IndicatorArray& indicators) {
        // Compute global smoothness indicator tau = |β₀ - β₂|
        // This detects higher regularity in the solution
        Set::Scalar tau = std::abs(indicators[0] - indicators[2]);
        
        std::array<Set::Scalar, NumCandidateStencils> weights;
        Set::Scalar alpha_sum = 0.0;
        
        // WENO-Z weight computation: α₍ᵣ₎ = d₍ᵣ₎(1 + τ/(β₍ᵣ₎ + ε))^p
        for (int r = 0; r < NumCandidateStencils; ++r) {
            Set::Scalar alpha = OptimalWeights[r] * 
                std::pow(1.0 + tau / (indicators[r] + Epsilon), Power);
            
            weights[r] = alpha;
            alpha_sum += alpha;
        }
        
        // Normalize weights
        const Set::Scalar inv_alpha_sum = 1.0 / alpha_sum;
        for (int r = 0; r < NumCandidateStencils; ++r) {
            weights[r] *= inv_alpha_sum;
        }
        
        return weights;
    }
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
     * @return Array of smoothness indicators for each candidate stencil
     */
    template <typename StencilArray>
    AMREX_GPU_DEVICE
    static auto ComputeSmoothnessIndicators(const StencilArray& stencil) {
        std::array<Set::Scalar, Traits::NumCandidateStencils> indicators;
        
        // Loop over each candidate stencil
        for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
            indicators[r] = 0.0;
            
            // Loop over each difference term in the smoothness indicator
            for (int term = 0; term < static_cast<int>(Traits::DifferenceOperators::TermWeights.size()); ++term) {
                // Compute the weighted squared difference for this term
                Set::Scalar diff_value = 0.0;
                
                // Apply the difference operator for this term
                for (int p = 0; p < Traits::PointsPerStencil; ++p) {
                    int stencil_idx = Traits::SubStencilOffsets[r][p] + Traits::TotalPoints/2;
                    diff_value += Traits::DifferenceOperators::Coefficients[r][term][p] * stencil[stencil_idx];
                }
                
                // Add the weighted squared difference to the indicator
                indicators[r] += Traits::DifferenceOperators::TermWeights[term] * diff_value * diff_value;
            }
        }
        
        return indicators;
    }
    
    /**
     * Compute candidate fluxes for each stencil
     * 
     * @param stencil Array of stencil values
     * @return Array of candidate reconstructed values for each stencil
     */
    template <typename StencilArray>
    AMREX_GPU_DEVICE
    static auto ComputeCandidateFluxes(const StencilArray& stencil) {
        std::array<Set::Scalar, Traits::NumCandidateStencils> fluxes;
        
        // Compute reconstruction for each candidate stencil
        for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
            fluxes[r] = 0.0;
            
            // Apply reconstruction formula using the coefficients
            for (int p = 0; p < Traits::PointsPerStencil; ++p) {
                int stencil_idx = Traits::SubStencilOffsets[r][p] + Traits::TotalPoints/2;
                fluxes[r] += Traits::CandidateStencilCoefficients[r][p] * stencil[stencil_idx];
            }
        }
        
        return fluxes;
    }
    
    /**
     * Perform WENO reconstruction for a single point
     * 
     * @param stencil Array of stencil values
     * @return Reconstructed value at the interface
     */
    template <typename StencilArray>
    AMREX_GPU_DEVICE
    static Set::Scalar ReconstructValue(
        const StencilArray& stencil,
        int i, int j, int k,  // Cell indices for debugging
        bool debugMode = true
    ) {
        // 1. Compute smoothness indicators
        auto indicators = ComputeSmoothnessIndicators(stencil);
        
        // Debug smoothness indicators if needed
        if (debugMode && Util::ScimitarX_Util::Debug::IsTargetLocation(i, j, k)) {
            Util::Message(INFO, "WENO Smoothness Indicators at (" + 
                std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + "):");
            for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
                Util::Message(INFO, "  β[" + std::to_string(r) + "] = " + std::to_string(indicators[r]));
                // Check for NaN or Inf in smoothness indicators
                Util::ScimitarX_Util::Debug::CheckNaNAndAbort(i, j, k, indicators[r], 
                    "WENO smoothness indicator", "Reconstruction");
            }
        }
        
        // 2. Compute candidate reconstructions
        auto candidate_fluxes = ComputeCandidateFluxes(stencil);
        
        // Debug candidate fluxes if needed
        if (debugMode && Util::ScimitarX_Util::Debug::IsTargetLocation(i, j, k)) {
            Util::Message(INFO, "WENO Candidate Fluxes at (" + 
                std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + "):");
            for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
                Util::Message(INFO, "  flux[" + std::to_string(r) + "] = " + std::to_string(candidate_fluxes[r]));
                // Check for NaN or Inf in candidate fluxes
                Util::ScimitarX_Util::Debug::CheckNaNAndAbort(i, j, k, candidate_fluxes[r], 
                    "WENO candidate flux", "Reconstruction");
            }
        }
        
        // 3. Compute non-linear weights
        auto weights = Traits::ComputeNonlinearWeights(indicators);
        
        // Debug weights if needed
        if (debugMode && Util::ScimitarX_Util::Debug::IsTargetLocation(i, j, k)) {
            Util::Message(INFO, "WENO Weights at (" + 
                std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + "):");
            Set::Scalar sum = 0.0;
            for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
                Util::Message(INFO, "  ω[" + std::to_string(r) + "] = " + std::to_string(weights[r]));
                sum += weights[r];
                // Check for NaN or Inf in weights
                Util::ScimitarX_Util::Debug::CheckNaNAndAbort(i, j, k, weights[r], 
                    "WENO weight", "Reconstruction");
            }
            Util::Message(INFO, "  Sum of weights: " + std::to_string(sum));
            
            // Check if weights sum to 1.0 (allowing for small numerical errors)
            if (std::abs(sum - 1.0) > 1.0e-14) {
                Util::Warning(INFO, "WENO weights do not sum to 1.0 at (" + 
                    std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + 
                    "). Sum = " + std::to_string(sum));
            }
        }
        
        // 4. Form convex combination for final reconstruction
        Set::Scalar reconstructed_value = 0.0;
        for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
            reconstructed_value += weights[r] * candidate_fluxes[r];
        }
        
        // Debug reconstructed value if needed
        if (debugMode && Util::ScimitarX_Util::Debug::IsTargetLocation(i, j, k)) {
            Util::Message(INFO, "WENO Reconstructed Value at (" + 
                std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + 
                "): " + std::to_string(reconstructed_value));
            // Check for NaN, Inf, or negative value
            Util::ScimitarX_Util::Debug::CheckNaNAndAbort(i, j, k, reconstructed_value, 
                "WENO reconstructed value", "Final Reconstruction");
        }
        
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

        Util::ScimitarX_Util::Debug::SetTargetDebugLocationIndices(100,5,0);     

        // Process each grid patch
        for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
            const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);
            
            auto const& in_arr = working_buffer.const_array(mfi);
            auto const& ql_arr = QL_stencil.array(mfi);
            
            // Find the min and max stencil offsets
            int min_stencil_offset = 0;
            int max_stencil_offset = 0;
            
            // Analyze each sub-stencil to find the overall stencil extents
            for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
                for (int p = 0; p < Traits::PointsPerStencil; ++p) {
                    int offset = Traits::SubStencilOffsets[r][p];
                    min_stencil_offset = std::min(min_stencil_offset, offset);
                    max_stencil_offset = std::max(max_stencil_offset, offset);
                }
            }
            
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
                    
                    // Debug stencil center if enabled
                    if (enableDebug && Util::ScimitarX_Util::Debug::IsTargetLocation(iface, jface, kface)) {
                        Util::Message(INFO, "WENO Stencil Center at (" + 
                            std::to_string(center_idx[0]) + ", " + 
                            std::to_string(center_idx[1]) + ", " + 
                            std::to_string(center_idx[2]) + ") for component " + 
                            std::to_string(n));
                    }
                    
                    // Extract the full stencil using the min and max offsets
                    for (int s = 0; s < Traits::TotalPoints; ++s) {
                        // Calculate stencil point index with proper offset
                        int stencil_idx[3] = {center_idx[0], center_idx[1], center_idx[2]};
                        
                        // Apply offset in the reconstruction direction
                        // Map s to range [min_offset, max_offset]
                        int offset = min_stencil_offset + s;
                        stencil_idx[direction] = center_idx[direction] + offset;
                        
                        // Clamp indices to ensure they're within bounds
                        ClampIndices(stencil_idx, lower_bounds, upper_bounds);
                        
                        // Extract data value for this stencil point
                        stencil[s] = in_arr(stencil_idx[0], stencil_idx[1], stencil_idx[2], n);

                        // Debug stencil center if enabled
                        if (enableDebug && Util::ScimitarX_Util::Debug::IsTargetLocation(iface, jface, kface)) {
                            Util::Message(INFO, "WENO Stencil Center at (" + 
                            std::to_string(stencil_idx[0]) + ", " + 
                            std::to_string(stencil_idx[1]) + ", " + 
                            std::to_string(stencil_idx[2]) + ") for component " + 
                            std::to_string(n));
                        }

                        // Debug stencil point values if enabled
                        if (enableDebug && Util::ScimitarX_Util::Debug::IsTargetLocation(iface, jface, kface)) {
                            Util::ScimitarX_Util::Debug::DebugValuesIfTarget(
                                iface, jface, kface,
                                stencil[s],
                                "Stencil Point " + std::to_string(s),
                                "WENO Reconstruction", 
                                true,  // Abort if NaN
                                n      // Component
                            );
                        }
                    }                    
                    // Reconstruct left state (QL) at interface using WENO
                    ql_arr(index[0], index[1], index[2], n) = 
                        Policy::ReconstructValue(stencil, iface, jface, kface, enableDebug);                    
                }
            });
        }

        // Process each grid patch
        for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
            const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);
            
            auto const& in_arr = working_buffer.const_array(mfi);
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
                    
                    // Check for negative density or internal energy after reconstruction
                   /* if (n == solver->variableIndex.DENS || n == solver->variableIndex.IE) {
                        Set::Scalar reconstructed_value = ql_arr(index[0], index[1], index[2], n);
                        
                        // Apply a limiter to ensure physical values - prevent negative density/energy
                        if (reconstructed_value < 0.0) {
                            if (enableDebug) {
                                Util::Warning(INFO, "Negative value detected in WENO reconstruction at (" +
                                    std::to_string(iface) + ", " + std::to_string(jface) + ", " + 
                                    std::to_string(kface) + ") for component " + std::to_string(n) + 
                                    ". Value: " + std::to_string(reconstructed_value) + 
                                    ". Applying positivity-preserving limiter.");
                            }
                            
                            // Apply positivity-preserving limiter - use first-order reconstruction
                            // This is a simple fallback to first-order upwind in problematic regions
                               const Set::Scalar QL = in_arr(index[0], index[1], index[2], n);
                               const Set::Scalar QR = in_arr(right_index[0], right_index[1], right_index[2], n);
                                
                                qr_arr(index[0], index[1], index[2], n) = QR;
                                ql_arr(index[0], index[1], index[2], n) = QL;
                        }
                    }*/
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

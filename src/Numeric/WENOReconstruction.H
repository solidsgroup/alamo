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
    static constexpr Set::Scalar Epsilon = 1e-40;  // Small parameter to avoid division by zero
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
    static constexpr Set::Scalar Epsilon = 1e-40;  // Small parameter to avoid division by zero
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
        bool enableDebug
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
            Traits::NumCandidateStencils,
            Traits::TotalPoints,
            Traits::NumDifferenceTerms,
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
        bool enableDebug
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
            Traits::NumCandidateStencils,
            Traits::PointsPerStencil,
            Traits::TotalPoints,
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
        bool enableDebug
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
            Traits::NumCandidateStencils,
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
    template <typename StencilMatrixType>
    AMREX_GPU_DEVICE
    static Set::MultiVector ReconstructValue(
        const StencilMatrixType& input_stencil,
        int i, int j, int k,
        int direction [[maybe_unused]],
        bool enableDebug)
    {

        int num_components = input_stencil.cols();
        // Apply WENO reconstruction to each component (no branching in the core algorithm)
        Set::MultiVector reconstructed_values(num_components);
        reconstructed_values.setZero();

        for (int comp = 0; comp < num_components; ++comp) {
            // Extract stencil values for this component
            std::array<Set::Scalar, Traits::TotalPoints> component_stencil;
            for (int s = 0; s < Traits::TotalPoints; ++s) {
                component_stencil[s] = input_stencil(s, comp);
            }

            // Apply core WENO algorithm to this component (no branching)
            auto indicators = ComputeSmoothnessIndicators(component_stencil, i, j, k, enableDebug);
            auto candidate_fluxes = ComputeCandidateFluxes(component_stencil, i, j, k, enableDebug);
            auto weights = ComputeNonlinearWeights(indicators, i, j, k, enableDebug);

            // Form convex combination of candidate fluxes (no branching)
            Set::Scalar reconstructed_component = 0.0;
            for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
                reconstructed_component += weights[r] * candidate_fluxes[r];
            }

            // Store the reconstructed value
            reconstructed_values(comp) = reconstructed_component;
        }

        return reconstructed_values;
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
        bool enableDebug
    ) {
        return WENOReconstructionPolicy<WENOJS5Traits>::ComputeSmoothnessIndicators(
            stencil, i, j, k, enableDebug);
    }

    template <typename StencilArray>
    AMREX_GPU_DEVICE
    static auto ComputeCandidateFluxes(
        const StencilArray& stencil,
        int i, int j, int k,
        bool enableDebug
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
        bool enableDebug
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
            Traits::NumCandidateStencils, 
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
    template <typename StencilMatrixType>
    AMREX_GPU_DEVICE
    static Set::MultiVector ReconstructValue(
        const StencilMatrixType& input_stencil,
        int i, int j, int k,
        int direction [[maybe_unused]],
        bool enableDebug)
    {
        int num_components = input_stencil.cols();

        // Apply WENO reconstruction to each component (no branching in the core algorithm)
        Set::MultiVector reconstructed_values(num_components);
        reconstructed_values.setZero();

        for (int comp = 0; comp < num_components; ++comp) {
            // Extract stencil values for this component
            std::array<Set::Scalar, Traits::TotalPoints> component_stencil;
            for (int s = 0; s < Traits::TotalPoints; ++s) {
                component_stencil[s] = input_stencil(s, comp);
            }

            // Apply core WENO algorithm to this component (no branching)
            auto indicators = ComputeSmoothnessIndicators(component_stencil, i, j, k, enableDebug);
            auto candidate_fluxes = ComputeCandidateFluxes(component_stencil, i, j, k, enableDebug);
            auto weights = ComputeNonlinearWeights(indicators, i, j, k, enableDebug);

            // Form convex combination of candidate fluxes (no branching)
            Set::Scalar reconstructed_component = 0.0;
            for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
                reconstructed_component += weights[r] * candidate_fluxes[r];
            }

            // Store the reconstructed value
            reconstructed_values(comp) = reconstructed_component;

            // Debug call for the final reconstruction
            Util::ScimitarX_Util::Debug::DebugWENOReconstruction(
                i, j, k, candidate_fluxes, weights, 
                Traits::GetOptimalWeights(), reconstructed_component,
                Traits::NumCandidateStencils,
                enableDebug);

        }

        return reconstructed_values;

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

    // Helper function to get variant name, only used when debug is enabled
    std::string GetVariantName() const {
        if constexpr (std::is_same_v<VariantTraits, WENOJS3Traits>) {
            return "WENOJS3";
        } else if constexpr (std::is_same_v<VariantTraits, WENOJS5Traits> && 
                            !std::is_same_v<VariantTraits, WENOZ5Traits>) {
            return "WENOJS5";
        } else if constexpr (std::is_same_v<VariantTraits, WENOZ5Traits>) {
            return "WENOZ5";
        } else {
            return "Unknown";
        }
    }


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
    void ReconstructFluxes(
        int direction, 
        int /*lev*/, 
        IntegratorDerivedClass* solver,
        const amrex::MultiFab& Flux_Negative,
        const amrex::MultiFab& Flux_Positive, 
        amrex::MultiFab& SummedFlux,
        const std::shared_ptr<GenericVariableAccessor>& variable_accessor,
        const amrex::MultiFab* AverageStateBuffer
    ) const override {
        const int nghosts = solver->number_of_ghost_cells;
        const int num_components = solver->number_of_components;
        const bool enableDebug = false;  // Enable detailed debugging

        // Get Reconstruction Mode 
        ReconstructionMode mode = variable_accessor->getCurrentReconstructionMode();

        // Create a bool to check for ReconstructionMode and if AverageStateBuffer is not null pointer
        bool useCharacteristic = (mode == ReconstructionMode::Characteristic && AverageStateBuffer != nullptr);

        //Util::ScimitarX_Util::Debug::SetTargetDebugLocationIndices(101, 10, 0);  // Set target to cell (10,10,0)
        // Only log the variant when debug is enabled
        if (enableDebug) {
            Util::Message(INFO, "ReconstructFluxes using WENO variant: " + GetVariantName());
        }

        // Process each grid patch
        for (amrex::MFIter mfi(SummedFlux, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);

            auto const& FluxNeg_arr = Flux_Negative.const_array(mfi);
            auto const& FluxPos_arr = Flux_Positive.const_array(mfi);
            auto const& TotalFlux_arr = SummedFlux.array(mfi);

            // Get average state array if available
            amrex::Array4<const Set::Scalar> Avg_arr;
            if (useCharacteristic) {
                Avg_arr = AverageStateBuffer->const_array(mfi);
            }

            // Process each cell in the valid box
            amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {
                // Construct index arrays with correct number of components
                int index[3];
                int lower_bounds[3];
                int upper_bounds[3];

                // Set up indices and bounds
                index[0] = iface; index[1] = jface; index[2] = kface;

                // Define bounds including ghost cells
                lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); 
                lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); 
                lower_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.smallEnd(2) : 0;
                upper_bounds[0] = face_bx_with_ghosts.bigEnd(0) -1;
                upper_bounds[1] = face_bx_with_ghosts.bigEnd(1) -1;
                upper_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.bigEnd(2) - 1 : 0;

                // Clamp the main index to ensure it's within bounds
                ClampIndices(index, lower_bounds, upper_bounds);

                // Store average state buffer in Vector Form
                Set::MultiVector avg_state(num_components);
                avg_state.setZero();

                // Store Fpos and Fneg Reconstructed
                Set::MultiVector Fpos_reconstructed(num_components);
                Set::MultiVector Fneg_reconstructed(num_components);
                Set::MultiVector Ftotal_reconstructed(num_components);
                Fpos_reconstructed.setZero();
                Fneg_reconstructed.setZero();
                Ftotal_reconstructed.setZero();

                if (useCharacteristic) {
                avg_state = Numeric::FieldToMultiVector(Avg_arr, index[0], index[1], index[2], num_components);
                }    

                //Create Stencil Matrix where each row represents a full state vector
                Set::MultiMatrix stencilpos_matrix(Traits::TotalPoints, num_components);
                Set::MultiMatrix stencilneg_matrix(Traits::TotalPoints, num_components);
                stencilpos_matrix.setZero();
                stencilneg_matrix.setZero();

                //Fill stencil matrix with state vectors from each stencil point 
                for (int s = 0; s < Traits::TotalPoints; ++s) {
                    // Calculate stencil point index with proper offset
                    int stencil_idx[3] = {index[0], index[1], index[2]};

                    // Apply offset in the reconstruction direction
                    // Map s to range [-(TotalPoints/2), +(TotalPoints/2)]
                    int offset = s - Traits::TotalPoints/2;
                    stencil_idx[direction] += offset;

                    // Clamp indices to ensure they're within bounds
                    ClampIndices(stencil_idx, lower_bounds, upper_bounds);

                    for (int n = 0; n < num_components; ++n) {
                    // Extract data value for this stencil point
                    stencilneg_matrix(s, n) = FluxNeg_arr(stencil_idx[0], stencil_idx[1], stencil_idx[2], n);
                    stencilpos_matrix(s, n) = FluxPos_arr(stencil_idx[0], stencil_idx[1], stencil_idx[2], n);

                    }
                }

                // Prepare stencil data for reconstruction
                // Create a copy regardless of whether we're using characteristic variables
                Set::MultiMatrix stencilneg = stencilneg_matrix;
                Set::MultiMatrix stencilpos = stencilpos_matrix;

                // Conditionally transform to characteristic space (branch happens once for all threads)
                if (useCharacteristic && variable_accessor) {
                    stencilneg = variable_accessor->TransformStencilToCharacteristic(    
                        index[0], index[1], index[2], stencilneg_matrix, direction, avg_state);
                    stencilpos = variable_accessor->TransformStencilToCharacteristic(    
                        index[0], index[1], index[2], stencilpos_matrix, direction, avg_state);
                }

                // Reconstruct Flux Positive and Negative at interface using WENO
                Fpos_reconstructed = Policy::ReconstructValue(stencilpos, index[0], index[1], index[2], direction, enableDebug);
                Fneg_reconstructed = Policy::ReconstructValue(stencilneg, index[0], index[1], index[2], direction, enableDebug);

                // Compute final summed flux at face center
                for (int n = 0; n < solver->number_of_components; ++n) {
                    Ftotal_reconstructed(n) = Fpos_reconstructed(n) + Fneg_reconstructed(n);
                }

                // Conditionally transform back to physical space (branch happens once for all threads)
                if (useCharacteristic && variable_accessor) {
                    Ftotal_reconstructed = variable_accessor->TransformFromCharacteristic(
                        Ftotal_reconstructed, direction, avg_state);
                }

                // Compute final summed flux at face center
                for (int n = 0; n < solver->number_of_components; ++n) {
                    TotalFlux_arr(index[0], index[1], index[2], n) = Ftotal_reconstructed(n);
                }

            });
        }

    SummedFlux.FillBoundary();

    }    

void ReconstructVariables(
    int direction, 
    int /*lev*/, 
    IntegratorDerivedClass* solver,
    const amrex::MultiFab& VariableBuffer,
    amrex::MultiFab& LeftStates, 
    amrex::MultiFab& RightStates,
    const std::shared_ptr<GenericVariableAccessor>& variable_accessor,
    const amrex::MultiFab* AverageStateBuffer
) const override {
    const int nghosts = solver->number_of_ghost_cells;
    const int num_components = solver->number_of_components;
    const bool enableDebug = false;  // Enable detailed debugging

    ReconstructionMode mode = variable_accessor->getCurrentReconstructionMode();
    // Create a bool to check for ReconstructionMode and if AverageStateBuffer is not null pointer
    bool useCharacteristic = (mode == ReconstructionMode::Characteristic && AverageStateBuffer != nullptr);

    // Create average state vector once outside the loop if needed
    Set::MultiVector avg_state;
    if (useCharacteristic) {
        //Initialize the average state outside the component loop
        avg_state.resize(num_components);
    }    

    //Util::ScimitarX_Util::Debug::SetTargetDebugLocationIndices(101, 10, 0);  // Set target to cell (10,10,0)
    // Only log the variant when debug is enabled
    if (enableDebug) {
        Util::Message(INFO, "ReconstructFluxes using WENO variant: " + GetVariantName());
    }

    // First pass: Compute left states using WENO reconstruction
    for (amrex::MFIter mfi(VariableBuffer, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);

        auto const& Var_arr = VariableBuffer.const_array(mfi);
        auto const& LeftStates_arr = LeftStates.array(mfi);
        auto const& RightStates_arr = RightStates.array(mfi);

        // Get average state array if available
        amrex::Array4<const Set::Scalar> Avg_arr;
        if (useCharacteristic) {
            Avg_arr = AverageStateBuffer->const_array(mfi);
        }

        amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {
            // Construct index arrays
            int index[3];
            int lower_bounds[3];
            int upper_bounds[3];

            // Set up indices
            index[0] = iface; index[1] = jface; index[2] = kface;

            // Define bounds including ghost cells
            lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); 
            lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); 
            lower_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.smallEnd(2) : 0;
            upper_bounds[0] = face_bx_with_ghosts.bigEnd(0) - 1;
            upper_bounds[1] = face_bx_with_ghosts.bigEnd(1) - 1;
            upper_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.bigEnd(2) - 1 : 0;

            // Clamp the main index
            ClampIndices(index, lower_bounds, upper_bounds);

            // store average state buffer in vector form  
            Set::MultiVector avg_state(num_components);
            avg_state.setZero();

            //Store left face reconstructed values in an array
            Set::MultiVector leftvalue_reconstructed(num_components);
            leftvalue_reconstructed.setZero();

            if (useCharacteristic) {
                avg_state = Numeric::FieldToMultiVector(Avg_arr, index[0], index[1], index[2], num_components);    
            }   

            // Create Stencil Matrix Where each row represents a full state vector
            Set::MultiMatrix stencil_matrix(Traits::TotalPoints, num_components);

            //Fill stencil matrix with state vectors from each stencil point 
            for (int s = 0; s < Traits::TotalPoints; ++s) {
                    // Calculate stencil point index with proper offset
                    int stencil_idx[3] = {index[0], index[1], index[2]};

                    // Apply offset in the reconstruction direction
                    // Map s to range [-(TotalPoints/2), +(TotalPoints/2)]
                    int offset = s - Traits::TotalPoints/2;
                    stencil_idx[direction] += offset;

                    // Clamp indices to ensure they're within bounds
                    ClampIndices(stencil_idx, lower_bounds, upper_bounds);

                    for (int n = 0; n < num_components; ++n) {
                    // Extract data value for this stencil point
                    stencil_matrix(s, n) = Var_arr(stencil_idx[0], stencil_idx[1], stencil_idx[2], n);

                    }
            }

            // Prepare stencil data for reconstruction
            // Create a copy regardless of whether we're using characteristic variables
            Set::MultiMatrix input_stencil = stencil_matrix;

            // Conditionally transform to characteristic space (branch happens once for all threads)
            if (useCharacteristic && variable_accessor) {
                input_stencil = variable_accessor->TransformStencilToCharacteristic(
                index[0], index[1], index[2], stencil_matrix, direction, avg_state);
            }

            leftvalue_reconstructed = Policy::ReconstructValue(input_stencil, index[0], index[1], index[2], direction, enableDebug);

            // Conditionally transform back to physical space (branch happens once for all threads)
            if (useCharacteristic && variable_accessor) {
                        leftvalue_reconstructed = variable_accessor->TransformFromCharacteristic(
                        leftvalue_reconstructed, direction, avg_state);
            }

            // Process each component
            for (int n = 0; n < num_components; ++n) {                
                // Store reconstructed left state
                LeftStates_arr(index[0], index[1], index[2], n) = leftvalue_reconstructed(n);
            }

            // ********** CHANGE - BEGIN **************

            //  Compute Right State (q_R) at face 'iface' 
            Set::MultiMatrix stencil_R_matrix(Traits::TotalPoints, num_components);
            for (int s = 0; s < Traits::TotalPoints; ++s) {
                int stencil_idx[3] = {iface, jface, kface};
                // Stencil for right state at face 'iface' is {i-1, i, i+1, i+2, i+3} where i = iface
                int offset = s - 1; 
                stencil_idx[direction] += offset;
                ClampIndices(stencil_idx, lower_bounds, upper_bounds);
                for (int n = 0; n < num_components; ++n) {
                    stencil_R_matrix(s, n) = Var_arr(stencil_idx[0], stencil_idx[1], stencil_idx[2], n);
                }
            }

            // To compute the right state with a left-biased policy, we reverse the stencil data.
            Set::MultiMatrix input_stencil_R = stencil_R_matrix.colwise().reverse();

            if (useCharacteristic) {
                // The transform should be on the original, non-reversed stencil data for physical correctness.
                Set::MultiMatrix temp_transform_stencil = stencil_R_matrix;
                input_stencil_R = variable_accessor->TransformStencilToCharacteristic(index[0], index[1], index[2], temp_transform_stencil, direction, avg_state);
                input_stencil_R.colwise().reverseInPlace();
            }

            Set::MultiVector right_reconstructed = Policy::ReconstructValue(input_stencil_R, index[0], index[1], index[2], direction, enableDebug);
            if (useCharacteristic) {
                right_reconstructed = variable_accessor->TransformFromCharacteristic(right_reconstructed, direction, avg_state);
            }
            for (int n = 0; n < num_components; ++n) {
                RightStates_arr(index[0], index[1], index[2], n) = right_reconstructed(n);
            }

        });
    }


    LeftStates.FillBoundary();
    RightStates.FillBoundary();
}

};          // ********** CHANGE - END ************

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

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
#include "Numeric/NumericTypes.H"
#include "Numeric/IntegratorVariableAccessLayer.H"
#include "Numeric/FluxHandler.H"
#include "Numeric/NumericFactory.H"

namespace Numeric {

// Forward declarations
template<typename VariantTraits> class WENOReconstructionPolicy;

//==============================================================================
// Base traits template that defines common interface for all WENO variants
//==============================================================================
template <int Order>
struct WENOBaseTraits {
    static constexpr int OrderValue = Order;
    
    // These constants and methods should be defined by specific variants
    static constexpr Set::Scalar Epsilon = 1.0e-6;
    static constexpr int NumCandidateStencils = 0;
    static constexpr int PointsPerStencil = 0;
    static constexpr int TotalPoints = 0;
};

//==============================================================================
// WENO-JS5 Specific Traits
//==============================================================================
struct WENOJS5Traits : public WENOBaseTraits<5> {
    static constexpr Set::Scalar Epsilon = 1.0e-6;
    static constexpr int NumCandidateStencils = 3;
    static constexpr int PointsPerStencil = 3;
    static constexpr int TotalPoints = 5;
    static constexpr int Power = 2;
    
    // Optimal linear weights for JS variant
    static constexpr std::array<Set::Scalar, NumCandidateStencils> OptimalWeights = {
        0.1, 0.6, 0.3
    };
    
    // Coefficients for candidate stencils (3 stencils with 3 points each)
    static constexpr std::array<std::array<Set::Scalar, PointsPerStencil>, NumCandidateStencils> 
    CandidateStencilCoefficients = {{
        {1.0/3.0, -7.0/6.0, 11.0/6.0},   // Left stencil S₀
        {-1.0/6.0, 5.0/6.0, 1.0/3.0},    // Center stencil S₁
        {1.0/3.0, 5.0/6.0, -1.0/6.0}     // Right stencil S₂
    }};
    
    // Stencil offsets define which points to use for each candidate stencil
    static constexpr std::array<std::array<int, PointsPerStencil>, NumCandidateStencils> 
    SubStencilOffsets = {{
        {0, 1, 2},   // Left stencil S₀ uses points i-2, i-1, i
        {1, 2, 3},   // Center stencil S₁ uses points i-1, i, i+1
        {2, 3, 4}    // Right stencil S₂ uses points i, i+1, i+2
    }};
    
    // Smoothness indicator coefficients for WENO5
    struct SmoothnessCoefficients {

        static constexpr std::array<Set::Scalar, 2> Stencil0 = {13.0/12.0, 1.0/4.0};
        
        static constexpr std::array<Set::Scalar, 2> Stencil1 = {13.0/12.0, 1.0/4.0};
        
        static constexpr std::array<Set::Scalar, 2> Stencil2 = {13.0/12.0, 1.0/4.0};
    };
    
    // Nonlinear weights calculation method - specific to JS
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
// WENO-Z5 Specific Traits - inherits most from JS5 but modifies weight calculation
//==============================================================================
struct WENOZ5Traits : public WENOJS5Traits {
    // Override nonlinear weights calculation method - Z-variant specific
    template <typename IndicatorArray>
    AMREX_GPU_DEVICE
    static auto ComputeNonlinearWeights(const IndicatorArray& indicators) {
        // Compute global smoothness indicator tau = |β₀ - β₂|
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
template<typename VariantTraits>
class WENOReconstructionPolicy {
public:
    using Traits = VariantTraits;
    
    // Compute smoothness indicators - can be specialized for efficiency depending on the variant
    template <typename StencilArray>
    AMREX_GPU_DEVICE
    static auto ComputeSmoothnessIndicators(const StencilArray& stencil) {
        std::array<Set::Scalar, Traits::NumCandidateStencils> indicators;
        
        if constexpr (Traits::OrderValue == 5) {
            // Efficient WENO5 smoothness indicator calculation
            // β₀ = (13/12)(f₍ᵢ₋₂₎ - 2f₍ᵢ₋₁₎ + fᵢ)² + (1/4)(f₍ᵢ₋₂₎ - 4f₍ᵢ₋₁₎ + 3fᵢ)²
            Set::Scalar diff1 = stencil[0] - 2.0*stencil[1] + stencil[2];
            Set::Scalar diff2 = stencil[0] - 4.0*stencil[1] + 3.0*stencil[2];
            indicators[0] = Traits::SmoothnessCoefficients::Stencil0[0] * diff1*diff1 + 
                           Traits::SmoothnessCoefficients::Stencil0[1] * diff2*diff2;
            
            // β₁ = (13/12)(f₍ᵢ₋₁₎ - 2fᵢ + f₍ᵢ₊₁₎)² + (1/4)(f₍ᵢ₋₁₎ - f₍ᵢ₊₁₎)²
            diff1 = stencil[1] - 2.0*stencil[2] + stencil[3];
            diff2 = stencil[1] - stencil[3];
            indicators[1] = Traits::SmoothnessCoefficients::Stencil1[0] * diff1*diff1 + 
                           Traits::SmoothnessCoefficients::Stencil1[1] * diff2*diff2;
            
            // β₂ = (13/12)(fᵢ - 2f₍ᵢ₊₁₎ + f₍ᵢ₊₂₎)² + (1/4)(3fᵢ - 4f₍ᵢ₊₁₎ + f₍ᵢ₊₂₎)²
            diff1 = stencil[2] - 2.0*stencil[3] + stencil[4];
            diff2 = 3.0*stencil[2] - 4.0*stencil[3] + stencil[4];
            indicators[2] = Traits::SmoothnessCoefficients::Stencil2[0] * diff1*diff1 + 
                           Traits::SmoothnessCoefficients::Stencil2[1] * diff2*diff2;
        } else if constexpr (Traits::OrderValue == 3) {
            // WENO3 smoothness indicators (if needed)
            indicators[0] = std::pow(stencil[0] - stencil[1], 2);
            indicators[1] = std::pow(stencil[1] - stencil[2], 2);
        } else {
            // Generic fallback for other orders (could be optimized further)
            for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
                indicators[r] = 0.0;
                for (int p = 0; p < Traits::PointsPerStencil; ++p) {
                    for (int q = p + 1; q < Traits::PointsPerStencil; ++q) {
                        int idx_p = Traits::SubStencilOffsets[r][p];
                        int idx_q = Traits::SubStencilOffsets[r][q];
                        Set::Scalar diff = stencil[idx_p] - stencil[idx_q];
                        indicators[r] += diff * diff;
                    }
                }
            }
        }
        
        return indicators;
    }
    
    // Compute candidate fluxes for each stencil
    template <typename StencilArray>
    AMREX_GPU_DEVICE
    static auto ComputeCandidateFluxes(const StencilArray& stencil) {
        std::array<Set::Scalar, Traits::NumCandidateStencils> fluxes;
        
        for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
            fluxes[r] = 0.0;
            for (int p = 0; p < Traits::PointsPerStencil; ++p) {
                int point_idx = Traits::SubStencilOffsets[r][p];
                fluxes[r] += Traits::CandidateStencilCoefficients[r][p] * stencil[point_idx];
            }
        }
        
        return fluxes;
    }
    
    // Main reconstruction method
    template <typename StencilArray>
    AMREX_GPU_DEVICE
    static Set::Scalar ReconstructValue(const StencilArray& stencil) {
        // 1. Compute candidate fluxes for each candidate stencil
        auto candidate_fluxes = ComputeCandidateFluxes(stencil);
        
        // 2. Compute smoothness indicators
        auto indicators = ComputeSmoothnessIndicators(stencil);
        
        // 3. Compute nonlinear weights (using variant-specific method)
        auto weights = Traits::ComputeNonlinearWeights(indicators);
        
        // 4. Form convex combination for final reconstruction
        Set::Scalar reconstructed_value = 0.0;
        for (int r = 0; r < Traits::NumCandidateStencils; ++r) {
            reconstructed_value += weights[r] * candidate_fluxes[r];
        }
        
        return reconstructed_value;
    }
};

//==============================================================================
// WENO Reconstruction class (implements FluxReconstruction interface)
//==============================================================================
template<typename IntegratorDerivedClass, typename VariantTraits>
class WENOReconstruction : public FluxReconstruction<IntegratorDerivedClass> {
public:
    using Traits = VariantTraits;
    using Policy = WENOReconstructionPolicy<Traits>;
    
    WENOReconstruction() {}

    
    void ReconstructFlux(
        int direction, 
        int lev, 
        IntegratorDerivedClass* solver,
        amrex::MultiFab& working_buffer,
        amrex::MultiFab& QL_stencil, 
        amrex::MultiFab& QR_stencil
    ) const override {

        const int nghosts = solver->number_of_ghost_cells;

        // Process each grid patch
        for (amrex::MFIter mfi(working_buffer, false); mfi.isValid(); ++mfi) {
            const amrex::Box& face_bx_with_ghosts = mfi.grownnodaltilebox(direction, nghosts);
            
            auto const& in_arr = working_buffer.const_array(mfi);
            auto const& ql_arr = QL_stencil.array(mfi);
            auto const& qr_arr = QR_stencil.array(mfi);
            
            const amrex::Box& domain = solver->geom[lev].Domain();
            
            amrex::ParallelFor(face_bx_with_ghosts, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {
                // Construct index arrays with correct number of components
                int index[3];
                int left_index[3];
                int right_index[3];
                int lower_bounds[3];
                int upper_bounds[3];
                
                // Set up indices and bounds
                index[0] = iface; index[1] = jface; index[2] = kface;
                left_index[0] = iface; left_index[1] = jface; left_index[2] = kface;
                right_index[0] = iface; right_index[1] = jface; right_index[2] = kface;

                lower_bounds[0] = face_bx_with_ghosts.smallEnd(0); 
                lower_bounds[1] = face_bx_with_ghosts.smallEnd(1); 
                lower_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.smallEnd(2) : 0;
                upper_bounds[0] = face_bx_with_ghosts.bigEnd(0);
                upper_bounds[1] = face_bx_with_ghosts.bigEnd(1);
                upper_bounds[2] = AMREX_SPACEDIM == 3 ? face_bx_with_ghosts.bigEnd(2) : 0;
                
                // Clamp the main index
                ClampIndices(index, lower_bounds, upper_bounds);
                
                // Shift and clamp for left and right indices
                int left_offset[3] = {0, 0, 0};
                int right_offset[3] = {0, 0, 0};
                left_offset[direction] = 0;  // Shift left
                right_offset[direction] = 1;  // Shift right
    
                ShiftAndClampIndices(left_index, left_offset, lower_bounds, upper_bounds, direction);
                ShiftAndClampIndices(right_index, right_offset, lower_bounds, upper_bounds, direction);

                // Process each component
                for (int n = 0; n < solver->number_of_components; ++n) {
                    // Extract stencil values based on direction
                    std::array<Set::Scalar, Traits::TotalPoints> stencil;
                    
                    // Stencil center point
                    int half_points = Traits::TotalPoints / 2;
                    
                    // Extract stencil in the specified direction
                    for (int s = 0; s < Traits::TotalPoints; ++s) {
                        // Create stencil index (offsetting in the correct direction)
                        int stencil_index[3] = {left_index[0], left_index[1], left_index[2]};
                        stencil_index[direction] = left_index[direction] - half_points + s;
                        
                        // Clamp stencil index to domain boundaries
                        ClampIndices(stencil_index, lower_bounds, upper_bounds);
                        
                        // Extract value
                        stencil[s] = in_arr(stencil_index[0], stencil_index[1], stencil_index[2], n);
                    }
                    
                    // Reconstruct left state (QL) at interface
                    ql_arr(left_index[0], left_index[1], left_index[2], n) = Policy::ReconstructValue(stencil);
                    
                    // Only set QR if we're not at the domain boundary
                    qr_arr(right_index[0], right_index[1], right_index[2], n) = ql_arr(left_index[0], left_index[1], left_index[2], n);

                }
            });
        }
    }    
        
};

//==============================================================================
// Convenience type aliases for common WENO variants
//==============================================================================
template<typename IntegratorDerivedClass>
using WENOJS5 = WENOReconstruction<IntegratorDerivedClass, WENOJS5Traits>;

template<typename IntegratorDerivedClass>
using WENOZ5 = WENOReconstruction<IntegratorDerivedClass, WENOZ5Traits>;

//==============================================================================
// Template for creating new WENO variants easily
//==============================================================================

// Example: How to create a new WENO variant (WENO-M5)
/*
struct WENOM5Traits : public WENOJS5Traits {
    // Override only what's different - in this case, the nonlinear weight calculation
    template <typename IndicatorArray>
    AMREX_GPU_DEVICE
    static auto ComputeNonlinearWeights(const IndicatorArray& indicators) {
        // WENO-M specific implementation
        // ...
    }
};

template<typename IntegratorDerivedClass>
using WENOM5 = WENOReconstruction<IntegratorDerivedClass, WENOM5Traits>;
*/

} // namespace Numeric

#endif // NUMERIC_WENO_RECONSTRUCTION_H

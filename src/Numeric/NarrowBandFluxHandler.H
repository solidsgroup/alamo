#ifndef NUMERIC_NARROWBANDFLUXHANDLER_H
#define NUMERIC_NARROWBANDFLUXHANDLER_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Set/Set.H"
#include "Integrator/NarrowBandLevelset.H"
#include "FluxHandler.H"
#include "Util/Util.H"
#include "Numeric/Stencil.H"                    
#include <memory>
#include <cmath>
#include <algorithm>

// Include foward declare for NarrowBandLevelset.H
namespace Integrator {
    class NarrowBandLevelset;
}

namespace Numeric {
// Base class for flux reconstruction
template <>
class FluxReconstruction<Integrator::NarrowBandLevelset> {
public:
    virtual ~FluxReconstruction() = default;

    /**
    * @brief Perform flux reconstruction for the specified direction.
    */
    virtual void ReconstructFlux(int direction, const amrex::iMultiFab& mask,
        const amrex::MultiFab& narrowband_LS, amrex::MultiFab& Flux_mf,
        const Set::Scalar* DX) const = 0; 
};

// NarrowBand specific flux class
template <>
class FluxMethod<Integrator::NarrowBandLevelset> {
public:
    virtual ~FluxMethod() = default;

    /**
    * @brief Compute the flux for the specified direction.
    */
    virtual void ComputeFlux(int direction, const amrex::iMultiFab& mask, 
        const amrex::MultiFab& narrowband_LS, amrex::MultiFab& Flux_mf) const = 0;
};

// FluxHandler class for Integrator::NarrowBandLevelset to correctly assess distribution 
// mapping for narrowband
template <>
class FluxHandler<Integrator::NarrowBandLevelset> {
public:
    enum Directions { Xdir = 0, Ydir = 1, Zdir = 2 };

    FluxHandler() {}

    void SetReconstruction(std::shared_ptr<FluxReconstruction<Integrator::NarrowBandLevelset>> recon) {
        reconstruction = std::move(recon);
    }

    void SetFluxMethod(std::shared_ptr<FluxMethod<Integrator::NarrowBandLevelset>> method) {
        flux_method = std::move(method);
    }

    void ConstructFluxes(int lev, Integrator::NarrowBandLevelset* solver) const {
        if (!reconstruction || !flux_method) {
            amrex::Abort("FluxHandler: Reconstruction or flux method is not set.");
        }

        // Define LSData structure
        auto& LSData = solver->LSData;
        const int num_objs = LSData.num_objs;

        // Loop through objects
        for (int obj=0; obj<num_objs; obj++) {
            // Get object structure
            auto& object = LSData.objects[lev][obj];
            const auto& DX = object.FlowData->geom.CellSize();

            // Define object geometry
            const auto& Tube_ba = object.Tube_ba;
            const auto& Tube_dm = object.Tube_dm;
            const int ngrow     = object.FlowData->GhostCells();

            // Define (i)MultiFabs
            const auto& LS_mf = *object.TubeLS_mf;
            const auto& Tube_imf = *object.TubeMask_imf;
            auto& Flux_mf = *object.TubeFlux_mf;
            Flux_mf.setVal(0.0);

            // Loop through directions
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                // Perform reconstruction/flux computation
                reconstruction->ReconstructFlux(dir, Tube_imf, LS_mf, Flux_mf, DX);
                flux_method->ComputeFlux(dir, Tube_imf, LS_mf, Flux_mf);
            }

            // Fill Ghosts
            Flux_mf.FillBoundary();
        }
    }

private:
    std::shared_ptr<FluxReconstruction<Integrator::NarrowBandLevelset>> reconstruction;  // Reconstruction method
    std::shared_ptr<FluxMethod<Integrator::NarrowBandLevelset>> flux_method;             // Flux method
};

template <typename Props>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void GatherStencil(Set::MultiVector& stencil,
                    const amrex::Array4<const amrex::Real>& LS_arr,
                    const amrex::Array4<const int>& Tube_arr,
                    int i, int j, int k,
                    const amrex::IntVect debug,
                    int direction, int SignV, 
                    int LScomp, int Tubecomp,
                    const int OUTSIDETUBE,
                    const Set::Scalar* DX)
{
    for (int s = 0; s < Props::TotalPoints; ++s) {
        // Define offsets
        Set::MultiMatrix indices = Props::GetStencilIndices();

        // Get the current stencil point using the offsets matrix
        int idx[3] = { i, j, k };
        idx[direction] += indices(s, SignV);
        amrex::IntVect i_plus(AMREX_D_DECL(idx[0], idx[1], idx[2]));

        // Get the negative stencil point for difference
        amrex::IntVect i_minus = i_plus;
        i_minus[direction] -= 1;

        // Perform flux computation
        const bool in_tube = Tube_arr(i_plus,Tubecomp) != OUTSIDETUBE &&
                            Tube_arr(i_minus,Tubecomp) != OUTSIDETUBE;
        const auto& LS_plus = LS_arr(i_plus,LScomp);
        const auto& LS_minus = LS_arr(i_minus,LScomp);
        
        stencil[s] = in_tube ? (LS_plus - LS_minus) / DX[direction] : 0.0;
        /*if (amrex::IntVect(AMREX_D_DECL(i,j,k)) == debug) {
            amrex::Print() << "Debugging cell " << debug << " in direction " << direction << std::endl;
            amrex::Print() << "i_plus: " << i_plus << ", i_minus: " <<i_minus << std::endl;
            amrex::Print() << "in_tube: " << in_tube << std::endl;
            amrex::Print() << "LS_plus: " << LS_plus << ", LS_minus: " << LS_minus << std::endl;
            amrex::Print() << "flux: " << stencil[s] << std::endl;
        }*/
    }
}

template <typename Props>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Set::Scalar ComputeSmoothness(const Set::MultiVector& stencil,
                                const Set::MultiMatrix& diff_ops,
                                const std::array<Set::Scalar, Props::NumDifferenceTerms>& weights)
{
    Set::Scalar beta = 0.0;
    for (int t = 0; t < Props::NumDifferenceTerms; ++t) {
        Set::Scalar val = diff_ops.row(t).dot(stencil);
        beta += weights[t] * val * val;
    }
    return beta;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Set::Scalar ComputeCandidateFlux(const Set::MultiVector& stencil,
                                const Set::MultiVector& coeffs)
{
    return coeffs.dot(stencil);
}

template <typename Props>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void ComputeWeights(const std::array<Set::Scalar, Props::NumCandidateStencils>& betas,
                    const Set::MultiVector& optimal,
                    std::array<Set::Scalar, Props::NumCandidateStencils>& weights)
{
    std::array<Set::Scalar, Props::NumCandidateStencils> alphas;
    Set::Scalar sum_alpha = 0.0;
    for (int s = 0; s < Props::NumCandidateStencils; ++s) {
        alphas[s] = optimal(s) / std::pow(betas[s] + Props::Epsilon, Props::Power);
        sum_alpha += alphas[s];
    }
    for (int s = 0; s < Props::NumCandidateStencils; ++s) {
        weights[s] = alphas[s] / sum_alpha;
    }
}

template <typename Props>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Set::Scalar ApplyReconstruction(const std::array<Set::Scalar, Props::NumCandidateStencils>& weights,
                                const std::array<Set::Scalar, Props::NumCandidateStencils>& fluxes)
{
    Set::Scalar result = 0.0;
    for (int s = 0; s < Props::NumCandidateStencils; ++s)
        result += weights[s] * fluxes[s];
    return result;
}

//==============================================================================
// FirstOrder Traits - 1st Order Scheme
//==============================================================================

/**
* Implements a 1st order reconstruction scheme.
*/
struct FirstOrderProps {
    static constexpr int OrderValue = 1;
    static constexpr int GhostCells = 1;
    static constexpr Set::Scalar Epsilon = 1e-6;  // Small parameter to avoid division by zero
    static constexpr int NumCandidateStencils = 1;  // Number of candidate stencils
    static constexpr int PointsPerStencil = 1;      // Points per candidate stencil
    static constexpr int TotalPoints = 1;           // Total points in the full stencil
    static constexpr int NumDifferenceTerms = 1;    // Number of difference terms in smoothness indicator
    static constexpr int Power = 1;                 // Power parameter in weights formula

    // Optimal (linear) weights for WENO-JS5 as a MultiVector
    static Set::MultiVector GetOptimalWeights() {
        Set::MultiVector weights(NumCandidateStencils);
        weights << 1.0;
        return weights;
    }

    // Coefficients for candidate stencils reconstruction as a MultiMatrix
    static Set::MultiMatrix GetCandidateStencilCoefficients() {
        Set::MultiMatrix coeffs(NumCandidateStencils, PointsPerStencil);
        coeffs << 1.0;    // Right stencil S₂ coefficients
        return coeffs;
    }

    // Stencil Indices too use in flux computations
    static Set::MultiMatrix GetStencilIndices() {
        Set::MultiMatrix Indices(TotalPoints, 2);
        Indices << 0, 1;
        return Indices;
    }

    // Stencil offsets define which points to use for each candidate stencil
    static Set::MultiMatrix GetSubStencilOffsets() {
        Set::MultiMatrix offsets(NumCandidateStencils, PointsPerStencil);
        offsets << 0;
        return offsets;
    }

    /**
    * @struct FirstOrderDifferenceOperator
    * @brief Difference operators for First Order smoothness indicators
    */
    struct DifferenceOperatorNB {

        // Define term weights for the smoothness indicator calculation
        static constexpr std::array<Set::Scalar, NumDifferenceTerms> TermWeights = {
            1.0
        };

        // Function to get coefficient matrices
        static std::vector<Set::MultiMatrix> GetCoefficients() {
            std::vector<Set::MultiMatrix> coeffs(NumCandidateStencils);

            // Stencil 0 coefficients
            coeffs[0] = (Set::MultiMatrix(1, 1) << 1.0).finished();

            return coeffs;
        }
    };
};

//==============================================================================
// WENO-JS5 Traits - 5th order Jiang-Shu WENO scheme
//==============================================================================

/**
* Implements the 5th order WENO scheme from Jiang & Shu (1996).
* This scheme uses 3 candidate stencils with 3 points each.
*/
struct WENOJS5Props {
    static constexpr int OrderValue = 5;
    static constexpr int GhostCells = 3;
    static constexpr Set::Scalar Epsilon = 1e-6;  // Small parameter to avoid division by zero
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

    // Stencil Indices to use in flux computations
    static Set::MultiMatrix GetStencilIndices() {
        Set::MultiMatrix Indices(TotalPoints, 2);
        Indices << -2, 3,
                    -1, 2,
                    0, 1,
                    1, 0,
                    2, -1;
        return Indices;
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
    struct DifferenceOperatorNB {

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
};

template <typename Props, typename IntegratorT>
class GenericFluxReconstruction : public FluxReconstruction<IntegratorT> {
public:
    void ReconstructFlux(int direction,
                         const amrex::iMultiFab& Tube,
                         const amrex::MultiFab& LS,
                         amrex::MultiFab& Flux,
                         const Set::Scalar* DX) const override
    {
        constexpr int nghosts = Props::GhostCells; 
        constexpr int OUTSIDETUBE = NarrowBandTubeType::OUTSIDETUBE;

        const LevelSetData::ObjectMeta Data;
        const int LScomp   = Data.Comp(LevelSetData::ObjectMeta::TubeFieldComp::LS);
        const int Tubecomp = Data.Comp(LevelSetData::ObjectMeta::TubeIntComp::TUBE_MASK);
        constexpr int Velcomp[AMREX_SPACEDIM] = {
            Data.Comp(LevelSetData::ObjectMeta::TubeFieldComp::VELOCITY_X)
#if AMREX_SPACEDIM >= 2
            , Data.Comp(LevelSetData::ObjectMeta::TubeFieldComp::VELOCITY_Y)
#endif
#if AMREX_SPACEDIM == 3
            , Data.Comp(LevelSetData::ObjectMeta::TubeFieldComp::VELOCITY_Z)
#endif
        };

        const auto subOffsets      = Props::GetSubStencilOffsets();
        const auto candidateCoeffs = Props::GetCandidateStencilCoefficients();
        const auto smoothCoeffs    = Props::DifferenceOperatorNB::GetCoefficients();
        const auto termWeights     = Props::DifferenceOperatorNB::TermWeights;
        const auto optimalWeights  = Props::GetOptimalWeights();

        const amrex::IntVect debug(AMREX_D_DECL(37,66,0));

        for (amrex::MFIter mfi(Flux, false); mfi.isValid(); ++mfi) {
            const auto& ghost_bx = mfi.growntilebox(nghosts);
            const auto& valid_bx = mfi.validbox();
        
            auto const& LS_arr   = LS.const_array(mfi);
            auto const& Tube_arr = Tube.const_array(mfi);
            auto const& Flux_arr = Flux.array(mfi);

            amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if (!valid_bx.contains(i,j,k)) return;
                if (Tube_arr(i,j,k,Tubecomp) == OUTSIDETUBE) return;

                const Set::Scalar V = LS_arr(i,j,k,Velcomp[direction]);
                const int SignV = (V > 0) ? 0 : 1;

                Set::MultiVector stencil(Props::TotalPoints);
                GatherStencil<Props>(stencil, LS_arr, Tube_arr, i, j, k,
                                     debug, direction, SignV, LScomp, Tubecomp,
                                     OUTSIDETUBE, DX);

                std::array<Set::Scalar, Props::NumCandidateStencils> beta;
                std::array<Set::Scalar, Props::NumCandidateStencils> flux;
                std::array<Set::Scalar, Props::NumCandidateStencils> weights;

                for (int s = 0; s < Props::NumCandidateStencils; ++s) {
                    Set::MultiVector vals(Props::PointsPerStencil);
                    for (int p = 0; p < Props::PointsPerStencil; ++p) {
                        int stencil_idx = subOffsets(s, p) + Props::TotalPoints / 2;
                        vals(p) = stencil[stencil_idx];
                    }
                    beta[s] = ComputeSmoothness<Props>(vals, smoothCoeffs[s], termWeights);
                    flux[s] = ComputeCandidateFlux(vals, candidateCoeffs.row(s));
                }

                ComputeWeights<Props>(beta, optimalWeights, weights);
                Flux_arr(i,j,k,direction) = V * ApplyReconstruction<Props>(weights, flux);
            });
        }
    }
};

// First-order specific to NarrowBandLevelset
template <>
class FirstOrderReconstruction<Integrator::NarrowBandLevelset>
    : public GenericFluxReconstruction<FirstOrderProps, Integrator::NarrowBandLevelset> {};

// Fifth-order WENO JS specific to NarrowBandLevelset
template <>
class FifthOrderWENOReconstruction<Integrator::NarrowBandLevelset>
    : public GenericFluxReconstruction<WENOJS5Props, Integrator::NarrowBandLevelset> {};


// WILL DELETE THIS
template <>
class NarrowbandFluxComputation<Integrator::NarrowBandLevelset> : public FluxMethod<Integrator::NarrowBandLevelset> {
public:
    void ComputeFlux(int direction,
        const amrex::iMultiFab& Tube,
        const amrex::MultiFab& LS, 
        amrex::MultiFab& Flux) const override {
        const int nghosts = 1; 

        // Define OutsideNarrowband
        const int OutsideTube = NarrowBandTubeType::OUTSIDETUBE;

        // Define components
        const LevelSetData::ObjectMeta Data;
        const int mask_comp = Data.Comp(LevelSetData::ObjectMeta::TubeIntComp::TUBE_MASK);
        constexpr int vel_comp[AMREX_SPACEDIM] = {
            Data.Comp(LevelSetData::ObjectMeta::TubeFieldComp::VELOCITY_X)
#if AMREX_SPACEDIM >= 2
            , Data.Comp(LevelSetData::ObjectMeta::TubeFieldComp::VELOCITY_Y)
#endif
#if AMREX_SPACEDIM == 3
            , Data.Comp(LevelSetData::ObjectMeta::TubeFieldComp::VELOCITY_Z)
#endif
        };

        const amrex::IntVect debug(AMREX_D_DECL(60,64,0));
    }
};

}  // namespace Numeric

#endif  // NUMERIC_NARROWBANDFLUXHANDLER_H

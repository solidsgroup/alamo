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
        const amrex::MultiFab& narrowband_LS, amrex::MultiFab& QL_stencil, 
        amrex::MultiFab& QR_stencil) const = 0; 
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
        const amrex::MultiFab& narrowband_LS, const amrex::MultiFab& QL_stencil, const amrex::MultiFab& QR_stencil, 
        amrex::MultiFab& flux_mf) const = 0;
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

            // Define object geometry
            const auto& Tube_ba     = object.Tube_ba;
            const auto& Tube_dm     = object.Tube_dm;
            const int ngrow         = object.FlowData->GhostCells();

            // Define(i)MultiFabs
            const auto& LS_mf = *object.TubeLS_mf;
            const auto& Tube_imf = *object.TubeMask_imf;

            // Loop through directions
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                // Create temporary multifabs to store levelset reconstruction
                amrex::BoxArray Face_ba = amrex::convert(Tube_ba, amrex::IntVect::TheDimensionVector(dir));
                amrex::MultiFab QL(Face_ba, Tube_dm, 1, ngrow);
                amrex::MultiFab QR(Face_ba, Tube_dm, 1, ngrow);

                // Select proper flux multifab based on direction
                amrex::MultiFab& flux_mf = 
                    [&, dir]() -> amrex::MultiFab& {
                        switch (dir) {
                            case 0: return *object.TubeXFlux_mf;
                #if AMREX_SPACEDIM >= 2
                            case 1: return *object.TubeYFlux_mf;
                #endif
                #if AMREX_SPACEDIM == 3
                            case 2: return *object.TubeZFlux_mf;
                #endif
                            default: amrex::Abort("Invalid direction"); return *object.TubeXFlux_mf; // dummy return for compiler
                        }
                    }();

                    // Set flux to 0
                    flux_mf.setVal(0.0);

                    // Perform reconstruction/flux computation
                    reconstruction->ReconstructFlux(dir, Tube_imf, LS_mf, QL, QR);
                    flux_method->ComputeFlux(dir, Tube_imf, LS_mf, QL, QR, flux_mf);
            }
        }

    }

private:
    std::shared_ptr<FluxReconstruction<Integrator::NarrowBandLevelset>> reconstruction;  // Reconstruction method
    std::shared_ptr<FluxMethod<Integrator::NarrowBandLevelset>> flux_method;             // Flux method
};

// First-order specific to NarrowBandLevelset integrator class
template <>
class FirstOrderReconstruction<Integrator::NarrowBandLevelset> : public FluxReconstruction<Integrator::NarrowBandLevelset> {
public:
    void ReconstructFlux(int direction,
                        const amrex::iMultiFab& mask,
                        const amrex::MultiFab& narrowband_LS,
                        amrex::MultiFab& QL_stencil,
                        amrex::MultiFab& QR_stencil) const override {

        constexpr int nghosts = 1;
        constexpr int OutsideTube = NarrowBandTubeType::OUTSIDETUBE;

        const LevelSetData::ObjectMeta Data;
        const int LScomp    = Data.Comp(LevelSetData::ObjectMeta::TubeFieldComp::LS_OLD);
        const int mask_comp = Data.Comp(LevelSetData::ObjectMeta::TubeIntComp::TUBE_MASK);

        for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
            const amrex::Box& face_bx = mfi.grownnodaltilebox(direction, nghosts);

            int lo[3], hi[3];
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                lo[d] = face_bx.smallEnd(d);
                hi[d] = face_bx.bigEnd(d) - 1;
            }

            auto const& LS_arr   = narrowband_LS.const_array(mfi);
            auto const& mask_arr = mask.const_array(mfi);
            auto const& QL_arr   = QL_stencil.array(mfi);
            auto const& QR_arr   = QR_stencil.array(mfi);

            amrex::ParallelFor(face_bx, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {
                // Clamp main face index
                int idxFace[3] = {iface, jface, kface};
                ClampIndices(idxFace, lo, hi);
                const amrex::IntVect face_iv(AMREX_D_DECL(idxFace[0], idxFace[1], idxFace[2]));

                // Cell-centered: left and center
                int idxC[3] = {iface, jface, kface};
                int idxL[3] = {iface, jface, kface};
                int idxL_offsets[3] = {0, 0, 0};
                idxL_offsets[direction] -= 1;

                ClampIndices(idxC, lo, hi);
                ShiftAndClampIndices(idxL, idxL_offsets, lo, hi, direction);

                const amrex::IntVect ivC(AMREX_D_DECL(idxC[0], idxC[1], idxC[2]));
                const amrex::IntVect ivL(AMREX_D_DECL(idxL[0], idxL[1], idxL[2]));

                if (mask_arr(ivC, mask_comp) == OutsideTube && mask_arr(ivL, mask_comp) == OutsideTube)
                    return;

                QL_arr(face_iv, 0) = LS_arr(ivL, LScomp);
                QR_arr(face_iv, 0) = LS_arr(ivC, LScomp);
            });
        }

        QL_stencil.FillBoundary();
        QR_stencil.FillBoundary();
    }
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

    // Offsets for the left stencil (for QL at i+1/2): i-3 to i+1 - cell-centered offsets!
    static constexpr std::array<int, TotalPoints> LeftStencilOffsets = {-3, -2, -1, 0, +1};

    // Offsets for the right stencil (for QR at i+1/2): i-2 to i+2 - cell-centered offsets!
    static constexpr std::array<int, TotalPoints> RightStencilOffsets = {-2, -1, 0, +1, +2};

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
    struct WENOJS5DifferenceOperatorNB {

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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void GatherStencil(Set::MultiVector& stencil,
                    const amrex::Array4<const amrex::Real>& LS_arr,
                    const int* idxC,
                    int dir, int LScomp,
                    const int* lo, const int* hi,
                    const std::array<int, WENOJS5Props::TotalPoints>& offsets)
{
    for (int s = 0; s < WENOJS5Props::TotalPoints; ++s) {
        int idx[3] = {idxC[0], idxC[1], idxC[2]};
        int idx_offset[3] = {0, 0, 0};
        idx_offset[dir] += offsets[s];
        ShiftAndClampIndices(idx, idx_offset, lo, hi, dir);
        stencil[s] = LS_arr(amrex::IntVect(AMREX_D_DECL(idx[0], idx[1], idx[2])), LScomp);
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Set::Scalar ComputeSmoothness(const Set::MultiVector& stencil,
                                const Set::MultiMatrix& diff_ops,
                                const std::array<Set::Scalar, WENOJS5Props::NumDifferenceTerms>& weights)
{
    Set::Scalar beta = 0.0;
    for (int t = 0; t < WENOJS5Props::NumDifferenceTerms; ++t) {
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void ComputeWeights(const std::array<Set::Scalar, WENOJS5Props::NumCandidateStencils>& betas,
                    const Set::MultiVector& optimal,
                    std::array<Set::Scalar, WENOJS5Props::NumCandidateStencils>& weights)
{
    std::array<Set::Scalar, WENOJS5Props::NumCandidateStencils> alphas;
    Set::Scalar sum_alpha = 0.0;
    for (int s = 0; s < WENOJS5Props::NumCandidateStencils; ++s) {
        alphas[s] = optimal(s) / std::pow(betas[s] + WENOJS5Props::Epsilon, WENOJS5Props::Power);
        sum_alpha += alphas[s];
    }
    for (int s = 0; s < WENOJS5Props::NumCandidateStencils; ++s) {
        weights[s] = alphas[s] / sum_alpha;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Set::Scalar ApplyReconstruction(const std::array<Set::Scalar, WENOJS5Props::NumCandidateStencils>& weights,
                                const std::array<Set::Scalar, WENOJS5Props::NumCandidateStencils>& fluxes)
{
    Set::Scalar result = 0.0;
    for (int s = 0; s < WENOJS5Props::NumCandidateStencils; ++s)
        result += weights[s] * fluxes[s];
    return result;
}

// Fifth-order WENO JS specific to NarrowBandLevelset integrator class
template <>
class FifthOrderWENOReconstruction<Integrator::NarrowBandLevelset>
    : public FluxReconstruction<Integrator::NarrowBandLevelset> {
public:
    void ReconstructFlux(int direction,
                        const amrex::iMultiFab& mask,
                        const amrex::MultiFab& narrowband_LS,
                        amrex::MultiFab& QL_stencil,
                        amrex::MultiFab& QR_stencil) const override
    {
        const int nghosts = 3;
        const int OutsideTube = NarrowBandTubeType::OUTSIDETUBE;

        const LevelSetData::ObjectMeta Data;
        const int LScomp    = Data.Comp(LevelSetData::ObjectMeta::TubeFieldComp::LS_OLD);
        const int mask_comp = Data.Comp(LevelSetData::ObjectMeta::TubeIntComp::TUBE_MASK);

        const auto offsetsL = WENOJS5Props::LeftStencilOffsets;
        const auto offsetsR = WENOJS5Props::RightStencilOffsets;
        const auto subOffsets = WENOJS5Props::GetSubStencilOffsets();
        const auto candidateCoeffs = WENOJS5Props::GetCandidateStencilCoefficients();
        const auto smoothCoeffs = WENOJS5Props::WENOJS5DifferenceOperatorNB::GetCoefficients();
        const auto termWeights = WENOJS5Props::WENOJS5DifferenceOperatorNB::TermWeights;
        const auto optimalWeights = WENOJS5Props::GetOptimalWeights();

        for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
            const amrex::Box& face_bx = mfi.grownnodaltilebox(direction, nghosts);

            int lo[3], hi[3];
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                lo[d] = face_bx.smallEnd(d);
                hi[d] = face_bx.bigEnd(d) - 1;
            }

            auto const& LS_arr   = narrowband_LS.const_array(mfi);
            auto const& mask_arr = mask.const_array(mfi);
            auto const& QL_arr   = QL_stencil.array(mfi);
            auto const& QR_arr   = QR_stencil.array(mfi);

            amrex::ParallelFor(face_bx, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {
                // Clamp main face index
                int idxFace[3] = {iface, jface, kface};
                ClampIndices(idxFace, lo, hi);
                const amrex::IntVect face_iv(AMREX_D_DECL(idxFace[0], idxFace[1], idxFace[2]));

                // Cell-centered: left and center
                int idxC[3] = {iface, jface, kface};
                int idxL[3] = {iface, jface, kface};
                int idxL_offsets[3] = {0, 0, 0};
                idxL_offsets[direction] -= 1;

                ClampIndices(idxC, lo, hi);
                ShiftAndClampIndices(idxL, idxL_offsets, lo, hi, direction);

                const amrex::IntVect ivC(AMREX_D_DECL(idxC[0], idxC[1], idxC[2]));
                const amrex::IntVect ivL(AMREX_D_DECL(idxL[0], idxL[1], idxL[2]));

                if (mask_arr(ivC, mask_comp) == OutsideTube && mask_arr(ivL, mask_comp) == OutsideTube)
                    return;

                Set::MultiVector stencilL(WENOJS5Props::TotalPoints);
                Set::MultiVector stencilR(WENOJS5Props::TotalPoints);

                GatherStencil(stencilL, LS_arr, idxC, direction, LScomp, lo, hi, offsetsL);
                GatherStencil(stencilR, LS_arr, idxC, direction, LScomp, lo, hi, offsetsR);

                std::array<Set::Scalar, WENOJS5Props::NumCandidateStencils> betaL, betaR;
                std::array<Set::Scalar, WENOJS5Props::NumCandidateStencils> fluxL, fluxR;
                std::array<Set::Scalar, WENOJS5Props::NumCandidateStencils> weightsL, weightsR;

                for (int s = 0; s < WENOJS5Props::NumCandidateStencils; ++s) {
                    Set::MultiVector valsL(WENOJS5Props::PointsPerStencil);
                    Set::MultiVector valsR(WENOJS5Props::PointsPerStencil);

                    for (int p = 0; p < WENOJS5Props::PointsPerStencil; ++p) {
                        int stencil_idx = subOffsets(s, p) + WENOJS5Props::TotalPoints / 2;
                        valsL(p) = stencilL[stencil_idx];
                        valsR(p) = stencilR[stencil_idx];
                    }

                    betaL[s] = ComputeSmoothness(valsL, smoothCoeffs[s], termWeights);
                    betaR[s] = ComputeSmoothness(valsR, smoothCoeffs[s], termWeights);
                    fluxL[s] = ComputeCandidateFlux(valsL, candidateCoeffs.row(s));
                    fluxR[s] = ComputeCandidateFlux(valsR, candidateCoeffs.row(s));
                }

                ComputeWeights(betaL, optimalWeights, weightsL);
                ComputeWeights(betaR, optimalWeights, weightsR);

                QL_arr(face_iv, 0) = ApplyReconstruction(weightsL, fluxL);
                QR_arr(face_iv, 0) = ApplyReconstruction(weightsR, fluxR);
            });
        }

        QL_stencil.FillBoundary();
        QR_stencil.FillBoundary();
    }
};

template <>
class NarrowbandFluxComputation<Integrator::NarrowBandLevelset> : public FluxMethod<Integrator::NarrowBandLevelset> {
public:
    void ComputeFlux(int direction,
        const amrex::iMultiFab& mask,
        const amrex::MultiFab& narrowband_LS, 
        const amrex::MultiFab& QL_stencil, 
        const amrex::MultiFab& QR_stencil, 
        amrex::MultiFab& flux_mf) const override {
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

        // Loop through face-centered cells
        for (amrex::MFIter mfi(QL_stencil, false); mfi.isValid(); ++mfi) {
            // Use a face-centered box to compute fluxes
            const amrex::Box& face_bx = mfi.grownnodaltilebox(direction, nghosts);

            int lo[3], hi[3];
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                lo[d] = face_bx.smallEnd(d);
                hi[d] = face_bx.bigEnd(d) - 1;
            }

            // Define multifab arrays
            auto const& mask_arr = mask.const_array(mfi);
            auto const& LS_arr = narrowband_LS.const_array(mfi);
            auto const& QL_arr = QL_stencil.const_array(mfi);
            auto const& QR_arr = QR_stencil.const_array(mfi);
            auto const& flux_arr = flux_mf.array(mfi);

            // Loop through face indices
            amrex::ParallelFor(face_bx, [=] AMREX_GPU_DEVICE(int iface, int jface, int kface) noexcept {
                // Clamp main face index
                int idxFace[3] = {iface, jface, kface};
                ClampIndices(idxFace, lo, hi);
                const amrex::IntVect face_iv(AMREX_D_DECL(idxFace[0], idxFace[1], idxFace[2]));

                // Cell-centered: left and center
                int idxC[3] = {iface, jface, kface};
                int idxL[3] = {iface, jface, kface};
                int idxL_offsets[3] = {0, 0, 0};
                idxL_offsets[direction] -= 1;

                ClampIndices(idxC, lo, hi);
                ShiftAndClampIndices(idxL, idxL_offsets, lo, hi, direction);

                const amrex::IntVect ivC(AMREX_D_DECL(idxC[0], idxC[1], idxC[2]));
                const amrex::IntVect ivL(AMREX_D_DECL(idxL[0], idxL[1], idxL[2]));

                if (mask_arr(ivC, mask_comp) == OutsideTube && mask_arr(ivL, mask_comp) == OutsideTube)
                    return;

                // Get velocity of cell center
                Set::Scalar V = LS_arr(ivC, vel_comp[direction]);

                // Fetch reconstructed level-set values
                Set::Scalar QL = QL_arr(face_iv, 0);
                Set::Scalar QR = QR_arr(face_iv, 0);

                // Get proper flux based on velocity
                Set::Scalar flux = (V > 0) ? V * QL : V * QR;

                // Store computed flux
                flux_arr(face_iv, 0) = flux;
            });
        }

        // Fill ghost cells
        flux_mf.FillBoundary();
    }
};

}  // namespace Numeric

#endif  // NUMERIC_NARROWBANDFLUXHANDLER_H

#ifndef INTEGRATOR_NARROWBANDLEVELSET_H // Include guards
#define INTEGRATOR_NARROWBANDLEVELSET_H // 

// Standard library includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <vector>
#include <stdexcept>
#include <memory>

// Alamo includes
#include "IO/ParmParse.H"
#include "Integrator/Integrator.H"
#include "BC/BC.H"
#include "BC/Nothing.H"
#include "IC/IC.H"

#include "Numeric/Stencil.H"

namespace Numeric {
    template <typename T> class FluxReconstruction;
    template <typename T> class FluxMethod;
    template <typename T> class FluxHandler;  
    template <typename T> class FirstOrderReconstruction;
    template <typename T> class FifthOrderWENOReconstruction;
    template <typename T> class NarrowbandFluxComputation;
    //template <typename T> class LocalLaxFriedrichsMethod;
    template <typename T> class TimeStepper;
    template <typename T> class TimeSteppingScheme;
    template <typename T> class EulerForwardScheme;  
    template <typename T> class RK3Scheme;
}

// Define a structure for the values of the narrowband
struct NarrowBandTubeType
{
    static constexpr int INTERFACE   = 0;           // 0 level set
    static constexpr int INNERBAND   = 1;           // |phi| < 4.0 * min_DX
    static constexpr int INSIDETUBE  = 2;           // |phi| < 6.0 * min_DX
    static constexpr int EDGEPOINT   = 3;           // |phi| >= 6.0 * min_DX & nbr <= INSIDETUBE
    static constexpr int OUTSIDETUBE = 4;           // |phi| >= 6.0 * min_DX
};

struct FlagType
{
    static constexpr int NEVER   = -1;
    static constexpr int NO      =  0;
    static constexpr int YES     =  1;
};

// Quick structure for single connectivity (2, 4, 6) to loop through neighbor values
struct Neighbors {
    static constexpr int num_neighbors =
    #if AMREX_SPACEDIM == 1
        2;
    #elif AMREX_SPACEDIM == 2
        4;
    #elif AMREX_SPACEDIM == 3
        6;
    #endif

    static constexpr int offsets[num_neighbors][3] =
    #if AMREX_SPACEDIM == 1
    {
        {-1, 0, 0},
        { 1, 0, 0}
    };
    #elif AMREX_SPACEDIM == 2
    {
        {-1, 0, 0},
        { 1, 0, 0},
        { 0,-1, 0},
        { 0, 1, 0}
    };
    #elif AMREX_SPACEDIM == 3
    {
        {-1, 0, 0},
        { 1, 0, 0},
        { 0,-1, 0},
        { 0, 1, 0},
        { 0, 0,-1},
        { 0, 0, 1}
    };
    #endif
};

struct LevelSetData {
    
    // Flux Reconstruction variables and required ghosts
    enum class FluxReconstruction {
        FirstOrder,
        ThirdOrderENO,
        FifthOrderWENO,
    };

    static constexpr int RequiredGhosts(FluxReconstruction scheme) {
        using FR = FluxReconstruction;
        switch (scheme) {
            case FR::FirstOrder:     return 1;
            case FR::ThirdOrderENO:  return 2;
            case FR::FifthOrderWENO: return 3;
            default:                 return 3;  // fallback
        }
    }

    int GhostCells() const {
        return RequiredGhosts(reconstruction_method);
    }

    // FluxReconstruction 
    FluxReconstruction reconstruction_method = FluxReconstruction::FifthOrderWENO;

    // Flow domain information
    amrex::Geometry geom;
    amrex::Box flow_domain;

    // Helper function to return flow grid spacing for stencil derivatives
    Set::Scalar min_DX() const {
        const Set::Scalar* DX = geom.CellSize(); 
        return *std::min_element(DX, DX + AMREX_SPACEDIM);
    }
 
    // Band width information
    Set::Scalar Inner_Tube_width;
    Set::Scalar Tube_width;

    // Helper function to compute band and inner band widths
    void ComputeTubeWidths() {
        Inner_Tube_width =  4.0 * min_DX();
        Tube_width = 6.0 * min_DX();
    }

    // Structure to store object data. Assumes each object has its own levelset field
    struct ObjectMeta {
        // Pointer to parent for physical_domain/ghost cell count
        LevelSetData* FlowData = nullptr;  

        // Metadata
        int object_id = 0; 
        int ls_id = -1;           
        int material_id = 0;
        bool has_narrowband = false;
        bool reinitialize = false;

        // Narrowband geometry information
        amrex::Box Tube_domain;
        amrex::BoxArray Tube_ba;
        amrex::DistributionMapping Tube_dm;

        // Field Data
        std::unique_ptr<amrex::MultiFab>  TubeLS_mf;
        std::unique_ptr<amrex::iMultiFab> TubeMask_imf;

        std::unique_ptr<amrex::MultiFab>  TubeXFlux_mf;
    #if AMREX_SPACEDIM >= 2
        std::unique_ptr<amrex::MultiFab>  TubeYFlux_mf;
    #endif
    #if AMREX_SPACEDIM == 3
        std::unique_ptr<amrex::MultiFab>  TubeZFlux_mf;
    #endif

        // Enum classes to store components for LS and mask (i)MultiFabs
        enum class TubeFieldComp : int {
            LS = 0,
            LS_OLD,
            VELOCITY_X,
    #if AMREX_SPACEDIM >= 2
            VELOCITY_Y,
    #endif
    #if AMREX_SPACEDIM == 3
            VELOCITY_Z,
    #endif
            NORMAL_X,
    #if AMREX_SPACEDIM >= 2
            NORMAL_Y,
    #endif
    #if AMREX_SPACEDIM == 3
            NORMAL_Z,
    #endif
            CURVATURE,
            NUM_REAL_COMPONENTS
        };

        enum class TubeIntComp : int {
            TUBE_MASK = 0,
            ZERO_LS,
            CPT_ID,
            NUM_INT_COMPONENTS
        };

        // Helper to get integer of components in each Tube (i)mf
        template <typename CompType>
        static constexpr int Comp(CompType c) {
            return static_cast<int>(c);
        }

        // Structures to hold LS Boundary 
        // Enum for user-configurable BC types
        enum class FieldBCType {
            Neumann,
            Dirichlet,
            Periodic,
            Reflective,
            ShockPulse
        };

        struct VelocityBCConfig {
            FieldBCType bc_type = FieldBCType::Neumann;

            // Dirichlet velocity values per component
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dirichlet_vals = {0.0};

            // Add other shock velocity parameters here
        };

        // User-configurable BC settings
        struct TubeBCConfig {
            FieldBCType normal_bc = FieldBCType::Neumann;
            VelocityBCConfig velocity_bc;
        };

        // ObjectMeta TubeBCConfig variable
        TubeBCConfig bc_config;

        // Manual BC custom functions
        // Template Neumann function
        template <typename MFType>
        void ApplyNeumannBCs(MFType& mf, int comp)
        {
            static_assert(std::is_same<MFType, amrex::MultiFab>::value ||
                        std::is_same<MFType, amrex::iMultiFab>::value,
                        "ApplyNeumannBCs only supports MultiFab or iMultiFab");

            // Get flow domain info
            const auto& flow_domain = FlowData->flow_domain;
            const auto& flow_lo = flow_domain.smallEnd();
            const auto& flow_hi = flow_domain.bigEnd();

            for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const amrex::Box& full = mfi.growntilebox();
                const auto& fab = mf.array(mfi);

                amrex::ParallelFor(full, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    if (!flow_domain.contains(i, j, k)) {
                        int ii = amrex::min(amrex::max(i, flow_lo[0]), flow_hi[0]);
        #if AMREX_SPACEDIM >= 2
                        int jj = amrex::min(amrex::max(j, flow_lo[1]), flow_hi[1]);
        #else
                        int jj = j;
        #endif
        #if AMREX_SPACEDIM == 3
                        int kk = amrex::min(amrex::max(k, flow_lo[2]), flow_hi[2]);
        #else
                        int kk = k;
        #endif
                        fab(i,j,k,comp) = fab(ii,jj,kk,comp);
                    }
                });
            }
        }

        void ApplyLSBCs() {
            using TFC = TubeFieldComp;
            AMREX_ALWAYS_ASSERT(TubeLS_mf != nullptr);

            for (auto comp : {TFC::LS, TFC::LS_OLD}) {
                int c = Comp(comp);

                ApplyNeumannBCs(*TubeLS_mf, c);
            }
        }

        void ApplyVelBCs(VelocityBCConfig vel_bc) {
            using TFC = TubeFieldComp;
            AMREX_ALWAYS_ASSERT(TubeLS_mf != nullptr);

            std::vector<int> comps = {
                Comp(TFC::VELOCITY_X)
        #if AMREX_SPACEDIM >= 2
                , Comp(TFC::VELOCITY_Y)
        #endif
        #if AMREX_SPACEDIM == 3
                , Comp(TFC::VELOCITY_Z)
        #endif
            };

            const FieldBCType bc_type = vel_bc.bc_type;

            // Get flow domain info
            const auto& flow_domain = FlowData->flow_domain;
            const auto& flow_lo = flow_domain.smallEnd();
            const auto& flow_hi = flow_domain.bigEnd();

            for (amrex::MFIter mfi(*TubeLS_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const auto& full  = mfi.growntilebox();
                auto fab = TubeLS_mf->array(mfi);

                amrex::ParallelFor(full, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    if (!flow_domain.contains(i,j,k)) {
                        int ii = amrex::min(amrex::max(i, flow_lo[0]), flow_hi[0]);
        #if AMREX_SPACEDIM >= 2
                        int jj = amrex::min(amrex::max(j, flow_lo[1]), flow_hi[1]);
        #else
                        int jj = j;
        #endif
        #if AMREX_SPACEDIM == 3
                        int kk = amrex::min(amrex::max(k, flow_lo[2]), flow_hi[2]);
        #else
                        int kk = k;
        #endif

                        for (std::size_t d = 0; d < comps.size(); ++d) {
                            int comp = comps[d];

                            if (bc_type == FieldBCType::Neumann) {
                                fab(i,j,k,comp) = fab(ii,jj,kk,comp);
                            } else if (bc_type == FieldBCType::Reflective) {
                                fab(i,j,k,comp) = -fab(ii,jj,kk,comp);
                            } else if (bc_type == FieldBCType::Dirichlet) {
                                fab(i,j,k,comp) = vel_bc.dirichlet_vals[d];
                            } 
                        }
                    }
                });
            }
        }

        void ApplyNormBCs(FieldBCType bc_type) {
            using TFC = TubeFieldComp;
            AMREX_ALWAYS_ASSERT(TubeLS_mf != nullptr);

            std::vector<int> comps;
            comps.push_back(Comp(TFC::NORMAL_X));
        #if AMREX_SPACEDIM >= 2
            comps.push_back(Comp(TFC::NORMAL_Y));
        #endif
        #if AMREX_SPACEDIM == 3
            comps.push_back(Comp(TFC::NORMAL_Z));
        #endif

            // Get flow domain info
            const auto& flow_domain = FlowData->flow_domain;
            const auto& flow_lo = flow_domain.smallEnd();
            const auto& flow_hi = flow_domain.bigEnd();

            for (amrex::MFIter mfi(*TubeLS_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const auto& full  = mfi.growntilebox();
                auto fab = TubeLS_mf->array(mfi);

                amrex::ParallelFor(full, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    if (!flow_domain.contains(i,j,k)) {
                        int ii = amrex::min(amrex::max(i, flow_lo[0]), flow_hi[0]);
        #if AMREX_SPACEDIM >= 2
                        int jj = amrex::min(amrex::max(j, flow_lo[1]), flow_hi[1]);
        #else
                        int jj = j;
        #endif
        #if AMREX_SPACEDIM == 3
                        int kk = amrex::min(amrex::max(k, flow_lo[2]), flow_hi[2]);
        #else
                        int kk = k;
        #endif

                        // Mirror all
                        for (int idx : comps) {
                            fab(i,j,k,idx) = fab(ii,jj,kk,idx);
                        }

                        if (bc_type == FieldBCType::Reflective) {
                            if (i < flow_lo[0] || i > flow_hi[0]) fab(i,j,k,Comp(TFC::NORMAL_X)) *= -1;
        #if AMREX_SPACEDIM >= 2
                            if (j < flow_lo[1] || j > flow_hi[1]) fab(i,j,k,Comp(TFC::NORMAL_Y)) *= -1;
        #endif
        #if AMREX_SPACEDIM == 3
                            if (k < flow_lo[2] || k > flow_hi[2]) fab(i,j,k,Comp(TFC::NORMAL_Z)) *= -1;
        #endif
                        }
                    }
                });
            }
        }

        void ApplyCurveBCs() {
            using TFC = TubeFieldComp;
            AMREX_ALWAYS_ASSERT(TubeLS_mf != nullptr);

            int comp = Comp(TFC::CURVATURE);

            ApplyNeumannBCs(*TubeLS_mf, comp);
        }

        void ApplyAllFieldBCs() {
            ApplyLSBCs();
            ApplyNormBCs(bc_config.normal_bc);
            ApplyVelBCs(bc_config.velocity_bc);
            ApplyCurveBCs();
        }

        void ApplyMaskBCs() {
            using TIC = TubeIntComp;
            AMREX_ALWAYS_ASSERT(TubeMask_imf != nullptr);

            int comp = Comp(TIC::TUBE_MASK);

            ApplyNeumannBCs(*TubeMask_imf, comp);
        }

        void ApplyCPTBCs() {
            using TIC = TubeIntComp;
            AMREX_ALWAYS_ASSERT(TubeMask_imf != nullptr);

            int comp = Comp(TIC::CPT_ID);

            ApplyNeumannBCs(*TubeMask_imf, comp);
        }

        void ApplyAllIntBCs() {
            ApplyMaskBCs();
            ApplyCPTBCs();
        }

        template <typename FabType>
        void InitLayout(const FabType& mask, int comp, int tile_size)
        {
            amrex::BoxArray ba(Tube_domain);
            ba.maxSize(tile_size);

            std::vector<amrex::Box> kept_boxes;

            using value_type = typename FabType::value_type;

            // Define OUTSIDETUBE based on FabType
            value_type OUTSIDETUBE;
            if constexpr (std::is_same_v<value_type, int>) {
                OUTSIDETUBE = NarrowBandTubeType::OUTSIDETUBE;
            } else {
                OUTSIDETUBE = static_cast<value_type>(NarrowBandTubeType::OUTSIDETUBE);
            }

            for (int i = 0; i < ba.size(); ++i)
            {
                const amrex::Box& tile = ba[i];
                bool has_band = false;

                for (amrex::MFIter mfi(mask, false); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& valid_bx = mfi.validbox();
                    const amrex::Box overlap = tile & valid_bx;
                    if (!overlap.ok()) continue;

                    auto const& arr = mask.const_array(mfi);

                    for (amrex::IntVect iv = overlap.smallEnd(); overlap.contains(iv); overlap.next(iv))
                    {
                        if (arr(iv, comp) != OUTSIDETUBE)
                        {
                            has_band = true;
                            break;
                        }
                    }

                    if (has_band) break;
                }

                if (has_band)
                {
                    kept_boxes.push_back(tile);
                }
            }

            amrex::BoxList bl;
            for (const auto& box : kept_boxes) {
                bl.push_back(box);
            }

            Tube_ba.define(bl);
            Tube_dm = amrex::DistributionMapping(Tube_ba);
        }

        template <typename FabType>
        void AllocateFabs(const FabType& mask, int comp, int tile_size = 32)
        {
            // Build Tube_ba and Tube_dm by checking where mask is not OUTSIDETUBE
            InitLayout(mask, comp, tile_size);

            const int n_real_comps = Comp(TubeFieldComp::NUM_REAL_COMPONENTS);
            const int n_int_comps  = Comp(TubeIntComp::NUM_INT_COMPONENTS);
            const int ngrow = FlowData->GhostCells();

            // Always allocate real field (MultiFab)
            TubeLS_mf = std::make_unique<amrex::MultiFab>(Tube_ba, Tube_dm, n_real_comps, ngrow);
            TubeLS_mf->setVal(0.0);

            // Always allocate integer mask (iMultiFab)
            TubeMask_imf = std::make_unique<amrex::iMultiFab>(Tube_ba, Tube_dm, n_int_comps, ngrow);

            TubeMask_imf->setVal(NarrowBandTubeType::OUTSIDETUBE, Comp(TubeIntComp::TUBE_MASK), 1, ngrow);
            TubeMask_imf->setVal(0, Comp(TubeIntComp::ZERO_LS), 1, ngrow);
            TubeMask_imf->setVal(0, Comp(TubeIntComp::CPT_ID), 1, ngrow);

            // Allocate fluxes (always real-valued)
            amrex::BoxArray xflux_ba = amrex::convert(Tube_ba, amrex::IntVect::TheDimensionVector(0));
            TubeXFlux_mf = std::make_unique<amrex::MultiFab>(xflux_ba, Tube_dm, 1, ngrow);
            TubeXFlux_mf->setVal(0.0);

        #if AMREX_SPACEDIM >= 2
            amrex::BoxArray yflux_ba = amrex::convert(Tube_ba, amrex::IntVect::TheDimensionVector(1));
            TubeYFlux_mf = std::make_unique<amrex::MultiFab>(yflux_ba, Tube_dm, 1, ngrow);
            TubeYFlux_mf->setVal(0.0);
        #endif

        #if AMREX_SPACEDIM == 3
            amrex::BoxArray zflux_ba = amrex::convert(Tube_ba, amrex::IntVect::TheDimensionVector(2));
            TubeZFlux_mf = std::make_unique<amrex::MultiFab>(zflux_ba, Tube_dm, 1, ngrow);
            TubeZFlux_mf->setVal(0.0);
        #endif
        }

        void SwapLSandLSold()
        {
            const int LScomp = Comp(TubeFieldComp::LS);
            const int LSoldcomp = Comp(TubeFieldComp::LS_OLD);

            for (amrex::MFIter mfi(*TubeLS_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& ghost_bx = mfi.growntilebox(FlowData->GhostCells());
                auto arr = TubeLS_mf->array(mfi);

                amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    using std::swap;
                    swap(arr(i,j,k,LScomp), arr(i,j,k,LSoldcomp));
                });
            }
        }

        // Deallocate memory
        void DeallocateFabs() {
            if (TubeLS_mf) TubeLS_mf.reset();
            if (TubeMask_imf) TubeMask_imf.reset();
            if (TubeXFlux_mf) TubeXFlux_mf.reset();

        #if AMREX_SPACEDIM >= 2
            if (TubeYFlux_mf) TubeYFlux_mf.reset();
        #endif

        #if AMREX_SPACEDIM == 3
            if (TubeZFlux_mf) TubeZFlux_mf.reset();
        #endif
        }

        // Parallel copy flow domain variables to narrowband real fields
        void LoadLSFields(const amrex::MultiFab& Flow_LS, const amrex::MultiFab& Flow_LSold)
        {
            const int ngrow = FlowData->GhostCells();

            const int LScomp    = Comp(TubeFieldComp::LS);
            const int LSoldcomp = Comp(TubeFieldComp::LS_OLD);

            // Parallel copy LS
            TubeLS_mf->ParallelCopy(Flow_LS, ls_id, LScomp, 1, ngrow, ngrow);
            TubeLS_mf->ParallelCopy(Flow_LSold, ls_id, LSoldcomp, 1, ngrow, ngrow);

            // Apply boundary conditions
            ApplyLSBCs();
        }

        // Load flow domain masks to Tube iMultiFab
        void LoadIntFields(const amrex::MultiFab& Tubeflow,
                            const amrex::MultiFab& Zeroflow,
                            const amrex::MultiFab& CPTflow)
        {
            // Aliases to target component indices
            const int Tubecomp = Comp(TubeIntComp::TUBE_MASK);
            const int Zerocomp = Comp(TubeIntComp::ZERO_LS);
            const int CPTcomp  = Comp(TubeIntComp::CPT_ID);
            const int num_comp = Comp(TubeIntComp::NUM_INT_COMPONENTS);
            const int ngrow    = FlowData->GhostCells();

            // Define temp multifab to store amrex::real data
            amrex::MultiFab tmp_Tube_mf(Tube_ba, Tube_dm, num_comp, ngrow);

            // Parallel copy data from flow MultiFabs
            tmp_Tube_mf.ParallelCopy(Tubeflow, ls_id, Tubecomp, 1, ngrow, ngrow);
            tmp_Tube_mf.ParallelCopy(Zeroflow, ls_id, Zerocomp, 1, ngrow, ngrow);
            tmp_Tube_mf.ParallelCopy(CPTflow, ls_id, CPTcomp, 1, ngrow, ngrow);

            // Fill ghosts
            tmp_Tube_mf.FillBoundary();

            // Loop through and copy 
            for (amrex::MFIter mfi(*TubeMask_imf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const amrex::Box& ghost_bx = mfi.growntilebox(ngrow);

                auto const& src = tmp_Tube_mf.array(mfi);
                auto const& dst = (*TubeMask_imf)[mfi].array();

                amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    for (int n=0; n<num_comp; n++){
                        dst(i, j, k, n) = static_cast<int>(src(i, j, k, n));
                    }
                });
            }

            // Fill ghosts
            ApplyAllIntBCs();
        }

        void LoadVelField(const amrex::MultiFab& Flow_vel)
        {
            constexpr int vel_comp[AMREX_SPACEDIM] = {
                Comp(TubeFieldComp::VELOCITY_X)
        #if AMREX_SPACEDIM >= 2
                , Comp(TubeFieldComp::VELOCITY_Y)
        #endif
        #if AMREX_SPACEDIM == 3
                , Comp(TubeFieldComp::VELOCITY_Z)
        #endif
            };

            const int ngrow = FlowData->GhostCells();

            // Copy
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                TubeLS_mf->ParallelCopy(Flow_vel, d, vel_comp[d], 1, ngrow, ngrow);
            }

            // ApplyBCs
            ApplyVelBCs(bc_config.velocity_bc);
        }

    #if AMREX_SPACEDIM == 1
        void CopyGeometryToFlowFields(amrex::MultiFab& NormalXflow) const
        {
            const int Xcomp = Comp(TubeFieldComp::NORMAL_X);

            NormalXflow.ParallelCopy(*TubeLS_mf, Xcomp, ls_id, 1);

            NormalXflow.FillBoundary(ls_id, 1);
        }
    #endif

    #if AMREX_SPACEDIM == 2
        void CopyGeometryToFlowFields(amrex::MultiFab& NormalXflow,
                                        amrex::MultiFab& NormalYflow) const
        {
            const int Xcomp = Comp(TubeFieldComp::NORMAL_X);
            const int Ycomp = Comp(TubeFieldComp::NORMAL_Y);

            NormalXflow.ParallelCopy(*TubeLS_mf, Xcomp, ls_id, 1);
            NormalYflow.ParallelCopy(*TubeLS_mf, Ycomp, ls_id, 1);

            NormalXflow.FillBoundary(ls_id, 1);
            NormalYflow.FillBoundary(ls_id, 1);
        }
    #endif

    #if AMREX_SPACEDIM == 3
        void CopyGeometryToFlowFields(amrex::MultiFab& NormalXflow,
                                        amrex::MultiFab& NormalYflow,
                                        amrex::MultiFab& NormalZflow) const
        {
            const int Xcomp = Comp(TubeFieldComp::NORMAL_X);
            const int Ycomp = Comp(TubeFieldComp::NORMAL_Y);
            const int Zcomp = Comp(TubeFieldComp::NORMAL_Z);

            NormalXflow.ParallelCopy(*TubeLS_mf, Xcomp, ls_id, 1);
            NormalYflow.ParallelCopy(*TubeLS_mf, Ycomp, ls_id, 1);
            NormalZflow.ParallelCopy(*TubeLS_mf, Zcomp, ls_id, 1);

            NormalXflow.FillBoundary(ls_id, 1);
            NormalYflow.FillBoundary(ls_id, 1);
            NormalZflow.FillBoundary(ls_id, 1);
        }
    #endif

        void ZeroTheTube(amrex::MultiFab& Tubeflow, amrex::MultiFab& Zeroflow) {
            const Set::Scalar EDGEPOINT  = NarrowBandTubeType::EDGEPOINT;
            const Set::Scalar INSIDETUBE = NarrowBandTubeType::INSIDETUBE;

            const int ngrow = FlowData->GhostCells();

            for (amrex::MFIter mfi(Tubeflow, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // Determine the overlapping region between flow_domain and Tube_domain
                const auto& flow_ghost_bx = mfi.growntilebox(ngrow);
                const amrex::Box overlap = Tube_domain.grow(ngrow) & flow_ghost_bx;
                if (overlap.isEmpty()) continue;

                const auto& Tube_arr = Tubeflow.array(mfi);
                const auto& Zero_arr = Zeroflow.array(mfi);

                amrex::ParallelFor(overlap, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    if (Tube_arr(i,j,k,ls_id) == EDGEPOINT) Tube_arr(i,j,k,ls_id) = INSIDETUBE;
                    if (Zero_arr(i,j,k,ls_id) == object_id) Zero_arr(i,j,k,ls_id) = 0;
                });
            }

            // Save ghosts
            Tubeflow.FillBoundary(ls_id, 1);
            Zeroflow.FillBoundary(ls_id, 1);
        }

        void UpdateNarrowbandTubeandMapping(amrex::MultiFab& LSflow, 
            const amrex::MultiFab& LSoldflow,
            amrex::MultiFab& Tubeflow,
            amrex::MultiFab& Zeroflow,
            amrex::MultiFab& CPTflow)
        {
            // Constants
            const Set::Scalar INTERFACE   = NarrowBandTubeType::INTERFACE;
            const Set::Scalar INNERBAND   = NarrowBandTubeType::INNERBAND;
            const Set::Scalar INSIDETUBE  = NarrowBandTubeType::INSIDETUBE;
            const Set::Scalar EDGEPOINT   = NarrowBandTubeType::EDGEPOINT;
            const Set::Scalar OUTSIDETUBE = NarrowBandTubeType::OUTSIDETUBE;

            const auto& physical_domain = FlowData->flow_domain;
            const Set::Scalar INNERTUBE = FlowData->Inner_Tube_width;
            const Set::Scalar OUTERTUBE = FlowData->Tube_width;
            const int ngrow = FlowData->GhostCells();

            // Set the EDGEPOINT cells to OUTSIDETUBE
            ZeroTheTube(Tubeflow, Zeroflow);

            // Track tile-wise min/max for bounding box
            amrex::IntVect global_min(amrex::IntVect::TheMaxVector());
            amrex::IntVect global_max(amrex::IntVect::TheMinVector());

            // Operate directly on flow domain using domain from previous timestep
            for (amrex::MFIter mfi(LSflow, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // Determine the overlapping region between flow_domain and Tube_domain
                const auto& flow_ghost_bx = mfi.growntilebox(ngrow);
                const amrex::Box overlap = Tube_domain.grow(ngrow) & flow_ghost_bx;
                if (overlap.isEmpty()) continue;

                // Get array views
                const auto& LS_arr = LSflow.array(mfi);
                const auto& Tube_arr = Tubeflow.array(mfi);
                const auto& Zero_arr = Zeroflow.array(mfi);
                const auto& CPT_arr  = CPTflow.array(mfi);

                // Local tile-wise bounds
                amrex::IntVect local_min(amrex::IntVect::TheMaxVector());
                amrex::IntVect local_max(amrex::IntVect::TheMinVector());

                // ------------------------------------------------------------------
                // Pass 1: Compute initial tube classification
                // ------------------------------------------------------------------
                amrex::ParallelFor(overlap, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Skip cells outside the physical domain
                    if (!physical_domain.contains(i,j,k)) return;

                    // Skip cells that were classified as OUTSIDETUBE the past timestep
                    Set::Scalar& Tubeval = Tube_arr(i,j,k,ls_id);
                    if (Tubeval == OUTSIDETUBE) return;

                    // Begin as OUTSIDETUBE
                    Tubeval = OUTSIDETUBE;

                    // Clamp LS and get abs value
                    auto& LS = LS_arr(i,j,k,ls_id);
                    LS = std::clamp(LS, -OUTERTUBE, OUTERTUBE);
                    Set::Scalar absLS = std::abs(LS);

                    // Get CPT
                    CPT_arr(i,j,k,ls_id) = (LS <= 0.0) ? object_id : 0;

                    // Check if cell is in the tube
                    if (absLS < OUTERTUBE) {
                        // Mark as INSIDETUBE
                        Tubeval = INSIDETUBE;

                        // If cell is in the outer half of tube, it's done
                        if (absLS >= INNERTUBE) return;

                        // Mark as INNERBAND
                        Tubeval = INNERBAND;

                        // Loop through all of the neighbor directions
                        for (int d = 0; d < Neighbors::num_neighbors; ++d) {
                            int ni = i + Neighbors::offsets[d][0];
                            int nj = j + Neighbors::offsets[d][1];
                            int nk = k + Neighbors::offsets[d][2];

                            // Skip nbr cells outside ghost range
                            if (!overlap.contains(ni,nj,nk)) continue;

                            // Skip nbr cells outside domain boundary
                            if (!physical_domain.contains(ni,nj,nk)) continue; 

                            // Skip if neighbor is not in the LS tube
                            if (Tube_arr(ni,nj,nk,ls_id) == OUTSIDETUBE) continue;

                            // Get LS value of neighbor in current direction
                            const Set::Scalar LSnbr = LS_arr(ni,nj,nk,ls_id);

                            // Check if neighbor is on the opposite side of interface
                            if (LS * LSnbr > 0.0 || Zero_arr(i,j,k,ls_id) != 0) continue;
                            
                            // Mark INTERFACE cells
                            Tubeval = INTERFACE;
                            Zero_arr(i,j,k,ls_id) = object_id;
                        }
                    }
                });

                // ----------------------------------
                // Pass 2: Determine EDGEPOINTs
                // ----------------------------------
                amrex::ParallelFor(overlap, [=, &local_min, &local_max] AMREX_GPU_DEVICE(int i, int j, int k) {
                    // Skip cells outside domain
                    if (!physical_domain.contains(i,j,k)) return;

                    // Skip cells not classified as INSIDETUBE
                    const auto& Tubeval = Tube_arr(i,j,k,ls_id);
                    if (Tubeval != INSIDETUBE) return;

                    // Loop through nbrs and assign EDGEPOINTS
                    for (int d = 0; d < Neighbors::num_neighbors; ++d) {
                        int ni = i + Neighbors::offsets[d][0];
                        int nj = j + Neighbors::offsets[d][1];
                        int nk = k + Neighbors::offsets[d][2];

                        // Skip cells outside ghost range
                        if (!overlap.contains(ni,nj,nk)) continue;

                        // Skip nbrs outside domain boundary
                        if(!physical_domain.contains(ni,nj,nk)) continue;

                        if (Tube_arr(ni,nj,nk,ls_id) == OUTSIDETUBE) {
                            Tube_arr(ni,nj,nk,ls_id) = EDGEPOINT;
                        }
                    }
                });

                // ----------------------------------
                // Pass 3: Bounding box for EDGEPOINTs
                // ----------------------------------
                amrex::ParallelFor(overlap, [=, &local_min, &local_max] AMREX_GPU_DEVICE(int i, int j, int k) {
                    if (Tube_arr(i,j,k,ls_id) == EDGEPOINT) {
                        AMREX_D_TERM(
                            amrex::Gpu::Atomic::Min(&local_min[0], i);
                            amrex::Gpu::Atomic::Max(&local_max[0], i);,
                            amrex::Gpu::Atomic::Min(&local_min[1], j);
                            amrex::Gpu::Atomic::Max(&local_max[1], j);,
                            amrex::Gpu::Atomic::Min(&local_min[2], k);
                            amrex::Gpu::Atomic::Max(&local_max[2], k);
                        )
                    }
                });

                // Sync tile min/max across CPU
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    global_min[d] = std::min(global_min[d], local_min[d]);
                    global_max[d] = std::max(global_max[d], local_max[d]);
                }
            }

            // Fill ghost cells
            LSflow.FillBoundary(ls_id, 1);
            Tubeflow.FillBoundary(ls_id, 1);
            Zeroflow.FillBoundary(ls_id, 1);
            CPTflow.FillBoundary(ls_id, 1);

            // Check if narrowband exists
            has_narrowband = (global_min <= global_max);

            // Domain update and buffer swap
            if (has_narrowband)
            {
                // Update Tube_domain using new coordinates
                Tube_domain = amrex::Box(global_min, global_max) & physical_domain;

                // Copy velocities to tmp velocity field
                constexpr int vel_comp[AMREX_SPACEDIM] = {
                    Comp(TubeFieldComp::VELOCITY_X)
        #if AMREX_SPACEDIM >= 2
                    , Comp(TubeFieldComp::VELOCITY_Y)
        #endif
        #if AMREX_SPACEDIM == 3
                    , Comp(TubeFieldComp::VELOCITY_Z)
        #endif
                };

                amrex::MultiFab tmpvel(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), AMREX_SPACEDIM, ngrow);

                // Copy
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    tmpvel.ParallelCopy(*TubeLS_mf, vel_comp[d], d, 1, ngrow, ngrow);
                }

                // Reallocate pointers
                AllocateFabs(Tubeflow, ls_id);

                // Copy back to pointers
                LoadLSFields(LSflow, LSoldflow);
                LoadIntFields(Tubeflow, Zeroflow, CPTflow);

                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    TubeLS_mf->ParallelCopy(tmpvel, d, vel_comp[d], 1, ngrow, ngrow);
                }
                ApplyVelBCs(bc_config.velocity_bc);
            }
        }

        void ComputeGeometryQuantities() 
        {
            // Define DX for gradient
            const amrex::Real* DX = FlowData->geom.CellSize();

            // Define OutsideNarrowband constant
            const int OUTSIDETUBE = NarrowBandTubeType::OUTSIDETUBE;

            // Get number of ghosts
            const int ngrow = FlowData->GhostCells();

            // Define component values
            const int LScomp =  Comp(TubeFieldComp::LS);
            const int Tubecomp = Comp(TubeIntComp::TUBE_MASK);

            // Define normal components and set to 0
            const int normalXcomp = Comp(TubeFieldComp::NORMAL_X);
        #if AMREX_SPACEDIM >= 2
            const int normalYcomp = Comp(TubeFieldComp::NORMAL_Y);
        #endif
        #if AMREX_SPACEDIM == 3
            const int normalZcomp = Comp(TubeFieldComp::NORMAL_Z);
        #endif

            // Create array mapping dimension to normal component
            constexpr int norm_comp[AMREX_SPACEDIM] = {
                normalXcomp
        #if AMREX_SPACEDIM >= 2
                , normalYcomp
        #endif
        #if AMREX_SPACEDIM == 3
                , normalZccomp
        #endif
            };

            const int Curvecomp = Comp(TubeFieldComp::CURVATURE);

            // Iterate over only the narrowband cells
            for (amrex::MFIter mfi(*TubeLS_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const amrex::Box& ghost_bx = mfi.growntilebox(1);

                const auto& LS_arr = TubeLS_mf->array(mfi);
                const auto& Tube_arr = TubeMask_imf->const_array(mfi);
                
                amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    // Skip non-narrowband cells
                    amrex::IntVect coord(AMREX_D_DECL(i,j,k));
                    if (Tube_arr(i, j, k, Tubecomp) == OUTSIDETUBE) return; 

                    // Utilize Numeric::Gradient for simplified computations
                    std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil = Numeric::GetStencil(i,j,k,ghost_bx);
                    Set::Vector gradient = Numeric::Gradient(LS_arr, i, j, k, LScomp, DX, stencil);

                    Set::Scalar grad_mag = gradient.norm();

                    for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
                        LS_arr(i, j, k, norm_comp[dim]) = (grad_mag > 0.0) ? (gradient[dim] / grad_mag) : 0.0;
                    }
                });
            }  

            // Fill ghost cells
            ApplyNormBCs(bc_config.normal_bc);

            /*for (amrex::MFIter mfi(*LS_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const amrex::Box& valid_bx = mfi.validbox();

                const auto& LS_arr = LS_mf->array(mfi);
                const auto& Tube_arr = Tube_imf->const_array(mfi);

                amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    // skip if not within narowband
                    amrex::IntVect coord(AMREX_D_DECL(i,j,k));
                    if (Tube_arr(coord, maskcomp) == OUTSIDETUBE) return; 
                    
                    // Get curvature as divergence of normal
                    std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil = Numeric::GetStencil(i, j, k, physical_domain);
                    LS_arr(coord, curvecomp) = Numeric::Divergence(LS_arr, i, j, k, normalXcomp, DX, stencil);
                });
            }*/
        }

        // Assuming first-order gradients. Will change later to higher order
        void ReinitializeLS(amrex::MultiFab& LSflow, Set::Scalar cflReinit)
        {
            // Define flag constants
            const int NEVER   = FlagType::NEVER;
            const int NO      = FlagType::NO;
            const int YES     = FlagType::YES;

            // Define (i)MultiFab components
            const int LScomp = Comp(TubeFieldComp::LS);

            const int Tubecomp = Comp(TubeIntComp::TUBE_MASK);
            const int Zerocomp = Comp(TubeIntComp::ZERO_LS);
            const int CPTcomp  = Comp(TubeIntComp::CPT_ID);

            // Define Tube constants
            const int OUTSIDETUBE = NarrowBandTubeType::OUTSIDETUBE;
            const Set::Scalar INNERTUBE = FlowData->Inner_Tube_width;
            const Set::Scalar OUTERTUBE = FlowData->Tube_width;

            // Define flow constants
            const auto& physical_domain = FlowData->flow_domain;
            const int ngrow = FlowData->GhostCells();
            const amrex::Real* DX = FlowData->geom.CellSize();
            const amrex::Dim3 lo = amrex::lbound(physical_domain);
            const amrex::Dim3 hi = amrex::ubound(physical_domain);

            // Allocate and initialize temporary multifabs
            amrex::MultiFab LSreinit(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, ngrow);
            amrex::MultiFab LSreinitold(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, ngrow);
            amrex::iMultiFab compflag(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, ngrow);

            LSreinit.setVal(10.0 * OUTERTUBE);
            compflag.setVal(NO);

            // Set LSreinit to narrowband values and find INTERFACE cells
            for (amrex::MFIter mfi(*TubeLS_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // Determine the overlapping region between flow_domain and Tube_domain
                const auto& ghost_bx = mfi.growntilebox(ngrow);

                // Get array views
                const auto& LS_arr = TubeLS_mf->const_array(mfi);
                const auto& Tube_arr = TubeMask_imf->const_array(mfi);
                const auto& reinit_arr = LSreinit.array(mfi);
                const auto& flag_arr = compflag.array(mfi);

                amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Skip cells outside flow domain
                    if(!physical_domain.contains(i,j,k)) return;

                    // Skip non-narrowband cells
                    if (Tube_arr(i,j,k,Tubecomp) == OUTSIDETUBE) return;

                    // Set the narrowband value in the LSreinit MultiFab
                    Set::Scalar LS = LS_arr(i,j,k,LScomp);
                    Set::Scalar signedLS = LS;
                    if (Tube_arr(i,j,k,CPTcomp) == object_id) signedLS = -std::abs(LS);
                    reinit_arr(i,j,k) = signedLS;

                    // Set the comp flag for interface cells to NEVER
                    if (Tube_arr(i,j,k,Zerocomp) != 0) {
                        flag_arr(i,j,k) = NEVER;

                        // Check nbrs
                        for (int d = 0; d < Neighbors::num_neighbors; ++d) {
                            int ni = i + Neighbors::offsets[d][0];
                            int nj = j + Neighbors::offsets[d][1];
                            int nk = k + Neighbors::offsets[d][2];

                            // Skip if outside ghost box
                            if (!ghost_bx.contains(ni,nj,nk)) continue;

                            // Skip if nbr is outside domain
                            if (!physical_domain.contains(ni,nj,nk)) continue;

                            // Check for sign change and that both LS values are within tube
                            Set::Scalar LSnbr = LS_arr(ni,nj,nk,LScomp);

                            if (LS * LSnbr <= 0.0 && std::max(LS, LSnbr) < INNERTUBE) {
                                flag_arr(ni,nj,nk) = NEVER;
                            }
                        }
                    }
                });
            }

            // Fill ghost cells
            LSreinit.FillBoundary();
            compflag.FillBoundary();

            // Set initial expanding tube points
            for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // Determine the overlapping region between flow_domain and Tube_domain
                const auto& ghost_bx = mfi.growntilebox(ngrow);

                // Get array views
                const auto& Tube_arr = TubeMask_imf->const_array(mfi);
                const auto& reinit_arr = LSreinit.array(mfi);
                const auto& flag_arr = compflag.array(mfi);

                amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Skip non-narrowband cells
                    if (Tube_arr(i,j,k,Tubecomp) == OUTSIDETUBE) return;

                    // Skip cells not in flow domain
                    if (!physical_domain.contains(i,j,k)) return;

                    // Skip interface cells
                    if (flag_arr(i,j,k) == NEVER) return;

                    // Get current LS value
                    const Set::Scalar LS = reinit_arr(i,j,k);

                    // Loop through nbrs
                    for (int d = 0; d < Neighbors::num_neighbors; ++d) {
                        int ni = i + Neighbors::offsets[d][0];
                        int nj = j + Neighbors::offsets[d][1];
                        int nk = k + Neighbors::offsets[d][2];

                        // Skip cells outside ghost cells
                        if (!ghost_bx.contains(ni,nj,nk)) continue;

                        // Skip cells outside flow domain
                        if (!physical_domain.contains(ni,nj,nk)) continue;

                        // Skip cells if nbr flag is NEVER
                        if (flag_arr(ni,nj,nk) == NEVER) continue;

                        // Set LS sign of nbr to that of current cell
                        auto& LSnbr = reinit_arr(ni,nj,nk);
                        LSnbr = std::copysign(LSnbr, LS);
                    }

                    // Set compflag to yes
                    flag_arr(i,j,k) = YES;
                });
            }

            // Fill ghosts
            LSreinit.FillBoundary();
            compflag.FillBoundary();

            // Copy LSreinit to LSreinitold
            LSreinitold.ParallelCopy(LSreinit, 0, 0, 1, ngrow, ngrow);
            LSreinitold.FillBoundary();

            // Define variables for timestep loop
            amrex::MultiFab error(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, 0);
            error.setVal(0.0);

            const Set::Scalar toler_error_reinit = 1e-3;
            const int LSREINIT_MAX_ITERS = 100;
            const Set::Scalar DTreinit = cflReinit * FlowData->min_DX();

            Set::Scalar error_max = 2.0 * toler_error_reinit;
            int icount = 0;

            // Begin artifical timestep loop
            while (error_max > toler_error_reinit && icount <= LSREINIT_MAX_ITERS) {
                // Increase iteration count
                icount += 1;

                // Swap
                std::swap(LSreinit, LSreinitold);

                // Fill ghosts
                LSreinitold.FillBoundary();

                // Set error array to 0.0
                if (icount % 5 == 0) error.setVal(0.0);

                // Main PDE loop
                for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    // Determine the overlapping region between flow_domain and Tube_domain
                    const auto& valid_bx = mfi.validbox();

                    // Get array views
                    const auto& Tube_arr = TubeMask_imf->const_array(mfi);
                    const auto& reinit_arr = LSreinit.array(mfi);
                    const auto& reinitold_arr = LSreinitold.array(mfi);
                    const auto& flag_arr = compflag.array(mfi);
                    const auto& error_arr = error.array(mfi);

                    amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        // Skip non-narrowband cells
                        if (Tube_arr(i,j,k,Tubecomp) == OUTSIDETUBE) return;

                        // ----------------------------------
                        // We may not need the following flag variables. Keeping for now to 
                        // match Fortran code
                        // ----------------------------------
                        // Cycle if compflag is NO or NEVER
                        if (flag_arr(i,j,k) <= NO) return;
                        // ----------------------------------

                        // Define sign of ls for upwinding
                        const Set::Scalar LS = reinitold_arr(i,j,k);
                        const Set::Scalar signLS = (LS >= 0.0 ? 1.0 : -1.0);

                        // Compute first order differences
                        // x-direction
                        Set::Scalar dxm = (i > lo.x)
                            ? Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(reinitold_arr, i, j, k, 0, DX, Numeric::XLo)
                            : 0.0;

                        Set::Scalar dxp = (i < hi.x)
                            ? Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(reinitold_arr, i, j, k, 0, DX, Numeric::XHi)
                            : 0.0;

                        Set::Scalar dfx = std::max(
                            std::pow(std::max(signLS * dxm, 0.0), 2),
                            std::pow(std::min(signLS * dxp, 0.0), 2)
                        );

                        // y-direction
                        Set::Scalar dfy = 0.0;
                        #if AMREX_SPACEDIM >= 2
                        Set::Scalar dym = (j > lo.y)
                            ? Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(reinitold_arr, i, j, k, 0, DX, Numeric::YLo)
                            : 0.0;

                        Set::Scalar dyp = (j < hi.y)
                            ? Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(reinitold_arr, i, j, k, 0, DX, Numeric::YHi)
                            : 0.0;

                        dfy = std::max(
                            std::pow(std::max(signLS * dym, 0.0), 2),
                            std::pow(std::min(signLS * dyp, 0.0), 2)
                        );
                        #endif

                        // z-direction
                        Set::Scalar dfz = 0.0;
                        #if AMREX_SPACEDIM == 3
                        Set::Scalar dzm = (k > lo.z)
                            ? Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(reinitold_arr, i, j, k, 0, DX, Numeric::ZLo)
                            : 0.0;

                        Set::Scalar dzp = (k < hi.z)
                            ? Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(reinitold_arr, i, j, k, 0, DX, Numeric::ZHi)
                            : 0.0;

                        dfz = std::max(
                            std::pow(std::max(signLS * dzm, 0.0), 2),
                            std::pow(std::min(signLS * dzp, 0.0), 2)
                        ); 
                        #endif

                        // Compute gradient and update
                        Set::Scalar nablaG = std::sqrt(dfx + dfy + dfz);
                        reinit_arr(i,j,k) = LS - DTreinit * signLS * (nablaG - 1.0);
                        reinit_arr(i,j,k) = std::copysign(reinit_arr(i,j,k),signLS);

                        // Check for slope error
                        if (std::abs(reinit_arr(i,j,k)) <= OUTERTUBE && flag_arr(i,j,k) > NO) {
                            error_arr(i,j,k) = std::abs(nablaG - 1.0);
                        }
                    });
                }

                // Find max error
                error_max = 10.0;
                if (icount % 5 == 0) {
                    error_max = error.max(0);            
                } 
            }

            // Fill ghost cells
            LSreinit.FillBoundary();

            // Clamp levelset values
            for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const amrex::Box& ghost_bx = mfi.growntilebox();

                const auto& reinit_arr = LSreinit.array(mfi);
                const auto& Tube_arr = TubeMask_imf->const_array(mfi);

                amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Skip cells outside domain
                    if (!physical_domain.contains(i,j,k)) return;

                    // Loop through all cells in domain including those outside tube
                    Set::Scalar LS = reinit_arr(i,j,k);
                    if (Tube_arr(i,j,k,CPTcomp) == object_id && LS >= 0.0) LS = -std::abs(LS);
                    reinit_arr(i,j,k) = std::clamp(LS, -OUTERTUBE, OUTERTUBE);
                });
            }

            LSreinit.FillBoundary();

            // Copy to flow
            LSflow.ParallelCopy(LSreinit, 0, ls_id, 1, ngrow, ngrow);
            LSflow.FillBoundary();
        }
    };

    // Object-level metadata and levelset id mapping
    int num_objs = 0;
    std::vector<std::vector<ObjectMeta>> objects;                  // One entry per object
    std::vector<std::unordered_map<int, std::vector<int>>> levelset_object_map; // Used to map object id's to levelset
};

namespace Integrator
{
class NarrowBandLevelset : virtual public Integrator
{
public:
    static constexpr const char* name = "narrowbandlevelset";

    // Constructors
    NarrowBandLevelset(int a_nghost = 3); // Default constructor
    NarrowBandLevelset(IO::ParmParse& pp); // Constructor that triggers parse
    static void Parse(NarrowBandLevelset& value, IO::ParmParse& pp);

    // Destructor
    virtual ~NarrowBandLevelset()
    {
        delete ic_ls;
        delete ic_velocity;
        delete bc_ls;
        delete bc_velocity;
        
    }
    
    // Flux Reconstruction Schemes
    enum class FluxReconstruction {
        FirstOrder,               // Basic first-order reconstruction
        ThirdOrderENO,            // Third-order Essentially Non-Oscillatory (ENO)
        FifthOrderWENO,           // Fifth-order Weighted Essentially Non-Oscillatory (WENO)
    };

    // Flux Schemes
    enum class FluxScheme {
        NarrowbandFluxComputation
    };

    enum class TimeSteppingScheme {
        ForwardEuler,           // First-order explicit Euler method
        RK3                    // Third-order Runge-Kutta method
    };

    // Subclass to manage feature maps
    class FeatureMaps {
    public:
        // Singleton access to ensure one instance
        static FeatureMaps& getInstance() {
            static FeatureMaps instance;
            return instance;
        }

        // Accessor methods for the maps
        const std::map<std::string, FluxReconstruction>& getFluxReconstructionMap() const {
            return fluxReconstructionMap;
        }

        const std::map<std::string, FluxScheme>& getFluxSchemeMap() const {
            return fluxSchemeMap;
        }

        const std::map<std::string, TimeSteppingScheme>& getTimeSteppingSchemeMap() const {
            return timeSteppingSchemeMap;
        }

    private:
        // Private constructor to initialize the maps
        FeatureMaps() {
            fluxReconstructionMap = {
                {"FirstOrder", FluxReconstruction::FirstOrder},
                {"ThirdOrderENO", FluxReconstruction::ThirdOrderENO},
                {"FifthOrderWENO", FluxReconstruction::FifthOrderWENO},
            };

            fluxSchemeMap = {
                {"NarrowbandFluxComputation", FluxScheme::NarrowbandFluxComputation}
            };

            timeSteppingSchemeMap = {
                {"ForwardEuler", TimeSteppingScheme::ForwardEuler},
                {"RK3", TimeSteppingScheme::RK3}
            };
        }

        // Disable copy and assignment
        FeatureMaps(const FeatureMaps&) = delete;
        FeatureMaps& operator=(const FeatureMaps&) = delete;

        // The maps for different schemes
        std::map<std::string, FluxReconstruction> fluxReconstructionMap;
        std::map<std::string, FluxScheme> fluxSchemeMap;
        std::map<std::string, TimeSteppingScheme> timeSteppingSchemeMap;
    };

    // Accessors to the singleton
    static const FeatureMaps& getFeatureMaps() {
        return FeatureMaps::getInstance();
    }

protected:
    // Variables (MultiFabs)
    Set::Field<Set::Scalar>    LS_mf;                   // Levelset field variable (current timestep)
    Set::Field<Set::Scalar>    LSold_mf;                // Levelset field variable (previous timestep)
    Set::Field<Set::Scalar>    NormalX_mf;              // Normal field in x-direction
#if AMREX_SPACEDIM >= 2
    Set::Field<Set::Scalar>    NormalY_mf;              // Normal field in y-direction
#endif
#if AMREX_SPACEDIM == 3
    Set::Field<Set::Scalar>    NormalZ_mf;              // Normal field in x-direction
#endif
    Set::Field<Set::Scalar>    Curvature_mf;            // Curvature field
    Set::Field<Set::Scalar>    LSvel_mf;                // Velocity field in each dimension

    // Will try to convert to iMultiFabs in future
    Set::Field<Set::Scalar>    Tube_mf;                 // Narrowband Tube field variable for plotting
    Set::Field<Set::Scalar>    Zero_mf;                 // Zero levelset field variable for plotting
    Set::Field<Set::Scalar>    CPT_mf;                  // CPT field variable for plotting

    std::string velocity_method;                        // Determines the velocity field based on validation case

    // Solver setup and control functions (override)
    void Initialize(int lev) override;
    void TimeStepBegin(Set::Scalar time, int iter) override;
    void TimeStepComplete(Set::Scalar time, int lev) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;
    void Regrid(int lev, Set::Scalar time) override;

    // Timestepping
    void ComputeAndSetNewTimeStep();
    Set::Scalar GetTimeStep();
    
    // Advance
    void Advect(int lev, Set::Scalar time, Set::Scalar dt);
    void UpdateInterfaceVelocity(int lev);
    
    std::shared_ptr<Numeric::FluxHandler<NarrowBandLevelset>> fluxHandler;
    std::shared_ptr<Numeric::TimeStepper<NarrowBandLevelset>> timeStepper;

    template <typename T>
    friend class Numeric::FluxReconstruction;

    template <typename T>
    friend class Numeric::FluxMethod;    

    template <typename T>
    friend class Numeric::FluxHandler;    

    template <typename T>
    friend class Numeric::FirstOrderReconstruction;

    template <typename T>
    friend class Numeric::FifthOrderWENOReconstruction;

    template <typename T>
    friend class Numeric::NarrowbandFluxComputation;

    template <typename T>
    friend class Numeric::TimeSteppingScheme;

    template <typename T>
    friend class Numeric::EulerForwardScheme;

    template <typename T>
    friend class Numeric::RK3Scheme;

private:
    //
    // Definition of parameters set only at instantiation by
    // constructors. 
    //
    int current_timestep            = 0;
    int number_of_ghost_cells       = -1;     // Number of ghost cells
    int number_of_levelsets         = -1;
    
    // Define LevelSetData structure object
    LevelSetData LSData;

    //
    // Definition of user-determined pointer variables.
    //

    IC::IC<Set::Scalar>* ic_ls = nullptr;
    IC::IC<Set::Scalar>* ic_velocity = nullptr;

    BC::BC<Set::Scalar>* bc_ls = nullptr;
    BC::BC<Set::Scalar>* bc_velocity = nullptr;
    BC::Nothing bc_nothing;
    
    FluxReconstruction reconstruction_method;    // Selected flux reconstruction scheme
    FluxScheme flux_scheme;                      // Selected flux scheme
    TimeSteppingScheme temporal_scheme;          // Selected time-stepping scheme

    Set::Scalar cflNumber = NAN;                 // CFL number for flow (read from input)
    Set::Scalar cflReinit = 0.4;                 // CFL for reinitialization only

}; //NarrowBandLevelset Class

} // namespace Integrator
#endif

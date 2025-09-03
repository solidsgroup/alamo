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
#include <type_traits>  // For std::integral_constant

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
            case FR::FirstOrder:     return 2;
            case FR::ThirdOrderENO:  return 2;
            case FR::FifthOrderWENO: return 3;
            default:                 return 3;  // fallback
        }
    }

    static constexpr int RequiredReinitInterval(FluxReconstruction scheme) {
        using FR = FluxReconstruction;
        switch (scheme) {
            case FR::FirstOrder:     return 1;
            case FR::ThirdOrderENO:  return 5;
            case FR::FifthOrderWENO: return 1;
            default:                 return 5;  // fallback
        }
    }

    int GhostCells() const {
        return RequiredGhosts(reconstruction_method);
    }

    int GetReinitInterval() const {
        return RequiredReinitInterval(reconstruction_method);
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
        int material_id = -1;
        bool has_narrowband = false;
        bool Reinitialize = false;
        int ReinitOrder = 1;
        int ReinitInterval = -1;
        bool IS_MOVING = true;
        
        // Narrowband geometry information
        amrex::Box Tube_domain;
        amrex::BoxArray Tube_ba;
        amrex::DistributionMapping Tube_dm;

        // Field Data
        std::unique_ptr<amrex::MultiFab>  TubeLS_mf;
        std::unique_ptr<amrex::iMultiFab> TubeMask_imf;
        std::unique_ptr<amrex::MultiFab>  TubeFlux_mf;

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

        enum class TubeFluxComp : int {
            FLUX_X = 0,
    #if AMREX_SPACEDIM >= 2
            FLUX_Y,
    #endif
    #if AMREX_SPACEDIM == 3
            FLUX_Z,
    #endif
            NUM_FLUX_COMPONENTS
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
            amrex::GpuArray<Set::Scalar, AMREX_SPACEDIM> dirichlet_vals = {0.0};

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

            const int ngrow = FlowData->GhostCells();

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
                    // Grow the tile by ngrow to safely ensure tiles that contain necessary
                    // ghost cells are included even if tile does not contain narrowband cells
                    const amrex::Box overlap = amrex::grow(tile, ngrow) & valid_bx;
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

            // Define constants
            const int n_real_comps = Comp(TubeFieldComp::NUM_REAL_COMPONENTS);
            const int n_int_comps  = Comp(TubeIntComp::NUM_INT_COMPONENTS);
            const int n_flux_comps = Comp(TubeFluxComp::NUM_FLUX_COMPONENTS);
            const int ngrow = FlowData->GhostCells();

            // Always allocate real field (MultiFab)
            TubeLS_mf = std::make_unique<amrex::MultiFab>(Tube_ba, Tube_dm, n_real_comps, ngrow);

            // Always allocate integer mask (iMultiFab)
            TubeMask_imf = std::make_unique<amrex::iMultiFab>(Tube_ba, Tube_dm, n_int_comps, ngrow);

            // Allocate fluxes (always real-valued)
            TubeFlux_mf = std::make_unique<amrex::MultiFab>(Tube_ba, Tube_dm, n_flux_comps, ngrow);
        }

        void CopyLSToLSOld()
        {
            const int LScomp = Comp(TubeFieldComp::LS);
            const int LSoldcomp = Comp(TubeFieldComp::LS_OLD);

            for (amrex::MFIter mfi(*TubeLS_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& ghost_bx = mfi.growntilebox();
                auto arr = TubeLS_mf->array(mfi);

                amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    arr(i,j,k,LSoldcomp) = arr(i,j,k,LScomp);
                });
            } 
        }

        // Deallocate memory
        void DeallocateFabs() {
            if (TubeLS_mf) TubeLS_mf.reset();
            if (TubeMask_imf) TubeMask_imf.reset();
            if (TubeFlux_mf) TubeFlux_mf.reset();
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

            // Ensure ghost cells are up to date
            TubeLS_mf->FillBoundary(LScomp, 1);
            TubeLS_mf->FillBoundary(LSoldcomp, 1);

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
            TubeMask_imf->FillBoundary();
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

            // Apply ghosts
            TubeLS_mf->FillBoundary(vel_comp[0], AMREX_SPACEDIM);

            // ApplyBCs
            ApplyVelBCs(bc_config.velocity_bc);
        }

    #if AMREX_SPACEDIM == 1
        void CopyGeometryToFlowFields(amrex::MultiFab& NormalXflow) const
        {
            const int Xcomp = Comp(TubeFieldComp::NORMAL_X);

            NormalXflow.ParallelCopy(*TubeLS_mf, Xcomp, ls_id, 1);
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
        }
    #endif

        // Function to ensure new zero crossing is within tube
        void CheckZeroCrossings(amrex::MultiFab& LSflow, amrex::MultiFab& LSOldflow) {
            // Constants
            const int LScomp      = Comp(TubeFieldComp::LS);
            const int LSoldcomp   = Comp(TubeFieldComp::LS_OLD);

            const int OUTSIDETUBE = NarrowBandTubeType::OUTSIDETUBE;
            const int Tubecomp    = Comp(TubeIntComp::TUBE_MASK);
            const int Zerocomp    = Comp(TubeIntComp::ZERO_LS);

            // Loop through TubeLS_mf
            for (amrex::MFIter mfi(*TubeLS_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // Use validbox
                const auto valid_bx = mfi.validbox();

                auto const& LS_arr    = TubeLS_mf->array(mfi);
                auto const& Tube_arr  = TubeMask_imf->const_array(mfi);

                amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // SKip non-narrowband cells
                    if (Tube_arr(i,j,k,Tubecomp) == OUTSIDETUBE) return;

                    // Find new crossings
                    Set::Scalar& LS = LS_arr(i,j,k,LScomp);
                    Set::Scalar& LSold = LS_arr(i,j,k,LSoldcomp);
                    // New crossing but not original interface cell
                    if (LS * LSold < 0.0 && Tube_arr(i,j,k,Zerocomp) == 0) LS = LSold;
                });
            }

            // Update internal ghost cells
            TubeLS_mf->FillBoundary(LScomp, 1);

            // Copy back to flow domain
            LSflow.ParallelCopy(*TubeLS_mf, LScomp, ls_id, 1);
            LSOldflow.ParallelCopy(*TubeLS_mf, LSoldcomp, ls_id, 1);
            LSflow.FillBoundary();
            LSOldflow.FillBoundary();
        }

        // Function to update CPT flags after advection
        void UpdateCPTFlagAfterAdvection() {
            // Constants
            const int LScomp      = Comp(TubeFieldComp::LS);

            const int OUTSIDETUBE = NarrowBandTubeType::OUTSIDETUBE;
            const int Tubecomp    = Comp(TubeIntComp::TUBE_MASK);
            const int Zerocomp    = Comp(TubeIntComp::ZERO_LS);
            const int CPTcomp     = Comp(TubeIntComp::CPT_ID);

            // Loop through TubeLS_mf
            for (amrex::MFIter mfi(*TubeLS_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // Use validbox
                const auto valid_bx = mfi.validbox();

                auto const& LS_arr    = TubeLS_mf->const_array(mfi);
                auto const& Tube_arr  = TubeMask_imf->array(mfi);

                amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Skip non-narrowband cells
                    if (Tube_arr(i,j,k,Tubecomp) == OUTSIDETUBE) return;

                    // Skip non-interface cells
                    if (Tube_arr(i,j,k,Zerocomp) == 0) return;

                    // Set CPT to 0
                    const Set::Scalar LS = LS_arr(i,j,k,LScomp);
                    Tube_arr(i,j,k,CPTcomp) = (LS <= 0.0) ? object_id : 0;
                });
            }

            // Fill internal ghosts
            TubeMask_imf->FillBoundary(CPTcomp, 1);
        }

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
                // Pass 4: Bounding box for EDGEPOINTs
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

                tmpvel.FillBoundary();

                // Reallocate pointers
                AllocateFabs(Tubeflow, ls_id);

                // Copy back to pointers
                LoadLSFields(LSflow, LSoldflow);
                LoadIntFields(Tubeflow, Zeroflow, CPTflow);

                /*for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    TubeLS_mf->ParallelCopy(tmpvel, d, vel_comp[d], 1, ngrow, ngrow);
                }
                TubeLS_mf->FillBoundary(vel_comp[0], AMREX_SPACEDIM);*/
            }
        }

        void ComputeGeometryQuantities() 
        {
            // Define DX for gradient
            const amrex::Real* DX = FlowData->geom.CellSize();

            // Define OutsideNarrowband constant
            const int OUTSIDETUBE = NarrowBandTubeType::OUTSIDETUBE;

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
                , normalZcomp
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

        // --- Reinitialization structures/functions ---
        // Struct to hold the common parameters
        struct StencilParams {
            const amrex::Array4<const Set::Scalar>& LS_arr;
            const amrex::Array4<const int>& Tube_arr;
            const int i, j, k;
            const int LScomp, Tubecomp, OUTSIDETUBE;
            const amrex::Dim3& lo;
            const amrex::Dim3& hi;
        };

        // Template to build stencil based on Order
        template<int Order>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void BuildStencil(
            Set::Scalar (&stencil)[AMREX_SPACEDIM][(Order == 1 ? 3 : 5)],
            const StencilParams& params,
            const amrex::IntVect debug)
        {
            // Compile-time dispatch based on Order template parameter
            if constexpr (Order == 1) {
                BuildFirstOrderStencil(stencil,
                                        params.LS_arr, params.Tube_arr,
                                        params.i, params.j, params.k,
                                        debug,
                                        params.LScomp, params.Tubecomp,
                                        params.OUTSIDETUBE,
                                        params.lo, params.hi);
            } else { // Only Order==1 or Order==2 supported
                BuildSecondOrderStencil(stencil,
                                        params.LS_arr, params.Tube_arr,
                                        params.i, params.j, params.k,
                                        debug,
                                        params.LScomp, params.Tubecomp,
                                        params.OUTSIDETUBE,
                                        params.lo, params.hi);
            }
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void BuildFirstOrderStencil(Set::Scalar (&stencil)[AMREX_SPACEDIM][3],
                            const amrex::Array4<const Set::Scalar>& LS_arr,
                            const amrex::Array4<const int>& Tube_arr,
                            int i, int j, int k,
                            const amrex::IntVect debug,
                            const int LScomp,
                            const int Tubecomp,
                            const int OUTSIDETUBE,
                            const amrex::Dim3& lo,
                            const amrex::Dim3& hi)
        {
            // Define index and lo/hi vectors
            int idx[3] = {i, j, k};
            int lo_v[3] = {lo.x, lo.y, lo.z};
            int hi_v[3] = {hi.x, hi.y, hi.z};

            const int is_debug = amrex::IntVect(AMREX_D_DECL(i,j,k)) == debug;

            // Loop through dimensions
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                // Neighbor offsets in (i,j,k)
                int di = (d==0);
                int dj = (d==1);
                int dk = (d==2);

                // Index shortcuts
                const int id = idx[d];
                const int lo_d = lo_v[d];
                const int hi_d = hi_v[d];

                const int is_lo = (id == lo_d);
                const int is_hi = (id == hi_d);

                // Load stencil values (raw, uncorrected)
                stencil[d][0] = LS_arr(i-1*di, j-1*dj, k-1*dk, LScomp);
                stencil[d][1] = LS_arr(i,      j,      k,      LScomp);
                stencil[d][2] = LS_arr(i+1*di, j+1*dj, k+1*dk, LScomp);

                // Mask checks
                const int out_m = Tube_arr(i-1*di, j-1*dj, k-1*dk, Tubecomp) == OUTSIDETUBE;
                const int out_p = Tube_arr(i+1*di, j+1*dj, k+1*dk, Tubecomp) == OUTSIDETUBE;

                // Apply corrections branchlessly
                auto &s = stencil[d]; // alias
                s[0] = (is_lo || out_m) ? s[1] : s[0];
                s[2] = (is_hi || out_p) ? s[1] : s[2];

                /*if (is_debug) {
                    amrex::Print() << "Debugging point " << debug << " in direction " << d << std::endl;
                    amrex::Print() << "di: " << di << ", dj: "<< dj << ", dk: " << dk << std::endl;
                    amrex::Print() << "LSLeft: "<< LS_arr(i-1*di, j-1*dj, k-1*dk, LScomp) << std::endl;
                    amrex::Print() << "LS: "<< LS_arr(i,      j,      k,      LScomp) << std::endl;
                    amrex::Print() << "LSRight: "<< LS_arr(i+1*di, j+1*dj, k+1*dk, LScomp) << std::endl;
                    amrex::Print() << "is_lo: " << is_lo << ", is_hi: " << is_hi << std::endl;
                    amrex::Print() << "out_m: " << out_m << ", out_p: " << out_p << std::endl;
                    amrex::Print() << "s[0]: " << s[0] << ", s[1]: " << s[1] << ", s[2]: " << s[2] << std::endl;
                }*/
            }
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void BuildSecondOrderStencil(Set::Scalar (&stencil)[AMREX_SPACEDIM][5],
                            const amrex::Array4<const Set::Scalar>& LS_arr,
                            const amrex::Array4<const int>& Tube_arr,
                            int i, int j, int k,
                            const amrex::IntVect debug,
                            const int LScomp,
                            const int Tubecomp,
                            const int OUTSIDETUBE,
                            const amrex::Dim3& lo,
                            const amrex::Dim3& hi)
        {
            // Define index and lo/hi vectors
            int idx[3] = {i, j, k};
            int lo_v[3] = {lo.x, lo.y, lo.z};
            int hi_v[3] = {hi.x, hi.y, hi.z};

            const int is_debug = amrex::IntVect(AMREX_D_DECL(i,j,k)) == debug;

            // Loop through dimensions
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                // Neighbor offsets in (i,j,k)
                int di = (d==0);
                int dj = (d==1);
                int dk = (d==2);

                // Index shortcuts
                const int id = idx[d];
                const int lo_d = lo_v[d];
                const int hi_d = hi_v[d];

                const int is_lo  = (id == lo_d);
                const int is_hi  = (id == hi_d);
                const int is_lo1 = (id == lo_d+1);
                const int is_hi1 = (id == hi_d-1);

                // Load stencil values (raw, uncorrected)
                stencil[d][0] = LS_arr(i-2*di, j-2*dj, k-2*dk, LScomp); // West/South/Bottom order 2
                stencil[d][1] = LS_arr(i-1*di, j-1*dj, k-1*dk, LScomp); // West/South/Bottom
                stencil[d][2] = LS_arr(i,      j,      k,      LScomp); // LP (point)
                stencil[d][3] = LS_arr(i+1*di, j+1*dj, k+1*dk, LScomp); // East/North/Top
                stencil[d][4] = LS_arr(i+2*di, j+2*dj, k+2*dk, LScomp); // East/North/Top order 2

                // Mask checks
                const int out_m2  = (Tube_arr(i-2*di, j-2*dj, k-2*dk, Tubecomp) == OUTSIDETUBE);
                const int out_m1  = (Tube_arr(i-1*di, j-1*dj, k-1*dk, Tubecomp) == OUTSIDETUBE);
                const int out_p1  = (Tube_arr(i+1*di, j+1*dj, k+1*dk, Tubecomp) == OUTSIDETUBE);
                const int out_p2  = (Tube_arr(i+2*di, j+2*dj, k+2*dk, Tubecomp) == OUTSIDETUBE);
                const int out_neg = out_m2 || out_m1;
                const int out_pos = out_p2 || out_p1;

                // Apply corrections branchlessly
                auto &s = stencil[d]; // alias
                s[0] = (is_lo || out_neg) ? s[2] : s[0];
                s[1] = (is_lo || out_neg) ? s[2] : s[1];
                s[3] = (is_hi || out_pos) ? s[2] : s[3];
                s[4] = (is_hi || out_pos) ? s[2] : s[4];

                // Edge cases (mirror neighbors)
                s[0] = is_lo1 ? s[1] : s[0];
                s[4] = is_hi1 ? s[3] : s[4];

                /*if (is_debug) {
                    amrex::Print() << "Debugging point " << debug << " in direction " << d << std::endl;
                    amrex::Print() << "di: " << di << ", dj: "<< dj << ", dk: " << dk << std::endl;
                    amrex::Print() << "LSLeft2: "<< LS_arr(i-2*di, j-2*dj, k-2*dk, LScomp) << std::endl;
                    amrex::Print() << "LSLeft: "<< LS_arr(i-1*di, j-1*dj, k-1*dk, LScomp) << std::endl;
                    amrex::Print() << "LS: "<< LS_arr(i,      j,      k,      LScomp) << std::endl;
                    amrex::Print() << "LSRight: "<< LS_arr(i+1*di, j+1*dj, k+1*dk, LScomp) << std::endl;
                    amrex::Print() << "LSRight2: "<< LS_arr(i+2*di, j+2*dj, k+2*dk, LScomp) << std::endl;
                    amrex::Print() << "is_lo: " << is_lo << ", is_hi: " << is_hi << std::endl;
                    amrex::Print() << "out_neg: " << out_neg << ", out_pos: " << out_pos << std::endl;
                    amrex::Print() << "s[0]: " << s[0] << ", s[1]: " << s[1] << std::endl;
                    amrex::Print() << "s[2]: " << s[2] << ", s[3]: " << s[3] << ", s[4]: " << s[4] << std::endl;
                }*/
            }
        }

        template<int Order>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void DiffSquared(const StencilParams& params,
                        const Set::Scalar* DX,
                        const Set::Scalar signLS,
                        const amrex::IntVect debug,
                        Set::Scalar& dfx, Set::Scalar& dfy, Set::Scalar& dfz)
        {
            constexpr int stencil_size = (Order == 1 ? 3 : 5);
            Set::Scalar stencil[AMREX_SPACEDIM][stencil_size];

            BuildStencil<Order>(stencil, params, debug);

            if constexpr (Order == 1) {
                FirstOrderDiffSquared(stencil, params, DX, signLS, debug, dfx, dfy, dfz);
            } else {
                SecondOrderDiffSquared(stencil, params, DX, signLS, debug, dfx, dfy, dfz);
            }
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void FirstOrderDiffSquared(
            const Set::Scalar (&stencil)[AMREX_SPACEDIM][3],
            const StencilParams params,
            const Set::Scalar* DX,
            const Set::Scalar signLS,
            const amrex::IntVect debug,
            Set::Scalar& dfx,
            Set::Scalar& dfy,
            Set::Scalar& dfz)
        {
            // Initialize fluxes
            dfx = 0.0;
            dfy = 0.0;
            dfz = 0.0;

            const int is_debug = amrex::IntVect(AMREX_D_DECL(params.i, params.j, params.k)) == debug;

            // Loop through dimensions and compute flux
            Set::Scalar flux[AMREX_SPACEDIM];
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                // Get stencil alias
                const Set::Scalar* s = stencil[d];

                // Define DX component
                const Set::Scalar idx = 1.0 / DX[d];

                // Compute derivatives
                const Set::Scalar dxm = (s[1] - s[0]) * idx;
                const Set::Scalar dxp = (s[2] - s[1]) * idx;

                // Compute a+/- & b+/-
                const Set::Scalar a = amrex::max(signLS * dxm, 0.0); 
                const Set::Scalar b = amrex::min(signLS * dxp, 0.0);

                // Compute directional flux (dfx/dfy/dfz)
                flux[d] = amrex::max(a*a, b*b);

                /*if (is_debug) {
                    amrex::Print() << "Debugging point " << debug << " in direction " << d << std::endl;
                    amrex::Print() << "dxm: " << dxm << ", dxp: "<< dxp << std::endl;
                    amrex::Print() << "signLS: " << signLS << std::endl;
                    amrex::Print() << "a: " << a << ", b: "<< b << std::endl;
                    amrex::Print() << "Flux[" << d << "]: " << flux[d] << std::endl;
                }*/
            }

            // Assign corresponding flux components
            dfx = flux[0];
        #if (AMREX_SPACEDIM >= 2)
            dfy = flux[1];
        #endif
        #if (AMREX_SPACEDIM == 3)
            dfz = flux[2];
        #endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void SecondOrderDiffSquared(
            const Set::Scalar (&stencil)[AMREX_SPACEDIM][5],
            const StencilParams params,
            const Set::Scalar* DX,
            const Set::Scalar signLS,
            const amrex::IntVect debug,
            Set::Scalar& dfx,
            Set::Scalar& dfy,
            Set::Scalar& dfz)
        {
            // Initialize fluxes
            dfx = 0.0;
            dfy = 0.0;
            dfz = 0.0;

            const int is_debug = amrex::IntVect(AMREX_D_DECL(params.i, params.j, params.k)) == debug;

            // Loop through dimensions and compute flux
            Set::Scalar flux[AMREX_SPACEDIM];
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                // Get stencil alias
                const Set::Scalar* s = stencil[d];

                // Define DX component
                const Set::Scalar dx = DX[d];
                const Set::Scalar IDX = 1.0 / dx;
                const Set::Scalar IDXSQR = IDX*IDX;

                // Compute first order derivatives
                Set::Scalar dxm  = (s[2] - s[1]) * IDX; // LP - (West/South/Bottom)
                Set::Scalar dxpp = (s[3] - s[2]) * IDX; // (East/North/Top) - LP

                // Compute second order derivatives
                const Set::Scalar dx2Cent  = (s[1] - 2.0*s[2] + s[3]) * IDXSQR; // Central 2nd order
                const Set::Scalar dx2Left  = (s[0] - 2.0*s[1] + s[2]) * IDXSQR;
                const Set::Scalar dx2Right = (s[2] - 2.0*s[3] + s[4]) * IDXSQR;

                // Inline branchless minmod logic
                auto minmod = [](Set::Scalar a, Set::Scalar b) -> Set::Scalar {
                    const Set::Scalar prod = a * b;
                    const Set::Scalar abs_a = amrex::Math::abs(a);
                    const Set::Scalar abs_b = amrex::Math::abs(b);
                    return (prod <= 0.0) ? 0.0 : ((abs_a < abs_b) ? a : b);
                };

                // Adjust dxm/dxpp
                dxm  += 0.5 * dx * minmod(dx2Cent, dx2Left);
                dxpp -= 0.5 * dx * minmod(dx2Cent, dx2Right);

                // Directional Upwinding
                const Set::Scalar dxmS = dxm * signLS, dxpS = dxpp * signLS;

                // Compute directional flux (dfx/dfy/dfz)
                const Set::Scalar flux_d = (dxpS < 0.0 && dxmS < -dxpS) ? dxpp :
                    (dxmS >= 0.0 && dxpS >= -dxmS) ? dxm : 0.0;
                flux[d] = flux_d * flux_d;

                /*if (is_debug) {
                    amrex::Print() << "Debugging point " << debug << " in direction " << d << std::endl;
                    amrex::Print() << "dxm: " << dxm << ", dxpp: "<< dxpp << std::endl;
                    amrex::Print() << "dx2Left: " << dx2Left << ", dx2Cent: "<< dx2Cent <<", dx2Right: "<< dx2Right << std::endl;
                    amrex::Print() << "signLS: " << signLS << std::endl;
                    amrex::Print() << "dxmS: " << dxmS << ", dxps: "<< dxpS << std::endl;
                    amrex::Print() << "flux_d: " << flux_d << std::endl;
                    amrex::Print() << "flux[" << d << "]: " << flux[d] << std::endl;
                }*/
            }

            // Assign corresponding flux components
            dfx = flux[0];
        #if (AMREX_SPACEDIM >= 2)
            dfy = flux[1];
        #endif
        #if (AMREX_SPACEDIM == 3)
            dfz = flux[2];
        #endif
        }

        /*AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void FirstOrderDiffSquared(
            const amrex::Array4<const Set::Scalar>& LS_arr,
            const amrex::Array4<const int>& Tube_arr,
            int i, int j, int k,
            const amrex::Dim3& lo,
            const amrex::Dim3& hi,
            const amrex::Real* DX,
            const int LScomp,
            const int Tubecomp,
            const int OUTSIDETUBE,
            Set::Scalar signLS,
            Set::Scalar& dfx,
            Set::Scalar& dfy,
            Set::Scalar& dfz)
        {
            dfx = 0.0;
            dfy = 0.0;
            dfz = 0.0;

            // x-direction
            if (i > lo.x) {
                Set::Scalar dxm = (Tube_arr(i-1,j,k,Tubecomp) == OUTSIDETUBE) ? 0.0 : 
                    Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(LS_arr, i, j, k, LScomp, DX, Numeric::XLo);
                const Set::Scalar a = amrex::max(signLS * dxm, 0.0);
                dfx = a*a;
            }
            if (i < hi.x) {
                Set::Scalar dxp = (Tube_arr(i+1,j,k,Tubecomp) == OUTSIDETUBE) ? 0.0 : 
                    Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(LS_arr, i, j, k, LScomp, DX, Numeric::XHi);
                const Set::Scalar b = amrex::min(signLS * dxp, 0.0);
                dfx = amrex::max(dfx, b*b);
            }

        #if AMREX_SPACEDIM >= 2
            // y-direction
            if (j > lo.y) {
                Set::Scalar dym = (Tube_arr(i,j-1,k,Tubecomp) == OUTSIDETUBE) ? 0.0 :
                    Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(LS_arr, i, j, k, LScomp, DX, Numeric::YLo);
                const Set::Scalar c = amrex::max(signLS * dym, 0.0);
                dfy = c*c;
            }
            if (j < hi.y) {
                Set::Scalar dyp = (Tube_arr(i,j+1,k,Tubecomp) == OUTSIDETUBE) ? 0.0 :
                    Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(LS_arr, i, j, k, LScomp, DX, Numeric::YHi);
                const Set::Scalar d = amrex::min(signLS * dyp, 0.0);
                dfy = amrex::max(dfy, d*d);
            }
        #endif

        #if AMREX_SPACEDIM == 3
            // z-direction
            if (k > lo.z) {
                Set::Scalar dzm = (Tube_arr(i,j,k-1,Tubecomp) == OUTSIDETUBE) ? 0.0 :
                    Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(LS_arr, i, j, k, LScomp, DX, Numeric::ZLo);
                const Set::Scalar e = amrex::max(signLS * dzm, 0.0);
                dfz = e*e;
            }
            if (k < hi.z) {
                Set::Scalar dzp = (Tube_arr(i,j,k+1,Tubecomp) == OUTSIDETUBE) ? 0.0 :
                    Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(LS_arr, i, j, k, LScomp, DX, Numeric::ZHi);
                const Set::Scalar f = amrex::min(signLS * dzp, 0.0);
                dfz = amrex::max(dfz, f*f);
            }
        #endif
        }*/

        // Switch b/w First/Second Order based on Order value
        void ReinitializeLSFirstOrder(amrex::MultiFab& LSflow, amrex::MultiFab& Signflow)
        {
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
            const int ngrow = 1;
            const amrex::Real* DX = FlowData->geom.CellSize();
            const amrex::Dim3 lo = amrex::lbound(physical_domain);
            const amrex::Dim3 hi = amrex::ubound(physical_domain);
            const amrex::IntVect debug(37,67);

            // Allocate and initialize temporary multifabs
            amrex::MultiFab LSreinit(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, ngrow);
            amrex::MultiFab LSreinitold(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, ngrow);
            amrex::MultiFab SignLS(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, ngrow);
            LSreinit.setVal(10.0 * OUTERTUBE);
            SignLS.setVal(-2.0);

            // Set LS reinit and SignLS values
            for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // Loop through ghosts
                const auto& ghost_bx = mfi.growntilebox(); 

                // Get array views
                const auto& LS_arr = TubeLS_mf->const_array(mfi);
                const auto& Reinit_arr = LSreinit.array(mfi);
                const auto& Sign_arr = SignLS.array(mfi);

                const auto& Tube_arr = TubeMask_imf->const_array(mfi);
                
                // First loop - set LSreinit and freeze INTERFACE cells
                amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Skip cells outside physical domain
                    if (!physical_domain.contains(i,j,k)) return;

                    // Define LS, mask, CPT, and initial val
                    const Set::Scalar LS = LS_arr(i,j,k,LScomp);
                    const int mask = Tube_arr(i,j,k,Tubecomp);
                    const int CPT = Tube_arr(i,j,k,CPTcomp);
                    const int Zero = Tube_arr(i,j,k,Zerocomp);

                    // Initialize LSreinit to 10.0 * OUTERTUBE with correct sign
                    Set::Scalar& Reinitval = Reinit_arr(i,j,k,0);
                    Reinitval = (CPT == object_id) ? -Reinitval : Reinitval;

                    // Loop through narrowband cells
                    if (mask == OUTSIDETUBE) return;

                    // Set the LSreinit cell to LS after advection
                    Reinitval = (CPT == object_id) ? -amrex::Math::abs(LS) : LS; 

                    // Freeze the old and new INTERFACE cells
                    auto& S = Sign_arr(i,j,k,0);
                    if (Zero != 0) {
                        S = 0.0;

                        // Loop through neighbors
                        for (int d = 0; d < Neighbors::num_neighbors; ++d) {
                            int ni = i + Neighbors::offsets[d][0];
                            int nj = j + Neighbors::offsets[d][1];
                            int nk = k + Neighbors::offsets[d][2];

                            // Skip if outside ghost box
                            if (!ghost_bx.contains(ni,nj,nk)) continue;

                            // Skip if nbr is outside domain
                            if (!physical_domain.contains(ni,nj,nk)) continue;

                            // Check for sign change and that both LS values are within tube
                            const Set::Scalar LSnbr = LS_arr(ni,nj,nk,LScomp);

                            if (LS * LSnbr <= 0.0 && std::max(LS, LSnbr) < INNERTUBE) {
                                S = 0.0;
                                Sign_arr(ni,nj,nk) = 0.0;
                            }
                        } 
                    }
                });

                // Second loop - Tube expansion
                amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Skip cells outside physical domain
                    if (!physical_domain.contains(i,j,k)) return;

                    // Skip non-narrowband cells
                    if (Tube_arr(i,j,k,Tubecomp) == OUTSIDETUBE) return;

                    // Skip INTERFACE cells
                    auto& S = Sign_arr(i,j,k,0);
                    if (S == 0.0) return;

                    // Get current LS value
                    const Set::Scalar LS = Reinit_arr(i,j,k,0);

                    // Loop through neighbors
                    for (int d = 0; d < Neighbors::num_neighbors; ++d) {
                        int ni = i + Neighbors::offsets[d][0];
                        int nj = j + Neighbors::offsets[d][1];
                        int nk = k + Neighbors::offsets[d][2];

                        // Skip if outside ghost box
                        if (!ghost_bx.contains(ni,nj,nk)) continue;

                        // Skip if nbr is outside domain
                        if (!physical_domain.contains(ni,nj,nk)) continue;

                        // skip if nbr is outside tube
                        if (Tube_arr(ni,nj,nk,Tubecomp) == OUTSIDETUBE) continue;

                        // Skip if nbr is INTERFACE
                        if (Sign_arr(ni,nj,nk,0) == 0.0) continue;

                        // Set sign of cell to its nbr
                        auto& LSnbr = Reinit_arr(ni,nj,nk,0); 
                        LSnbr = std::copysign(LSnbr, LS);
                    } 

                    // Update SignLS
                    S = (LS > 0.0) ? 1.0 : -1.0;
                });
            }

            // Fill ghost cells
            LSreinit.FillBoundary();
            SignLS.FillBoundary();

            // Plot signLS
            Signflow.ParallelCopy(SignLS, 0, ls_id, 1);
            Signflow.FillBoundary();

            // Copy LSreinit to LSreinitold
            LSreinitold.ParallelCopy(LSreinit, 0, 0, 1, ngrow, ngrow);
            LSreinitold.FillBoundary();

            // Define variables for timestep loop
            amrex::MultiFab error(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, 0);
            error.setVal(0.0);

            const Set::Scalar toler_error_reinit = 1e-3;
            const int LSREINIT_MAX_ITERS = 100;
            const Set::Scalar cflReinit = 0.4;
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

                // Main PDE loop
                for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    // Determine the overlapping region between flow_domain and Tube_domain
                    const auto& valid_bx = mfi.validbox();

                    // Get array views
                    const auto& Tube_arr = TubeMask_imf->const_array(mfi);
                    const auto& Sign_arr = SignLS.const_array(mfi);
                    const auto& reinit_arr = LSreinit.array(mfi);
                    const auto& reinitold_arr = LSreinitold.array(mfi);
                    
                    const auto& error_arr = error.array(mfi);

                    amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        // Skip non-narrowband or INTERFACE cells
                        const int& mask = Tube_arr(i,j,k,Tubecomp);
                        const Set::Scalar& signLS = Sign_arr(i,j,k,0);
                        if (mask == OUTSIDETUBE || signLS == 0.0) return;
                        
                        // Compute gradients
                        Set::Scalar dfx, dfy, dfz;
                        const StencilParams params = {reinitold_arr, Tube_arr, i, j, k, 
                            LScomp, Tubecomp, OUTSIDETUBE, lo, hi};
                        DiffSquared<1>(params, DX, signLS, debug, dfx, dfy, dfz);

                        // Compute gradient and update
                        const Set::Scalar nablaG = std::sqrt(dfx + dfy + dfz);
                        const Set::Scalar LSold = reinitold_arr(i,j,k);
                        Set::Scalar& LS = reinit_arr(i,j,k);
                        LS = LSold - DTreinit * signLS * (nablaG - 1.0);
                        LS = std::copysign(LS,signLS);

                        // Check for slope error
                        if (std::abs(LS) <= OUTERTUBE) {
                            error_arr(i,j,k) = std::abs(nablaG - 1.0);
                        }
                    });
                }

                // Find max error
                error_max = 10.0;
                if (icount % 5 == 0) {
                    error_max = error.max(0); 
                    error.setVal(0.0);  
                    amrex::Print() << "Iter: " << icount << std::endl;
                    amrex::Print() << "error_max: " << error_max << std::endl;        
                } 
            }

            // Clamp OUTSIDETUBE to +/- OUTERTUBE
            for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const amrex::Box& valid_bx = mfi.validbox();

                const auto& reinit_arr = LSreinit.array(mfi);
                const auto& Tube_arr = TubeMask_imf->const_array(mfi);

                amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    Set::Scalar& LS = reinit_arr(i,j,k);
                    LS = std::clamp(LS, -OUTERTUBE, OUTERTUBE);
                });
            }

            // Copy to flow
            LSflow.ParallelCopy(LSreinit, 0, ls_id, 1);
            LSflow.FillBoundary();
        }

        void ReinitializeLSSecondOrder(amrex::MultiFab& LSflow, amrex::MultiFab& Signflow)
        {
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
            const int ngrow = 2;
            const amrex::Real* DX = FlowData->geom.CellSize();
            const amrex::Dim3 lo = amrex::lbound(physical_domain);
            const amrex::Dim3 hi = amrex::ubound(physical_domain);
            const amrex::IntVect debug(-4,-4);

            // Allocate and initialize temporary multifabs
            amrex::MultiFab LSreinit(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, ngrow);
            amrex::MultiFab LSreinitold(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, ngrow);
            amrex::MultiFab SignLS(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, ngrow);
            LSreinit.setVal(10.0 * OUTERTUBE);
            SignLS.setVal(-2.0);

            // Set LS reinit and SignLS values
            for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // Loop through ghosts
                const auto& ghost_bx = mfi.growntilebox(); 

                // Get array views
                const auto& LS_arr = TubeLS_mf->const_array(mfi);
                const auto& Reinit_arr = LSreinit.array(mfi);
                const auto& Sign_arr = SignLS.array(mfi);

                const auto& Tube_arr = TubeMask_imf->const_array(mfi);
                
                // First loop - set LSreinit and freeze INTERFACE cells
                amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Skip cells outside physical domain
                    if (!physical_domain.contains(i,j,k)) return;

                    // Define LS, mask, CPT, and initial val
                    const Set::Scalar LS = LS_arr(i,j,k,LScomp);
                    const int mask = Tube_arr(i,j,k,Tubecomp);
                    const int CPT = Tube_arr(i,j,k,CPTcomp);
                    const int Zero = Tube_arr(i,j,k,Zerocomp);

                    // Initialize LSreinit to 10.0 * OUTERTUBE with correct sign
                    Set::Scalar& Reinitval = Reinit_arr(i,j,k,0);
                    Reinitval = (CPT == object_id) ? -Reinitval : Reinitval;

                    // Loop through narrowband cells
                    if (mask == OUTSIDETUBE) return;

                    // Set the LSreinit cell to LS after advection
                    Reinitval = (CPT == object_id) ? -amrex::Math::abs(LS) : LS; 

                    // Freeze the old and new INTERFACE cells
                    auto& S = Sign_arr(i,j,k,0);
                    if (Zero != 0) {
                        S = 0.0;

                        // Loop through neighbors
                        for (int d = 0; d < Neighbors::num_neighbors; ++d) {
                            int ni = i + Neighbors::offsets[d][0];
                            int nj = j + Neighbors::offsets[d][1];
                            int nk = k + Neighbors::offsets[d][2];

                            // Skip if outside ghost box
                            if (!ghost_bx.contains(ni,nj,nk)) continue;

                            // Skip if nbr is outside domain
                            if (!physical_domain.contains(ni,nj,nk)) continue;

                            // Check for sign change and that both LS values are within tube
                            const Set::Scalar LSnbr = LS_arr(ni,nj,nk,LScomp);

                            if (LS * LSnbr <= 0.0 && std::max(LS, LSnbr) < INNERTUBE) {
                                S = 0.0;
                                Sign_arr(ni,nj,nk) = 0.0;
                            }
                        } 
                    }
                });

                // Second loop - Tube expansion
                amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Skip cells outside physical domain
                    if (!physical_domain.contains(i,j,k)) return;

                    // Skip non-narrowband cells
                    if (Tube_arr(i,j,k,Tubecomp) == OUTSIDETUBE) return;

                    // Skip INTERFACE cells
                    auto& S = Sign_arr(i,j,k,0);
                    if (S == 0.0) return;

                    // Get current LS value
                    const Set::Scalar LS = Reinit_arr(i,j,k,0);

                    // Loop through neighbors
                    for (int d = 0; d < Neighbors::num_neighbors; ++d) {
                        int ni = i + Neighbors::offsets[d][0];
                        int nj = j + Neighbors::offsets[d][1];
                        int nk = k + Neighbors::offsets[d][2];

                        // Skip if outside ghost box
                        if (!ghost_bx.contains(ni,nj,nk)) continue;

                        // Skip if nbr is outside domain
                        if (!physical_domain.contains(ni,nj,nk)) continue;

                        // skip if nbr is outside tube
                        if (Tube_arr(ni,nj,nk,Tubecomp) == OUTSIDETUBE) continue;

                        // Skip if nbr is INTERFACE
                        if (Sign_arr(ni,nj,nk,0) == 0.0) continue;

                        // Set sign of cell to its nbr
                        auto& LSnbr = Reinit_arr(ni,nj,nk,0); 
                        LSnbr = std::copysign(LSnbr, LS);
                    } 

                    // Update SignLS
                    S = (LS > 0.0) ? 1.0 : -1.0;
                });
            }

            // Fill ghost cells
            LSreinit.FillBoundary();
            SignLS.FillBoundary();

            // Plot signLS
            Signflow.ParallelCopy(SignLS, 0, ls_id, 1);
            Signflow.FillBoundary();

            // Copy LSreinit to LSreinitold
            LSreinitold.ParallelCopy(LSreinit, 0, 0, 1, ngrow, ngrow);
            LSreinitold.FillBoundary();

            // Define variables for timestep loop
            amrex::MultiFab error(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, 0);
            error.setVal(0.0);

            const Set::Scalar toler_error_reinit = 1e-2;
            const int LSREINIT_MAX_ITERS = 100;
            const Set::Scalar cflReinit = 0.4;
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

                // Main PDE loop
                for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    // Determine the overlapping region between flow_domain and Tube_domain
                    const auto& valid_bx = mfi.validbox();

                    // Get array views
                    const auto& Tube_arr = TubeMask_imf->const_array(mfi);
                    const auto& Sign_arr = SignLS.const_array(mfi);
                    const auto& reinit_arr = LSreinit.array(mfi);
                    const auto& reinitold_arr = LSreinitold.array(mfi);
                    
                    const auto& error_arr = error.array(mfi);

                    amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        // Skip non-narrowband or INTERFACE cells
                        const int& mask = Tube_arr(i,j,k,Tubecomp);
                        const Set::Scalar& signLS = Sign_arr(i,j,k,0);
                        if (mask == OUTSIDETUBE || signLS == 0.0) return;

                        const int is_debug = amrex::IntVect(AMREX_D_DECL(i,j,k)) == debug;

                        
                        // Compute gradients
                        Set::Scalar dfx2, dfy2, dfz2;
                        const StencilParams params = {reinitold_arr, Tube_arr, i, j, k, 
                            LScomp, Tubecomp, OUTSIDETUBE, lo, hi};
                        DiffSquared<2>(params, DX, signLS, debug, dfx2, dfy2, dfz2);

                        // Compute gradient and update
                        const Set::Scalar nablaG = std::sqrt(dfx2 + dfy2 + dfz2);
                        const Set::Scalar LSold = reinitold_arr(i,j,k);
                        Set::Scalar& LS = reinit_arr(i,j,k);
                        LS = LSold - DTreinit * signLS * (nablaG - 1.0);
                        LS = std::copysign(LS,signLS);
                        if (is_debug) {
                            //amrex::Print() << "nablaG: " << nablaG << std::endl;
                            //amrex::Print() << "LS: " << LS << std::endl;
                        }

                        // Check for slope error
                        if (std::abs(LS) <= OUTERTUBE) {
                            error_arr(i,j,k) = amrex::Math::abs(nablaG - 1.0);
                        }
                    });
                }

                // Find max error
                error_max = 10.0;
                if (icount % 5 == 0) {
                    error_max = error.max(0); 

                    /*// Loop through error
                    Set::Scalar max_val = -std::numeric_limits<Set::Scalar>::infinity();
                    amrex::IntVect max_loc(AMREX_D_DECL(-1, -1, -1));
                    for (amrex::MFIter mfi(error, false); mfi.isValid(); ++mfi)
                    {
                        const amrex::FArrayBox& fab = error[mfi];
                        const amrex::Box& bx = mfi.validbox();
                        const auto& arr = fab.array();

                        for (amrex::BoxIterator bit(bx); bit.ok(); ++bit)
                        {
                            const amrex::IntVect& iv = bit();
                            Set::Scalar val = arr(iv, 0);
                            if (val > max_val)
                            {
                                max_val = val;
                                max_loc = iv;
                            }
                        }
                    }*/
                    error.setVal(0.0);   
                    amrex::Print() << "Iter: " << icount << std::endl;
                    amrex::Print() << "error_max: " << error_max << std::endl;  
                    //amrex::Print() << "Max error at " << max_loc << std::endl;      
                } 
            }

            // Clamp OUTSIDETUBE to +/- OUTERTUBE
            for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const amrex::Box& valid_bx = mfi.validbox();

                const auto& reinit_arr = LSreinit.array(mfi);
                const auto& Tube_arr = TubeMask_imf->const_array(mfi);

                amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    Set::Scalar& LS = reinit_arr(i,j,k);
                    LS = std::clamp(LS, -OUTERTUBE, OUTERTUBE);
                });
            }

            // Copy to flow
            LSflow.ParallelCopy(LSreinit, 0, ls_id, 1);
            LSflow.FillBoundary();
        }

        /*AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Set::Scalar minmod(Set::Scalar a, Set::Scalar b)
        {
            if (a*b <= 0.0) return 0.0;
            return (amrex::Math::abs(a) < amrex::Math::abs(b)) ? a : b;
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void SecondOrderUpwindDerivatives(
            const amrex::Array4<const Set::Scalar>& LS_arr,
            const amrex::Array4<const Set::Scalar>& Tube_arr,
            int i, int j, int k,
            const amrex::IntVect debug,
            const amrex::Dim3& lo,
            const amrex::Dim3& hi,
            const Set::Scalar* DX,
            const int LScomp,
            const int Tubecomp,
            const int OUTSIDETUBE,
            Set::Scalar signLS,
            Set::Scalar& dfx,
            Set::Scalar& dfy,
            Set::Scalar& dfz)
        {
            using namespace Numeric;

            // --- X ---
            // Define boolean coonditioons
            const bool iminus2_out_tube = Tube_arr(i-2,j,k,Tubecomp) == OUTSIDETUBE;
            const bool iminus_out_tube = Tube_arr(i-1,j,k,Tubecomp) == OUTSIDETUBE;
            const bool iplus2_out_tube = Tube_arr(i+2,j,k,Tubecomp) == OUTSIDETUBE;
            const bool iplus_out_tube = Tube_arr(i+1,j,k,Tubecomp) == OUTSIDETUBE;
            const bool minus_nbrs_out_tube = iminus2_out_tube || iminus_out_tube;
            const bool plus_nbrs_out_tube = iplus2_out_tube || iplus_out_tube;

            const bool is_lo_x  = (i == lo.x);
            const bool is_lo1_x = (i == lo.x + 1);
            const bool is_hi_x  = (i == hi.x);
            const bool is_hi1_x = (i == hi.x - 1);

            // Get point to check if debug point
            const amrex::IntVect point = amrex::IntVect(AMREX_D_DECL(i,j,k));

            // Define stencil points (i-2, i-1, i, i+1, i+2) for flux computations
            Set::Scalar iminus2 = LS_arr(i-2,j,k,LScomp);
            Set::Scalar iminus = LS_arr(i-1,j,k,LScomp);
            const Set::Scalar p = LS_arr(i-1,j,k,LScomp);
            Set::Scalar iplus = LS_arr(i+1,j,k,LScomp);
            Set::Scalar iplus2 = LS_arr(i+2,j,k,LScomp);

            // Correct stencil points at physical/tube boundaries
            if (is_lo_x || minus_nbrs_out_tube) {
                iminus2 = p;
                iminus = p;
            }
            if (is_hi_x || plus_nbrs_out_tube) {
                iplus2 = p;
                iplus = p;
            }
            if (is_lo1_x) iminus2 = iminus;
            if (is_hi1_x) iplus2 = iplus;

            // Define denominator terms
            const Set::Scalar idx = 1.0 / DX[0];
            const Set::Scalar idx2 = 1.0 / (DX[0] * DX[0]);

            // Compute fluxes
            Set::Scalar dxm = (p - iminus) * idx;
            Set::Scalar dxpp = (iplus - p) * idx2;
            if (point == debug) {
                amrex::Print()<< "LS: " << p << std::endl;
                amrex::Print() << "LS_i-2: " << iminus2 << ", LS_i-1: " << iminus << std::endl;
                amrex::Print() << "LS_i+1: " << iplus << ", LS_i+2: " << iplus2 << std::endl;
                amrex::Print() << "dxm: " << dxm <<", dxpp: " << dxpp << std::endl;
            }

            const Set::Scalar dx2Cent = (iminus - 2.0*p + iplus) * idx2;
            const Set::Scalar dx2Left = (iminus2 - 2.0*iminus + p) * idx2;
            const Set::Scalar dx2Right = (p - 2.0*iplus + iplus2) * idx2;

            dxm  += 0.5 * DX[0] * minmod(dx2Cent, dx2Left);
            dxpp -= 0.5 * DX[0] * minmod(dx2Cent, dx2Right);

            if (point == debug) {
                amrex::Print()<< "dx2Cent: " << dx2Cent << std::endl;
                amrex::Print()<< "dx2Left: " << dx2Left << ", dx2Right: " << dx2Right << std::endl;
                amrex::Print() << "dxm: " << dxm <<", dxpp: " << dxpp << std::endl;
            }

            // --- Y ---
            dfy = 0.0;
        #if AMREX_SPACEDIM >= 2
            const bool is_lo_y  = (j == lo.y);
            const bool is_hi_y  = (j == hi.y);
            const bool is_lo1_y = (j == lo.y + 1);
            const bool is_hi1_y = (j == hi.y - 1);

            Set::Scalar dym = is_lo_y ? 0.0 : Stencil<Set::Scalar,0,1,0>::D(arr,i,j,k,0,DX,YLo);
            Set::Scalar dypp = is_hi_y ? 0.0 : Stencil<Set::Scalar,0,1,0>::D(arr,i,j,k,0,DX,YHi);
            if (point == debug) {
                amrex::Print() << "LS_j-2: " << arr(i,j-2,k) << ", LS_j-1: " << arr(i,j-1,k) << std::endl;
                amrex::Print() << "LS_j+1: " << arr(i,j+1,k) << ", LS_j+2: " << arr(i,j+2,k) << std::endl;
                amrex::Print() << "dym: " << dym <<", dypp: " << dypp << std::endl;
            }

            Set::Scalar dy2Cent = (!is_lo_y && !is_hi_y)
                ? Stencil<Set::Scalar,0,2,0>::D(arr,i,j,k,0,DX,DefaultType)
                : (is_lo_y ? (arr(i,j+1,k) - arr(i,j,k)) : (arr(i,j-1,k) - arr(i,j,k))) / (DX[1] * DX[1]);

            Set::Scalar dy2Left = (!is_lo_y && !is_lo1_y)
                ? Stencil<Set::Scalar,0,2,0>::D(arr,i,j-1,k,0,DX,DefaultType)
                : (is_lo1_y ? (arr(i,j,k) - arr(i,j-1,k)) / (DX[1] * DX[1]) : 0.0);

            Set::Scalar dy2Right = (!is_hi_y && !is_hi1_y)
                ? Stencil<Set::Scalar,0,2,0>::D(arr,i,j+1,k,0,DX,DefaultType)
                : (is_hi1_y ? (arr(i,j,k) - arr(i,j+1,k)) / (DX[1] * DX[1]) : 0.0);

            dym  += 0.5 * DX[1] * minmod(dy2Cent, dy2Left);
            dypp -= 0.5 * DX[1] * minmod(dy2Cent, dy2Right);
            if (point == debug) {
                amrex::Print()<< "dy2Cent: " << dy2Cent << std::endl;
                amrex::Print()<< "dy2Left: " << dy2Left << ", dy2Right: " << dy2Right << std::endl;
                amrex::Print() << "dym: " << dym <<", dypp: " << dypp << std::endl;
            }
        #endif

            // --- Z ---
            dfz = 0.0;
        #if AMREX_SPACEDIM == 3
            const bool is_lo_z  = (k == lo.z);
            const bool is_hi_z  = (k == hi.z);
            const bool is_lo1_z = (k == lo.z + 1);
            const bool is_hi1_z = (k == hi.z - 1);

            Set::Scalar dzm  = is_lo_z ? 0.0 : Stencil<Set::Scalar,0,0,1>::D(arr,i,j,k,0,DX,ZLo);
            Set::Scalar dzpp = is_hi_z ? 0.0 : Stencil<Set::Scalar,0,0,1>::D(arr,i,j,k,0,DX,ZHi);

            Set::Scalar dz2Cent = (!is_lo_z && !is_hi_z)
                ? Stencil<Set::Scalar,0,0,2>::D(arr,i,j,k,0,DX,DefaultType)
                : (is_lo_z ? (arr(i,j,k+1) - arr(i,j,k)) : (arr(i,j,k-1) - arr(i,j,k))) / (DX[2] * DX[2]);

            Set::Scalar dz2Left = (!is_lo_z && !is_lo1_z)
                ? Stencil<Set::Scalar,0,0,2>::D(arr,i,j,k-1,0,DX,DefaultType)
                : (is_lo1_z ? (arr(i,j,k) - arr(i,j,k-1)) / (DX[2] * DX[2]) : 0.0);

            Set::Scalar dz2Right = (!is_hi_z && !is_hi1_z)
                ? Stencil<Set::Scalar,0,0,2>::D(arr,i,j,k+1,0,DX,DefaultType)
                : (is_hi1_z ? (arr(i,j,k) - arr(i,j,k+1)) / (DX[2] * DX[2]) : 0.0);

            dzm  += 0.5 * DX[2] * minmod(dz2Cent, dz2Left);
            dzpp -= 0.5 * DX[2] * minmod(dz2Cent, dz2Right);
        #endif

            // --- Directional Upwinding ---
            const Set::Scalar dxmS = dxm * signLS, dxpS = dxpp * signLS;

            #if AMREX_SPACEDIM >= 2
            const Set::Scalar dymS = dym * signLS, dypS = dypp * signLS;
            #else
            const Set::Scalar dymS = 0.0, dypS = 0.0;
            #endif

            #if AMREX_SPACEDIM == 3
            const Set::Scalar dzmS = dzm * signLS, dzpS = dzpp * signLS;
            #else
            const Set::Scalar dzmS = 0.0, dzpS = 0.0;
            #endif


            dfx = (dxpS < 0.0 && dxmS < -dxpS) ? dxpp :
                (dxmS >= 0.0 && dxpS >= -dxmS) ? dxm : 0.0;

        #if AMREX_SPACEDIM >= 2
            dfy = (dypS < 0.0 && dymS < -dypS) ? dypp :
                (dymS >= 0.0 && dypS >= -dymS) ? dym : 0.0;
        #endif

        #if AMREX_SPACEDIM == 3
            dfz = (dzpS < 0.0 && dzmS < -dzpS) ? dzpp :
                (dzmS >= 0.0 && dzpS >= -dzmS) ? dzm : 0.0;
        #endif
        }*/

        /*AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Set::Scalar minmod(Set::Scalar a, Set::Scalar b)
        {
            if (a*b <= 0.0) return 0.0;
            return (std::abs(a) < std::abs(b)) ? a : b;
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void SecondOrderUpwindDerivatives(
            const amrex::Array4<const Set::Scalar>& arr,
            int i, int j, int k,
            const amrex::IntVect debug,
            const amrex::Dim3& lo,
            const amrex::Dim3& hi,
            const Set::Scalar* DX,
            Set::Scalar signLS,
            Set::Scalar& dfx,
            Set::Scalar& dfy,
            Set::Scalar& dfz)
        {
            using namespace Numeric;

            // --- X ---
            const bool is_lo_x  = (i == lo.x);
            const bool is_hi_x  = (i == hi.x);
            const bool is_lo1_x = (i == lo.x + 1);
            const bool is_hi1_x = (i == hi.x - 1);
            const amrex::IntVect point = amrex::IntVect(AMREX_D_DECL(i,j,k));

            Set::Scalar dxm = is_lo_x ? 0.0 : Stencil<Set::Scalar,1,0,0>::D(arr,i,j,k,0,DX,XLo);
            Set::Scalar dxpp = is_hi_x ? 0.0 : Stencil<Set::Scalar,1,0,0>::D(arr,i,j,k,0,DX,XHi);
            if (point == debug) {
                amrex::Print()<< "LS: " << arr(i,j,k) << std::endl;
                amrex::Print() << "LS_i-2: " << arr(i-2,j,k) << ", LS_i-1: " << arr(i-1,j,k) << std::endl;
                amrex::Print() << "LS_i+1: " << arr(i+1,j,k) << ", LS_i+2: " << arr(i+2,j,k) << std::endl;
                amrex::Print() << "dxm: " << dxm <<", dxpp: " << dxpp << std::endl;
            }

            Set::Scalar dx2Cent = (!is_lo_x && !is_hi_x)
                ? Stencil<Set::Scalar,2,0,0>::D(arr,i,j,k,0,DX,DefaultType)
                : (is_lo_x ? (arr(i+1,j,k) - arr(i,j,k)) : (arr(i-1,j,k) - arr(i,j,k))) / (DX[0] * DX[0]);

            Set::Scalar dx2Left = (!is_lo_x && !is_lo1_x)
                ? Stencil<Set::Scalar,2,0,0>::D(arr,i-1,j,k,0,DX,DefaultType)
                : (is_lo1_x ? (arr(i,j,k) - arr(i-1,j,k)) / (DX[0] * DX[0]) : 0.0);

            Set::Scalar dx2Right = (!is_hi_x && !is_hi1_x)
                ? Stencil<Set::Scalar,2,0,0>::D(arr,i+1,j,k,0,DX,DefaultType)
                : (is_hi1_x ? (arr(i,j,k) - arr(i+1,j,k)) / (DX[0] * DX[0]) : 0.0);

            dxm  += 0.5 * DX[0] * minmod(dx2Cent, dx2Left);
            dxpp -= 0.5 * DX[0] * minmod(dx2Cent, dx2Right);

            if (point == debug) {
                amrex::Print()<< "dx2Cent: " << dx2Cent << std::endl;
                amrex::Print()<< "dx2Left: " << dx2Left << ", dx2Right: " << dx2Right << std::endl;
                amrex::Print() << "dxm: " << dxm <<", dxpp: " << dxpp << std::endl;
            }

            // --- Y ---
            dfy = 0.0;
        #if AMREX_SPACEDIM >= 2
            const bool is_lo_y  = (j == lo.y);
            const bool is_hi_y  = (j == hi.y);
            const bool is_lo1_y = (j == lo.y + 1);
            const bool is_hi1_y = (j == hi.y - 1);

            Set::Scalar dym = is_lo_y ? 0.0 : Stencil<Set::Scalar,0,1,0>::D(arr,i,j,k,0,DX,YLo);
            Set::Scalar dypp = is_hi_y ? 0.0 : Stencil<Set::Scalar,0,1,0>::D(arr,i,j,k,0,DX,YHi);
            if (point == debug) {
                amrex::Print() << "LS_j-2: " << arr(i,j-2,k) << ", LS_j-1: " << arr(i,j-1,k) << std::endl;
                amrex::Print() << "LS_j+1: " << arr(i,j+1,k) << ", LS_j+2: " << arr(i,j+2,k) << std::endl;
                amrex::Print() << "dym: " << dym <<", dypp: " << dypp << std::endl;
            }

            Set::Scalar dy2Cent = (!is_lo_y && !is_hi_y)
                ? Stencil<Set::Scalar,0,2,0>::D(arr,i,j,k,0,DX,DefaultType)
                : (is_lo_y ? (arr(i,j+1,k) - arr(i,j,k)) : (arr(i,j-1,k) - arr(i,j,k))) / (DX[1] * DX[1]);

            Set::Scalar dy2Left = (!is_lo_y && !is_lo1_y)
                ? Stencil<Set::Scalar,0,2,0>::D(arr,i,j-1,k,0,DX,DefaultType)
                : (is_lo1_y ? (arr(i,j,k) - arr(i,j-1,k)) / (DX[1] * DX[1]) : 0.0);

            Set::Scalar dy2Right = (!is_hi_y && !is_hi1_y)
                ? Stencil<Set::Scalar,0,2,0>::D(arr,i,j+1,k,0,DX,DefaultType)
                : (is_hi1_y ? (arr(i,j,k) - arr(i,j+1,k)) / (DX[1] * DX[1]) : 0.0);

            dym  += 0.5 * DX[1] * minmod(dy2Cent, dy2Left);
            dypp -= 0.5 * DX[1] * minmod(dy2Cent, dy2Right);
            if (point == debug) {
                amrex::Print()<< "dy2Cent: " << dy2Cent << std::endl;
                amrex::Print()<< "dy2Left: " << dy2Left << ", dy2Right: " << dy2Right << std::endl;
                amrex::Print() << "dym: " << dym <<", dypp: " << dypp << std::endl;
            }
        #endif

            // --- Z ---
            dfz = 0.0;
        #if AMREX_SPACEDIM == 3
            const bool is_lo_z  = (k == lo.z);
            const bool is_hi_z  = (k == hi.z);
            const bool is_lo1_z = (k == lo.z + 1);
            const bool is_hi1_z = (k == hi.z - 1);

            Set::Scalar dzm  = is_lo_z ? 0.0 : Stencil<Set::Scalar,0,0,1>::D(arr,i,j,k,0,DX,ZLo);
            Set::Scalar dzpp = is_hi_z ? 0.0 : Stencil<Set::Scalar,0,0,1>::D(arr,i,j,k,0,DX,ZHi);

            Set::Scalar dz2Cent = (!is_lo_z && !is_hi_z)
                ? Stencil<Set::Scalar,0,0,2>::D(arr,i,j,k,0,DX,DefaultType)
                : (is_lo_z ? (arr(i,j,k+1) - arr(i,j,k)) : (arr(i,j,k-1) - arr(i,j,k))) / (DX[2] * DX[2]);

            Set::Scalar dz2Left = (!is_lo_z && !is_lo1_z)
                ? Stencil<Set::Scalar,0,0,2>::D(arr,i,j,k-1,0,DX,DefaultType)
                : (is_lo1_z ? (arr(i,j,k) - arr(i,j,k-1)) / (DX[2] * DX[2]) : 0.0);

            Set::Scalar dz2Right = (!is_hi_z && !is_hi1_z)
                ? Stencil<Set::Scalar,0,0,2>::D(arr,i,j,k+1,0,DX,DefaultType)
                : (is_hi1_z ? (arr(i,j,k) - arr(i,j,k+1)) / (DX[2] * DX[2]) : 0.0);

            dzm  += 0.5 * DX[2] * minmod(dz2Cent, dz2Left);
            dzpp -= 0.5 * DX[2] * minmod(dz2Cent, dz2Right);
        #endif

            // --- Directional Upwinding ---
            const Set::Scalar dxmS = dxm * signLS, dxpS = dxpp * signLS;

            #if AMREX_SPACEDIM >= 2
            const Set::Scalar dymS = dym * signLS, dypS = dypp * signLS;
            #else
            const Set::Scalar dymS = 0.0, dypS = 0.0;
            #endif

            #if AMREX_SPACEDIM == 3
            const Set::Scalar dzmS = dzm * signLS, dzpS = dzpp * signLS;
            #else
            const Set::Scalar dzmS = 0.0, dzpS = 0.0;
            #endif


            dfx = (dxpS < 0.0 && dxmS < -dxpS) ? dxpp :
                (dxmS >= 0.0 && dxpS >= -dxmS) ? dxm : 0.0;

        #if AMREX_SPACEDIM >= 2
            dfy = (dypS < 0.0 && dymS < -dypS) ? dypp :
                (dymS >= 0.0 && dypS >= -dymS) ? dym : 0.0;
        #endif

        #if AMREX_SPACEDIM == 3
            dfz = (dzpS < 0.0 && dzmS < -dzpS) ? dzpp :
                (dzmS >= 0.0 && dzpS >= -dzmS) ? dzm : 0.0;
        #endif
        }*/

        /*void ReinitializeLSSecondOrder(amrex::MultiFab& LSflow, amrex::MultiFab& Signflow)
        {
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
            const int ngrow = 2;
            const amrex::Real* DX = FlowData->geom.CellSize();
            const amrex::Dim3 lo = amrex::lbound(physical_domain);
            const amrex::Dim3 hi = amrex::ubound(physical_domain);

            amrex::IntVect debug(AMREX_D_DECL(0,0,0));

            // Allocate and initialize temporary multifabs
            amrex::MultiFab LSreinit(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, ngrow);
            amrex::MultiFab LSreinitold(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, ngrow);
            amrex::MultiFab SignLS(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, ngrow);

            // Set LS reinit and SignLS values
            for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // Loop through ghosts
                const auto& ghost_bx = mfi.growntilebox(); 

                // Get array views
                const auto& LS_arr = TubeLS_mf->const_array(mfi);
                const auto& Reinit_arr = LSreinit.array(mfi);
                const auto& Sign_arr = SignLS.array(mfi);

                const auto& Tube_arr = TubeMask_imf->const_array(mfi);
                
                // First loop - set LSreinit and initialize SignLS
                amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Skip cells outside physical domain
                    if (!physical_domain.contains(i,j,k)) return;

                    // Define LS, mask, CPT, and initial val
                    const Set::Scalar LS = LS_arr(i,j,k,LScomp);
                    const int mask = Tube_arr(i,j,k,Tubecomp);
                    const int CPT = Tube_arr(i,j,k,CPTcomp);
                    const Set::Scalar val = 10.0 * OUTERTUBE;

                    // Initialize LSreinit to 10.0 * OUTERTUBE with correct sign
                    Set::Scalar& Reinitval = Reinit_arr(i,j,k,0);
                    Reinitval = (CPT == object_id) ? -val : val;

                    // Initialize Sign
                    Set::Scalar& S = Sign_arr(i,j,k,0);
                    S = (LS > 0) ? 1.0 : -1.0;

                    // Exit loop if non-narrowband cell
                    if (mask == OUTSIDETUBE) return;

                    // Set LSreinit to LS
                    Reinitval = LS;
                });

                // Second loop - freeze INTERFACE cells
                amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Skip non-INTERFACE cells
                    if (Tube_arr(i,j,k,Zerocomp) == 0) return;

                    // Get current LS value
                    const Set::Scalar LS = LS_arr(i,j,k,LScomp);

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
                        const Set::Scalar LSnbr = LS_arr(ni,nj,nk,LScomp);

                        if (LS * LSnbr <= 0.0 && std::max(LS, LSnbr) < INNERTUBE) {
                            Sign_arr(i,j,k) = 0.0;
                            Sign_arr(ni,nj,nk) = 0.0;
                        }
                    }
                });
            }

            // Fill ghost cells
            LSreinit.FillBoundary();
            SignLS.FillBoundary();

            // Plot signLS
            Signflow.ParallelCopy(SignLS, 0, ls_id, 1);
            Signflow.FillBoundary();

            // Copy LSreinit to LSreinitold
            LSreinitold.ParallelCopy(LSreinit, 0, 0, 1, ngrow, ngrow);
            LSreinitold.FillBoundary();

            // Define variables for timestep loop
            amrex::MultiFab error(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, 0);
            error.setVal(0.0);

            const Set::Scalar toler_error_reinit = 1e-3;
            const int LSREINIT_MAX_ITERS = 100;
            const Set::Scalar cflReinit = 0.4;
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

                // Main PDE loop
                for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    // Determine the overlapping region between flow_domain and Tube_domain
                    const auto& ghost_bx = mfi.growntilebox();
                    const auto& valid_bx = mfi.validbox();

                    // Get array views
                    const auto& Tube_arr = TubeMask_imf->const_array(mfi);
                    const auto& Sign_arr = SignLS.const_array(mfi);
                    const auto& reinit_arr = LSreinit.array(mfi);
                    const auto& reinitold_arr = LSreinitold.array(mfi);
                    
                    const auto& error_arr = error.array(mfi);

                    amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        // Skip non-valid cells
                        if (!valid_bx.contains(i,j,k)) return;
                        
                        // Skip non-narrowband or INTERFACE cells
                        const int& mask = Tube_arr(i,j,k,Tubecomp);
                        const Set::Scalar& signLS = Sign_arr(i,j,k,0);
                        if (mask == OUTSIDETUBE || signLS == 0.0) return;

                        // Compute second order differences
                        Set::Scalar dfx, dfy, dfz;
                        //SecondOrderUpwindDerivatives(reinitold_arr, i, j, k, debug, lo, hi, DX, signLS, dfx, dfy, dfz);

                        // Compute gradient and update
                        const Set::Scalar nablaG = std::sqrt(dfx*dfx + dfy*dfy + dfz*dfz);
                        const Set::Scalar LSold = reinitold_arr(i,j,k);
                        Set::Scalar& LS = reinit_arr(i,j,k);
                        LS = LSold - DTreinit * signLS * (nablaG - 1.0);
                        LS = std::copysign(LS,signLS);

                        // Check for slope error
                        if (std::abs(LS) <= OUTERTUBE) {
                            error_arr(i,j,k) = std::abs(nablaG - 1.0);
                        }
                    });
                }

                // Find max error
                error_max = 10.0;
                if (icount % 5 == 0) {
                    error_max = error.max(0);  
                    error.setVal(0.0);
                    amrex::Print() << "Iter: " << icount << std::endl;
                    amrex::Print() << "error_max: " << error_max << std::endl;         
                } 
            }

            // Loop through error
            Set::Scalar max_val = -std::numeric_limits<Set::Scalar>::infinity();
            amrex::IntVect max_loc(AMREX_D_DECL(-1, -1, -1));
            for (amrex::MFIter mfi(error, false); mfi.isValid(); ++mfi)
            {
                const amrex::FArrayBox& fab = error[mfi];
                const amrex::Box& bx = mfi.validbox();
                const auto& arr = fab.array();

                for (amrex::BoxIterator bit(bx); bit.ok(); ++bit)
                {
                    const amrex::IntVect& iv = bit();
                    Set::Scalar val = arr(iv, 0);
                    if (val > max_val)
                    {
                        max_val = val;
                        max_loc = iv;
                    }
                }
            }
            amrex::Print() << "Max error at " << max_loc << std::endl;

            // Clamp OUTSIDETUBE to +/- OUTERTUBE
            for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const amrex::Box& valid_bx = mfi.validbox();

                const auto& reinit_arr = LSreinit.array(mfi);
                const auto& Tube_arr = TubeMask_imf->const_array(mfi);

                amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    Set::Scalar& LS = reinit_arr(i,j,k);
                    LS = std::clamp(LS, -OUTERTUBE, OUTERTUBE);
                });
            }

            // Copy to flow
            LSflow.ParallelCopy(LSreinit, 0, ls_id, 1);
            LSflow.FillBoundary();
        }*/

        /*void ReinitializeLSSecondOrder(amrex::MultiFab& LSflow, amrex::MultiFab& Signflow)
        {
            // Define (i)MultiFab components
            const int LScomp = Comp(TubeFieldComp::LS);
            const int Signcomp = Comp(TubeFieldComp::SIGN_PHI);

            const int Tubecomp = Comp(TubeIntComp::TUBE_MASK);
            const int CPTcomp  = Comp(TubeIntComp::CPT_ID);

            // Define Tube constants
            const int OUTSIDETUBE = NarrowBandTubeType::OUTSIDETUBE;
            const Set::Scalar INNERTUBE = FlowData->Inner_Tube_width;
            const Set::Scalar OUTERTUBE = FlowData->Tube_width;

            // Define flow constants
            const auto& physical_domain = FlowData->flow_domain;
            const int ngrow = 2;
            const amrex::Real* DX = FlowData->geom.CellSize();
            const amrex::Dim3 lo = amrex::lbound(physical_domain);
            const amrex::Dim3 hi = amrex::ubound(physical_domain);

            // Allocate and initialize temporary multifabs
            amrex::MultiFab LSreinit(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, ngrow);
            amrex::MultiFab LSreinitold(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, ngrow);

            amrex::IntVect debug(AMREX_D_DECL(37,50,0));
            //amrex::Print() << "Debugging point: " << debug << std::endl;
            
            // Properly set LSreinit to 10.0 * OUTERTUBE where the
            // sign depends on the sign of the levelset
            for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // Loop through ghosts
                const auto& ghost_bx = mfi.growntilebox(ngrow);

                // Get array views
                const auto& LS_arr = TubeLS_mf->const_array(mfi);
                const auto& reinit_arr = LSreinit.array(mfi);

                amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const Set::Scalar val = 10.0 * OUTERTUBE;
                    const Set::Scalar LS = LS_arr(i,j,k,LScomp);

                    reinit_arr(i,j,k,0) = (LS > 0) ? val : -val;
                });
            }

            // Fill ghost cells
            LSreinit.FillBoundary();

            // Set LSreinit to narrowband values and properly save signLS
            for (amrex::MFIter mfi(*TubeLS_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // Determine the overlapping region between flow_domain and Tube_domain
                const auto& valid_bx = mfi.validbox();

                // Get array views
                const auto& LS_arr = TubeLS_mf->array(mfi);
                const auto& Tube_arr = TubeMask_imf->const_array(mfi);
                const auto& reinit_arr = LSreinit.array(mfi);

                amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Skip cells outside flow domain
                    if(!physical_domain.contains(i,j,k)) return;

                    // Skip non-narrowband cells
                    if (Tube_arr(i,j,k,Tubecomp) == OUTSIDETUBE) return;

                    // Set the narrowband value in the LSreinit MultiFab
                    const Set::Scalar LS = LS_arr(i,j,k,LScomp);
                    reinit_arr(i,j,k,0) = (Tube_arr(i,j,k,CPTcomp) == object_id) ? 
                        -std::abs(LS) : LS;

                    // Update the sign
                    Set::Scalar& S = LS_arr(i,j,k,Signcomp);
                    S = (LS < 0) ? -1.0 : 1.0;

                    // Set the sign of INTERFACE cells to 0
                    if (std::abs(LS) < INNERTUBE) {
                        // Check nbrs
                        for (int d = 0; d < Neighbors::num_neighbors; ++d) {
                            int ni = i + Neighbors::offsets[d][0];
                            int nj = j + Neighbors::offsets[d][1];
                            int nk = k + Neighbors::offsets[d][2];

                            // Skip nbr cells outside domain boundary
                            if (!physical_domain.contains(ni,nj,nk)) continue; 

                            // Skip if neighbor is not in the LS tube
                            if (Tube_arr(ni,nj,nk,Tubecomp) == OUTSIDETUBE) continue;

                            // Get LS value of neighbor in current direction
                            const Set::Scalar LSnbr = LS_arr(ni,nj,nk,LScomp);

                            // Check if neighbor is on the opposite side of interface
                            if (LS * LSnbr > 0.0 || S == 0.0) continue;
                            
                            // Set sign of INTERFACE cells to 0
                            S = 0.0;
                        }
                    }
                });
            }

            // Fill ghost cells
            LSreinit.FillBoundary();
            TubeLS_mf->FillBoundary(Signcomp, 1);

            // Copy LSreinit to LSreinitold
            LSreinitold.ParallelCopy(LSreinit, 0, 0, 1, ngrow, ngrow);
            LSreinitold.FillBoundary();

            // Plot signLS
            Signflow.ParallelCopy(*TubeLS_mf, Signcomp, ls_id, 1);
            Signflow.FillBoundary();

            // Define variables for timestep loop
            amrex::MultiFab error(TubeLS_mf->boxArray(), TubeLS_mf->DistributionMap(), 1, 0);
            error.setVal(0.0);

            const Set::Scalar toler_error_reinit = 1e-3;
            const int LSREINIT_MAX_ITERS = 100;
            const Set::Scalar cflReinit = 0.4;
            const Set::Scalar DTreinit = cflReinit * FlowData->min_DX();

            Set::Scalar error_max = 2.0 * toler_error_reinit;
            int icount = 0;

            // Begin artifical timestep loop
            while (error_max > toler_error_reinit && icount <= LSREINIT_MAX_ITERS) {
                // Increase iteration count
                icount += 1;

                // Main PDE loop
                for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    // Loop through valid cells
                    const auto& valid_bx = mfi.validbox();

                    // Get array views
                    const auto& LS_arr = TubeLS_mf->const_array(mfi);
                    const auto& Tube_arr = TubeMask_imf->const_array(mfi);
                    const auto& reinit_arr = LSreinit.array(mfi);
                    const auto& reinitold_arr = LSreinitold.const_array(mfi);
                    const auto& error_arr = error.array(mfi);

                    amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        // Skip non-narrowband cells
                        if (Tube_arr(i,j,k,Tubecomp) == OUTSIDETUBE) return;

                        // Skip INTERFACE cells
                        const Set::Scalar& signLS = LS_arr(i,j,k,Signcomp);
                        if (signLS == 0.0) return;

                        // Apply neumann-like BCs at boundaries
                        if (i == lo.x) {
                            reinit_arr(i,j,k,0) = reinitold_arr(i+1,j,k,0);
                            return;
                        } 
                        if (i == hi.x) {
                            reinit_arr(i,j,k,0) = reinitold_arr(i-1,j,k,0);
                            return;
                        } 
    #if AMREX_SPACEDIM >= 2
                        if (j == lo.y) {
                            reinit_arr(i,j,k,0) = reinitold_arr(i,j+1,k,0);
                            return;
                        } 
                        if (j == hi.y) {
                            reinit_arr(i,j,k,0) = reinitold_arr(i,j-1,k,0);
                            return;
                        } 
    #endif
    #if AMREX_SPACEDIM == 3
                        if (k == lo.z) {
                            reinit_arr(i,j,k,0) = reinitold_arr(i,j,k+1,0);
                            return;
                        } 
                        if (k == hi.z) {
                            reinit_arr(i,j,k,0) = reinitold_arr(i,j,k-1,0);
                            return;
                        } 
    #endif

                        // Compute per-direction gradients
                        Set::Scalar dxm = Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(reinitold_arr, i, j, k, 0, DX, Numeric::XLo);
                        Set::Scalar dxpp = Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(reinitold_arr, i, j, k, 0, DX, Numeric::XHi);
                        Set::Scalar dx2Cent = Numeric::Stencil<Set::Scalar, 2, 0, 0>::D(reinitold_arr, i, j, k, 0, DX, Numeric::DefaultType);

                        Set::Scalar dx2Left = Numeric::Stencil<Set::Scalar, 2, 0, 0>::D(reinitold_arr, i-1, j, k, 0, DX, Numeric::DefaultType);
                        Set::Scalar dfx2 = minmod(dx2Cent, dx2Left);
                        dxm += 0.5 * DX[0] * dfx2;
                        /*if (i==debug[0] && j==debug[1]) {
                            amrex::Print() << "phi2L: " << reinitold_arr(i-2,j,k) << std::endl;
                            amrex::Print() << "phiL: " << reinitold_arr(i-1,j,k) << std::endl;
                            if (reinitold_arr(i-1,j,k) > 0)  amrex::Print() << "PhiL > 0" << std::endl;
                            amrex::Print() << "phi: " << reinitold_arr(i,j,k) << std::endl;
                            amrex::Print() << "phiR: " << reinitold_arr(i+1,j,k) << std::endl;
                            amrex::Print() << "phi2R: " << reinitold_arr(i+2,j,k) << std::endl;
                            amrex::Print() << "dx2a: " << dx2a << std::endl;
                            amrex::Print() << "dx2b: " << dx2b << std::endl;
                            amrex::Print() << "dfx2: " << dfx2 << std::endl;
                            amrex::Print() << "dxm: " << dxm << std::endl;
                        }

                        Set::Scalar dx2Right = Numeric::Stencil<Set::Scalar, 2, 0, 0>::D(reinitold_arr, i+1, j, k, 0, DX, Numeric::DefaultType);
                        dfx2 = minmod(dx2Cent, dx2Right);
                        dxpp -= 0.5 * DX[0] * dfx2;
                        /*if (i==debug[0] && j==debug[1]) {
                            amrex::Print() << "dx2b: " << dx2b << std::endl;
                            amrex::Print() << "dfx2: " << dfx2 << std::endl;
                            amrex::Print() << "dxpp: " << dxpp << std::endl;
                        }

                        Set::Scalar dym = 0.0;
                        Set::Scalar dypp = 0.0;
    #if AMREX_SPACEDIM >= 2
                        dym = Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(reinitold_arr, i, j, k, 0, DX, Numeric::YLo);
                        dypp = Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(reinitold_arr, i, j, k, 0, DX, Numeric::YHi);
                        Set::Scalar dy2Cent = Numeric::Stencil<Set::Scalar, 0, 2, 0>::D(reinitold_arr, i, j, k, 0, DX, Numeric::DefaultType);

                        Set::Scalar dy2Left = Numeric::Stencil<Set::Scalar, 0, 2, 0>::D(reinitold_arr, i, j-1, k, 0, DX, Numeric::DefaultType);
                        Set::Scalar dfy2 = minmod(dy2Cent, dy2Left);
                        dym += 0.5 * DX[1] * dfy2;
                        /*if (i==debug[0] && j==debug[1]) {
                            amrex::Print() << "phi2S: " << reinitold_arr(i,j-2,k) << std::endl;
                            amrex::Print() << "phiS: " << reinitold_arr(i,j-1,k) << std::endl;
                            amrex::Print() << "phiN: " << reinitold_arr(i,j+1,k) << std::endl;
                            amrex::Print() << "phi2N: " << reinitold_arr(i,j+2,k) << std::endl;
                            amrex::Print() << "dy2a: " << dy2a << std::endl;
                            amrex::Print() << "dy2b: " << dy2b << std::endl;
                            amrex::Print() << "dfy2: " << dfy2 << std::endl;
                            amrex::Print() << "dym: " << dym << std::endl;
                        }

                        Set::Scalar dy2Right = Numeric::Stencil<Set::Scalar, 0, 2, 0>::D(reinitold_arr, i, j+1, k, 0, DX, Numeric::DefaultType);
                        dfy2 = minmod(dy2Cent, dy2Right);
                        dypp -= 0.5 * DX[1] * dfy2;
                        /*if (i==debug[0] && j==debug[1]) {
                            amrex::Print() << "dy2b: " << dy2b << std::endl;
                            amrex::Print() << "dfy2: " << dfy2 << std::endl;
                            amrex::Print() << "dypp: " << dypp << std::endl;
                        }
    #endif
                        Set::Scalar dzm = 0.0;
                        Set::Scalar dzpp = 0.0;
    #if AMREX_SPACEDIM == 3
                        dzm = Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(reinitold_arr, i, j, k, 0, DX, Numeric::ZLo);
                        dzpp = Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(reinitold_arr, i, j, k, 0, DX, Numeric::ZHi);
                        Set::Scalar dz2a = Numeric::Stencil<Set::Scalar, 0, 0, 2>::D(reinitold_arr, i, j, k, 0, DX, Numeric::DefaultType);

                        Set::Scalar dz2b = Numeric::Stencil<Set::Scalar, 0, 0, 2>::D(reinitold_arr, i, j, k-1, 0, DX, Numeric::DefaultType);
                        Set::Scalar dfz2 = minmod(dz2a, dz2b);
                        dzm += 0.5 * DX[2] * dfz2;

                        dz2b = Numeric::Stencil<Set::Scalar, 0, 0, 2>::D(reinitold_arr, i, j, k+1, 0, DX, Numeric::DefaultType);
                        dfz2 = minmod(dz2a, dz2b);
                        dzpp -= 0.5 * DX[2] * dfz2;
    #endif
                        // Compute fluxes
                        Set::Scalar dxmS = dxm * signLS; Set::Scalar dxpS = dxpp * signLS;
                        Set::Scalar dymS = dym * signLS; Set::Scalar dypS = dypp * signLS;
                        Set::Scalar dzmS = dzm * signLS; Set::Scalar dzpS = dzpp * signLS;
                        /*if (i==debug[0] && j==debug[1]) {
                            amrex::Print() << "dxmS: " << dxmS << std::endl;
                            amrex::Print() << "dxpS: " << dxpS << std::endl;
                            amrex::Print() << "dymS: " << dymS << std::endl;
                            amrex::Print() << "dypS: " << dypS << std::endl;
                            //amrex::Print() << "dzmS: " << dzmS << std::endl;
                            //amrex::Print() << "dzpS: " << dzpS << std::endl;
                        }

                        // Calculate gradient nablaG
                        Set::Scalar dfx = 0.0; Set::Scalar dfy = 0.0; Set::Scalar dfz = 0.0;

                        if (dxpS < 0.0 && dxmS < -dxpS) 
                        {
                            dfx = dxpp;
                        } else if (dxmS >= 0.0 && dxpS >= -dxmS)
                        {
                            dfx = dxm;
                        }

                        if (dypS < 0.0 && dymS < -dypS) 
                        {
                            dfy = dypp;
                        } else if (dymS >= 0.0 && dypS >= -dymS)
                        {
                            dfy = dym;
                        }

                        if (dzpS < 0.0 && dzmS < -dzpS) 
                        {
                            dfz = dzpp;
                        } else if (dzmS >= 0.0 && dzpS >= -dzmS)
                        {
                            dfz = dzm;
                        }
                        /*if (i==debug[0] && j==debug[1]) {
                            amrex::Print() << "dfx: " << dfx << std::endl;
                            amrex::Print() << "dfy: " << dfy << std::endl;
                            //amrex::Print() << "dfz: " << dfz << std::endl;
                        }
                        
                        Set::Scalar nablaG = std::sqrt(dfx*dfx + dfy*dfy + dfz*dfz);
                        if (i==debug[0] && j==debug[1]) {
                            //amrex::Print() << "nablaG: " << nablaG << std::endl;
                        }

                        // Compute distance function
                        reinit_arr(i,j,k,0) = reinitold_arr(i,j,k,0) - DTreinit * signLS * (nablaG - 1.0);
                        reinit_arr(i,j,k,0) = std::copysign(reinit_arr(i,j,k,0), signLS);
                        if (i==debug[0] && j==debug[1]) {
                            //amrex::Print() << "LSreinit(" << debug << "): " << reinit_arr(i,j,k) << std::endl;
                        }

                        // Check for slope error
                        if (std::abs(reinit_arr(i,j,k,0)) <= OUTERTUBE && signLS != 0.0) {
                            error_arr(i,j,k,0) = std::abs(nablaG - 1.0);
                            if (i==debug[0] && j==debug[1]) {
                                //amrex::Print() << "error: " << error_arr(i,j,k,0) << std::endl;
                            }
                        }
                    });
                }

                // Find max error
                error_max = 10.0;
                if (icount % 5 == 0) {
                    //amrex::Print() << "Iter: " << icount << std::endl;
                    error_max = error.max(0);  
                    error.setVal(0.0);
                    //amrex::Print() << "error_max: " << error_max << std::endl;        
                } 

                // Copy LSreinit into LSreinitold for next iteration
                LSreinitold.ParallelCopy(LSreinit, 0, 0, 1, ngrow, ngrow);
                LSreinitold.FillBoundary();
            }

            // Fill ghost cells
            LSreinit.FillBoundary();

            // Loop through error
            Set::Scalar max_val = -std::numeric_limits<Set::Scalar>::infinity();
            amrex::IntVect max_loc(AMREX_D_DECL(-1, -1, -1));
            for (amrex::MFIter mfi(error, false); mfi.isValid(); ++mfi)
            {
                const amrex::FArrayBox& fab = error[mfi];
                const amrex::Box& bx = mfi.validbox();
                const auto& arr = fab.array();

                for (amrex::BoxIterator bit(bx); bit.ok(); ++bit)
                {
                    const amrex::IntVect& iv = bit();
                    Set::Scalar val = arr(iv, 0);
                    if (val > max_val)
                    {
                        max_val = val;
                        max_loc = iv;
                    }
                }
            }
            amrex::Print() << "Max error at " << max_loc << std::endl;
            amrex::Print() << "Max error: " << max_val << std::endl;

            // Clamp levelset values
            for (amrex::MFIter mfi(LSreinit, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const amrex::Box& ghost_bx = mfi.growntilebox();

                const auto& reinit_arr = LSreinit.array(mfi);
                const auto& Tube_arr = TubeMask_imf->const_array(mfi);

                amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Skip cells outside domain
                    if (!physical_domain.contains(i,j,k)) return;

                    Set::Scalar& LS = reinit_arr(i,j,k,0);
                    if (Tube_arr(i,j,k,CPTcomp) == object_id && LS >= 0.0) LS = -std::abs(LS);
                    LS = std::clamp(LS, -OUTERTUBE, OUTERTUBE);
                });
            }

            LSreinit.FillBoundary(); 

            // Copy to flow
            LSflow.ParallelCopy(LSreinit, 0, ls_id, 1);
            LSflow.FillBoundary();
        }*/
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
    Set::Field<Set::Scalar>    SignLS_mf;                  // CPT field variable for plotting

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
    int ReinitOrder                 = -1;
    
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

}; //NarrowBandLevelset Class

} // namespace Integrator
#endif

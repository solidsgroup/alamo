#include <array>

#include "NarrowBandLevelset.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "BC/Constant.H"
#include "BC/Nothing.H"
#include "BC/Expression.H"
#include "Numeric/NBStencil.H"

// Inline helper functions
AMREX_GPU_DEVICE 
AMREX_FORCE_INLINE
bool is_interface (int i, int j, int k,
                   amrex::Array4<const amrex::Real> const& phi)
{
    const amrex::Real phi0 = phi(i,j,k);

    if (phi0 == amrex::Real(0.0)) return true;

    const bool pos = (phi0 > 0.0);

    if ((phi(i-1,j,k) > 0.0) != pos) return true;
    if ((phi(i+1,j,k) > 0.0) != pos) return true;

#if (AMREX_SPACEDIM >= 2)
    if ((phi(i,j-1,k) > 0.0) != pos) return true;
    if ((phi(i,j+1,k) > 0.0) != pos) return true;
#endif

#if (AMREX_SPACEDIM >= 3)
    if ((phi(i,j,k-1) > 0.0) != pos) return true;
    if ((phi(i,j,k+1) > 0.0) != pos) return true;
#endif

    return false;
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
bool has_tube_neighbor (int i, int j, int k,
                        amrex::Array4<const amrex::Real> const& phi,
                        amrex::Real tube_cutoff)
{
    using amrex::Math::abs;

    if (abs(phi(i-1,j,k)) < tube_cutoff) return true;
    if (abs(phi(i+1,j,k)) < tube_cutoff) return true;

#if (AMREX_SPACEDIM >= 2)
    if (abs(phi(i,j-1,k)) < tube_cutoff) return true;
    if (abs(phi(i,j+1,k)) < tube_cutoff) return true;
#endif

#if (AMREX_SPACEDIM == 3)
    if (abs(phi(i,j,k-1)) < tube_cutoff) return true;
    if (abs(phi(i,j,k+1)) < tube_cutoff) return true;
#endif

    return false;
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Scalar to_scalar(int v)
{
    return static_cast<Set::Scalar>(v);
}


namespace NarrowBandData
{
    struct TubeIdx
    {
        static constexpr int NB_MASK   = 0;  // 0/1: in narrow band
        static constexpr int Zero_ID   = 1;  // zero level set id
        static constexpr int OBJ_ID    = 2;  // object id
        static constexpr int TUBE_TYPE = 3;  // enum below

        static constexpr int NCOMP     = 4;
    };

    struct TubeType
    {
        static constexpr int INTERFACE   = 0;
        static constexpr int INNERBAND   = 1;
        static constexpr int INSIDETUBE  = 2;
        static constexpr int EDGEPOINT   = 3;
        static constexpr int OUTSIDETUBE = 4;
    };
};

// =====================================
// Per-object per-level metadata
// =====================================
struct ObjectMetaLevel
{
    amrex::Box bbox;  // bounding box for this object at this AMR level
    amrex::Box prev_bbox; // bounding box from previous time step (for tracking movement in tube update)
};

// =====================================
// Object metadata (global)
// =====================================
struct ObjectMeta
{
    int id = -1;

    amrex::Vector<ObjectMetaLevel> leveldata; // size = number of AMR levels

    // --- Setup helpers ---
    void resize_leveldata(int finest_level)
    {
        leveldata.resize(finest_level + 1);
    }

    bool active = true;
    bool moving = true;
    bool needs_reinit = false;

    amrex::Real volume = 0.0;
    amrex::Real centroid[AMREX_SPACEDIM] = {0.0};
};

// =====================================
// Level set object
// =====================================
struct LevelSetObject
{
    ObjectMeta meta;
};

// -------------------------------
// Level set field (AMR-ready)
// -------------------------------
struct LevelSetField
{
    // Define Field data
    Set::Field<Set::Scalar> LS;
    Set::Field<Set::Scalar> LSold;

    Set::Field<Set::Vector> LSvel;
    Set::Field<Set::Vector> LSnormal;

    Set::Field<Set::Scalar> LSkappa;

    Set::Field<Set::Scalar> TUBE;

    // Store collection of objects in this level set field
    amrex::Vector<LevelSetObject> objects; // Number of objects per field (global)

    // Store IC/BC for level set field
    std::unique_ptr<IC::IC<Set::Scalar>> ic_LS = nullptr;
    std::unique_ptr<BC::BC<Set::Scalar>> bc_LS = nullptr;

    // Store bandwidth constants
    const Set::Scalar HALFINNERTUBE = 4.0;
    const Set::Scalar TUBECUTOFF    = 6.0;

    // Save boolean to force reinitialization of entire level set field
    bool force_reinit_all = false;

    // Helper function for connected component analysis (Placeholder)
    int GetNumObjects() {
        return 1;
    }

    // Helper function to update narrowband tube
    void update_Tube(int lev, const amrex::Geometry& geom, bool init=false)
    {
        // BL_PROFILE("LevelSetField::update_Tube");

        // ----------------------------------------------------------
        // Cutoffs
        // ----------------------------------------------------------
        const Set::Scalar* DX = geom.CellSize();
        Set::Scalar min_dx = *std::min_element(DX, DX + AMREX_SPACEDIM);

        const Set::Scalar inner_cutoff = HALFINNERTUBE * min_dx;
        const Set::Scalar outer_cutoff = TUBECUTOFF * min_dx;

        const amrex::Box& domain = geom.Domain();
        const int Ng = LS[lev]->nGrow();

        // ----------------------------------------------------------
        // Loop objects (you planned for multi-object)
        // ----------------------------------------------------------
        for (LevelSetObject& obj : objects)
        {
            ObjectMetaLevel &meta = obj.meta.leveldata[lev];

            // Save previous bbox
            meta.prev_bbox = meta.bbox;
            meta.bbox = amrex::Box();

            // ------------------------------------------------------
            // Build search region
            // ------------------------------------------------------
            amrex::Box search_box = meta.prev_bbox;

            if (!search_box.ok())
                search_box = domain;

            search_box.grow(Ng);
            search_box &= domain;

            // ------------------------------------------------------
            // Initialize bbox reduction
            // ------------------------------------------------------
            amrex::IntVect lo(AMREX_D_DECL(INT_MAX, INT_MAX, INT_MAX));
            amrex::IntVect hi(AMREX_D_DECL(INT_MIN, INT_MIN, INT_MIN));

            // ------------------------------------------------------
            // Tile loop
            // ------------------------------------------------------
            for (amrex::MFIter mfi(*LS[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& tilebox = mfi.tilebox();

                // Skip irrelevant tiles
                if (!tilebox.intersects(search_box)) continue;

                const amrex::Box workbox = tilebox & search_box;
                if (!workbox.ok()) continue;

                Set::Patch<Set::Scalar> phi_arr  = LS.Patch(lev, mfi);
                Set::Patch<Set::Scalar> tube_arr = TUBE.Patch(lev, mfi);

                // Local bbox
                amrex::IntVect local_lo = lo;
                amrex::IntVect local_hi = hi;

                // --------------------------------------------------
                // GPU kernel
                // --------------------------------------------------
                amrex::ParallelFor(workbox,
                [=, &local_lo, &local_hi] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Set::Scalar& phi = phi_arr(i,j,k);

                    if (init)
                        phi = amrex::Clamp(phi, -outer_cutoff, outer_cutoff);

                    const Set::Scalar abs_phi = amrex::Math::abs(phi);

                    // --------------------------------------------------
                    // Classification (same logic as your original)
                    // --------------------------------------------------
                    bool inside = (abs_phi < outer_cutoff);
                    bool inner  = (abs_phi < inner_cutoff);
                    bool inside_tube = inside && (!inner);

                    bool edge = false;
                    if (!inside)
                        edge = has_tube_neighbor(i,j,k,phi_arr,outer_cutoff);

                    bool interface = inner && is_interface(i,j,k,phi_arr);

                    int tube_type = NarrowBandData::TubeType::OUTSIDETUBE;

                    if (edge)         tube_type = NarrowBandData::TubeType::EDGEPOINT;
                    if (inside_tube)  tube_type = NarrowBandData::TubeType::INSIDETUBE;
                    if (inner)        tube_type = NarrowBandData::TubeType::INNERBAND;
                    if (interface)    tube_type = NarrowBandData::TubeType::INTERFACE;

                    // --------------------------------------------------
                    // Write data
                    // --------------------------------------------------
                    tube_arr(i,j,k,NarrowBandData::TubeIdx::TUBE_TYPE) = to_scalar(tube_type);

                    int nb_mask = (tube_type != NarrowBandData::TubeType::OUTSIDETUBE);
                    tube_arr(i,j,k,NarrowBandData::TubeIdx::NB_MASK) = to_scalar(nb_mask);

                    int obj_id = (phi < 0.0) ? obj.meta.id : -1;
                    tube_arr(i,j,k,NarrowBandData::TubeIdx::OBJ_ID) = to_scalar(obj_id);

                    int zero_id = (tube_type == NarrowBandData::TubeType::INTERFACE) ? obj.meta.id : -1;
                    tube_arr(i,j,k,NarrowBandData::TubeIdx::Zero_ID) = to_scalar(zero_id);

                    // --------------------------------------------------
                    // Update bbox (only if in narrowband)
                    // --------------------------------------------------
                    if (nb_mask)
                    {
                        amrex::Gpu::Atomic::Min(&local_lo[0], i);

    #if AMREX_SPACEDIM >= 2
                        amrex::Gpu::Atomic::Min(&local_lo[1], j);
    #endif

    #if AMREX_SPACEDIM == 3
                        amrex::Gpu::Atomic::Min(&local_lo[2], k);
    #endif

                        amrex::Gpu::Atomic::Max(&local_hi[0], i);

    #if AMREX_SPACEDIM >= 2
                        amrex::Gpu::Atomic::Max(&local_hi[1], j);
    #endif

    #if AMREX_SPACEDIM == 3
                        amrex::Gpu::Atomic::Max(&local_hi[2], k);
    #endif
                    }
                });

                lo = amrex::min(lo, local_lo);
                hi = amrex::max(hi, local_hi);
            }

            // ------------------------------------------------------
            // MPI reduction
            // ------------------------------------------------------
            amrex::ParallelDescriptor::ReduceIntMin(lo.getVect(), AMREX_SPACEDIM);
            amrex::ParallelDescriptor::ReduceIntMax(hi.getVect(), AMREX_SPACEDIM);

            // ------------------------------------------------------
            // Final bbox update
            // ------------------------------------------------------
            if (hi[0] >= lo[0])
            {
                meta.bbox = amrex::Box(lo, hi);
            }
            else
            {
                meta.bbox = amrex::Box(); // empty
            }
        }
    }

    void compute_normal(int lev, const amrex::Geometry& geom) {
        // Set normal to 0.0
        LSnormal[lev]->setVal(Set::Vector::Zero());

        // Get dx per dimension outside loop
        auto const dx = geom.CellSizeArray();

        for (LevelSetObject& obj : objects)
        {
            ObjectMetaLevel &meta = obj.meta.leveldata[lev];

            // Build search region
            const amrex::Box& bbox = meta.bbox;

            for (amrex::MFIter mfi(*LSnormal[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& tilebox = mfi.tilebox();

                // Skip irrelevant tiles
                if (!tilebox.intersects(bbox)) continue;

                const amrex::Box workbox = tilebox & bbox;
                if (!workbox.ok()) continue;

                Set::Patch<const Set::Scalar> phi_arr  = LS.Patch(lev, mfi);
                Set::Patch<const Set::Scalar> tube_arr = TUBE.Patch(lev, mfi);
                Set::Patch<Set::Vector> norm_arr = LSnormal.Patch(lev, mfi);

                amrex::ParallelFor(workbox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    // Skip non-narrowband cells
                    const int bandval = static_cast<int>(tube_arr(i, j, k, NarrowBandData::TubeIdx::NB_MASK));
                    if (bandval == 0) return;

                    // Get the normal gradient
                    Set::Vector grad = NBStencil::norm_gradient(i, j, k, phi_arr, tube_arr, dx);

                    // Normalize
                    Set::Scalar mag = grad.norm();
                    norm_arr(i, j, k) = grad / (mag + 1e-12);
                });
            }
        }
    }
};

// ======================================
// Simulation Domain (AMR FIXED)
// ======================================
struct SimulationDomain
{
    // Store reference to integrator
    Integrator::Integrator* integrator = nullptr;

    // Save ghost cell width (3 for 5th order WENO)
    int nghost = 3;

    // Store collection of level set fields
    amrex::Vector<std::unique_ptr<LevelSetField>> fields;

    // Store velocity field (will copy from flow domain in future)
    Set::Field<Set::Vector> Flowvel;

    std::unique_ptr<IC::IC<Set::Vector>> ic_Vel;

    // Define function to attach integrator to retrieve geometry information
    void Define(Integrator::Integrator* _integrator)
    {
        integrator = _integrator;
    }

    // Define helper functions to retrieve geometry information from intergrator
    const amrex::Geometry& Geom(int lev) const
    {
        return integrator->Geom(lev);
    }

    const amrex::BoxArray& BA(int lev) const
    {
        return integrator->boxArray(lev);
    }

    const amrex::DistributionMapping& DM(int lev) const
    {
        return integrator->DistributionMap(lev);
    }

    int FinestLevel() const
    {
        return integrator->finestLevel();
    }
};

// Define Integrator namespace
namespace Integrator
{
    // Define Constructor/Destructor
    NarrowBandLevelset::NarrowBandLevelset() = default;
    NarrowBandLevelset::NarrowBandLevelset(IO::ParmParse& pp)
        : Integrator()
    {
        Parse(*this, pp);
    }

    NarrowBandLevelset::~NarrowBandLevelset() = default;

    // Define Parser
    void NarrowBandLevelset::Parse(
        NarrowBandLevelset& value,
        IO::ParmParse& pp)
    {
        // -------------------------
        // Allocate domain
        // -------------------------
        value.domain_ = std::make_unique<SimulationDomain>();

        // Attach integrator
        value.domain_->Define(&value);

        // -------------------------
        // Initialize velocity field (Will pull from flow domain in future)
        // -------------------------
        IC::IC<Set::Vector>* ic_vel_tmp = nullptr;

        pp.select_default<IC::Constant, IC::Expression>(
            "velocity.ic",
            ic_vel_tmp,
            value.geom
        );

        value.domain_->ic_Vel =
            std::unique_ptr<IC::IC<Set::Vector>>(ic_vel_tmp);

        value.AddField<Set::Vector, Set::Hypercube::Cell>(
            value.domain_->Flowvel,          // Set::Field<Set::Vector>
            nullptr,                         // BC::BC<Set::Vector>*
            /* ncomp   */ 1,                 // ONE vector per cell
            /* nghost  */ 0,                 // No ghost communication
            /* name    */ "Flowvel",         // plotfile base name
            /* writeout*/ false,             // include in plotfiles
            /* evolving*/ false              // evolves in time

        );
        
        // -------------------------
        // Number of level sets
        // -------------------------
        int nls = 1;
        pp_query_default("ls.number_of_levelsets", nls, 1);

        // -------------------------
        // Allocate fields
        // -------------------------
        value.domain_->fields.resize(nls);

        for (int i = 0; i < nls; ++i)
        {
            value.domain_->fields[i] = std::make_unique<LevelSetField>();
            LevelSetField& ls = *value.domain_->fields[i];

            std::string prefix = "ls" + std::to_string(i);
            IO::ParmParse pp_ls(prefix.c_str());

            ls.ic_LS = nullptr;
            ls.bc_LS = nullptr;

            IC::IC<Set::Scalar>* ic_tmp = nullptr;
            BC::BC<Set::Scalar>* bc_tmp = nullptr;

            // -------------------------
            // IC (AMR FIXED)
            // -------------------------
            pp_ls.select_default<IC::Constant, IC::Expression>(
                "ic",
                ic_tmp,
                value.geom
            );
            ls.ic_LS = std::unique_ptr<IC::IC<Set::Scalar>>(ic_tmp);

            // -------------------------
            // BC
            // -------------------------
            pp_ls.select_default<BC::Constant, BC::Expression>(
                "bc",
                bc_tmp,
                1
            );
            ls.bc_LS = std::unique_ptr<BC::BC<Set::Scalar>>(bc_tmp);

            // -------------------------
            // Register fields (AMR handled internally)
            // -------------------------
            value.RegisterOneField(ls, prefix);
        }
    }

    void NarrowBandLevelset::RegisterOneField(
        LevelSetField& ls,
        const std::string& prefix)
    {
        // ====================================================
        // Evolving level set fields (PDE state)
        // ====================================================

        RegisterNewFab(
            ls.LS,
            ls.bc_LS.get(),
            1,
            domain_->nghost,
            prefix + "_LS",
            true
        );

        RegisterNewFab(
            ls.LSold,
            ls.bc_LS.get(),
            1,
            domain_->nghost,
            prefix + "_LSold",
            true
        );

        // ====================================================
        // Derived / diagnostic fields (non-evolving, no BCs)
        // ====================================================

        // Level-set velocity (vector)
        AddField<Set::Vector, Set::Hypercube::Cell>(
            ls.LSvel,
            nullptr,          // no BCs
            1,                // one vector per cell
            domain_->nghost,                // 
            prefix + "_LSvel",
            true,             // plot
            false             // NOT evolving
        );

        // Level-set normal (vector)
        AddField<Set::Vector, Set::Hypercube::Cell>(
            ls.LSnormal,
            nullptr,
            1,
            domain_->nghost,
            prefix + "_LSnormal",
            true,
            false
        );

        // Curvature (scalar)
        AddField<Set::Scalar, Set::Hypercube::Cell>(
            ls.LSkappa,
            nullptr,
            1,
            0,
            prefix + "_LSkappa",
            true,
            false
        );

        // Tube data (int flags packed as scalars)
        AddField<Set::Scalar, Set::Hypercube::Cell>(
            ls.TUBE,
            nullptr,
            NarrowBandData::TubeIdx::NCOMP,                
            domain_->nghost,
            prefix + "_TUBE",
            true,
            false
        );
    }


    // Placeholder override functions
    void NarrowBandLevelset::Initialize(int lev) {
        // Initialize global velocity field (will use flow domain in future)
        domain_->ic_Vel->Initialize(lev, domain_->Flowvel);

        // Loop over all level set fields and initialize
        for (auto& ls_ptr : domain_->fields)
        {
            // Get the pointer to the current level set field
            LevelSetField& ls = *ls_ptr;
            const auto& geom = domain_->Geom(lev);

            // If first level, size all vector data
            if (lev == 0) {
                // Get number of objects and resize
                int num_obj = ls.GetNumObjects();
                ls.objects.resize(num_obj);

                // Resize metalevel data
                for (int i = 0; i < num_obj; ++i) {
                    ObjectMeta& meta = ls.objects[lev].meta;

                    meta.id = i;
                    meta.resize_leveldata(domain_->FinestLevel());
                }
            }

            // Initialize LS from IC
            ls.ic_LS->Initialize(lev, ls.LS);

            // Build the tube
            ls.update_Tube(lev, geom, true);

            // Copy LS to LSold
            amrex::MultiFab::Copy(*ls.LSold[lev],
                *ls.LS[lev],
                0,
                0,
                ls.LS[lev]->nComp(),
                ls.LS[lev]->nGrow()
            );

            // Copy domain velocity to level set velocity
            amrex::Copy(
                *ls.LSvel[lev],
                *domain_->Flowvel[lev],
                0, 0, 1, 0
            );

            // Set all real fields to 0
            ls.LSkappa[lev]->setVal(0.0);

            // COmpute normal
            ls.compute_normal(lev, geom);
        }
    }

    void NarrowBandLevelset::TimeStepBegin(Set::Scalar time, int iter) {
        return;
    }

    void NarrowBandLevelset::TimeStepComplete(Set::Scalar time, int lev) {
        return;
    }

    void NarrowBandLevelset::Advance(int lev, Set::Scalar time, Set::Scalar dt) {
        return;
    }

    void NarrowBandLevelset::TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) {
        return;
    }

    void NarrowBandLevelset::Regrid(int lev, Set::Scalar time) {
        return;
    }

}

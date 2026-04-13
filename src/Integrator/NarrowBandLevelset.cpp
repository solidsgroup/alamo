#include <array>

#include "NarrowBandLevelset.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "BC/Constant.H"
#include "BC/Nothing.H"
#include "BC/Expression.H"

// Structures
struct NarrowBandData
{
    // ===============================
    // Owns dense stencil compute data
    // ===============================
    amrex::Box bbox;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;

    amrex::MultiFab real_mf;
    std::array<amrex::MultiFab, AMREX_SPACEDIM> flux_mf;

    amrex::iMultiFab int_mf;

    struct NarrowBandTubeType
    {
        static constexpr int INTERFACE   = 0;  // 0 level set
        static constexpr int INNERBAND   = 1;  // |phi| <  4.0 * min_DX
        static constexpr int INSIDETUBE  = 2;  // |phi| <  6.0 * min_DX
        static constexpr int EDGEPOINT   = 3;  // |phi| >= 6.0 * min_DX & nbr == INSIDETUBE
        static constexpr int OUTSIDETUBE = 4;  // |phi| >= 6.0 * min_DX
    };

    struct NBRealIdx
    {
        static constexpr int LS = 0;

        static constexpr int VEL_START = 1;
        static constexpr int VEL_SIZE  = AMREX_SPACEDIM;

        static constexpr int ACC_START = VEL_START + VEL_SIZE;
        static constexpr int ACC_SIZE  = AMREX_SPACEDIM;

        static constexpr int NORM_START = ACC_START + ACC_SIZE;
        static constexpr int NORM_SIZE  = AMREX_SPACEDIM;

        static constexpr int KAPPA = NORM_START + NORM_SIZE;

        static constexpr int NCOMP = KAPPA + 1;
    };

    struct NBIntIdx
    {
        enum {
            NB_MASK = 0,     // narrowband flag
            TUBE_TYPE,       // uses NarrowBandTubeType
            NCOMP
        };
    };
};

struct ObjectMeta
{
    // =====================
    // Owns bookkeeping only
    // =====================
    int id = -1;
    int material_id = -1;

    bool active = true;  // Narrowband exists. Delete if false.
    bool moving = true;  // Object is moving. Reinitialize if false.
    bool needs_reinit = false;

    amrex::Box bbox;

    Set::Scalar volume = 0.0;
    Set::Scalar centroid[AMREX_SPACEDIM] = {0.0};
};

struct LevelSetObject
{
    // =========================
    // Owns one object’s physics
    // =========================
    ObjectMeta meta;
    NarrowBandData nb;

    void mark_for_reinit() { meta.needs_reinit = true; }
};

struct LevelSetField
{
    // ===============================
    // Owns φ and object decomposition
    // ===============================

    // Full domain level set mfs
    Set::Field<Set::Scalar> LS;
    Set::Field<Set::Scalar> LSold;

    // Full domain geometry mfs
    Set::Field<Set::Scalar> LSvel; // Tube extended velocity
    Set::Field<Set::Scalar> LSnormal;
    Set::Field<Set::Scalar> LSkappa;

    // Full domain int multifabs (for debugging)
    Set::Field<Set::Scalar> TUBE_int;

    std::vector<LevelSetObject> objects;

    // Define IC/BC for this field
    std::unique_ptr<IC::IC<Set::Scalar>> ic_LS;
    std::unique_ptr<BC::BC<Set::Scalar>> bc_LS;

    std::unique_ptr<BC::BC<Set::Scalar>> bc_nothing; // BC::Nothing for all fields that are not LS/LSold

    bool force_reinit_all = false;  // Reinitialize all objects if true
};

struct SimulationDomain
{
    // ========================
    // Owns the grid and fields
    // ========================
    const amrex::Geometry* geom = nullptr;
    const amrex::BoxArray* ba = nullptr;
    const amrex::DistributionMapping* dm = nullptr;
    
    int nghost = 3;

    std::vector<std::unique_ptr<LevelSetField>> fields;

    // Global velocity field from flow domain
    Set::Field<Set::Scalar> velocity;

    // Define function to get global geom info from Integrator
    void DefineFromIntegrator(Integrator::Integrator* integrator, int lev)
        {
            geom = &integrator->Geom(lev);
            ba   = &integrator->boxArray(lev);
            dm   = &integrator->DistributionMap(lev);
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

        // -------------------------
        // Number of level sets
        // -------------------------
        int nls = 1;
        pp_query_default("ls.number_of_levelsets", nls, 1);

        // -------------------------
        // Allocate fields (no geometry yet!)
        // -------------------------
        value.domain_->fields.resize(nls);

        for (int i = 0; i < nls; ++i)
        {
            value.domain_->fields[i] = std::make_unique<LevelSetField>();

            // Define LevelSetField object reference
            LevelSetField& ls = *value.domain_->fields[i];

            // Create per-levelset ParmParse
            std::string prefix = "ls" + std::to_string(i);
            IO::ParmParse pp_ls(prefix.c_str());

            // Initialize pointers to nullptr
            ls.ic_LS = nullptr;
            ls.bc_LS = nullptr;

            ls.bc_nothing = std::make_unique<BC::Nothing>();

            // Define temporary raw pointers for IC/BC selection
            IC::IC<Set::Scalar>* ic_tmp = nullptr;
            BC::BC<Set::Scalar>* bc_tmp = nullptr;

            // Parse and transfer ownership back to smart pointers
            pp_ls.select_default<IC::Constant, IC::Expression>(
                "ic",
                ic_tmp,
                value.geom
            );
            ls.ic_LS = std::unique_ptr<IC::IC<Set::Scalar>>(ic_tmp);

            pp_ls.select_default<BC::Constant, BC::Expression>(
                "bc",
                bc_tmp,
                1
            );
            ls.bc_LS = std::unique_ptr<BC::BC<Set::Scalar>>(bc_tmp);

            // Register fields for this level set
            value.RegisterOneField(ls, prefix);
        }
    }

    // Wrapper function for registering fields
    void NarrowBandLevelset::RegisterOneField(
        LevelSetField& ls,
        const std::string& prefix)
    {
        // ============================
        // Register each field ========
        // ============================

        // Register LS and LSold
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

        // Register geometry fields (nothing BC)
        RegisterNewFab(
            ls.LSvel,
            ls.bc_nothing.get(),
            AMREX_SPACEDIM,
            domain_->nghost,
            prefix + "_LSvel",
            true
        );

        RegisterNewFab(
            ls.LSnormal,
            ls.bc_nothing.get(),
            AMREX_SPACEDIM,
            domain_->nghost,
            prefix + "_LSnormal",
            true
        );

        RegisterNewFab(
            ls.LSkappa,
            ls.bc_nothing.get(),
            1,
            domain_->nghost,
            prefix + "_LSkappa",
            true
        );

        // Register int field
        RegisterNewFab(
            ls.TUBE_int,
            ls.bc_nothing.get(),
            3, // narrowband flag, tube type, and object id
            domain_->nghost,
            prefix + "_TUBE_int",
            true
        );
    }

    // Placeholder override functions
    void NarrowBandLevelset::Initialize(int lev) {
        return;
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

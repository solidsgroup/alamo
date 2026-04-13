#include <array>

#include "NarrowBandLevelset.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "BC/Constant.H"
#include "BC/Nothing.H"
#include "BC/Expression.H"

// Structures

// -------------------------------
// Per-level narrowband container
// -------------------------------
struct NarrowBandLevel
{
    amrex::Box bbox;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;

    amrex::MultiFab real_mf;
    std::array<amrex::MultiFab, AMREX_SPACEDIM> flux_mf;

    amrex::iMultiFab int_mf;
};

// -------------------------------
// Narrowband (AMR-aware)
// -------------------------------
struct NarrowBandData
{
    // One narrowband per AMR level
    amrex::Vector<NarrowBandLevel> leveldata;

    struct NarrowBandTubeType
    {
        static constexpr int INTERFACE   = 0;
        static constexpr int INNERBAND   = 1;
        static constexpr int INSIDETUBE  = 2;
        static constexpr int EDGEPOINT   = 3;
        static constexpr int OUTSIDETUBE = 4;
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
            NB_MASK = 0,
            TUBE_TYPE,
            NCOMP
        };
    };
};

// -------------------------------
// Object metadata
// -------------------------------
struct ObjectMeta
{
    int id = -1;
    int material_id = -1;

    bool active = true;
    bool moving = true;
    bool needs_reinit = false;

    amrex::Box bbox;

    Set::Scalar volume = 0.0;
    Set::Scalar centroid[AMREX_SPACEDIM] = {0.0};
};

// -------------------------------
struct LevelSetObject
{
    ObjectMeta meta;
    NarrowBandData nb;

    void mark_for_reinit() { meta.needs_reinit = true; }
};

// -------------------------------
// Level set field (AMR-ready)
// -------------------------------
struct LevelSetField
{
    Set::Field<Set::Scalar> LS;
    Set::Field<Set::Scalar> LSold;

    Set::Field<Set::Scalar> LSvel;
    Set::Field<Set::Scalar> LSnormal;
    Set::Field<Set::Scalar> LSkappa;

    Set::Field<Set::Scalar> TUBE_int;

    std::vector<LevelSetObject> objects;

    std::unique_ptr<IC::IC<Set::Scalar>> ic_LS;
    std::unique_ptr<BC::BC<Set::Scalar>> bc_LS;
    std::unique_ptr<BC::BC<Set::Scalar>> bc_nothing;

    bool force_reinit_all = false;

    void clip_phi(int lev, const amrex::Geometry& geom) {
        // Define the value to clip to based on 6.0*min_dx
        const Set::Scalar* DX = geom.CellSize();
        Set::Scalar max_phi = 6.0 * std::min({DX[0], DX[1], DX[2]});

        for (amrex::MFIter mfi(*LS[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& gbox = mfi.growntilebox();

            Set::Patch<Set::Scalar> phi_arr = LS.Patch(lev, mfi);

            amrex::ParallelFor(gbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Scalar& phi = phi_arr(i, j, k, 0);

                phi = std::clamp(phi, -max_phi, max_phi);
            });
        }
    }
};

// ======================================
// Simulation Domain (AMR FIXED)
// ======================================
struct SimulationDomain
{
    Integrator::Integrator* integrator = nullptr;

    int nghost = 3;

    std::vector<std::unique_ptr<LevelSetField>> fields;

    Set::Field<Set::Scalar> velocity;

    std::unique_ptr<IC::IC<Set::Scalar>> ic_Vel;
    std::unique_ptr<BC::BC<Set::Scalar>> bc_Vel = std::make_unique<BC::Nothing>();

    void Define(Integrator::Integrator* _integrator)
    {
        integrator = _integrator;
    }

    // -------- AMR-safe accessors --------
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
            ls.bc_nothing = std::make_unique<BC::Nothing>();

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

        // Initialize velocity field (Will pull from flow domain in future)
        IC::IC<Set::Scalar>* ic_vel_tmp = nullptr;

        pp.select_default<IC::Constant, IC::Expression>(
            "velocity.ic",
            ic_vel_tmp,
            value.geom
        );

        value.domain_->ic_Vel =
            std::unique_ptr<IC::IC<Set::Scalar>>(ic_vel_tmp);

        value.RegisterNewFab(value.domain_->velocity, value.domain_->bc_Vel.get(), AMREX_SPACEDIM, 0, "LSvel", true);
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
        // Initialize global velocity field (will use flow domain in future)
        domain_->ic_Vel->Initialize(lev, domain_->velocity);

        // Loop over all level set fields and initialize
        for (auto& ls_ptr : domain_->fields)
        {
            // Get the pointer to the current level set field
            LevelSetField& ls = *ls_ptr;

            // Initialize LS from IC
            if (ls.ic_LS)
            {
                ls.ic_LS->Initialize(lev, ls.LS);
            }
            else
            {
                amrex::Abort("IC for level set is not defined!");
            }

            // Clip phi to |phi| <= 6.0 * min(dx) to track interface
            ls.clip_phi(lev, domain_->Geom(lev));

            amrex::MultiFab::Copy(*ls.LSold[lev],
                *ls.LS[lev],
                0,
                0,
                ls.LS[lev]->nComp(),
                ls.LS[lev]->nGrow()
            );

            // Copy domain velocity to level set velocity
            amrex::MultiFab::Copy(*ls.LSvel[lev], *domain_->velocity[lev], 0, 0, AMREX_SPACEDIM, 0);

            // Set all real fields to 0
            ls.LSnormal[lev]->setVal(0.0);
            ls.LSkappa[lev]->setVal(0.0);

            // Initialize int fields
            ls.TUBE_int[lev]->setVal(0.0, 0, 1, 0); // narrowband flag
            ls.TUBE_int[lev]->setVal(static_cast<Set::Scalar>(NarrowBandData::NarrowBandTubeType::OUTSIDETUBE), 1, 1, 0); // tube type
            ls.TUBE_int[lev]->setVal(-1.0, 2, 1, 0); // object id
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

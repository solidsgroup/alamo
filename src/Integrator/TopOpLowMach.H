#ifndef INTEGRATOR_TOPOPLOWMACH_H
#define INTEGRATOR_TOPOPLOWMACH_H

#include "Integrator/TopOp.H"
#include "Integrator/LowMach.H"

namespace Integrator
{
// Couples Integrator::TopOp (phase-field topology optimization) with
// Integrator::LowMach (reacting low-Mach flow, including its own elastic
// solve) so a structure can be optimized under a load supplied by a fluid.
//
// Both bases derive `virtual public Base::Mechanics<MODEL>`; instantiating
// TopOp on LowMach's own elastic model (NeoHookeanPredeformed) collapses the
// diamond to a single Base::Mechanics subobject -- one model_mf, one
// rhs_mf, one stress_mf/strain_mf, one solver -- following the same pattern
// as SFI<PF> (SFI.H) and ThermoElastic (ThermoElastic.H).
//
// LowMach is not robust against voids in the field it solves the flow
// around. Topology optimization exists to create voids. So this class
// keeps two phase fields:
//   - psi_TO (TopOp's own psi_mf): the topology-optimization density.
//     Develops voids; masks the elastic solve (via a soft ersatz void
//     model blended by psi_TO in UpdateModel -- see there).
//   - eta_LM (this class's own eta_LM_mf): the geometry the fluid sees.
//     In a later commit this is relaxed toward psi_TO by a
//     gradient-limited relaxation (so it develops no interior voids) and
//     pushed into LowMach's component_density_mf for the rigid species,
//     which is where LowMach reconstructs rigid_eta_mf/rigid_species_eta_mf
//     from every step (LowMach::UpdateComponentState).
//
// UpdateModel below reads eta_LM_mf, not rigid_eta_mf, for the interfacial
// pressure traction: eta_LM is the controlled target geometry, while
// rigid_eta_mf can drift slightly from it between TOPOP-phase solves as the
// flow's own advection/diffusion/chemistry perturb component_density_mf.
// The elastic mask (model_mf's stiffness blend) uses psi_TO, not eta_LM, so
// the two fields being slightly out of step is an inherent, accepted
// approximation of this design -- not a bug to chase.
class TopOpLowMach :
    virtual public TopOp<Model::Solid::Finite::NeoHookeanPredeformed>,
    virtual public LowMach
{
    using elastic_model_type = Model::Solid::Finite::NeoHookeanPredeformed;
    using TopOpBase = TopOp<elastic_model_type>;

public:
    TopOpLowMach(): TopOpBase(), LowMach() {}
    TopOpLowMach(IO::ParmParse& pp): TopOpLowMach()
    {
        Parse(*this, pp);
    }

    ~TopOpLowMach()
    {
        delete eta_LM_bc;
        delete ic_eta_LM;
    }

    static void Parse(TopOpLowMach& value, IO::ParmParse& pp)
    {
        // TopOp::ParseCommon must run before LowMach::Parse: LowMach::Parse
        // ends with an unused-input check across the whole (unprefixed)
        // input scope, and would flag every "topop.*" key as unused if they
        // had not already been consumed by the time it runs.
        //
        // ParseCommon (not Parse) is used here so
        // Base::Mechanics<elastic_model_type>::Parse -- which registers
        // disp_mf/rhs_mf/stress_mf/strain_mf/model_mf and must run exactly
        // once -- is not invoked a second time; LowMach::Parse below is the
        // sole caller of it (under its own "elastic." prefix).
        pp.pushPrefix("topop");
        TopOpBase::ParseCommon(value, pp);
        pp.popPrefix();

        // eta_LM: the geometry the fluid sees (see class comment). Parsed
        // and registered here, before LowMach::Parse, for the same
        // unused-input-check-ordering reason as topop.* above. Zero-Neumann
        // like LowMach's own rigid_eta_bc, so physical-wall ghost cells
        // extrapolate from the interior rather than defaulting to 0 (void)
        // -- see the fix note on rigid_eta_bc in LowMach.H.
        value.eta_LM_bc = new BC::Constant::ZeroNeumann(1);
        value.AddField<Set::Scalar, Set::HC::Cell>(
            value.eta_LM_mf, value.eta_LM_bc, 1, 2, "eta_LM", true, false);
        pp.select<IC::Constant, IC::Expression>(
            "eta_lm.ic", value.ic_eta_LM, pp.forward_args(value.geom));

        LowMach::Parse(value, pp);

        if (value.m_type != Base::Mechanics<elastic_model_type>::Type::Disable &&
            !value.elastic.void_model_on)
        {
            Util::Abort(INFO,
                "TopOpLowMach requires elastic.void.model.mu/.kappa (or "
                ".E/.nu) to be set. The topology void is masked by a soft "
                "ersatz model blended with psi_TO in UpdateModel, not by "
                "masking the elastic operator with psi directly -- see the "
                "class comment in TopOpLowMach.H.");
        }
    }

protected:
    void Initialize(int lev) override
    {
        // Diamond ambiguity: both TopOp and LowMach override Initialize
        // (each ultimately reaching the shared Base::Mechanics::Initialize,
        // called twice here -- harmless, it only zeroes disp_mf/rhs_mf).
        LowMach::Initialize(lev);
        TopOpBase::Initialize(lev);
        ic_eta_LM->Initialize(lev, eta_LM_mf);
    }

    void TimeStepBegin(Set::Scalar time, int iter) override
    {
        // Base::Mechanics::TimeStepBegin (reached via LowMach::TimeStepBegin)
        // calls solver.setPsi(psi_mf) whenever psi_on is true, multiplying
        // the elastic operator by psi_TO on top of the ersatz void blend
        // UpdateModel already built into model_mf. That double-masks the
        // void: psi's 1e-8 floor (Operator::Elastic::m_psi_small) would
        // reintroduce, on top of the soft void modulus, exactly the
        // near-singular operator the ersatz model exists to avoid (see the
        // same reasoning LowMach::Parse documents for its own
        // elastic.void_model_on branch). Clear psi_on only around this call;
        // TopOp::Initialize and TagCellsForRefinement need it back on
        // afterward, since it also gates whether psi_mf exists at all.
        bool psi_on_saved = psi_on;
        psi_on = false;
        LowMach::TimeStepBegin(time, iter);
        psi_on = psi_on_saved;
    }

    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override
    {
        // No phase logic yet: psi_TO and eta_LM stay at their ICs until a
        // later commit wires in AdvancePsi / eta_LM relaxation on their own
        // pseudo-time schedule.
        LowMach::Advance(lev, time, dt);
    }

    void TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar a_time, int a_ngrow) override
    {
        // Same harmless double-call to Base::Mechanics::TagCellsForRefinement
        // as Initialize above; TopOpBase's own psi_TO-gradient criterion is
        // gated correctly since psi_on is true here (only cleared inside
        // TimeStepBegin's solve window).
        LowMach::TagCellsForRefinement(lev, a_tags, a_time, a_ngrow);
        TopOpBase::TagCellsForRefinement(lev, a_tags, a_time, a_ngrow);
    }

    // Builds model_mf (TopOp's uniform solid model blended against the soft
    // ersatz void model by psi_TO) and rhs_mf (a constant traction plus,
    // optionally, the local fluid pressure acting normal to the eta_LM
    // interface -- still one-way, not fed back into the momentum equation).
    // A hybrid of TopOp::UpdateModel (the uniform solid model) and
    // LowMach::UpdateModel (the pressure-traction kernel and the void
    // blend, LowMach.cpp:2245-2352) -- see the class comment above for why
    // psi_TO drives the model blend while eta_LM drives the traction.
    void UpdateModel(int /*a_step*/, Set::Scalar /*a_time*/) override
    {
        if (m_type == Base::Mechanics<elastic_model_type>::Type::Disable) return;

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            amrex::Box domain = geom[lev].Domain();
            domain.convert(amrex::IntVect::TheNodeVector());
            const Set::Scalar* DX = geom[lev].CellSize();

            // Ghost fills before this function's grad(eta_LM) RHS forcing
            // and cell-to-node averaging read them -- mirrors
            // LowMach::UpdateModel's preamble (LowMach.cpp:2270-2279).
            eta_LM_bc->define(geom[lev]);
            eta_LM_bc->FillBoundary(*eta_LM_mf[lev], 0, 1, 0.0, 0);
            eta_LM_mf[lev]->FillBoundary(geom[lev].periodicity());
            Util::RealFillBoundary(*psi_mf[lev], geom[lev]);
            pressure_mf[lev]->FillBoundary(geom[lev].periodicity());

            const elastic_model_type solid_model = model;
            const elastic_model_type void_model = elastic.void_model;
            const Set::Scalar traction = elastic.traction;
            const Set::Scalar etacutoff = elastic.etacutoff;
            const bool use_fluid_pressure = elastic.apply_fluid_pressure != 0.0;

            for (amrex::MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box smallbox = mfi.nodaltilebox();
                amrex::Box bx       = mfi.grownnodaltilebox() & domain;

                Set::Patch<elastic_model_type> model_arr = model_mf.Patch(lev, mfi);
                Set::Patch<Set::Vector> rhs = rhs_mf.Patch(lev, mfi);
                Set::Patch<const Set::Scalar> eta_LM = eta_LM_mf.Patch(lev, mfi);
                Set::Patch<const Set::Scalar> psi = psi_mf.Patch(lev, mfi);
                Set::Patch<const Set::Scalar> pressure = pressure_mf.Patch(lev, mfi);

                // Interfacial body force: a constant traction plus,
                // optionally, the local fluid pressure. grad_eta points from
                // fluid into the solid, so a positive (compressive) fluid
                // pressure must push back against that gradient -- hence the
                // minus sign, matching LowMach::UpdateModel.
                amrex::ParallelFor(smallbox, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Vector grad_eta = Numeric::CellGradientOnNode(eta_LM, i, j, k, 0, DX);
                    Set::Scalar eta_node = Numeric::Interpolate::CellToNodeAverage(eta_LM, i, j, k, 0);
                    if (eta_node <= etacutoff) { rhs(i,j,k) = Set::Vector::Zero(); return; }
                    Set::Scalar p_node = use_fluid_pressure
                        ? Numeric::Interpolate::CellToNodeAverage(pressure, i, j, k, 0) : 0.0;
                    rhs(i,j,k) = (traction - p_node) * grad_eta;
                });

                // Soft ersatz void blend, weighted by psi_TO. The clamp is
                // load-bearing, not defensive: NeoHookeanPredeformed's
                // operator* scales F0, so this affine combination is only
                // valid while the two weights sum to exactly 1 (cf.
                // LowMach.cpp:2338-2342, the same reasoning for its
                // rigid-fraction blend).
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar psi_node = Util::Clamp(
                        Numeric::Interpolate::CellToNodeAverage(psi, i, j, k, 0), 0.0, 1.0);
                    model_arr(i,j,k) = solid_model * psi_node + void_model * (1.0 - psi_node);
                });
            }
            Util::RealFillBoundary(*model_mf[lev], geom[lev]);
        }
    }

protected:
    Set::Field<Set::Scalar> eta_LM_mf;
    BC::BC<Set::Scalar>* eta_LM_bc = nullptr;
    IC::IC<Set::Scalar>* ic_eta_LM = nullptr;
};
} // namespace Integrator
#endif

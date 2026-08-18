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
//   - eta_LM (LowMach's rigid_eta_mf, reconstructed each step from
//     component_density_mf): the geometry the fluid sees. In later
//     commits this is relaxed toward psi_TO by a gradient-limited
//     relaxation so it develops no interior voids.
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

    void UpdateModel(int a_step, Set::Scalar a_time) override
    {
        // Placeholder: reuse LowMach's own species-mixture UpdateModel
        // unmodified. A later commit replaces this with the psi_TO/ersatz
        // void blend plus the (traction - pressure)*grad(eta_LM) body
        // force described in the class comment above.
        LowMach::UpdateModel(a_step, a_time);
    }
};
} // namespace Integrator
#endif

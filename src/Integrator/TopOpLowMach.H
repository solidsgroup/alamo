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
// UpdateModel below reads eta_LM_mf, not LowMach's own rigid_eta_mf (the
// reconstruction off component_density_mf), for the interfacial pressure
// traction: eta_LM is the controlled target geometry, while LowMach's
// rigid_eta_mf can drift slightly from it between TOPOP-phase solves as the
// flow's own advection/diffusion/chemistry perturb component_density_mf.
// The elastic mask (model_mf's stiffness blend) uses psi_TO, not eta_LM, so
// the two fields being slightly out of step is an inherent, accepted
// approximation of this design -- not a bug to chase.
//
// A third, optional field, TopOp::rigid_eta_mf, marks geometry that must
// stay solid unconditionally (e.g. a fixed support wall) -- see the class
// comment in TopOp.H. It is never written to by AdvancePsi, RelaxEtaLM, or
// PushEtaToComponentDensity below; wherever this class needs "what's solid
// here", it computes a union of rigid_eta with psi_TO (for the model blend)
// or eta_LM (for the traction and the LowMach density push) *at the point
// of use*, rather than folding rigid_eta into eta_LM_mf itself. Storing the
// union back into eta_LM_mf would expose it to RelaxEtaLM's gradient-based
// relaxation and to the interior-node BC fill in UpdateModel's preamble --
// boundary/interpolation effects that were exactly the mechanism behind the
// original mask leakage this design replaces. Combining only at the point
// of use means rigid_eta itself is never touched by any of that.
//
// Scheduling: one global alamo step is one phase. A cycle is
// topop.n_flow FLOW steps, then 1 TOPOP step, then topop.n_relax RELAX
// steps, repeating (see Phase/PhaseForStep below); TOPOP is placed at the
// start of each cycle (not after FLOW) so it lines up with
// Base::Mechanics's own elastic.interval gate, which fires whenever
// (global step) % elastic.interval == 0 -- Parse requires
// elastic.interval == topop.n_flow + 1 + topop.n_relax so the solve fires
// exactly once per cycle. Consequence: a TOPOP step uses the flow state
// left over from the FLOW phase at the *end of the previous* cycle (fully
// stabilized, for every cycle but the first -- the first TOPOP, at global
// step 0, necessarily uses the IC pressure field, same as how
// Base::Mechanics always performs one solve at step 0 regardless of
// elastic.interval; this is an accepted edge case, not a bug).
// FLOW runs LowMach::Advance/TimeStepComplete normally; RELAX runs
// RelaxEtaLM + PushEtaToComponentDensity and freezes the rest of the flow
// (LowMach::Advance/TimeStepComplete are not called) so eta_LM alone
// moves; TOPOP runs TopOp::AdvancePsi on its own topop.dt pseudo-time
// clock (m_topop_time), not the physical clock, since the flow's dt and
// the phase-field's natural pseudo-timestep are unrelated scales.
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
        value.AddField<Set::Scalar, Set::HC::Cell>(
            value.eta_LM_old_mf, value.eta_LM_bc, 1, 2, "eta_LM_old", false, false);
        pp.select<IC::Constant, IC::Expression>(
            "eta_lm.ic", value.ic_eta_LM, pp.forward_args(value.geom));

        // Phase schedule -- see the class comment above.
        pp.pushPrefix("topop");
        pp_query_default("n_flow", value.topop_n_flow, 0);   // FLOW steps per cycle
        pp_query_default("n_relax", value.topop_n_relax, 0); // RELAX steps per cycle
        // TOPOP phase-field pseudo-dt: parsed only so existing inputs that
        // set it don't trip LowMach::Parse's unused-input check. The actual
        // value used each TOPOP step is overwritten dynamically in
        // TimeStepBegin -- see TopOp::ComputePsiDt/topop.max_dpsidt.
        pp_query_default("dt", value.topop_dt, 0.0);
        pp_query_default("relax.dt", value.topop_relax_dt, 0.0); // RELAX pseudo-dt
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

        const int cycle_length = value.topop_n_flow + 1 + value.topop_n_relax;
        if (value.m_type != Base::Mechanics<elastic_model_type>::Type::Disable &&
            value.m_interval != cycle_length)
        {
            Util::Abort(INFO,
                "TopOpLowMach requires elastic.interval == topop.n_flow + 1 "
                "+ topop.n_relax (the phase cycle length), so the elastic "
                "solve fires exactly once per cycle, at the TOPOP step. Got "
                "elastic.interval=", value.m_interval, " but topop.n_flow=",
                value.topop_n_flow, " + 1 + topop.n_relax=", value.topop_n_relax,
                " = ", cycle_length, ".");
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

    enum class Phase { Flow, TopOpSolve, Relax };

    // See the class comment above for the schedule and why TOPOP sits at
    // the start of each cycle rather than after FLOW.
    Phase PhaseForStep(int step) const
    {
        const int cycle_length = topop_n_flow + 1 + topop_n_relax;
        if (cycle_length <= 0) return Phase::Flow; // degenerate/unconfigured
        const int s = step % cycle_length;
        if (s == 0) return Phase::TopOpSolve;
        if (s <= topop_n_relax) return Phase::Relax;
        return Phase::Flow;
    }

    void TimeStepBegin(Set::Scalar time, int iter) override
    {
        m_phase = PhaseForStep(iter);

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
        //
        // Called every step regardless of phase (UpdateModel is cheap, and
        // the elastic.interval gate inside this call is what actually
        // restricts the solve to the TOPOP step -- see Parse's check that
        // elastic.interval equals the cycle length).
        bool psi_on_saved = psi_on;
        psi_on = false;
        LowMach::TimeStepBegin(time, iter);
        psi_on = psi_on_saved;

        // Dynamic pseudo-dt for the upcoming AdvancePsi call (see
        // TopOp::ComputePsiDt): scans the driving force -- now current,
        // since the elastic solve above just refreshed stress_mf/strain_mf
        // -- across every level and picks a single dt so max(|d(psi)/dt|)
        // stays under topop.max_dpsidt. Only meaningful on the TOPOP step;
        // computed here (once per cycle) rather than per-level in Advance,
        // since it must be uniform across the whole domain.
        if (m_phase == Phase::TopOpSolve)
            topop_dt = TopOpBase::ComputePsiDt(m_topop_time);
    }

    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override
    {
        switch (m_phase)
        {
        case Phase::Flow:
            LowMach::Advance(lev, time, dt);
            break;
        case Phase::TopOpSolve:
            // Reads stress_mf/strain_mf from the solve TimeStepBegin just
            // performed (Evolve() runs TimeStepBegin, then IntegrateVariables
            // to refresh `volume`, then Advance -- exactly the order
            // AdvancePsi's driving force needs).
            TopOpBase::AdvancePsi(lev, m_topop_time, topop_dt);
            break;
        case Phase::Relax:
            RelaxEtaLM(lev, topop_relax_dt);
            PushEtaToComponentDensity(lev);
            break;
        }
    }

    void TimeStepComplete(Set::Scalar time, int iter) override
    {
        // FLOW's own postprocessing (implicit diffusion/phase-change,
        // pressure projection, diagnostics) only runs during FLOW: RELAX's
        // entire point is that everything except eta_LM stays frozen (the
        // "old flow domain" the next FLOW phase restarts from), and running
        // a projection against RELAX's un-backfilled component_density_mf
        // (see PushEtaToComponentDensity) is exactly the instability
        // tests/TopOpLowMachRelax's earlier, pre-phase-counter version hit.
        if (m_phase == Phase::Flow)
            LowMach::TimeStepComplete(time, iter);
        else if (m_phase == Phase::TopOpSolve)
            m_topop_time += topop_dt;
    }

    // d(eta_LM)/dt = (psi_TO - eta_LM) * |grad(eta_LM)|, explicit Euler.
    // The |grad(eta_LM)| factor means eta_LM can only move where it
    // already has an interface: a cell interior to a uniform eta_LM region
    // (grad = 0) stays exactly put regardless of what psi_TO does there,
    // so psi_TO developing an interior void cannot nucleate one in eta_LM
    // -- see tests/TopOpLowMachRelax, which checks exactly this.
    void RelaxEtaLM(int lev, Set::Scalar dt)
    {
        BL_PROFILE("TopOpLowMach::RelaxEtaLM");
        std::swap(eta_LM_old_mf[lev], eta_LM_mf[lev]);
        const Set::Scalar* DX = geom[lev].CellSize();
        amrex::Box domain = geom[lev].Domain();

        for (amrex::MFIter mfi(*eta_LM_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            bx.grow(1);
            bx = bx & domain;
            amrex::Array4<const Set::Scalar> const& eta_old = (*eta_LM_old_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& psi = (*psi_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& eta_new = (*eta_LM_mf[lev]).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar grad_mag = Numeric::Gradient(eta_old, i, j, k, 0, DX).lpNorm<2>();
                Set::Scalar val = eta_old(i, j, k) + dt * (psi(i, j, k) - eta_old(i, j, k)) * grad_mag;
                eta_new(i, j, k) = Util::Clamp(val, 0.0, 1.0);
            });
        }
        eta_LM_bc->define(geom[lev]);
        eta_LM_bc->FillBoundary(*eta_LM_mf[lev], 0, 1, 0.0, 0);
        eta_LM_mf[lev]->FillBoundary(geom[lev].periodicity());
    }

    // Writes union(rigid_eta, eta_LM) into component_density_mf for each
    // rigid_solid species -- the exact inverse of
    // LowMach::UpdateComponentState's reconstruction (LowMach.cpp:630-641:
    // rigid_species_eta = component_density / reference_density). LowMach
    // then rebuilds rigid_eta_mf/rigid_species_eta_mf from this every step
    // with no LowMach changes. rigid_eta is combined in here (not folded
    // into eta_LM_mf) so it never passes through RelaxEtaLM -- see the
    // class comment above.
    //
    // This does not conserve mass: it removes/adds condensed mass without
    // touching any other species, so the volume fractions LowMach's inputs
    // assume sum to 1 will not, wherever the combined field changes.
    // Deliberately deferred -- see the "Mass backfill" deferral in the
    // coupling plan.
    void PushEtaToComponentDensity(int lev)
    {
        BL_PROFILE("TopOpLowMach::PushEtaToComponentDensity");
        const bool use_rigid = rigid_eta_on;
        for (amrex::MFIter mfi(*eta_LM_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.growntilebox(eta_LM_mf[lev]->nGrow());
            amrex::Array4<const Set::Scalar> const& eta_LM = (*eta_LM_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& rigid_eta =
                use_rigid ? (*TopOpBase::rigid_eta_mf[lev]).array(mfi) : amrex::Array4<const Set::Scalar>();
            amrex::Array4<Set::Scalar> const& component_density = (*component_density_mf[lev]).array(mfi);
            for (const int n : rigid_solid_species)
            {
                const Set::Scalar rho_ref = reference_density[n];
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar eta_val = eta_LM(i, j, k);
                    if (use_rigid) eta_val = rigid_eta(i, j, k) + eta_val * (1.0 - rigid_eta(i, j, k));
                    component_density(i, j, k, n) = eta_val * rho_ref;
                });
            }
        }
        component_density_bc->define(geom[lev]);
        component_density_bc->FillBoundary(*component_density_mf[lev], 0, nspecies, 0.0, 0);
        component_density_mf[lev]->FillBoundary(geom[lev].periodicity());
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
            if (rigid_eta_on) Util::RealFillBoundary(*TopOpBase::rigid_eta_mf[lev], geom[lev]);

            const elastic_model_type solid_model = model;
            const elastic_model_type void_model = elastic.void_model;
            const Set::Scalar traction = elastic.traction;
            const Set::Scalar etacutoff = elastic.etacutoff;
            const bool use_fluid_pressure = elastic.apply_fluid_pressure != 0.0;
            const bool use_rigid = rigid_eta_on;

            for (amrex::MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box smallbox = mfi.nodaltilebox();
                amrex::Box bx       = mfi.grownnodaltilebox() & domain;

                Set::Patch<elastic_model_type> model_arr = model_mf.Patch(lev, mfi);
                Set::Patch<Set::Vector> rhs = rhs_mf.Patch(lev, mfi);
                Set::Patch<const Set::Scalar> eta_LM = eta_LM_mf.Patch(lev, mfi);
                Set::Patch<const Set::Scalar> psi = psi_mf.Patch(lev, mfi);
                Set::Patch<const Set::Scalar> pressure = pressure_mf.Patch(lev, mfi);
                amrex::Array4<const Set::Scalar> rigid_eta =
                    use_rigid ? (*TopOpBase::rigid_eta_mf[lev]).array(mfi) : amrex::Array4<const Set::Scalar>();

                // Interfacial body force: a constant traction plus,
                // optionally, the local fluid pressure. grad_eta points from
                // fluid into the solid, so a positive (compressive) fluid
                // pressure must push back against that gradient -- hence the
                // minus sign, matching LowMach::UpdateModel.
                //
                // "eta" here is union(rigid_eta, eta_LM), not eta_LM alone,
                // so a fixed rigid wall gets the same interfacial pressure
                // traction as the design boundary. The union's node value
                // and gradient are formed analytically via the product rule
                // on combined = rigid + eta_LM*(1-rigid) -- this reads both
                // fields directly rather than materializing a combined
                // field, so rigid_eta itself is never written to (see the
                // class comment above).
                amrex::ParallelFor(smallbox, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Vector grad_eta_lm = Numeric::CellGradientOnNode(eta_LM, i, j, k, 0, DX);
                    Set::Scalar eta_lm_node = Numeric::Interpolate::CellToNodeAverage(eta_LM, i, j, k, 0);
                    Set::Vector grad_eta = grad_eta_lm;
                    Set::Scalar eta_node = eta_lm_node;
                    if (use_rigid)
                    {
                        Set::Vector grad_rigid = Numeric::CellGradientOnNode(rigid_eta, i, j, k, 0, DX);
                        Set::Scalar rigid_node = Numeric::Interpolate::CellToNodeAverage(rigid_eta, i, j, k, 0);
                        eta_node = rigid_node + eta_lm_node * (1.0 - rigid_node);
                        grad_eta = grad_rigid * (1.0 - eta_lm_node) + grad_eta_lm * (1.0 - rigid_node);
                    }
                    eta_node = Util::Clamp(eta_node, 0.0, 1.0);
                    if (eta_node <= etacutoff) { rhs(i,j,k) = Set::Vector::Zero(); return; }
                    Set::Scalar p_node = use_fluid_pressure
                        ? Numeric::Interpolate::CellToNodeAverage(pressure, i, j, k, 0) : 0.0;
                    rhs(i,j,k) = (traction - p_node) * grad_eta;
                });

                // Soft ersatz void blend, weighted by union(rigid_eta,
                // psi_TO) -- a fixed rigid region reads as solid regardless
                // of what psi_TO does underneath it (see the class comment
                // above). The clamp is load-bearing, not defensive:
                // NeoHookeanPredeformed's operator* scales F0, so this
                // affine combination is only valid while the two weights
                // sum to exactly 1 (cf. LowMach.cpp:2338-2342, the same
                // reasoning for its rigid-fraction blend).
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar psi_node = Util::Clamp(
                        Numeric::Interpolate::CellToNodeAverage(psi, i, j, k, 0), 0.0, 1.0);
                    if (use_rigid)
                    {
                        Set::Scalar rigid_node = Util::Clamp(
                            Numeric::Interpolate::CellToNodeAverage(rigid_eta, i, j, k, 0), 0.0, 1.0);
                        psi_node = Util::Clamp(rigid_node + psi_node * (1.0 - rigid_node), 0.0, 1.0);
                    }
                    model_arr(i,j,k) = solid_model * psi_node + void_model * (1.0 - psi_node);
                });
            }
            Util::RealFillBoundary(*model_mf[lev], geom[lev]);
        }
    }

protected:
    Set::Field<Set::Scalar> eta_LM_mf;
    Set::Field<Set::Scalar> eta_LM_old_mf;
    BC::BC<Set::Scalar>* eta_LM_bc = nullptr;
    IC::IC<Set::Scalar>* ic_eta_LM = nullptr;

    int topop_n_flow = 0;
    int topop_n_relax = 0;
    Set::Scalar topop_dt = 0.0;
    Set::Scalar topop_relax_dt = 0.0;
    Set::Scalar m_topop_time = 0.0;
    Phase m_phase = Phase::Flow;
};
} // namespace Integrator
#endif

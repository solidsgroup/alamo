#ifndef INTEGRATOR_TOPOP_H
#define INTEGRATOR_TOPOP_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <cmath>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "IO/ParmParse.H"
#include "Integrator/Base/Mechanics.H"


#include "IC/IC.H"
#include "BC/BC.H"
#include "BC/Operator/Elastic/Constant.H"
#include "BC/Operator/Elastic/TensionTest.H"
#include "BC/Operator/Elastic/Expression.H"

#include "IC/Ellipse.H"
#include "IC/Voronoi.H"
#include "IC/Constant.H"
#include "IC/BMP.H"
#include "IC/Expression.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"

#include "Model/Solid/Solid.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"

#include "Operator/Operator.H"


namespace Integrator
{
template<class MODEL>
class TopOp: virtual public Base::Mechanics<MODEL>
{
public:

    TopOp(): Base::Mechanics<MODEL>() {}
    TopOp(IO::ParmParse& pp): Base::Mechanics<MODEL>()
    {
        Parse(*this, pp);
    }

    ~TopOp()
    {
        delete ic_psi;
        delete bc_psi;
    }

    // The mechanics integrator manages the solution of an elastic
    // solve using the MLMG solver.
    static void Parse(TopOp& value, IO::ParmParse& pp)
    {
        Base::Mechanics<MODEL>::Parse(value, pp);
        ParseCommon(value, pp);
    }

    // Everything TopOp needs beyond the elastic solve itself (the psi field,
    // the phase-field parameters, the volume constraint). Split out from
    // Parse so a coupled integrator (e.g. TopOpLowMach) can call
    // Base::Mechanics<MODEL>::Parse exactly once -- from whichever base owns
    // it -- and then bring in this half under its own prefix.
    static void ParseCommon(TopOp& value, IO::ParmParse& pp)
    {
        pp.queryclass<MODEL>("model", value.model);

        // Read in IC for psi
        if (pp.contains("psi.ic.type")) 
        {
            value.psi_on = true;
            value.bc_psi = new BC::Constant(1, pp, "psi.bc");
            value.RegisterNewFab(value.psi_mf, value.bc_psi, 1, 2, "psi", true);
            value.RegisterNewFab(value.psi_old_mf, value.bc_psi, 1, 2, "psiold", false);
        }

        // Initial condition for psi field
        pp.select_default<IC::Constant::One,IC::Ellipse,IC::Constant,IC::Expression>("psi.ic",value.ic_psi,pp.forward_args(value.geom));

        // Target maximum |d(psi)/dt|, opt-in via presence of this key.
        // When set, standalone TopOp::Advance drives AdvancePsi with a
        // dynamically computed pseudo-dt (see ComputePsiDt) instead of the
        // physical dt, picking a single domain-wide dt each step so this
        // bound holds everywhere, uniformly in space (never per-cell).
        // Guards against the checkerboarding a fixed dt produces once the
        // interface's length scale shrinks enough that the reaction/
        // diffusion driving force outruns it. Off by default so existing
        // TopOp inputs/references (e.g. tests/TopOp) are unaffected;
        // TopOpLowMach always computes it, regardless of this flag, since
        // its TOPOP phase never shared a timescale with the physical dt to
        // begin with.
        value.m_dynamic_dt = pp.contains("max_dpsidt");
        pp_query_default("max_dpsidt", value.m_max_dpsidt, 0.05);

        pp_query_default("eta_ref_threshold", value.m_eta_ref_threshold, 0.01); // Refinement threshold based on eta
        pp_query_default("alpha", value.alpha, 1.0); // :math:`\alpha` parameter
        pp_query_default("beta", value.beta, 1.0); // :math:`\beta` parameter
        pp_query_default("gamma", value.gamma, 1.0); // :math:`\gamma` parameter
        pp_queryclass("L", value.L); // Mobility
        if (pp.contains("volume0frac"))
        {
            if (pp.contains("volume0")) Util::Abort(INFO, "Cannot specify volume0frac and volume0");
            Set::Scalar volumefrac;
            pp_query_default("volume0frac", volumefrac, 0.5); // Prescribed volume fraction
            value.volume0 = volumefrac *
                AMREX_D_TERM((value.geom[0].ProbHi()[0] - value.geom[0].ProbLo()[0]),
                    *(value.geom[0].ProbHi()[1] - value.geom[0].ProbLo()[1]),
                    *(value.geom[0].ProbHi()[2] - value.geom[0].ProbLo()[2]));
        }
        else
            pp_query_default("volume0", value.volume0, 0.5); // Prescribed total vlume

        pp_queryclass("lambda", value.lambda); // Lagrange multiplier (can be interplated)


        value.RegisterIntegratedVariable(&value.volume, "volume");
        value.RegisterIntegratedVariable(&value.w_chem_potential, "chem_potential");
        value.RegisterIntegratedVariable(&value.w_bndry, "bndry");
        value.RegisterIntegratedVariable(&value.w_elastic, "elastic");
    }

    void Initialize(int lev) override
    {
        Base::Mechanics<MODEL>::Initialize(lev);
        if (psi_on) ic_psi->Initialize(lev, psi_mf);
        if (psi_on) ic_psi->Initialize(lev, psi_old_mf);
    }

    virtual void UpdateModel(int a_step, Set::Scalar /*a_time*/) override
    {
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;

        if (a_step > 0) return;

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            model_mf[lev]->setVal(model);
            Util::RealFillBoundary(*model_mf[lev], geom[lev]);
            Util::RealFillBoundary(*psi_mf[lev], geom[lev]);
            Util::RealFillBoundary(*psi_old_mf[lev], geom[lev]);
        }

    }


    // Computes the driving force once per (global, not per-level) step, via
    // TimeStepBegin below -- see the class comment on Advance for why this
    // must not be recomputed per level.
    void TimeStepBegin(Set::Scalar time, int iter) override
    {
        Base::Mechanics<MODEL>::TimeStepBegin(time, iter);
        if (psi_on && m_dynamic_dt) ComputePsiDt(time);
    }

    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override
    {
        BL_PROFILE("TopOp::Advance");
        Base::Mechanics<MODEL>::Advance(lev, time, dt);
        // Dynamic mode (topop.max_dpsidt set): drive AdvancePsi with the
        // computed pseudo-dt (see ComputePsiDt) instead of the physical dt
        // argument -- psi's stability limit and the flow/mechanics
        // timestep are unrelated scales. Otherwise, unchanged legacy
        // behavior: psi advances on the physical dt.
        AdvancePsi(lev, time, m_dynamic_dt ? m_psi_dt : dt);
    }

    // Scans the driving force (the same formula AdvancePsi applies) across
    // every level and cell and returns a single, domain-wide pseudo-dt such
    // that max(|L * driving_force|) * dt == max_dpsidt. Must be called once
    // per step, before any level's psi has been advanced -- AdvancePsi is
    // applied independently per level (AMR subcycling), and computing this
    // per-level would mix already-updated psi from earlier levels into the
    // scan for later ones. TimeStepBegin above is the call site for
    // standalone TopOp; TopOpLowMach calls it directly from its own
    // TimeStepBegin, once per cycle, on the TOPOP step.
    Set::Scalar ComputePsiDt(Set::Scalar time)
    {
        Set::Scalar Lnow = L(time);
        Set::Scalar lambdaT = lambda(time);
        Set::Scalar alpha_ = alpha, beta_ = beta, gamma_ = gamma;
        Set::Scalar vol = volume, vol0 = volume0;

        Set::Scalar max_rate = 0.0;
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            const Set::Scalar* DX = geom[lev].CellSize();
            amrex::Box domain = geom[lev].Domain();
            for (amrex::MFIter mfi(*psi_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.tilebox();
                bx.grow(1);
                bx = bx & domain;
                amrex::Array4<const Set::Matrix> const& sig = (*stress_mf[lev]).array(mfi);
                amrex::Array4<const Set::Matrix> const& eps = (*strain_mf[lev]).array(mfi);
                amrex::Array4<const MODEL> const& model_arr = (*model_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const& psi = (*psi_mf[lev]).array(mfi);

                amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
                amrex::ReduceData<Set::Scalar> reduce_data(reduce_op);
                using ReduceTuple = typename decltype(reduce_data)::Type;
                reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
                {
                    Set::Scalar driving_force = 0.0;
                    driving_force += alpha_ * 2.0 * psi(i, j, k) * (2.0 * psi(i, j, k) * psi(i, j, k) - 3.0 * psi(i, j, k) + 1.0);
                    driving_force += -beta_ * Numeric::Laplacian(psi, i, j, k, 0, DX);
                    driving_force += -gamma_ * StrainEnergy(model_arr, sig, eps, i, j, k) * psi(i, j, k);
                    driving_force += lambdaT * (vol - vol0);
                    return { std::abs(driving_force) };
                });
                ReduceTuple hv = reduce_data.value();
                max_rate = std::max(max_rate, amrex::get<0>(hv));
            }
        }
        amrex::ParallelDescriptor::ReduceRealMax(max_rate);

        Set::Scalar rate = Lnow * max_rate;
        // rate ~ 0: the driving force is uniformly zero everywhere (e.g.
        // L.str = 0, a frozen psi, or a fully converged field) -- psi will
        // not move this step regardless of dt, so leave the last computed
        // dt in place (or seed it with max_dpsidt on the very first call).
        if (rate > 1.0e-12) m_psi_dt = m_max_dpsidt / rate;
        else if (std::isnan(m_psi_dt)) m_psi_dt = m_max_dpsidt;
        return m_psi_dt;
    }

    // The phase-field (psi) update, factored out of Advance so a coupled
    // integrator can drive it on its own pseudo-time/dt rather than the
    // physical clock (see TopOpLowMach).
    void AdvancePsi(int lev, Set::Scalar time, Set::Scalar dt)
    {
        BL_PROFILE("TopOp::AdvancePsi");
        std::swap(psi_old_mf[lev], psi_mf[lev]);
        const Set::Scalar* DX = geom[lev].CellSize();
        amrex::Box domain = geom[lev].Domain();

        Set::Scalar Lnow = L(time);
        Set::Scalar lambdaT = lambda(time);

        for (amrex::MFIter mfi(*psi_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            bx.grow(1);
            bx = bx & domain;
            amrex::Array4<const Set::Matrix> const& sig = (*stress_mf[lev]).array(mfi);
            amrex::Array4<const Set::Matrix> const& eps = (*strain_mf[lev]).array(mfi);
            amrex::Array4<const MODEL> const& model_arr = (*model_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& psi = (*psi_old_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& psinew = (*psi_mf[lev]).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar driving_force = 0.0;

                driving_force += alpha * 2.0 * psi(i, j, k) * (2.0 * psi(i, j, k) * psi(i, j, k) - 3.0 * psi(i, j, k) + 1.0);
                driving_force += -beta * Numeric::Laplacian(psi, i, j, k, 0, DX);

                driving_force += -gamma * StrainEnergy(model_arr, sig, eps, i, j, k) * psi(i, j, k);

                driving_force += lambdaT * (volume - volume0);

                psinew(i, j, k) = psi(i, j, k) - Lnow * dt * driving_force;
                if (psinew(i, j, k) < 0.0) psinew(i, j, k) = 0.0;
                if (psinew(i, j, k) > 1.0) psinew(i, j, k) = 1.0;
            });
        }
    }

    // Cell-averaged strain energy at (i,j,k), used both for the psi driving
    // force and for the integrated w_elastic diagnostic.
    //
    // For small-strain models (kinvar != F), stress/strain are the Cauchy
    // stress / infinitesimal strain and 0.5*(sig^T*eps).trace() is twice the
    // strain energy density (the pre-existing behavior, unchanged here).
    //
    // For finite-strain models (kinvar == F), Base::Mechanics stores
    // strain_mf = F (the deformation gradient, = I at zero displacement) and
    // stress_mf = DW(F) (the first Piola-Kirchhoff stress). 0.5*(sig^T*eps)
    // then evaluates to a spurious nonzero trace(P) even with no
    // displacement, so use the model's own strain energy density W(F)
    // instead. This branches on MODEL::kinvar, a static constexpr/const, so
    // it compiles away to whichever branch is live for MODEL.
    AMREX_FORCE_INLINE static AMREX_GPU_DEVICE
    Set::Scalar StrainEnergy(amrex::Array4<const MODEL> const& model_arr,
        amrex::Array4<const Set::Matrix> const& sig,
        amrex::Array4<const Set::Matrix> const& eps,
        int i, int j, int k)
    {
        if constexpr (MODEL::kinvar == Model::Solid::KinematicVariable::F)
        {
            Set::Matrix F_avg = Numeric::Interpolate::NodeToCellAverage(eps, i, j, k, 0);
            return model_arr(i, j, k).W(F_avg);
        }
        else
        {
            Set::Matrix sig_avg = Numeric::Interpolate::NodeToCellAverage(sig, i, j, k, 0);
            Set::Matrix eps_avg = Numeric::Interpolate::NodeToCellAverage(eps, i, j, k, 0);
            return 0.5 * (sig_avg.transpose() * eps_avg).trace();
        }
    }

    void Integrate(int amrlev, Set::Scalar time, int step,
        const amrex::MFIter& mfi, const amrex::Box& box) override
    {
        BL_PROFILE("TopOp::Integrate");
        Base::Mechanics<MODEL>::Integrate(amrlev, time, step, mfi, box);

        const amrex::Real* DX = geom[amrlev].CellSize();
        Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);

        amrex::Array4<amrex::Real> const& psi = (*psi_mf[amrlev]).array(mfi);
        amrex::Array4<const Set::Matrix> const& sig = (*stress_mf[amrlev]).array(mfi);
        amrex::Array4<const Set::Matrix> const& eps = (*strain_mf[amrlev]).array(mfi);
        amrex::Array4<const MODEL> const& model_arr = (*model_mf[amrlev]).array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            volume += psi(i, j, k, 0) * dv;
            w_chem_potential += alpha * psi(i, j, k, 0) * psi(i, j, k, 0) * (1. - psi(i, j, k, 0) * psi(i, j, k, 0)) * dv;
            w_bndry += beta * 0.5 * Numeric::Gradient(psi, i, j, k, 0, DX).squaredNorm() * dv;

            w_elastic += gamma * StrainEnergy(model_arr, sig, eps, i, j, k) * psi(i, j, k) * dv;
        });
    }


    void TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar a_time, int a_ngrow) override
    {
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;
        Base::Mechanics<MODEL>::TagCellsForRefinement(lev, a_tags, a_time, a_ngrow);

        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();
        for (amrex::MFIter mfi(*model_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            amrex::Array4<char> const& tags = a_tags.array(mfi);
            if (psi_on)
            {
                amrex::Array4<Set::Scalar> const& psi = psi_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    auto sten = Numeric::GetStencil(i, j, k, bx);
                    {
                        Set::Vector gradpsi = Numeric::Gradient(psi, i, j, k, 0, DX.data(), sten);
                        if (gradpsi.lpNorm<2>() * DXnorm > m_eta_ref_threshold)
                            tags(i, j, k) = amrex::TagBox::SET;
                    }
                });
            }
        }
    }


protected:
    MODEL model;
    IC::IC<Set::Scalar>* ic_psi = nullptr;
    BC::BC<Set::Scalar>* bc_psi = nullptr;
    Set::Scalar m_eta_ref_threshold = NAN;
    Set::Field<Set::Scalar> psi_old_mf;
    bool m_dynamic_dt = false;
    Set::Scalar m_max_dpsidt = NAN;
    Set::Scalar m_psi_dt = NAN;


    using Base::Mechanics<MODEL>::m_type;
    using Base::Mechanics<MODEL>::finest_level;
    using Base::Mechanics<MODEL>::geom;
    using Base::Mechanics<MODEL>::model_mf;
    using Base::Mechanics<MODEL>::psi_mf;
    using Base::Mechanics<MODEL>::psi_on;
    using Base::Mechanics<MODEL>::stress_mf;
    using Base::Mechanics<MODEL>::strain_mf;


    Set::Scalar alpha = NAN;
    Set::Scalar beta = NAN;
    Set::Scalar gamma = NAN;
    //Set::Scalar L = 1.0;
    Numeric::Interpolator::Linear<Set::Scalar> L;
    Set::Scalar volume0 = NAN;
    //Set::Scalar lambda = 1.0;
    Numeric::Interpolator::Linear<Set::Scalar> lambda;

    Set::Scalar volume = NAN;
    Set::Scalar w_chem_potential = NAN;
    Set::Scalar w_bndry = NAN;
    Set::Scalar w_elastic = NAN;
};










} // namespace Integrator
#endif

#ifndef SOLVER_NONLOCAL_LINEAR
#define SOLVER_NONLOCAL_LINEAR
#include "Operator/Operator.H"
#include <AMReX_MLMG.H>
#include <cmath>

namespace Solver
{
namespace Nonlocal
{
/// \brief Multigrid Linear solver for multicomponent, multi-level operators
/// 
/// This class is a thin wrapper for the `amrex::MLMG` solver.
/// It exists to set a range of default MLMG settings automatically, for instance,
/// `setCFStrategy`, which may not be obvious to the user.
///
/// It also exists as a compatibility layer so that future fixes for compatibility
/// with AMReX can be implemented here.
class Linear  // : public amrex::MLMG
{
public:

    Linear() //Operator::Operator<Grid::Node>& a_lp) : MLMG(a_lp), linop(a_lp)
    {
    }

    Linear(Operator::Operator<Grid::Node>& a_lp)
    {
        this->Define(a_lp);
    }

    ~Linear()
    {
        if (m_defined) Clear();
    }

    void Define(Operator::Operator<Grid::Node>& a_lp)
    {
        if (m_defined) Util::Abort(INFO, "Solver cannot be re-defined");
        this->linop = &a_lp;
        this->mlmg = new amrex::MLMG(a_lp);
        m_defined = true;
        PrepareMLMG(*mlmg);
    }
    void Clear()
    {
        if (!m_defined) Util::Abort(INFO, "Solver cannot be cleared if not defined");
        this->linop = nullptr;
        if (this->mlmg) delete this->mlmg;
        this->mlmg = nullptr;
        m_defined = false;
    }

    // Iteration/residual diagnostics (C7). Purely additive read-outs of numbers
    // amrex::MLMG already computes internally but that nothing in this repo
    // previously consumed -- iteration counts have historically been recovered
    // by regexing stdout instead.
    [[nodiscard]] int getNumIters() const { return mlmg->getNumIters(); }
    [[nodiscard]] Set::Scalar getFinalResidual() const { return mlmg->getFinalResidual(); }
    [[nodiscard]] Set::Scalar getInitResidual() const { return mlmg->getInitResidual(); }
    [[nodiscard]] Set::Scalar getInitRHS() const { return mlmg->getInitRHS(); }
    [[nodiscard]] const amrex::Vector<Set::Scalar>& getResidualHistory() const { return mlmg->getResidualHistory(); }

    // Geometric-mean convergence factor over the residual history. This is the
    // number that actually distinguishes "slow" from "broken" -- an iteration
    // count alone conflates the two.
    [[nodiscard]] Set::Scalar getConvergenceFactor() const
    {
        const amrex::Vector<Set::Scalar>& hist = mlmg->getResidualHistory();
        if (hist.size() < 2) return NAN;
        Set::Scalar logsum = 0.0;
        int n = 0;
        for (int i = 1; i < (int)hist.size(); i++)
        {
            if (hist[i - 1] <= 0.0 || hist[i] <= 0.0) continue;
            logsum += std::log(hist[i] / hist[i - 1]);
            n++;
        }
        if (n == 0) return NAN;
        return std::exp(logsum / n);
    }

    enum class SolveStatus { Converged, Stalled, Diverged, Unknown };

    [[nodiscard]] SolveStatus getSolveStatus() const { return m_last_status; }

    Set::Scalar solveaffine(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_sol,
        amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_rhs,
        Real a_tol_rel, Real a_tol_abs, bool copyrhs = false,
        const char* checkpoint_file = nullptr)
    {
        if (!m_defined) Util::Abort(INFO, "Solver not defined");
        amrex::Vector<amrex::MultiFab*> rhs_tmp(a_rhs.size());
        amrex::Vector<amrex::MultiFab*> zero_tmp(a_rhs.size());
        for (int i = 0; i < rhs_tmp.size(); i++)
        {
            rhs_tmp[i] = new amrex::MultiFab(a_rhs[i]->boxArray(), a_rhs[i]->DistributionMap(), a_rhs[i]->nComp(), a_rhs[i]->nGrow());
            zero_tmp[i] = new amrex::MultiFab(a_rhs[i]->boxArray(), a_rhs[i]->DistributionMap(), a_rhs[i]->nComp(), a_rhs[i]->nGrow());
            rhs_tmp[i]->setVal(0.0);
            zero_tmp[i]->setVal(0.0);
            Util::Message(INFO, rhs_tmp[i]->norm0());
        }

        linop->SetHomogeneous(false);
        mlmg->apply(rhs_tmp, zero_tmp);

        for (int lev = 0; lev < rhs_tmp.size(); lev++)
        {
            amrex::Box domain = linop->Geom(lev).Domain();
            domain.convert(amrex::IntVect::TheNodeVector());
            const Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);
            for (MFIter mfi(*rhs_tmp[lev], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
            {
                amrex::Box bx = mfi.growntilebox(rhs_tmp[lev]->nGrow());
                bx = bx & domain;
                amrex::Array4<amrex::Real> const& rhstmp = rhs_tmp[lev]->array(mfi);
                for (int n = 0; n < rhs_tmp[lev]->nComp(); n++)
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        bool    AMREX_D_DECL(xmin = (i == lo.x), ymin = (j == lo.y), zmin = (k == lo.z)),
                            AMREX_D_DECL(xmax = (i == hi.x), ymax = (j == hi.y), zmax = (k == hi.z));
                        if (AMREX_D_TERM(xmax || xmin, || ymax || ymin, || zmax || zmin))
                            rhstmp(i, j, k, n) = 0.0;
                        else
                            rhstmp(i, j, k, n) *= -1.0;
                    });
                }
            }
            Util::Message(INFO, rhs_tmp[lev]->norm0());
            rhs_tmp[lev]->setMultiGhost(true);
            rhs_tmp[lev]->FillBoundaryAndSync(linop->Geom(lev).periodicity());
        }

        for (int lev = 0; lev < rhs_tmp.size(); lev++)
        {
            Util::Message(INFO, rhs_tmp[lev]->norm0());
            amrex::Add(*rhs_tmp[lev], *a_rhs[lev], 0, 0, rhs_tmp[lev]->nComp(), rhs_tmp[lev]->nGrow());
            if (copyrhs)
                amrex::Copy(*a_rhs[lev], *rhs_tmp[lev], 0, 0, rhs_tmp[lev]->nComp(), rhs_tmp[lev]->nGrow());
            Util::Message(INFO, rhs_tmp[lev]->norm0());
        }

        linop->SetHomogeneous(true);
        PrepareMLMG(*mlmg);
        Set::Scalar retval = NAN;
        try
        {
            retval = mlmg->solve(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(rhs_tmp), a_tol_rel, a_tol_abs, checkpoint_file);
            reportSolverStats(false, nullptr);
        }
        catch (const std::exception& e)
        {
            reportSolverStats(true, e.what());
            if (m_dump_on_fail || m_dump_always) dumpOnConvergenceFail(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(rhs_tmp));
            if (m_abort_on_fail) Util::Abort(INFO, e.what());
        }
        if (m_dump_always && m_last_status == SolveStatus::Converged)
            dumpOnConvergenceFail(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(rhs_tmp));
        if (a_sol[0]->contains_nan())
        {
            dumpOnConvergenceFail(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(rhs_tmp));
            Util::Abort(INFO);
        }
        runDiagnosticChecks(a_sol);

        return retval;
    };

    Set::Scalar solve(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_sol,
        amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_rhs,
        Real a_tol_rel, Real a_tol_abs, const char* checkpoint_file = nullptr)
    {
        PrepareMLMG(*mlmg);
        Set::Scalar retval = NAN;
        try
        {
            retval = mlmg->solve(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(a_rhs), a_tol_rel, a_tol_abs, checkpoint_file);
            reportSolverStats(false, nullptr);
        }
        catch (const std::exception& e)
        {
            reportSolverStats(true, e.what());
            if (m_dump_on_fail || m_dump_always) dumpOnConvergenceFail(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(a_rhs));
            if (m_abort_on_fail) Util::Abort(INFO, e.what());
        }
        if (m_dump_always && m_last_status == SolveStatus::Converged)
            dumpOnConvergenceFail(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(a_rhs));
        runDiagnosticChecks(a_sol);
        return retval;
    };
    Set::Scalar solve(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_sol,
        amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_rhs)
    {
        PrepareMLMG(*mlmg);
        Set::Scalar retval = NAN;
        try
        {
            retval = mlmg->solve(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(a_rhs), tol_rel, tol_abs);
            reportSolverStats(false, nullptr);
        }
        catch (const std::exception& e)
        {
            reportSolverStats(true, e.what());
            if (m_dump_on_fail || m_dump_always) dumpOnConvergenceFail(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(a_rhs));
            if (m_abort_on_fail) Util::Abort(INFO, e.what());
        }
        if (m_dump_always && m_last_status == SolveStatus::Converged)
            dumpOnConvergenceFail(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(a_rhs));
        if (a_sol[0]->contains_nan())
        {
            dumpOnConvergenceFail(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(a_rhs));
            Util::Abort(INFO);
        }
        runDiagnosticChecks(a_sol);
        return retval;
    };
    void apply(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_rhs,
        amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_sol)
    {
        PrepareMLMG(*mlmg);
        mlmg->apply(GetVecOfPtrs(a_rhs), GetVecOfPtrs(a_sol));
    };



    void setMaxIter(const int a_max_iter) { max_iter = a_max_iter; }
    void setBottomMaxIter(const int a_bottom_max_iter) { bottom_max_iter = a_bottom_max_iter; }
    void setMaxFmgIter(const int a_max_fmg_iter) { max_fmg_iter = a_max_fmg_iter; }
    void setFixedIter(const int a_fixed_iter) { fixed_iter = a_fixed_iter; }
    void setVerbose(const int a_verbose) { verbose = a_verbose; }
    void setPreSmooth(const int a_pre_smooth) { pre_smooth = a_pre_smooth; }
    void setPostSmooth(const int a_post_smooth) { post_smooth = a_post_smooth; }

    void dumpOnConvergenceFail(const amrex::Vector<amrex::MultiFab*>& a_sol_mf,
        const amrex::Vector<amrex::MultiFab const*>& a_rhs_mf)
    {
        int nlevs = a_sol_mf.size();
        int ncomps = a_sol_mf[0]->nComp();

        amrex::Vector<amrex::Geometry> geom;
        amrex::Vector<int> iter;
        amrex::Vector<amrex::IntVect> refratio;
        amrex::Vector<std::string> names;
        for (int i = 0; i < nlevs; i++)
        {
            geom.push_back(linop->Geom(i));
            iter.push_back(0);
            if (i > 0) refratio.push_back(amrex::IntVect(2));
        }
        for (int n = 0; n < ncomps; n++)
        {
            names.push_back("var" + std::to_string(n));
        }

        std::string outputdir = Util::GetFileName();
        WriteMultiLevelPlotfile(outputdir + "/mlmg_sol", nlevs,
            amrex::GetVecOfConstPtrs(a_sol_mf),
            names, geom, 0, iter, refratio);
        WriteMultiLevelPlotfile(outputdir + "/mlmg_rhs", nlevs,
            amrex::GetVecOfConstPtrs(a_rhs_mf),
            names, geom, 0, iter, refratio);

        Set::Field<Set::Scalar> res_mf(nlevs);
        for (int lev = 0; lev < nlevs; lev++)
        {
            res_mf.Define(lev, a_sol_mf[lev]->boxArray(), a_sol_mf[lev]->DistributionMap(),
                ncomps, a_sol_mf[lev]->nGrow());
        }

        mlmg->compResidual(amrex::GetVecOfPtrs(res_mf), a_sol_mf, a_rhs_mf);

        WriteMultiLevelPlotfile(outputdir + "/mlmg_res", nlevs,
            amrex::GetVecOfConstPtrs(res_mf),
            names, geom, 0, iter, refratio);

    }

    // A5: is the assembled operator symmetric?
    //
    // Tests x.(Ay) - y.(Ax) on random vectors through this class's own
    // apply() (note apply()'s argument order is (out, in)), on the actual
    // configuration -- real BCs, real contrast, real AMR levels -- rather
    // than a synthetic operator.
    //
    // Domain-boundary nodes are zeroed out of both x and y in every
    // component before applying A. Displacement-type BC rows are *replaced*
    // (BC::Operator::Elastic::Constant::set, f(i)=u(i)) while neighboring
    // interior rows keep their original coupling to that boundary node --
    // that row-replacement is expected, known-asymmetric bookkeeping, not a
    // property of the interior discretization. Masking it out isolates the
    // question this test actually exists to answer: is the interior/traction
    // stencil symmetric. (SetHomogeneous is not touched here: m_homogeneous
    // on Operator::Elastic is set but never read anywhere in Fapply, so
    // Fapply is already a purely linear map independent of that flag --
    // see improve-MLMG-solver-ideas.md R5.)
    //
    // Separately: this only exercises Fapply, not Fsmooth. The
    // relaxCoarseFineGhostRows() ghost-row zeroing in Operator::Fsmooth
    // (active on the non-conservative path) is invisible to this test --
    // a clean result here does not exonerate that suspect.
    //
    // Checks single-level dot products per AMR level (no composite
    // correction for overlapping fine/coarse regions) -- adequate for a
    // diagnostic, not for a quantitative norm.
    Set::Scalar checkSymmetry(const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_template,
        int n_samples = 3)
    {
        if (!m_defined) Util::Abort(INFO, "Solver not defined");
        int nlevs = a_template.size();
        Set::Scalar max_relative_asymmetry = 0.0;

        for (int s = 0; s < n_samples; s++)
        {
            amrex::Vector<std::unique_ptr<amrex::MultiFab> > x(nlevs), y(nlevs), Ax(nlevs), Ay(nlevs);
            for (int lev = 0; lev < nlevs; lev++)
            {
                const amrex::BoxArray& ba = a_template[lev]->boxArray();
                const amrex::DistributionMapping& dm = a_template[lev]->DistributionMap();
                int ncomp = a_template[lev]->nComp();
                amrex::IntVect ng = a_template[lev]->nGrowVect();
                x[lev] = std::make_unique<amrex::MultiFab>(ba, dm, ncomp, ng);
                y[lev] = std::make_unique<amrex::MultiFab>(ba, dm, ncomp, ng);
                Ax[lev] = std::make_unique<amrex::MultiFab>(ba, dm, ncomp, ng);
                Ay[lev] = std::make_unique<amrex::MultiFab>(ba, dm, ncomp, ng);

                amrex::Box domain(linop->Geom(lev).Domain());
                domain.convert(amrex::IntVect::TheNodeVector());
                const Dim3 dlo = amrex::lbound(domain), dhi = amrex::ubound(domain);

                for (MFIter mfi(*x[lev], false); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& bx = mfi.validbox();
                    amrex::Array4<Set::Scalar> const& xfab = x[lev]->array(mfi);
                    amrex::Array4<Set::Scalar> const& yfab = y[lev]->array(mfi);
                    amrex::LoopOnCpu(bx, ncomp, [=](int i, int j, int k, int n)
                    {
                        bool AMREX_D_DECL(xmin = (i == dlo.x), ymin = (j == dlo.y), zmin = (k == dlo.z)),
                             AMREX_D_DECL(xmax = (i == dhi.x), ymax = (j == dhi.y), zmax = (k == dhi.z));
                        bool on_boundary = AMREX_D_TERM(xmax || xmin, || ymax || ymin, || zmax || zmin);
                        xfab(i, j, k, n) = on_boundary ? 0.0 : Util::Random() - 0.5;
                        yfab(i, j, k, n) = on_boundary ? 0.0 : Util::Random() - 0.5;
                    });
                }
                x[lev]->FillBoundary(); y[lev]->FillBoundary();
            }

            apply(Ax, x);
            apply(Ay, y);

            Set::Scalar xAy = 0.0, yAx = 0.0, normAx_normy = 0.0;
            for (int lev = 0; lev < nlevs; lev++)
            {
                int ncomp = a_template[lev]->nComp();
                xAy += amrex::MultiFab::Dot(*x[lev], 0, *Ay[lev], 0, ncomp, 0);
                yAx += amrex::MultiFab::Dot(*y[lev], 0, *Ax[lev], 0, ncomp, 0);
                normAx_normy += Ax[lev]->norm0() * y[lev]->norm0();
            }
            Set::Scalar relative_asymmetry =
                (normAx_normy > 0.0) ? std::abs(xAy - yAx) / normAx_normy : std::abs(xAy - yAx);
            max_relative_asymmetry = std::max(max_relative_asymmetry, relative_asymmetry);

            Util::Message(INFO, "checkSymmetry: sample ", s,
                " x.Ay=", xAy, " y.Ax=", yAx,
                " relative asymmetry=", relative_asymmetry);
        }
        return max_relative_asymmetry;
    }

    // A1: does the analytic diagonal (Elastic::Diagonal, a closed form) agree
    // with the probed diagonal (the base Operator<Grid::Node>::Diagonal,
    // which finds the diagonal by hitting Fapply with unit vectors)?
    //
    // Checked only at (amrlev=0, mglev=0) -- the finest level, where
    // divergence would be explained outright with no coarse-grid story at
    // all. Expected to agree to roundoff on the conservative path (the
    // conservative diagonal is built from the same DDW(...,face+1) the
    // conservative Fapply uses for its flux), so this is cheap insurance
    // more than a leading suspect -- but it must pass before any offline
    // matrix assembly (A4) can be trusted.
    Set::Scalar checkDiagonalConsistency(const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_template)
    {
        if (!m_defined) Util::Abort(INFO, "Solver not defined");
        const amrex::BoxArray& ba = a_template[0]->boxArray();
        const amrex::DistributionMapping& dm = a_template[0]->DistributionMap();
        int ncomp = a_template[0]->nComp();

        amrex::MultiFab analytic(ba, dm, ncomp, 0);
        amrex::MultiFab probed(ba, dm, ncomp, 0);
        linop->Diagonal(0, 0, analytic);
        linop->Operator::Operator<Grid::Node>::Diagonal(0, 0, probed);

        amrex::MultiFab diff(ba, dm, ncomp, 0);
        amrex::MultiFab::Copy(diff, analytic, 0, 0, ncomp, 0);
        amrex::MultiFab::Subtract(diff, probed, 0, 0, ncomp, 0);
        Set::Scalar max_abs_diff = diff.norm0();
        Set::Scalar scale = std::max(analytic.norm0(), probed.norm0());
        Set::Scalar relative_diff = (scale > 0.0) ? max_abs_diff / scale : max_abs_diff;

        Util::Message(INFO, "checkDiagonalConsistency: max abs diff=", max_abs_diff,
            " relative diff=", relative_diff);
        return relative_diff;
    }

    // Runs whichever of A5/A1 the user opted into (check_symmetry,
    // check_diagonal), using the just-solved (or just-attempted) solution
    // vector purely as a layout template -- its values are not read.
    void runDiagnosticChecks(const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_template)
    {
        if (m_check_symmetry_trials > 0) checkSymmetry(a_template, m_check_symmetry_trials);
        if (m_check_diagonal) checkDiagonalConsistency(a_template);
    }

    //using MLMG::solve;
protected:
    int max_iter = -1;
    int bottom_max_iter = -1;
    int max_fmg_iter = -1;
    int fixed_iter = -1;
    int verbose = -1;
    int pre_smooth = -1;
    int post_smooth = -1;
    int final_smooth = -1;
    int bottom_smooth = -1;
    std::string bottom_solver;
    Set::Scalar cg_tol_rel = -1.0;
    Set::Scalar cg_tol_abs = -1.0;
    Set::Scalar bottom_tol_rel = -1.0;
    Set::Scalar bottom_tol_abs = -1.0;
    Set::Scalar tol_rel = -1.0;
    Set::Scalar tol_abs = -1.0;
    Set::Scalar omega = -1.0;
    bool average_down_coeffs = false;
    bool normalize_ddw = false;

    Operator::Operator<Grid::Node>* linop;
    amrex::MLMG* mlmg;

    SolveStatus m_last_status = SolveStatus::Unknown;
    bool m_dump_always = false;
    int m_check_symmetry_trials = 0;
    bool m_check_diagonal = false;

    // Classifies the outcome of the mlmg->solve() call just made and, if
    // verbose, prints a single machine-readable line -- SOLVER_STATS -- so
    // scripts (e.g. scripts/solverlib.py) can read iteration counts and
    // residuals without regexing AMReX's prose. Must be called before any
    // Util::Abort so the numbers are still recoverable from stdout.
    void reportSolverStats(bool threw, const char* what)
    {
        if (threw)
        {
            // amrex::MLMG only throws these specific, distinguishable
            // messages ("MLMG blew up." / "MLMG failed to converge.") when
            // its own throw_exception flag is set -- PrepareMLMG does that
            // below whenever dump_on_fail/abort_on_fail=false is in play.
            // Without it, MLMG instead calls the generic amrex::Abort("MLMG
            // failed."), which -- once amrex.signal_handling=0 and
            // amrex.throw_exception=1 are set, as dump_on_fail/
            // abort_on_fail already require -- surfaces here as a caught
            // exception with that generic text, indistinguishable between
            // a stall and a blowup. Default that ambiguous case to Diverged
            // (the worse-case reading) rather than silently mislabeling it
            // Stalled.
            std::string msg = what ? what : "";
            if (msg.find("blew up") != std::string::npos)
                m_last_status = SolveStatus::Diverged;
            else if (msg.find("failed to converge") != std::string::npos)
                m_last_status = SolveStatus::Stalled;
            else
                m_last_status = SolveStatus::Diverged;
        }
        else
        {
            m_last_status = SolveStatus::Converged;
        }

        if (verbose >= 1)
        {
            const char* status_str =
                m_last_status == SolveStatus::Converged ? "converged" :
                m_last_status == SolveStatus::Stalled ? "stalled" :
                m_last_status == SolveStatus::Diverged ? "diverged" : "unknown";
            Util::Message(INFO, "SOLVER_STATS iters=", getNumIters(),
                " resid0=", getInitResidual(),
                " resid=", getFinalResidual(),
                " rhs0=", getInitRHS(),
                " convfactor=", getConvergenceFactor(),
                " status=", status_str);
        }
    }

    void PrepareMLMG(amrex::MLMG& mlmg)
    {
        if (!m_defined) Util::Message(INFO, "Solver not defined");
        mlmg.setBottomSolver(MLMG::BottomSolver::bicgstab);
        mlmg.setCFStrategy(MLMG::CFStrategy::ghostnodes);
        mlmg.setFinalFillBC(false);
        mlmg.setMaxFmgIter(100000000);

        // Whenever the user has opted into non-aborting failure handling
        // (dump_on_fail or abort_on_fail=false -- which already forces
        // amrex.signal_handling=0/amrex.throw_exception=1 in Parse below),
        // also ask MLMG itself to throw. Without this, MLMG's own
        // convergence-failure paths call the generic amrex::Abort("MLMG
        // failed.") instead of throwing their specific "MLMG blew up." /
        // "MLMG failed to converge." messages -- the global amrex settings
        // still convert that Abort into a caught exception, but with
        // indistinguishable text, which defeats reportSolverStats's
        // stalled-vs-diverged classification.
        if (m_dump_on_fail || !m_abort_on_fail) mlmg.setThrowException(true);

        if (max_iter >= 0)        mlmg.setMaxIter(max_iter);
        if (bottom_max_iter >= 0) mlmg.setBottomMaxIter(bottom_max_iter);
        if (max_fmg_iter >= 0)    mlmg.setMaxFmgIter(max_fmg_iter);
        if (fixed_iter >= 0)      mlmg.setFixedIter(fixed_iter);
        if (verbose >= 0)
        {
            mlmg.setVerbose(verbose - 1);
            if (verbose > 4)      mlmg.setBottomVerbose(verbose);
            else                  mlmg.setBottomVerbose(0);
        }

        if (pre_smooth >= 0)      mlmg.setPreSmooth(pre_smooth);
        if (post_smooth >= 0)     mlmg.setPostSmooth(post_smooth);
        if (final_smooth >= 0)    mlmg.setFinalSmooth(final_smooth);
        if (bottom_smooth >= 0)   mlmg.setBottomSmooth(bottom_smooth);

        if (bottom_solver == "cg")       mlmg.setBottomSolver(MLMG::BottomSolver::cg);
        else if (bottom_solver == "bicgstab") mlmg.setBottomSolver(MLMG::BottomSolver::bicgstab);
        else if (bottom_solver == "smoother") mlmg.setBottomSolver(MLMG::BottomSolver::smoother);
        else if (!bottom_solver.empty())
            // amrex::MLLinOp::BottomSolver also has Default, bicgcg, cgbicg,
            // hypre, petsc -- none of those are wired up here, and an
            // unrecognized string was previously ignored silently, leaving
            // whatever the default (bicgstab) already was.
            Util::Warning(INFO, "Unrecognized elastic.solver.bottom_solver = '",
                bottom_solver, "' -- ignored. Recognized: cg, bicgstab, smoother.");

        if (bottom_tol_rel >= 0) mlmg.setBottomTolerance(bottom_tol_rel);
        if (bottom_tol_abs >= 0) mlmg.setBottomToleranceAbs(bottom_tol_abs);

        if (omega >= 0) this->linop->SetOmega(omega);
        if (average_down_coeffs) this->linop->SetAverageDownCoeffs(true);
        if (normalize_ddw) this->linop->SetNormalizeDDW(true);
    }


public:
    // These are the parameters that are read in for a standard 
    // multigrid linear solve.
    static void Parse(Linear& value, amrex::ParmParse& pp)
    {
        // Max number of iterations to perform before erroring out
        pp_query("max_iter", value.max_iter);

        // Max number of iterations on the bottom solver
        pp_query("bottom_max_iter", value.bottom_max_iter);

        // Max number of F-cycle iterations to perform
        pp_query("max_fmg_iter", value.max_fmg_iter);

        // DEPRICATED - do not use
        if (pp.contains("max_fixed_iter"))
            Util::Abort(INFO, "max_fixed_iter is depricated. Use fixed_iter instead.");

        // Number of fixed iterations to perform before exiting gracefully
        pp_query("fixed_iter", value.fixed_iter);

        // Verbosity of the solver (1-5)
        pp_query("verbose", value.verbose);

        // Number of smoothing operations before bottom solve (2)
        pp_query("pre_smooth", value.pre_smooth);

        // Number of smoothing operations after bottom solve (2)
        pp_query("post_smooth", value.post_smooth);

        // Number of final smoothing operations when smoother is used as bottom solver (8)
        pp_query("final_smooth", value.final_smooth);

        // Additional smoothing after bottom CG solver (0)
        pp_query("bottom_smooth", value.bottom_smooth);

        // The method that is used for the multigrid bottom solve (cg, bicgstab, smoother)
        pp_query("bottom_solver", value.bottom_solver);

        if (pp.contains("cg_tol_rel"))
            Util::Abort(INFO, "cg_tol_rel is depricated. Use bottom_tol_rel instead.");
        if (pp.contains("cg_tol_abs"))
            Util::Abort(INFO, "cg_tol_abs is depricated. Use bottom_tol_abs instead.");

        // Relative tolerance on bottom solver
        pp_query("bottom_tol_rel", value.bottom_tol_rel);

        // Absolute tolerance on bottom solver
        pp_query("bottom_tol_abs", value.bottom_tol_abs);

        // Relative tolerance
        pp_query("tol_rel", value.tol_rel);

        // Absolute tolerance
        pp_query("tol_abs", value.tol_abs);

        // Omega (used in gauss-seidel solver)
        pp_query("omega", value.omega);

        // Whether to average down coefficients or use the ones given.
        // (Setting this to true is important for fracture.)
        pp_query("average_down_coeffs", value.average_down_coeffs);

        // Whether to normalize DDW when calculating the diagonal.
        // This is primarily used when DDW is near-singular - like when there
        // is a "void" region or when doing phase field fracture.
        pp_query("normalize_ddw", value.normalize_ddw);

        // [false] 
        // If set to true, output diagnostic multifab information 
        // whenever the MLMG solver fails to converge.
        // (Note: you must also set :code:`amrex.signalhandling=0`
        // and :code:`amrex.throw_exception=1` for this to work.)
        pp_query("dump_on_fail", value.m_dump_on_fail);

        // [false]
        // A3 diagnostic. If set to true, output the same diagnostic
        // multifab information as dump_on_fail after *every* solve, whether
        // it converged or not. Combine with fixed_iter (or a small max_iter
        // plus abort_on_fail=0) to localize where the residual lives after
        // a fixed number of iterations on a case that is merely slow, not
        // one that outright fails.
        pp_query("dump_always", value.m_dump_always);

        // [0]
        // A5 diagnostic. If > 0, run this many random-vector symmetry
        // trials (x.(Ay) vs y.(Ax), boundary DOFs masked out) through
        // Linear::apply after every solve and report the relative
        // asymmetry. See Linear::checkSymmetry for what this can and
        // cannot see.
        pp_query("check_symmetry", value.m_check_symmetry_trials);

        // [false]
        // A1 diagnostic. If true, after every solve compare the analytic
        // diagonal (e.g. Elastic::Diagonal's closed form) against the
        // probed diagonal (the base Operator<Grid::Node>::Diagonal, which
        // hits Fapply with unit vectors) at the finest level, and report
        // the relative difference. See Linear::checkDiagonalConsistency.
        pp_query("check_diagonal", value.m_check_diagonal);

        // [true]
        // If set to false, MLMG will not die if convergence criterion
        // is not reached.
        // (Note: you must also set :code:`amrex.signalhandling=0`
        // and :code:`amrex.throw_exception=1` for this to work.)
        pp_query("abort_on_fail", value.m_abort_on_fail);

        if (value.m_dump_on_fail || !value.m_abort_on_fail)
        {
            IO::ParmParse pp("amrex");
            pp.add("signal_handling", 0);
            pp.add("throw_exception", 1);
        }

    }
protected:
    bool m_defined = false;

    bool m_dump_on_fail = false;
    bool m_abort_on_fail = true;

};
}
}
#endif

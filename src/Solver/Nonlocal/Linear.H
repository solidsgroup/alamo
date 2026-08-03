#ifndef SOLVER_NONLOCAL_LINEAR
#define SOLVER_NONLOCAL_LINEAR
#include "Operator/Operator.H"
#include <AMReX_MLMG.H>

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
        }
        catch (const std::exception& e)
        {
            if (m_dump_on_fail) dumpOnConvergenceFail(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(rhs_tmp));
            if (m_abort_on_fail) Util::Abort(INFO, e.what());
        }
        if (a_sol[0]->contains_nan()) 
        {
            dumpOnConvergenceFail(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(rhs_tmp));
            Util::Abort(INFO);
        }

        return retval;
    };

    Set::Scalar solve(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_sol,
        amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_rhs,
        Real a_tol_rel, Real a_tol_abs, const char* checkpoint_file = nullptr)
    {
        if (m_use_fgmres) return solveFGMRES(a_sol, a_rhs, a_tol_rel, a_tol_abs);
        PrepareMLMG(*mlmg);
        Set::Scalar retval = NAN;
        try
        {
            retval = mlmg->solve(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(a_rhs), a_tol_rel, a_tol_abs, checkpoint_file);
        }
        catch (const std::exception& e)
        {
            if (m_dump_on_fail) dumpOnConvergenceFail(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(a_rhs));
            if (m_abort_on_fail) Util::Abort(INFO, e.what());
        }
        return retval;
    };
    Set::Scalar solve(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_sol,
        amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_rhs)
    {
        if (m_use_fgmres) return solveFGMRES(a_sol, a_rhs, tol_rel, tol_abs);
        PrepareMLMG(*mlmg);
        Set::Scalar retval = NAN;
        try
        {
            retval = mlmg->solve(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(a_rhs), tol_rel, tol_abs);
        }
        catch (const std::exception& e)
        {
            if (m_dump_on_fail) dumpOnConvergenceFail(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(a_rhs));
            if (m_abort_on_fail) Util::Abort(INFO, e.what());
        }
        if (a_sol[0]->contains_nan())
        {
            dumpOnConvergenceFail(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(a_rhs));
            Util::Abort(INFO);
        }
        return retval;
    };
    void apply(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_rhs,
        amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_sol)
    {
        PrepareMLMG(*mlmg);
        mlmg->apply(GetVecOfPtrs(a_rhs), GetVecOfPtrs(a_sol));
    };

    /// \brief Restarted, right-preconditioned FGMRES with MLMG as the preconditioner.
    ///
    /// Outer Krylov acceleration so the solve can survive a poorly-scaled/converging
    /// multigrid operator (near-singular coarse levels, high modulus contrast) instead
    /// of stalling. mlmg->apply() is used as the (exact, linear) matrix-vector product,
    /// and mlmg->precond() -- a fixed number of V/F-cycles from a zero initial guess --
    /// is used as the preconditioner. See getFinalResidual()/getNumIters() after calling
    /// this for diagnostics.
    ///
    /// CAVEAT: mlmg->apply() only computes a correct composite matrix-vector product for
    /// single-level problems. Operator<Grid::Node>::solutionResidual() ignores its
    /// crse_bcdata argument, so coarse/fine ghost coupling is not re-established for an
    /// arbitrary (non-solution) vector the way it is inside a real MLMG V-cycle. On
    /// amr.max_level>0 problems this makes FGMRES diverge; only use solver.fgmres=1 (or
    /// elastic.solver.fgmres=1) with amr.max_level=0 until that gap is fixed.
    Set::Scalar solveFGMRES(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_sol,
        amrex::Vector<std::unique_ptr<amrex::MultiFab> >& a_rhs,
        Real a_tol_rel, Real a_tol_abs)
    {
        if (!m_defined) Util::Abort(INFO, "Solver not defined");
        PrepareMLMG(*mlmg);

        using AMRVec = amrex::Vector<std::unique_ptr<amrex::MultiFab> >;
        const int nlevs = a_sol.size();
        if (nlevs > 1)
            Util::Abort(INFO,
                "solver.fgmres=1 (FGMRES) is only valid for single-level (amr.max_level=0) "
                "problems: Operator<Grid::Node>::solutionResidual() ignores its crse_bcdata "
                "argument, so mlmg->apply() does not re-establish coarse/fine ghost coupling "
                "for an arbitrary vector and FGMRES will diverge on amr.max_level>0 problems. "
                "Set amr.max_level=0 to use FGMRES, or unset solver.fgmres to use plain MLMG.");
        const int ncomp = a_sol[0]->nComp();
        const int restart = fgmres_restart > 0 ? fgmres_restart : 30;
        const int max_outer = fgmres_max_outer > 0 ? fgmres_max_outer : 10;
        const int precond_iters = fgmres_precond_iter > 0 ? fgmres_precond_iter : 1;

        auto makeLike = [&](AMRVec& ref) {
            AMRVec v(nlevs);
            for (int lev = 0; lev < nlevs; ++lev)
                v[lev] = std::make_unique<amrex::MultiFab>(ref[lev]->boxArray(), ref[lev]->DistributionMap(), ncomp, ref[lev]->nGrow());
            return v;
        };
        auto setValAll = [&](AMRVec& v, Set::Scalar val) { for (auto& mf : v) mf->setVal(val); };
        auto dot = [&](AMRVec& x, AMRVec& y) {
            Set::Scalar s = 0.0;
            for (int lev = 0; lev < nlevs; ++lev)
                s += amrex::MultiFab::Dot(*x[lev], 0, *y[lev], 0, ncomp, 0);
            return s;
        };
        auto vnorm = [&](AMRVec& x) { return std::sqrt(dot(x, x)); };
        auto axpy = [&](AMRVec& y, Set::Scalar a, AMRVec& x) {
            for (int lev = 0; lev < nlevs; ++lev)
                amrex::MultiFab::Saxpy(*y[lev], a, *x[lev], 0, 0, ncomp, 0);
        };
        auto copyv = [&](AMRVec& dst, AMRVec& src) {
            for (int lev = 0; lev < nlevs; ++lev)
                amrex::MultiFab::Copy(*dst[lev], *src[lev], 0, 0, ncomp, src[lev]->nGrow());
        };
        auto scalev = [&](AMRVec& x, Set::Scalar a) {
            for (int lev = 0; lev < nlevs; ++lev) x[lev]->mult(a, 0, ncomp, 0);
        };
        auto matvec = [&](AMRVec& v, AMRVec& Av) {
            mlmg->apply(GetVecOfPtrs(Av), GetVecOfPtrs(v));
        };
        auto precond = [&](AMRVec& r, AMRVec& z) {
            setValAll(z, 0.0);
            mlmg->setPrecondIter(precond_iters);
            try { mlmg->precond(GetVecOfPtrs(z), GetVecOfConstPtrs(r), 1.0e-16, 1.0e-16); }
            catch (const std::exception&) {}
        };

        AMRVec r = makeLike(a_sol);
        AMRVec Ax = makeLike(a_sol);
        matvec(a_sol, Ax);
        copyv(r, a_rhs);
        axpy(r, -1.0, Ax);

        Set::Scalar bnorm = vnorm(a_rhs);
        if (bnorm < 1.0e-300) bnorm = 1.0;
        Set::Scalar rnorm = vnorm(r);
        Set::Scalar tol = std::max(a_tol_rel * bnorm, a_tol_abs);

        if (verbose > 0) Util::Message(INFO, "FGMRES: initial residual = ", rnorm, " (target = ", tol, ")");

        for (int outer = 0; outer < max_outer && rnorm > tol; ++outer)
        {
            std::vector<AMRVec> V(restart + 1);
            std::vector<AMRVec> Z(restart);
            std::vector<std::vector<Set::Scalar> > H(restart + 1, std::vector<Set::Scalar>(restart, 0.0));
            std::vector<Set::Scalar> cs(restart, 0.0), sn(restart, 0.0), g(restart + 1, 0.0);

            Set::Scalar beta = rnorm;
            V[0] = makeLike(a_sol);
            copyv(V[0], r);
            scalev(V[0], 1.0 / beta);
            g[0] = beta;

            int kk = -1;
            for (int j = 0; j < restart; ++j)
            {
                Z[j] = makeLike(a_sol);
                precond(V[j], Z[j]);

                V[j + 1] = makeLike(a_sol);
                matvec(Z[j], V[j + 1]);

                for (int i = 0; i <= j; ++i)
                {
                    H[i][j] = dot(V[j + 1], V[i]);
                    axpy(V[j + 1], -H[i][j], V[i]);
                }
                H[j + 1][j] = vnorm(V[j + 1]);
                if (H[j + 1][j] > 1.0e-300) scalev(V[j + 1], 1.0 / H[j + 1][j]);

                for (int i = 0; i < j; ++i)
                {
                    Set::Scalar temp = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
                    H[i + 1][j] = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
                    H[i][j] = temp;
                }
                Set::Scalar denom = std::sqrt(H[j][j] * H[j][j] + H[j + 1][j] * H[j + 1][j]);
                if (denom < 1.0e-300) denom = 1.0e-300;
                cs[j] = H[j][j] / denom;
                sn[j] = H[j + 1][j] / denom;
                H[j][j] = cs[j] * H[j][j] + sn[j] * H[j + 1][j];
                H[j + 1][j] = 0.0;
                g[j + 1] = -sn[j] * g[j];
                g[j] = cs[j] * g[j];

                kk = j;
                rnorm = std::abs(g[j + 1]);
                if (verbose > 0) Util::Message(INFO, "FGMRES: outer ", outer, " inner ", j, " residual = ", rnorm);
                if (rnorm <= tol) break;
            }

            std::vector<Set::Scalar> y(kk + 1);
            for (int i = kk; i >= 0; --i)
            {
                y[i] = g[i];
                for (int k = i + 1; k <= kk; ++k) y[i] -= H[i][k] * y[k];
                y[i] /= H[i][i];
            }
            for (int i = 0; i <= kk; ++i)
                axpy(a_sol, y[i], Z[i]);

            matvec(a_sol, Ax);
            copyv(r, a_rhs);
            axpy(r, -1.0, Ax);
            rnorm = vnorm(r);
            if (verbose > 0) Util::Message(INFO, "FGMRES: outer ", outer, " true residual = ", rnorm);
        }

        if (a_sol[0]->contains_nan())
        {
            dumpOnConvergenceFail(GetVecOfPtrs(a_sol), GetVecOfConstPtrs(a_rhs));
            Util::Abort(INFO);
        }

        return rnorm;
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

    bool m_use_fgmres = false;
    int fgmres_restart = -1;
    int fgmres_max_outer = -1;
    int fgmres_precond_iter = -1;

    Operator::Operator<Grid::Node>* linop;
    amrex::MLMG* mlmg;

    void PrepareMLMG(amrex::MLMG& mlmg)
    {
        if (!m_defined) Util::Message(INFO, "Solver not defined");
        mlmg.setBottomSolver(MLMG::BottomSolver::bicgstab);
        mlmg.setCFStrategy(MLMG::CFStrategy::ghostnodes);
        mlmg.setFinalFillBC(false);
        mlmg.setMaxFmgIter(100000000);


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
        // Use a hand-rolled restarted FGMRES (MLMG as preconditioner) instead of
        // calling MLMG::solve() directly. Intended for cases where MLMG alone
        // stalls/diverges (e.g. high modulus contrast); more robust but slower
        // per-iteration than a plain V/F-cycle. CAVEAT: only valid for single-level
        // (amr.max_level=0) problems -- see solveFGMRES() doc comment for why.
        pp_query("fgmres", value.m_use_fgmres);

        // Restart length (number of Krylov vectors kept before restarting) [30]
        pp_query("fgmres_restart", value.fgmres_restart);

        // Max number of outer restart cycles [10]
        pp_query("fgmres_max_outer", value.fgmres_max_outer);

        // Number of MLMG V/F-cycles used per FGMRES preconditioner application [1]
        pp_query("fgmres_precond_iter", value.fgmres_precond_iter);

        // [false] 
        // If set to true, output diagnostic multifab information 
        // whenever the MLMG solver fails to converge.
        // (Note: you must also set :code:`amrex.signalhandling=0`
        // and :code:`amrex.throw_exception=1` for this to work.)
        pp_query("dump_on_fail", value.m_dump_on_fail);

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

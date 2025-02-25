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
            linop->realFillBoundary(*rhs_tmp[lev], linop->Geom(lev));
            Util::Message(INFO, rhs_tmp[lev]->norm0());
            //rhs_tmp[lev]->FillBoundary();
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

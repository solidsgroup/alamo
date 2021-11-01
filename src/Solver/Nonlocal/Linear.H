#ifndef SOLVER_NONLOCAL_LINEAR
#define SOLVER_NONLOCAL_LINEAR
#include "Operator/Operator.H"
#include <AMReX_MLMG.H>
#include "IC/Trig.H"

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
class Linear : public amrex::MLMG
{
public:
    /// This constructor sets several defaults for NODAL solves
    Linear (Operator::Operator<Grid::Node>& a_lp) : MLMG(a_lp), linop(a_lp)
    {
        MLMG::setBottomSolver(MLMG::BottomSolver::bicgstab);
        MLMG::setCFStrategy(MLMG::CFStrategy::ghostnodes);
        MLMG::setFinalFillBC(false);
        MLMG::setMaxFmgIter(100000000);
    }
    void setVerbose(int verbosity)
    {
        m_verbose = verbosity;
        MLMG::setVerbose(m_verbose-1);
        if (m_verbose > 4) MLMG::setBottomVerbose(m_verbose);
        else
        {
            MLMG::setBottomVerbose(0);
        }
    }

    Set::Scalar solveaffine (amrex::Vector<std::unique_ptr<amrex::MultiFab> > & a_sol, 
                        amrex::Vector<std::unique_ptr<amrex::MultiFab> > & a_rhs,
                        Real a_tol_rel, Real a_tol_abs, bool copyrhs = false, 
                        const char* checkpoint_file = nullptr)
    {
        amrex::Vector<amrex::MultiFab *> rhs_tmp(a_rhs.size());
        amrex::Vector<amrex::MultiFab *> zero_tmp(a_rhs.size());
        for (int i = 0; i < rhs_tmp.size(); i++)
        {
            rhs_tmp[i]  = new amrex::MultiFab(a_rhs[i]->boxArray(),a_rhs[i]->DistributionMap(),a_rhs[i]->nComp(),a_rhs[i]->nGrow());
            zero_tmp[i] = new amrex::MultiFab(a_rhs[i]->boxArray(),a_rhs[i]->DistributionMap(),a_rhs[i]->nComp(),a_rhs[i]->nGrow());
            rhs_tmp[i]->setVal(0.0);
            zero_tmp[i]->setVal(0.0);
            Util::Message(INFO,rhs_tmp[i]->norm0());
        }

        linop.SetHomogeneous(false);
        MLMG::apply(rhs_tmp,zero_tmp);

        for (int lev = 0; lev < rhs_tmp.size(); lev++)
        {
            amrex::Box domain = linop.Geom(lev).Domain();
            domain.convert(amrex::IntVect::TheNodeVector());
            const Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
            for (MFIter mfi(*rhs_tmp[lev],amrex::TilingIfNotGPU());mfi.isValid();++mfi)
            {
                amrex::Box bx = mfi.growntilebox(rhs_tmp[lev]->nGrow());
                bx = bx & domain;
                amrex::Array4<amrex::Real> const& rhstmp = rhs_tmp[lev]->array(mfi);
                for (int n = 0; n < rhs_tmp[lev]->nComp(); n++)
                {
                    amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) 
                    {
                        bool    AMREX_D_DECL(xmin = (i == lo.x), ymin = (j==lo.y), zmin = (k==lo.z)),
                                AMREX_D_DECL(xmax = (i == hi.x), ymax = (j==hi.y), zmax = (k==hi.z));
                        if (AMREX_D_TERM(xmax || xmin, || ymax || ymin, || zmax || zmin))
                            rhstmp(i,j,k,n) = 0.0;
                        else
                            rhstmp(i,j,k,n) *= -1.0;
                    });
                }
            }
            Util::Message(INFO,rhs_tmp[lev]->norm0());
            linop.realFillBoundary(*rhs_tmp[lev],linop.Geom(lev));
            Util::Message(INFO,rhs_tmp[lev]->norm0());
            //rhs_tmp[lev]->FillBoundary();
        }

        for (int lev = 0; lev < rhs_tmp.size(); lev++)
        {
            Util::Message(INFO,rhs_tmp[lev]->norm0());
            amrex::Add(*rhs_tmp[lev],*a_rhs[lev],0,0,rhs_tmp[lev]->nComp(),rhs_tmp[lev]->nGrow());
            if (copyrhs)
                amrex::Copy(*a_rhs[lev],*rhs_tmp[lev],0,0,rhs_tmp[lev]->nComp(),rhs_tmp[lev]->nGrow());
            Util::Message(INFO,rhs_tmp[lev]->norm0());
        }

        linop.SetHomogeneous(true);
        return MLMG::solve(GetVecOfPtrs(a_sol),GetVecOfConstPtrs(rhs_tmp),a_tol_rel,a_tol_abs,checkpoint_file);
    };

    Set::Scalar solve (amrex::Vector<std::unique_ptr<amrex::MultiFab> > & a_sol, 
                        amrex::Vector<std::unique_ptr<amrex::MultiFab> > & a_rhs,
                        Real a_tol_rel, Real a_tol_abs, const char* checkpoint_file = nullptr)
    {
        return MLMG::solve(GetVecOfPtrs(a_sol),GetVecOfConstPtrs(a_rhs),a_tol_rel,a_tol_abs,checkpoint_file);
    };
    Set::Scalar solve (amrex::Vector<std::unique_ptr<amrex::MultiFab> > & a_sol, 
                        amrex::Vector<std::unique_ptr<amrex::MultiFab> > & a_rhs)
    {
        return MLMG::solve(GetVecOfPtrs(a_sol),GetVecOfConstPtrs(a_rhs),m_tol_rel,m_tol_abs);
    };

    using MLMG::solve;
protected:
    Operator::Operator<Grid::Node> &linop;
    int m_verbose = 0;
    Set::Scalar m_tol_rel = 1E-8, m_tol_abs = 1E-8;

public:
    // These are the parameters that are read in for a standard 
    // multigrid linear solve.
    static void Parse(Linear & value, amrex::ParmParse & pp)
    {
        if (pp.contains("max_iter"))
        // Max number of iterations to perform before erroring out
        { int max_iter; pp.query("max_iter",max_iter);value.setMaxIter(max_iter);}

        if (pp.contains("bottom_max_iter"))
        // Max number of iterations on the bottom solver
        { int bottom_max_iter; pp.query("bottom_max_iter",bottom_max_iter);value.setBottomMaxIter(bottom_max_iter);}

        if (pp.contains("max_fmg_iter"))
        // Max number of F-cycle iterations to perform
        { int max_fmg_iter; pp.query("max_fmg_iter",max_fmg_iter);value.setMaxFmgIter(max_fmg_iter);}

        // DEPRICATED - do not use
        if (pp.contains("max_fixed_iter"))
            Util::Abort(INFO, "max_fixed_iter is depricated. Use fixed_iter instead." );

        if (pp.contains("fixed_iter"))
        // Number of fixed iterations to perform before exiting gracefully
        { int fixed_iter; pp.query("fixed_iter",fixed_iter);value.setFixedIter(fixed_iter);}

        if (pp.contains("verbose"))
        // Verbosity of the solver (1-5)
        { int verbose; pp.query("verbose",verbose);value.setVerbose(verbose);}

        if (pp.contains("pre_smooth"))
        // Number of smoothing operations before bottom solve (2)
        { int pre_smooth; pp.query("pre_smooth",pre_smooth);value.setPreSmooth(pre_smooth);}

        if (pp.contains("post_smooth"))
        // Number of smoothing operations after bottom solve (2)
        { int post_smooth; pp.query("post_smooth",post_smooth);value.setPostSmooth(post_smooth);}
    
        if (pp.contains("bottom_solver"))
        {
            std::string bottom_solver;
            // Type of bottom solver to use [cg, bicgstab, smoother]
            pp.query("bottom_solver",bottom_solver);
            if (bottom_solver == "cg")
                value.setBottomSolver(MLMG::BottomSolver::cg);
            if (bottom_solver == "bicgstab")
                value.setBottomSolver(MLMG::BottomSolver::bicgstab);
            if (bottom_solver == "smoother")
                value.setBottomSolver(MLMG::BottomSolver::smoother);
            else
                Util::Abort(INFO,"Unsupported bottom_solver " + bottom_solver);
        }
        
        if (pp.contains("cg_tol_rel"))
            Util::Abort(INFO,"cg_tol_rel is depricated. Use bottom_tol_rel instead.");
        if (pp.contains("cg_tol_abs"))
            Util::Abort(INFO,"cg_tol_abs is depricated. Use bottom_tol_abs instead.");
        
        if (pp.contains("bottom_tol_rel"))
        {
            Set::Scalar bottom_tol_rel;
            // Relative tolerance on bottom solver
            pp.query("bottom_tol_rel",bottom_tol_rel);
            value.setBottomTolerance(bottom_tol_rel);
        }
        if (pp.contains("bottom_tol_abs"))
        {
            Set::Scalar bottom_tol_abs;
            // Absolute tolerance on bottom solver
            pp.query("bottom_tol_abs",bottom_tol_abs);
            value.setBottomToleranceAbs(bottom_tol_abs);
        }

        // Relative tolerance
        pp.query("tol_rel",value.m_tol_rel);
        // Absolute tolerance
        pp.query("tol_abs",value.m_tol_abs);
    }

    
};
}
}
#endif

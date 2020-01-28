#ifndef SOLVER_NONLOCAL_NEWTON
#define SOLVER_NONLOCAL_NEWTON

#include "Set/Set.H"
#include "Operator/Elastic.H"
#include "Solver/Nonlocal/Linear.H"
#include "IO/ParmParse.H"
#include "Model/Solid/Elastic/NeoHookean.H"
#include "Numeric/Stencil.H"

namespace Solver
{
namespace Nonlocal
{
template <class T>
class Newton : public Linear
{
public:
    Newton(Operator::Elastic<T> &a_op) : Linear(a_op){};

    void setNRIters(int a_nriters) { m_nriters = a_nriters; }

    Set::Scalar solve (const amrex::Vector<MultiFab*>& a_sol_mf, 
                       const amrex::Vector<MultiFab*>& a_rhs_mf,
                       amrex::Vector<amrex::FabArray<amrex::BaseFab<T> > > &a_model,
                       Real a_tol_rel, Real a_tol_abs, bool copyrhs = false, const char* checkpoint_file = nullptr)
    {
        Operator::Elastic<T> &elastic = dynamic_cast<Operator::Elastic<T> &>(linop);

        amrex::Vector<std::unique_ptr<amrex::MultiFab>> dsol_mf, rhs_mf;
        amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<Set::Matrix> > > > dw_mf;
        dsol_mf.resize(a_sol_mf.size());
        dw_mf.resize(a_sol_mf.size());
        rhs_mf.resize(a_sol_mf.size());
        for (int lev = 0; lev < a_sol_mf.size(); lev++)
        {
            dsol_mf[lev].reset(new amrex::MultiFab(a_sol_mf[lev]->boxArray(),
                                                   a_sol_mf[lev]->DistributionMap(),
                                                   a_sol_mf[lev]->nComp(), a_sol_mf[lev]->nGrow()));
            dw_mf[lev].reset(new amrex::FabArray<amrex::BaseFab<Set::Matrix> >(a_rhs_mf[lev]->boxArray(),
                                                 a_rhs_mf[lev]->DistributionMap(),
                                                 1, a_rhs_mf[lev]->nGrow()));
            rhs_mf[lev].reset(new amrex::MultiFab(a_rhs_mf[lev]->boxArray(),
                                                  a_rhs_mf[lev]->DistributionMap(),
                                                  a_rhs_mf[lev]->nComp(), a_rhs_mf[lev]->nGrow()));

            
            dsol_mf[lev]->setVal(0.0);
            dw_mf[lev]->setVal(Set::Matrix::Zero());
            
            amrex::MultiFab::Copy(*rhs_mf[lev], *a_rhs_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
        }

        for (int nriter = 0; nriter < (copyrhs ? m_nriters+1 : m_nriters); nriter++)
        {
            if (m_verbose > 0 && nriter < m_nriters) Util::Message(INFO, "Newton Iteration ", nriter+1, " of ", m_nriters);

            //

            BL_PROFILE_VAR_START(newton_initialization);
            for (int lev = 0; lev < a_rhs_mf.size(); ++lev)
            {
                amrex::Box domain(linop.Geom(lev).Domain());
                domain.convert(amrex::IntVect::TheNodeVector());
                const Set::Scalar *dx = linop.Geom(lev).CellSize();
                Set::Vector DX(linop.Geom(lev).CellSize());

                //
                //  TODO: implement a general method for handling homogeneous / inhomogeneous BCs
                //        Right now, only Dirichlet works.
                //
                //if (nriter > 0) elastic.GetBC().Init((rhs_mf[lev]).get(),elastic.Geom(lev),true);
                
                for (MFIter mfi(a_model[lev], amrex::TilingIfNotGPU); mfi.isValid(); ++mfi)
                {
                    amrex::Box bx = mfi.growntilebox(2);
                    amrex::Box domainint = domain;
                    domainint.grow(-1); // shrink so that we are only on the interior of the domain
                    bx = bx & domain;

                    amrex::Array4<Set::Scalar> const &u = a_sol_mf[lev]->array(mfi);
                    amrex::Array4<Set::Matrix> const &dw = dw_mf[lev]->array(mfi);
                    amrex::Array4<Set::Scalar> const &graddw = rhs_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const &b = a_rhs_mf[lev]->array(mfi);
                    amrex::Array4<T> const &model = a_model[lev].array(mfi);

                    // Set model internal dw and ddw.
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        std::array<Numeric::StencilType, AMREX_SPACEDIM> sten = Numeric::GetStencil(i, j, k, bx);

                        Set::Matrix gradu = Numeric::Gradient(u, i, j, k, dx, sten);

                        if (model(i,j,k).kinvar == Model::Solid::KinematicVariable::gradu)
                        {
                            dw(i,j,k) = model(i, j, k).DW(gradu);
                            model(i, j, k).ddw = model(i, j, k).DDW(gradu);
                        }
                        else if (model(i,j,k).kinvar == Model::Solid::KinematicVariable::epsilon)
                        {
                            Set::Matrix eps = 0.5 * (gradu + gradu.transpose());
                            dw(i,j,k) = model(i, j, k).DW(eps);
                            model(i, j, k).ddw = model(i, j, k).DDW(eps);
                        }
                        else if (model(i,j,k).kinvar == Model::Solid::KinematicVariable::F)
                        {
                            Set::Matrix F = gradu + Set::Matrix::Identity();
                            dw(i,j,k) = model(i, j, k).DW(F);
                            model(i, j, k).ddw = model(i, j, k).DDW(F);
                        }
                    });

                    if (nriter > 0)
                    rhs_mf[lev]->setVal(0.0);

                    bx = mfi.validbox() & domainint;
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        for (int p = 0; p < AMREX_SPACEDIM; p++)
                        {
                            #if AMREX_SPACEDIM == 2
                            graddw(i,j,k,p) = b(i,j,k,p) - 
                                ((dw(i + 1, j, k)(p, 0) - dw(i - 1, j, k)(p, 0)) / 2. / DX(0) +
                                 (dw(i, j + 1, k)(p, 1) - dw(i, j - 1, k)(p, 1)) / 2. / DX(1));
                            #elif AMREX_SPACEDIM == 3
                            graddw(i,j,k,p) = b(i,j,k,p) - 
                                ((dw(i + 1, j, k)(p, 0) - dw(i - 1, j, k)(p, 0)) / 2. / DX(0) +
                                 (dw(i, j + 1, k)(p, 1) - dw(i, j - 1, k)(p, 1)) / 2. / DX(1) +
                                 (dw(i, j, k + 1)(p, 2) - dw(i, j, k - 1)(p, 2)) / 2. / DX(2));
                            #endif
                        }
                    });
                }

                amrex::Geometry geom = elastic.Geom(lev);
                for (int i = 0; i < 2; i++)
                {
                    amrex::FabArray<amrex::BaseFab<T>> &mf = a_model[lev];
                    mf.FillBoundary(geom.periodicity());
                    const int ncomp = mf.nComp();
                    const int ng1 = 1;
                    const int ng2 = 2;
                    amrex::FabArray<amrex::BaseFab<T>> tmpmf(mf.boxArray(), mf.DistributionMap(), ncomp, ng1);
                    amrex::Copy(tmpmf, mf, 0, 0, ncomp, ng1);
                    mf.ParallelCopy(tmpmf, 0, 0, ncomp, ng1, ng2, geom.periodicity());
                }
            }
            BL_PROFILE_VAR_STOP(newton_initialization);
            elastic.SetModel(a_model);

            if (nriter == m_nriters) break;
            
            Solver::Nonlocal::Linear::solve(dsol_mf,rhs_mf, a_tol_rel, a_tol_abs,checkpoint_file);

            Set::Scalar cornorm = 0, solnorm = 0;
            for (int lev = 0; lev < dsol_mf.size(); ++lev)
            {
                for (int comp = 0; comp < AMREX_SPACEDIM; comp++)
                {
                    Set::Scalar tmpcornorm = dsol_mf[lev]->norm0(comp,0);
                    if (tmpcornorm > cornorm) cornorm = tmpcornorm;

                    Set::Scalar tmpsolnorm = a_sol_mf[lev]->norm0(comp,0);
                    if (tmpsolnorm > solnorm) solnorm = tmpsolnorm;
                }

            }
            Set::Scalar relnorm;
            if (solnorm == 0) relnorm = cornorm;
            else relnorm = cornorm / solnorm;
            if (m_verbose > 1) Util::Message(INFO,"NR relative norm(ddisp) = ",relnorm);

            for (int lev = 0; lev < dsol_mf.size(); ++lev)
                amrex::MultiFab::Add(*a_sol_mf[lev], *dsol_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
        }

        if (copyrhs)
        {
            for (int lev = 0; lev < a_sol_mf.size(); lev++)
                amrex::MultiFab::Copy(*a_rhs_mf[lev], *rhs_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
        }
        return 0.0;
    }

    Set::Scalar solve(amrex::Vector<std::unique_ptr<amrex::MultiFab> > & a_sol_mf,
                      amrex::Vector<std::unique_ptr<amrex::MultiFab> > & a_rhs_mf,
                      amrex::Vector<amrex::FabArray<amrex::BaseFab<T> > > &a_model,
                      Set::Scalar a_tol_rel, Set::Scalar a_tol_abs,
                      bool copyrhs = false, const char *checkpoint_file = nullptr)
    {
        return solve(GetVecOfPtrs(a_sol_mf),
                     GetVecOfPtrs(a_rhs_mf),a_model,a_tol_rel,a_tol_abs,copyrhs,checkpoint_file);
    }


private:
    int m_nriters = 1;

public:
    static void Parse(Newton<T> & value, amrex::ParmParse & pp)
    {
        Linear::Parse(value,pp);
        
        pp.query("nriters",value.m_nriters);
    }

};
} // namespace Nonlocal
} // namespace Solver


#endif

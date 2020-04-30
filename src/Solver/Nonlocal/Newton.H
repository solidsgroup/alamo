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
    Newton(Operator::Elastic<T> &a_op) : Linear(a_op), m_elastic(dynamic_cast<Operator::Elastic<T> &>(linop))
    {};

    void setNRIters(int a_nriters) { m_nriters = a_nriters; }


private:

    void prepareForSolve(const Set::Field<Set::Scalar>& a_u_mf, 
                         const Set::Field<Set::Scalar>& a_b_mf,
                         const Set::Field<Set::Scalar>& a_rhs_mf,
                         Set::Field<Set::Matrix> &a_dw_mf,
                         Set::Field<T> &a_model_mf)
    {
            for (int lev = 0; lev < a_b_mf.size(); ++lev)
            {
                amrex::Box domain(linop.Geom(lev).Domain());
                domain.convert(amrex::IntVect::TheNodeVector());
                const Set::Scalar *dx = linop.Geom(lev).CellSize();
                Set::Vector DX(linop.Geom(lev).CellSize());
                
                for (MFIter mfi(*a_model_mf[lev], amrex::TilingIfNotGPU); mfi.isValid(); ++mfi)
                {
                    amrex::Box bx = mfi.growntilebox(2);
                    amrex::Box domainint = domain;
                    domainint.grow(-1); // shrink so that we are only on the interior of the domain
                    bx = bx & domain;

                    amrex::Array4<const Set::Scalar> const &u     = a_u_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const &b     = a_b_mf[lev]->array(mfi);
                    amrex::Array4<Set::Scalar>       const &rhs   = a_rhs_mf[lev]->array(mfi);
                    amrex::Array4<Set::Matrix>       const &dw    = a_dw_mf[lev]->array(mfi);
                    amrex::Array4<T>                 const &model = a_model_mf[lev]->array(mfi);

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


                    bx = mfi.validbox() & domainint;
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        for (int p = 0; p < AMREX_SPACEDIM; p++)
                        {
                            #if AMREX_SPACEDIM == 2
                            rhs(i,j,k,p) = b(i,j,k,p) - 
                                ((dw(i + 1, j, k)(p, 0) - dw(i - 1, j, k)(p, 0)) / 2. / DX(0) +
                                 (dw(i, j + 1, k)(p, 1) - dw(i, j - 1, k)(p, 1)) / 2. / DX(1));
                            #elif AMREX_SPACEDIM == 3
                            rhs(i,j,k,p) = b(i,j,k,p) - 
                                ((dw(i + 1, j, k)(p, 0) - dw(i - 1, j, k)(p, 0)) / 2. / DX(0) +
                                 (dw(i, j + 1, k)(p, 1) - dw(i, j - 1, k)(p, 1)) / 2. / DX(1) +
                                 (dw(i, j, k + 1)(p, 2) - dw(i, j, k - 1)(p, 2)) / 2. / DX(2));
                            #endif
                        }
                    });
                }
                Util::RealFillBoundary(*a_model_mf[lev],m_elastic.Geom(lev));
            }
    }

public:
    Set::Scalar solve (const Set::Field<Set::Scalar> & a_u_mf, 
                       const Set::Field<Set::Scalar> & a_b_mf,
                       Set::Field<T> &a_model_mf,
                       Real a_tol_rel, Real a_tol_abs, const char* checkpoint_file = nullptr)
    {
        Set::Field<Set::Scalar> dsol_mf, rhs_mf;
        Set::Field<Set::Matrix> dw_mf;

        dsol_mf.resize(a_u_mf.size());
        dw_mf.resize(a_u_mf.size());
        rhs_mf.resize(a_u_mf.size());
        for (int lev = 0; lev < a_u_mf.size(); lev++)
        {
            dsol_mf.Define(lev, a_u_mf[lev]->boxArray(),
                                a_u_mf[lev]->DistributionMap(),
                                a_u_mf[lev]->nComp(), 
                                a_u_mf[lev]->nGrow());
            dw_mf.Define(lev,   a_b_mf[lev]->boxArray(),
                                a_b_mf[lev]->DistributionMap(),
                                1, 
                                a_b_mf[lev]->nGrow());
            rhs_mf.Define(lev,  a_b_mf[lev]->boxArray(),
                                a_b_mf[lev]->DistributionMap(),
                                a_b_mf[lev]->nComp(), 
                                a_b_mf[lev]->nGrow());
            
            dsol_mf[lev]->setVal(0.0);
            dw_mf[lev]->setVal(Set::Matrix::Zero());
            
            amrex::MultiFab::Copy(*rhs_mf[lev], *a_b_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
        }

        for (int lev = 0; lev < a_b_mf.size(); ++lev)
            m_elastic.GetBC().Init(rhs_mf[lev].get(),m_elastic.Geom(lev),false);

        for (int nriter = 0; nriter < m_nriters; nriter++)
        {
            if (m_verbose > 0 && nriter < m_nriters) Util::Message(INFO, "Newton Iteration ", nriter+1, " of ", m_nriters);

            for (int lev = 0; lev < a_b_mf.size(); ++lev)
                m_elastic.GetBC().Init(rhs_mf[lev].get(),m_elastic.Geom(lev),nriter>0);

            prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf, a_model_mf);

            m_elastic.SetModel(a_model_mf);

            if (nriter == m_nriters) break;
            
            Solver::Nonlocal::Linear::solve(GetVecOfPtrs(dsol_mf), GetVecOfConstPtrs(rhs_mf), a_tol_rel, a_tol_abs,checkpoint_file);

            Set::Scalar cornorm = 0, solnorm = 0;
            for (int lev = 0; lev < dsol_mf.size(); ++lev)
            {
                for (int comp = 0; comp < AMREX_SPACEDIM; comp++)
                {
                    Set::Scalar tmpcornorm = dsol_mf[lev]->norm0(comp,0);
                    if (tmpcornorm > cornorm) cornorm = tmpcornorm;

                    Set::Scalar tmpsolnorm = a_u_mf[lev]->norm0(comp,0);
                    if (tmpsolnorm > solnorm) solnorm = tmpsolnorm;
                }

            }
            Set::Scalar relnorm;
            if (solnorm == 0) relnorm = cornorm;
            else relnorm = cornorm / solnorm;
            if (m_verbose > 1) Util::Message(INFO,"NR relative norm(ddisp) = ",relnorm);

            for (int lev = 0; lev < dsol_mf.size(); ++lev)
                amrex::MultiFab::Add(*a_u_mf[lev], *dsol_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
        }

        return 0.0;
    }

//    Set::Scalar solve(amrex::Vector<std::unique_ptr<amrex::MultiFab> > & a_u_mf,
//                      amrex::Vector<std::unique_ptr<amrex::MultiFab> > & a_b_mf,
//                      amrex::Vector<amrex::FabArray<amrex::BaseFab<T> > > &a_model,
//                      Set::Scalar a_tol_rel, Set::Scalar a_tol_abs,
//                      const char *checkpoint_file = nullptr)
//    {
//        return solve(GetVecOfPtrs(a_u_mf),
//                     GetVecOfPtrs(a_b_mf),a_model,a_tol_rel,a_tol_abs,checkpoint_file);
//    }

    void compResidual(Set::Field<Set::Scalar> & a_res_mf,
                      Set::Field<Set::Scalar> & a_u_mf,
                      Set::Field<Set::Scalar> & a_b_mf,
                      Set::Field<T> &a_model_mf)
    {
        Set::Field<Set::Matrix> dw_mf;
        dw_mf.resize(a_u_mf.size());
        for (int lev = 0; lev < a_u_mf.size(); lev++)
        {
            dw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                              a_b_mf[lev]->DistributionMap(),
                              1, a_b_mf[lev]->nGrow());
            dw_mf[lev]->setVal(Set::Matrix::Zero());
        }
        
        for (int lev = 0; lev < a_b_mf.size(); ++lev)
            m_elastic.GetBC().Init(a_b_mf[lev].get(),m_elastic.Geom(lev),true);
            
        prepareForSolve(a_u_mf, a_b_mf, a_res_mf, dw_mf, a_model_mf);
        
    }

private:
    int m_nriters = 1;
    Operator::Elastic<T> &m_elastic;

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
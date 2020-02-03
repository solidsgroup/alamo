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
    Newton(Operator::Elastic<T> &a_op) : 
        Linear(a_op), 
        m_elastic(dynamic_cast<Operator::Elastic<T> &>(linop)),
        m_bc(m_elastic.GetBC())
    {};

    void setNRIters(int a_nriters) { m_nriters = a_nriters; }


private:

#define SQRT3INV 0.57735026919
#define SQRT2INV 0.70710678118

    void prepareForSolve(const amrex::Vector<MultiFab*>& a_u_mf, 
                         const amrex::Vector<MultiFab*>& a_b_mf,
                         const amrex::Vector<MultiFab*>& a_rhs_mf,
                         amrex::Vector<amrex::FabArray<amrex::BaseFab<Set::Matrix> > > &a_dw_mf,
                         amrex::Vector<amrex::FabArray<amrex::BaseFab<T> > > &a_model_mf)
    {
            for (int lev = 0; lev < a_b_mf.size(); ++lev)
            {
                amrex::Box domain(linop.Geom(lev).Domain());
                domain.convert(amrex::IntVect::TheNodeVector());
                const Set::Scalar *dx = linop.Geom(lev).CellSize();
                Set::Vector DX(linop.Geom(lev).CellSize());

                const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);
                
                for (MFIter mfi(a_model_mf[lev], amrex::TilingIfNotGPU); mfi.isValid(); ++mfi)
                {
                    amrex::Box bx = mfi.growntilebox(2);
                    bx = bx & domain;

                    amrex::Array4<const Set::Scalar> const &u     = a_u_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const &b     = a_b_mf[lev]->array(mfi);
                    amrex::Array4<Set::Scalar>       const &rhs   = a_rhs_mf[lev]->array(mfi);
                    amrex::Array4<Set::Matrix>       const &dw    = a_dw_mf[lev].array(mfi);
                    amrex::Array4<T>                 const &model = a_model_mf[lev].array(mfi);

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


                    amrex::Box domainint = domain; domainint.grow(-1); // shrink so that we are only on the interior of the domain
                    
                    bx               = mfi.validbox() & domain;
                    //amrex::Box bxint = mfi.validbox() & domainint;
                    //Util::Message(INFO,bx);
                    //Util::Message(INFO,bxint);

                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        #if AMREX_SPACEDIM == 2
                        if (i == lo.x || i == hi.x || j == lo.y || j == hi.y)
                        {
                            std::array<Numeric::StencilType, AMREX_SPACEDIM> sten = Numeric::GetStencil(i, j, k, bx);
                            Set::Matrix gradu = Numeric::Gradient(u, i, j, k, dx, sten);

                            Set::Vector U(u(i,j,k,0),u(i,j,k,1));
                            Set::Vector ret = m_bc(U, gradu, dw(i,j,k), i, j, k, bx);

                            for (int p = 0; p < 2; p++)
                                rhs(i,j,k,p) = b(i,j,k,p) - ret(p);
                        }
                        else
                        {
                            for (int p = 0; p < 2; p++)
                                rhs(i,j,k,p) = b(i,j,k,p) - 
                                    ((dw(i + 1, j, k)(p, 0) - dw(i - 1, j, k)(p, 0)) / 2. / DX(0) +
                                     (dw(i, j + 1, k)(p, 1) - dw(i, j - 1, k)(p, 1)) / 2. / DX(1));
                        }

                        #elif AMREX_SPACEDIM == 3

                        if (i == lo.x || i == hi.x || j == lo.y || j == hi.y || k == lo.z || k == hi.z)
                        {
//                            std::array<Numeric::StencilType, AMREX_SPACEDIM> sten = Numeric::GetStencil(i, j, k, bx);
//                            Set::Matrix gradu = Numeric::Gradient(u, i, j, k, dx, sten);
//
//                            Set::Vector U(u(i,j,k,0),u(i,j,k,1),u(i,j,k,2));
//                            Set::Vector ret = m_bc(U, gradu, dw(i,j,k), i, j, k, bx);
//
//                            for (int p = 0; p < 3; p++)
//                                rhs(i,j,k,p) = b(i,j,k,p) - ret(p);
                        }
                        else
                        {
                            for (int p = 0; p < 3; p++)
                                rhs(i,j,k,p) = b(i,j,k,p) - 
                                    ((dw(i + 1, j, k)(p, 0) - dw(i - 1, j, k)(p, 0)) / 2. / DX(0) +
                                     (dw(i, j + 1, k)(p, 1) - dw(i, j - 1, k)(p, 1)) / 2. / DX(1) +
                                     (dw(i, j, k + 1)(p, 2) - dw(i, j, k - 1)(p, 2)) / 2. / DX(2));
                        }
                        #endif
                    });
                }
                Util::RealFillBoundary(a_model_mf[lev],m_elastic.Geom(lev));
            }
    }

public:
    Set::Scalar solve (const amrex::Vector<MultiFab*>& a_u_mf, 
                       const amrex::Vector<MultiFab*>& a_b_mf,
                       amrex::Vector<amrex::FabArray<amrex::BaseFab<T> > > &a_model_mf,
                       Real a_tol_rel, Real a_tol_abs, const char* checkpoint_file = nullptr)
    {
        //Operator::Elastic<T> &elastic = dynamic_cast<Operator::Elastic<T> &>(linop);

        amrex::Vector<amrex::MultiFab> dsol_mf, rhs_mf;
        amrex::Vector<amrex::FabArray<amrex::BaseFab<Set::Matrix> > > dw_mf;
        dsol_mf.resize(a_u_mf.size());
        dw_mf.resize(a_u_mf.size());
        rhs_mf.resize(a_u_mf.size());
        for (int lev = 0; lev < a_u_mf.size(); lev++)
        {
            dsol_mf[lev].define(a_u_mf[lev]->boxArray(),
                                a_u_mf[lev]->DistributionMap(),
                                a_u_mf[lev]->nComp(), a_u_mf[lev]->nGrow());
            dw_mf[lev].define(a_b_mf[lev]->boxArray(),
                              a_b_mf[lev]->DistributionMap(),
                              1, a_b_mf[lev]->nGrow());
            rhs_mf[lev].define(a_b_mf[lev]->boxArray(),
                               a_b_mf[lev]->DistributionMap(),
                               a_b_mf[lev]->nComp(), a_b_mf[lev]->nGrow());
            
            dsol_mf[lev].setVal(0.0);
            dw_mf[lev].setVal(Set::Matrix::Zero());
            
            amrex::MultiFab::Copy(rhs_mf[lev], *a_b_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
        }

//        for (int lev = 0; lev < a_b_mf.size(); ++lev)
//            m_elastic.GetBC().Init(&(rhs_mf[lev]),m_elastic.Geom(lev),false);

        for (int nriter = 0; nriter < m_nriters; nriter++)
        {
            if (m_verbose > 0 && nriter < m_nriters) Util::Message(INFO, "Newton Iteration ", nriter+1, " of ", m_nriters);

            for (int lev = 0; lev < a_b_mf.size(); ++lev)
                m_elastic.GetBC().Init(&(rhs_mf[lev]),m_elastic.Geom(lev),nriter>0);

            prepareForSolve(a_u_mf, a_b_mf, GetVecOfPtrs(rhs_mf), dw_mf, a_model_mf);

            m_elastic.SetModel(a_model_mf);

            if (nriter == m_nriters) break;
            
            Solver::Nonlocal::Linear::solve(GetVecOfPtrs(dsol_mf), GetVecOfConstPtrs(rhs_mf), a_tol_rel, a_tol_abs,checkpoint_file);

            Set::Scalar cornorm = 0, solnorm = 0;
            for (int lev = 0; lev < dsol_mf.size(); ++lev)
            {
                for (int comp = 0; comp < AMREX_SPACEDIM; comp++)
                {
                    Set::Scalar tmpcornorm = dsol_mf[lev].norm0(comp,0);
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
                amrex::MultiFab::Add(*a_u_mf[lev], dsol_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
        }

        return 0.0;
    }

    Set::Scalar solve(amrex::Vector<std::unique_ptr<amrex::MultiFab> > & a_u_mf,
                      amrex::Vector<std::unique_ptr<amrex::MultiFab> > & a_b_mf,
                      amrex::Vector<amrex::FabArray<amrex::BaseFab<T> > > &a_model,
                      Set::Scalar a_tol_rel, Set::Scalar a_tol_abs,
                      const char *checkpoint_file = nullptr)
    {
        return solve(GetVecOfPtrs(a_u_mf),
                     GetVecOfPtrs(a_b_mf),a_model,a_tol_rel,a_tol_abs,checkpoint_file);
    }

    void compResidual (const amrex::Vector<MultiFab*>& a_res_mf, 
                       const amrex::Vector<MultiFab*>& a_u_mf,
                       const amrex::Vector<MultiFab*>& a_b_mf,
                       amrex::Vector<amrex::FabArray<amrex::BaseFab<T> > > &a_model_mf)
    {
        amrex::Vector<amrex::FabArray<amrex::BaseFab<Set::Matrix> > > dw_mf;
        dw_mf.resize(a_u_mf.size());
        for (int lev = 0; lev < a_u_mf.size(); lev++)
        {
            dw_mf[lev].define(a_b_mf[lev]->boxArray(),
                              a_b_mf[lev]->DistributionMap(),
                              1, a_b_mf[lev]->nGrow());
            dw_mf[lev].setVal(Set::Matrix::Zero());
        }
        
        //for (int lev = 0; lev < a_b_mf.size(); ++lev)
        //    m_elastic.GetBC().Init(a_b_mf[lev],m_elastic.Geom(lev),true);
            
        prepareForSolve(a_u_mf, a_b_mf, a_res_mf, dw_mf, a_model_mf);
        
    }
    void compResidual(amrex::Vector<std::unique_ptr<amrex::MultiFab> > & a_res_mf,
                      amrex::Vector<std::unique_ptr<amrex::MultiFab> > & a_u_mf,
                      amrex::Vector<std::unique_ptr<amrex::MultiFab> > & a_b_mf,
                      amrex::Vector<amrex::FabArray<amrex::BaseFab<T> > > &a_model_mf)
    {
        compResidual(GetVecOfPtrs(a_res_mf),GetVecOfPtrs(a_u_mf),GetVecOfPtrs(a_b_mf),a_model_mf);
    }                      

private:
    int m_nriters = 1;
    Operator::Elastic<T> &m_elastic;
    BC::Operator::Elastic<T> &m_bc;

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

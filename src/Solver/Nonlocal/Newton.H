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
        dsol_mf.resize(a_sol_mf.size());
        rhs_mf.resize(a_sol_mf.size());
        for (int lev = 0; lev < a_sol_mf.size(); lev++)
        {
            dsol_mf[lev].reset(new amrex::MultiFab(a_sol_mf[lev]->boxArray(),
                                                   a_sol_mf[lev]->DistributionMap(),
                                                   a_sol_mf[lev]->nComp(), a_sol_mf[lev]->nGrow()));
            rhs_mf[lev].reset(new amrex::MultiFab(a_rhs_mf[lev]->boxArray(),
                                                  a_rhs_mf[lev]->DistributionMap(),
                                                  a_rhs_mf[lev]->nComp(), a_rhs_mf[lev]->nGrow()));

            dsol_mf[lev]->setVal(0.0);
            
            //rhs_mf[lev]->setVal(0.0); Util::Warning(INFO,"Should not need to do this");
            
            amrex::MultiFab::Copy(*rhs_mf[lev], *a_rhs_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
        }

        for (int nriter = 0; nriter < (copyrhs ? m_nriters+1 : m_nriters); nriter++)
        {
            //Util::Message(INFO, "Newton Iteration ", nriter+1, " of ", m_nriters);

            for (int lev = 0; lev < a_rhs_mf.size(); ++lev)
            {
                amrex::Box domain(linop.Geom(lev).Domain());
                domain.convert(amrex::IntVect::TheNodeVector());
                const Set::Scalar *dx = linop.Geom(lev).CellSize();
                Set::Vector DX(linop.Geom(lev).CellSize());

//
//  TODO: re-enable this. Don't need it now since BCs are already homogenized.
//
//                // homogenize boundary conditions for subsequent iterations (not the first)
//                if (nriter > 0)
//                {
//                    for (MFIter mfi(a_model[lev], amrex::TilingIfNotGPU); mfi.isValid(); ++mfi)
//                    {
//                        amrex::Box bx = mfi.tilebox();
//
//                        amrex::Array4<Set::Scalar> const &rhs = rhs_mf[lev]->array(mfi);
//
//                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
//                            if (i == domain.loVect()[0] || j == domain.loVect()[1] || k == domain.loVect()[2] ||
//                                i == domain.hiVect()[0] || j == domain.hiVect()[1] || k == domain.hiVect()[2])
//                            {
//                                rhs(i, j, k, 0) = 0.0;
//                                rhs(i, j, k, 1) = 0.0;
//                                rhs(i, j, k, 2) = 0.0;
//                            }
//                        });
//                    }
//                }

                for (MFIter mfi(a_model[lev], amrex::TilingIfNotGPU); mfi.isValid(); ++mfi)
                {
                    amrex::Box bx = mfi.growntilebox(2);
                    amrex::Box domainint = domain;
                    domainint.grow(-1); // shrink so that we are only on the interior of the domain
                    bx = bx & domain;

                    amrex::Array4<Set::Scalar> const &u = a_sol_mf[lev]->array(mfi);
                    //amrex::Array4<Set::Scalar> const &f = a_rhs_mf[lev]->array(mfi);
                    amrex::Array4<Set::Scalar> const &graddw = rhs_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const &b = a_rhs_mf[lev]->array(mfi);
                    amrex::Array4<T> const &model = a_model[lev].array(mfi);

                    // Set model internal dw and ddw.
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        std::array<Numeric::StencilType, AMREX_SPACEDIM> sten = Numeric::GetStencil(i, j, k, bx);

                        Set::Matrix gradu = Numeric::Gradient(u, i, j, k, dx, sten);

                        if (model(i,j,k).kinvar == Model::Solid::KinematicVariable::gradu)
                        {
                            model(i, j, k).dw = model(i, j, k).DW(gradu);
                            model(i, j, k).ddw = model(i, j, k).DDW(gradu);
                        }
                        else if (model(i,j,k).kinvar == Model::Solid::KinematicVariable::epsilon)
                        {
                            Set::Matrix eps = 0.5 * (gradu + gradu.transpose());
                            model(i, j, k).dw = model(i, j, k).DW(eps);
                            model(i, j, k).ddw = model(i, j, k).DDW(eps);
                        }
                        else if (model(i,j,k).kinvar == Model::Solid::KinematicVariable::F)
                        {
                            Set::Matrix F = gradu + Set::Matrix::Identity();
                            model(i, j, k).dw = model(i, j, k).DW(F);
                            model(i, j, k).ddw = model(i, j, k).DDW(F);
                        }
                    });

                    if (nriter > 0)
                    rhs_mf[lev]->setVal(0.0);

                    //bx = bx & domainint;
                    //bx = mfi.growntilebox(1);
                    bx = mfi.validbox() & domainint;
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        for (int p = 0; p < AMREX_SPACEDIM; p++)
                        {
                            #if AMREX_SPACEDIM == 2
                            graddw(i,j,k,p) = b(i,j,k,p) -   // exact same as above
                                ((model(i + 1, j, k).dw(p, 0) - model(i - 1, j, k).dw(p, 0)) / 2. / DX(0) +
                                 (model(i, j + 1, k).dw(p, 1) - model(i, j - 1, k).dw(p, 1)) / 2. / DX(1));
                            #elif AMREX_SPACEDIM == 3
                            graddw(i,j,k,p) = b(i,j,k,p) -   // exact same as above
                                ((model(i + 1, j, k).dw(p, 0) - model(i - 1, j, k).dw(p, 0)) / 2. / DX(0) +
                                 (model(i, j + 1, k).dw(p, 1) - model(i, j - 1, k).dw(p, 1)) / 2. / DX(1) +
                                 (model(i, j, k + 1).dw(p, 2) - model(i, j, k - 1).dw(p, 2)) / 2. / DX(2));
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
            elastic.SetModel(a_model);


            if (nriter == m_nriters) break;
            
            Solver::Nonlocal::Linear::solve(dsol_mf,rhs_mf, a_tol_rel, a_tol_abs,checkpoint_file);

            for (int lev = 0; lev < dsol_mf.size(); ++lev)
            {
                //amrex::IntVect::TheUnitVector()
                //Util::Message(INFO,"DSol Norm Level ", lev, " = (",dsol_mf[lev]->norm2(0,elastic.Geom(lev).periodicity()),",", dsol_mf[lev]->norm2(1,elastic.Geom(lev).periodicity()) ,",", dsol_mf[lev]->norm2(2,elastic.Geom(lev).periodicity()),")");
                //Util::Message(INFO,"RHS  Norm Level ", lev, " = (",rhs_mf[lev]->norm2(0) ,",", rhs_mf[lev]->norm2(1)  ,",", rhs_mf[lev]->norm2(2),")");
            }

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
        
         //value.m_nriters
    }

};
} // namespace Nonlocal
} // namespace Solver


#endif